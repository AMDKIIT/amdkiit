MODULE wfn_initialize

  CONTAINS
  SUBROUTINE wfn_initialize_orb
  USE kinds
  USE mympi, ONLY : MPI_GlobSumC,MPI_GlobSumR2s,MPI_RBcast,MPI_GlobSumR2
  USE system_data_types
  USE pseudopotential
  IMPLICIT NONE
  REAL(KIND=dp) :: sfc
  INTEGER :: iii,il,istate,ig,is,lxx,natst,ish,l,ll,iv,ly,lpp(5),I,j,IAORB,IAT,ia,ixx,ORB,IST,NN,IOPT,ik,q,NN1
  REAL(KIND=dp) :: norm
  COMPLEX*16, DIMENSION(:,:),ALLOCATABLE :: work
  COMPLEX*16::ci,CATOM(ngpw_l,NATTOT)
  DATA       LPP /0,1,4,9,16/ 
  REAL(kind=dp) vol!,cc 
  REAL(kind=dp) DDOT,FC(20)!,SCR(73083)
  real(kind=dp), DIMENSION(:),ALLOCATABLE :: ggng,n_cc
  
  EXTERNAL DDOT
  REAL(kind=dp) OVLAP_MAT(NATTOT,NATTOT),XMATAT(NATTOT,NATTOT),XSMAT(NATTOT,NATTOT),WORKAT(3*NATTOT),EIVAT(NATTOT)

  IF(nstate<=0.OR.ngpw_l<=0)STOP 'Error! nstate/ngw is not initialized'
  !IF(.NOT.ASSOCIATED(psi))THEN
  !ALLOCATE(psi(ngpw_l,nstate)) 
  !ELSE
  !IF(SIZE(psi,1)/=ngpw_l)STOP 'Error! psi - array size mismatch'
  !IF(SIZE(psi,2)/=nstate)STOP 'Error! psi - array size mismatch'
  !END IF

  ALLOCATE(C_0(NGPW_l,NSTATE),C2(NGPW_l,nstate),ggng(nspln),n_cc(ngpw_l))    

  DO IL =1,nspln
    ggng(il)=(il-1)*(Gcutoff%pw)/DBLE(nspln-1)
  ENDDO
  
    C_0(:,:)=(0._dp,0._dp)
   ! psi(:,:)=(0._dp,0._dp)
    OVLAP_MAT(:,:)=0._dp
    CATOM(:,:)=(0.0d0,0.0d0)
    IAORB=1
    IAT=0
    ORB=0
      DO IS=1,sp_t
        DO IA=1,atom_p_sp(IS)
           IAT=IAT+1
           VOL=1.D0/SQRT(cell_volume)
           LXX=0
          DO ISH=1,NSHELL(IS)
             L=LSHELL(ISH,IS)
             CI=(0.0D0,1.0D0)**L
             LL=2*L+1
             CALL spline_inter(nspln,ggng,cat(1,1,ish,is),hg(1),n_cc(1),NGPW_l)
            DO IV=1,LL
               LXX=LXX+1
               ORB=ORB+1
               LY=LPP(L+1)+IV
              DO IG=1,NGPW_l
                 CATOM(IG,orb)=CI*YLG(LY,IG,gvec(1:3,1))*n_CC(ig)*EIGR(IG,IAT)*vol
              ENDDO
               FOC(LXX)=OC(ISH,IS)/DBLE(LL)
            ENDDO
          ENDDO

          NATST=LXX
          DO IXX=IAORB,IAORB+NATST-1
        
            SFC=DOTP(NGPW_l,CATOM(1,IXX),CATOM(1,IXX))
            CALL MPI_GlobSumR2s(SFC)
            IF(SFC.EQ.0.D0) THEN
              WRITE(6,'(A,A,I5,A)') ' ATRHO|',&
              &' THE NORM OF ATOMIC ORBITAL (',IXX,') IS NULL'
               STOP
            ELSE
               SFC=1.D0/SQRT(SFC)
            ENDIF
            CALL DSCAL(2*ngpw_l,SFC,CATOM(1,IXX),1)
       
          ENDDO
          IAORB=IAORB+NATST
        ENDDO
      ENDDO
      DO I=1,20
      FC(I)=1.0d0
      ENDDO
      CALL GS_DISORTHO(NATTOT,CATOM) 
      
      IST=1
     if(l_upf)then
      DO IS=1,sp_t
        DO IA=1,atom_p_sp(IS)
   
           CALL UPF_fnl(catom(:,ist:ist+numaor(is)-1),NUMAOR(is)) 
           C2(1:NGPW_L:1,1:NUMAOR(IS):1)=(0._dp,0._dp)
           CALL VPSI(CATOM(:,ist:ist+numaor(is)-1),C2,FC,NUMAOR(IS))
    
           CALL UPF_FNONLOC(C2,NUMAOR(IS),FC)
     
           CALL OVLAP_D(NGPW_l,NATTOT,NUMAOR(IS),OVLAP_MAT(1,IST),CATOM,C2)
           IST=IST+NUMAOR(IS)
        ENDDO
      ENDDO
      else
      DO IS=1,sp_t
        DO IA=1,atom_p_sp(IS)
           CALL fnl(catom(:,ist:ist+numaor(is)-1),NUMAOR(is))
           C2(1:NGPW_L:1,1:NUMAOR(IS):1)=(0._dp,0._dp)
           CALL VPSI(CATOM(:,ist:ist+numaor(is)-1),C2,FC,NUMAOR(IS))
           CALL FNONLOC(C2,NUMAOR(IS),FC)
           CALL OVLAP_D(NGPW_l,NATTOT,NUMAOR(IS),OVLAP_MAT(1,IST),CATOM,C2)!!IST),CATOM,C2)
           IST=IST+NUMAOR(IS)
        ENDDO
        ENDdo
      endif
    CALL DSCAL(NATTOT*NATTOT,-1.D0,OVLAP_MAT(1,1),1)

    CALL OVLAP_D(NGPW_l,NATTOT,NATTOT,XSMAT,CATOM,CATOM)
    CALL MPI_GlobSumR2(OVLAP_MAT,NATTOT*NATTOT)
    CALL MPI_GlobSumR2(XSMAT,NATTOT*NATTOT)
    XMATAT(:,:)=0._dp
   IOPT=1
   IF(icpu==1)CALL DSYGVX(IOPT,OVLAP_MAT,NATTOT,XSMAT,NATTOT,EIVAT,XMATAT,NATTOT,NATTOT,WORKAT,3*NATTOT)
   CALL MPI_RBcast(XMATAT,NATTOT*NATTOT)
!================================================
   IF(LOPEN_SHELL)THEN
     NN=MIN(NEL_UP,NATTOT)
     IF(NN.GT.0 .AND. NGPW_l.GT.0)THEN
       CALL DGEMM('N','N',2*NGPW_l,NN,NATTOT,1.0D0,CATOM(1,1),2*NGPW_l,XMATAT,NATTOT,0.0D0,C_0(1,1),2*NGPW_l)
     ENDIF
     NN=MIN(NEL_DOWN,NATTOT)
     NN1=NEL_UP+1
     IF(NN.GT.0 .AND. NGPW_l.GT.0)THEN
       CALL DGEMM('N','N',2*NGPW_l,NN,NATTOT,1.0D0,CATOM(1,1),2*NGPW_l,XMATAT,NATTOT,0.0D0,C_0(1,NN1),2*NGPW_l)
     ENDIF
   ELSE
     NN=MIN(NSTATE,NATTOT)

     IF(NN.GT.0 .AND. NGPW_l.GT.0)THEN
      CALL DGEMM('N','N',2*NGPW_l,NN,NATTOT,1.0D0,CATOM(1,1),2*NGPW_l,XMATAT,NATTOT,0.0D0,C_0(1,1),2*NGPW_l)
     ENDIF
   ENDIF
!================================================
DEALLOCATE(GGNG,n_cc,cat,oc,NUMAOR)
RETURN
END SUBROUTINE wfn_initialize_orb


      SUBROUTINE DSYGVX(IOPT,A,LDA,B,LDB,W,Z,LDZ,N,AUX,NAUX)
!     ==--------------------------------------------------------------==
!     == GENERALIZED DIAGONALIZATION ROUTINE:                         ==
!     ==--------------------------------------------------------------==
!     == GENERALIZED DIAGONALIZATION ROUTINE:                         ==
!     == GENERALIZED DIAGONALIZATION ROUTINE:                         ==
!     == FOLLOW THE ESSL CONVENTION                                   ==
!     ==--------------------------------------------------------------==
      USE KINDS
      IMPLICIT NONE
!     Arguments
      INTEGER  IOPT,LDA,LDB,LDZ,N,NAUX
      REAL(kind=dp)   A(LDA,N),B(LDB,N),W(N),Z(LDZ,N),AUX(NAUX)
!     Variables
      INTEGER  INFO,I
      !WRITE(*,*)LDA,LDB,LDZ,N,NAUX
      !write(*,*)A(LDA,N),B(LDB,N),W(N),Z(LDZ,N),AUX(NAUX)
!     ==--------------------------------------------------------------==
      IF(IOPT.EQ.1) THEN
        INFO=0
        CALL DSYGV(1,'V','U',N,A,LDA,B,LDB,W,AUX,NAUX,INFO)
        DO I=1,N
          CALL DCOPY(N,A(1,I),1,Z(1,I),1)
        ENDDO
      ELSE
        WRITE(6,*)"DSYGVX MISSING LIBRARY OPTION"
        STOP
      ENDIF
      IF(INFO.NE.0) THEN
        WRITE(*,'(/," DSYGVX| INFO=",I5)') INFO
        STOP
      ENDIF
!     ==--------------------------------------------------------------==
      RETURN
      END



!      SUBROUTINE SUMMAT(A,NSTAT)!,SCR,LSCR)
!      USE system_data_types
!      USE mympi
!      !USE utils
!      USE kinds
!      IMPLICIT NONE
!      INTEGER NSTAT,LSCR
!      REAL(kind=dp)  A(NSTAT,NSTAT)!,SCR(LSCR)
!      REAL(kind=dp), DIMENSION(:),ALLOCATABLE :: SCR
!      INTEGER N2,NTR,K,I,J
!      ==--------------------------------------------------------------==
!      == GLOBAL SUMMATION OF A SYMMETRIC MATRIX                       ==
!      ==--------------------------------------------------------------==
!      IF(ncpu.LE.1) RETURN
!      N2 = NSTAT*NSTAT
!      NTR = (NSTAT*(NSTAT+1))/2
!      ALLOCATE(SCR(NTR))
!      SCR=0.0
!      DO I=1,NSTAT
!        DO J=I,NSTAT
!          K=I+(J*(J-1))/2
!          SCR(K)=A(I,J)
!          !if(icpu==1)write(6,*)i,j,"kk",I,J,A(I,J)
!        ENDDO
!      ENDDO
!      CALL MPI_GlobSumR2(A,N2)      
!      RETURN
!      END SUBROUTINE



SUBROUTINE OVLAP(NGPW_l,N1,N2,A,C1,C2)
USE kinds
!    ==--------------------------------------------------------------==
!    ==         COMPUTES THE OVERLAP MATRIX A = < C1 | C2 >          ==
!    ==--------------------------------------------------------------==
      IMPLICIT NONE
      INTEGER NGPW_l,N1,N2
      COMPLEX*16 C1(NGPW_l,*),C2(NGPW_l,*)
      REAL(kind=dp) A(N1,N2)

      A(:,:)=0._dp
    
      CALL DGEMM('T','N',N1,N2,2*NGPW_l,2.0D0,C1(1,1),2*NGPW_l,C2(1,1),2*NGPW_l,0.0D0,A(1,1),N1)
      CALL DGER(N1,N2,-1.0D0,C1(1,1),2*NGPW_l,C2(1,1),2*NGPW_l,A(1,1),N1)
      
END SUBROUTINE

SUBROUTINE OVLAP_D(NG_PW,N1,N2,A,CA,CB)
USE KINDS
use SYSTEM_DATA_TYPES!gvectors, ONLY: G0_STAT
!    ==--------------------------------------------------------------==
!    ==         COMPUTES THE OVERLAP MATRIX A = < C1 | C2 >          ==
!    ==--------------------------------------------------------------==
      IMPLICIT NONE
      INTEGER NG_PW,N1,N2,NN1
      COMPLEX*16 CA(NG_PW,*),CB(NG_PW,*)
      REAL(kind=dp) A(N1,N2)
!====================================================
      IF(LOPEN_SHELL.AND.NLSD.EQ.2)THEN
        A(:,:)=0._dp
        NN1=NEL_UP+1
!alpha electron
        CALL DGEMM('T','N',NEL_UP,NEL_UP,2*NG_PW,2.0D0,CA(1,1),2*NG_PW,CB(1,1),2*NG_PW,0.0D0,A(1,1),N1)
        if(g0_stat)CALL DGER(NEL_UP,NEL_UP,-1.0D0,CA(1,1),2*NG_PW,CB(1,1),2*NG_PW,A(1,1),N1)
!beta electron
        CALL DGEMM('T','N',NEL_DOWN,NEL_DOWN,2*NG_PW,2.0D0,CA(1,NN1),2*NG_PW,CB(1,NN1),2*NG_PW,0.0D0,A(NN1,NN1),N1)
        IF(g0_stat)CALL DGER(NEL_DOWN,NEL_DOWN,-1.0D0,CA(1,NN1),2*NG_PW,CB(1,NN1),2*NG_PW,A(NN1,NN1),N1)
      ELSE
        A(:,:)=0._dp
        CALL DGEMM('T','N',N1,N2,2*NG_PW,2.0D0,CA(1,1),2*NG_PW,CB(1,1),2*NG_PW,0.0D0,A(1,1),N1)
        if(g0_stat)CALL DGER(N1,N2,-1.0D0,CA(1,1),2*NG_PW,CB(1,1),2*NG_PW,A(1,1),N1)
      ENDIF
!====================================================
      !A(:,:)=0._dp
      !CALL DGEMM('T','N',N1,N2,2*NG_PW,2.0D0,CA(1,1),2*NG_PW,CB(1,1),2*NG_PW,0.0D0,A(1,1),N1)
      !if(g0_stat) CALL DGER(N1,N2,-1.0D0,CA(1,1),2*NG_PW,CB(1,1),2*NG_PW,A(1,1),N1)
END SUBROUTINE




SUBROUTINE vpsi(C_ATOM,C2_nl,FC,NOSTATE)
USE system_data_types
USE constants
USE kinds
USE fft_interface, ONLY: prepare_fft, fft_forward,fft_backward
IMPLICIT NONE

COMPLEX*16 C_ATOM(NGPW_l,nostate)!,C2(NGPW,1)
COMPLEX*16,INTENT(INOUT):: C2_nl(NGPW_l,nostate)
complex*16, pointer  :: PSY(:)
INTEGER IA,IG,IR,I,J,IB,IS1,IS2,IBB,NSTA,NOSTATE
COMPLEX*16 FP,FM
REAL(kind=dp)     FI,FIP1,FC(NOSTATE)
ALLOCATE(PSY(MAX_FFT))

            dO IA=1,NOSTATE,2
                psy=(0._dp,0._dp)
                NSTA=MIN((NOSTATE-IA+2)/2,1)!NGROUP)
               DO IB=1,NSTA
                I=IA+2*(IB-1)
                IBB=(IB-1)*(Ggridinfo(6,icpu))*nrlead(1)!
                IS1 = I
                IS2 = I+1
                
                   IF(IS2.GT.NOSTATE) THEN
                    DO IG=1,NGPW_l
                    PSY(map_grid1d_p_pw(IG)+IBB)=C_ATOM(IG,IS1)
                    PSY(map_grid1d_m_pw(IG)+IBB)=DCONJG(C_ATOM(IG,IS1))
                    
                    ENDDO
                    IF(G0_stat)PSY(map_grid1d_p_pw(1)+IBB)=C_ATOM(1,IS1)
                ELSE 
                    DO IG=1,NGPW_l
                    PSY(map_grid1d_p_pw(IG)+IBB)=C_ATOM(IG,IS1)+UIMAG*C_ATOM(IG,IS2)
                    PSY(map_grid1d_m_pw(IG)+IBB)=DCONJG(C_ATOM(IG,IS1))+UIMAG*DCONJG(C_ATOM(IG,IS2))
                    
                    

                    ENDDO
                    if(g0_stat)PSY(map_grid1d_p_pw(1)+IBB)=C_ATOM(1,IS1)+UIMAG*C_ATOM(1,IS2)
                ENDIF
               ENDDO !IB
!   Transform the wavefunction to real space
                CALL fft_backward(PSY)
       !Sudhir DBG
        IF(LOPEN_SHELL.AND.NLSD.EQ.2)THEN
          IF(IA.EQ.NEL_UP)THEN
            DO IG=1,NNR1
              PSY(IG)=V_POT(IG,1)*DREAL(PSY(IG))+UIMAG*V_POT(IG,2)*DIMAG(PSY(IG))
            ENDDO
          ELSEIF(IA.GT.NEL_UP)THEN
            DO IG=1,NNR1
              PSY(IG)=V_POT(IG,2)*PSY(IG)
            ENDDO
          ELSE
            DO IG=1,NNR1
              PSY(IG)=V_POT(IG,1)*PSY(IG)
            ENDDO
          ENDIF
        ELSE
          DO IG=1,NNR1
            PSY(IG)=V_POT(IG,NLSD)*PSY(IG)
           !write(104,*)v_pot(ig,nlsd),ig
          ENDDO
        ENDIF
      !Sudhir DBG

       !         DO IG=1,NNR1
       !         PSY(IG)=V_POT(IG)*PSY(IG)
       !         ENDDO
              !if(icpu==1)write(3335,*)"VPOT",V_POT(1)!,PSY(1),NNR1
!   Back transform to reciprocal space the product V.
                CALL fft_forward(PSY)
              !if(icpu==1)write(3335,*)"#",PSY(1)
              DO IB=1,NSTA
                I=IA+2*(IB-1)
                IBB=(IB-1)*(IB-1)*(Ggridinfo(6,icpu))*nrlead(1)
                IS1 = I
                IS2 = I+1
                FI=FC(IS1)*0.5D0
                IF(FI.EQ.0.D0) FI=1.D0
                FIP1=0.D0
                IF(IS2.LE.NOSTATE) FIP1=FC(IS2)*0.5D0
                IF(FIP1.EQ.0.D0) FIP1=1.D0


                DO IG=1,NGPW_l
                FP=PSY(map_grid1d_p_pw(IG)+IBB)+PSY(map_grid1d_m_pw(IG)+IBB)
                FM=PSY(map_grid1d_p_pw(IG)+IBB)-PSY(map_grid1d_m_pw(IG)+IBB)
                    C2_nl(IG,IS1)=-FI*((twopibya2*HG(IG))*C_ATOM(IG,IS1)+DCMPLX(DREAL(FP),DIMAG(FM)))
                IF(IS2.LE.NOSTATE)C2_nl(IG,IS2)=-FIP1*((twopibya2*HG(IG))*C_ATOM(IG,IS2)+DCMPLX(DIMAG(FP),-DREAL(FM)))
               ! write(103,*)c2_nl(ig,is1),ig,is1
                ENDDO
              ENDDO
           ENDDO!IA
          DEALLOCATE(PSY)

END SUBROUTINE vpsi

      FUNCTION YLG(L,IG,GK)
      USE system_data_types
      USE constants
      USE kinds
      IMPLICIT NONE

      INTEGER L,IG
      REAL(kind=dp)  GK(3,*)
      REAL(kind=dp)  YLG
      REAL(kind=dp)  YLM
        CALL compute_yL(L,1,HG(ig),GK(1,IG),YLM)
        YLG=YLM
      RETURN
      END

      FUNCTION DOTP(N,CA,CB)
!C     ==--------------------------------------------------------------==
      USE kinds
      USE system_data_types

      IMPLICIT NONE
!C     Arguments
      INTEGER N
      COMPLEX*16 CA(N),CB(N)
      REAL(kind=dp) DOTP,A(2*N),B(2*N)
!C     Variables
      REAL(kind=dp) DDOT
      INTEGER I
      EXTERNAL DDOT
!C     ==--------------------------------------------------------------==
      !write(17,*)CA(1),CB(1),N
      IF(N.EQ.0) THEN
        DOTP=0.D0
        RETURN
      ELSE
        DO I=1,2*N,2
        A(I)=DREAL(CA((I/2)+1))
        A(I+1)=DIMAG(CA((I/2)+1))
        B(I)=DREAL(CB((I/2)+1))
        B(I+1)=DIMAG(CB((I/2)+1))
        ENDDO
        !write(35,*)N,A(1),A(2),A(3),A(4)
      ENDIF

      IF(g0_stat) THEN
        DOTP=A(1)*B(1)
      ELSE
        DOTP=2.0D0*(A(1)*B(1)+A(2)*B(2))
      ENDIF

      IF(N.GT.1) THEN
      DOTP=DOTP+2.0D0*DDOT(2*N-2,A(3),1,B(3),1)
      ENDIF
      DO I=1,2*N
      !write(35,*)N,A(1),A(2),A(3),A(4)
      ENDDO

!C     ==--------------------------------------------------------------==
      RETURN
      END

    subroutine update_wf
    USE kinds
!    USE mympi, ONLY : MPI_GlobSumC
    USE system_data_types
    USE interpolgrid
    USE constants
!    USE utils
    IMPLICIT NONE
    INTEGER K,IG

        DO K=1,NSTATE
          DO IG=1,NGPW_L
           !write(33,*)IG,C_0(IG,K),C2(IG,K)
            C_0(IG,K)=C_0(IG,K)+LAMBDA*C2(IG,K)
          ENDDO
        ENDDO
   

      end subroutine update_wf

      SUBROUTINE ortho_gs(CP)
      USE kinds
      USE system_data_types
      USE mympi
      IMPLICIT NONE
      REAL(kind=dp),DIMENSION(:,:),ALLOCATABLE :: SMAT
      COMPLEX*16 CP(NGPW_L,*)
      INTEGER    ISUB
      IF(NSTATE.LE.0) RETURN
      !===========================================================
            IF(LOPEN_SHELL.AND.NLSD.EQ.2)THEN
!alpha electron
        ALLOCATE(SMAT(NEL_UP,NEL_UP))
        SMAT=0.0
        IF(NGPW_L.GT.0)CALL DSYRK('U','T',NEL_UP,2*NGPW_L,2.D0,CP(1,1),2*NGPW_L,0.D0,SMAT,NEL_UP)
        IF(g0_stat)CALL DGER(NEL_UP,NEL_UP,-1.0D0,CP(1,1),2*NGPW_L,CP(1,1),2*NGPW_L,SMAT,NEL_UP)
        CALL MPI_GlobSumR2(SMAT,NEL_UP*NEL_UP)
        CALL UINV('U',SMAT,NEL_UP,NEL_UP)
        IF(NGPW_L.GT.0)CALL DTRMM('R','U','N','N',2*NGPW_L,NEL_UP,1.D0,SMAT,NEL_UP,CP(1,1),2*NGPW_L)
        IF(g0_stat) CALL ZCLEAN(CP(1,1),NEL_UP,NGPW_L)!todo
        DEALLOCATE(SMAT)
!beta electron
        ALLOCATE(SMAT(NEL_DOWN,NEL_DOWN))
        SMAT=0.0
        IF(NGPW_L.GT.0)CALL DSYRK('U','T',NEL_DOWN,2*NGPW_L,2.D0,CP(1,NEL_UP+1),2*NGPW_L,0.D0,SMAT,NEL_DOWN)
        IF(g0_stat)CALL DGER(NEL_DOWN,NEL_DOWN,-1.0D0,CP(1,NEL_UP+1),2*NGPW_L,CP(1,NEL_UP+1),2*NGPW_L,SMAT,NEL_DOWN)
        CALL MPI_GlobSumR2(SMAT,NEL_DOWN*NEL_DOWN)
        CALL UINV('U',SMAT,NEL_DOWN,NEL_DOWN)
        IF(NGPW_L.GT.0)CALL DTRMM('R','U','N','N',2*NGPW_L,NEL_DOWN,1.D0,SMAT,NEL_DOWN,CP(1,NEL_UP+1),2*NGPW_L)
        IF(g0_stat) CALL ZCLEAN(CP(1,NEL_UP+1),NEL_DOWN,NGPW_L)!todo
      ELSE
        ALLOCATE(SMAT(NSTATE,NSTATE))
        SMAT=0.0
        IF(NGPW_L.GT.0)CALL DSYRK('U','T',NSTATE,2*NGPW_L,2.D0,CP,2*NGPW_L,0.D0,SMAT,NSTATE)
        IF(g0_stat)CALL DGER(NSTATE,NSTATE,-1.0D0,CP,2*NGPW_L,CP,2*NGPW_L,SMAT,NSTATE)
        CALL MPI_GlobSumR2(SMAT,NSTATE*NSTATE)
        CALL UINV('U',SMAT,NSTATE,NSTATE)
        IF(NGPW_L.GT.0)CALL DTRMM('R','U','N','N',2*NGPW_L,NSTATE,1.D0,SMAT,NSTATE,CP,2*NGPW_L)
        IF(g0_stat) CALL ZCLEAN(CP,NSTATE,NGPW_L)!todo
      ENDIF
      !===========================================================
      !SMAT=0._dp
      !!
      !IF(NGPW_L.GT.0)CALL DSYRK('U','T',NSTATE,2*NGPW_L,2.D0,CP,2*NGPW_L,0.D0,SMAT,NSTATE)
      !IF(g0_stat)CALL DGER(NSTATE,NSTATE,-1.0D0,CP,2*NGPW_L,CP,2*NGPW_L,SMAT,NSTATE)
      !CALL MPI_GlobSumR2(SMAT,NSTATE*NSTATE)
      !CALL UINV('U',SMAT,NSTATE,NSTATE)
      !IF(NGPW_L.GT.0)CALL DTRMM('R','U','N','N',2*NGPW_L,NSTATE,1.D0,SMAT,NSTATE,CP,2*NGPW_L)
      !IF(g0_stat) CALL ZCLEAN(CP,NSTATE,NGPW_L)!todo
      DEALLOCATE(SMAT)
      RETURN
      END




      SUBROUTINE UINV(UPLO,SMAT,LDA,N)
!C     == Inversion of a positive definite matrix                      =
!C     == Use for orthogonalisation                                    ==
      USE kinds
      IMPLICIT NONE
      CHARACTER UPLO*1
      INTEGER   LDA,N
      REAL(KIND=DP)    SMAT(LDA,*)
      INTEGER   INFO
      CALL DPOTRF(UPLO,N,SMAT,LDA,INFO)
      IF(INFO.NE.0) THEN
        IF(INFO.LT.0) THEN
          WRITE(*,'(A,I5,A)')' UNIV| THE I-TH (',-INFO,') ARGUMENT HAD AN ILLEGAL VALUE'
        ELSE
          WRITE(*,'(A,I5,A)')' UNIV| THE LEADING MINOR OF ORDER',INFO,' IS NOT POSITIVE DEFINITE, '
          WRITE(*,'(A)')' UNIV| AND THE FACTORIZATION COULD NOT BE COMPLETED.'
        ENDIF
        write(6,*)"UINV','ILLEGAL RESULTS DPOTRF"
        STOP
      ENDIF
      CALL DTRTRI(UPLO,'N',N,SMAT,LDA,INFO)
      IF(INFO.NE.0) THEN
      write(6,*)"UINV','ILLEGAL RESULTS DPOTRF"
      STOP
      ENDIF
      RETURN
      END
      SUBROUTINE ZCLEAN(CA,N,NGPW_L)
      USE kinds
      IMPLICIT NONE
      INTEGER N,NGPW_L
      COMPLEX*16 CA(NGPW_L,N)
      REAL(KIND=DP)  A(2,NGPW_L,N)
      INTEGER I,J,K
        DO K=1,N
        DO J=1,NGPW_L
        A(1,J,K)=DREAL(CA(J,K))
        A(2,J,K)=DIMAG(CA(J,K))
        ENDDO
        ENDDO
      IF(NGPW_l.GT.0) THEN
              DO I=1,N
           A(2,1,I)=0.D0
          ! DIMAG(CA(1,I))=0.0d0
        ENDDO
      ENDIF

        DO K=1,N
        DO J=1,NGPW_L
          CA(J,K)=CMPLX(A(1,J,K) ,A(2,J,K))
        ENDDO
        ENDDO

      RETURN
      END

    SUBROUTINE FNONLOC(CA,NOSTATE,FC)
    USE kinds
    USE system_data_types
    USE constants
    USE pseudopotential
    IMPLICIT NONE
    INTEGER NOSTATE
    COMPLEX*16, INTENT(INOUT):: CA(NGPW_L,NOSTATE)
    REAL(kind=dp) FC(NOSTATE)

    INTEGER IA,ISA,ISA0,IS,IV,I,IG,J
    COMPLEX*16 AUXC(NGPW_L,NOSTATE),C2_l(NGPW_L,NOSTATE)
    REAL(kind=dp) DDIA(NA_MAX,NOSTATE),FFI,T
    COMPLEX*16 CTM

    AUXC=(0.0D0,0.0D0)
    c2_l=(0.0_dp,0.0_dp)
    ISA0=0

      DO IS=1,sp_t
        DO IV=1,NGH(IS)
          DDIA=0.0D0
          DO I=1,NOSTATE
             CALL DCOPY(atom_p_sp(IS),NL(1,ISA0+1,IV,I),1,DDIA(1,I),1)
          ENDDO

          IF(NGPW_L.GT.0) THEN
            CALL DGEMM('N','N',2*NGPW_L,NOSTATE,atom_p_sp(IS),1.D0,EIGR_pw(1,ISA0+1),2*NGPW_L,DDIA,NA_MAX,0.0D0,AUXC,2*NGPW_L)
            DO I=1,NOSTATE
              CTM=(0.0D0,-1.0D0)**NGHTOL(IV,IS)
              FFI=FC(I)
              IF(FFI.LT.1.D-5) FFI=1.0D0
              T=-FFI*WSG(IS,IV)
              CTM=CTM*T
              DO IG=1,NGPW_L
                AUXC(IG,I)=CTM*AUXC(IG,I)
                C2_L(IG,I)=C2_L(IG,I)+TWNL(IG,IV,IS)*AUXC(IG,I)                   
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ISA0=ISA0+atom_p_sp(IS)
      ENDDO

      IF (NGPW_L.GT.0.AND.NGPW.GT.0) THEN
        DO i=1,nostate
          DO j=1,ngpw_l
            CA(j,i) = CA(j,i) + C2_l(j,i)
          ENDDO
        ENDDO
      ENDIF
    RETURN
    END

      SUBROUTINE UPF_FNONLOC(CA,NOSTATE,FC)
      USE kinds
      USE system_data_types
      USE constants
      USE pseudopotential
      IMPLICIT NONE
      INTEGER IV,L,KI,LI,JV,L2,LJ,KJ,NOSTATE,IG,ISA0,IS,I,J,ISA,ia
      REAL(kind=dp)  DDIA(NA_MAX,NOSTATE),FC(NOSTATE),FFI,T,EI,ER
      COMPLEX*16,intent(inout):: CA(NGPW_l,NOSTATE)
      COMPLEX*16 AUXC(NGPW_L,NOSTATE),CTM
      COMPLEX*16::c2_local(NGPW_L,NOSTATE)
      REAL(KIND=DP)::  t1, ti, tr
      c2_local=(0._dp,0._dp)
      ISA0=0
      DO IS=1,sp_t
          DO iv=1,upf_ngh(is)
             ddia=0.0D0
             !CALL zeroing(ddia)!,imagp*maxsys%nax*nstate)
             l=upf_nghtol(iv,is)+1
             ki=lfval(iv,is)
             li=lpval(iv,is)
             DO jv=1,UPF_ngh(is)
                l2=UPF_nghtol(jv,is)+1
                lj=lpval(jv,is)
                IF (l.EQ.l2.AND.li.EQ.lj) THEN
                   kj=lfval(jv,is)
                   DO i=1,NOSTATE
                      DO ia=1,atom_p_sp(is)
                         isa=isa0+ia
                         ddia(ia,i)=ddia(ia,i)+nl(1,isa,jv,i)*hlsg(ki,kj,l,is)
                      ENDDO
                   ENDDO
                ENDIF
             ENDDO
             IF(NGPW_L.GT.0) THEN
                ! vw here we can build only the submatrix
         CALL DGEMM('N','N',2*NGPW_L,NOSTATE,atom_p_sp(IS),1.D0,EIGR_pw(1,ISA0+1),2*NGPW_L,DDIA,NA_MAX,0.0D0,AUXC,2*NGPW_L)
                DO i=1,nostate
                 CTM=(0.0D0,-1.0D0)**upf_NGHTOL(IV,IS)
                 FFI=FC(I)
                 IF(FFI.LT.1.D-5) FFI=1.0D0
                 T=-FFI!*WSG(IS,IV)
                 CTM=CTM*T
                   ! MAKE USE OF THE SPECIAL STRUCTURE OF CTM
                   IF (ABS(REAL(ctm)).GT.0.01) THEN
                      ! CTM IS REAL
                      DO ig=1,NGpw_l
                         t1=twnl(ig,iv,is)*REAL(ctm)
                         tr=REAL(auxc(ig,i))
                         ti=AIMAG(auxc(ig,i))
                        
                         
                         C2_local(ig,i)=c2_local(ig,i)+CMPLX(t1*tr,t1*ti,kind=dp)
                       
                      
                     
                    
                      ENDDO
                   ELSE
                      ! CTM IS IMAGINARY
                      DO ig=1,NGpw_l
                         t1=twnl(ig,iv,is)*AIMAG(ctm)
                         tr=REAL(auxc(ig,i))
                         ti=AIMAG(auxc(ig,i))
                         C2_local(ig,i)=c2_local(ig,i)+CMPLX(-t1*ti,t1*tr,kind=dp)
                   
                      ENDDO
                   ENDIF
                   !IF (G0_stat)C2_local(ig,i)=0.0D0
                ENDDO        ! NSTATE
             ENDIF          ! NGWK.GT.0
          ENDDO              ! IV
     ENDDO
IF (NGPW_L.GT.0.AND.NGPW.GT.0) THEN
     ! copy back to the right C0 block
     DO i=1,nostate
        DO j=1,ngpw_l
           CA(j,i) = CA(j,i) + C2_local(j,i)
        ENDDO
     ENDDO
ENDIF
END SUBROUTINE
!     ==================================================================
      SUBROUTINE GS_DISORTHO(N_STATE,C0)
!     ==--------------------------------------------------------------==
!     ==     ORTHOGONALIZE A SET OF WAVEFUNCTIONS C0                  ==
!     ==     USING A DISTRIBUTED OVERLAP MATRIX ALGORITHM             ==
!     ==--------------------------------------------------------------==
      USE kinds
      USE mympi, ONLY : MPI_GlobSumC,MPI_GlobSumR2s,MPI_RBcast,MPI_GlobSumR
      USE system_data_types

      IMPLICIT NONE
      INTEGER N_STATE
      COMPLEX*16 C0(NGPW_L,*)
      INTEGER NORB,NORBX,MMX,MSGLEN,IP,NX,IPJ,NPJ,IPI,NPI
      INTEGER I,J,I1,J1,IJ,N1,N2,INFO,I2,MMX1
      INTEGER ISUB
      INTEGER  IERR, INDEX1, INDEX0
      REAL(KIND=DP)   OPW(N_STATE), LOCAL_OPW(N_STATE), LOCAL_NRM, NRM
      REAL(KIND=DP)   DDOT
      EXTERNAL DDOT
!     ==--------------------------------------------------------------==
     LOCAL_NRM = 2 * DDOT(2*NGpw_l, C0(1,1), 1, C0(1,1), 1)
      IF (G0_STAT) THEN
         LOCAL_NRM = LOCAL_NRM - REAL(C0(1,1))**2
      ENDIF
      MSGLEN=1*8
      !CALL MY_COMBINE(LOCAL_NRM,NRM,MSGLEN,1,ALLGRP)
      
                  CALL MPI_GlobSumR2s(LOCAL_NRM)
            NRM = DSQRT(LOCAL_NRM)
      
     
      CALL DSCAL(2*NGPW_L, 1/NRM, C0(1,1), 1 )
      DO INDEX0 = 2, N_STATE
         DO INDEX1=1,INDEX0-1 
            LOCAL_OPW(INDEX1) = 2* DDOT(2*NGPW_L, C0(1,INDEX1),1, C0(1,INDEX0),1)
            IF (G0_STAT) THEN
               LOCAL_OPW(INDEX1) = LOCAL_OPW(INDEX1) -REAL(C0(1,INDEX1))*REAL(C0(1,INDEX0))
            ENDIF
         ENDDO
         MSGLEN=(INDEX0-1)*8
         CALL MPI_GlobSumR(LOCAL_OPW,OPW,N_STATE)
         DO INDEX1=1, INDEX0-1 
            CALL DAXPY(2*NGPW_L,-OPW(INDEX1),C0(1,INDEX1), 1, C0(1,INDEX0), 1)
         ENDDO
!        NORMALIZE TO ONE
         LOCAL_NRM = 2 * DDOT(2*NGPW_L,C0(1,INDEX0),1,C0(1,INDEX0),1)
         IF (G0_STAT) THEN
            LOCAL_NRM = LOCAL_NRM -REAL(C0(1,INDEX0))*REAL(C0(1,INDEX0))
         ENDIF
         !MSGLEN=1*8
       !  CALL MY_COMBINE(LOCAL_NRM,NRM,MSGLEN,1,ALLGRP       
       CALL  MPI_GlobSumR2s(LOCAL_NRM)
            NRM = DSQRT(LOCAL_NRM)

        
         CALL DSCAL(2*NGPW_L, 1/NRM, C0(1,INDEX0), 1 )
      ENDDO

      
!     ==--------------------------------------------------------------==
      RETURN
      END

END MODULE wfn_initialize
