     MODULE nuclear_grad

     CONTAINS
     SUBROUTINE POTFOR
     
     USE system_data_types
     USE constants
     USE gradient
     USE kinds
     IMPLICIT NONE
 
     COMPLEX*16 RP,RHET,RHOG,RHETS,RHOGS,GX,GY,GZ,VCGS,TXX,TYY,TZZ,EI123
     REAL(KIND=DP)    OMTP
     INTEGER    IG1,IG,ISA,IS,IA
      OMTP=2.D0*cell_volume*twopibya
      IG1=1
      IF(g0_stat) IG1=2
        DO IG=IG1,ngrho_l
          rp=eigrxrhos(ig)
          rhet=rho_g(ig,1)
          rhog=rhet+rp
          rhets=dconjg(rhet)
          rhogs=dconjg(rhog)
          !write(49,*)"rp",rp
          !write(49,*)"rh",rhet
          gx=dcmplx(0.d0,gvec(1,ig))
          gy=dcmplx(0.d0,gvec(2,ig))
          gz=dcmplx(0.d0,gvec(3,ig))
          vcgs=(DCMPLX(fourpi/(twopibya2*hg(ig)),0.D0))*rhogs
          isa=0
          do is=1,sp_t
            txx=(rhos(is,ig)*vcgs+vps(is,ig)*rhets)*gx
            tyy=(rhos(is,ig)*vcgs+vps(is,ig)*rhets)*gy
            tzz=(rhos(is,ig)*vcgs+vps(is,ig)*rhets)*gz
            !write(47,*)gz!(DCMPLX(fourpi/(twopibya2*hg(ig)),0.D0))
            do ia=1,atom_p_sp(is)
              isa=isa+1
              !write(46,*)"tzz=",tzz,"eigr=",eigr(ig,isa)!)!*tzz)*omtp
              !write(46,*) "tzz*eigr=",eigr(ig,isa)*tzz
              force(1,ia,is)=force(1,ia,is)+dreal(eigr(ig,isa)*txx)*omtp
              force(2,ia,is)=force(2,ia,is)+dreal(eigr(ig,isa)*tyy)*omtp
              force(3,ia,is)=force(3,ia,is)+real(eigr(ig,isa)*tzz)*omtp
              !write(46,*)"*",force(3,ia,is)
            enddo
          enddo
        enddo


      RETURN
      END
       SUBROUTINE NLFOR
         USE system_data_types
         USE constants
         USE gradient
         USE pseudopotential
         USE kinds
         IMPLICIT NONE

      REAL(KIND=DP)        WEIGHT,TEMP!,TT,TDBL
      INTEGER ISUB,IK,ISA0,IS,IV,JV,I,II,ISPIN,ISA,L,KI,LI,L2,LJ,KJ,IA
      real(KIND=DP)  wk1_1,wk1_2,wk1_3,wk2_1,wk2_2,wk2_3
      CALL D_FNL
          ISA0=0
          DO IS=1,sp_t
              DO IV=1,NGH(IS)
                TEMP=2.D0*WSG(IS,IV)
               ! WRITE(6,*)"**",TEMP,IS,IV
                DO I=NST12(icpu,1),NST12(icpu,2)! do this part loadpa.f NST12(MEPOS,1),NST12(MEPOS,2)
                  WEIGHT=OCCUPATION(I) !WK(IK)*F(I)
                  IF(ABS(WEIGHT).GT.1.D-12) THEN
                    II=I-NST12(icpu,1)+1
                        DO IA=1,atom_p_sp(IS)
                          ISA=ISA0+IA
                          !write(6,*)"fnl",nl(1,1,1,1)
                          wk1_1=dfnl(1,isa,iv,1,ii)*nl(1,isa,iv,   i)
                          wk1_2=dfnl(1,isa,iv,2,ii)*nl(1,isa,iv,   i)
                          wk1_3=dfnl(1,isa,iv,3,ii)*nl(1,isa,iv,   i)
                          force(1,ia,is)=force(1,ia,is)-temp*weight*wk1_1
                          force(2,ia,is)=force(2,ia,is)-temp*weight*wk1_2
                          force(3,ia,is)=force(3,ia,is)-temp*weight*wk1_3
                        ENDDO
                  ENDIF    
                ENDDO
              ENDDO
            
            ISA0 = ISA0 + atom_p_sp(IS)
          ENDDO
       DEALLOCATE(DFNL)
       END SUBROUTINE
       SUBROUTINE UPF_D_FNL
       USE system_data_types
       USE constants
       USE gradient
       USE pseudopotential
       USE mympi
       USE kinds
       IMPLICIT NONE

       INTEGER   k,isa0,is,iv,ia,ig,isa,ikind,nat,ik,i,nstat,II
       REAL(kind=dp)    CIR,CII,ARG,ER,EI,TFAC
       COMPLEX*16 ci,eiscr(ngpw_l,NA_MAX)!,C0(NGPW,NSTAT)
       REAL(kind=dp)     DAI(1,NA_MAX,NSTATE)
      IF(.NOT.ASSOCIATED(DFNL)) ALLOCATE(DFNL(1,atom_t,UPF_NGHMAX,3,NORBPE))!to do
      TFAC=2.D0*twopibya
      DAI=0.0
       DO K=1,3
        ISA0=0
        DO IS=1,sp_t
          DO IV=1,upf_NGH(IS)
            CI=(0.0D0,-1.0D0)**(upf_NGHTOL(IV,IS)+1)
            CIR=DREAL(CI)
            CII=DIMAG(CI)
              DO IA=1,atom_p_sp(IS)
                ISA=ISA0+IA
                IF(DABS(CIR).GT.0.5D0) THEN
                 !CI is real
                  DO IG=1,NGPW_l
                    ARG=GVEC(K,IG)*TWNL(IG,IV,IS)*CIR
                    ER=DREAL(EIGR(IG,ISA))
                    EI=DIMAG(EIGR(IG,ISA))
                    EISCR(IG,IA) = DCMPLX(ARG*ER,ARG*EI)
                  ENDDO
                ELSE
                 !CI is imaginary
                  DO IG=1,NGPW_l
                    ARG=GVEC(K,IG)*TWNL(IG,IV,IS)*CII
                    ER=DREAL(EIGR(IG,ISA))
                    EI=DIMAG(EIGR(IG,ISA))
                    EISCR(IG,IA) = DCMPLX(-ARG*EI,ARG*ER)
                  ENDDO
                ENDIF
                IF(g0_stat) EISCR(1,IA)=0.5D0*EISCR(1,IA)
              ENDDO
              IF (atom_p_sp(IS).GT.1) THEN
            CALL DGEMM('T','N',atom_p_sp(IS),NSTATE,2*NGPW_l,TFAC,EISCR(1,1),&
                           2*NGPW_l,C_0(1,1),2*NGPW_l,0.0D0, DAI(1,1,1),NA_MAX)
              ELSE
                CALL DGEMV('T',2*NGPW_l,NSTATE,TFAC,C_0(1,1),2*NGPW_l,EISCR(1,1),1,0.0D0,DAI(1,1,1),NA_MAX)
              ENDIF
            CALL MPI_GlobSumR2(DAI,NA_MAX*NSTATE)
            DO I=NST12(ICPU,1),NST12(ICPU,2)
              II=I-NST12(ICPU,1)+1
              CALL DCOPY(1*atom_p_sp(IS),DAI(1,1,I),1,DFNL(1,ISA0+1,IV,K,II),1)
            ENDDO
          ENDDO
          ISA0=ISA0+atom_p_sp(IS)
        ENDDO
       ENDDO
    END


!-----------------------------------------------------------------------------------------------------------
SUBROUTINE UPF_NLFOR
USE system_data_types
USE constants
USE gradient
USE pseudopotential
USE kinds

IMPLICIT NONE

REAL(KIND=DP)        WEIGHT,TEMP,TT
INTEGER ISUB,IK,ISA0,IS,IV,JV,I,II,ISPIN,ISA,L,KI,LI,L2,LJ,KJ,IA
!REAL(KIND=DP)  
CALL UPF_D_FNL
       isa0=0
       DO is=1,sp_t
         ! Stefan Goedecker pp
           DO iv=1,upf_ngh(is)
                l=upf_nghtol(iv,is)+1
                ki=lfval(iv,is)
                li=lpval(iv,is)
                DO jv=1,upf_ngh(is)
                   l2=upf_nghtol(jv,is)+1
                   lj=lpval(jv,is)
                   IF (l2.EQ.l.AND.li.EQ.lj) THEN
                      kj=lfval(jv,is)
                      !DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                      DO I=NST12(ICPU,1),NST12(ICPU,2)
                         WEIGHT=OCCUPATION(I)
                         IF(ABS(WEIGHT).GT.1.D-12) THEN
                         II=I-NST12(icpu,1)+1
                            tt=2.0*weight*hlsg(ki,kj,l,is)
                                  DO ia=1,atom_p_sp(is)
                                     isa=isa0+ia
                                     force(1,ia,is)=force(1,ia,is)-tt*nl(1,isa,jv,i)*dfnl(1,isa,iv,1,ii)
                                     force(2,ia,is)=force(2,ia,is)-tt*nl(1,isa,jv,i)*dfnl(1,isa,iv,2,ii)
                                     force(3,ia,is)=force(3,ia,is)-tt*nl(1,isa,jv,i)*dfnl(1,isa,iv,3,ii)
                                   ENDDO
                                        
                         ENDIF
                       
                      ENDDO
                   ENDIF
                ENDDO
           ENDDO
       ENDDO

END SUBROUTINE


!-----------------------------------------------------------------------------------------------------------

















    subroutine convert3dto1d(a,m,n,o,b,mno)
    USE kinds
    USE system_data_types

    implicit none

    integer, intent(in) :: m,n,o
    REAL(KIND=DP), intent(in) :: a(m,n,o) 
    integer, intent(in) :: mno
    REAL(KIND=DP), intent(inout) :: b(mno)  
    INTEGER I,J,K,count
   
   count=1
   do k = 1, o
    do j = 1, atom_p_sp(k)
      do i = 1, m
        b(count)=a(I,J,K)
        count=count+1
      end do
    end do
  enddo
    end subroutine

    subroutine convert1dto3d(b,mno,a,m,n,o)
    USE kinds
    USE system_data_types

    implicit none

    integer, intent(in) :: m,n,o
    REAL(KIND=DP), intent(inout) :: a(m,n,o)
    integer, intent(in) :: mno
    REAL(KIND=DP), intent(in) :: b(mno)
    INTEGER I,J,K,count

   count=1
   do k = 1, o
    do j = 1, atom_p_sp(k)
      do i = 1, m
        a(I,J,K)=b(count)
        count=count+1
      end do
    end do
  enddo
    end subroutine




    subroutine PRINT_COORDINATE(a,forc,n,o,totatom,totenergy,gopt_step)
    USE kinds
    USE system_data_types

    implicit none

    integer, intent(in) :: n,o,totatom,gopt_step
    REAL(KIND=DP), intent(in) :: a(3,n,o),forc(3,n,o),totenergy

    INTEGER I,J,K
    
    if(ionode)then
    !file_name = 'GEOOPT.xyz'!// trim(adjustl(file_id)) // '.xyz'
!!    open(12,file = 'atom_coord.xyz')!trim(file_name))
    OPEN(12,FILE="atom_coord.xyz",STATUS="unknown",position="APPEND")

    WRITE(12,*)totatom
    write(12,*)"ENERGY=",totenergy, gopt_step!!"GRADIENT=",gnorm
         do k = 1, o
          do j = 1, atom_p_sp(k)
           write(12,*)SYMBOL(z(k)),((a(i,j,k)*0.529177),i=1,3)!,(force(IS,K,J),is=1,3)!0.529177

           end do
         end do
!    write(12,*)"   " ## RK edit
    CLOSE(12)
    endif
    
    if(ionode)then
    !open(13,file = 'atom_force.dat')!trim(file_name))
    OPEN(13,FILE="atom_force.dat",STATUS="unknown",position="APPEND")
    WRITE(13,*)totatom
    write(13,*)"ENERGY=",totenergy, gopt_step!!"GRADIENT=",gnorm
         do k = 1, o
          do j = 1, atom_p_sp(k)
            write(13,*)SYMBOL(z(k)),((forc(i,j,k)*1),i=1,3)!,(force(IS,K,J),is=1,3)!0.529177

           end do
         end do
!    write(13,*)"   " 
    CLOSE(13)
    endif

    end subroutine

       SUBROUTINE D_FNL
       USE system_data_types
       USE constants
       USE gradient
       USE pseudopotential
       USE mympi
       USE kinds
       IMPLICIT NONE

       INTEGER   k,isa0,is,iv,ia,ig,isa,ikind,nat,ik,i,nstat,II
       REAL(kind=dp)    CIR,CII,ARG,ER,EI,TFAC
       COMPLEX*16 ci,eiscr(ngpw_l,NA_MAX)!,C0(NGPW,NSTAT)
       REAL(kind=dp)     DAI(1,NA_MAX,NSTATE)
      IF(.NOT.ASSOCIATED(DFNL)) ALLOCATE(DFNL(1,atom_t,NGHMAX,3,NORBPE))!to do
      TFAC=2.D0*twopibya
      DAI=0.0
       DO K=1,3
        ISA0=0
        DO IS=1,sp_t
          DO IV=1,NGH(IS)
            CI=(0.0D0,-1.0D0)**(NGHTOL(IV,IS)+1)
            CIR=DREAL(CI)
            CII=DIMAG(CI)
              DO IA=1,atom_p_sp(IS)
                ISA=ISA0+IA
                IF(DABS(CIR).GT.0.5D0) THEN
                 !CI is real
                  DO IG=1,NGPW_l
                    ARG=GVEC(K,IG)*TWNL(IG,IV,IS)*CIR
                    ER=DREAL(EIGR(IG,ISA))
                    EI=DIMAG(EIGR(IG,ISA))
                    EISCR(IG,IA) = DCMPLX(ARG*ER,ARG*EI)
                  ENDDO
                ELSE
!CI is imaginary
                  DO IG=1,NGPW_l
                    ARG=GVEC(K,IG)*TWNL(IG,IV,IS)*CII
                    ER=DREAL(EIGR(IG,ISA))
                    EI=DIMAG(EIGR(IG,ISA))
                    EISCR(IG,IA) = DCMPLX(-ARG*EI,ARG*ER)
                  ENDDO
                ENDIF
                IF(g0_stat) EISCR(1,IA)=0.5D0*EISCR(1,IA)
              ENDDO
              IF (atom_p_sp(IS).GT.1) THEN
                CALL DGEMM('T','N',atom_p_sp(IS),NSTATE,2*NGPW_l,TFAC,EISCR(1,1),2*NGPW_l,C_0(1,1),2*NGPW_l,0.0D0,DAI(1,1,1),NA_MAX)
              ELSE
                CALL DGEMV('T',2*NGPW_l,NSTATE,TFAC,C_0(1,1),2*NGPW_l,EISCR(1,1),1,0.0D0,DAI(1,1,1),NA_MAX)
              ENDIF
            CALL MPI_GlobSumR2(DAI,NA_MAX*NSTATE)
            DO I=NST12(ICPU,1),NST12(ICPU,2)
              II=I-NST12(ICPU,1)+1
              CALL DCOPY(1*atom_p_sp(IS),DAI(1,1,I),1,DFNL(1,ISA0+1,IV,K,II),1)
            ENDDO
          ENDDO
          ISA0=ISA0+atom_p_sp(IS)
        ENDDO
       ENDDO
    END

END MODULE
