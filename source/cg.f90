MODULE cg!onjugate_gradient
CONTAINS
SUBROUTINE hess(VPP)
USE system_data_types
USE kinds
USE constants
USE mympi
IMPLICIT NONE
REAL(kind=dp) VPP_LOC,g0_vp
REAL(kind=dp), intent(inout) ::VPP(:)
INTEGER IG,IS,J,IV,ISA0
PRECONDITION=.TRUE.
IF(PRECONDITION) THEN
        IF(g0_stat)then 
        VPP_LOC=eigrxvps(1)
        ELSE
        VPP_LOC=0.0
        ENDIF
        CALL MPI_RBcasts(VPP_LOC,g0_ip)
        if(l_upf) then
        DO IG=1,NGPW_L
          VPP(IG)=VPP_LOC
          ISA0=0
          DO IS=1,sp_t
              DO IV=1,upf_NGH(IS)
                VPP(IG)=VPP(IG)+atom_p_sp(IS)*WSG(IS,IV)*TWNL(IG,IV,IS)*TWNL(IG,IV,IS)
                !if(icpu==2)write(1111,*)vpp(ig),ig,is,eigrxvps(1)
              ENDDO
              ISA0=ISA0+atom_p_sp(IS)
          ENDDO
        ENDDO
        else
         DO IG=1,NGPW_L
          VPP(IG)=VPP_LOC
          ISA0=0
          DO IS=1,sp_t
              DO IV=1,NGH(IS)
                VPP(IG)=VPP(IG)+atom_p_sp(IS)*WSG(IS,IV)*TWNL(IG,IV,IS)*TWNL(IG,IV,IS)
                !if(icpu==2)write(1111,*)vpp(ig),ig,is,eigrxvps(1)
              ENDDO
              ISA0=ISA0+atom_p_sp(IS)
          ENDDO
        ENDDO

        endif

        DO J=1,NGPW_l
          VPP(J)=0.5D0*TWOPIBYA2*HG(J)+VPP(J)
          IF(VPP(J).LT.HCUT) VPP(J)=HCUT
        ENDDO
        DO J=1,NGPW_l
          VPP(J)=1.0D0/VPP(J)
        ENDDO
ELSE
        DO J=1,NGPW_l
          VPP(J)=5.D0*5.D0/400.D0!!DT2BYE=DELT_ELEC*DELT_ELEC/EMASS
        ENDDO

ENDIF

END SUBROUTINE hess

SUBROUTINE PCG(istep)
USE system_data_types
USE kinds
USE constants
USE wfn_initialize
USE mympi
IMPLICIT NONE
INTEGER NOCC,I,ISTEP,NSTEP,IHIST,ILSR,IG
REAL(KIND=DP) ::VPP(NGPW_l),GGNORM,ALAM,DE
REAL(KIND=DP) ::GHIST(10),FHIST(10),DETOT,GAMMA,FORGET
SAVE     GHIST,FHIST
COMPLEX*16 C_2(NGPW_l,NSTATE)!,HNM1(NGPW_l,NSTATE)
!COMPLEX(KIND=dp), DIMENSION(:,:),       POINTER ::
LOGICAL PCGMIN,PREC !todo
REAL(KIND=DP)  DT2BYE
SAVE    ILSR
DT2BYE=5.D0*5.D0/400.D0
!      DELT_ELEC=5.D0
!      EMASS=400.D0

PREC=.TRUE.
PCGMIN=.TRUE.
ILSR=0




IF(.NOT.ASSOCIATED(HNM1))ALLOCATE(HNM1(NGPW_l,NSTATE))
CALL HESS(VPP)
!HNM1(NGPW_l,NSTATE)=(0._dp,0._dp)

      NOCC=0
      DO I=1,NSTATE
        IF(OCCUPATION(I).GT.1.D-5) THEN
          NOCC=NOCC+1
        ENDIF
      ENDDO

      IF(PRECONDITION)CALL PRECON(VPP,C_2)
      GGNORM=0.0D0
      DO I=1,NOCC
        GGNORM=GGNORM+DOTP(NGPW_l,C_2(1,I),C_2(1,I))
      ENDDO

      CALL MPI_GlobSumR2s(GGNORM)
       
      if(istep.eq.1)then
          FHIST(1)=5.D0*5.D0/400.D0
          GHIST(1)=GGNORM
          CALL DCOPY(2*NGPW_l*NSTATE,C_2(1,1),1,HNM1(1,1),1)
          CALL DAXPY(2*NGPW_l*NSTATE,FHIST(1),C_2(1,1),1,C_0(1,1),1)
      else
          IF(ISTEP.LE.10) THEN
!       Almost steepest descent
          FORGET=0.60D0
          ELSEIF(ISTEP.LE.20) THEN
          FORGET=0.91D0
          ELSE
          FORGET=0.95D0
          ENDIF

          IF(ISTEP.GT.10) THEN
          IHIST=10
          DO I=1,10-1
            FHIST(I)=FHIST(I+1)
            GHIST(I)=GHIST(I+1)
          ENDDO
          ELSE
          IHIST=ISTEP
          ENDIF
        GHIST(IHIST)=GGNORM
        FHIST(IHIST)=FHIST(IHIST-1)
        GAMMA=GHIST(IHIST)/GHIST(IHIST-1)

!       Multiply by a forget factor
        GAMMA=FORGET*GAMMA
        CALL DAXPY(2*NGPW_l*NSTATE,GAMMA,HNM1(1,1),1,C_2(1,1),1)
        CALL DCOPY(2*NGPW_l*NSTATE,C_2(1,1),1,HNM1(1,1),1)
        CALL DAXPY(2*NGPW_l*NSTATE,FHIST(IHIST),HNM1(1,1),1,C_0(1,1),1)
      !print*,"GGNORM",GGNORM
      endif
      !DEALlOCATE(HNM1)
END SUBROUTINE

       SUBROUTINE LINESR(P,FAC,X,Cc2,ALAM,DE,ILSR)
       USE system_data_types
       USE kinds
       USE constants
       USE wfn_initialize
       USE mympi
       IMPLICIT NONE
       COMPLEX*16 P(NGPW_l,NSTATE),X(:,:),Cc2(:,:)
       REAL(KIND=DP)  FAC,DD,ALAM,ALAM0
       REAL(KIND=DP)  EEE,E0,E1,E2,DX,X0,A,DE
       integer I,ILSR
       CALL DSCAL(2*NGPW_l*NSTATE,FAC,P(1,1),1)
       E0=ETOTAL
       !if
        ALAM=1.0D0
        DD=0.0D0
        DO I=1,NSTATE
          DD=DD+DOTP(NGPW_l,P(1,I),P(1,I))
        ENDDO
        CALL MPI_GlobSumR2s(DD)
       !CALL GLOSUM(1,DD)
!       First point
        CALL DCOPY(2*NGPW_l*NSTATE,X(1,1),1,Cc2(1,1),1)
        CALL DAXPY(2*NGPW_l*NSTATE,ALAM,P(1,1),1,X(1,1),1)
        ALAM0=ALAM
        CALL XETOT(X)!(NSTATE,X,SC0,TAU0,TSCR,RHOE,PSI,SCR,LSCR)
        IF(ETOTAL.GT.E0) THEN
          E2=ETOTAL
          ALAM=0.5D0
        ELSE
          E1=ETOTAL
          ALAM=2.0D0
        ENDIF
!C       Second point
        CALL DCOPY(2*NGPW_l*NSTATE,Cc2(1,1),1,X(1,1),1)
        CALL DAXPY(2*NGPW_l*NSTATE,ALAM,P(1,1),1,X(1,1),1)
        ALAM0=ALAM
!        Calculate the new point
        CALL XETOT(X)!,SC0,TAU0,TSCR,RHOE,PSI,SCR,LSCR)
        IF(ALAM.GT.1.D0) THEN
          E2=ETOTAL
          DX=1.0D0
        ELSE
          E1=ETOTAL
          DX=0.5D0
        ENDIF
!C       Extrapolation
        IF((2*E0-4*E1+2*E2).EQ.0.D0) THEN
!C         Very rare (T.D.)
          X0=DX
        ELSE
          X0=DX*(4*E1-E2-3*E0)/(2*E0-4*E1+2*E2)
        ENDIF
        IF((2*X0+DX).EQ.0.D0) THEN
          A=0.D0
        ELSE
          A=(E1-E0)/DX/(2*X0+DX)
        ENDIF

        EEE=E0-A*X0*X0
        ALAM=-X0
        IF(ALAM.GT.0.D0) THEN
!C         Limit the move
          ILSR=0
          IF(ALAM.GT.3.D0) THEN
!C           New minimisation
            ALAM=3.0D0
          ENDIF
        ELSE
!C         ALAM is negative, take the lowest energy between E0, E1 and E2
          IF(E2.LT.E0) THEN
            ALAM=2*DX
            ILSR=0
          ELSEIF(E1.LT.E0) THEN
            ALAM=DX
            ILSR=0
          ELSE
!C           E0 < E1 and E2
            ALAM=0.05D0
!C           Reset the direction
            ILSR=-2
          ENDIF
        ENDIF
        DE=EEE+A*(X0+ALAM)**2
        A=A/DD
!C       New point
        CALL DCOPY(2*NGPW_l*NSTATE,CC2(1,1),1,X(1,1),1)
        CALL DAXPY(2*NGPW_l*NSTATE,ALAM,P(1,1),1,X(1,1),1)
      !ENDIF

      CALL DSCAL(2*NGPW_l*NSTATE,1.D0/FAC,P(1,1),1)
      IF(ALAM.LT.0.1D0) THEN
        FAC=FAC*0.25D0
      ELSEIF(ALAM.GT.1.0D0) THEN
        FAC=FAC*1.50D0
      ELSEIF(ALAM.LT.0.2D0) THEN
        FAC=FAC*0.75D0
      ENDIF
        
       END SUBROUTINE
       
       SUBROUTINE XETOT(X)
       USE system_data_types
       USE kinds
       USE constants
       USE wfn_initialize
       USE density
       USE potential
       use xc
       USE pseudopotential
       USE total_energy

       IMPLICIT NONE
       COMPLEX*16 X(:,:)!,Cc2(:,:)
       !REAL*8  FAC,DD,ALAM,ALAM0
       integer IS,NOCC,NP
       NOCC=0
        DO IS=1,NSTATE
        !  IF(OCCUPATION(IS).GT.1.D-5) NOCC=NOCC+1
        ENDDO
        IF(NOCC.LT.NSTATE) THEN
       !   NP=NSTATE-NOCC
          !CALL GS_ORTHO(C0,NOCC,C0(1,NOCC+1),NP,SCR(NSTATE+1),SCR)
        ENDIF
        CALL ortho_gs(X)!,ngpw_l,nstate)
        CALL eval_density
        CALL cal_vpot_ex
        CALL pp_energy
        CALL e_total

       END SUBROUTINE



      SUBROUTINE PRECON(VPP,C_2)
       USE system_data_types
       USE kinds
       USE constants

      IMPLICIT NONE
           
      
      REAL(kind=dp) ::VPP(:)
      COMPLEX*16 C_2(:,:)
      REAL(KIND=DP)     FI(NSTATE)
      INTEGER    I,J

      DO I=1,NSTATE
        IF(OCCUPATION(I).GT.0.01D0) THEN
          FI(I)=1.0D0/OCCUPATION(I)
        ELSE
          FI(I)=1.0D0
        ENDIF
      ENDDO

        DO J=1,NGPW_l
          DO I=1,NSTATE

            C_2(J,I)=FI(I)*VPP(J)*C2(J,I)


          ENDDO
        ENDDO
      RETURN
      END

END MODULE
