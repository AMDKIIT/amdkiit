     MODULE gradient

     CONTAINS

     subroutine forces
     USE system_data_types
     USE constants
     USE wfn_initialize
     USE mympi
     USE kinds
     IMPLICIT NONE
     INTEGER I,J
     REAL(kind=dp) GAM(NSTATE,NSTATE),FSC(NSTATE)

            DO I=1,NSTATE
              IF(OCCUPATION(I).EQ.0.D0) THEN
                FSC(I)=0.D0
              ELSE
                FSC(I)=1.D0
              ENDIF
            ENDDO
         
            CALL OVLAP_D(NGPW_l,NSTATE,NSTATE,GAM,C2,C_0)             
            CALL MPI_GlobSumR2(GAM,NSTATE*NSTATE)

             DO J=1,NSTATE
               DO I=1,NSTATE
                 GAM(I,J)=FSC(I)*FSC(J)*GAM(I,J)
               ENDDO
             ENDDO
            CALL DGEMM('N','N',2*NGPW_L,NSTATE,NSTATE,-1.0D0,C_0(1,1),2*NGPW_L,GAM(1,1),NSTATE,1.0D0,C2(1,1),2*NGPW_L)
            IF(g0_stat) CALL ZCLEAN(C2,NSTATE,NGPW_L)
            CALL eval_GMAX(C2)
     end subroutine

     SUBROUTINE eval_GMAX(C_2g)
     USE kinds
     USE system_data_types
     USE mympi
     USE wfn_initialize
     IMPLICIT NONE

       INTEGER NOCC,I,IIABS
       complex*16 C_2g(NGPW_L,nstate) !//1)
          NOCC=0
          DO I=1,NSTATE
            IF(OCCUPATION(I).GT.1.D-5) NOCC=NOCC+1
          ENDDO
          IIABS=IZAMAX(NOCC*NGPW_L,C_2g,1)
          GEMAX=CDABS(ZGIVE(C_2g(1,1),IIABS))
          CNORM=0.0D0
          DO I=1,NOCC
              CNORM=CNORM+DOTP(NGPW_L,C_2g(1,I),C_2g(1,I))
          ENDDO

          CALL MPI_GlobSumR2s(CNORM)
          CALL MPI_GlobMaxR2s(GEMAX) 
          CALL CHK_CNVGRAD(CONVWF,GEMAX)
       CNORM=DSQRT(CNORM/DBLE(NOCC*NGpw))
        
      END SUBROUTINE eval_GMAX

      FUNCTION ZGIVE(A,N)
      IMPLICIT NONE
      COMPLEX*16 ZGIVE
      INTEGER    N
      COMPLEX*16 A(N)
      ZGIVE=A(N)
      RETURN
      END

      FUNCTION IZAMAX(N,A,INC)
      USE kinds
      IMPLICIT NONE
      INTEGER IZAMAX
      INTEGER    N,INC
      COMPLEX*16 A(*)
      INTEGER    I,II
      REAL(kind=dp)     XMAX,Y
      IZAMAX=0
      XMAX=-1.0D0
      DO I=1,N
        II=(I-1)*INC+1
        Y=ABS(A(II))
        IF(Y.GT.XMAX) THEN
          XMAX=Y
          IZAMAX=I
        ENDIF
      ENDDO
      RETURN
      END

     !***********************************************************************
     SUBROUTINE CHK_CNVGRAD(CONV_WF, GMAX)
     !    == CHECK IF WAVEFUNCTION IS CONVERGED BASED ON GRADIENT         ==
     USE kinds
     USE system_data_types
     USE constants

     IMPLICIT NONE     
     LOGICAL,INTENT(INOUT) ::CONV_WF
     REAL(kind=dp) ::GMAX

     IF(GMAX.LT.wf_opt%g_cutoff) THEN
      CONV_WF=.TRUE.
     ELSE
      CONV_WF=.FALSE.
     ENDIF

     RETURN
     END
     !***********************************************************************
     SUBROUTINE CHK_TFOR(WFSTEP,GMAX)
     USE kinds
     USE system_data_types
     USE constants
     IMPLICIT NONE
     LOGICAL CALSTE,TNOFOR
     REAL(kind=dp)  GMAX,d_etotal
     INTEGER WFSTEP

      IF (WFSTEP.EQ.1) THEN
         E_OLD=0.0
         !CALSTE = .FALSE.
         TFOR   = .FALSE.
         TNOFOR = .FALSE.       ! Based on energy change
      ELSE if(GMAX.LT.(3*wf_opt%g_cutoff)) THEN!.AND. .NOT.(ABS(etotal-E_old).GT.0.0000)) then!TOLREL*TOLFOR*DELAST)!TOLREL*TOLFOR*TOLOG .AND. .NOT.TNOFOR) THEN
         !CALSTE = TPRCP
         TFOR   = .TRUE.
         E_OLD = etotal
      !write(6,*)tfor,e_old
      ENDIF
     ! write(6,*)"#",tfor,GMAX,3*chk_wf,(ABS(etotal-E_old).GT.0.0000)
      !`E_OLD = etotal
     RETURN
     END 

   END MODULE GRADIENT
