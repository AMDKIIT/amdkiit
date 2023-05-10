      SUBROUTINE PBC(X1,Y1,Z1,X2,Y2,Z2,M,A,IBRAV)
      USE system_data_types
      USE kinds
!     ==--------------------------------------------------------------==
!     == CALCULATES THE PERIODIC BOUNDARY CONDITIONS FOR SOME         ==
!     == STRUCTURES:                                                  ==
!     ==             SIMPLE CUBIC          (IBRAV=1)                  ==
!     ==             FACE CENTERED CUBIC   (IBRAV=2)                  ==
!     ==             ORTHOROMBIC           (IBRAV=8)                  ==
!     ==             TRIGONAL              (IBRAV=12)                 ==
!     ==                                                              ==
!     == OTHER THAN THESE, WE USE A GENERAL ROUTINE                   ==
!     ==                                                              ==
!     == FOR FCC THE ROUTINE IS TAKEN FROM                            ==
!     == CCP5 INF. QUAT. #10,SEPT(1983)                               ==
!     ==--------------------------------------------------------------==
      IMPLICIT NONE

!     Arguments
      REAL(kind=dp)  X1,Y1,Z1,X2,Y2,Z2,A(4)
      INTEGER M,IBRAV
!     Variables
      REAL(kind=dp)  B(4),C(4),DX,R(3),S(3)
      REAL(kind=dp)  TA,AA,S2,X2T,Y2T,XT,X0,Y0,X3,Y3,XIC
      INTEGER I
      REAL(kind=dp)   V1CRYS(4),V2CRYS(4),V3CRYS(4)

      IF(SYMMETRY.ne.12)THEN
       DO I=1,4
          V1CRYS(I)=0.D0
          V2CRYS(I)=0.D0
          V3CRYS(I)=0.D0
        ENDDO 
      ENDIF
!     ==--------------------------------------------------------------==
      TA=2.D0*A(1)*M
      AA=A(1)*M
      !GO TO (5) IBRAV
      GO TO (5,10,99,25,99,99,99,15,99,99,99,20,99,25) IBRAV
!     ==--------------------------------------------------------------==
!     == Simple Cubic                                                 ==
!     ==--------------------------------------------------------------==
    5 CONTINUE
      X2=DMOD(X1,TA)
      Y2=DMOD(Y1,TA)
      Z2=DMOD(Z1,TA)
      X2=X2-TA*DINT(X2/AA)
      Y2=Y2-TA*DINT(Y2/AA)
      Z2=Z2-TA*DINT(Z2/AA)
      RETURN
!     ==--------------------------------------------------------------==
!     == Face Centered Cubic                                          ==
!     ==--------------------------------------------------------------==
   10 CONTINUE
1     S2=DSQRT(2.0D0)
      X2T=( X1+Y1)/S2
      Y2T=(-X1+Y1)/S2
      X2=X2T
      Y2=Y2T
      X2=DMOD(X2,TA)
      Y2=DMOD(Y2,TA)
      Z2=DMOD(Z1,TA*S2)
      X2=X2-TA*DINT(X2/AA)
      Y2=Y2-TA*DINT(Y2/AA)
      Z2=Z2-TA*S2*DINT(Z2/(S2*AA))
      IF(DABS(X2)+DABS(Y2)+S2*DABS(Z2).LT.TA) GOTO 11
      X2=X2-DSIGN(AA,X2)
      Y2=Y2-DSIGN(AA,Y2)
      Z2=Z2-DSIGN(S2*AA,Z2)
   11 CONTINUE
      XT=(X2-Y2)/S2
      Y2=(X2+Y2)/S2
      X2=XT
      RETURN
!     ==--------------------------------------------------------------==
!     == Orthorhombic                                                 ==
!     ==--------------------------------------------------------------==
   15 CONTINUE
      DO I=1,3
        B(I)=2*A(I)*M
        C(I)=A(I)*M
      ENDDO
      X2=DMOD(X1,B(1))
      Y2=DMOD(Y1,B(2))
      Z2=DMOD(Z1,B(3))
      X2=X2-B(1)*DINT(X2/C(1))
      Y2=Y2-B(2)*DINT(Y2/C(2))
      Z2=Z2-B(3)*DINT(Z2/C(3))
      RETURN
!     ==--------------------------------------------------------------==
!     == Hexagonal and special case of triclinic: monocl., gamma<>90  ==
!     ==--------------------------------------------------------------==
 25   CONTINUE
      IF(A(2).EQ.0.0D0) GOTO 99
      DO I=1,3
        C(I)=A(I)*M
        B(I)=2*C(I)
      ENDDO
      Y0=Y1
      Y2=DMOD(Y1,B(2))
      Z2=DMOD(Z1,B(3))
      Y2=Y2-B(2)*DINT(Y2/C(2))
      Z2=Z2-B(3)*DINT(Z2/C(3))
      DX=(Y2-Y0)*A(4)           ! A(4)=SQRT(3)/3
      X2=DMOD(X1-DX,B(1))
      X2=X2-B(1)*DINT(X2/C(1))
      RETURN
!     ==--------------------------------------------------------------==
!     == Trigonal                                                     ==
!     ==--------------------------------------------------------------==
   20 CONTINUE
      DO 22 I=1,4
      B(I)=2*A(I)*M
      C(I)=A(I)*M
   22 CONTINUE
      Z2=DMOD(Z1,B(3))
      Z2=Z2-B(3)*DINT(Z2/C(3))
      CALL ROT(X1,Y1,X0,Y0,V1CRYS,1)
      X0=DMOD(X0,B(1))
      X0=X0-B(1)*DINT(X0/C(1))
      CALL ROT(X0,Y0,X3,Y3,V1CRYS,-1)
      CALL ROT(X3,Y3,X0,Y0,V2CRYS,1)
      X0=DMOD(X0,B(2))
      X0=X0-B(2)*DINT(X0/C(2))
      CALL ROT(X0,Y0,X3,Y3,V2CRYS,-1)
      CALL ROT(X3,Y3,X0,Y0,V3CRYS,1)
      X0=DMOD(X0,B(4))
      X0=X0-B(4)*DINT(X0/C(4))
      CALL ROT(X0,Y0,X3,Y3,V3CRYS,-1)
      CALL ROT(X3,Y3,X0,Y0,V1CRYS,1)
      X0=DMOD(X0,B(1))
      X0=X0-B(1)*DINT(X0/C(1))
      CALL ROT(X0,Y0,X2,Y2,V1CRYS,-1)
      RETURN
!     ==--------------------------------------------------------------==
!     == The rest : use a general procedure                           ==
!     ==--------------------------------------------------------------==
   99 CONTINUE
      R(1)=X1
      R(2)=Y1
      R(3)=Z1
      CALL DGEMV('T',3,3,1.0D0,b_cell,3,R,1,0.0D0,S,1)
                               !HTM1
      XIC=DBLE(M)
      S(1)=S(1)-NINT(S(1)/XIC)*XIC
      S(2)=S(2)-NINT(S(2)/XIC)*XIC
      S(3)=S(3)-NINT(S(3)/XIC)*XIC
      CALL DGEMV('T',3,3,1.0D0,a_cell,3,S,1,0.0D0,R,1)
                              !HT
      X2=R(1)
      Y2=R(2)
      Z2=R(3)
!     ==--------------------------------------------------------------==
      RETURN
      END
!     ==================================================================
      SUBROUTINE ROT(X1,Y1,X0,Y0,VV,ISIGN)
!     ==--------------------------------------------------------------==
!     == THIS ROUTINE TRANSFORMS THE POINT (X1,Y1) INTO THE POINT     ==
!     == (X0,Y0) CORRESPONDING TO A ROTATION OF AN ANGLE WHICH COS    ==
!     == AND SIN ARE GIVEN BY VV(3) AND VV(4) RESPECTIVELY. THE       ==
!     == SENSE OF THE ROTATION IS GIVEN BY ISIGN (+ FOR CLOCKWISE,    ==
!     == AND - FOR ANTICLOCKWISE)                                     ==
!     ==--------------------------------------------------------------==
      IMPLICIT NONE
!     Arguments
      INTEGER ISIGN
      REAL*8  X1,Y1,X0,Y0,VV(4)
      X0=X1*VV(3)+Y1*VV(4)*DBLE(ISIGN)
      Y0=-X1*VV(4)*DBLE(ISIGN)+Y1*VV(3)
!     ==--------------------------------------------------------------==
      RETURN
      END
!     ==================================================================
