MODULE INTERPOLGRID  
CONTAINS 
      SUBROUTINE INTRP_GRID(RP,MP,RW,MW,F,MMAXX)
!     ==--------------------------------------------------------------==
!     == Interpolate F(1:MP) from RP(1:MP) grid to RW(1:MW)           ==
!     == Output: F(1:MW)                                              ==
!     ==--------------------------------------------------------------==
      USE kinds
      IMPLICIT NONE
!     Arguments
      INTEGER MP,MW,MMAXX
      REAL(kind=dp)  RP(MP),RW(MW),F(MMAXX),CC
      REAL(KIND=DP)  FC(MMAXX)
      INTEGER I,IERR
!     ==--------------------------------------------------------------==
      FC=0.0D0 
      CALL dCOPY(mp,F(1),1,FC(1),1)
     
      DO I=1,MW
        call spline_inter(MP,RP,FC,RW(I),CC,1)!
        f(i)=cc
      ENDDO

      RETURN
      END
      
END
