MODULE random
USE kinds

IMPLICIT NONE

PUBLIC :: uniform_rnd

CONTAINS

SUBROUTINE uniform_rnd(n,a)
   INTEGER :: n
   REAL(KIND=DP), DIMENSION(n) :: a
   REAL(KIND=DP) :: rnd
   INTEGER :: i
   integer :: seed = 1

   CALL random_seed(seed)
  
   DO i=1,n
      CALL RANDOM_NUMBER(rnd)
      a(i)=rnd
   ENDDO
END SUBROUTINE uniform_rnd

END MODULE random
