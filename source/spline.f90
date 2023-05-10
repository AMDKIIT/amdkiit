MODULE spline
CONTAINS

SUBROUTINE spline_params(n,xp,yp,h,v,u)
USE kinds
  IMPLICIT NONE
  INTEGER, INTENT(IN):: n
  REAL(kind=dp), POINTER :: xp(:), yp(:)
  REAL(kind=dp), POINTER :: h(:),v(:),u(:)
  INTEGER :: i

  REAL(kind=dp), ALLOCATABLE :: b(:)

  ALLOCATE(b(n))
  DO i=1,n-1
     h(i)=xp(i+1)-xp(i)
     b(i)=(yp(i+1)-yp(i))/h(i)
  END DO
  v(1)=0.d0
  u(1)=0.d0
  DO i=2,n-1
    v(i)=2.d0*(h(i-1)+h(i))
    u(i)=6.d0*(b(i)-b(i-1))
  END DO
  DEALLOCATE(b)
END SUBROUTINE spline_params


SUBROUTINE spline_thomas(n,h,v,u,sdd)
USE kinds
!Ref: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  REAL(kind=dp), POINTER     :: h(:),v(:),u(:)
  REAL(kind=dp), POINTER     :: sdd(:) 

  REAL(kind=dp), ALLOCATABLE :: d(:),a(:),b(:),c(:),x(:)
  REAL(kind=dp) :: s
  INTEGER :: i, m

!Sisdde of the matrix is m
  m=n-2

  ALLOCATE(d(m))
  ALLOCATE(a(m))
  ALLOCATE(b(m))
  ALLOCATE(c(m))
  ALLOCATE(x(m))

  !Setup the dmatrix
  DO i=1,m
   d(i)=v(i+1)
   a(i)=h(i+1)
   c(i)=h(i+1)
   b(i)=u(i+1)
  END DO

  !Elimination
  DO i=2,m
    s=a(i-1)/d(i-1)
    d(i)=d(i)-s*c(i-1)
    b(i)=b(i)-s*b(i-1) 
  END DO

  !Backward Substitution
  x(m)=b(m)/d(m)
  DO i=m-1,1,-1
    x(i)=(b(i)-c(i)*x(i+1))/d(i)
  END DO
  sdd(1)=0.d0
  sdd(n)=0.d0
  DO i=2,n-1
     sdd(i)=x(i-1)
  END DO
  DEALLOCATE(d,a,b,c,x)
END SUBROUTINE spline_thomas

SUBROUTINE spline_interpolate(n,xp,yp,h,v,u,sdd,ig,x,y)
  USE kinds
  IMPLICIT NONE
  INTEGER :: n, ig
  REAL(kind=dp), POINTER :: xp(:),yp(:),h(:),v(:),u(:),sdd(:)
  REAL(kind=dp) :: x, y

  !if(ig==n)ig=ig-1

  y=sdd(ig+1)/6.d0/h(ig)*(x-xp(ig))**3 + &
    sdd(ig)/  6.d0/h(ig)*(xp(ig+1)-x)**3 + &
    (yp(ig+1)/h(ig) - sdd(ig+1)/6.d0*h(ig) )*(x-xp(ig)) + &
    (yp(ig)/h(ig) - h(ig)/6.d0*sdd(ig))*(xp(ig+1)-x)
END SUBROUTINE spline_interpolate

SUBROUTINE spline_hunt(n,x,xp,ig)
  USE kinds
  IMPLICIT NONE
  INTEGER :: ig,n
  REAL(kind=dp), POINTER :: xp(:)
  REAL(kind=dp) :: x
  INTEGER :: i

  REAL(kind=dp) :: dist
!  ig=2 !!DBG BY PG 
!         S1          S2         S3
!     t1        t2          t3        t4
!     t0        t1          t2        t3
!    *1*        *2*        *3*       *4*  
  ig=1
  outer:DO i=2,n
     IF(x<xp(i))THEN
       ig=i-1
       EXIT outer  
     END IF
  END DO outer
  !print *, "|debug|", i, x, xp(i), ig!!DBG BY PG
  IF(ig > n)STOP "error in spline grids"
END SUBROUTINE spline_hunt

END MODULE spline



