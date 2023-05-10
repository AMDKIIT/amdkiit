SUBROUTINE spline_inter(n,xpp,ypp,intp_x1,intp_y,dx,intp_n)
USE spline
USE system_data_types
USE kinds

IMPLICIT NONE
  integer :: n
  REAL(kind=dp):: xpp(n),ypp(n),intp_x1
  REAL(kind=dp),intent(out):: intp_y(intp_n)
  REAL(kind=dp), POINTER :: xp(:),yp(:)
  REAL(kind=dp), POINTER :: h(:),v(:),u(:),sdd(:)
  REAL(kind=dp) :: x, y, dx

  INTEGER ::  i, j, ig, intp_n
  ALLOCATE(xp(n),yp(n))
  ALLOCATE(h(n),v(n),u(n),sdd(n))
!
  DO i=1,n
  xp(i)=xpp(i)
  yp(i)=ypp(i)
  END DO
!Setup the spline parameters for the given data
  CALL spline_params(n,xp,yp,h,v,u)

!Perform Gausselimination of tridiagonal method (Thomas Method)
  CALL spline_thomas(n,h,v,u,sdd)
  x=intp_x1
  grid: DO i=1,intp_n
    IF(x>xp(n))EXIT grid
    CALL spline_hunt(n,x,xp,ig)
    CALL spline_interpolate(n,xp,yp,h,v,u,sdd,ig,x,y)
    intp_y(i)=y
    x=x+dx
  END DO grid
END SUBROUTINE spline_inter

