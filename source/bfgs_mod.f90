MODULE bfgs_mod

  CONTAINS

  SUBROUTINE GEOOPT_bfgs(opt_step,x,g,b,n,chk_cnv)
  USE kinds
  USE system_data_types
  use mympi
  USE backtracking 
  IMPLICIT NONE
  INTEGER, INTENT(IN):: opt_step,n
  REAL(KIND=dp),pointer,intent(inout) :: x(:),g(:),b(:,:)
  REAL(KIND=dp),SAVE,pointer :: xold(:),gold(:)
  INTEGER I 
  REAL(KIND=dp)::gnrm,energy_old,dx,dxmax
  REAL(KIND=dp), EXTERNAL :: DDOT
  REAL(KIND=dp), PARAMETER ::step_max=0.5_dp
  LOGICAL, intent(inout) :: chk_cnv
 
  
  IF(.NOT.ASSOCIATED(gold)) ALLOCATE(gold(3*atom_t))
  IF(.NOT.ASSOCIATED(xold)) ALLOCATE(xold(3*atom_t))
  
  chk_cnv=.false.
  if(opt_step==0)then
    !Initialize Hessian Inverse
  b(:,:)=0._dp
  DO i=1,n
    b(i,i)=1._dp
  END DO
  endif
  
  GNMAX=0.0
  DO i=1,n
     IF(ABS(g(i)).GT.GNMAX) GNMAX=ABS(g(i))
  ENDDO
     
   !   BFGS update of inverse Hessian
   IF(opt_step>0) CALL bfgs(x,xold,g,gold,b)
   !   Find the magnitude of gradient vector
   
   !print*,"o",x(1),x(2),x(3)
   gnrm=DSQRT(DDOT(n,g,1,g,1))/DFLOAT(atom_t) !normalized with number of atoms

           IF(gnrm<geo_opt%g_cutoff)then
              chk_cnv=.true.
              IF(IONODE)write(6,"(A)")repeat("*", 93)
              if(ionode) write(6,"(3A)")repeat(" ",35),"OPTIMIZATION COMPLETED",repeat(" ",35) 
              goto 10
           endif
  !setting the search direction for linesearch
   xold(:)=x(:)
   gold(:)=g(:)
   energy_old=etotal
   !:wqprint*,"1gold",gold(1)
     
!  print*,"old=",x
  !Linesearch
  CALL line_search_backtrack(xold,energy_old,g,b,step_max,x)
!  print*,"new=",x
           !DXMAX=0.0D0
           !DO I=1,n
            !  DX=ABS(x(i)-xold(i))
            !  DXMAX=MAX(DX,DXMAX)
           !ENDDO
            !IF(DXMAX<1.0E-6_dp)then
            !    chk_cnv=.TRUE.
            !goto 10
           !end if


  10 continue
  RETURN 
  END SUBROUTINE GEOOPT_bfgs

  SUBROUTINE bfgs(x,xold,g,gold,b)
    USE kinds
    USE system_data_types
    IMPLICIT NONE
    
    REAL(KIND=dp), POINTER :: x(:),xold(:), g(:), gold(:), b(:,:)

    REAL(KIND=dp), ALLOCATABLE :: s(:), y(:), scratch(:,:), bb(:,:), by(:)
    REAL(KIND=dp) :: sy,fac
    REAL(KIND=dp), EXTERNAL :: DDOT
   integer n
   n=atom_t
!    return
    ALLOCATE(s(3*n),y(3*n),scratch(3*n,3*n),bb(3*n,3*n),by(3*n))

    s(:)=x(:)-xold(:)
    y(:)=g(:)-gold(:)
!    y(:)=gold(:)-g(:)

!copy of B matrix
    bb(:,:)=b(:,:)
!s^Ty
    sy=DDOT(3*n,s,1,y,1)
    sy=1._dp/sy
!sy^T/s^Ty
    CALL DGEMM('N','T',3*n,3*n,1  ,sy     ,s      ,3*n,y,3*n,0._dp,scratch,3*n)
!B_k = B_k- sy^T B_k/s^Ty
    CALL DGEMM('N','N',3*n,3*n,3*n,-1._dp  ,scratch,3*n,bb,3*n,1._dp,b     ,3*n)
!ys^T/s^Ty
    CALL DGEMM('N','T',3*n,3*n,1  ,sy     ,y      ,3*n,s,3*n,0._dp,scratch,3*n)
!B_k = B_k - B_k ys^T/s^Ty
    CALL DGEMM('N','N',3*n,3*n,3*n,-1._dp  ,bb,3*n,scratch,3*n,1._dp,b     ,3*n)
!B_k y/s^Ty
    CALL DGEMV('N',3*n,3*n,sy,bb,3*n,y,1,0._dp,by,1)
!1+y^T B_k y/s^Ty
    fac=DDOT(3*n,y,1,by,1)
    fac=1._dp+fac
!ss^T/s^Ty
    CALL DGEMM('N','T',3*n,3*n,1  ,sy     ,s      ,3*n,s,3*n,0._dp,scratch,3*n)
!Bk+(1+y^T B_k y/S^Ty)ss^T/s^Ty
    CALL DAXPY(3*n*3*n,fac,scratch,1,b,1)
    DEALLOCATE(s,y,bb,scratch,by)
  END SUBROUTINE bfgs
  

END MODULE bfgs_mod
