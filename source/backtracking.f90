MODULE backtracking 

CONTAINS
  
  SUBROUTINE line_search_backtrack(x_old,f_old,g,b,step_max,x)
    !CALL line_search_backtrack(xold,energy_old,g,b,step_max,x,converged)
!
    USE kinds, ONLY : dp
    USE system_data_types
    USE function_val, ONLY : func
    IMPLICIT NONE
!
!   x_old: current positions
!   g_old: current gradients
!   p: direction vector
!   x: new positions with new lambda
!   step_max: maximum step size allowed
!
    REAL(KIND=dp), POINTER :: x_old(:),x(:),b(:,:),g(:)
    REAL (KIND=dp), INTENT (IN) :: f_old,step_max
    LOGICAL :: converged

    REAL(KIND=dp), PARAMETER :: alpha=1.0E-3, lambda_min=1.0E-6
    REAL(KIND=dp) :: slope,lambda,lambda2,f,f2,a1,b1,lambda_new
    REAL(KIND=dp) :: dum1,dum2,dum3
    INTEGER :: niter,n
    INTEGER, PARAMETER :: maxiter = 9999
    LOGICAL :: done_qudratic
    REAL(KIND=dp), EXTERNAL :: DDOT
    REAL(KIND=dp), ALLOCATABLE :: p(:),g_old(:)
    n=atom_t
    ALLOCATE(p(3*n),g_old(3*n))
!    write(6,*)"#8",g!,converged


   
    p(:)=-g(:)
    CALL DGEMV('N',3*n,3*n,1._dp,b,3*n,p,1,0._dp,g_old,1)
!   CALL DGEMV('N',3*n,3*n,-1._dp,b,3*n,g,1,0._dp,g_old,1)
!   dum1=DDOT(3*n,p,1,p,1)
    dum1=DDOT(3*n,g_old,1,g_old,1)
    dum1=SQRT(dum1)

    !normalizing the step size if it is large
    !IF(dum1>step_max) p(:) = p(:)*step_max/dum1  
    IF(dum1>step_max) g_old(:) = g_old(:)*step_max/dum1  

!   slope=DDOT(3*n,p,1,g_old,1)
    slope=DDOT(3*n,p,1,g_old,1)

    IF(slope>0._dp)PRINT *, "Warning! slope > 0 in linesearch"

    lambda=1._dp !update with the Netwon step (full)
    niter=0
    converged=.FALSE.
    done_qudratic=.FALSE.
    iteration: DO 
       niter=niter+1
       !update x 
       !x(:)=x_old(:)+lambda*p(:)
       x(:)=x_old(:)+lambda*g_old(:)
       write(43,*)"#",x(1),lambda,g_old(1),func(x),f_old+alpha*lambda*slope
       !get function value at the new x
       f=func(x)
       !Doing checks
       IF(f <= f_old+alpha*lambda*slope)THEN
!          PRINT *, "           |bt| doing full newton; return"
          EXIT iteration  !There is sufficient decrease in the functional value for lambda=1 
       ELSE ! backtrack
          IF(.NOT.done_qudratic)THEN !Quadratic approx for the first backtrack
!            print *, "          |bt| trying quadratic "
            done_qudratic=.TRUE.
            dum1=f-f_old-slope 
            lambda_new = -slope/2._dp/dum1
!            print *, "          |bt| qudratic| lambda_new=", lambda_new
          ELSE !Cubic approx for the later ones
            !print *, "trying cubic "
            dum1=f-f_old-slope*lambda
            dum2=f2-f_old-slope*lambda2 
            dum3=1._dp/(lambda-lambda2)

            a1=dum1/lambda**2 - dum2/lambda2**2
            a1=a1*dum3

            b1=-lambda2/lambda**2*dum1+lambda/lambda2**2*dum2
            b1=b1*dum3          

            dum1=-b1**2-3.0_dp*a1*slope
            IF(dum1<0._dp)THEN
               lambda_new=0.5d0*lambda !brute-force to avoid complex numbers
            ELSE 
               dum2=-b1+SQRT(dum1)
               lambda_new=dum2/3._dp/a1
            END IF
!            print *, "                 |bt| cubic| lambda_new1=", lambda_new
            IF(lambda_new>0.5_dp*lambda)lambda_new=lambda   !lambda <= 0.5 lambda1
!            print *, "                 |bt| cubic| lambda_new changed =", lambda_new
          END IF
       END IF
       lambda2=lambda
       f2=f
       lambda=MAX(lambda_new,0.1_dp*lambda)

       IF(niter>maxiter)THEN
         PRINT *, "Warning! Linesearch has reached maximum iteration",niter,maxiter 
         EXIT iteration
       END IF


    END DO iteration
!       WRITE(111,"(i10,4f16.8)")niter,x(1:2),f,lambda 
    DEALLOCATE(p,g_old)

    RETURN
  END SUBROUTINE line_search_backtrack

    SUBROUTINE CPMD_TRACK(X_PAR,DX_PAR,HESSI,N_DIM,STP_FAC,g_step)
    USE kinds, ONLY : dp
    USE system_data_types
!    USE function_val, ONLY : func
    IMPLICIT NONE
    
      INTEGER N_DIM
      REAL(KIND=DP),INTENT(INOUT)::  X_PAR(N_DIM)
      REAL(KIND=DP) DX_PAR(N_DIM),HESSI(N_DIM,N_DIM),STP_FAC!,SCR(9),dtbm(*)

      REAL(kind=dp)    EPS
      PARAMETER (EPS=1.D-6)
     REAL(kind=dp) HESSVD(N_DIM,N_DIM),S_FION(N_DIM),SVAL(N_DIM),SCR(N_DIM*N_DIM)!DTBYM(N_DIM*N_DIM)
     integer g_step
      !CHARACTER TAG*30
      INTEGER   I,J,K,L,NRANK,INFO
      REAL(kind=dp)    SDLEN,DDOT,QNRLEN,TR,SS
      EXTERNAL  DDOT
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(g_step==0)then
      L=0
      DO i=1,sp_t
        DO J=1,atom_p_sp(i)
          DO K=1,3
              L=L+1
             ! DTBYM(L)=5.D0*5.D0/(atwt(z(i))*1822.888485D0)!DT2BYM(IS)!*LSKCOR(K,IAT)
          ENDDO
        ENDDO
      ENDDO           
                    endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      CALL DSCAL(N_DIM,STP_FAC,DTBYM(1),1)
      CALL DCOPY(N_DIM*N_DIM,HESSI(1,1),1,HESSVD(1,1),1)
      
      DO I=1,N_DIM
        S_FION(I)=-DX_PAR(I)
        SCR(I)=DTbYM(I)*DX_PAR(I)
        WRITE(*,*)DTBYM(I),stp_fac
      ENDDO
  
      SDLEN=SQRT(DDOT(N_DIM,SCR,1,SCR,1))
      CALL DGELSS(N_DIM,N_DIM,1,HESSVD,N_DIM,S_FION,N_DIM,SVAL,EPS,NRANK,SCR,6*N_DIM,INFO)
      QNRLEN=SQRT(DDOT(N_DIM,S_FION,1,S_FION,1))
      TR=100.0D0 * SDLEN
      IF(TR.LT.0.01D0) TR=0.01D0
      TR=MIN(TR,0.3D0)
      IF(QNRLEN.GT.TR) THEN
        SS=SQRT(TR/QNRLEN)
        CALL DSCAL(N_DIM,SS,S_FION,1)
      ENDIF
      CALL DAXPY(N_DIM,1.0D0,S_FION,1,X_PAR(1),1)
      RETURN
   END SUBROUTINE 
   END MODULE backtracking
