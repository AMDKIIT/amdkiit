    SUBROUTINE apply_pbc(X1,Y1,Z1,X2,Y2,Z2,M,A,IBRAV)
      USE system_data_types
      USE kinds
      USE CONSTANTS
      USE math,  ONLY : inverse_matrix
      IMPLICIT NONE


      REAL(kind=dp)  X1,Y1,Z1,X2,Y2,Z2,A,s2,x2t,y2t,xt,r(3),xic,s(3),m_a,a_pbc(4),dx,y0
      INTEGER M,IBRAV,i
      REAL(kind=dp)  aby2,theta,ALPHA_a,BETA_a,GAMMA_a
      REAL(KIND=dp),    DIMENSION(:,:),       POINTER :: test_1
      allocate(test_1(3,3))
      theta=90.0d0
      ALPHA_a=COS((theta)*pi/180.0d0)
      theta=90.0d0
      BETA_a=COS((theta)*pi/180.0d0)
      theta=90.0D0
      GAMMA_a=COS((theta)*pi/180.0d0)

      !write(*,*)"**",x1,y1,z1
      aby2=(A*M)/2.d0
      !Simple Cubic 
      select case (IBRAV)
      case(1) 
       X2=DMOD(X1,A)
       Y2=DMOD(Y1,A)
       Z2=DMOD(Z1,A)
       X2=X2-A*DINT(X2/aby2)
       Y2=Y2-A*DINT(Y2/aby2)
       Z2=Z2-A*DINT(Z2/aby2)
       RETURN
       case(2)
       m_a=m*a*SQRT(2.0d0)/4.0d0
       s2=SQRT(2.0d0)
       x2t=( x1+y1)/s2
       y2t=(-x1+y1)/s2
       x2=x2t
       y2=y2t
       x2=dmod(x2,2.d0*m_a)
       y2=dmod(y2,2.d0*m_a)
       z2=dmod(z1,2.d0*m_a*s2)
       x2=x2-2.d0*m_a*dint(x2/m_a)
       y2=y2-2.d0*m_a*dint(y2/m_a)
       z2=z2-2.d0*m_a*s2*dint(z2/(s2*m_a))
       IF (ABS(x2)+ABS(y2)+s2*ABS(z2).LT.2.d0*m_a) GOTO 11
       x2=x2-dsign(m_a,x2)
       y2=y2-dsign(m_a,y2)
       z2=z2-dsign(s2*m_a,z2)
       11  CONTINUE
       xt=(x2-y2)/s2
       y2=(x2+y2)/s2
       x2=xt
       RETURN
       case(3,5,6,7,9,10,11,13)
       r(1)=x1
       r(2)=y1
       r(3)=z1
       test_1=(b_cell/a_lattice)
       print*,"#",a_cell(1,1),a_cell(1,2),a_cell(1,3)
       print*,"#",test_1(1,1),test_1(1,2),test_1(1,3)
       CALL dgemv('T',3,3,1.0d0,test_1,3,r,1,0.0d0,s,1)
       xic=REAL(m,kind=dp)
       s(1)=s(1)-NINT(s(1)/xic)*xic
       s(2)=s(2)-NINT(s(2)/xic)*xic
       s(3)=s(3)-NINT(s(3)/xic)*xic
       CALL dgemv('T',3,3,1.0d0,TRANSPOSE(a_cell),3,s,1,0.0d0,r,1)
       x2=r(1)
       y2=r(2)
       z2=r(3)
       print*,"xyz",x2,y2,z2
       return 
       case(8)
        a_pbc(1)=a/2.d0
        a_pbc(2)=a/2.d0*primitive_vec(2)
        a_pbc(3)=a/2.d0*primitive_vec(3)

       DO i=1,3
       r(i)=a_pbc(i)*m
       ENDDO
          x2=dmod(x1,2*r(1))
          y2=dmod(y1,2*r(2))
          z2=dmod(z1,2*r(3))
        x2=x2-2*r(1)*dint(x2/r(1))
        y2=y2-2*r(2)*dint(y2/r(2))
        z2=z2-2*r(3)*dint(z2/r(3))
       RETURN
     case(4,14)
          a_pbc(1)=a/2.0d0
          a_pbc(2)=a*SQRT(3.d0)/4.d0
          a_pbc(3)=a/2.d0*primitive_vec(3)
          a_pbc(4)=SQRT(3.0d0)/3.0d0
    
          IF(ABS(AlpHA_A).LT.1.e-10 .AND.ABS(BETA_a).LT.1.e-10) THEN
             a_pbc(1)=a/2.0D0
             a_pbc(2)=a/2.D0*primitive_vec(2)*SQRT(1.d0-GAMMA_A**2)
             a_pbc(3)=a/2.D0*primitive_vec(3)
             a_pbc(4)=-GAMMA_A/SQRT(1.D0-GAMMA_A**2)
          ELSE
             a_pbc(2)=0.0d0
          ENDIF

           DO i=1,3
           r(i)=a_pbc(i)*m
           ENDDO
      y0=y1
      y2=dmod(y1,2*r(2))
      z2=dmod(z1,2*r(3))
      y2=y2-2*r(2)*dint(y2/r(2))
      z2=z2-2*r(3)*dint(z2/r(3))
      dx=(y2-y0)*a_pbc(4)           
      x2=dmod(x1-dx,2*r(1))
      x2=x2-2*r(1)*dint(x2/r(1))
      !write(*,*)"#",x2,y2,z2
      RETURN

!                DO i=1,4
!       b(i)=2*a(i)*m
!       c(i)=a(i)*m
!    END DO
!    z2=dmod(z1,b(3))
!    z2=z2-b(3)*dint(z2/c(3))
!    CALL rot(x1,y1,x0,y0,bc_com%v1crys,1)
!    x0=dmod(x0,b(1))
!    x0=x0-b(1)*dint(x0/c(1))
!    CALL rot(x0,y0,x3,y3,bc_com%v1crys,-1)
!    CALL rot(x3,y3,x0,y0,bc_com%v2crys,1)
!    x0=dmod(x0,b(2))
!    x0=x0-b(2)*dint(x0/c(2))
!    CALL rot(x0,y0,x3,y3,bc_com%v2crys,-1)
!    CALL rot(x3,y3,x0,y0,bc_com%v3crys,1)
!    x0=dmod(x0,b(4))
!    x0=x0-b(4)*dint(x0/c(4))
!    CALL rot(x0,y0,x3,y3,bc_com%v3crys,-1)
!    CALL rot(x3,y3,x0,y0,bc_com%v1crys,1)
!    x0=dmod(x0,b(1))
!    x0=x0-b(1)*dint(x0/c(1))
!    CALL rot(x0,y0,x2,y2,bc_com%v1crys,-1)
!    RETURN
!       case(3)
!         xIC = REAL( M, DP )
!        X2 = DNINT(X1/xIC)*xIC-x1
!        Y2 = DNINT(Y1/xIC)*xIC-y1
!        Z2 = DNINT(Z1/xIC)*xIC-z1
!        write(222,*)"#",DNINT(X1/xIC)*xIC
       case default
       IF(IONODE)WRITE(6,*)"P.B.C. is not implemented for"
       STOP
       END select
    deallocate(test_1)
    END SUBROUTINE
