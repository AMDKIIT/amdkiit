    SUBROUTINE apply_pbc(X1,Y1,Z1,X2,Y2,Z2,M,A,IBRAV)
      USE system_data_types
      USE kinds
      USE CONSTANTS
      USE math,  ONLY : inverse_matrix,sort_array2
      IMPLICIT NONE


      REAL(kind=dp)  X1,Y1,Z1,X2,Y2,Z2,A,s2,x2t,y2t,xt,r(4),xic,s(3),m_a,a_pbc(4),dx,y0
      INTEGER M,IBRAV,i,kk
      REAL(kind=dp)  aby2!,theta,ALPHA_a,BETA_a,GAMMA_a
      REAL(KIND=dp),    DIMENSION(:,:),       POINTER :: test_1
      REAL(kind=dp) bbbb,cccc,xl(4), v1crys(4),v2crys(4),v3crys(4),aux
      INTEGER indx(4),xl_int(4)
      REAL(kind=dp) x0,x3,y3
       
      allocate(test_1(3,3))
      a_pbc=0.0d0
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
       CALL dgemv('T',3,3,1.0d0,test_1,3,r,1,0.0d0,s,1)
       xic=REAL(m,kind=dp)
       s(1)=s(1)-NINT(s(1)/xic)*xic
       s(2)=s(2)-NINT(s(2)/xic)*xic
       s(3)=s(3)-NINT(s(3)/xic)*xic
       CALL dgemv('T',3,3,1.0d0,TRANSPOSE(a_cell),3,s,1,0.0d0,r,1)
       x2=r(1)
       y2=r(2)
       z2=r(3)
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
     case(4)
          a_pbc(1)=a/2.0d0
          a_pbc(2)=a*SQRT(3.d0)/4.d0
          a_pbc(3)=a/2.d0*primitive_vec(3)/primitive_vec(1)
          a_pbc(4)=SQRT(3.0d0)/3.0d0
          IF(a_pbc(2)==0)THEN
               r(1)=x1
               r(2)=y1
               r(3)=z1
               test_1=(b_cell/a_lattice)
               CALL dgemv('T',3,3,1.0d0,test_1,3,r,1,0.0d0,s,1)
               xic=REAL(m,kind=dp)
               s(1)=s(1)-NINT(s(1)/xic)*xic
               s(2)=s(2)-NINT(s(2)/xic)*xic
               s(3)=s(3)-NINT(s(3)/xic)*xic
               CALL dgemv('T',3,3,1.0d0,TRANSPOSE(a_cell),3,s,1,0.0d0,r,1)
               x2=r(1)
               y2=r(2)
               z2=r(3)
               RETURN
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
       RETURN

      case(14)
          IF(ABS(COS((a_theta(1))*pi/180.0d0)).LT.1.e-10.AND.ABS(COS((a_theta(2))*pi/180.0d0)).LT.1.e-10) THEN
             a_pbc(1)=a/2.0D0
             a_pbc(2)=a/2.D0*primitive_vec(2)/primitive_vec(1)*SQRT(1.d0-(COS((a_theta(2))*pi/180.0d0))**2)
             a_pbc(3)=a/2.D0*primitive_vec(3)/primitive_vec(1)
             a_pbc(4)=-(COS((a_theta(2))*pi/180.0d0))/SQRT(1.D0-(COS((a_theta(2))*pi/180.0d0))**2)
          ELSE
             a_pbc(2)=0.0d0
               r(1)=x1
               r(2)=y1
               r(3)=z1
               test_1=(b_cell/a_lattice)
               CALL dgemv('T',3,3,1.0d0,test_1,3,r,1,0.0d0,s,1)
               xic=REAL(m,kind=dp)
               s(1)=s(1)-NINT(s(1)/xic)*xic
               s(2)=s(2)-NINT(s(2)/xic)*xic
               s(3)=s(3)-NINT(s(3)/xic)*xic
               CALL dgemv('T',3,3,1.0d0,TRANSPOSE(a_cell),3,s,1,0.0d0,r,1)
               x2=r(1)
               y2=r(2)
               z2=r(3)
               RETURN
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
       RETURN
     case(12)
          !      NEW GEOMETRY PROGRAMMED                  
          bbbb=a*(primitive_vec(2)/primitive_vec(1))*COS((a_theta(1))*pi/180.0d0)
          cccc=a*(primitive_vec(2)/primitive_vec(1))*SQRT(1.0d0-(COS((a_theta(1))*pi/180.0d0))**2)
                  
          ! LOOK FOR THE CLOSEST IMAGES TO BUILD THE WIGNER-SEITZ CELL.
          DO kk=1,4
             xl(kk)=(a-REAL(kk-1,kind=dp)*bbbb)**2+(REAL(kk-1,kind=dp)*cccc)**2
          ENDDO

          ! ORDER THE IMAGES WITH INCREASING DISTANCES
          CALL sort_array2(4,xl,indx)
                  
          ! HERE WE DEFINE THE DIMENSIONS OF THE W-S CELL. APBC(1),APBC(2)
          ! AND APBC(4) ARE IN THE (X,Y) PLANE, WHERE THE CELL IS A NON-
          ! REGULAR HEXAGON; WHILE APBC(3) CORRESPONDS TO THE Z-DIRECTION.
          ! THE HEXAGON IS DEFINED BY THE 3 CLOSEST IMAGES IN THE (X,Y)
          ! PLANE. TWO OF THEM WERE CALCULATED BEFORE, AND THE 3RD (THE
          ! SMALLEST ONE) IS THAT IN THE DIRECTION OF THE LATTICE VECTOR
          ! NUMBER 2 (CELLDM(2)).
                  
          a_pbc(1)=SQRT(xl(1))/2.0d0
          a_pbc(2)=a/2.0d0*primitive_vec(2)/primitive_vec(1)
          a_pbc(3)=a/2.0d0*primitive_vec(3)/primitive_vec(1)
          a_pbc(4)=SQRT(xl(2))/2.0d0
          
          !============================================================
           v1crys(1)=a_cell(1,1)-REAL(indx(1)-1,kind=DP)*a_cell(1,2)
           v1crys(2)=-REAL(indx(1)-1,kind=DP)*a_cell(2,2)
            aux=SQRT(v1crys(1)**2+v1crys(2)**2)
               v1crys(3)=v1crys(1)/aux
               v1crys(4)=v1crys(2)/aux
         v2crys(1)=a_cell(1,2)
         v2crys(2)=a_cell(2,2)
         aux=SQRT(v2crys(1)**2+v2crys(2)**2)
         v2crys(3)=v2crys(1)/aux
        v2crys(4)=v2crys(2)/aux
        v3crys(1)=a_cell(1,1)-REAL(indx(2)-1,kind=DP)*a_cell(1,2)
        v3crys(2)=-REAL(indx(2)-1,kind=DP)*a_cell(2,2)
        aux=SQRT(v3crys(1)**2+v3crys(2)**2)
        v3crys(3)=v3crys(1)/aux
                v3crys(4)=v3crys(2)/aux
          !============================================================
   
          DO i=1,4
          r(i)=a_pbc(i)*m
          END DO
    z2=dmod(z1,2*r(3))
    z2=z2-2*r(3)*dint(z2/r(3))
    CALL rot(x1,y1,x0,y0,v1crys,1)
    x0=dmod(x0,2*r(1))
    x0=x0-2*r(1)*dint(x0/r(1))
    CALL rot(x0,y0,x3,y3,v1crys,-1)
    CALL rot(x3,y3,x0,y0,v2crys,1)
    x0=dmod(x0,2*r(2))
    x0=x0-2*r(2)*dint(x0/r(2))
    CALL rot(x0,y0,x3,y3,v2crys,-1)
    CALL rot(x3,y3,x0,y0,v3crys,1)
    x0=dmod(x0,2*r(4))
    x0=x0-2*r(4)*dint(x0/r(4))
    CALL rot(x0,y0,x3,y3,v3crys,-1)
    CALL rot(x3,y3,x0,y0,v1crys,1)
    x0=dmod(x0,2*r(1))
    x0=x0-2*r(1)*dint(x0/r(1))
    CALL rot(x0,y0,x2,y2,v1crys,-1)
    RETURN
       case default
       IF(IONODE)WRITE(6,*)"P.B.C. is not implemented for"
       STOP
       END select
    deallocate(test_1)
    END SUBROUTINE

  SUBROUTINE rot(x1,y1,x0,y0,vv,isign)
    ! ==--------------------------------------------------------------==
    ! == THIS ROUTINE TRANSFORMS THE POINT (X1,Y1) INTO THE POINT     ==
    ! == (X0,Y0) CORRESPONDING TO A ROTATION OF AN ANGLE WHICH COS    ==
    ! == AND SIN ARE GIVEN BY VV(3) AND VV(4) RESPECTIVELY. THE       ==
    ! == SENSE OF THE ROTATION IS GIVEN BY ISIGN (+ FOR CLOCKWISE,    ==
    ! == AND - FOR ANTICLOCKWISE)                                     ==
    ! ==--------------------------------------------------------------==
    USE kinds
    REAL(kind=dp)                             :: x1, y1, x0, y0, vv(4)
    INTEGER                                  :: isign

    x0=x1*vv(3)+y1*vv(4)*REAL(isign,kind=dp)
    y0=-x1*vv(4)*REAL(isign,kind=dp)+y1*vv(3)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rot

