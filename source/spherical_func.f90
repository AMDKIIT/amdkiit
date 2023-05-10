      SUBROUTINE compute_yL(L,ng,H_G,G,y_L)
      !     REAL SPHERICAL HARMONICS,  
      !     (L=1,2...10), ORDER:  S, P_X, P_Z, P_Y, D_XY, D_XZ, D_Z^2, D_YZ, D_X^2-Y^2  ....
      !     FOR IRREDUCIBLE REPRESENTATIONS OF THE GROUP O

      USE system_data_types
      USE kinds
      USE constants

      IMPLICIT NONE
      INTEGER,intent(in) ::L,ng
      REAL(kind=dp)  h_g(ng),G(3,*),ylm(10)
      REAL(KIND=DP),INTENT(OUT)::Y_l(*)
      INTEGER ig
      REAL(kind=dp), PARAMETER                  :: thr = 1.0e-6_dp
      DO ig=1,ng
         IF (ABS(h_g(ig)).GE.thr) THEN
            SELECT CASE (L)
          CASE (1)
            ylm(1) = SQRT(1.0_dp/fourpi) 
          CASE (2)
            ylm(2) = SQRT(3._dp/fourpi)*g(1,ig)/SQRT(h_g(ig))
          CASE (3)
            ylm(3) = SQRT(3._dp/fourpi)*g(3,ig)/SQRT(h_g(ig))
          CASE (4)
            ylm(4) = SQRT(3._dp/fourpi)*g(2,ig)/SQRT(h_g(ig))
          CASE (5)
            ylm(5) = SQRT(15._dp/fourpi)*g(1,ig)*g(2,ig)/h_g(ig) 
          CASE (6)
            ylm(6) = SQRT(15._dp/fourpi)*g(1,ig)*g(3,ig)/h_g(ig)
          CASE (7)
            ylm(7) = SQRT(5._dp/fourpi/4._dp)*(3._dp*g(3,ig)**2/h_g(ig)-1._dp)
          CASE (8)
            ylm(8) = SQRT(15._dp/fourpi)*g(2,ig)*g(3,ig)/h_g(ig)
          CASE (9)
            ylm(9)=  SQRT(15._dp/fourpi/4._dp)*(g(1,ig)**2-g(2,ig)**2)/h_g(ig)
          CASE(10)
            ylm(10)= SQRT(7._dp/fourpi)*5._dp/2._dp*g(1,ig)* &
                      (g(1,ig)**2-0.6_dp*h_g(ig))/(h_g(ig)*SQRT(h_g(ig)))
          END SELECT 
             y_l(ig) =  ylm(L)
         ELSE
             y_l(ig) = 0.0_dp
             h_g(ig) = 0.0_dp
             if(l==1)y_l(ig)=SQRT(1.0_dp/fourpi)
         ENDIF
       ENDDO
       RETURN
      END SUBROUTINE compute_yL

