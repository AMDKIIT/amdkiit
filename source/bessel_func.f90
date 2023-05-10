MODULE bessel_func
CONTAINS
 FUNCTION sph_bes(n,xx) RESULT(sph_besl)!(bes,besp,bespp,l,xx)

 USE system_data_types
 USE kinds
 USE constants

 implicit none
 integer :: n
 real(kind=dp) :: xx
 real(kind=dp) :: sph_besl

 integer,parameter :: kmax=40
 real(kind=dp),parameter :: toler=1.d-15

 real(kind=dp) :: fact,jl,jlp,jr
 integer :: k,id

 if (abs(xx)<toler) then
   sph_besl=0
   if (n==0) sph_besl=1
   return
 end if

 jlp=0
 
 if (xx<1) then
   fact=1
   do id=1,n
     fact=fact*xx/dble(2*id+1)
   end do
   jl=1
   jr=1
   k=0
   do while(abs(jr)>=toler.and.k<kmax)
     k=k+1
     jr=-jr*(0.5*xx*xx)/dble(k*(2*(n+k)+1))
     jl=jl+jr
   end do
   sph_besl=jl*fact
   if (abs(jr)>toler)  print*,'Bessel function converge!'
 else
   jl =sin(xx)/xx
   jlp=(-cos(xx)+jl)/xx
   do id=2,n+1
     jr=-jl+dble(2*id-1)*jlp/xx
     jl=jlp
     jlp=jr
   end do
   sph_besl=jl
 end if
 END FUNCTION sph_bes

 SUBROUTINE ov_sph_bes_f(RG,dlog_c,NSIZE,FUN,L,G,RGMAX,OVLP)
 USE system_data_types
 USE kinds

 IMPLICIT NONE
 REAL(KIND=DP), INTENT(IN) :: RG, dlog_c, RGMAX, G
 INTEGER, INTENT(IN) :: L,NSIZE
 REAL(KIND=DP), INTENT(IN) :: FUN(NSIZE)
 REAL(KIND=DP), INTENT(OUT) :: OVLP

 INTEGER, PARAMETER ::  tot_n_min = 400
 REAL(KIND=DP), PARAMETER :: conf(4) = [17.D0, 59.D0, 43.D0,49.D0]/48.D0
 INTEGER :: tot_n, i
 REAL(KIND=DP) :: R, X, FVAL

      tot_n = MAX(tot_n_min, INT(20.D0 * G / 6.D0))
      OVLP=0
      R=rgmax

      DO I=1,tot_n-3
        IF(I.NE.1)R=RGMAX*DBLE(I-1)/tot_n
        X=1.D0+DLOG(R/RG) /dlog_c
        CALL get_interpolated_value(NSIZE,FUN,X,FVAL)
        IF(I.LE.4) THEN
          OVLP=OVLP+RGMAX*conf(I)*FVAL/tot_n*sph_bes(L,R*G)*R**2
                  IF(I.NE.1)THEN
                     R=RGMAX-R
                     X=1.D0+DLOG(R/RG) /dlog_c
                     CALL get_interpolated_value(NSIZE,FUN,X,FVAL)
                     OVLP=OVLP+RGMAX*conf(I)*FVAL/tot_n*sph_bes(L,R*G)*R**2
                  ENDIF
        ELSE
          OVLP=OVLP+RGMAX*FVAL/tot_n*sph_bes(L,R*G)*R**2
        ENDIF
      ENDDO

 RETURN
 END 

 SUBROUTINE get_interpolated_value(nval,func,at_x,rval)
 USE system_data_types
 USE kinds

 IMPLICIT NONE

 INTEGER, INTENT(IN) :: nval
 REAL(KIND=DP), INTENT(IN) :: FUNc(Nval),at_x
 REAL(KIND=DP), INTENT(OUT) :: rval

 REAL(KIND=DP) HX(4),R(3),POLY(2)
 INTEGER INCR,I

  INCR= MIN(MAX(2, INT(at_x)), nval-2)

  DO I=1,4
     HX(I)=at_X-DBLE(INCR-1)-DBLE(I-1)
  ENDDO
 R(1:3) = -HX(2:4) * Func(INCR-1:INCR+1) + HX(1:3) * Func(INCR:INCR+2)
 poly(1:2) = 0.5D0 * (-HX(3:4) * R(1:2) + HX(1:2) * R(2:3))
 rval=1.D0/3.D0*(-HX(4)*POLY(1)+HX(1)*POLY(2))

 RETURN
 END

END MODULE BESSEL_func
