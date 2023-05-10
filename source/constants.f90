MODULE constants
  USE kinds
  REAL(KIND=dp), PARAMETER :: pi=4.0_dp*ATAN(1.0_dp)
  REAL(KIND=dp), PARAMETER :: twopi=8.0_dp*ATAN(1.0_dp)
  REAL(KIND=dp), PARAMETER :: fourpi=16.0_dp*ATAN(1.0_dp)
  REAL(KIND=dp), PARAMETER :: eps8=1.0E-8_dp
  REAL(KIND=dp), PARAMETER :: gc_cutoff=1.0e-8_dp !Density cutoff for the calculation of the gradient correction
  REAL(KIND=DP), PARAMETER :: e_2 = 2.0_DP      ! the square of the electron charge
  REAL(KIND=dp), PARAMETER :: chk_wf=1.0E-4_dp
  REAL(KIND=dp), PARAMETER :: lambda=.400
  REAL(KIND=dp), PARAMETER :: hcut=0.50
     COMPLEX*16, PARAMETER :: uimag=DCMPLX(0.D0,1.D0)
  REAL(KIND=dp), PARAMETER :: dmax=0.3D0 
  
END MODULE constants
