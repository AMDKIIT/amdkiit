MODULE fft_support 
  CONTAINS
!=----------------------------------------------------------------------=!
!
!         FFT support Functions/Subroutines
!
!=----------------------------------------------------------------------=!

  FUNCTION optimal_fft_grid_fft(n) RESULT(n_out)
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: n
    INTEGER :: n_out,i,ng(3)
    INTEGER, PARAMETER :: limit_prefer2=10 !Prefer 2^n if n_out is close to it 
                                           !as powers of 2 is fast with FFTW
!   TODO: Find a number close to 2^a 3^b 5^c 7^d

    ng(1)=close_to_power_m(n,2) !2^m
    ng(2)=close_to_power_m(n,3) !3^m
    ng(3)=close_to_power_m(n,5) !7^m

    n_out=MINVAL(ng)
    IF((n_out-ng(1))<limit_prefer2)n_out=ng(1)
    !print *, "ng =", ng(1:3), " nout", n_out

  END FUNCTION optimal_fft_grid_fft

  FUNCTION close_to_power_m(n,m) RESULT(nn)
    IMPLICIT NONE
    INTEGER :: n, m
    INTEGER :: nn
    INTEGER :: i
    INTEGER, PARAMETER :: imax=999

    i=1
    power_m: DO
      IF(n == m**i)THEN
        nn=m**(i)
        EXIT power_m
      ELSE IF(n > m**i .and. n < m**(i+1) )THEN
        nn=m**(i+1)
        EXIT power_m
      ELSE
        i=i+1
        IF(i>imax)STOP 'ERROR IN CLOSE_TO_POWER_M'
      END IF
    END DO power_m

  END FUNCTION close_to_power_m

  FUNCTION permitted (n_val)
  ! find if the fft dimension is a good one
  ! a "bad one" is either not implemented (as on IBM with ESSL)
  ! or implemented but with awful performances (most other cases)
  implicit none
  integer :: n_val

  logical :: permitted
  integer :: power (5)
  integer :: m_val, i, fac, p, max_power
  integer :: factors( 5 ) = (/ 2, 3, 5, 7, 11 /)

  ! find the factors of the fft dimension

  m_val  = n_val
  power = 0
  factors_loop: do i = 1, 5
     fac = factors (i)
     max_power = NINT ( LOG( DBLE (m_val) ) / LOG( DBLE (fac) ) ) + 1
     do p = 1, max_power
        if ( m_val == 1 ) EXIT factors_loop
        if ( MOD (m_val, fac) == 0 ) then
           m_val = m_val / fac
           power (i) = power (i) + 1
        endif
     enddo
  end do factors_loop

  IF ( n_val /= ( m_val * 2**power (1) * 3**power (2) * 5**power (3) * 7**power(4) * 11**power (5) )) write(6,*)"ERROR ORDER"

  IF ( m_val /= 1 ) then
     permitted = .false.
  ELSE
     permitted = (( power(4) == 0) .and. (power(5) == 0))
  ENDIF
  RETURN
  END FUNCTION permitted 

   INTEGER FUNCTION optimal_fft_order( n_old, n_order )
  !
  !    This function find a "good" fft order value greater or equal to "n_old"
  !
  !    n_old  (input) tentative order n of a fft
  !
  !    n_order  (optional input) if present restrict the search of the order
  !        in the ensemble of multiples of n_order
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: n_old
  INTEGER, OPTIONAL, INTENT(IN) :: n_order
  INTEGER :: n_new

  n_new = n_old
  IF( PRESENT( n_order ) ) THEN
    IF (n_order <= 0 .OR. n_order > n_old) WRITE(6,*)"ERROR FFT ORDER "
    DO WHILE( ( ( .NOT. permitted( n_new ) ) .OR. ( MOD( n_new, n_order ) /= 0 )) .AND. ( n_new <= 2049) )
      n_new = n_new + 1
    END DO
  ELSE
    DO WHILE( ( .NOT. permitted( n_new ) ) .AND. ( n_new <= 2049 ) )
      n_new = n_new + 1
    END DO
  END IF

  IF( n_new > 2049 ) WRITE(6,*)"ERROR FFT ORDER "
  optimal_fft_order = n_new

  RETURN
  END FUNCTION optimal_fft_order
END MODULE fft_support
