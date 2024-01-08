MODULE set_cell

CONTAINS
  SUBROUTINE set_vec(a,b,c,a_cell,lattice_str)
  USE kinds
  USE CONSTANTS
  USE system_data_types, only : ionode,symmetry,a_theta
  IMPLICIT NONE
  REAL(kind=dp), POINTER :: a_cell(:,:)
  REAL(kind=dp):: a,b,c
  CHARACTER (LEN=30) ::lattice_str 

  a_cell=0.0d0
  
  ! Identify the symmetry number based on the lattice type
    
  select case (trim(adjustl(lattice_str)))
  case ("CUBIC")
     symmetry = 1
  case ("FCC")
     symmetry = 2
  case ("BCC")
     symmetry = 3
  case ("HEXAGONAL")
     symmetry = 4
  !case ("TRIGONAL")
  CASE ("TRIGONAL", "RHOMBOHEDRAL")
     symmetry = 5
  case ("TETRAGONAL")
     symmetry = 6
  case ("BODY CENTERED TETRAGONAL")
     symmetry = 7
  case ("ORTHORHOMBIC")
     symmetry = 8
  case ("MONOCLINIC")
     symmetry = 12
  case ("TRICLINIC")
     symmetry = 14
  case default
     if(ionode)write(*,*) "Invalid lattice type entered! See amdkiit_manual"
     stop
  end select

    IF (symmetry.EQ.1) THEN !Simple Cubic
      a_cell(1,1)=a
      a_cell(2,2)=b
      a_cell(3,3)=c

    ELSEIF(symmetry.EQ.2) THEN !Face Centered Cubic
      a_cell = reshape([ -a/2.D0,   0.D0, a/2.D0, &
                            0.D0, a/2.D0, a/2.D0, &
                         -a/2.D0, a/2.D0,   0.D0] ,[3,3])

    ELSEIF(symmetry.EQ.3) THEN !Body Centered Cubic
      a_cell = reshape([ a/2.D0, a/2.D0, a/2.D0, &
                        -a/2.D0, a/2.D0, a/2.D0, &
                        -a/2.D0,-a/2.D0, a/2.D0 ], [3,3])

    ELSEIF(symmetry.EQ.4) THEN !Hexagonal
      IF(c.EQ.0.D0) print*,'C IS NULL'
      a_cell = reshape([    a,              0.0D0, 0.0D0,&
                     -a/2.0D0,a*sqrt(3.0d0)/2.0D0, 0.0D0,&
                        0.0D0,              0.0D0, c], [3,3])

    ELSEIF(symmetry.EQ.5) THEN !Trigonal or rhombohedral
    a_cell = reshape([a*SQRT(1.D0-(COS((a_theta(1))*pi/180.0d0)))/SQRT(2.D0),&
                     -a*SQRT(1.D0-(COS((a_theta(1))*pi/180.0d0)))/SQRT(2.D0)/SQRT(3.0d0),&
                     (a*SQRT(1.d0+2.d0*(COS((a_theta(1))*pi/180.0d0)))/SQRT(3.0d0)), &
                0.d0,   SQRT(2.d0)*a*SQRT(1.D0-(COS((a_theta(1))*pi/180.0d0)))/SQRT(3.0d0),&
                     (a*SQRT(1.d0+2.d0*(COS((a_theta(1))*pi/180.0d0)))/SQRT(3.0d0)), &
                     -a*SQRT(1.D0-(COS((a_theta(1))*pi/180.0d0)))/SQRT(2.D0),&
                     -a*SQRT(1.D0-(COS((a_theta(1))*pi/180.0d0)))/SQRT(2.D0)/SQRT(3.0d0),&
                     (a*SQRT(1.d0+2.d0*(COS((a_theta(1))*pi/180.0d0)))/SQRT(3.0d0))],[3,3])
    ELSEIF(symmetry.EQ.6) THEN !Tetragonal
      IF(c.EQ.0.D0) PRINT*,'C IS NULL'
       a_cell = reshape([    a, 0.D0, 0.D0, &
                          0.D0,    a, 0.D0, &
                          0.D0, 0.D0,    c], [3,3])

    ELSEIF(symmetry.EQ.7) THEN !Body Centred Tetragonal
    a_cell = reshape([ a/2.D0,-a/2.D0, 1/2.D0*c, &
                       a/2.D0, a/2.D0, 1/2.D0*c, &
                      -a/2.D0,-a/2.D0, 1/2.D0*c ], [3,3])

    ELSEIF(symmetry.EQ.8) THEN !Orthorhombic
       IF(b.EQ.0.D0)PRINT*,'B IS NULL'
       IF(c.EQ.0.D0)PRINT*,'C IS NULL'
     a_cell = reshape([   a, 0.D0, 0.D0, &
                          0.D0,    b, 0.D0, &
                          0.D0, 0.D0,    c], [3,3])

    ELSEIF(symmetry.EQ.12) THEN !Monoclinic
       IF(b.EQ.0.D0)PRINT*,'B IS NULL'
       IF(c.EQ.0.D0)PRINT*,'C IS NULL'
       IF(a_theta(1).EQ.0.D0) PRINT*,'ANGLE IS NULL'
      
       a_cell = reshape([ a,0.0D0, 0.0d0, &
                       b*COS((a_theta(1))*pi/180.0d0), b*SQRT(1.D0-(COS((a_theta(1))*pi/180.0d0))**2), 0.0D0, &
                      0.0D0,0.0d0, c ], [3,3])
    ELSEIF(symmetry.EQ.14) THEN !Triclinic
      IF (b.EQ.0.D0)print*,'THE SECOND ARGUMENT IS NULL'
      IF (c.EQ.0.D0)print*,'THE THIRD ARGUMENT IS NULL'
      IF (a_theta(1).EQ.1.D0)print*,'THE FIRST ANGLE IS NULL! '
      IF (a_theta(2).EQ.1.D0)print*,'THE SECOND ANGLE IS NULL! '
      IF (a_theta(3).EQ.1.D0)print*,'THE THIRD ANGLE IS NULL! '

       IF (ACOS(a_theta(1))+ACOS(a_theta(2)).LT.ACOS(a_theta(3)))print*,'LATGEN','ALPHA + BETA >= GAMMA'
       IF (ACOS(a_theta(2))+ACOS(a_theta(3)).LT.ACOS(a_theta(1)))print*,'LATGEN','BETA + GAMMA >= ALPHA'
       IF (ACOS(a_theta(3))+ACOS(a_theta(1)).LT.ACOS(a_theta(2)))print*,'LATGEN','GAMMA + ALPHA >= BETA'

       
    a_cell = reshape([ a,0.D0, 0.D0, &
                       b*COS((a_theta(3))*pi/180.0d0), b*(SQRT(1.D0-(COS((a_theta(3))*pi/180.0d0))**2)), 0.D0, &
    c*(COS((a_theta(2))*pi/180.0d0)),  &
    c*((COS((a_theta(1))*pi/180.0d0))-(COS((a_theta(2))*pi/180.0d0))*(COS((a_theta(3))*pi/180.0d0))) &
       /SQRT(1.D0-(COS((a_theta(3))*pi/180.0d0))**2), &
    c*(SQRT((1.D0+2.D0*(COS((a_theta(1))*pi/180.0d0))*(COS((a_theta(2))*pi/180.0d0))*(COS((a_theta(3))*pi/180.0d0)) &
            -(COS((a_theta(1))*pi/180.0d0))**2-(COS((a_theta(2))*pi/180.0d0))**2 &
            -(COS((a_theta(3))*pi/180.0d0))**2)/(1.D0-(COS((a_theta(3))*pi/180.0d0))**2))) &
                     ] , [3,3])
                      
    ELSEIF(symmetry.EQ.12.OR.symmetry.EQ.14.OR.symmetry.EQ.9.OR. symmetry.EQ.10.OR.&
            symmetry.EQ.11.OR.symmetry.EQ.13) THEN
      IF(IONODE) WRITE(*,'(A,I3,A)') ' BRAVAIS LATTICE',symmetry,' NOT PROGRAMMED'
    STOP
    ENDIF
    RETURN
  
  END SUBROUTINE

  FUNCTION get_cell_volume(a_cell) RESULT(vol)
    !
    ! Calculating volume (Omega) of the cell by
    ! Omega= det(a), where a=[a1,a2,a3]
    !
    USE kinds
    USE math, ONLY : determinant
    IMPLICIT NONE

    REAL(KIND=dp), POINTER :: a_cell(:,:)
    REAL(KIND=dp) :: vol
    
!    vol=determinant(a_cell)
        vol   = a_cell(1,1)*(a_cell(2,2)*a_cell(3,3)-a_cell(2,3)*a_cell(3,2)) +&
            a_cell(1,2)*(a_cell(2,3)*a_cell(3,1)-a_cell(2,1)*a_cell(3,3)) + &
            a_cell(1,3)*(a_cell(2,1)*a_cell(3,2)-a_cell(2,2)*a_cell(3,1))

    vol=ABS(vol)
  END FUNCTION get_cell_volume

  SUBROUTINE set_b_cell(a_lattice,a_cell,b_cell)
    !
    ! Calculating b matrix from a
    ! Here a = [a1,a2,a3]
    ! here b = [b1,b2,b3] = 2pi (a^T)^-1
    ! Internally, b is stored in the units of 2pi/|a1_vector|. Thus we don't multiply by 2pi, and we divide by a_lattice
    !
    USE kinds
    USE math,  ONLY : inverse_matrix
    USE constants
    IMPLICIT NONE

    REAL(KIND=dp)          :: a_lattice
    REAL(kind=dp), POINTER :: a_cell(:,:), b_cell(:,:)

    REAL(KIND=dp), POINTER :: work(:,:)
    ALLOCATE(work(3,3))
    work = TRANSPOSE(a_cell)

    CALL inverse_matrix(work,b_cell)
    b_cell(:,:)=b_cell(:,:)*a_lattice

    DEALLOCATE(work)
  END SUBROUTINE set_b_cell


  SUBROUTINE set_r_grids(a_cell,ecut_rho,nrgrids)
    !Using Sampling Theorem to find optimal real space grids
    !Finally, the number will be adjusted for an optimal FFT
    !Sampling Theorem: Length(L) / number_of_points (N)
    !Nyquist critical frequency f_c = 1/(2Delta)
    !For a given plane-wave cutoff (frequency), there is a minimum number of
    !equidistant real space grid points needed for the same accuracy
    !  N = (L/Pi)*\sqrt{E_cut}
    USE kinds
    USE constants
    USE fft_support
    IMPLICIT NONE

    REAL(KIND=dp), DIMENSION(:,:), POINTER :: a_cell
    REAL(KIND=dp) :: ecut_rho
    INTEGER, DIMENSION(:), POINTER :: nrgrids

    REAL(KIND=dp) :: ax,ay,az
    INTEGER :: nx,ny,nz

    !Length along a direction
    ax=SQRT(DOT_PRODUCT(a_cell(1:3,1),a_cell(1:3,1)))
    ay=SQRT(DOT_PRODUCT(a_cell(1:3,2),a_cell(1:3,2)))
    az=SQRT(DOT_PRODUCT(a_cell(1:3,3),a_cell(1:3,3)))

    !N=(L/pi)*E_cut^{1/2)  (Sampling Theorm)
    nx=NINT(ax/pi*SQRT(ecut_rho)+.5D0)
    ny=NINT(ay/pi*SQRT(ecut_rho)+.5D0)
    nz=NINT(az/pi*SQRT(ecut_rho)+.5D0)


    nrgrids(1)=optimal_fft_order(nx,2) 
    nrgrids(2)=optimal_fft_order(ny,2)!optimal_fft_grid_fft(ny) 
    nrgrids(3)=optimal_fft_order(nz,2)
    
  END SUBROUTINE set_r_grids
END MODULE set_cell



