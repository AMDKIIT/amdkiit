MODULE init_module
CONTAINS
SUBROUTINE initialize
  USE kinds
  USE mympi
  USE max_parameter_pp
  USE constants
  USE system_data_types
  USE yamlread
  USE set_cell
  USE gvectors
  USE fft_interface, ONLY: prepare_fft
  USE exp_igr
  USE atomic_basis
  USE potential
  USE density
  USE xc
  USE wfn_initialize
  USE pseudopotential

  IMPLICIT NONE
  INTEGER ::i,is,ia
  IF(ionode)WRITE(6,'(A,I10)')'Number of MPI processes = ', ncpu

  ! READ INPUTS
  PATHOFINPUT='./'
  if(IONODE)CALL read_and_print
  CALL reading_input
  CALL atom_info
  !CALL readcoordinate(coordinate_file)
  
  nlp=0!TODO
  gamma_point=.TRUE.
  NLSD = 1 !dfault
  is_lsd_but_not=.FALSE.

! Cutoffs
  cutoff%rho  = cutoff%pw*4.d0

! Setting Lattice Vectors
! a_cell is consists of three columns; a1_vector, a2_vector, a3_vector as columns

  ALLOCATE(a_cell(3,3))
  CALL set_vec(primitive_vec(1),primitive_vec(2),primitive_vec(3),a_cell,symmetry_str)
! Calculate Cell Volume
  cell_volume=get_cell_volume(a_cell)

! a_1 vector mod is set as this constant
  a_lattice =  primitive_vec(1)!DSQRT(a_cell(1,1)**2+a_cell(2,1)**2+a_cell(3,1)**2)
  twopibya = 2.0_dp*pi/a_lattice
  twopibya2 = twopibya*twopibya

! G-cutoffs
  Gcutoff%pw  = cutoff%pw/twopibya2
  Gcutoff%rho = cutoff%rho/twopibya2

! Calculate Reciprocal Lattice vectors
  ALLOCATE(b_cell(3,3))
  CALL set_b_cell(a_lattice,a_cell,b_cell) !b_cell in 2pi/|a1} units
  IF(IONODE)write(6, "(3A)")repeat(" ",35),"SIMULATION PARAMETERS",repeat(" ", 27)
  IF(IONODE)write(6,"(A)")repeat("*", 93)

  IF(IONODE)WRITE(*,"(A40,F16.6)") "Wavefunction Cutoff (Ry) = ", cutoff%pw
  IF(IONODE)WRITE(*,"(A40,F16.6)") "Density Cutoff (Ry)     = ", cutoff%rho
  IF(IONODE)WRITE(*,"(A40,3A)") "Exchange Corelation funational"
  IF(IONODE)WRITE(*,"(A40,A)")"Exchange : ",func_string(1),"Corelation : ",func_string(2), "Exchange-Corelation : ",func_string(3)
  IF(IONODE)write(6,"(A)")repeat("*", 93)
  IF(IONODE)write(6, "(3A)")repeat(" ",42),"SUPERCELL",repeat(" ", 42)
  IF(IONODE)write(6,"(A)")repeat("*", 93) 
  IF(ionode)THEN
    WRITE(*,"(2A)") "Lattice Type = ",symmetry_str
    WRITE(*,"(A,F16.6)") "Lattice Constant (Bohr) = ", a_lattice
    WRITE(*,"(A,F16.6)") "Cell Volume (Bohr^3)    = ", cell_volume
    write(*,"(A,F16.6)") "Reciprocal Lattice Vector (2pi Bohr^-1)"
    write(*,"(A,3F16.6)") 'B1 =', b_cell(1:3,1)/a_lattice
    write(*,"(A,3F16.6)") 'B2 =', b_cell(1:3,2)/a_lattice
    write(*,"(A,3F16.6)") 'B3 =', b_cell(1:3,3)/a_lattice
  END IF
  IF(IONODE)write(6,"(A)")repeat("*", 93)
  IF(IONODE)write(6, "(3A)")repeat(" ",42),"SPECIES",repeat(" ", 42)
  IF(IONODE)write(6,"(A)")repeat("*", 93)
  IF(IONODE)write(6, "(2A)")repeat(" ",30),coordinate_unit
  IF(IONODE)write(6,"(A5,3A16,A20)")"TYPE","X    ","Y    ","Z    ","PSEUDOPOTENTIAL"
  IF(IONODE)write(6,"(A)")repeat("-", 93)
!  DO IS=1,sp_t
!   DO IA=1,atom_p_sp(IS)
!     IF(IONODE)write(6,"(A5,3F16.7,A33)")symbol(z(IS)),(ATCO(I,IA,IS),I=1,3),PPFILE(IS)
!   ENDDO
!  ENDDO 
!  IF(IONODE)write(6,"(A)")repeat("*", 93)
  
  ! Set real space grids (tuned to FFT)
  ALLOCATE(nrgrids(3))
  CALL set_r_grids(a_cell,cutoff%rho,nrgrids)
  CALL read_pp_file
    DO IS=1,sp_t
   DO IA=1,atom_p_sp(IS)
     atmass(z(is))=atom_mass(z(is))* amu_au
     IF(IONODE)write(6,"(A5,3F16.7,10x,A33)")symbol(z(IS)),(ATCO(I,IA,IS),I=1,3),PPFILE(IS)  ! PG debug
   ENDDO
  ENDDO
  IF(IONODE)write(6,"(A)")repeat("*", 93)


  CALL set_ngvectors
  CALL set_reciprocal
  CALL prepare_fft

  if(l_upf)then
    CALL form_factor_upf
  else
     CALL form_factor
  endif
  CALL exp_igrxfactor !!    STRUCTURE FACTOR CALCULATION
  CALL setbasis
  IF(IONODE)WRITE(*,"(A40,I5)") "Total Number of Valence Electrons = ", tn_velec
  IF(IONODE)WRITE(*,"(A40,I5)") "Total Number of Valence Orbital = ",nstate
  !IF(IONODE)WRITE(*,"(A40,13F5.1)") "Occupation = ", (occupation(I),I=1,nstate)
  IF(IONODE)WRITE(*,"(A40)") "Occupation = "
  IF(IONODE)WRITE(*,"(13F5.1)")(occupation(I),I=1,nstate)

  IF(LOPEN_SHELL)THEN
    NSPIN = (MULTIPLICITY-1) ! no of unpaired electron
    NEL_UP   = INT((NSTATE-NSPIN)/2.0d0) ! no of alpha electron
    NEL_DOWN = INT((NSTATE+NSPIN)/2.0d0) ! no of alpha electron
    IF(ionode)WRITE(*,*)"No. of alpha electrons =", NEL_UP
    IF(ionode)WRITE(*,*)"No. of beta electrons =", NEL_DOWN
    IF(NSTATE.NE.NEL_UP+NEL_DOWN.OR.NEL_UP.LT.0.OR.NEL_DOWN.LT.0)THEN
        WRITE(*,*)"Incorrect multiplicity"
        STOP
    ENDIF
  ENDIF

  call distribute_orb

  call distribute_atom

  CALL eval_pot

  CALL rhog2r
  CALL cal_vpot_ex
  if(l_upf) then
      DO i=1,sp_t
      CALL upf_nlppgrid(i)
      ENDDO
      CALL upf_nlpro
  else
      DO i=1,sp_t
      CALL nlppgrid(i)
      ENDDO
      CALL nlpro
  endif

  CALL wfn_initialize_orb
    IF(LOPEN_SHELL)THEN
    is_lsd_but_not=.FALSE.
    NLSD = 2
    DEALLOCATE(RHO,Rho_G,V_POT)
    ALLOCATE(RHO(NNR1,NLSD),Rho_G(NGRHO_l,NLSD),V_POT(NNR1,NLSD))
    ENDIF
END SUBROUTINE initialize
END MODULE init_module
