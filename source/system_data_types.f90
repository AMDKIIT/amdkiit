  MODULE system_data_types
  USE kinds
  USE max_parameter_pp
  INTEGER       :: ncpu, & !Number of CPUs
                   icpu    ! id of each CPU
   
  CHARACTER (LEN=80):: task_option

  TYPE :: unit_data
   CHARACTER (LEN=80) :: en      !< Energy unit
   CHARACTER (LEN=80) :: coordi  !< Co-ordinate unit
  END TYPE unit_data

  TYPE :: initialization
   CHARACTER (LEN=80) :: wfunc      !< TODO 
   CHARACTER (LEN=80) :: coordi     !<
   CHARACTER (LEN=80) :: hessian    !<
  END TYPE initialization

  TYPE :: dynamics
   INTEGER:: max_step
   REAL(KIND=DP) :: temperature, time_step
  END TYPE dynamics

  TYPE :: optimization
   CHARACTER (LEN=80) ::minimizer
   INTEGER::max_step
   REAL(KIND=DP) ::g_cutoff,e_cutoff
  END TYPE  optimization

  TYPE:: cut_off
   REAL(KIND=DP) :: pw  !> Wavefunction cutoff
   REAL(KIND=DP) :: rho !> Density cutoff
  END TYPE cut_off
  
  REAL(KIND=DP) :: order,extrapolate,number,frequency
  CHARACTER (LEN=80)::cell_unit
  REAL(KIND=DP), DIMENSION(3)::primitive_vec,a_theta
  
  CHARACTER (LEN=80):: coordinate_file,coordinate_unit
  !CHARACTER,DIMENSION(atommax)                :: LABEL
  INTEGER:: sp_t,tn_velec
  CHARACTER (LEN=80), DIMENSION(:), POINTER :: ppfile
  CHARACTER (LEN=80)::input_filename
  LOGICAL, DIMENSION(:), POINTER :: tkb   !>True if Kleinman-Bylander separable form
  INTEGER, DIMENSION(:), POINTER :: skip,l_local,atom_p_sp,lmax
  !REAL(KIND=DP), DIMENSION(:), POINTER ::atom_mass
  

  TYPE(unit_data):: units
  TYPE(initialization):: init
  TYPE(dynamics):: md !RITAMA
  TYPE(optimization):: wf_opt
  TYPE(optimization):: geo_opt
  TYPE(cut_off)::cutoff  !> Real space
  TYPE(cut_off)::Gcutoff !> Reciprocal space
  

  CHARACTER     :: PATHOFINPUT*42,PATHOFOUTPUT*42
  !REAL(KIND=dp) :: Gcutoff%pw, Gcutoff%rho
  REAL(KIND=dp) :: a_cell_scale       !Parameter to scale the cell 

  REAL(KIND=dp) :: cell_volume, &        !Cell Volume
                   a_lattice, &
                   twopibya, &
                   twopibya2
  INTEGER       :: symmetry
  CHARACTER (LEN=30)::symmetry_str 
  REAL(KIND=dp) :: rr(splnmax,atommax)
  INTEGER       :: atom_t,na_max
  REAL(KIND=dp) :: rw(splnmax,spmax)
  REAL(KIND=dp), PARAMETER :: raggio=1.2d0
  INTEGER,       DIMENSION(:,:),   POINTER ::iatpt,NST12
  INTEGER,       DIMENSION(:),   POINTER ::iatpe
  REAL(KIND=dp) :: etxc
  COMPLEX*16    :: eloc,eself,eh,ee,eion,elpp
  REAL(KIND=dp) :: exc,estat,enl,etotal,esr
  REAL(KIND=dp) :: wsg(spmax,maxngh)
  INTEGER       ::  nghtol(maxngh,spmax),nghcom(maxngh,spmax)! HOW TO CALCULATE NHX
  INTEGER       :: ngh(spmax),nl_angmom !TOTAL NUMBER OF NON-LOCAL ANGULAR MOMENTUM
  REAL(KIND=dp) :: foc(20),E_OLD
  REAL(KIND=dp) :: wns(nspln,5,spmax)!,wns(nspln,2,5,spmax)!define in pseudopotential file
  INTEGER(8)    :: NNR1
  INTEGER       :: ngpw, &   !Number of planewaves (< Gcutoff%pw)
                   ngrho,&   !Number of planewaves (< Gcutoff%rho)
                   nstate,&  !Number of electronic states
                   !ngroup,&
                   nghmax,&
                   ngpw_l,&  !Number of planewaves(< gcut_wf) on a processor iproc
                   ngrho_l,& !Number of planewaves(< Gcutoff%rho) on a processor iproc
                   !nnr1,&      !Total Number of grids (nx * ny * nz without distributing over processors)
                   g0_ip     !Processor ID holding g0
  LOGICAL       :: gamma_point, &
                   ionode,debug,print_density,&
                   g0_stat,&  !TRUE if the processor is holding the g0 component 
                   tfor
  !REAL(KIND=dp) :: primitive_vec(3)
  REAL(KIND=dp) :: eke      !Kinetic Energy
!######################################### RITAMA ##########################################
  REAL(KIND=dp) :: atke, atpe, atte!, com(3)      !Atomic Kinetic Energy 
  REAL(KIND=dp) :: temperature,temp_inst      !Temperature (in K) 
  REAL(KIND=dp) :: time_step_fs, time_step ! Time Step (in fs) 
  INTEGER       :: maxmd !maximum number of MD step
  INTEGER       :: ndof !Degrees of Freedom
!######################################### RITAMA ##########################################
  INTEGER       :: nlp !Number of lattice point
  REAL(KIND=dp) :: atwt(99)
      DATA atwt /1.00797d0,  4.0026d0,    6.939d0,   9.0122d0,   10.811d0,&
     & 12.01115d0,14.0067d0,  15.9994d0,  18.9984d0,   20.183d0,&
     & 22.9898d0,  24.312d0,  26.9815d0,   28.086d0,  30.9738d0,&
     & 32.064d0,   35.453d0,   39.948d0,   39.102d0,   40.080d0,&
     & 44.956d0,   47.900d0,   50.942d0,   51.996d0,   54.938d0,&
     & 55.847d0,   58.933d0,   58.710d0,   63.540d0,   65.370d0,&
     & 69.720d0,   72.590d0,   74.922d0,   78.960d0,   79.909d0,&
     & 83.800d0,   85.470d0,   87.620d0,   88.905d0,   91.220d0,&
     & 92.906d0,   95.940d0,   98.000d0,  101.070d0,  102.905d0,&
     &106.400d0,  107.870d0,  112.400d0,  114.820d0,  118.690d0,&
     &121.750d0,  127.600d0,  126.904d0,  131.300d0,  132.905d0,&
     &137.340d0,  138.910d0,  140.120d0,  140.907d0,  144.240d0,&
     &147.000d0,  150.350d0,  151.960d0,  157.250d0,  158.924d0,&
     &162.500d0,  164.930d0,  167.260d0,  168.934d0,  173.040d0,&
     &174.970d0,  178.490d0,  180.948d0,  183.850d0,  186.200d0,&
     &190.200d0,  192.200d0,  195.090d0,  196.967d0,  200.590d0,&
     &204.370d0,  207.190d0,  208.980d0,  210.000d0,  210.000d0,&
     &222.000d0, 13*250.0d0/

  REAL(KIND=dp),    DIMENSION(:,:),       POINTER :: a_cell    !Lattice Vector
  REAL(KIND=dp),    DIMENSION(:,:),       POINTER :: b_cell    !Reciprocal Lattice Vector
  INTEGER,          DIMENSION(:),         POINTER :: nrgrids   !Number of real space grids
  INTEGER,          DIMENSION(:),         POINTER :: nrgrids_l !Number of real space grids (in each processor iproc)
  INTEGER,          DIMENSION(:),         POINTER :: nrlead    !Number of real space grids, leading dimension
  INTEGER,          DIMENSION(:),         POINTER :: nrlead_l  !Number of real space grids (in each processor iproc)
  INTEGER,          DIMENSION(:,:),       POINTER :: nrxplane  !realspace x-grid is distributed over processors (min,max indices)
  INTEGER,          DIMENSION(:,:),       POINTER :: inyh      !Stores grid index of a G-vector with G^2<Gcut_rho
  INTEGER,          DIMENSION(:,:),       POINTER :: Ggridinfo !Information on #grids,#rays on each processor
  INTEGER,          DIMENSION(:),         POINTER :: mapgp     !pointer to G vectors
  INTEGER,          DIMENSION(:),         POINTER :: map_grid1d_p !map3d to 1d grid
  INTEGER,          DIMENSION(:),         POINTER :: map_grid1d_m !map3d to 1d grid
  REAL(KIND=dp),    DIMENSION(:),         POINTER :: hg        !G.G values are stored
  REAL(KIND=dp),    DIMENSION(:,:),       POINTER :: gvec      !G-vecors
  REAL(KIND=dp),    DIMENSION(:,:,:),     POINTER :: twnl
  REAL(KIND=dp),    DIMENSION(:,:,:,:),   POINTER :: nl
  REAL(KIND=dp),    DIMENSION(:,:,:,:,:), POINTER :: dfnl
  COMPLEX(KIND=dp), DIMENSION(:,:),       POINTER :: psi,C_0   !Wavefunction coefficients (G space)
  COMPLEX(KIND=dp), DIMENSION(:,:),       POINTER :: C2,eigr,eigr_pw,HNM1
  COMPLEX(KIND=dp), DIMENSION(:),         POINTER :: eigrxvps,eigrxrhos
  COMPLEX(KIND=dp), DIMENSION(:,:),         POINTER :: rho_g
  REAL(KIND=dp),    DIMENSION(:,:),       POINTER :: vps,rhos
  !REAL(KIND=dp),    DIMENSION(:),         POINTER :: v_pot,rho
  REAL(KIND=dp),    DIMENSION(:,:),       POINTER :: v_pot,rho !Sudhir DBG 
  REAL(KIND=dp),    DIMENSION(:,:,:),     POINTER :: force
  INTEGER,          DIMENSION(:,:),       POINTER :: xrays_pw ! count of G-vec x-rays
  INTEGER,          DIMENSION(:),         POINTER :: ngvy       ! count of G-vec along y
  INTEGER,          DIMENSION(:),         POINTER :: ngvz       ! count of G-vec along z
  INTEGER                                         :: gstart_ygrid, gstart_zgrid,&
                                                     gend_ygrid,gend_zgrid,NORBPE
  INTEGER                                         :: max_nrgrid,max_ngrays,max_nhrays,max_nrlead
!FFT related
  INTEGER                                         :: max_fft
  INTEGER,          DIMENSION(:),         POINTER :: xscatter_fft_cat !gather-scatter array !MSP
  INTEGER,          DIMENSION(:),         POINTER :: map_grid1d_p_pw,&
                                                     map_grid1d_m_pw
  COMPLEX*16,       DIMENSION(:),         POINTER :: work1_fft, work2_fft
  REAL(KIND=dp)                                   :: GEMAX,CNORM
!CHARACTER(LEN=20),  DIMENSION(:),         POINTER :: ppfile
  REAL(KIND=dp)                                   :: charge 
  REAL(KIND=dp),    DIMENSION(:,:,:),     POINTER :: atco  !atomic coordinates
  REAL(KIND=dp),    DIMENSION(:),         POINTER :: atcharge  !atomic charge
  REAL(KIND=dp),    DIMENSION(:,:,:),     POINTER :: atvel  !atomic velocities !RITAMA
  !LOGICAL,          DIMENSION(:),         POINTER :: tkb   !True if Kleinman-Bylander separable form
  !INTEGER,          DIMENSION(:),         POINTER :: l_max!,skip,l_local,atom_p_sp
  INTEGER,          DIMENSION(:),         POINTER :: z,zv,xc_fun,meshv,meshw
  REAL(KIND=dp),    DIMENSION(:),         POINTER :: slatr_ex
  REAL(KIND=dp) :: vr(splnmax,spmax,lmaxx),rps(splnmax,spmax,lmaxx) !potential and wavefunction grid
          LOGICAL,  DIMENSION(:),         POINTER :: tnum!  Use of numerical values for PP
!!upf
  !INTEGER,          DIMENSION(:),         POINTER :: z,zv,xc_fun,meshv,meshw

!!!!!!need to modify from here
LOGICAL CONVGEO,CONVWF,PRECONDITION,GOPT
INTEGER NSHELL(2)     !Number of shell n per species          
INTEGER LSHELL(20,2)  !Number of orbital per species 
INTEGER M1SHL,e_config(4,7,99)
CHARACTER(LEN=2) SYMBOL(99)
REAL(dp), DIMENSION(99) :: atom_mass
REAL(dp), DIMENSION(99) :: atmass
REAL(kind=dp), DIMENSION(:),  POINTER ::OCCUPATION !NSTATE !TODO
INTEGER nattot   !Total number of  orbitals
INTEGER, DIMENSION(:),   POINTER:: NUMAOR !Number of  orbitals per species
INTEGER:: NUMAORMAX
 real(kind=dp), dimension(:,:,:,:), ALLOCATABLE::cat
 real(kind=dp), dimension(:,:), ALLOCATABLE::oc
 !real(kind=dp), dimension(:), pointer :: XPAR,DXPAR,YPAR,DYPAR
 real(kind=dp), dimension(:,:), pointer :: HESSIAN
Integer NODIM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
CHARACTER(LEN=20):: species,gminimizr,func_string(3)!L_XC,L_C,L_X
INTEGER :: maxwopt,maxgopt,hess_unit,init_wave,init_coor,wminimizr
REAL(KIND=dp), dimension(:),ALLOCATABLE :: mass
CHARACTER(LEN=20),dimension(:),ALLOCATABLE::ppf,lmx,lcl
REAL(KIND=dp) :: g_cutoff,wfn_convg,worder,geo_convg,gnmax
REAL(KIND=DP)::    gnl(splnmax,spmax,lmaxx)
REAL(KIND=DP)::    upf_gnl(splnmax,spmax,lmaxx)

!----------------------------------------------------------------
! SG TODO 
!REAL(KIND=DP)   RCSG(spmax)!,CLSG(MCFUN,spmax),
!REAL(KIND=DP)      RCNL(LMAXX,spmax)
REAL(KIND=DP)       HLSG(MPRO,MPRO,LMAXX,spmax)
INTEGER       NPRO(LMAXX,spmax),LPVAL(LMAXX*LMAXX*MPRO,spmax)
INTEGER             LFVAL(LMAXX*LMAXX*MPRO,spmax)
INTEGER ::             PPLIST(LMAXX*4,2,spmax)
LOGICAL :: L_UPF

  INTEGER       :: upf_nghtol(maxngh,spmax),upf_nghcom(maxngh,spmax)! HOW TO CALCULATE NHX
  INTEGER       :: upf_ngh(spmax)!,nl_angmom !TO
  INTEGER       :: upf_nghmax
!------------------------------------------------------------------
                    LOGICAL, DIMENSION(:),    POINTER :: upf_l_coulomb
          CHARACTER (len=2), DIMENSION(:),    POINTER :: upf_element
          CHARACTER (len=2), DIMENSION(:),    POINTER :: upf_typ
         CHARACTER (len=20), DIMENSION(:),    POINTER :: upf_rel
        CHARACTER (len=100), DIMENSION(:),    POINTER :: upf_dft
             INTEGER,        DIMENSION(:),    POINTER :: upf_meshv,upf_lmax,upf_local,upf_nwfc,upf_nbeta !z,zv,xc_fun,meshv,meshw
       REAL(kind=dp),        DIMENSION(:,:),  POINTER :: upf_r,upf_rab,upf_vloc
       REAL(kind=dp),        DIMENSION(:),    POINTER :: upf_zp
             INTEGER,        DIMENSION(:,:),  POINTER :: upf_els_beta,upf_lll,upf_kbeta
       REAL(kind=dp),        DIMENSION(:,:),  POINTER :: upf_rcut,upf_rcutus,upf_oc,upf_epseu,upf_rcut_chi,rho_atom
       REAL(kind=dp),        DIMENSION(:,:,:),POINTER :: upf_beta,upf_chi
!Sudhir DBG 
LOGICAL LOPEN_SHELL,is_lsd_but_not
INTEGER MULTIPLICITY, NSPIN, NEL_UP, NEL_DOWN, NLSD
!Sudhir DBG 
REAL(KINd=DP) DTBYM(3*spmax*atommax)
END MODULE system_data_types
