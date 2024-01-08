MODULE yamlread
  CONTAINS
  SUBROUTINE reading_input
  USE constants
  USE, INTRINSIC :: iso_fortran_env, ONLY:  output_unit
  USE yaml, ONLY: parse, error_length
  USE yaml_types, ONLY: type_node, type_dictionary, type_error,  &
                        type_list, type_list_item, type_scalar
  USE system_data_types 
 
  CLASS(type_node), POINTER :: root
  CHARACTER(LEN=error_length) :: error 
  CLASS(type_dictionary), POINTER :: dict,dict1
  CLASS (type_list), POINTER :: list
  CLASS (type_list_item), POINTER :: item
  TYPE (type_error), POINTER :: io_err
    
  CHARACTER(LEN=:), ALLOCATABLE :: string,sms 
  INTEGER::is,i
  CHARACTER(LEN=80)::id
  CHARACTER(LEN=2), ALLOCATABLE :: LABEL(:)  
  REAL(KIND=DP):: M_FACT

  atom_t=0
  root => parse(trim(PATHOFINPUT)//input_filename,error = error)
  IF (error/='') then
    IF(IONODE)print*,trim(error)   
  END IF
  
  SELECT TYPE (root)
  CLASS is (type_dictionary)
  
       task_option = root%get_string('TASK',error=io_err)                          !TASK
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
       IF (task_option.NE.'GEOMETRY_OPTIMIZATION' .AND. task_option.NE.'WAVEFUNCTION_OPTIMIZATION' &
        .AND. task_option.NE.'MOLECULAR_DYNAMICS')THEN
              IF(IONODE)write(6,*)"Invaild TASK !! Choose 'GEOMETRY_OPTIMIZATION' or 'WAVEFUNCTION_OPTIMIZATION' &
	      or 'MOLECULAR_DYNAMICS'" 
       !IF (task_option.NE.'GEOMETRY_OPTIMIZATION' .AND. task_option.NE.'WAVEFUNCTION_OPTIMIZATION')THEN
       !IF(IONODE)write(6,*)"Invaild TASK !! Choose 'GEOMETRY_OPTIMIZATION' or 'WAVEFUNCTION_OPTIMIZATION'" 
              STOP
       ENDIF      
       !#dict => root%get_dictionary('UNITS',required=.TRUE.,error=io_err)           !UNITS
       !#IF (associated(io_err)) THEN 
       !  PRINT* ,trim(io_err%message)
       !  units%en="RYDBERG"
       !  units%coordi="BOHR"
       !ELSE 
       !  units%en = trim(dict%get_string("ENERGY","RYDBERG",error = io_err))
       !  IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
       !  units%coordi = trim(dict%get_string("COORDINATES","BOHR",error = io_err))
       !  IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
       !END IF 


       dict => root%get_dictionary('INITIALIZATION',required=.TRUE.,error=io_err)    !INITIALIZATION
       IF (associated(io_err)) THEN 
         IF(IONODE)PRINT* ,trim(io_err%message)
         init%wfunc="ATOMIC"
         init%hessian="UNIT"
         Lopen_shell=.FALSE.
         multiplicity=1
       ELSE 
         init%wfunc = trim(dict%get_string("WAVEFUNCTION","ATOMIC",error = io_err))
         IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
         !init%coordi = trim(dict%get_string("COORDINATES",error = io_err))
         !IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
         init%hessian = trim(dict%get_string("HESSIAN","UNIT",error = io_err))
         IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
         Lopen_shell = (dict%get_logical("OPEN_SHELL",.FALSE.,error = io_err))
         IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
         multiplicity = (dict%get_integer("MULTIPLICITY",1,error = io_err))
         IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
       ENDIF
      
!############################################ RITAMA ####################################### 

         IF(TASK_option=="MOLECULAR_DYNAMICS") THEN
           dict => root%get_dictionary('DYNAMICS',required=.TRUE.,error=io_err)       !MOLECULAR DYNAMICS
           IF (associated(io_err)) THEN 
                IF(IONODE)PRINT*,(trim(io_err%message))
                md%temperature=300.d0  ! in kelvin
                md%max_step=999
                md%time_step=1.0d0     ! in fs
           ELSE
                IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
                md%temperature = (dict%get_real("TEMPERATURE",error = io_err))
                IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
                md%max_step =  (dict%get_integer("MAX_STEPS",error = io_err))
                IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
                md%time_step = (dict%get_real("TIMESTEP",error = io_err))
                IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
           ENDIF
         ENDIF  
       
!############################################ RITAMA ####################################### 
       dict => root%get_dictionary('OPTIMIZATION',required=.TRUE.,error=io_err)    !WAVEFUNCTION_OPTIMIZATION
       IF (associated(io_err)) THEN 
         IF(IONODE)PRINT*,(trim(io_err%message))
         wf_opt%minimizer="CONJUGATE_GRADIENT"
         wf_opt%max_step=800
         wf_opt%g_cutoff=0.0000001d0
         IF(TASK_option=="GEOMETRY_OPTIMIZATION") THEN
         geo_opt%minimizer="LBFGS"
         geo_opt%max_step=100
         geo_opt%g_cutoff=0.0001d0
         ENDIF
       ELSE
         dict1 => dict%get_dictionary('WAVEFUNCTION',required=.TRUE.,error=io_err)    
         IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
         wf_opt%minimizer = trim(dict1%get_string("MINIMIZER","CONJUGATE_GRADIENT",error = io_err))
         IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
         wf_opt%max_step = (dict1%get_integer("MAX_STEPS",800,error = io_err))
         IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
         wf_opt%g_cutoff = (dict1%get_real("GRAD_CUTOFF",0.0000001d0,error = io_err))
         IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
         IF(TASK_option=="GEOMETRY_OPTIMIZATION") THEN
         dict1 => dict%get_dictionary('GEOMETRY',required=.TRUE.,error=io_err)       !GEOMETRY OPTIMIZATION
         IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
         geo_opt%minimizer = trim(dict1%get_string("MINIMIZER","BFGS",error = io_err))
         IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
         geo_opt%max_step = (dict1%get_integer("MAX_STEPS",100,error = io_err))
         IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
         geo_opt%g_cutoff = (dict1%get_real("GRAD_CUTOFF",0.001d0,error = io_err))
         IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
         ENDIF
       ENDIF

       symmetry_str = root%get_string('SYMMETRY',"CUBIC",error=io_err)                     !SYMMETRY
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))

       dict => root%get_dictionary('CUTOFF',required=.TRUE.,error=io_err)  !WAVEFUNCTION CUTOFF
       IF (associated(io_err)) THEN
         IF(IONODE) PRINT*,trim(io_err%message)
         cutoff%pw= 80
         units%en="RYDBERG"
       ELSE
         cutoff%pw = (dict%get_real("VALUE",80.0d0,error = io_err))
         IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
         units%en = trim(dict%get_string("UNIT","RYDBERG",error = io_err))
         IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
       END IF

      
       dict => root%get_dictionary('CELL',required=.TRUE.,error=io_err)
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
       cell_unit = trim(dict%get_string("UNIT",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
       list => dict%get_list('DIMENSION',required=.TRUE.,error=io_err)
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))

 
    item => list%first
    i=1
    DO WHILE(associated(item))
      SELECT TYPE (element => item%node)
      CLASS is (type_scalar)
        string =(element%string)
        READ(string,*)primitive_vec(i)
        item => item%next
      END SELECT
      i=i+1
    ENDDO
     !-------------------------------------------------------------
     list => dict%get_list('ANGLE',required=.TRUE.,error=io_err)
     IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))


      item => list%first
    i=1
    DO WHILE(associated(item))
      SELECT TYPE (element => item%node)
      CLASS is (type_scalar)
        string =(element%string)
        READ(string,*)a_theta(i)
        item => item%next
      END SELECT
      i=i+1
    ENDDO

  select case (trim(adjustl(cell_unit)))
  case ("ANGSTROM")
    primitive_vec = primitive_vec*angs2bohr
  case ("BOHR")
    primitive_vec = primitive_vec*1.0_dp
  case default
     if(ionode)write(*,*) "Invalid unit. See amdkiit_manual"
     stop
  end select


     !-------------------------------------------------------------
    dict => root%get_dictionary('INPUT',required=.TRUE.,error=io_err)
    IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
      
     coordinate_file = trim(dict%get_string("COORDINATES_FILE",error = io_err))
     IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
     
     coordinate_unit = trim(dict%get_string("UNIT",error =io_err))
     IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))

  select case (trim(adjustl(coordinate_unit)))
  case ("ANGSTROM")
     m_fact = angs2bohr 
  case ("BOHR")
     m_fact = 1.0_dp
  case default
     if(ionode)write(*,*) "Invalid unit. See amdkiit_manual"
     stop
  end select
    


     !sp_t = root%get_integer('NUMBER_OF_SPECIES',error=io_err)                          !TASK
     !IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
   call readcoordinate(coordinate_file,LABEL,m_fact)
!---------------------------------------------------------------------------------------
    ALLOCATE(PPFILE(sp_t),TKB(sp_t),LMAX(sp_t),SKIP(sp_t),l_local(sp_t))!ATOM_MASS(sp_t))!atom_p_sp(sp_t),ATOM_MASS(sp_t))
!-------------------------------------------------------------------------------------- 
  DO IS=1,sp_t
    write(id,*)is    
     dict1 => root%get_dictionary('SPECIES',required=.TRUE.,error=io_err)      !ATOMS
      
      dict => dict1%get_dictionary('LABEL'//trim(adjustl(id)),required=.TRUE.,error=io_err)      !ATOMS
      IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
      

       string = trim(dict%get_string("ATOM",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
       if(string.ne.label(is)) then
        if(ionode)print*,"Species name different !!!!"
        STOP
       endif       
       ppfile(is)= trim(dict%get_string("PP_FILE",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
      
       string = trim(dict%get_string("PP_SCHEME",'NONE',error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))     
           IF(INDEX(string,'KLEINMAN').NE.0)THEN
           tKB(Is)=.TRUE.
           ELSE
           tKB(Is)=.FALSE.
           ENDIF

        STRING =trim(dict%get_string("LMAX",error = io_err))
        IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
      
        LMAX(IS)=LVAL(STRING)   

        !STRING = (dict%get_STRING("LOCAL",error = io_err))
        !IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
     
        L_LOCAL(IS)=LVAL(STRING)
 
        !atom_mass(is) = (dict%get_real("MASS",error = io_err))
        !IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
      
        !atom_p_sp(is) = (dict%get_integer("NUMBER",error = io_err))
        !IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
        !atom_t=atom_t+atom_p_sp(is)
  END DO


    dict => root%get_dictionary('OUTPUT',required=.TRUE.,error=io_err)         !OUTPUT
           IF (associated(io_err)) THEN
             IF(IONODE)PRINT*,(trim(io_err%message))
             debug=.FALSE.
             print_density=.FALSE.
           ELSE
            debug = (dict%get_logical("DEBUG",.FALSE.,error = io_err))
            IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))

            print_density = (dict%get_logical("PRINT_DENSITY",.FALSE.,error = io_err))
            IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
          endif

       !print_density = trim(dict%get_string("PRINT_DENSITY",.FALSE.,error = io_err))
       !IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
       
       !frequency = (dict%get_real("FREQUENCY",error = io_err))
       !IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
      
    dict => root%get_dictionary('DFT',required=.TRUE.,error=io_err)                 !!DFT
    IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
      
       func_string(1)= trim(dict%get_string("EXCHANGE",'NONE',error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
      
       func_string(2) = trim(dict%get_string("CORRELATION",'NONE',error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
     
       func_string(3) = trim(dict%get_string("XC",'NONE',error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
     
  END SELECT
  CALL root%finalize()
  DEALLOCATE(root)
  

  DO Is=1,sp_t
    SKIP(Is)=-1
    IF(TKB(IS).AND.L_LOCAL(IS).GE.0.AND.L_LOCAL(IS).LE.LMAX(IS))THEN
    L_LOCAL(IS)=L_LOCAL(IS)+1
    ELSE
    L_LOCAL(IS)=LMAX(IS)+1
    ENDIF

    IF(TKB(IS).AND.SKIP(IS).GE.0.AND.SKIP(IS).LE.LMAX(IS))THEN
    SKIP(IS)=SKIP(IS)+1
    ELSE
    SKIP(IS)=LMAX(IS)+2
    ENDIF
    LMAX(IS)=LMAX(IS)+1
  ENDDO



END SUBROUTINE reading_input



SUBROUTINE CHECK_ERROR(msg)
  USE yaml
  USE yaml_types
  USE system_data_types
IMPLICIT NONE
CHARACTER(LEN=*)              :: msg
IF(IONODE)PRINT* ,msg
STOP 1
END SUBROUTINE



FUNCTION LVAL(L)
 IMPLICIT NONE
 CHARACTER,INTENT(IN)::L
 INTEGER  :: LVAL

 IF(L.EQ.'S')THEN
 LVAL=0
 ELSE IF(L.EQ.'P')THEN
 LVAL=1
 ELSE IF(L.EQ.'D')THEN
 LVAL=2
 ELSE IF(L.EQ.'F')THEN
 LVAL=3
 ELSE IF(L.EQ.'G')THEN
 LVAL=4
 ELSE IF(L.EQ.'H')THEN
 LVAL=5
 ELSE
 LVAL=-1
 ENDIF
END FUNCTION



   SUBROUTINE readcoordinate(filename,LABEL,fact_m)
   USE system_data_types

   IMPLICIT NONE
   CHARACTER(LEN=*),intent(IN)               :: filename
   CHARACTER(LEN=*),ALLOCATABLE, intent(OUT)              :: LABEL(:)
   REAL(KIND=DP),INTENT(IN)::fact_m
   INTEGER                                   :: J,I,IA,numLines,iunit
   CHARACTER LINE*80, PATH*42

    character(len=2) :: current_element, previous_element
    integer :: num_elements,count
    logical :: first_element
    character(len=2), allocatable :: elements(:)
    integer, allocatable :: scr_counts(:)

   ALLOCATE(scr_counts(10))
   !ALLOCATE(LABEL(SP_T))
    previous_element="  "
    atom_t = 0
    sp_t=0
    scr_counts=1

    OPEN(UNIT=11,FILE=trim(PATHOFINPUT)//filename, ACTION = "READ" ,STATUS='OLD' )
    read(11,*,iostat=i) 
    read(11,*,iostat=i)
 
    DO
        read(11,*,iostat=i) current_element
        if (i /= 0) exit ! Exit loop when end of file is reached
        if(current_element.ne.previous_element)then
        sp_t=sp_t+1
        else 
        scr_counts(sp_t)= scr_counts(sp_t) +1
        endif
        previous_element=current_element
        atom_t = atom_t + 1
    END DO
   CLOSE(11)
   
    ALLOCATE(atom_p_sp(sp_t))
    DO I=1,sp_t
    atom_p_sp(I)=scr_counts(I)  
    ENDDO

    DEALLOCATE(scr_counts)

   DO I=1,sp_t
    NA_MAX=MAX(atom_p_sp(I),NA_MAX)
   ENDDO

  
   ALLOCATE(ATCO(3,NA_MAX,sp_t))
   ALLOCATE(LABEL(SP_T))

   OPEN(UNIT=12,FILE=trim(PATHOFINPUT)//filename, ACTION = "READ" , STATUS='OLD' )
   READ(12,*)
   READ(12,*)
   DO J=1,sp_t
     DO IA=1,atom_p_sp(J)
     READ(12,*)LABEL(J),(ATCO(I,IA,J),I=1,3)
     !K=K+1
     ENDDO
   ENDDO
   ATCO=ATCO*Fact_m
   CLOSE(12)
 END SUBROUTINE readcoordinate
 
 SUBROUTINE read_and_print
 USE system_data_types
 IMPLICIT NONE 
 INTEGER IOS
 CHARACTER(LEN=100) :: line ! declare a character variable to store each line of text
   write(6,"(A)") repeat("*", 93)
   write(6, "(3A)")repeat(" ", 40),"INPUT DETAILS",repeat(" ", 40)
   write(6,"(A)") repeat("*", 93)
   open(10, file=trim(PATHOFINPUT)//input_filename, status='old', action='read') ! open the input file for reading
   do 
     read(10, '(A)', iostat=ios) line ! read a line of text from the file and check the status
     if (ios < 0) exit ! if end of file is reached, exit the loop
     if (ios > 0) stop 'Error reading file' ! if an error occurs, stop the program
     write(*, '(A)') line ! write the line of text to the standard output
   end do 
   close(unit=10) ! close the input file
   write(6,"(A)") repeat("*", 93)
 END SUBROUTINE read_and_print
END MODULE


