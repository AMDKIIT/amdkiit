MODULE yamlread
  CONTAINS
  SUBROUTINE reading_input
  USE, INTRINSIC :: iso_fortran_env, ONLY:  output_unit
  USE yaml, ONLY: parse, error_length
  USE yaml_types, ONLY: type_node, type_dictionary, type_error, dp, &
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
  
  atom_t=0
  root => parse(trim(PATHOFINPUT)//"input.yaml",error = error)
  IF (error/='') then
    print*,trim(error)   
  END IF
  
  SELECT TYPE (root)
  CLASS is (type_dictionary)
  
    task_option = root%get_string('TASK',error=io_err)                          !TASK
    IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
    
    dict => root%get_dictionary('UNITS',required=.TRUE.,error=io_err)           !UNITS
    IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
       units%en = trim(dict%get_string("ENERGY",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
       units%coordi = trim(dict%get_string("COORDINATES",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
      
    dict => root%get_dictionary('INITIALIZATION',required=.TRUE.,error=io_err)    !INITIALIZATION
    IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
       init%wfunc = trim(dict%get_string("WAVEFUNCTION",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
       init%coordi = trim(dict%get_string("COORDINATES",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
       init%hessian = trim(dict%get_string("HESSIAN",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
          Lopen_shell = (dict%get_logical("OPEN_SHELL",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
        multiplicity = (dict%get_integer("MULTIPLICITY",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))

 
    dict => root%get_dictionary('OPTIMIZATION',required=.TRUE.,error=io_err)    !WAVEFUNCTION_OPTIMIZATION
    IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))

       dict1 => dict%get_dictionary('WAVEFUNCTION',required=.TRUE.,error=io_err)    
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
       wf_opt%minimizer = trim(dict1%get_string("MINIMIZER",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
       wf_opt%max_step = (dict1%get_real("MAX_STEPS",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
       wf_opt%g_cutoff = (dict1%get_real("GRAD_CUTOFF",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
       wf_opt%e_cutoff = (dict1%get_real("ENERGY_CUTOFF",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))

      dict1 => dict%get_dictionary('GEOMETRY',required=.TRUE.,error=io_err)       !GEOMETRY OPTIMIZATION
      IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
       geo_opt%minimizer = trim(dict1%get_string("MINIMIZER",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
       geo_opt%max_step = (dict1%get_real("MAX_STEPS",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
       geo_opt%g_cutoff = (dict1%get_real("GRAD_CUTOFF",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
 
    symmetry_str = root%get_string('SYMMETRY',error=io_err)                     !SYMMETRY
    IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
    !IF(INDEX(symmetry_str,"CUBIC").NE.0)SYMMETRY=1
    dict => root%get_dictionary('WAVEFUNCTION_CUTOFF',required=.TRUE.,error=io_err)  !WAVEFUNCTION CUTOFF
    IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
       cutoff%pw = (dict%get_real("VALUE",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
      
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

    dict => root%get_dictionary('INPUT',required=.TRUE.,error=io_err)
    IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
      
      coordinate_file = trim(dict%get_string("COORDINATES_FILE",error = io_err))
      IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
      
     coordinate_unit = trim(dict%get_string("UNIT",error =io_err))
      IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
      
      sp_t = root%get_integer('NUMBER_OF_SPECIES',error=io_err)                          !TASK
      IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
    
!---------------------------------------------------------------------------------------
    ALLOCATE(PPFILE(sp_t),TKB(sp_t),LMAX(sp_t),SKIP(sp_t),l_local(sp_t),atom_p_sp(sp_t),ATOM_MASS(sp_t))
!-------------------------------------------------------------------------------------- 
   DO IS=1,sp_t
    write(id,*)is    
     dict1 => root%get_dictionary('SPECIES',required=.TRUE.,error=io_err)      !ATOMS
      
      dict => dict1%get_dictionary('LABEL'//trim(adjustl(id)),required=.TRUE.,error=io_err)      !ATOMS
      IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
      

       string = trim(dict%get_string("ATOM",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
      
       ppfile(is)= trim(dict%get_string("PP_FILE",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
      
       string = trim(dict%get_string("PP_SCHEME",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))     
           IF(INDEX(string,'KLEINMAN').NE.0)THEN
           tKB(Is)=.TRUE.
           ELSE
           tKB(Is)=.FALSE.
           ENDIF

        STRING =trim(dict%get_string("LMAX",error = io_err))
        IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
      
        LMAX(IS)=LVAL(STRING)   

        STRING = (dict%get_STRING("LOCAL",error = io_err))
        IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
     
        L_LOCAL(IS)=LVAL(STRING)
 
        atom_mass(is) = (dict%get_real("MASS",error = io_err))
        IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
      
        atom_p_sp(is) = (dict%get_integer("NUMBER",error = io_err))
        IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
        atom_t=atom_t+atom_p_sp(is)
   END DO


    dict => root%get_dictionary('OUTPUT',required=.TRUE.,error=io_err)         !OUTPUT
    IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
      
       debug = (dict%get_logical("DEBUG",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))

       string = trim(dict%get_string("FILE",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
       
       frequency = (dict%get_real("FREQUENCY",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
      
    dict => root%get_dictionary('DFT',required=.TRUE.,error=io_err)                 !!DFT
    IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
      
       func_string(1)= trim(dict%get_string("EXCHANGE",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
      
       func_string(2) = trim(dict%get_string("CORRELATION",error = io_err))
       IF (associated(io_err)) CALL CHECK_ERROR(trim(io_err%message))
     
       func_string(3) = trim(dict%get_string("XC",error = io_err))
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


  DO Is=1,sp_t
   NA_MAX=MAX(atom_p_sp(Is),NA_MAX)
  ENDDO

END SUBROUTINE reading_input



SUBROUTINE CHECK_ERROR(msg)
  USE yaml
  USE yaml_types
IMPLICIT NONE
CHARACTER(LEN=*)              :: msg
PRINT* ,msg
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



 SUBROUTINE readcoordinate(filename)
   USE system_data_types
   IMPLICIT NONE
   CHARACTER(LEN=*),intent(IN)               :: filename
   INTEGER                                   :: J,I,IA,K
   CHARACTER LINE*80, PATH*42

   ALLOCATE(ATCO(3,NA_MAX,sp_t))
   OPEN(UNIT=12,FILE=trim(PATHOFINPUT)//filename, ACTION = "READ" , STATUS='OLD' )
   READ(12,*)
   READ(12,*)
   DO J=1,sp_t
     DO IA=1,atom_p_sp(J)
     READ(12,*),LABEL(J),(ATCO(I,IA,J),I=1,3)
     K=K+1
     ENDDO
   ENDDO
   CLOSE(12)
 END SUBROUTINE readcoordinate
 
 SUBROUTINE read_and_print
 USE system_data_types
 implicit none
 INTEGER IOS
 character(len=100) :: line ! declare a character variable to store each line of text
  write(6,"(A)") repeat("*", 93)
  write(6, "(3A)")repeat(" ", 40),"INPUT DETAILS",repeat(" ", 40)
  write(6,"(A)") repeat("*", 93)
  open(10, file=trim(PATHOFINPUT)//"input.yaml", status='old', action='read') ! open the input file for reading
  do ! start a loop
    read(10, '(A)', iostat=ios) line ! read a line of text from the file and check the status
    if (ios < 0) exit ! if end of file is reached, exit the loop
    if (ios > 0) stop 'Error reading file' ! if an error occurs, stop the program
    write(*, '(A)') line ! write the line of text to the standard output
  end do ! end the loop
  close(unit=10) ! close the input file
  write(6,"(A)") repeat("*", 93)
END SUBROUTINE read_and_print
END MODULE


