PROGRAM main_dftk 
  USE kinds
  USE mympi
  USE system_data_types
  USE set_environment
  USE init_module
  USE do_opt
  USE molecular_dynamics
  IMPLICIT NONE
  REAL(KIND=dp):: t_start

  !CHARACTER(80) :: arg1
       INTEGER :: len_arg1, i
  ! Retrieve the first command-line argument

  CALL GET_COMMAND_ARGUMENT(1, input_filename, len_arg1)
  !WRITE(*, *) "Length of first argument:", len_arg1,trim(input_filename),"***"
  CALL setting_env(t_start)
  CALL initialize
  IF(TASK_OPTION .eq.'MOLECULAR_DYNAMICS')THEN
    CALL bomd
  ELSE
    CALL optimize
  ENDIF
  IF(ionode)CALL eprint
  CALL closing_env(t_start)
END PROGRAM main_dftk

!!
