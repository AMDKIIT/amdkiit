PROGRAM main_dftk 
  USE kinds
  USE mympi
  USE system_data_types
  USE set_environment
  USE init_module
  USE do_opt

  IMPLICIT NONE
  REAL(KIND=dp):: t_start
  CALL setting_env(t_start)
  CALL initialize
  CALL optimize
  CALL eprint 
  CALL closing_env(t_start)
END PROGRAM main_dftk

!!
