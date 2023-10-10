MODULE set_environment
CONTAINS
  SUBROUTINE setting_env(start)
  USE kinds
  USE mympi
  USE system_data_types
  IMPLICIT NONE
  INTEGER :: values(8)
  REAL(KIND=dp),INTENT(INOUT):: start
  CHARACTER(len=8) :: date, time
  CHARACTER(10) :: date_str
  CHARACTER(8)  :: time_str


  CALL mpi_init
  CALL cpu_time(start)
  CALL date_and_time(values=values)
  IF(ionode)CALL print_logo_amdkiit
  ! Format date and time string
  write(date_str, "(I2.2,'-',I2.2,'-',I4.4)") values(3), values(2),values(1)
  write(time_str, "(I2.2,':',I2.2,':',I2.2)") values(5), values(6),values(7)
 
  ! Print formatted strings
  IF(ionode)write(6,"(5A)")repeat(" ",22),"Program started at : ",trim(date_str),", TIME : ",trim(time_str)
  IF(ionode)write(6,"(A)") repeat("*", 93)

  END SUBROUTINE setting_env

  SUBROUTINE closing_env(start)
  USE kinds
  USE mympi
  USE system_data_types

  IMPLICIT NONE
  INTEGER ::values(8)
  REAL(KIND=dp):: finish,start
  CHARACTER(len=8) :: date, time
  CHARACTER(10) :: date_str
  CHARACTER(8)  :: time_str
 
  CALL date_and_time(values=values)

  ! Format date and time string
  write(date_str, "(I2.2,'-',I2.2,'-',I4.4)") values(3), values(2),values(1)
  write(time_str, "(I2.2,':',I2.2,':',I2.2)") values(5), values(6),values(7)
  IF(ionode)write(6,"(5A)")repeat(" ",1),"Program ended at ",trim(date_str),", Time : ",trim(time_str)
  CALL mpi_final

  CALL cpu_time(finish)
  IF(IONODE) write(*,*) "Elapsed Time = ",(finish-start)/60,"min"
END SUBROUTINE closing_env

SUBROUTINE mpi_init
  USE kinds
  USE mympi
  USE system_data_types
  IMPLICIT NONE
    ! MPI Initialization
  CALL MPI_Start
  CALL MPI_get_ncpu(ncpu)
  CALL MPI_get_cpuid(icpu) !cpu id starts with 1 to ncpu

  ionode=.FALSE.
  IF(icpu==1)ionode=.TRUE.
END SUBROUTINE mpi_init

SUBROUTINE mpi_final
  USE kinds
  USE mympi
  USE system_data_types
  IMPLICIT NONE
  CALL mpi_stop
END SUBROUTINE mpi_final
END MODULE set_environment
