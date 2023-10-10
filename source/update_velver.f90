MODULE update_velver

USE kinds
USE mympi
USE max_parameter_pp
USE constants
USE system_data_types
USE yamlread
USE random 
USE md_tools

CONTAINS

SUBROUTINE velup_velver(vel,fion)
  IMPLICIT NONE
  !
  REAL(KIND=dp),INTENT(INOUT) :: vel(:,:,:)
  REAL(KIND=dp),INTENT(IN) :: fion(:,:,:)
  REAL(kind=dp) :: fact
  !
  INTEGER :: ia, is
  !
    DO is=1,sp_t
      fact=(0.5d0*time_step)/atmass(z(is))
      DO ia=1,atom_p_sp(is)
        vel(1,ia,is)=vel(1,ia,is)+fact*fion(1,ia,is)
        vel(2,ia,is)=vel(2,ia,is)+fact*fion(2,ia,is)
        vel(3,ia,is)=vel(3,ia,is)+fact*fion(3,ia,is)
      END DO
    ENDDO
END SUBROUTINE velup_velver

SUBROUTINE posup_velver(pos,vel)
  IMPLICIT NONE
  !
  REAL(KIND=dp),INTENT(INOUT) :: pos(:,:,:)
  REAL(KIND=dp),INTENT(IN) :: vel(:,:,:)
  !
  INTEGER :: ia,is
  !

  DO is=1,sp_t
     DO ia=1,atom_p_sp(is)
        pos(1,ia,is)=pos(1,ia,is)+vel(1,ia,is)*time_step
        pos(2,ia,is)=pos(2,ia,is)+vel(2,ia,is)*time_step
        pos(3,ia,is)=pos(3,ia,is)+vel(3,ia,is)*time_step
     ENDDO
  END DO
END SUBROUTINE posup_velver

END MODULE update_velver
