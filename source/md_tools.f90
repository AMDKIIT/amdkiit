MODULE md_tools

USE kinds
USE mympi
USE max_parameter_pp
USE constants
USE system_data_types
USE yamlread
USE random 

IMPLICIT NONE
REAL(KIND=DP) :: sigma_vel!, sigma_sinr, t_ke, lambda_sinr, sinr_mass_1, &
!                 sinr_mass_2, lbylp1, emgt, em2gt, sqe2gt, esinr

REAL(KIND=DP) :: com(3)

CONTAINS

SUBROUTINE init_vel(ndof,temperature,atvel,atke)

INTEGER :: is, ia, ndof
REAL(KIND=DP) :: ran(3), rnd1, rnd2, rnd3
REAL(KIND=DP) :: temperature,atvel(:,:,:),atke

atvel(:,:,:) = 0.0_dp

CALL uniform_rnd(3,ran)

DO is=1,sp_t
   sigma_vel=dsqrt((temperature*kboltz)/atmass(z(is)))
   DO ia=1,atom_p_sp(is)
      rnd1=2.0_dp*pi*ran(1)
      rnd2=2.0_dp*pi*ran(2)
      rnd3=2.0_dp*pi*ran(3)
      
      atvel(1,ia,is)=DSQRT(DLOG(ran(2))*(-2.0_dp))*DCOS(rnd1)*sigma_vel
      atvel(2,ia,is)=DSQRT(DLOG(ran(2))*(-2.0_dp))*DSIN(rnd1)*sigma_vel
      atvel(3,ia,is)=DSQRT(DLOG(ran(3))*(-2.0_dp))*DCOS(rnd2)*sigma_vel
   ENDDO
ENDDO

  CALL temp_ke(ndof,atvel,temp_inst,atke)

  CALL remove_lin_momentum(atvel)
  CALL rescale_atvel(ndof,temperature,atvel,atke) !scale the velocities
  CALL temp_ke(ndof,atvel,temp_inst,atke)
  !WRITE(*,*)"Print within subroutine: Step-3:   ", temp_inst, temperature
END SUBROUTINE init_vel

SUBROUTINE temp_ke(ndof,atvel,temp_inst,atke)
    REAL(KIND=DP) :: atvel(:,:,:), temp_inst, atke

    REAL(KIND=DP) :: mv2
    INTEGER :: is, ia, ndof

    mv2=0.0_dp
    DO is=1,sp_t
       DO ia=1,atom_p_sp(is)
          mv2=mv2+atmass(z(is))*(atvel(1,ia,is)**2+atvel(2,ia,is)**2+atvel(3,ia,is)**2) 
       ENDDO
    ENDDO  
    temp_inst=(mv2*au_to_K)/(real(ndof,kind=dp))
!    temp_inst=(mv2*au_to_K)/((real(ndof,kind=dp))*lbylp1) ! IF SINR
    
    atke=0.5_dp*mv2

END SUBROUTINE temp_ke

SUBROUTINE rescale_atvel(ndof,temperature,atvel,atke)
  REAL(KIND=dp) :: temperature,atvel(:,:,:),atke

  REAL(kind=dp) :: scal
  INTEGER :: is, ia, ndof

  CALL temp_ke(ndof,atvel,temp_inst,atke) !calculate instanteneous temp
  scal=DSQRT(temperature/temp_inst)
    DO is=1,sp_t
       DO ia=1,atom_p_sp(is)
          atvel(1:3,ia,is)=atvel(1:3,ia,is)*scal
       ENDDO
    END DO

  atke=atke*scal*scal

END SUBROUTINE rescale_atvel

SUBROUTINE remove_lin_momentum(atvel)

  REAL(KIND=8) :: atvel(:,:,:)
  !
  REAL(KIND=8) :: rlm(3)
  INTEGER :: ia, is
  REAL(KIND=8) :: totmass

  rlm(1:3)=0.d0
  totmass=0.d0
  DO is=1,sp_t
     DO ia=1,atom_p_sp(is)
        rlm(1)=rlm(1) + atvel(1,ia,is)*atmass(z(is))
        rlm(2)=rlm(2) + atvel(2,ia,is)*atmass(z(is))
        rlm(3)=rlm(3) + atvel(3,ia,is)*atmass(z(is))
        totmass=totmass+atmass(z(is))
     ENDDO
  END DO
  rlm(1:3)=rlm(1:3)/totmass
  DO is=1,sp_t
     DO ia=1,atom_p_sp(is)
        atvel(1,ia,is)=atvel(1,ia,is) - rlm(1)
        atvel(2,ia,is)=atvel(2,ia,is) - rlm(2)
        atvel(3,ia,is)=atvel(3,ia,is) - rlm(3)
     ENDDO
  END DO

END SUBROUTINE remove_lin_momentum

SUBROUTINE get_com(com,atco)
  IMPLICIT NONE
  REAL(KIND=8) :: atco(:,:,:), com(3)
  !
  INTEGER :: ia,is, k
  DO k=1,3
     com(k) = 0.d0
     DO is=1,sp_t
        DO ia=1,atom_p_sp(is)
           com(k)=com(k) + atco(k,ia,is)
        ENDDO
     END DO

     com(k)=com(k)/dble(na_max) 
  END DO
END SUBROUTINE get_com

SUBROUTINE shift_com(com,atco)
  IMPLICIT NONE
  REAL (KIND=8) :: com(3), atco(:,:,:)

  INTEGER :: ia, is, k
  REAL (KIND=8) :: newcom(3)

  CALL get_com(newcom,atco)

     DO is=1,sp_t
        DO ia=1,atom_p_sp(is)
           DO k=1,3
              atco(k,ia,is)=atco(k,ia,is)-(newcom(k)-com(k))
           ENDDO
     END DO
  END DO
END SUBROUTINE shift_com


SUBROUTINE write_md_trajectory(imd,atom_t,atco,atvel)

REAL(KIND=DP) :: atco(:,:,:),atvel(:,:,:)
INTEGER :: atom_t, imd

INTEGER :: is, ia

OPEN(101,FILE="MD_TRAJECTORY.xyz",STATUS="unknown",position="APPEND")
WRITE(101,*)atom_t
WRITE(101,*)"MD STEP=",imd

DO is=1,sp_t
   DO ia=1,atom_p_sp(is)
      WRITE(101,'(A,6f16.6)')&
        symbol(z(is)),atco(1:3,ia,is)/angs2bohr,atvel(1:3,ia,is)/angs2bohr
   ENDDO
ENDDO

CLOSE(101)

END SUBROUTINE write_md_trajectory

SUBROUTINE write_md_energy(imd,temp_inst,atpe,atte,tcpu)

REAL(KIND=DP) :: temp_inst, atpe, atte, tcpu
INTEGER :: imd

OPEN(201,FILE="MD_ENERGY.dat",STATUS="unknown",position="APPEND")
WRITE(201,*)imd,temp_inst,atpe,atte,tcpu

CLOSE(201)

END SUBROUTINE write_md_energy

SUBROUTINE write_md_force(atco,force)

REAL(KIND=dp) :: atco(:,:,:),force(:,:,:)
INTEGER :: iat, ia, is

WRITE(*,'(/,T4,A,T18,A,T48,A)') 'ATOM','COORDINATES (BOHR)',&
         'FORCES (-GRADIENTS)'

iat=0
DO is=1,sp_t
   DO ia=1,atom_p_sp(is)
      iat=iat+1
      WRITE(*,'(1X,I3,2X,A2,3(F10.5),5X,3(1PE11.3))')&
                iat,symbol(z(is)),atco(1:3,ia,is),force(1:3,ia,is)
   ENDDO
ENDDO

END SUBROUTINE write_md_force

SUBROUTINE write_md_info(imd,temp_inst,atpe,atte,tcpu)

REAL(kind=dp) :: temp_inst,atpe,tcpu,atte
INTEGER :: imd

WRITE(*,'(/,5X,A5,A15,A25,A25,A18)')"ISTEP","TEMP (K)", "KS ENERGY","TOT ENERGY","CPU TIME (s)"

WRITE(*,"(I10,5X,F8.1,2X,2(10X,F14.6),8X,F9.2)")imd,temp_inst,atpe,atte,tcpu

END SUBROUTINE write_md_info

END MODULE md_tools
