MODULE function_val

CONTAINS
!Computing Energy
  FUNCTION func(x) RESULT(f)
    USE kinds
    USE system_data_types
    USE density
    USE exp_igr
    USE potential
    USE pseudopotential
    USE total_energy
    USE gradient
    USE xc
    use wfn_initialize
      USE nuclear_grad
    IMPLICIT NONE
    INTEGER :: i, j,k,is,n
    REAL(KIND=dp), DIMENSION(:), POINTER :: x
    REAL(KIND=dp) :: f

    !CALL form_factor!exp_igrxrhos
   I=1
   DO J=1,sp_t
    DO K=1,atom_p_sp(J)
     DO IS=1,3
       ATCO(IS,K,J)=x(i)
       I=I+1
     ENDDO
    ENDDO
   ENDDO
 write(45,*)"#",eke!,"=",enl+eloc+eke+exc+estat
   CALL exp_igrxfactor
!   CALL ortho_gs(c_0)
!   CALL eval_density
!   call rhor2g
!   CALL eval_pot
!  CALL cal_vpot_ex
  gopt=.TRUE.
  CALL pp_energy
!  CALL e_total
   CALL e_hartree
   CALL e_sr
  estat= (DREAL(eh)*cell_volume) + esr-eself!/DBLE(ncpu)!nproc)
!  CALL e_ke
!  CALL e_xc
! if(icpu==1)then
  !etotal=enl+eloc+eke+exc+estat
 write(45,*)"#",eke!,exc,estat
 write(45,*)"ETOTAL",etotal
 write(45,*)"f=",enl+eloc+eke+exc+estat
    f=enl+eloc+eke+exc+estat
!    print*,x(1),x(2),x(3)

  END FUNCTION func

!Computing Gradient

END MODULE
