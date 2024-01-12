Module do_opt
contains
subroutine optimize
  USE kinds
  USE mympi
  USE max_parameter_pp
  USE constants
  USE system_data_types
  implicit NONE
  IF(IONODE)write(6, "(3A)")repeat(" ",35),"OPTIMIZATION DETAILS",repeat(" ",27)
  IF(IONODE)write(6,"(A)")repeat("*", 93)
  
  if(TASK_OPTION.eq.'WAVEFUNCTION_OPTIMIZATION')then
    IF(IONODE)WRITE(*,"(A25)") "Wavefunction Optimization"
    IF(IONODE)write(6,"(A)")repeat("-", 93)
    IF(IONODE)WRITE(*,"(A22,A20)") "Minimizer : ",wf_opt%minimizer 
    IF(IONODE)WRITE(*,"(A22,I18)") "Max Steps : ",wf_opt%max_step
    IF(IONODE)WRITE(*,"(A22,E18.6)") "Gradiant Cutoff : ",wf_opt%g_cutoff
    IF(IONODE)WRITE(*,"(A22,E18.6)") "Energy Cutoff : ",wf_opt%e_cutoff
    IF(IONODE)WRITE(*,"(3A)")repeat(" ",28),"Wavefunction Optimization Starts",repeat(" ",28)
    IF(IONODE)write(6,"(A)")repeat("*", 93)
    IF(IONODE)WRITE(*,'(A5,A19,A25,A25,A18)')"ITER","KS ENERGY","ENERGY CHANGE","MAX. GRAIDENT","CPU TIME (ms)"
    IF(IONODE)write(6,"(A)")repeat("*", 93)

  CALL wfn_opt(0)
      
  else
    IF(IONODE)WRITE(*,"(A25)") "Wavefunction Optimization"
    IF(IONODE)write(6,"(A)")repeat("-", 93)
    IF(IONODE)WRITE(*,"(A20,A20)") "Minimizer : ",wf_opt%minimizer
    IF(IONODE)WRITE(*,"(A20,I18)") "Max Steps : ",wf_opt%max_step
    IF(IONODE)WRITE(*,"(A20,E18.6)") "Gradiant Cutoff : ",wf_opt%g_cutoff
    IF(IONODE)WRITE(*,"(A20,E18.6)") "Energy Cutoff : ",wf_opt%e_cutoff

    IF(IONODE)write(6,"(A)")repeat(".", 93)
    IF(IONODE)WRITE(*,"(A23)") "Geometry Optimization"
    IF(IONODE)write(6,"(A)")repeat("-", 93)
    IF(IONODE)WRITE(*,"(A20,A20)") "Minimizer : ", geo_opt%minimizer
    IF(IONODE)WRITE(*,"(A20,I18)") "Max Steps : ",geo_opt%max_step
    IF(IONODE)WRITE(*,"(A20,E18.6)") "Gradiant Cutoff : ",geo_opt%g_cutoff
    IF(IONODE)write(6,"(A)")repeat("-", 93)
    IF(IONODE)WRITE(*,"(3A)")repeat(" ",32),"Geomtry Optimization Starts",repeat(" ",32)
    IF(IONODE)write(6,"(A)")repeat("*", 93)
    IF(IONODE)WRITE(*,'(A5,A19,A25,A25,A18)')"ITER","KS ENERGY","ENERGY CHANGE","MAX. GRAIDENT","CPU TIME (ms)"
    IF(IONODE)write(6,"(A)")repeat("*", 93)
 
    IF(geo_opt%minimizer=="LBFGS")CALL geom_opt_lbfgs
  endif
end subroutine optimize

subroutine wfn_opt(geo_step)
  USE kinds
  USE mympi
  USE max_parameter_pp
  USE constants
  USE system_data_types
  USE bessel_func
  USE atomic_basis
  USE density
  USE exp_igr
  USE potential
  USE xc
  USE wfn_initialize
  USE gvectors
  USE fft_interface, ONLY: prepare_fft
  USE pseudopotential
  USE total_energy
  USE gradient
  USE nuclear_grad
  USE cg
  USE bfgs_mod
  USE FUNCTION_VAL
  USE backtracking

  implicit NONE
  INTEGER :: i
  REAL(KIND=DP)::etotal_prev,w_start,w_end,detot
  integer,intent(in)::geo_step
  etotal_prev=0.0d0
  DO I=1,wf_opt%max_step
  CALL cpu_time(w_start)
  CALL eval_density
  call rhor2g
  CALL eval_pot
  CALL cal_vpot_ex
  gopt=.FALSE.
  CALL pp_energy
  CALL e_total
  C2(1:NGPW_L:1,1:NSTATE)=(0._dp,0._dp)
  CALL VPSI(C_0,C2,OCCUPATION,NSTATE)
  if(l_upf)then
  CALL UPF_FNONLOC(C2,NSTATE,OCCUPATION)
  else
  CALL FNONLOC(C2,NSTATE,OCCUPATION)
  endif
  CALL FORCES
  CALL CHK_TFOR(I,GEMAX)
  if(tfor) then
  CALL E_SR
  CALL POTFOR!(force)
   if(l_upf)then
     CALL upf_NLFOR
   else
     CALL NLFOR
   endif
  endif
  !DEALLOCATE(V_POT) !TODO
  CALL PCG(I)

  !!!CALL eval_GMAX(C2)!gradient calculation and convergence chk
  !!!call update_wf
  CALL ortho_gs(c_0)
!###RITAMA####
!  detot=etotal_prev-etotal 
!  IF(I==1)detot=0.0d0
!###RITAMA####
  DEALLOCATE(V_POT,twnl,nl)
  CALL cpu_time(w_end)
  IF(IONODE) write(6,'(A,I4,3X,F16.6,6X,F16.6,12X,E16.6,8X,F10.1)')"#",i,etotal,abs(etotal_prev-etotal),GEMAX,(w_end-w_start)*1000
    IF(IONODE.AND.TASK_OPTION .eq.'WAVEFUNCTION_OPTIMIZATION')THEN
  OPEN(7,FILE="SP_ENERGY.dat",STATUS="unknown",position="APPEND")
    write(7,'(I4,3X,F16.6,8X,F10.1)')i,etotal,(w_end-w_start)*1000
  CLOSE(7)
  ENDIF
  etotal_prev=etotal
  IF(CONVWF) THEN
          GOTO 101
  ENDIF
  ENDDO
  DEALlOCATE(HNM1)
  101 CONTINUE
  if(ncpu==1.and.Lopen_shell.and.print_density) then
    call cubefile
  else if(ionode.and.Lopen_shell.and.print_density) then
    write(6,*)"Run serial for rho_spin.cube !! "
  endif
  return
end subroutine wfN_opt
subroutine geom_opt
  USE kinds
  USE mympi
  USE max_parameter_pp
  USE constants
  USE system_data_types
  USE exp_igr
  USE potential
  USE xc
  USE wfn_initialize
  USE nuclear_grad
  USE bfgs_mod
  USE FUNCTION_VAL
  USE backtracking
  IMPLICIT NONE
  INTEGER iopt,i
  REAL(KIND=dp), SAVE,POINTER :: x0(:)
  REAL(KIND=dp), SAVE,POINTER :: g0(:)
  REAL(KIND=dp), SAVE,POINTER :: HESS_IN(:,:)
  REAL(KIND=dp) ::gnorm, energyold,g_start,g_end
  SAVE energyold
  ALLOCATE(x0(3*atom_t),g0(3*atom_t),HESS_IN(3*atom_t,3*atom_t))
  CONVGEO=.FALSE.
  energyold=0.0d0
  iopt =0
  opt_loop: DO
  IF(iopt>geo_opt%max_step)EXIT opt_loop
        CALL cpu_time(g_start)
        CALL wfn_opt(iopt)
        !energyold=etotal
        CALL MPI_GlobSumR2(force,3*NA_MAX*sp_t)
        CALL convert3dto1d(ATCO,3,NA_max,sp_t,x0,3*atom_t)
        CALL convert3dto1d(-force,3,NA_max,sp_t,g0,3*atom_t)
        write(90,*)x0
        write(91,*)g0
        CALL PRINT_COORDINATE(ATCO,-force,NA_max,sp_t,atom_t,etotal,iopt)
        CALL GEOOPT_bfgs(iopt,x0,g0,HESS_IN,3*atom_t,energyold,CONVGEO)
        CALL cpu_time(g_end)
        IF(IONODE)write(6,"(A)")repeat("*", 93)
        IF(IONODE)write(6,*)"GEOMETRY UPDATE:"!,etotal,energyold
        IF(IONODE)write(6,'(I5,3X,F16.6,6X,F16.6,12X,E16.6,8X,F10.1)')iopt,etotal,abs(energyold-etotal),GNMAX,(g_end-g_start)*1000
        IF(IONODE)write(6,"(A)")repeat("*", 93)
        IF(CONVGEO)EXIT opt_loop
        CALL convert1dto3d(x0,3*atom_t,ATCO,3,NA_max,sp_t)
        !CALL PRINT_COORDINATE(ATCO,NA_max,sp_t,atom_t,etotal,iopt)
        CALL exp_igrxfactor
        CALL ortho_gs(c_0)
      iopt=iopt+1
      !IF(IONODE)write(6,"(A)")repeat("*", 93)
      !IF(IONODE)write(6,*)"GEOMETRY UPDATE:"!,etotal,energyold
      !  CALL cpu_time(g_end)
     ! IF(iopt==1)energyold=etotal
      !IF(IONODE)write(6,'(I5,3X,F16.6,6X,F16.6,12X,E16.6,8X,F10.1)')iopt,etotal,abs(energyold-etotal),GNMAX,(g_end-g_start)*1000
      !IF(IONODE)write(6,"(A)")repeat("*", 93)
      energyold=etotal
   END DO opt_loop
   deallocate(hess_in)
end subroutine geom_opt
subroutine geom_opt_lbfgs
  USE kinds
  USE mympi
  USE max_parameter_pp
  USE constants
  USE system_data_types
  USE exp_igr
  USE potential
  USE xc
  USE wfn_initialize
  USE nuclear_grad
  USE backtracking
  IMPLICIT NONE
  REAL(KIND=dp), SAVE,POINTER :: x0(:)
  REAL(KIND=dp), SAVE,POINTER :: g0(:)
  REAL(KIND=dp), SAVE,POINTER :: DIAG(:)
  REAL(KIND=dp), SAVE,POINTER :: W(:)
  INTEGER NWORK,MSAVE
  REAL(KIND=dp):: XTOL,GTOL,T1,T2,STPMIN,STPMAX
  INTEGER IPRINT(2),IFLAG,IOPT,N,M,MP,LP,J
  LOGICAL DIAGCO,i_opt

  REAL(KIND=dp) ::gnorm, energyold,g_start,g_end
  SAVE energyold
  EXTERNAL LB2
  COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX
      i_opt=.TRUE.
      MSAVE=7
      M=5
      IPRINT(1)= 1
      IPRINT(2)= 1
      DIAGCO= .FALSE.
      XTOL= 1.0D-16

      IOPT=0
      IFLAG=0

  ALLOCATE(x0(3*atom_t),g0(3*atom_t),DIAG(3*atom_t))
  NWORK=(3*atom_t)*(2*MSAVE +1)+2*MSAVE
  ALLOCATE(W(NWORK))
  !CONVGEO=.FALSE.
  energyold=0.0d0

  20  CONTINUE
        CALL cpu_time(g_start)
        CALL wfn_opt(iopt)
        CALL MPI_GlobSumR2(force,3*NA_MAX*sp_t)
        CALL convert3dto1d(ATCO,3,NA_max,sp_t,x0,3*atom_t)
        force=-force
        CALL convert3dto1d(force,3,NA_max,sp_t,g0,3*atom_t)
        if(i_opt) then 
               CALL PRINT_COORDINATE(ATCO,force,NA_max,sp_t,atom_t,etotal,iopt)
               iopt=iopt+1
               IF(IONODE)Write(6,"(3A)")repeat("*", 38)," GEOMETRY UPDATE ",repeat("*",38)
        endif
        force=-force
        CALL LBFGS(3*atom_t,M,X0,etotal,G0,DIAGCO,DIAG,IPRINT,geo_opt%g_cutoff,XTOL,W,IFLAG,i_opt)
        IF(IFLAG.LE.0) GO TO 50
        CALL cpu_time(g_end)
        CALL convert1dto3d(x0,3*atom_t,ATCO,3,NA_max,sp_t)
        DEALLOCATE(eigr,eigrxrhos,eigrxvps,eigr_pw)
        CALL exp_igrxfactor
        CALL ortho_gs(c_0)
        energyold=etotal
        deallocate(force)
        IF(iopt.GT.geo_opt%max_step) GO TO 50
      GO TO 20
  50  CONTINUE
  deallocate(atco,x0,g0,diag,w)
end subroutine geom_opt_lbfgs
end module do_opt


