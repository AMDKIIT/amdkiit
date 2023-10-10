MODULE molecular_dynamics

CONTAINS

SUBROUTINE bomd
USE kinds
USE mympi
USE max_parameter_pp
USE constants
USE system_data_types
USE yamlread
USE random 
USE exp_igr
USE wfn_initialize
USE md_tools
USE do_opt
!USE update_sinr
USE update_velver

   IMPLICIT NONE
   INTEGER :: imd, ia, is, iat, minmd, loopimd
   REAL(KIND=dp) :: tcpu0,tcpu1


temperature= md%temperature 
maxmd= md%max_step
time_step_fs= md%time_step  

time_step=time_step_fs/au_to_fs !convert time step to au


   CALL cpu_time(tcpu0)

!   imd=0

   CALL get_com(com,atco)

!   ndof=3*atom_t   !IF SINR
   ndof=3*atom_t-3 
   ALLOCATE(atvel(3,na_max,sp_t)) 

!   CALL init_sinr  ! IF SINR
   
   CALL init_vel(ndof,temperature,atvel,atke)
!   CALL vel_sinr(atvel,v1_sinr,v2_sinr) !IF SINR

   CALL temp_ke(ndof,atvel,temp_inst,atke)

!   CALL sinr_cons(v1_sinr,atvel,esinr)  !IF SINR

   IF(IONODE)write(*,"(A)")repeat("-", 93)
   IF(ionode)WRITE(*,"(3A)")repeat(" ",36),"FORCES INITIALIZATION",repeat(" ",36)
   IF(IONODE)write(*,"(A)")repeat("-", 93)
   IF(IONODE)WRITE(*,'(A5,A19,A25,A25,A18)')"ITER","KS ENERGY","ENERGY CHANGE","MAX. GRAIDENT","CPU TIME (ms)"
   
   imd=0
!############################## FORCE CALCULATION ##############################
   CALL wfn_opt(imd)
   CALL MPI_GlobSumR2(force,3*NA_MAX*sp_t)
!############################## FORCE CALCULATION ##############################
  
   IF(IONODE)write(*,"(A)")repeat("-", 93)
   IF(ionode)WRITE(*,"(3A)")repeat(" ",32),"END OF FORCES INITIALIZATION",repeat(" ",32)
  
   atpe=etotal 
   atte=atke+atpe

   CALL cpu_time(tcpu1)

!   IF(ionode)THEN
 IF(IONODE)THEN
    CALL write_md_force(atco,force)
    write(6,"(A)")repeat("=", 93)
     CALL write_md_trajectory(imd,atom_t,atco,atvel)
     CALL write_md_energy(imd,temp_inst,etotal,atte,tcpu1-tcpu0)
!     CALL write_md_sinr(imd,v1_sinr,v2_sinr,esinr) !IF SINR
!   ENDIF

     WRITE(6,'(A,T50,F8.2,A8)') ' TIME FOR INITIALIZATION:',tcpu1-tcpu0,' SECONDS'
     write(6,"(A)")repeat("=", 93)
     WRITE(6,"(3A)")repeat(" ",32), "MOLECULAR DYNAMICS SIMULATION",repeat(" ",32)
     write(6,"(A)")repeat("=", 93)
     WRITE(6,'(A,T50,F8.2,A8)') "TEMPERATURE:", temp_inst, ' K'
     WRITE(6,'(A,T50,F8.2,A8)') "MD TIME STEP:", time_step, ' A.U.'
     WRITE(6,'(A,T50,I8)') "NUMBER OF MD STEP :", maxmd
     WRITE(6,'(A,T50,I8)') "DEGREES OF FREEDOM :", ndof
     write(6,"(A)")repeat("=", 93)
 ENDIF

     minmd=imd+1 ! START MD STEP

     DO loopimd=minmd,maxmd

     CALL cpu_time(tcpu0)

     imd=imd+1

     CALL velup_velver(atvel,force)
     
     CALL posup_velver(atco,atvel)

     CALL shift_com(com,atco)

!###################### FORCE CALCULATION #######################!
     IF(IONODE)WRITE(*,'(A5,A19,A25,A25,A18)')"ITER","KS ENERGY","ENERGY CHANGE","MAX. GRAIDENT","CPU TIME (ms)"
     DEALLOCATE(eigr,eigrxrhos,eigrxvps,eigr_pw)
     CALL exp_igrxfactor   !update structure factor
     CALL ortho_gs(c_0)
     CALL wfn_opt(imd)  
     CALL MPI_GlobSumR2(force,3*NA_MAX*sp_t)
!###################### FORCE CALCULATION #######################!

     CALL velup_velver(atvel,force)

     CALL temp_ke(ndof,atvel,temp_inst,atke)
    
     atpe=etotal
     atte=atke+atpe

     CALL cpu_time(tcpu1)
    

     IF(ionode)THEN
        CALL write_md_force(atco,force)
        CALL write_md_info(imd,temp_inst,atpe,atte,tcpu1-tcpu0)
        CALL write_md_trajectory(imd,atom_t,atco,atvel)
        CALL write_md_energy(imd,temp_inst,atpe,atte,tcpu1-tcpu0)
!        CALL write_md_sinr(imd,v1_sinr,v2_sinr,esinr) !IF SINR
     ENDIF
   ENDDO !MD_LOOP

   DEALLOCATE(atvel,force,atco) 

END SUBROUTINE bomd

END MODULE molecular_dynamics
