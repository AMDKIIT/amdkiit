SUBROUTINE read_pp_file
USE system_data_types
USE readstring
USE readupf 
USE readpsp

IMPLICIT NONE

INTEGER :: IS,count,count_p
LOGICAL :: l_PSP
ALLOCATE(z(sp_t),zv(sp_t),xc_fun(sp_t),meshv(sp_t),meshw(sp_t))!,slatr_ex(sp_t),tnum(sp_t))
count=0
count_p=0
L_UPF=.FALSE.
l_psP=.faLSE.

 DO IS=1,sp_t
   IF(index(PPFILE(is),"upf").ne.0)then
    IF(count==0)call upf
    count=count+1
    call read_upf(PPFILE(is),is)
    l_upf=.TRUE.
   ELSE
    IF(count_p==0)call psp
    count_p=count_p+1
    call read_psp(PPFILE(is),is)
    l_psp=.TRUE.
   ENDIF
 ENDDO 

 if(l_upf.AND.l_psp)then
   print*,"ERROR!! pp files are different"
   STOP
 endif

RETURN

END SUBROUTINE read_pp_file
