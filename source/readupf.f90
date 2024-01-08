      MODULE READupf
      CONTAINS

      SUBROUTINE upf
      USE system_data_types
      IMPLICIT NONE
      ALLOCATE(upf_element(sp_t),upf_typ(sp_t),upf_rel(sp_t),upf_dft(sp_t))
      ALLOCATE(upf_meshv(sp_t),upf_lmax(sp_t))
      ALLOCATE(upf_local(sp_t),upf_l_coulomb(sp_t),upf_zp(sp_t))
      ALLOCATE(upf_nwfc(sp_t),upf_nbeta(sp_t))
      ALLOCATE(upf_r(splnmax,sp_t),upf_rab(splnmax,sp_t))
      ALLOCATE(upf_vloc(splnmax,sp_t))
      ALLOCATE(rho_atom(splnmax,sp_t))
      
      upf_nbeta=0
      END SUBROUTINE upf
     
!     ===================================================================================
      SUBROUTINE read_upf(filename,id_sp)
      USE system_data_types
      IMPLICIT NONE
      CHARACTER, INTENT(in) :: filename*(*)
      integer, INTENT(IN):: ID_SP
      CHARACTER (len=80) line
      CHARACTER (len=256) string
      REAL(KIND=DP)::TEST
      INTEGER   IUNIT, IOS,i,is,nhxs
      PARAMETER (IUNIT=20)
      LOGICAL ERREAD,BACK
      
      OPEN(UNIT=IUNIT,FILE=trim(PATHOFINPUT)//filename,IOSTAT=ios)
      !OPEN(UNIT=IUNIT,FILE=trim(PATHOFINPUT)//'PSEUDOPOTENTIAL/'//filename,IOSTAT=ios)
      IF(ios.NE.0)THEN
        WRITE(*,*)"    Reading ",filename," is not successful"
        STOP
      ENDIF
      !!---------------------------------------------------------------------------------
      !------->Search for Header
      CALL search(iunit, "HEADER", .true.)
      CALL read_header (id_sp, iunit)
      !!---------------------------------------------------------------------------------
      !------->Search for Header
      CALL search(iunit, "MESH", .true.)
      CALL search(iunit, "R", .true.)
      CALL read_rr(id_sp,iunit,upf_r(:,id_sp))
      CALL search(iunit, "RAB", .true.)
      CALL read_rr(id_sp,iunit,upf_rab(:,id_sp))
      !CALL search_end (iunit, "HEADER")
      !!---------------------------------------------------------------------------------
      !------->Read local potential
      IF(.not. upf_l_coulomb(id_sp)) THEN
         CALL search(iunit, "LOCAL", .true.)
         CALL read_rr(id_sp, iunit,upf_vloc(:,id_sp))
      ENDIF
      !!---------------------------------------------------------------------------------
      !------->Read nonlocal potential
      CALL search(iunit, "NONLOCAL", .true.)
      CALL read_nonlocal(id_sp, iunit)
      !!---------------------------------------------------------------------------------
      !------->Search for atomic wavefunctions
      CALL search(iunit, "PSWFC", .true.)
      !CALL search(iunit, "CHI", .true.)
       
      !ALLOCATE(upf_oc(upf_nwfc(id_sp),sp_t),upf_epseu(upf_nwfc(id_sp),sp_t))
      !ALLOCATE(upf_rcut_chi(upf_nwfc(id_sp),sp_t))
      ALLOCATE(upf_chi(upf_meshv(id_sp),sp_t,upf_nwfc(id_sp)))
      CALL read_wfc(id_sp, iunit)
      !----------------------------------------------------------------------------------
      !-------->Search for atomic charge   
      CALL search(iunit, "RHOATOM", .true.)
      CALL read_rhoatom(id_sp, iunit)
      !----------------------------------------------------------------------------------
      CALL convert_upf_to_psp(id_sp)
      zv(id_sp)=upf_zp(id_sp)
      meshv(id_sp)=upf_meshv(id_sp)
      meshw(id_sp)=meshv(id_sp)
 
      call atom_info
      DO i=1,99
      IF(trim(SYMBOL(i)).EQ.upf_element(id_sp))Z(id_sp)=i
      ENDDO 
  
      close(iunit)
      RETURN
      END SUBROUTINE read_upf
!     ===================================================================================
      SUBROUTINE search (id, string, l_rewind)
      !----------------------------------------------------------------------------------
      implicit none
      integer :: id
      character (len=*) :: string  ! Label to be matched
      logical :: l_rewind  ! Flag: if .true. rewind the file
      character (len=80) :: read_string ! String read from file
      integer :: ios

      ios = 0
      if (l_rewind) rewind (id)
      do while (ios.eq.0)
      read (id,iostat = ios, err = 99,FMT='(A)') read_string
      if (INDEX(read_string,"<PP_"//string).NE.0 ) then
      !!print*,read_string
      return
      endif
      enddo
      99 print*, 'ERROR!!! NO ',string,' block found', abs (ios)
      end subroutine search
!     ===================================================================================    
      SUBROUTINE search_end (id, string, l_rewind)
      !----------------------------------------------------------------------------------
      implicit none
      integer :: id
      character (len=*) :: string  ! Label to be matched
      logical :: l_rewind  ! Flag: if .true. rewind the file
      character (len=80) :: read_string ! String read from file
      integer :: ios

      ios = 0
      if (l_rewind) rewind (id)
      do while (ios.eq.0)
      read (id,iostat = ios, err = 99,FMT='(A)') read_string
      if (INDEX(read_string,"</PP_"//string).NE.0 ) return
      enddo
      99 print*, 'ERROR!!! No ',string,read_string,' block found', abs (ios)

      end subroutine search_end
!     ===================================================================================   
      SUBROUTINE read_header(i_sp,id)
      USE system_data_types
      IMPLICIT NONE
      INTEGER i_sp,id
      CHARACTER (len=100) dummy_var
      logical l_header_end
      l_header_end=.false.
      CALL HEADER(ID,"author",dummy_var,l_header_end)
      CALL HEADER(ID,"date",dummy_var,l_header_end)
      CALL HEADER(ID,"comment",dummy_var,l_header_end)      
      CALL HEADER(ID,"element",upf_element(i_sp),l_header_end)    
      CALL HEADER(ID,"pseudo_type",upf_typ(i_sp),l_header_end)
      CALL HEADER(ID,"relativistic",upf_rel(i_sp),l_header_end)
      CALL HEADER(ID,"is_ultrasoft",dummy_var,l_header_end)
      CALL HEADER(ID,"is_paw",dummy_var,l_header_end)
      CALL HEADER(ID,"is_coulomb",dummy_var,l_header_end)
      READ(dummy_var,*)upf_l_coulomb(i_sp)
      CALL HEADER(ID,"has_so",dummy_var,l_header_end)
      CALL HEADER(ID,"has_wfc",dummy_var,l_header_end)
      CALL HEADER(ID,"has_gipaw",dummy_var,l_header_end)
      CALL HEADER(ID,"paw_as_gipaw",dummy_var,l_header_end)
      CALL HEADER(ID,"core_correction",dummy_var,l_header_end)
      CALL HEADER(ID,"functional",upf_dft(i_sp),l_header_end)
      CALL HEADER(ID,"z_valence",dummy_var,l_header_end)
      READ(dummy_var,*)upf_zp(i_sp)
      CALL HEADER(ID,"total_psenergy",dummy_var,l_header_end)
      CALL HEADER(ID,"wfc_cutoff",dummy_var,l_header_end)
      CALL HEADER(ID,"rho_cutoff",dummy_var,l_header_end)
      CALL HEADER(ID,"l_max",dummy_var,l_header_end)
      READ(dummy_var,*)upf_lmax(i_sp)
      CALL HEADER(ID,"l_max_rho",dummy_var,l_header_end)
      CALL HEADER(ID,"l_local",dummy_var,l_header_end)
      READ(dummy_var,*)upf_local(i_sp)
      CALL HEADER(ID,"mesh_size",dummy_var,l_header_end)
      READ(dummy_var,*)upf_meshv(i_sp)
      CALL HEADER(ID,"number_of_wfc",dummy_var,l_header_end)
      READ(dummy_var,*)upf_nwfc(i_sp) 
      CALL HEADER(ID,"number_of_proj",dummy_var,l_header_end)
      READ(dummy_var,*)upf_nbeta(i_sp)
       upf_lmax(i_sp)=upf_lmax(i_sp)+1
      END SUBROUTINE read_header
!     ===================================================================================     
      SUBROUTINE header(id,string,var,l_search_end)
      USE system_data_types
      IMPLICIT NONE
      CHARACTER, INTENT(in) :: string*(*)
      CHARACTER, INTENT(out):: var*(*)
      CHARACTER (len=100) :: read_string
      INTEGER IOS,ID
      logical l_search_end
      ios=0
      do while (ios.eq.0)
       READ(id,iostat = ios,FMT='(A)')read_string
         IF(INDEX(read_string,string).NE.0 ) then  
          var =trim(read_string(INDEX(read_string, '"')+1:INDEX(read_string, '"',BACK =.TRUE.)-1))
          if(INDEX(read_string,">").NE.0 ) l_search_end=.true.
          return
         ELSE 
           if(INDEX(read_string,"/>").NE.0 )then      
           !print*, 'ERROR!!! "',string,'"',"  not found"
           rewind(id)
           EXIT 
           endif
         endif
       enddo
      END SUBROUTINE header
!     ===================================================================================     

      SUBROUTINE header_test(idt,string,var,l_search_end)
      USE system_data_types
      IMPLICIT NONE
      CHARACTER, INTENT(in) :: string*(*)
      CHARACTER, INTENT(out):: var*(*)
      CHARACTER (len=100) :: read_string
      INTEGER IOS,IDt
      logical l_search_end
       READ(idt,iostat = ios,FMT='(A)')read_string
         if(INDEX(read_string,string).NE.0 ) then
          var =trim(read_string(INDEX(read_string,'"')+1:INDEX(read_string, '"',BACK =.TRUE.)-1))
         else
            !print*, 'ERROR!!! "',string,'"',"  not found",ios!,l_search_end,read_string
            var="null"
            backspace(idt)
         end if
         if(INDEX(read_string,">").NE.0 )l_search_end=.TRUE.
      END SUBROUTINE header_test
!     ===================================================================================  

      SUBROUTINE read_rr(i_sp,id,READ_VAR)
      USE system_data_types
      IMPLICIT NONE
      INTEGER i_sp,id
      REAL(KIND=DP),intent(out) ::read_var(upf_meshv(i_sp))
      INTEGER ir,IOS
      CHARACTER (len=100) :: var
      
      !WRITE(*,*)upf_meshv(i_sp),size(read_var)
      !if(i_sp==2)then
      !read(id,iostat = ios,FMT='(A)')var
      !print*,var,id,size(read_var)!(upf_meshv(i_sp),i_sp))
      !stop
      !else
      !print*,"*",read_var(1),read_var(upf_meshv(i_sp)),i_sp
      read (id, *, err = 100, iostat = ios)(read_var(ir), ir=1,upf_meshv(i_sp))
      
      !print*,read_var(1),read_var(upf_meshv(i_sp)),i_sp
      return
      100 print*,"ERROR",ios
      END SUBROUTINE read_rr
!     ===================================================================================        
      SUBROUTINE read_nonlocal(i_sp,id)
      USE system_data_types
      USE readstring
      IMPLICIT NONE
      INTEGER i_sp,id,IV,LP,K,L,M,j,jv,i,A,E
      INTEGER ir,nb,IOS,NHXS,ierr,nproje
      character (len=80)dummy_var,var
      logical l_header_end
      integer upf_size,upf_index,upf_l
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)  :: dij_tmp
      LOGICAL EREAD

      !upf_rcut=0._dp
      !upf_rcutus=0._dp
      DO nb=1,upf_nbeta(I_SP)
       l_header_end=.false.
       READ(id,iostat = ios,FMT='(A)')dummy_var!read_string
       IF(INDEX(dummy_var,"<PP_BETA").EQ.0 )THEN
            PRINT*,"PP_BETA NOT FOUND",dummy_var
            STOP
       ENDIF
 
       CALL HEADER_TEST(ID,"type",dummy_var,l_header_end)
       CALL HEADER_TEST(ID,"size",dummy_var,l_header_end)
       READ(dummy_var,*)upf_size
       CALL HEADER_TEST(ID,"columns",dummy_var,l_header_end)
       CALL HEADER_TEST(ID,"index",dummy_var,l_header_end)
       READ(dummy_var,*)upf_index
       CALL HEADER_TEST(ID,"label",dummy_var,l_header_end)
       CALL HEADER_TEST(ID,"angular_momentum",dummy_var,l_header_end)
       READ(dummy_var,*)upf_l
       !if(dummy_var.ne."null") READ(dummy_var,*)upf_lll(nb,i_sp)
       CALL HEADER_TEST(ID,"cutoff_radius_index",dummy_var,l_header_end)
       !if(dummy_var.ne."null")then
       !  READ(dummy_var,*)upf_kbeta(nb,i_sp)
       !else 
       !  upf_kbeta(nb,i_sp)=upf_meshv(i_sp)
       !endif
       CALL HEADER_TEST(ID,"cutoff_radius",dummy_var,l_header_end)
       !if(dummy_var.ne."null") READ(dummy_var,*)upf_rcut(nb,i_sp)
!       CALL HEADER_TEST(ID,"ultrasoft_cutoff_radius",dummy_var,l_header_end)
       !if(dummy_var.ne."null") READ(dummy_var,*)upf_rcutus(nb,i_sp)
       read (id, *, err = 104, iostat = ios)(upf_gnl(ir,i_sp,nb),ir=1,upf_meshv(i_sp))
            READ(id,iostat = ios,FMT='(A)')dummy_var!read_string
       if(INDEX(dummy_var,"</PP_BETA").EQ.0 )THEN
            PRINT*,"/PP_BETA NOT FOUND",dummy_var
            STOP
       ENDIF
      !  DO i=1,lmaxx
      !    npro(i,i_sp)=0
      ! ENDDO

      !==================================================
       PPLIST(upf_index,1,I_SP)=upf_l
       NPRO(upf_l+1,I_SP)=NPRO(upf_l+1,I_SP)+1
       !print*,"**",upf_l,NPRO(upf_l+1,I_SP)
       IF(NPRO(upf_l+1,I_SP).GT.4) THEN
         print*,'RECPUPF','TOO MANY PROJECTORS'
      ENDIF
       PPLIST(upf_index,2,I_SP)=NPRO(upf_l+1,I_SP)
      ENDDO
       !C..set up 
       IV=0
       LP=0
       DO L=1,upf_LMAX(I_SP)
         DO M=1,2*L-1
          LP=LP+1
           DO K=1,NPRO(L,I_SP)
            IV=IV+1
            UPF_NGHTOL(IV,I_SP)=L-1
             LPVAL(IV,I_SP)=LP
             LFVAL(IV,I_SP)=K
           ENDDO
         ENDDO
       ENDDO
       UPF_NGH(I_SP)=IV
      !==================================================
      READ(id,iostat = ios,FMT='(A)')dummy_var!read_string
      if(INDEX(dummy_var,"<PP_DIJ").EQ.0 )THEN
            PRINT*,"PP_DIJ NOT FOUND",dummy_var
            STOP
      endif
          DO l=1,lmax(i_sp)
          DO i=1,mpro
             DO j=1,mpro
               hlsg(i,j,l,i_sp)=0.0d0
             ENDDO
          ENDDO
       ENDDO
       !IF (upf2) THEN
       nproje=upf_nbeta(i_sp)
          ALLOCATE(dij_tmp(nproje*nproje),STAT=ierr)
          IF(ierr/=0) print*,'allocation problem'
          READ(id,*, err = 104, iostat = ios)(dij_tmp(ir),ir=1,nproje*nproje)
          ir=0
          DO i=1,nproje
             DO j=1,nproje
                ir=ir+1
                IF (i/=j) THEN
                   IF (dij_tmp(ir)/=0.0d0) THEN
                      print*,'RECPUPF','NON-ZERO DIJ IS INVALID!'
                   ENDIF
                   CYCLE
                ENDIF
                IF (pplist(i,1,i_sp).NE.pplist(j,1,i_sp)) THEN
                   print*,'INVALID DIJ'
                ENDIF
                l=pplist(i,1,i_sp)+1
                iv=pplist(i,2,i_sp)
                jv=pplist(j,2,i_sp)
                hlsg(iv,jv,l,i_sp)=.5*dij_tmp(ir)
                hlsg(jv,iv,l,i_sp)=.5*dij_tmp(ir)
                !print*,hlsg(iv,jv,l,i_sp),iv,jv,l,i_sp
                !print*,hlsg(jv,iv,l,i_sp),jv,iv
             ENDDO
          ENDDO
      !==================================================
      return
      104 print*,"error1"
      END SUBROUTINE read_nonlocal
!     ===================================================================================        
      SUBROUTINE read_wfc(i_sp,id)
      USE system_data_types
      IMPLICIT NONE
      INTEGER i_sp,id
      INTEGER ir,nb,IOS
      character (len=80)dummy_var
      character (len=2)indx
      logical l_header_end
      !upf_oc=0._dp
      !upf_epseu=0._dp
      DO nb=1,upf_nwfc(I_SP)
      write(indx, '(i0)') nb
!      CALL search(iD, "CHI."//indx, .FALSE.)
      !Read (id, FMT='(A)')dummy_var
      l_header_end=.false.
            READ(id,iostat = ios,FMT='(A)')dummy_var!read_string
      IF(INDEX(dummy_var,"<PP_CHI."//indx).EQ.0 )THEN
            PRINT*,"PP_CHI.",INDX," NOT FOUND",dummy_var
            STOP
      ENDIF

       
      CALL HEADER_TEST(ID,"type",dummy_var,l_header_end)
      CALL HEADER_TEST(ID,"size",dummy_var,l_header_end)
      CALL HEADER_TEST(ID,"columns",dummy_var,l_header_end)
      CALL HEADER_TEST(ID,"index",dummy_var,l_header_end)
      CALL HEADER_TEST(ID,"occupation",dummy_var,l_header_end)
      !READ(dummy_var,*)upf_oc(nb,i_sp)
      CALL HEADER_TEST(ID,"pseudo_energy",dummy_var,l_header_end)
      !READ(dummy_var,*)upf_epseu(nb,i_sp)
      CALL HEADER_TEST(ID,"label",dummy_var,l_header_end)
      read (id, FMT='(A)')dummy_var
      
      if(INDEX(dummy_var,">").NE.0 ) l_HEADER_end=.true.
      if(l_header_end)read (id, *, err = 105, iostat =ios)(upf_chi(ir,i_sp,nb),ir=1,upf_meshv(i_sp))
      Read (id, FMT='(A)')dummy_var
      !CALL search_end(iD, "CHI."//INDX, .FALSE.)
      !do ir = upf_kbeta(nb,i_sp) + 1, upf_meshv(i_sp)
      !  upf_beta(ir, nb, i_sp) = 0.d0
      !enddo
      ENDDO
      return
      105 print*,"error"
      END SUBROUTINE read_wfc

      SUBROUTINE read_rhoatom(i_sp,id)
      USE system_data_types
      IMPLICIT NONE
      integer i_sp,id,ir,ios
      read (id,*,err=106,iostat=ios)(rho_atom(ir,i_sp), ir=1,upf_meshv(i_sp))
      return

      106 print* ,"ERROR!"
      END SUBROUTINE read_rhoatom
      
      SUBROUTINE convert_upf_to_psp(i_sp)
      USE system_data_types
      IMPLICIT NONE
      integer i_sp,l,ir,i,iv
      rr(1:upf_meshv(i_sp),i_sp)=upf_r(1:upf_meshv(i_sp),i_sp)!change name rr to rgrid
      rw(1:upf_meshv(i_sp),i_sp)=upf_rab(1:upf_meshv(i_sp),i_sp)
      !upf%vloc(:) = vnl(1:upf%mesh,upf%lloc)*e2
      vr(1:upf_meshv(I_SP),i_sp,1)=upf_vloc(1:upf_meshv(i_sp),i_sp)*0.5
      !CALL DSCAL(upf_meshV(I_SP),0.5D0,VR(1,I_SP,1),1)
      !upf%chi(:,:) = chi(1:upf%mesh,1:upf%nwfc)
      RPS(1:upf_meshv(I_SP),i_sp,1:upf_nwfc(i_sp))=upf_chi(1:upf_meshv(i_sp),i_sp,1:upf_nwfc(i_sp))
      RETURN
      END SUBROUTINE convert_upf_to_psp

      END MODULE READUPF
