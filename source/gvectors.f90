MODULE gvectors
  CONTAINS
!-----------------------------------------------------------------------------------------!
  SUBROUTINE set_ngvectors
!PURPOSE: Find the number of gvectors for the density and wfn cutoffs
    USE kinds
    USE system_data_types!, ONLY :  ngpw,ngrho, &
                         !          nrgrids,b_cell,gcut_wfn,gcut_rho,gamma_point,&
                         !          ionode
    USE math, ONLY : sort_array
    IMPLICIT NONE

!    REAL(KIND=dp), ALLOCATABLE :: gsq_mat(:)
    INTEGER :: nr1,nr2,nr3,ix,iy,iz,ymin,ymax,zmin,zmax
    REAL(KIND=DP)  :: gsq
!Construction of g-vectors
    nr1=nrgrids(1) 
    nr2=nrgrids(2) 
    nr3=nrgrids(3)
!    ALLOCATE(gsq_mat((nr1-1)*(nr2-1)*(nr3-1)/2))

!Ref. CPMD, numpw.F
   !print *, 'nr1 =',nr1,' nr2 =',nr2, 'nr 3=',nr3, 'gcut_wfn',gcut_wfn,'gcut_rho',gcut_rho
!open(777,file='check_g2.dat')
    ngpw=0
    ngrho=0
    DO ix=0,nr1-1
      ymin=-nr2+1
      ymax=nr2-1  
      IF(ix==0)ymin=0
      DO iy=ymin,ymax
        zmin=-nr3+1
        zmax=nr3-1
        IF(ix==0.AND.iy==0)zmin=0
        DO iz=zmin,zmax
           gsq=Gsquare(ix,iy,iz,b_cell)
           IF(gsq<=Gcutoff%pw) ngpw=ngpw+1
           IF(gsq<=Gcutoff%rho) ngrho=ngrho+1
        END DO
      END DO 
    END DO
!    IF(ionode)THEN
!     print *, '#pw for density cutoff              =',ngrho
!     print *, '#pw for wfn cutoff                  =',ngpw
!    END IF
!   close(777)
   
  END SUBROUTINE set_ngvectors

!-----------------------------------------------------------------------------------------!
  SUBROUTINE set_reciprocal !Based on rggen.F CPMD

    USE kinds
    USE system_data_types, ONLY :  nrgrids,nrlead,nnr1
    IMPLICIT NONE
    INTEGER :: nr1,nr2,nr3,nh1,nh2,nh3
  
    ALLOCATE(nrlead(3)) !kr1,kr2,kr3

    !compute nrlead arrays (number of grids - leading dimensions)
    CALL get_lead_dim_nrgrids(nrgrids,nrlead) !TODO c

!Setting variables as in rggen.F
    nnr1=nrlead(1)*nrlead(2)*nrlead(3) !Compute total number of grids

!    print *,'reciprocal|  ngrids =',nrgrids(1:3)
!    print *,'reciprocal|  nrlead =',nrlead(1:3)
!    print *,'reciprocal|  nnr1   =',nnr1

!Parallelize x-grids over processors
    CALL distribute_YZplanes(nrgrids(1)) !TODO: check if we should pass nrlead(1) or nrgrids(1)

!Parallelize g-vectors over processors (as rays) ; G-vectors are packed here
    CALL distribute_gvectors

  END SUBROUTINE set_reciprocal


!-----------------------------------------------------------------------------------------!
  SUBROUTINE distribute_YZplanes(nr1) !Based on loapa.F CPMD
    USE kinds
    USE system_data_types, ONLY :  ncpu,nrxplane
    IMPLICIT NONE

    INTEGER :: nr1
    REAL(KIND=dp) :: nplanes,xp_min,xp_max,chunk
    INTEGER :: i

!Disribute YZ planes over processors; NR1 (or x-grids) is (are) distributed over the processors
!Sets the array nrxplane(processor_index,1) starting index of the X grid on processor_index (=me+1)
!               nrxplane(processor_index,2) ending   index of the X grid on processor_index (=me+1)
    ALLOCATE(nrxplane(ncpu,2))
    nrxplane(:,:)=0

    nplanes=DBLE(nr1)
    chunk=nplanes/ncpu
    xp_min=0.0D0
    DO i=ncpu,1,-1
      xp_max=xp_min+chunk
      nrxplane(i,1)=NINT(xp_min)+1 !Starting Index
      nrxplane(i,2)=NINT(xp_max)   !Ending Indiex
      IF(NINT(xp_max)>nr1.OR.i==1) nrxplane(i,2)=nr1 !Enforcing the maximum value
      xp_min = xp_max
    ENDDO
 !   print *, 'distribute_YZplane | nrxplane(1:ncpu,1) =',nrxplane(1:ncpu,1)
 !   print *, 'distribute_YZplanes| nrxplane(1:ncpu,2) =',nrxplane(1:ncpu,2)

  END SUBROUTINE distribute_YZplanes
!-----------------------------------------------------------------------------------------!

  SUBROUTINE distribute_gvectors() !Based on loapa.F CPMD
    USE kinds
    USE system_data_types!, ONLY :  ncpu,nrgrids,nrlead,nrlead_l,ngpw,ngrho,ngpw_l,ngrho_l,&
                         !          nrgrids_l,hg,inyh,nrxplane,Ggridinfo,b_cell,&
                         !          gcut_rho,gcut_wfn,icpu,ionode,nnr1
    USE mympi, ONLY : MPI_GlobSumI2
!
    IMPLICIT NONE
    INTEGER, POINTER :: ixray(:,:),ihray(:,:)
    INTEGER          :: kr2s,kr3s,ixrays,nr1,nr2,nr3,nh1,nh2,nh3,&
                        ix,iy,iz,ii,ip,ipp,icpu1,in1,in2,in3,ngrays,nhrays,ig,ihrays,img,&
                        ymin,ymax,zmin,zmax
    INTEGER,       DIMENSION(:,:), POINTER :: mgpa
    REAL(kind=dp) :: sign,gsq
    REAL(kind=dp), PARAMETER :: EPSGX=1.D-12

!Below sets the variable mapping for compatibility with CPMD routines <TODO:Change>
    nr1=nrgrids(1) ; nr2=nrgrids(2) ; nr3=nrgrids(3)
    nh1=nr1/2+1    ; nh2=nr2/2+1    ; nh3=nr3/2+1
    kr2s=nrlead(2) ; kr3s=nrlead(3) 

  !  print *, 'distribute_gvectors| kr2s=',nrlead(2)
  !  print *, 'distribute_gvectors| kr3s=',nrlead(3)
!.....................
    ALLOCATE(ixray(kr2s,kr3s)) ; ALLOCATE(ihray(kr2s,kr3s))

!Compute ixray(iy,iz) -> number of G vectors within a cutoff for (iy,iz), and any value of ix
    CALL xfft(ixray,Gcutoff%pw) !Number of G vectors within GCUT_WFN cutoff in the YZ plane
    CALL xfft(ihray,Gcutoff%rho) !-do- within the GCUT_DEN cutoff
!
!   How many rays along X have non-zero G component? ixrays store this information
    ixrays = 0
    DO iz=1,kr3s
      DO iy=1,kr2s
         IF(ixray(iy,iz)/=0)ixrays=ixrays+1
      ENDDO
    ENDDO
   ! print *, 'ixrays =',ixrays

!   Assign each ray to a processor based on the length of the non-zero G-components
!   Rays with large number of non-zero G components are assigned to different IP. 
!   Assignment to processors are done like this (in the case of 2 processors, for e.g.): 1 ; 2 ; 2 ; 1 
!                                               (in the increasing order of number of compoennts)
    ipp=0
    DO ii=1,ixrays
      ipp=ipp+1
      ip=MOD(ipp,2*ncpu)
      IF(ip==0)ip=2*ncpu
      IF(ip>ncpu)ip=2*ncpu+1-ip 
      CALL iraymax(kr2s,kr3s,ixray,iy,iz) !Find the maximum G among the non-assigned ray
      !Overwrite the x-ray arrays with IPs assigned to it
      if(iy<0.or.iz<0.or.iy>kr2s.or.iz>kr3s)STOP 'error in iraymax'
      ixray(iy,iz)=-ip  !Negative number is given so that the ray with maximum G components can be easily identified 
    ENDDO

    call cxfft(ixray,Gcutoff%pw) !Check the IP assignment and and reassign the missing ones (either -G & G)


!.....................

!   Compute x-rays and assign to processors based on density cutoff
    ihrays=0
    DO iz=1,kr3s
      DO iy=1,kr2s
        IF(ixray(iy,iz)/=0)ihray(iy,iz)=ixray(iy,iz) !Use assignments based on wfn cutoff 
        IF(ihray(iy,iz)>0)ihrays=ihrays+1        !Count number of x-rays based on density cutoff
      END DO
    END DO
    !print *, 'ihrays =',ihrays
    DO ii=1,ihrays
      ipp=ipp+1
      ip=MOD(ipp,2*ncpu)
      IF(ip==0)ip=2*ncpu
      IF(ip>ncpu)ip=2*ncpu+1-ip
      CALL iraymax(kr2s,kr3s,ihray,iy,iz)   ! Assign processor ID based on the number of G-vectors (considering density cutoff)
      ihray(iy,iz)=-ip
    END DO

    CALL cxfft(ihray,Gcutoff%rho)  !Check the IP assignments and reassign the missing ones (either -G & G)
!......................

!   Local number of g-vectors (for allocation, prior to its exact calculation, if required)
!   Exact calculation of local g-vectors is done below
    ngpw_l=INT(DBLE(ngpw)/DBLE(ncpu)*1.2d0)
    ngpw_l=MIN(ngpw,ngpw_l)
    ngrho_l=INT(DBLE(ngrho)/DBLE(ncpu)*1.2d0)
    ngrho_l=MIN(ngrho,ngrho_l)

    !print *, 'ngpw_l (upper est.)=',ngpw_l
    !print *, 'ngrho_l (upper est.)=',ngrho_l

    ALLOCATE(hg(ngrho_l))      ! Stores G^2 value of local G vectors
    ALLOCATE(inyh(3,ngrho_l)) ! Grid-indices of a local Gvector if G^2 < G_cut  (Local G vector means, G-vector part of a processor after the allocation processor wise)

!    PRINT *, "ngpw_l=",npgw_l, " ngrho_l=",ngrho_l
    

!   Reorder G-vectors

!   Calculate the G-vectors on each processors (after the above allocation of rays to processors)
!   We use the processor ID information in ixray or ihray which were computed above
    hg(:)=0.d0
    inyh(:,:)=0
    ngrho_l=0 !Recalculate (exactly) ngrho_l (processor dependent)
    ngpw_l=0  !Recalculate (exactly) ngpw_l  (processor dependent)
    !print *, 'icpu =', icpu
    ig=0
    DO ix=0,nr1-1
      !print *, 'loop =',ix
      ymin=-nr2+1
      ymax=nr2-1  
      IF(ix==0)ymin=0
      DO iy=ymin,ymax
        zmin=-nr3+1
        zmax=nr3-1
        IF(ix==0.AND.iy==0)zmin=0
        DO iz=zmin,zmax
!
          gsq=Gsquare(ix,iy,iz,b_cell)
!
          IF(gsq<Gcutoff%rho)THEN
            ig=ig+1
            in1=nh1+ix
            in2=nh2+iy
            in3=nh3+iz
            IF(gsq<Gcutoff%pw)THEN
              icpu1=ixray(in2,in3)
              !Count only if the G vector is part of the X-ray of this processor
              IF(-icpu1==icpu) ngpw_l=ngpw_l+1 
              sign=-1.d0
            ELSE
              sign=1.d0
            END IF 
            icpu1=ihray(in2,in3)
            !Count only if the G vector is part of the X-ray of this processor
            IF(-icpu1==icpu)THEN
              ngrho_l=ngrho_l+1
              hg(ngrho_l)=gsq*(1.d0+DSQRT(DFLOAT(ig-1))*EPSGX*sign)
              inyh(1,ngrho_l)=in1 !Storing the grid index of a G-vector belonging to a processor if G^2<Gcut
              inyh(2,ngrho_l)=in2
              inyh(3,ngrho_l)=in3
            END IF
          END IF
        END DO
      END DO
    END DO

    !print *, 'ngrho_l', ngrho_l,'ngpw_l=',ngpw_l
     
!Compute the number of G components of a ray (on each processor)
!and store the the ray indices on each processor 
   ALLOCATE(mgpa(kr2s,kr3s)) !Ray Index array
   mgpa=0
   ngrays=0
   nhrays=0
   DO iz=1,kr3s
     DO iy=1,kr2s
       IF(ixray(iy,iz)==-icpu)THEN
         ngrays=ngrays+1 !Total #non zero rays in icpu  (based on pw cutoff)
         mgpa(iy,iz)=ngrays  !Store the ray index (corresponding to iy,iz)
       END IF
       IF(ihray(iy,iz)==-icpu)THEN
         nhrays=nhrays+1 !Total #non-zero rays in icpu (based on density cutoff)
       END IF
     END DO
   END DO
   img=ngrays
   DO iz=1,kr3s
     DO iy=1,kr2s
        IF(ihray(iy,iz)==-icpu.and.mgpa(iy,iz)==0)THEN
          img=img+1  ! Index of the ray is augmented if there is a ray in this processor belong to density cutoff
          mgpa(iy,iz)=img !Store the ray index (based on denstiy cutoff)
        END IF
     END DO
   END DO
   !print *, 'ngrays=',ngrays,'nhrays=',nhrays

!Modifying the local value of the nrgrids (due to parallelization over x grids)
   ALLOCATE(nrgrids_l(3))
   nrgrids_l(1)=nrxplane(icpu,2)-nrxplane(icpu,1)+1 !xmax-xmin+1 = #x-grids in the current processor 
   nrgrids_l(2:3)=nrgrids(2:3)

   ALLOCATE(nrlead_l(3))
   CALL get_lead_dim_nrgrids(nrgrids_l,nrlead_l) !TODO check if this is OK

   nnr1=nrlead_l(1)*nrlead_l(2)*nrlead_l(3)
   
   ALLOCATE(Ggridinfo(8,ncpu))
   Ggridinfo=0
!Copy local G vector / grid information and communicate to all
   Ggridinfo(1,icpu)=ngpw_l !local value of number of planewaves (pw cutoff)
   Ggridinfo(2,icpu)=ngrho_l !local value of number of planewaves (den cutoff)
   Ggridinfo(3,icpu)=nrgrids_l(1) !local value of #x-grids
   Ggridinfo(4,icpu)=nrgrids_l(2) !local value of #y-grids
   Ggridinfo(5,icpu)=nrgrids_l(3) !local value of #z-grids
   Ggridinfo(6,icpu)=ngrays       !local value of #rays (pw cutoff)
   Ggridinfo(7,icpu)=nhrays       !local value of #rays (den cutoff)
   Ggridinfo(8,icpu)=nnr1
   !print *, 'Ggridinfo done..'
   
   CALL MPI_GlobSumI2(Ggridinfo,SIZE(Ggridinfo,1)*SIZE(Ggridinfo,2))

   !print *, 'MPI_GlosumI2  done..'

   IF(ionode)THEN
    WRITE(*,*)
    WRITE(*,'(A21,A50)')"#Plane waves","Local  Dimensions"
    WRITE(*,'(A6,A8,A8,3A10,2A12,A16)') 'CPU-ID','PW-CUT','DEN-CUT','#X-GRIDS','#Y-GRIDS','#Z-GRIDS','#RAYS(PW)','#RAYS(DEN)'&
             ,'#GRIDS(G-SPACE)'
     DO icpu1=1,ncpu
     WRITE(*,'(I6,I8,I8,3I10,2I12,I16)') icpu1,Ggridinfo(1:8,icpu1)
       
     END DO
     WRITE(*,*)
   END IF


   DEALLOCATE(ixray,ihray,mgpa)


   CALL sort_gvectors(hg,inyh,ngrho_l)
   ! print *, ' to-check-order gvectors ', icpu
   CALL check_order_gvectors
    !print *, ' to-assign ipg0 ', icpu
   CALL assign_ipg0 !Identify the process ID with G=0 component 
    !print *, ' done-assign ipg0 ', icpu
   !CALL aliasing
   CALL set_FFT_1d_mapping_indices !Set mapping indices for +G/-G componenets from 3D mesh to 1D mesh
   !Store gvector indices and G^2
   CALL compute_gvec_gvec2
  END SUBROUTINE distribute_gvectors
!-----------------------------------------------------------------------------------------!
  !Ref. xfft() within loadpa.F in CPMD
  !Counting the rays
  SUBROUTINE xfft(iray,gcut)
    USE kinds
    USE system_data_types, ONLY :  ncpu,nrlead,nrgrids,b_cell
    IMPLICIT NONE

    INTEGER,  INTENT(OUT) :: iray(:,:)
    REAL(KIND=DP), INTENT(IN)            :: gcut

    INTEGER :: nr1,nr2,nr3,nh1,nh2,nh3,&
               ix,iy,iz,xmax,xmin,ymax,ymin,zmax,zmin,&
               in2,in3,id2,id3,ir1,ir2
    REAL(KIND=DP)  :: gsq

    !print *, 'xfft| gcut=',gcut
    iray(:,:)=0
    nr1=nrgrids(1)
    nr2=nrgrids(2)
    nr3=nrgrids(3)

    nh1=nr1/2+1
    nh2=nr2/2+1
    nh3=nr3/2+1

    DO ix=0,nr1-1
      ymin=-nr2+1
      ymax=nr2-1  
      IF(ix==0)ymin=0
      DO iy=ymin,ymax
        zmin=-nr3+1
        zmax=nr3-1
        IF(ix==0.AND.iy==0)zmin=0
        DO iz=zmin,zmax
!
          gsq=Gsquare(ix,iy,iz,b_cell)
!
          IF(gsq<gcut)THEN
            in2=nh2+iy
            in3=nh3+iz
            id2=2*nh2-in2
            id3=2*nh3-in3
            ir1=iray(in2,in3) 
            ir2=iray(id2,id3)
!Count the rays (nr/2+1 to nr)
!iray(id2,id3)=iray(id2,id3)+1  !Count the rays (1 to nr/2)
            IF(ir1>0)THEN
              iray(in2,in3)=ir1+1
              IF(ir2>0)THEN
                IF((in2/=id2).or.(in3/=id3))&
                  STOP 'inconsistent mesh| stop in xfft'
              END IF
            ELSEIF(ir2>0)THEN
              iray(id2,id3)=ir2+1
            ELSE
              iray(in2,in3)=1
            END IF
          ENDIF
        ENDDO
      ENDDO
    ENDDO
!  OPEN(800,FILE='iray.log')
!  DO iy=1,nr2
!  DO iz=1,nr3
!    write(800,*)iy,iz,iray(iy,iz)
!  END DO
!  END DO
!  stop
  END SUBROUTINE xfft 

!-----------------------------------------------------------------------------------------!
  SUBROUTINE cxfft(iray,gcut)
    USE kinds
    USE system_data_types, ONLY :  ncpu,nrlead,nrgrids,b_cell
    IMPLICIT NONE

    INTEGER, INTENT(INOUT) :: iray(:,:)
    REAL(kind=dp)  :: gcut

    INTEGER :: nh1,nh2,nh3,nr1,nr2,nr3,&
               ix,iy,iz,xmin,xmax,ymin,ymax,zmin,zmax,&
               in1,in2,in3,id1,id2,id3,icpu1,icpu2

    REAL(kind=dp) :: gsq

    nr1=nrgrids(1)
    nr2=nrgrids(2)
    nr3=nrgrids(3)

    nh1=nr1/2+1
    nh2=nr2/2+1
    nh3=nr3/2+1

    DO ix=0,nr1-1
      ymin=-nr2+1
      ymax=nr2-1  
      IF(ix==0)ymin=0
      DO iy=ymin,ymax
        zmin=-nr3+1
        zmax=nr3-1
        IF(ix==0.AND.iy==0)zmin=0
        DO iz=zmin,zmax
!
          gsq=Gsquare(ix,iy,iz,b_cell)
!
          IF(gsq<gcut)THEN
            in2=nh2+iy
            in3=nh3+iz
            id2=2*nh2-in2
            id3=2*nh3-in3
            icpu1=iray(in2,in3)
            icpu2=iray(id2,id3)
            IF(icpu2==0) THEN
               iray(id2,id3)=icpu1
            ELSEIF(icpu1==0) THEN
               iray(in2,in3)=icpu2
            ELSEIF(icpu1/=icpu2) THEN
               print *,   'ERROR IN SETTING XRAY FIELDS',icpu1,icpu2
               STOP
            ENDIF
          ENDIF
        END DO
      END DO 
    END DO
  END SUBROUTINE cxfft
!-----------------------------------------------------------------------------------------!
  SUBROUTINE iraymax(n1,n2,ir,i1,i2) !based on loadpa.F CPMD
!     Find the entry ir(i1,i2) which is the highest value
!     It outputs the indices i1, and i2
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: n1,n2,ir(n1,n2)   !Input
      INTEGER, INTENT(OUT) :: i1,i2             !Output
      INTEGER :: im,j1,j2

      im=-2**30 ! Large negative number
      i1=0 ; i2=0
      DO j2=1,n2
        DO j1=1,n1
          IF(ir(j1,j2)>=im) THEN
            im=ir(j1,j2)
            i1=j1 ; i2=j2
          ENDIF
        ENDDO
      ENDDO
  END SUBROUTINE iraymax



!Reference: LEADIM subroutine of CPMD
  SUBROUTINE get_lead_dim_nrgrids(nrgrids,nrlead)
    IMPLICIT NONE
    INTEGER nrgrids(3), nrlead(3)

    nrlead(:)=nrgrids(:)+MOD(nrgrids(:)+1,2) !If nrgrids is even, then nrlead will be nrgrids+1 ; else the same
  END SUBROUTINE 

!Returns the value of G^2 for a given grid index (ix,iy,iz)
  FUNCTION Gsquare(ix,iy,iz,b_cell) RESULT(gsq)
    USE KINDS
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ix,iy,iz
    REAL(kind=dp), INTENT(IN)  :: b_cell(3,3)
    REAL(kind=dp)             :: gsq

    gsq= (DBLE(ix)*b_cell(1,1)+DBLE(iy)*b_cell(1,2)+DBLE(iz)*b_cell(1,3))**2 + &
         (DBLE(ix)*b_cell(2,1)+DBLE(iy)*b_cell(2,2)+DBLE(iz)*b_cell(2,3))**2 + &
         (DBLE(ix)*b_cell(3,1)+DBLE(iy)*b_cell(3,2)+DBLE(iz)*b_cell(3,3))**2

  END FUNCTION Gsquare


  SUBROUTINE sort_gvectors(hg,inyh,ngrho_l)
  USE math, ONLY : sort_array2
  USE kinds
  IMPLICIT NONE
!
  INTEGER          :: ngrho_l
  INTEGER, POINTER :: inyh(:,:)
  REAL(KIND=DP), POINTER  :: hg(:)
!
  INTEGER, ALLOCATABLE :: g_index(:)
  INTEGER ig, jg, ind(3), iii
!
!     reorder G vectors in the increasing value 
  ALLOCATE(g_index(ngrho_l))
  CALL sort_array2(ngrho_l,hg,g_index)
!  open(888,file='debug3.dat')
!  do iii=1,ngrho_l
!    write(888,*)g_index(iii)
!  end do
!  close(888)

! Rearrange the inyh array 
  sorting_loop1: DO ig=1,ngrho_l-1
     iii=0
     sorting_loop2: DO
       iii=iii+1
       !if(iii>99999)STOP 'error in sorting'
       IF(g_index(ig)/=ig) THEN 
!       map g_index in ascending order
        jg=g_index(ig)
        g_index(ig)=g_index(jg)
        g_index(jg)=jg

        ind(:)=inyh(:,g_index(ig)) !current inyh indices
        inyh(:,g_index(ig))=inyh(:,g_index(jg)) !map with correct indices
        inyh(:,jg)=ind(:) !swap

       ELSE
         EXIT sorting_loop2
       END IF
     END DO sorting_loop2
  END DO sorting_loop1
!     DO loop ends here
  DEALLOCATE(g_index)
!  open(888,FILE='debug.dat')
!  do iii=1,ngrho_l
!    write(888,*)inyh(1:3,iii)
!  end do
!  close(888)
!  stop 'stoopping debugging'
  END SUBROUTINE sort_gvectors 

  SUBROUTINE check_order_gvectors
!Sort across the the processors
    USE kinds
    USE system_data_types, ONLY :  ngrho_l,mapgp,ngrho,ncpu,icpu,hg,&
                                   Ggridinfo,ngpw_l,ngpw
    USE mympi,ONLY: mpi_shift, mpi_globsumr2s
    IMPLICIT NONE

    INTEGER :: ig,ngrho_l_max,ip,l,mho,i,j,ipp,ngrho_ip
    REAL(KIND=dp), ALLOCATABLE :: ho(:),hx(:)
    REAL(KIND=dp)  :: xt,xm
    logical :: located

    ALLOCATE(mapgp(ngrho_l))
    DO ig=1,ngrho_l
      mapgp(ig)=ig
    END DO

    ngrho_l_max=0
    DO ip=1,ncpu
      ngrho_l_max=MAX(ngrho_l_max,Ggridinfo(2,ip)) !Find the maximum number of ngrho_l across the processors
    END DO
    !print *, 'ngrho_l=',ngrho_l,'ngrho_l_max=',ngrho_l_max

    ALLOCATE(ho(ngrho_l_max)) ; ALLOCATE(hx(ngrho_l_max))
    ho(1:ngrho_l_max)=0.d0 ; hx(1:ngrho_l_max)=0.d0

    hx(1:ngrho_l)=hg(1:ngrho_l) !TODO use dcopy 
!The G^2 matrix of the current IP is stored in HX.
!For any processor ip, the G^2 value is fetched in HO through mpi_shift
!
!Thereafter, we check if any G vector has to be send to another processor

    DO ip=1,ncpu-1
      ipp=MOD(icpu+ip-1,ncpu)
      !print *, 'ip =', ip, ' ipp=', ipp
      ngrho_ip=Ggridinfo(2,ipp+1) 
      CALL MPI_shift(hx,ho,ngrho_l_max*8,icpu-1,-ip)
      l=1
      located=.false.
      DO i=1,ngrho_l
         IF(hg(i)<ho(1))THEN !If the first G^2 is high, then the current G can't be in that ip
            mho=0
         ELSE
           !It is possible that the current G can be in ip. Need to identify the
           !location
           !
           !Setting Maximum possible values of l and mho
           !or give avalue between 1 and nghrho_ip-1 if a suitable location is
           !is located with G^2(i) < G^2(new_location) < G^2(i+1)
           do_locate: DO j=l,ngrho_ip-1 !Loop over all the G vectors
             IF(hg(i)>ho(j).and.hg(i)<ho(j+1))THEN
               l=j
               mho=j
               located=.true.
               EXIT do_locate
             END IF
           END DO do_locate
           IF(.not.located)THEN
             l=ngrho_ip 
             mho=ngrho_ip
           END IF
         END IF
         mapgp(i)=mapgp(i)+mho
      END DO
    END DO
    !print *, 'done mapgp ', icpu

!   IF(icpu.eq.1)THEN
!      open(888,FILE='debug5-1.dat')
!      print *, 'ngrho_l, before priting mapgp',ngrho_l
!      DO i=1,ngrho_l
!       WRITE(888,*)i,mapgp(i)
!      END DO
!      CLOSE(888)
!   ELSEIF(icpu.eq.2)THEN
!      open(889,FILE='debug5-2.dat')
!      print *, 'ngrho_l, before priting mapgp',ngrho_l
!      DO i=1,ngrho_l
!       WRITE(889,*)i,mapgp(i)
!      END DO
!      CLOSE(889)
!   END IF

!Checking the ordering
   xm=0.d0
    DO i=1,ngpw_l
      xm=xm+mapgp(i)
    END DO
    !Summing across the processors
  !  print *, 'me (before sum)= ', xm, icpu
 !   print *, 'to glosumr2s ', icpu
    CALL MPI_GlobSumR2s(xm)
 !   print *, 'done glosumr2s ', icpu
  !  print *, 'me (after sum)= ', xm, icpu
  !  print *, 'ngrho ', ngpw,icpu
    xt=DBLE(ngpw)
    xt=0.5d0*xt*(xt+1.d0)-xm
    IF(ABS(xt).gt.0.1d0)THEN
        PRINT *, 'warning! xt > 0.1 ; xt and xm=', xt,xm
    ELSE
        !PRINT *, 'Gvector pointer mapping list mapgp is tested OK'
    END IF
    DEALLOCATE(ho)
    DEALLOCATE(hx)
    
    !print *, 'done dealloc ', icpu
  END SUBROUTINE check_order_gvectors


  SUBROUTINE  assign_ipg0 !Identify the process ID with G=0 component 
!
  USE kinds
  USE system_data_types, ONLY : g0_stat,g0_ip,ngpw_l,ionode,hg,icpu
  USE mympi, ONLY: MPI_GlobsumI2s
  IMPLICIT NONE
  REAL(KIND=dp), PARAMETER :: myzero=1.D-5
  INTEGER :: ig
!
  g0_ip=0
  g0_stat=.FALSE.

  DO ig=1,ngpw_l
    IF(hg(ig)<myzero)THEN
      g0_stat=.TRUE.
      g0_ip=icpu
    END IF
  END DO
  CALL MPI_GlobsumI2s(g0_ip)
  !IF(ionode) write(6,'(a,i5)') '  G=0 component on processor = ',g0_ip
!
  END SUBROUTINE assign_ipg0
 
  
  SUBROUTINE set_FFT_1d_mapping_indices()
    USE kinds
    USE system_data_types, ONLY : map_grid1d_p, &
                                  map_grid1d_m,nrlead,nrgrids,inyh,ngrho_l,icpu
  
    IMPLICIT NONE
    INTEGER :: i1,i2,i3,kr1s,nr1,nr2,nr3,nh1,nh2,nh3,ig

    !print *, ' in set_fft_1d...', icpu
    ALLOCATE(map_grid1d_p(ngrho_l))
    ALLOCATE(map_grid1d_m(ngrho_l))
    nr1=nrlead(1)
    nr2=nrlead(2)
    nr3=nrlead(3)
    nh1=nrgrids(1)/2+1
    nh2=nrgrids(2)/2+1
    nh3=nrgrids(3)/2+1
    DO ig=1,ngrho_l
      i1=inyh(1,ig)    !+G part
      i2=inyh(2,ig)
      i3=inyh(3,ig)
      map_grid1d_p(ig)=i1+(i2-1)*nr1+(i3-1)*nr2*nr3 !linear mapping of 3D grid
      i1=-i1+2*nh1 !-G part
      i2=-i2+2*nh2
      i3=-i3+2*nh3
      map_grid1d_m(ig)=i1+(i2-1)*nr1+(i3-1)*nr2*nr3
    ENDDO
      !if(icpu.eq.1)THEN
       ! open(880,file='debug6.dat')
        !do ig=1,ngrho_l
          !write(880,*)map_grid1d_p(ig),map_grid1d_m(ig)
        !end do
        !close(880)
      !elseif(icpu.eq.2)THEN
        !open(881,file='debug7.dat')
        !do ig=1,ngrho_l
          !write(881,*)map_grid1d_p(ig),map_grid1d_m(ig)
        !end do
        !close(881)
      !end if
  END SUBROUTINE set_FFT_1d_mapping_indices

  SUBROUTINE compute_gvec_gvec2

    USE kinds
    USE system_data_types, ONLY : nrgrids,inyh,ngrho_l,gvec,hg,b_cell

    IMPLICIT NONE
    INTEGER :: ix,iy,iz,ig,nh1,nh2,nh3
    REAL(KIND=DP) T1,T2,T3

    nh1=nrgrids(1)/2+1
    nh2=nrgrids(2)/2+1
    nh3=nrgrids(3)/2+1

    IF(.NOT.ASSOCIATED(gvec)) ALLOCATE(gvec(3,ngrho_l))
    IF(.NOT.ASSOCIATED(hg)) ALLOCATE(hg(ngrho_l))
    IF(.NOT.ASSOCIATED(inyh)) STOP 'error| compute_gvec_gvec2 | inyh not allocated'

    DO ig=1,ngrho_l

        ix=inyh(1,ig)-nh1
        iy=inyh(2,ig)-nh2
        iz=inyh(3,ig)-nh3

        gvec(1,ig)=(DBLE(ix)*b_cell(1,1)+DBLE(iy)*b_cell(1,2)+DBLE(iz)*b_cell(1,3)) 
        gvec(2,ig)=(DBLE(ix)*b_cell(2,1)+DBLE(iy)*b_cell(2,2)+DBLE(iz)*b_cell(2,3))
        gvec(3,ig)=(DBLE(ix)*b_cell(3,1)+DBLE(iy)*b_cell(3,2)+DBLE(iz)*b_cell(3,3))

        !t1=(DBLE(ix)*b_cell(1,1)+DBLE(iy)*b_cell(2,1)+DBLE(iz)*b_cell(3,1))
        !t2=(DBLE(ix)*b_cell(1,2)+DBLE(iy)*b_cell(2,2)+DBLE(iz)*b_cell(3,2))
        !t3=(DBLE(ix)*b_cell(1,3)+DBLE(iy)*b_cell(2,3)+DBLE(iz)*b_cell(3,3))

        !hg(ig)=t1*t1+t2*t2+t3*t3
         hg(ig)=DOT_PRODUCT(gvec(1:3,ig),gvec(1:3,ig))
    ENDDO
  END SUBROUTINE compute_gvec_gvec2
  SUBROUTINE DISTRIBUTE_ORB
  !     DISTRIBUTE ORBITALS
  !   ==--------------------------------------------------------------==
  USE KINDS
  USE SYSTEM_DATA_TYPES
  REAL(KIND=DP)    XSTATES,XSNOW,XSAIM
  ALLOCATE(NST12(NCPU,2))
      XSTATES=DBLE(NSTATE)
      XSNOW=0.0D0
      DO I=NCPU,1,-1
        XSAIM = XSNOW + XSTATES/NCPU
        NST12(I,1)=NINT(XSNOW)+1
        NST12(I,2)=NINT(XSAIM)
        IF(NINT(XSAIM).GT.NSTATE) THEN
           NST12(I,2)=NSTATE
        ENDIF
        IF(I.EQ.1) THEN
           NST12(I,2)=NSTATE
        ENDIF
        XSNOW = XSAIM
      ENDDO
      NORBPE=NST12(icpu,2)-NST12(icpu,1)+1
      !IF(IONODE)Print*,"NDFNL=",NORBPE
 !write(6,*)"#8888888",NST12(icpu,2),NST12(icpu,1),ICPU
  END SUBROUTINE DISTRIBUTE_ORB
  SUBROUTINE DISTRIBUTE_ATOM
  !     DISTRIBUTE ATOMS
  !   ==--------------------------------------------------------------==
  USE KINDS
  USE SYSTEM_DATA_TYPES
  REAL(KIND=DP)     XSTATES,XSNOW,XSAIM
  INTEGER  IPEPT(2,NCPU)!IATPE(atom_t)
  ALLOCATE(IATPE(atom_t))
      XSTATES=DBLE(atom_t)
      XSNOW=0.0D0
      DO I=NCPU,1,-1
        XSAIM = XSNOW + XSTATES/NCPU
        IPEPT(1,I)=NINT(XSNOW)+1
        IPEPT(2,I)=NINT(XSAIM)
        IF(NINT(XSAIM).GT.atom_t) THEN
           IPEPT(2,I)=atom_t
        ENDIF
        IF(I.EQ.1) THEN
           IPEPT(2,I)=atom_t
        ENDIF
        XSNOW = XSAIM
       !write(6,*)"#*****",i, IPEPT(1,I), IPEPT(2,I),ICPU
      ENDDO
      DO I=1,NCPU
        DO J=IPEPT(1,I),IPEPT(2,I)
          IATPE(J)=I
       ! write(6,*)"#*****",IATPE(J),J,i
        ENDDO
      ENDDO
      !write(6,*)"#*****",NCPU
  END SUBROUTINE DISTRIBUTE_ATOM

END MODULE gvectors


