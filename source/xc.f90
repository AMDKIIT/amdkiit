MODULE xc

  CONTAINS  
  SUBROUTINE cal_vpot_ex
  USE system_data_types
  USE constants
  USE gvectors
  USE kinds
  USE fft_interface
  use xc_f90_lib_m
  use mympi
  IMPLICIT NONE
  TYPE(xc_f90_func_t) :: xc_func
  TYPE(xc_f90_func_info_t) :: xc_info
  real(8) :: GRAD(nnr1,NLSD*4) 
  integer :: i, vmajor, vminor, vmicro, func_id 
  INTEGER ir,ig,nr1,nr2,nr3,n_func,ilsd
  REAL(kind=dp)vfac3,vfac2,vfac
  REAL(kind=dp)vxc,sxc
  REAL(kind=dp), DIMENSION(3)               :: sigma_, vsigma_
  REAL(kind=dp)::vsigma(max_fft),ex_tot(max_fft)
  REAL(kind=dp), DIMENSION(:,:), POINTER ::vpot_tot,vsigma_tot
  REAL(kind=dp) sc,scale, smoo,dsmoo,rhoa,rhob,rho_tmp, sx, texp, v1_, v1a_, v1b_, v1c, v1ca, v1cb, v1x, v2x, v2c, v2_,sgcx,sgcc
  REAL(kind=dp) v1xa, v1xb, v2a_, v2ab_, v2b_, v2ca, v2cab, v2cb,  v2xa, v2xab, v2xb
  COMPLEX(KIND=dp), DIMENSION(:), POINTER :: scr_v,scr_vtmp
  COMPLEX(KIND=dp), DIMENSION(:), POINTER :: WORK,vpotg,vpotrb,vpotr
  !COMPLEX(kind=dp), DIMENSION(:,:), POINTER :: VPTG
  REAL(kind=dp), DIMENSION(:), allocatable::EX,ro
  REAL(kind=dp), DIMENSION(2):: vpot
  REAL(kind=dp):: vpot_(nnr1)
  COMPLEX(KIND=DP) :: fa, fb, vg1, vg2
  REAL(kind=dp)vtmp(nnr1),vtmpb(nnr1),vtmpab(nnr1)
  REAL(kind=DP)fact,fact2,fact3
  integer*8 nn
  !LOGICAL is_lsd_but_not
 
  ALLOCATE(WORK(MAX_FFT),vpotg(ngrho_l))!,vptg(ngrho_l,nlsd))
  ALLOCATE(scr_v(MAX_FFT),scr_vtmp(ngrho),vpotr(MAX_FFT),vpotrb(MAX_FFT))!vpotr(nnr1),vpotrb(nnr1))
 
     scale=1.0_dp
    
     sgcx=0.0_dp
     sgcc=0.0_dp
     etxc=0.0_dp
!============================================================
if(nlsd==2) then
 DO ir=1,nnr1
             rho(ir,1)=rho(ir,1)-rho(ir,2)
 ENDDO
endif
!============================================================
CALL grad_cal(rho(:,1),grad(:,1))
if(nlsd==2)CALL grad_cal(rho(:,2),grad(:,5))
         fact = 1.0_dp
         fact2 = 1.0_dp
         fact3 = 2.0_dp
         IF( is_lsd_but_not ) THEN
            fact =0.5_dp
            fact2 = 0.25_dp
            fact3 = 1.0_dp
        ENDIF

! IF(.not.LOPEN_SHELL)THEN
IF(NLSD.EQ.1)THEN
     allocate(ex(max_fft),ro(max_fft),vsigma_tot(max_fft,1),vpot_tot(max_fft,1))
     do i = 1, nnr1
       ro(i) = rho(i,nlsd)
     end do
     ex=0.0
     vpot_=0.0
     vpot_tot=0.0
     ex_tot=0.0
     vsigma=0.0_dp
     vsigma_tot=0.0_dp
     DO Ir=1,3
       IF(func_string(ir).ne.'NONE')then
         func_id = xc_f90_functional_get_number ( func_string(ir) )
         call xc_f90_func_init(xc_func, func_id, XC_UNPOLARIZED)
         xc_info = xc_f90_func_get_info(xc_func)
         select case (xc_f90_func_info_get_family(xc_info))
         case(XC_FAMILY_LDA)
         call xc_f90_lda_exc_vxc(xc_func,nnr1, ro(1), ex(1),vpot_(1))
         case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
         call xc_f90_gga_exc_vxc(xc_func, nnr1, ro(1), grad(1,1), ex(1),vpot_(1),vsigma(1)) 
         vsigma_tot(:,1) = vsigma_tot(:,1) + scale * vsigma(:) *2.0_dp
         end select
          CALL xc_f90_func_end(xc_func)
         do i = 1, nnr1
           ex_tot(i) = ex_tot(i) + scale* ex(i) * ro(i)
           vpot_tot(i,1) = vpot_tot(i,1) +  scale*vpot_(i)
         end do

       ENDIF
     END DO
       
   do i = 1, nnr1
     !v_pot(i,nlsd)=v_pot(i,nlsd)+vpot_tot(i,1)
   end do
   !------------------------step 1
         DO i = 1,nnr1
            rho_tmp   = MAX(rho(i,nlsd),0.0_dp)
            IF (rho_tmp.GT.0.1_dp*gc_cutoff) THEN
               IF (rho_tmp.GT.4.0_dp*gc_cutoff) THEN
                  smoo  = 1.0_dp
                  dsmoo = 0.0_dp
               ELSE
                  texp=EXP(3.0_dp*(1.0_dp-rho_tmp/gc_cutoff))
                  smoo=1.0_dp/(1.0_dp+texp)
                  dsmoo=3.0_dp/gc_cutoff * texp*smoo*smoo
               ENDIF

               sx = ex_tot(i)
               sc = 0.0_dp
               v1x = vpot_tot(i,1)
               v1c = 0.0_dp
               v2x = vsigma_tot(i,1)
               v2c = 0.0_dp

               sgcx  = sgcx + smoo*sx
               sgcc  = sgcc + smoo*sc
               v1_    = dsmoo*(sx+sc) + smoo*(v1x+v1c)
               v2_    = smoo*(v2x+v2c)
            ELSE
               v1_    = 0.0_dp
               v2_    = 0.0_dp
            ENDIF
            !v1(i,1) = v1_
            vtmp(i) = v2_
            vpotr(i)=v1_
       !     vpotr(i),1 = vpotr(i,1) + CMPLX(0.0_dp,v2_*grad(i,2),kind=dp)
         ENDDO
   
    etxc=(sgcc+sgcx)!*parm%omega/REAL(nnrs,kind=real_8)

ELSEIF(NLSD.EQ.2)THEN
   allocate(ex(1),ro(2),vsigma_tot(max_fft,3),vpot_tot(max_fft,2))
   vsigma_tot(:,:)=0.0d0
   ex_tot=0.0d0
   vpot=0.0d0
   vpot_tot(:,:)=0.0d0
   vsigma_=0.0_dp
   !write(77,*)rho(1,1),rho(1,2),v_pot(1,1)
   DO Ir=1,3
   IF(func_string(ir).ne.'NONE')then
         func_id = xc_f90_functional_get_number ( func_string(ir) )
   do i=1,nnr1
         etxc=0.0
         nn=1
         ex(1)=0.0d0
         ro(1)=rho(i,1)
         if(is_lsd_but_not)then
                  ro(2)=0.0_dp!rho(i,2)
         else
             ro(2)=rho(i,2)
         endif

     call xc_f90_func_init(xc_func, func_id, XC_POLARIZED)
     xc_info = xc_f90_func_get_info(xc_func)
     select case (xc_f90_func_info_get_family(xc_info))
     case(XC_FAMILY_LDA)
     call xc_f90_lda_exc_vxc(xc_func,nn,fact*ro,ex(1),vpot(1))
     case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
          sigma_(1) = grad(i,1) ! aa
          !------------------------------------
            IF( is_lsd_but_not ) THEN
                     sigma_(3) = 0.0_dp
                     sigma_(2) = 0.0_dp
                     ro(2) = 0.0_dp
                  ELSE
                     sigma_(3) = grad(i,5) ! bb
                     sigma_(2) =(grad(i,2)*grad(i,6)+grad(i,3)*grad(i,7)+grad(i,4)*grad(i,8)) !ab
                  ENDIF

          !------------------------------------
     call xc_f90_gga_exc_vxc(xc_func,nn,fact*ro,fact2*sigma_,ex(1),vpot,vsigma_)
          vsigma_tot(i,1) = vsigma_tot(i,1) + scale * vsigma_(1) * fact3
          !vw vsigma_ = alpha/alpha
          vsigma_tot(i,2) = vsigma_tot(i,2) + scale * vsigma_(3) * fact3
          !vw vsigma_ = beta/beta
          vsigma_tot(i,3) = vsigma_tot(i,3) + scale * vsigma_(2) * fact3*0.5_dp !vw vsigma_ = alpha/beta
     end select
          ex_tot(i) = ex_tot(i) + scale* (ex(1)*(ro(1)+ro(2)))
          vpot_tot(i,:)=vpot_tot(i,:)+vpot(:)
     ENDDO
     ENDIF
     ENDDO
         
            DO i = 1,nnr1
            rhoa   = MAX(rho(i,1),0.0_dp)
            IF( is_lsd_but_not ) THEN
               rhob = 0.0_dp
            ELSE
               rhob   = MAX(rho(i,2),0.0_dp)
            ENDIF
            rho_tmp = rhoa + rhob
 
            !rho_tmp = MAX(rho(i,1),0.0_dp)+MAX(rho(i,2),0.0_dp)
            
            IF (rho_tmp.GT.0.1_dp*gc_cutoff) THEN
               IF (rho_tmp.GT.4.0_dp*gc_cutoff) THEN
                  smoo  = 1.0_dp
                  dsmoo = 0.0_dp
               ELSE
                  texp=EXP(3.0_dp*(1.0_dp-rho_tmp/gc_cutoff))
                  smoo=1.0_dp/(1.0_dp+texp)
                  dsmoo=3.0_dp/gc_cutoff * texp*smoo*smoo
               ENDIF

               sx = ex_tot(i)
               sc = 0.0_dp
               v1xa = vpot_tot(i,1)
               v1ca = 0.0_dp
               v1xb = vpot_tot(i,2)
               v1cb = 0.0_dp
               v2xa = vsigma_tot(i,1)
               v2ca = 0.0_dp
               v2xb = vsigma_tot(i,2)
               v2cb = 0.0_dp
               v2xab = vsigma_tot(i,3)
               v2cab = 0.0_dp

               sgcx  = sgcx + smoo*sx
               sgcc  = sgcc + smoo*sc

               v1a_     = dsmoo*(sx+sc) + smoo*(v1xa+v1ca)
               v2a_     = smoo*(v2xa+v2ca)
               v1b_     = dsmoo*(sx+sc) + smoo*(v1xb+v1cb)
               v2b_     = smoo*(v2xb+v2cb)
               v2ab_    = smoo*(v2xab+v2cab)

            ELSE
               v1a_     = 0.0_dp
               v1b_     = 0.0_dp
               v2a_     = 0.0_dp
               v2b_     = 0.0_dp
               v2ab_    = 0.0_dp
            ENDIF

            vpotr(i) = v1a_
            IF(.NOT.is_lsd_but_not)vpotrb(i) = v1b_
            vtmp(i) = v2a_
            IF(.NOT.is_lsd_but_not)vtmpb(i) = v2b_
            IF(.NOT.is_lsd_but_not)vtmpab(i) = v2ab_
        ENDDO

       etxc=(sgcc+sgcx)
ELSE
    WRITE(*,*)"!== -- ERROR IN NLSD VALUE -- ==!"
    STOP   
ENDIF
!-------------------------------------------------------------------------------------------------------------------
IF(NLSD==1)THEN
     WORK=(0._dp,0._dp)
      DO IG=1,NNR1
        WORK(IG)  = dcmplx( v_pot(IG,1),0.D0)
      ENDDO
        CALL fft_forward(work)
      DO IR=1,NGRHO_l
       vpotg(ir)=work(map_grid1d_p(ir))
      ENDDO
!------------------------------------------------------------------------------------------------------------------- 
   do i =1,nnr1
    vpotr(i) = vpotr(i) + CMPLX(0.0_dp,vtmp(i)*grad(i,2),kind=dp)
   enddo 
      CALL  fft_forward(vpotr)
    DO ig=1,ngrho_l
     vfac=twopibya*gvec(1,ig)
     FA  = vpotr(map_grid1d_p(ig)) + vpotr(map_grid1d_m(ig))
     FB  = vpotr(map_grid1d_p(ig)) - vpotr(map_grid1d_m(ig))
     vg1 = 0.5_dp*CMPLX(REAL(fa),AIMAG(fb),kind=dp)
     vg2 = 0.5_dp*CMPLX(AIMAG(fa),-REAL(fb),kind=dp)
     vpotg(ig) = vpotg(ig) + vg1 - vfac*uimag*vg2
    ENDDO

    ! FFT of V2*RHOy and V2*RHOz to G-Space
    vpotr=(0._dp,0._dp)
    DO ir=1,nnr1
        vpotr(ir) = CMPLX(vtmp(ir)*grad(ir,3),vtmp(ir)*grad(ir,4),kind=dp)
    ENDDO

    CALL  fft_forward(vpotr)
    DO ig=1,ngrho_l
      vfac2=twopibya*gvec(2,ig)
      vfac3=twopibya*gvec(3,ig)
      FA  = vpotr(map_grid1d_p(ig)) + vpotr(map_grid1d_m(ig))
      FB  = vpotr(map_grid1d_p(ig)) - vpotr(map_grid1d_m(ig))       
      vg1 = 0.5_dp*CMPLX(REAL(fa),AIMAG(fb),kind=dp)
      vg2 = 0.5_dp*CMPLX(AIMAG(fa),-REAL(fb),kind=dp)
      vpotg(ig) = vpotg(ig) - vfac2*uimag*vg1 - vfac3*uimag*vg2
      !write(802,*)vpotg(ig)
    ENDDO
           
    vpotr=(0._dp,0._dp)
    DO ig=1,ngrho_l
         vpotr(map_grid1d_p(ig)) = vpotg(ig)
         vpotr(map_grid1d_m(ig)) = CONJG(vpotg(ig))
    ENDDO
 
    CALL fft_backward(vpotr)
    
    v_pot(:,1)=vpotr 
  
    DO I=1,NNR1
      ! write(802,*)v_pot(i,1),I
    ENDDO
ELSE
!***************************************************************************
!ALPHA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     WORK=(0._dp,0._dp)
      DO IG=1,NNR1
        WORK(IG)  = dcmplx( v_pot(IG,1),0.D0)
      ENDDO
        CALL fft_forward(work)
      DO IR=1,NGRHO_l
       vpotg(ir)=work(map_grid1d_p(ir))
      ENDDO
!------------------------------------------------------------------------------------------------------------------- 
   do i =1,nnr1
 !   write(551,*)vpotr(i),i
   vpotr(i) = vpotr(i) + CMPLX(0.0_dp,vtmp(i)*grad(i,2),kind=dp)+CMPLX(0.0_dp,vtmpab(i)*grad(i,6),kind=dp)
   !write(551,*)vpotr(i),i
   enddo
      CALL  fft_forward(vpotr)
    DO ig=1,ngrho_l
    
     vfac=twopibya*gvec(1,ig)
     FA  = vpotr(map_grid1d_p(ig)) + vpotr(map_grid1d_m(ig))
     FB  = vpotr(map_grid1d_p(ig)) - vpotr(map_grid1d_m(ig))
     vg1 = 0.5_dp*CMPLX(REAL(fa),AIMAG(fb),kind=dp)
     vg2 = 0.5_dp*CMPLX(AIMAG(fa),-REAL(fb),kind=dp)
     vpotg(ig) = vpotg(ig) + vg1 - vfac*uimag*vg2
   
    ENDDO

    ! FFT of V2*RHOy and V2*RHOz to G-Space
    vpotr=(0._dp,0._dp)
    DO ir=1,nnr1
    vpotr(ir) = CMPLX(vtmp(ir)*grad(ir,3),vtmp(ir)*grad(ir,4),kind=dp)&
            +CMPLX(vtmpab(ir)*grad(ir,7),vtmpab(ir)*grad(ir,8),kind=dp)
  
    ENDDO

    CALL  fft_forward(vpotr)
    DO ig=1,ngrho_l
      vfac2=twopibya*gvec(2,ig)
      vfac3=twopibya*gvec(3,ig)
      FA  = vpotr(map_grid1d_p(ig)) + vpotr(map_grid1d_m(ig))
      FB  = vpotr(map_grid1d_p(ig)) - vpotr(map_grid1d_m(ig))
      vg1 = 0.5_dp*CMPLX(REAL(fa),AIMAG(fb),kind=dp)
      vg2 = 0.5_dp*CMPLX(AIMAG(fa),-REAL(fb),kind=dp)
      vpotg(ig) = vpotg(ig) - vfac2*uimag*vg1 - vfac3*uimag*vg2
      !write(992,*)vpotg(ig)
    ENDDO

    vpotr=(0._dp,0._dp)
    DO ig=1,ngrho_l
         vpotr(map_grid1d_p(ig)) = vpotg(ig)
         vpotr(map_grid1d_m(ig)) = CONJG(vpotg(ig))
    ENDDO

    CALL fft_backward(vpotr)

    v_pot(:,1)=vpotr

    DO I=1,NNR1
       !write(951,*)v_pot(i,1),I
    ENDDO

!***************************************************************************
     WORK=(0._dp,0._dp)
      DO IG=1,NNR1
        WORK(IG)  = dcmplx( v_pot(IG,2),0.D0)
      ENDDO
        CALL fft_forward(work)
     vpotg=(0.0_dp,0.0_dp)
      DO IR=1,NGRHO_l
       vpotg(ir)=work(map_grid1d_p(ir))
      ENDDO
!-----------------------------------------------------------
    DO ir=1,nnr1
       vpotrb(ir) = vpotrb(ir) + CMPLX(0.0_dp,vtmpb(ir)*grad(ir,6),kind=dp)+ CMPLX(0.0_dp,vtmpab(ir)*grad(ir,2),kind=dp)
    ENDDO
    CALL  fft_forward(vpotrb)
    DO ig=1,ngrho_l
     vfac=twopibya*gvec(1,ig)
     FA  = vpotrb(map_grid1d_p(ig)) + vpotrb(map_grid1d_m(ig))
     FB  = vpotrb(map_grid1d_p(ig)) - vpotrb(map_grid1d_m(ig))
     vg1 = 0.5_dp*CMPLX(REAL(fa),AIMAG(fb),kind=dp)
     vg2 = 0.5_dp*CMPLX(AIMAG(fa),-REAL(fb),kind=dp)
     vpotg(ig) = vpotg(ig) + vg1 - vfac*uimag*vg2
    ENDDO

    ! FFT of V2*RHOy and V2*RHOz to G-Space
    vpotrb=(0._dp,0._dp)
    DO ir=1,nnr1
   vpotrb(ir) = CMPLX(vtmpb(ir)*grad(ir,7),vtmpb(ir)*grad(ir,8),kind=dp)+&
   CMPLX(vtmpab(ir)*grad(ir,3),vtmpab(ir)*grad(ir,4),kind=dp)
   ENdDO 
   CALL  fft_forward(vpotrb)
    DO ig=1,ngrho_l
      vfac2=twopibya*gvec(2,ig)
      vfac3=twopibya*gvec(3,ig)
      FA  = vpotrb(map_grid1d_p(ig)) + vpotrb(map_grid1d_m(ig))
      FB  = vpotrb(map_grid1d_p(ig)) - vpotrb(map_grid1d_m(ig))
      vg1 = 0.5_dp*CMPLX(REAL(fa),AIMAG(fb),kind=dp)
      vg2 = 0.5_dp*CMPLX(AIMAG(fa),-REAL(fb),kind=dp)
      vpotg(ig) = vpotg(ig) - vfac2*uimag*vg1 - vfac3*uimag*vg2
      !write(992,*)vpotg(ig)
    ENDDO

    vpotrb=(0._dp,0._dp)
    DO ig=1,ngrho_l
         vpotrb(map_grid1d_p(ig)) = vpotg(ig)
         vpotrb(map_grid1d_m(ig)) = CONJG(vpotg(ig))
    ENDDO

    CALL fft_backward(vpotrb)
    v_pot(:,2)=vpotrb
      DO I=1,NNR1
       !write(952,*)v_pot(i,2),I
    ENDDO


ENDIF
    !v_pot(:,1)=vpotrk:,1) 
        !DEALLOCATE( ex, vrho, vsigma, e_tot, vrho_tot, vsigma_tot )
    CALL MPI_GlobSumR2s(etxc)   
    
!   CALL xc_f90_func_end(xc_func)

deallocate(work)
deallocate(scr_v,scr_vtmp,ex,ro,vsigma_tot,vpot_tot,vpotrb,vpotr,vpotg)
  END SUBROUTINE cal_vpot_ex

SUBROUTINE GRAD_CAL(ro_r,grad)
 USE system_data_types
 USE constants
 USE gvectors
 USE kinds
 USE fft_interface
 IMPLICIT NONE
 INTEGER I,IG,IR
 COMPLEX(KIND=dp), DIMENSION(:), POINTER :: scr_v,scr_vtmp,work
 REAL(KIND=DP)ro_r(nnr1),grad(nnr1,4)
     
ALLOCATE(scr_v(MAX_FFT),scr_vtmp(ngrho))
 !   write(*,*)ro_r(1),rho(1,1),RHO(1,2)
    



!-----------------------------------------------------------------------
  ALLOCATE(WORK(MAX_FFT))
    
     WORK=(0._dp,0._dp)
     DO IG=1,NNR1
       WORK(IG)  = DCMPLX(ro_r(ig),0.D0)
     ENDDO
     CALL fft_forward(WORK)
     DO IR=1,NGRHO_l
       SCR_VTMP(IR)= WORK(map_grid1d_p(IR))
     ENDDO
!-----------------------------------------------------------------------
     !DO I=1,Ngrho_l
     !scr_vtmp(I)=rho_g(I,1)
     !ENDDO
!     ==  FFT OF RHO AND NABLA(X)*RHOE                                ==
     scr_v=(0._dp,0._dp)

     DO IG=1,NGrho_l
        scr_v(map_grid1d_p(IG))=scr_VTMP(IG)-twopibya*gvec(1,IG)*scr_VTMP(IG)
        scr_v(map_grid1d_m(IG))=DCONJG(scr_VTMP(IG)+twopibya*Gvec(1,IG)*scr_VTMP(IG))
     ENDDO
     CALL fft_backward(scr_V)

      DO IR=1,NNR1
       ! RHO(IR)=DREAL(SCR_V(IR))
        GRAD(IR,1)=DIMAG(SCR_V(IR))*DIMAG(SCR_V(IR))
        GRAD(IR,2)=DIMAG(SCR_V(IR))
     ENDDO
!     ==  FFT OF NABLA(Y)*RHO AND NABLA(Z)*RHOE                       ==
      scr_v=(0._dp,0._dp)
      DO IG=1,NGRHO_L
      scr_v(map_grid1d_p(IG))=twopibya*(UIMAG*Gvec(2,IG)-Gvec(3,IG))*scr_VTMP(IG)
      scr_v(map_grid1d_m(IG))=twopibya*(-UIMAG*Gvec(2,IG)+gvec(3,IG))*DCONJG(scr_VTMP(IG))
      ENDDO
      CALL fft_backward(scr_V)

     DO IR=1,NNR1
        GRAD(IR,1)=GRAD(IR,1)+DREAL(scr_V(IR)*DCONJG(scr_V(IR)))
        GRAD(IR,3)=DREAL(scr_V(IR))
        GRAD(IR,4)=DIMAG(scr_V(IR))
     ENDDO
deallocate(work,scr_vtmp,scr_v)
    !write(77,*)"#",grad(1,1),grad(1,2),grad(1,3),grad(1,4),ro_r(1)
END SUBROUTINE



END MODULE
