MODULE fft_interface

 CONTAINS

! Routine to do FFT forward
   SUBROUTINE fft_forward(fz)
      USE kinds
      IMPLICIT NONE

      COMPLEX*16, POINTER :: fz(:)
      INTEGER, PARAMETER :: sign_fft=1 !+1 for forward FFT

!      print *, 'calling fft_forward'
      CALL do_fft(sign_fft,fz)
!      print *, 'done fft_forward'

   END SUBROUTINE fft_forward

! Routine to do FFT backward
   SUBROUTINE fft_backward(fz)
      USE kinds
      IMPLICIT NONE

      COMPLEX*16, POINTER :: fz(:)
      INTEGER, PARAMETER :: sign_fft=-1 !-1 for backward FFT

      CALL do_fft(sign_fft,fz)

   END SUBROUTINE fft_backward
!


   SUBROUTINE do_fft(sign_fft,fz)
      USE kinds
      USE system_data_types, ONLY: nrlead,nrgrids,nrlead_l,nrgrids_l,&
                                   max_nhrays,max_nrgrid,&
                                   Ggridinfo,nrxplane,&
                                   icpu,ncpu,&
                                   work1_fft,work2_fft, max_fft,xscatter_fft_cat
      USE math, ONLY : icopy
                                   

      IMPLICIT NONE

      INTEGER :: sign_fft
      COMPLEX*16, POINTER :: fz(:)
     
      REAL(kind=dp)  :: f_scale
      INTEGER :: qr1,qr1s,qr2s,qr3s,lr1,lr1s,lr2s,lr3s,lfrm,lr1m,lmsq,nhrays,&
                 lrxpl(ncpu,2),sp5(ncpu),sp8(ncpu),sp9(ncpu),&
                 ip,ipx,n1,n2,m,lda,mm
      INTEGER :: xscatter_fft_f(max_nhrays*ncpu)
      INTEGER :: xscatter_local(max_nhrays,2,ncpu) 
     
      f_scale=HUGE(0)
      qr1=nrlead_l(1) !KR1
      qr1s=nrlead(1) !KR1S
      qr2s=nrlead(2) !KR2S
      qr3s=nrlead(3) !KR3S
      lr1=nrgrids_l(1)
      lr1s=nrgrids(1)
      lr2s=nrgrids(2)
      lr3s=nrgrids(3)
      lfrm=max_nhrays
      lr1m=max_nrgrid
      lmsq=max_nhrays
      !all_to_all_single=.false. !TR4A2A
      nhrays=Ggridinfo(7,icpu) !MFRAYS ! NHRAYS?
      !print *, 'max_fft=', max_fft, ' me=', icpu
!      IF(.not.ASSOCIATED(work1_fft))THEN
        ALLOCATE(work1_fft(max_fft))
!      END IF
!      IF(.not.ASSOCIATED(work2_fft))THEN
        ALLOCATE(work2_fft(max_fft))
!      END IF

      xscatter_local = RESHAPE(xscatter_fft_cat,(/max_nhrays,2,ncpu/)) !MSP
!      print *, 'do_fft 3'

      DO ip=1,ncpu
        lrxpl(ip,1)=nrxplane(ip,1)
        lrxpl(ip,2)=nrxplane(ip,2)
        sp5(ip)=Ggridinfo(3,ip)
        sp8(ip)=Ggridinfo(7,ip)
        sp9(ip)=Ggridinfo(6,ip)
        ipx=lmsq*(ip-1)
        IF(max_nhrays>0)THEN
!TODO: array indices are having some problem here....
          CALL icopy(max_nhrays,xscatter_local(1,1,ip),1,xscatter_fft_f(ipx+1),1) ! MSQF
!TODO: array indices are having some problem here....
          !CALL icopy(max_nhrays,xscatter_fft_cat(1,2,ip),1,xscatter_fft_s(ipx),1) !MSQS
        ENDIF
      ENDDO         
      !print *, 'do_fft 4'


      IF(sign_fft.eq.1)THEN
         f_scale=1.d0

         n1=lrxpl(icpu,1) 
         n2=lrxpl(icpu,2) 

      !print *, 'do_fft 5'
         CALL phasen(fz,qr1,qr2s,qr3s,n1,n2,lr2s,lr3s)
      !print *, 'do_fft 6',fz(1)
         m=qr1*qr2s
      !print *, m,qr3s,lr3s,m
         CALL mltfft('T','N',fz,m,qr3s,work1_fft,qr3s,m,lr3s,m,sign_fft,f_scale)
      !print *, 'do_fft 7'
         m=qr1*qr3s
         CALL mltfft('T','N',work1_fft,m,qr2s,work2_fft,qr2s,m,lr2s,m,sign_fft,f_scale)
      !print *, 'do_fft 8'
         lda=lfrm*lr1m
         mm=qr2s*qr3s
         CALL pack_y2x(work1_fft,work2_fft,mm,lr1,lda,xscatter_fft_f,lmsq,sp8,max_fft,ncpu)
      !print *, 'do_fft 9'
         CALL fft_comm(work1_fft,work2_fft,lda) !all to all communication
      !print *, 'do_fft 10'
         CALL unpack_y2x(work1_fft,work2_fft,mm,nhrays,lda,lrxpl,sp5,max_fft,ncpu) !TILL HERE
      !print *, 'do_fft 11'
         f_scale=1.D0/DBLE(lr1s*lr2s*lr3s)
         m=nhrays
         CALL mltfft('T','N',work1_fft,m,qr1s,fz,qr1s,m,lr1s,m,sign_fft,f_scale)
      !print *, 'do_fft 12'
         !fz(:)=fz(:)*f_scale
      ELSE IF(sign_fft.eq.-1)THEN
!         f_scale=1.d0 sks dbg 
         f_scale=1.d0
         m=nhrays
         CALL MLTFFT('N','T',fz,qr1s,m,work1_fft,m,qr1s,lr1s,m,sign_fft,f_scale)
         lda=lfrm*lr1m
         mm=qr2s*qr3s
         CALL pack_x2y(work1_fft,work2_fft,nhrays,lda,lrxpl,sp5,max_fft,ncpu)
         call fft_comm(work2_fft,work1_fft,lda)
         work2_fft(:)=(0._dp,0._dp)
         call unpack_x2y(work1_fft,work2_fft,mm,lr1,lda,xscatter_fft_f,lmsq,sp8,max_fft,ncpu)
         m=qr1*qr3s
         CALL MLTFFT('N','T',work2_fft,QR2S,M,work1_fft,M,QR2S,LR2S,M,sign_fft,f_scale)
         m=qr1*qr2s
         CALL MLTFFT('N','T',work1_fft,QR3S,M,fz,M,QR3S,LR3S,M,sign_fft,f_scale)
         n1=LRXPL(icpu,1)
         n2=LRXPL(icpu,2)
         CALL PHASEN(fz,QR1,QR2S,QR3S,n1,n2,LR2S,LR3S)
      END IF
  deallocate(work1_fft)
  deallocate(work2_fft)

  END SUBROUTINE do_fft

!based on fftprp.F
  SUBROUTINE prepare_fft()

    USE kinds
    USE mympi, ONLY: mpi_globsumi2,mpi_concat,mpi_globsumi
    USE system_data_types, ONLY : map_grid1d_p,map_grid1d_m,inyh,&
                                  nrlead,nrlead_l,nrgrids,ngrho_l,ngpw_l,&
                                  ncpu,icpu,&
                                  xrays_pw,& !MG
                                  ngvy,& !MY
                                  ngvz,& !MZ
                                  gstart_ygrid,gstart_zgrid,&
                                  gend_ygrid,gend_zgrid,&
                                  max_nrgrid,max_ngrays,max_nhrays,max_nrlead,&
                                  max_fft,&
                                  xscatter_fft_cat,map_grid1d_p_pw,map_grid1d_m_pw,&
                                  Ggridinfo
                                   
    USE math, ONLY: icopy
    IMPLICIT NONE
    INTEGER :: nr1,nr2,nr3,nh1,nh2,nh3, &
               ig,kr1,kr2s,kr3s,ip,&
               i,j,ii,jj,i1,i2,i3,j1,j2,j3,&
               ir,ngrays,nhrays
    INTEGER, DIMENSION(:), POINTER :: xscatter_fft
               

    kr2s=nrlead_l(2) ; kr3s=nrlead_l(3)

    nr1=nrgrids(1) ; nr2=nrgrids(2) ; nr3=nrgrids(3)

    nh1=nr1/2+1 ; nh2=nr2/2+1 ; nh3=nr3/2+1

    !print *, 'kr2s,kr3s=', kr2s,kr3s
    !print *, 'nr2,nr3=', nr2,nr3
!
!allocation of arrays, MG, MZ, MY (gather array for fft along X)
!
    ALLOCATE(xrays_pw(kr2s,kr3s)) !kr2s, kr3s TODO: is kr2s equal to nr2s?
    ALLOCATE(ngvy(2*nr2)) !2 x nr3s
    ALLOCATE(ngvz(2*nr3)) !2 x nr3s
    xrays_pw(:,:)=0 ; ngvz(:)=0 ; ngvy(:)=0
    !print *, 'position 1'

    DO ig=1,ngpw_l !for G vectors less than the PW-cutoff
      i2=inyh(2,ig) !+G
      i3=inyh(3,ig)
      j2=-i2+2*nh2  !-G
      j3=-i3+2*nh3

!MG -> xrays_pw
!MY -> ngvy
!MZ -> ngvz
      !TODO: check if this is already there in gvectors.f90?
      xrays_pw(i2,i3)=xrays_pw(i2,i3)+1  !Count x-rays for any y or z grid
      xrays_pw(j2,j3)=xrays_pw(j2,j3)+1  !Count x-rays for any y or z grid
      
      ngvy(i2)=ngvy(i2)+1 !count G-vectors for a y-grid point
      ngvy(j2)=ngvy(j2)+1

      ngvz(i3)=ngvz(i3)+1 !count G-vectors for a z-grid point
      ngvz(j3)=ngvz(j3)+1

    END DO
    !print *, 'position 2'
    
    CALL MPI_GlobSumI2(ngvy,nr2) !Sum the number across all the processors
    CALL MPI_GlobSumI2(ngvz,nr3) 
    !print *, 'position 3'

!   Find the position of the grid which has the non-zero G component along Y
    gstart_ygrid=1 !KR2MIN
    loop1: DO i=1,kr2s
      IF(ngvy(i)/=0)THEN
        gstart_ygrid=i ; EXIT loop1
      END IF
    END DO loop1
   
    gend_ygrid=kr2s !KR2MAX
    loop2: DO i=kr2s,1,-1
      IF(ngvy(i)/=0)THEN
        gend_ygrid=i ; EXIT loop2
      END IF
    END DO loop2

!   Find the position of the grid which has the non-zero G component along Z
    gstart_zgrid=1  !KR3MIN
    loop3: DO i=1,kr3s
      IF(ngvz(i)/=0)THEN
        gstart_zgrid=i ; EXIT loop3
      END IF
    END DO loop3
   
    gend_zgrid=kr3s   !KR3MAX
    loop4: DO i=kr3s,1,-1
      IF(ngvz(i)/=0)THEN
        gend_zgrid=i ; EXIT loop4
      END IF
    END DO loop4
    !print *, 'position 4'
!   Find maximum number of grids, rays for pw/density cutoffs
    max_nrgrid=0; max_ngrays=0 ; max_nhrays=0 
    DO ip=1,ncpu
      max_nrgrid = MAX(max_nrgrid,Ggridinfo(3,ip))  !NR1M !SPARM 5
      max_ngrays = MAX(max_ngrays,Ggridinfo(6,ip))  !NGRM !SPARM 8
      max_nhrays = MAX(max_nhrays,Ggridinfo(7,ip))  !NHRM !SPARM 9
    END DO
    kr1=nrlead_l(1)
    max_nrlead=MAX(max_nrgrid+MOD(max_nrgrid+1,2),kr1) !KR1M
    !print *, 'kr1m|',max_nrgrid,kr1
    !print *, 'position 5'

    ii=0
    DO j=1,nrlead(3)
      DO i=1,nrlead(2)
        IF(xrays_pw(i,j)/=0)THEN
          ii=ii+1 
          xrays_pw(i,j)=ii
        END IF
      END DO
    END DO
    !print *, ' ii =', ii, 'ngrays =', Ggridinfo(6,icpu), ' icpu=', icpu
    if(ii.ne.Ggridinfo(6,icpu))STOP 'programming error|ngrays'
    DO ig=ngpw_l+1,ngrho_l
      i2=inyh(2,ig)
      i3=inyh(3,ig)
      j2=-i2+2*nh2
      j3=-i3+2*nh3
      jj=xrays_pw(i2,i3)
      IF(jj==0)xrays_pw(i2,i3)=-1
      jj=xrays_pw(j2,j3)
      IF(jj==0)xrays_pw(j2,j3)=-1
    END DO
    DO j=1,nrlead(3)
      DO i=1,nrlead(2)
        IF(xrays_pw(i,j)<0)THEN
         ii=ii+1
         xrays_pw(i,j)=ii
        END IF
      END DO
    END DO
    !nhray=ii !is it needed?
    !print *, ' ii =', ii, 'nhrays =', Ggridinfo(7,icpu), ' icpu=', icpu
    if(ii.ne.Ggridinfo(7,icpu))STOP 'programming error | nhrays'

    !Scatter arrays for FFT along X
    ALLOCATE(xscatter_fft(2*max_nhrays)) !MS
    xscatter_fft(:)=0 

    DO i=1,nrlead(2)
      DO j=1,nrlead(3)
        ii=xrays_pw(i,j)
        IF(ii>0)THEN
!if(ii>2*max_nhrays)STOP 'ERROR in programming xscatter_fft1'
!if(ii+max_nhrays>2*max_nhrays)THEN
!    print *, 'ii =', ii
!    print *, 'max_nhrays=',max_nhrays
!     STOP 'ERROR in programming xscatter_fft2'
!end if
          xscatter_fft(ii)=i
          xscatter_fft(max_nhrays+ii)=j
        END IF
       END DO
     END DO
          
    !concatenate gather-scatter array
    ALLOCATE(xscatter_fft_cat(2*max_nhrays*ncpu+1)) !MSP
    CALL MPI_concat(xscatter_fft,xscatter_fft_cat,2*max_nhrays)

    !translate i,j to a single g/s index
    DO ip=0,ncpu-1
      nhrays=Ggridinfo(7,ip+1) !nhrays 
      DO ir=1,nhrays
        jj=ir+ip*2*max_nhrays
         i=xscatter_fft_cat(jj)
         j=xscatter_fft_cat(jj+max_nhrays)
         xscatter_fft_cat(jj)=i+(j-1)*nrlead(2)
         IF(ir<=Ggridinfo(6,ip+1))xscatter_fft_cat(jj+max_nhrays)=i+(j-gstart_zgrid)*nrlead(2) 
      END DO
    END DO

    !redefine nzh and indz for compressed storage
    !nn2=1
    DO ig=1,ngrho_l
      i1=inyh(1,ig)
      i2=inyh(2,ig)
      i3=inyh(3,ig)
      j1=-i1+2*nh1
      j2=-i2+2*nh2
      j3=-i3+2*nh3
      map_grid1d_p(ig)=i1+(xrays_pw(i2,i3)-1)*nrlead(1)  !NZH
      map_grid1d_m(ig)=j1+(xrays_pw(j2,j3)-1)*nrlead(1)  !INDZS
      !write(16,*)"PSIcalculation_1",ig,map_grid1d_p(ig),map_grid1d_m(ig)
    END DO

    IF(.not.ASSOCIATED(map_grid1d_p_pw))ALLOCATE(map_grid1d_p_pw(ngpw_l)) !TODO check 
    IF(.not.ASSOCIATED(map_grid1d_m_pw))ALLOCATE(map_grid1d_m_pw(ngpw_l)) !TODO check 

    DO ig=1,ngpw_l
      map_grid1d_p_pw(ig)=map_grid1d_p(ig)
      map_grid1d_m_pw(ig)=map_grid1d_m(ig)
    END DO

    DEALLOCATE(xscatter_fft)

!   Array sizes
    max_fft=MAX(max_nrlead*nrlead(2)*nrlead(3),ncpu*max_nrgrid*max_nhrays)
    !print *, 'max_fft| ', max_nrlead,nrlead(2),nrlead(3),ncpu,max_nrgrid,max_nhrays
    ! print *, 'max_fft =', ncpu,max_nrgrid,max_nhrays
      !LMSQMAX = NHRM
      !LNZF    = NHG
      !LNZS    = NGW        
    DEALLOCATE(xrays_pw,ngvy,ngvz) 
  END SUBROUTINE prepare_fft

!based on mltfft
  SUBROUTINE mltfft(transa,transb,a,ldax,lday,b,ldbx,ldby,n,m,i_sign,scale_factor)
         !CALL mltfft('T','N',fz,m,qr3s,work1_fft,qr3s,m,lr3s,m,sign_ft,f_scale)
    USE KINDS
    IMPLICIT NONE
    CHARACTER(LEN=1) :: transa,transb
    INTEGER :: ldax,lday,ldbx,ldby,n,m,i_sign
    COMPLEX*16 :: a(ldax,*), b(ldbx,*)
    REAL(KIND=DP) :: scale_factor
    LOGICAL :: use_wisdom

    INTEGER*8 :: plan !TODO need to check this 
    
    INTEGER, PARAMETER :: fftw_b=1,fftw_f=-1
    INTEGER, PARAMETER :: fftw_me=0, fftw_es=64

    INTEGER :: fftw_dir, fftw_flags
    INTEGER :: i,j
    logical :: tscal

    IF(i_sign==1)THEN
      fftw_dir=fftw_f
    ELSE
      fftw_dir=fftw_b
    END IF
    use_wisdom=.false. !TOdO ACTIVIATE IT

    IF(use_wisdom)THEN
      fftw_flags=fftw_me
    ELSE
      fftw_flags=fftw_es
    END IF
      TSCAL=(ABS(scale_factor-1.D0).GT.1.D-12)
    !print *, 'mltfft-1'
    !print *, 'mlfft |, ldax,lday,ldbx,ldby,n,m,isign'
    !print *, 'mlfft |',ldax,lday,ldbx,ldby,n,m,i_sign
    !print *, 'mlfft | testing a ',  a(1,1),a(ldax,1)
    !print *, 'mlfft | testing b ',  b(1,1),b(ldbx,1)
 
    !print *, 'mlfft | testing2 a ',  a(ldax,lday)
    !print *, 'mlfft | testing2 b ',  b(ldbx,ldby)

    IF(transa.EQ.'N'.OR.transa.EQ.'n') THEN
        IF(transb.EQ.'N'.OR.transb.EQ.'n') THEN
         !WRITE(8,*)"1",m,n,ldax,ldbx
         CALL dfftw_plan_many_dft(plan,1,n,m,a,n,1,ldax,b,n,1,ldbx,fftw_dir,fftw_flags)
         CALL dfftw_execute_dft(plan,a,b)
          IF(TSCAL) THEN
            DO I = 1,M
              DO J = 1,N
                B(J,I)=scale_factor*B(J,I)
              ENDDO
            ENDDO
          ENDIF


        ELSE
!WRITE(8,*)"2",m,n,ldax,ldbx
          CALL dfftw_plan_many_dft(plan,1,n,m,a,n,1,ldax,b,n,ldbx,1,fftw_dir,fftw_flags)
          CALL dfftw_execute_dft(plan,a,b)
           !b(:,:)=scale_factor*b(:,:)

        IF(TSCAL) THEN
            DO I = 1,N
              DO J = 1,M
                B(J,I)=scale_factor*B(J,I)
              ENDDO
            ENDDO
          ENDIF


        END IF
    ELSE
        IF(transb.EQ.'N'.OR.transb.EQ.'n') THEN
    !print *, 'mltfft-2'
    !WRITE(8,*)"3",m,n,ldax,ldbx,b(1,1)

     !     print *, 'n=',n, 'm=',m, 'ldax=', ldax, 'ldbx=',ldbx
          CALL dfftw_plan_many_dft(plan,1,n,m,a,n,ldax,1,b,n,1,ldbx,fftw_dir,fftw_flags)
    !print *, 'mltfft-3'
          CALL dfftw_execute_dft(plan,a,b)
           !b=scale_factor*b
    !print *, 'mltfft-5',b(1,1)
          IF(TSCAL) THEN
            DO I = 1,M
              DO J = 1,N
                B(J,I)=scale_factor*B(J,I)
              ENDDO
            ENDDO
          ENDIF
        ELSE
!WRITE(8,*)"4",m,n,ldax,ldbx

          CALL dfftw_plan_many_dft(plan,1,n,m,a,n,ldax,1,b,n,ldbx,1,fftw_dir,fftw_flags)
          CALL dfftw_execute_dft(plan,a,b)
           !b=scale_factor*b
          IF(TSCAL) THEN
            DO I = 1,N
              DO J = 1,M
                B(J,I)=scale_factor*B(J,I)
              ENDDO
            ENDDO
          ENDIF
        END IF
    END IF
    !print *, 'mltfft-6'
    CALL dfftw_destroy_plan(plan)
    !print *, 'mltfft-7'

    IF(transb.EQ.'N'.OR.transb.EQ.'n') THEN
        DO j=m+1,ldby
          DO i=1,ldbx
            b(i,j)=DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
        DO j=1,m
          DO i=n+1,ldbx
            b(I,J)=DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
    ELSE
        DO j=n+1,ldby
          DO i=1,ldbx
            b(I,J)=DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
        DO j=1,n
          DO i=m+1,ldbx
            b(i,j)=DCMPLX(0.0D0,0.0D0)
          ENDDO
        ENDDO
    ENDIF
    !print *, 'mltfft-8'
  END SUBROUTINE mltfft

  SUBROUTINE FFT_COMM(XF,YF,LDA)
    USE kinds
    USE mympi, ONLY: MPI_trans
    IMPLICIT NONE

    COMPLEX*16 :: XF(*),YF(*)
    INTEGER    :: LDA

    INTEGER   BLKLEN,ISUB1
!..All to all communication
    CALL MPI_trans(XF,YF,LDA)
  END SUBROUTINE FFT_COMM



  SUBROUTINE phasen(f,qr1,qr2s,qr3s,n1,n2,lr2s,lr3s)
  USE KINDS
  IMPLICIT NONE
  INTEGER :: qr1,qr2s,qr3s,n1,n2,lr2s,lr3s
  COMPLEX*16 :: f(qr1,qr2s,qr3s)
  INTEGER :: k,j,i,ii,ijk,isub
  REAL(KIND=DP), PARAMETER, DIMENSION(2)  :: pf=[1.d0,-1.d0]

  DO k=1,qr3s
    DO j=1,qr2s
      DO i=n1,n2
        ii=i-n1+1
        ijk=MOD(k+j+i+1,2)+1
        f(ii,j,k)=f(ii,j,k)*pf(ijk)
      ENDDO
    ENDDO
  ENDDO
  END SUBROUTINE phasen

  SUBROUTINE unpack_y2x(xf,yf,m,nrays,lda,jrxpl,sp5,maxfft,mproc)
  IMPLICIT NONE
  INTEGER :: m,nrays,maxfft,mproc,lda
  COMPLEX*16 :: xf(*), yf(*)
  integer :: jrxpl(1:mproc), sp5(1:mproc)
  
  integer :: blklen,isub1,ip,mrxp,i,j,ii,jj,isub2,nrx,nrs,ipp,k
 
   DO ip=1,mproc
     nrx = sp5(ip)*nrays
     nrs = (jrxpl(ip)-1)*nrays + 1
     ipp = (ip-1)*lda + 1
      CALL DCOPY(2*nrx,yf(ipp),1,xf(nrs),1)
   ENDDO
  END SUBROUTINE unpack_y2x

  SUBROUTINE pack_y2x(xf,yf,m,lr1,lda,msp,lmsp,sp8,maxfft,mproc)
  IMPLICIT NONE
  INTEGER :: m,maxfft,mproc,lda,lmsp,lr1
  COMPLEX*16 :: xf(*), yf(*)
  integer :: msp(lmsp,*),sp8(*)
  integer :: blklen,isub1,ip,mxrp,i,j,ii,jj,isub2,nrx,nrs,ipp,k
  
  DO ip=1,mproc
     mxrp = sp8(ip)
     DO i=1,lr1
       ii = (ip-1)*lda + (i-1)*mxrp + 1
       jj = (i-1)*m + 1
       CALL zgthr(mxrp,yf(jj),xf(ii),msp(1,ip))
     ENDDO
  ENDDO
  END SUBROUTINE pack_y2x

  subroutine pack_x2y(xf,yf,nrays,lda,jrxpl,sp5,maxfft,mproc)
      IMPLICIT NONE

      complex*16 xf(*),yf(*)
      integer nrays,maxfft,mproc,lda
      integer jrxpl(*),sp5(*)

      integer   blklen,isub1,ip,mxrp,i,ii,jj,isub2,nrx,nrs,ipp,k

      DO IP=1,mproc
        nrx = sp5(ip)*nrays
        nrs = (jrxpl(ip)-1)*nrays + 1
        ipp = (ip-1)*lda + 1
        call dcopy(2*nrx,xf(nrs),1,yf(ipp),1)
      ENDDO

  END SUBROUTINE pack_x2y


  SUBROUTINE unpack_x2y(xf,yf,m,lr1,lda,msp,lmsp,sp8,maxfft,mproc)
      IMPLICIT NONE
      INTEGER m,maxfft,mproc,lda,lmsp,lr1
      COMPLEX*16 xf(*),yf(*)
      INTEGER msp(lmsp,*),sp8(*)

      INTEGER ::  BLKLEN,ISUB1,IP,MXRP,I,II,JJ,ISUB2,NRX,NRS,IPP,K
       DO IP=1,MPROC
         MXRP = SP8(IP)
         DO I=1,LR1
            ii = (ip-1)*lda + (i-1)*mxrp + 1
            jj = (i-1)*m + 1
            CALL zsctr(mxrp,xf(ii),msp(1,ip),yf(jj))
         ENDDO
      ENDDO
  END SUBROUTINE unpack_x2y


  SUBROUTINE zgthr(n,a,b,ind)
     IMPLICIT NONE
     integer    n,ind(n)
     complex*16 a(*),b(*)

     INTEGER I

      do i=1,n
        b(i)=a(ind(i))
      enddo
   END SUBROUTINE zgthr

    SUBROUTINE zsctr(n,a,ind,b)
      IMPLICIT NONE
      INTEGER   N,IND(N)
      COMPLEX*16 A(N),B(*)
      INTEGER I
      DO I=1,N
        B(IND(I))=A(I)
      ENDDO
    END SUBROUTINE zsctr


END MODULE fft_interface
