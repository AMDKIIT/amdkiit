MODULE pseudopotential
  CONTAINS
  SUBROUTINE pp_energy
  USE max_parameter_pp
  USE system_data_types
  USE kinds
  USE constants
  USE fft_interface
  USE mympi
  USE READupf 
  IMPLICIT NONE

  INTEGER I,ig,IG_s,IR,is,iv,ia,isa0,isa
  INTEGER li,ki,l2,lj,jv,kj,l
  REAL(KIND=DP):: sum
  COMPLEX*16 el
  !COMPLEX(KIND=dp), DIMENSION(:), POINTER :: work

  !ALLOCATE(WORK(MAX_FFT))
  !   WORK=(0._dp,0._dp)

     el=(0._dp,0._dp) 

     IF(g0_stat)then
     el=0.5d0*eigrxvps(1)*DCONJG(rho_g(1,1))
     ig_s=2
     else
     el=(0.0D0,0.0D0) 
     ig_s=1
     endif

     DO ig=ig_s,ngrho_l
      el=el+(eigrxvps(ig)*DCONJG(rho_g(ig,1)))
     ENDDO
     CALL MPI_GlobSumC2s(el)
     eloc=2.d0*DREAL(el)*cell_volume
     !------------------------*********************-----------------------!

     if(l_upf) then
     CALL upf_nlpro
     else
     CALL nlpro
     endif

!    OVERLAP OF PROJECTORS AND BANDS
     if(.not.gopt) then
        if(l_upf) then 
             CALL upf_fnl(C_0,nstate)
        else
             CALL fnl(C_0,nstate)
        endif
     endif

     enl=0.0
     isa0=0
    if(.not.l_upf)then
       DO is=1,sp_t
    DO iv=1,ngh(is)
      sum=0.0d0
        DO i=1,nstate
          DO ia=1,atom_p_sp(is)
            isa=isa0+ia
            sum=sum+OCCUPATION(i)*1*nl(1,isa,iv,i)*nl(1,isa,iv,i)
         
          ENDDO
        ENDDO
      enl=enl+wsg(is,iv)*sum
     
    ENDDO
  isa0=isa0+atom_p_sp(is)
  ENDDO
  else
     Do is =1,sp_t
             DO iv=1,upf_ngh(is)
                l=upf_nghtol(iv,is)+1
                ki=lfval(iv,is)
                li=lpval(iv,is)
                DO jv=iv,upf_ngh(is)
                   l2=upf_nghtol(jv,is)+1
                   lj=lpval(jv,is)
                    ! print*,"*",jv,l2,lj
                   IF (l2.EQ.l.AND.li.EQ.lj) THEN
                      kj=lfval(jv,is)
                      sum=0.D0
                      DO i=1,nstate
                       !weight=wk(ik)*crge%f(i,ik)
                         IF (occupation(i).EQ.0.D0) GOTO 2000
                               DO ia=1,atom_p_sp(is)
                                  isa=isa0+ia
                                  sum=sum+occupation(i)*nl(1,isa,iv,i)*nl(1,isa,jv,i)
                               ENDDO
                         2000                     CONTINUE
                     ENDDO
                      IF (iv.NE.jv) sum=2.D0*sum
                      enl=enl+sum*hlsg(ki,kj,l,is)
                   ENDIF
                ENDDO
             ENDDO
     END DO
endif
   
  return
  END SUBROUTINE pp_energy
  



  SUBROUTINE nlpro
  USE system_data_types

  IMPLICIT NONE

  INTEGER   lpmax,is,iv,ig,l,lp,il
  REAL(KIND=DP)::    vol,xx,tw
  REAL(KIND=DP)::    dasum
  external  dasum
  REAL(KIND=dp), DIMENSION(:,:,:),   POINTER ::ylmb
  REAL(KIND=dp), DIMENSION(:),   POINTER ::m_wns
  REAL(KIND=DP),DIMENSION(:),POINTER::  ggng
  IF(.NOT.ASSOCIATED(twnl))ALLOCATE(twnl(ngpw_l,nghmax,spmax))
  ALLOCATE(m_wns(ngpw_l),ggng(nspln))
  m_wns=0.0d0
  DO il=1,nspln
            ggng(il)=(il-1)*(Gcutoff%pw)/DBLE(nspln-1)
  enddo

  lpmax=0
  DO is=1,sp_t
    IF(ngh(is).gt.0) THEN
      IF(tkb(is)) THEN
        iv=ngh(is)
        lpmax=max(lpmax,nghcom(iv,is))
      ENDIF
    ENDIF
  ENDDO
  allocate(ylmb(ngrho_l,lpmax,1))
  ylmb=0.0d0
  vol=1.d0/sqrt(cell_volume)
  DO lp=1,lpmax
  !       cALL ylmr2(3,ngrho_l,hg,gvec,ylmb(1,lp,1))
    CALL compute_yL(lp,ngrho_l,hg,gvec,ylmb(1,lp,1))
  ENDDO
  !write(601,*)lp
  !enddo
  DO is=1,sp_t
    DO iv=1,ngh(is)
      IF(tkb(is)) THEN
        l=nghtol(iv,is)+1
        lp=nghcom(iv,is)
      ENDIF
      
      xx=dasum(nspln,wns(1,iv,is),1)
      CALL spline_inter(nspln,ggng,wns(1,iv,is),hg(1),m_wns(1),ngpw_l)
      IF(xx.gt.1.d-12) THEN
        DO ig=1,ngpw_l
          twnl(ig,iv,is)=ylmb(ig,lp,1)*vol*m_wns(ig)
        ENDDO
      ENDIF
    ENDDO
  ENDDO
  DEALLOCATE(m_wns,ggng,ylmb)
  END SUBROUTINE nlpro

  SUBROUTINE fnl(C0,nstat)
  USE max_parameter_pp
  USE system_data_types
  USE MYMPI
  USE constants
  USE math
  IMPLICIT NONE
  INTEGER   isa0,is,iv,ia,ig,isa,ikind,nat,ik,i,nstat,ii
  REAL(KIND=DP)    er,ei,check,xg,xrg,flx,t,sum
  COMPLEX*16 ci,scr(ngpw_l,na_max),C0(NGPW_l,NSTAT)
  ALLOCATE(nl(1,atom_t,nghmax,nstat))
  ikind=1
  nl=0.0D0
  isa0=0
  DO is=1,sp_t
    DO iv=1,ngh(is)
      ci=(0.0d0,-1.0d0)**nghtol(iv,is)
      !write(72,*)(DABS(DREAL(ci))),c0(1,1)
      IF(DABS(DREAL(ci)).gt.0.5d0) THEN
        DO ia=1,atom_p_sp(is)
          isa=isa0+ia
          DO ig=1,ngpw_l
            er=DREAL(eigr(ig,isa))
            ei=DIMAG(eigr(ig,isa))
            t=twnl(ig,iv,is)*(DREAL(ci))
            scr(ig,ia)=DCMPLX(t*er,t*ei)
            !if(icpu==1)write(72,*)scr(ig,ia),ig,ia
          ENDDO
        ENDDO
      ELSE
        DO ia=1,atom_p_sp(is)
          isa=isa0+ia
          DO ig=1,ngpw_l
            er=DREAL(eigr(ig,isa))
            ei=DIMAG(eigr(ig,isa))
            t=twnl(ig,iv,is)*DIMAG(ci)
            scr(ig,ia)=DCMPLX(-t*ei,t*er)
            !if(icpu==1)write(660,*)scr(ig,ia),ig,ia
          ENDDO
        ENDDO
      ENDIF 
      IF (G0_STAT)THEN
      DO ia=1,atom_p_sp(is)
        scr(1,ia)=0.5d0*scr(1,ia)
      ENDDO
      ENDIF

      IF(atom_p_sp(is).gt.1) THEN
        CALL dgemm('T','N',atom_p_sp(is),nstat,2*ngpw_l,2.d0,scr(1,1),&
      &          2*ngpw_l,c0(1,1),2*ngpw_l,0.0d0,&
      &          nl(1,isa0+1,iv,1),atom_t*nghmax)
      ELSE
        CALL dgemv('t',2*ngpw_l,nstat,2.d0,c0(1,1),2*ngpw_l,scr(1,1),&
      &          1,0.0d0,nl(1,isa0+1,iv,1),atom_t*nghmax)
      ENDIF
      !write(72,*)nl(1,isa0+1,iv,1),isa0
    ENDDO
    isa0=isa0+atom_p_sp(is)
   ENDDO

  CALL MPI_GlobSumR2(NL,nstat*atom_t*NGHMAX)
  RETURN
  END SUBROUTINE fnl
  !!-----------------------------------------------------------------------------





  !________-----------------------------------------------------------------------
  SUBROUTINE upf_nlppgrid(IS)
  USE max_parameter_pp
  USE system_data_types
  USE interpolgrid
  USE constants
  USE math
  USE bessel_func 
  IMPLICIT NONE

  INTEGER   is,il,l,j,iv,ir,lm,k,lold,ierr
  REAL(KIND=DP)::    gmax(splnmax),wsgtmp(splnmax)
  REAL(KIND=DP)::   wc(splnmax),wfc(splnmax),wtemp(splnmax)
  !REAL(KIND=DP)::    ggng(1:nspln)      
  REAL(KIND=DP)::    jl(meshw(is))!,wns(nspln,2,upf_nghmax,sp_t)
  REAL(KIND=DP)::   xg,tmp,xgr
  REAL(KIND=DP)::    wfint(meshw(is)),fnt(splnmax),loggrid
  REAL(KIND=DP),DIMENSION(:),POINTER::  ggng 
  ALLOCATE(ggng(nspln))
fnt = 0.0D0
wfint = 0.0D0
gmax = 0.0D0
wsgtmp = 0.0D0
wc = 0.0D0
wfc = 0.0D0
wtemp = 0.0D0


      DO il=1,nspln
        ggng(il)=(il-1)*(Gcutoff%pw)/DBLE(nspln-1)
     ENDDO

    lold=-1
    DO iv=1,upf_ngh(is)
      l=upf_nghtol(iv,is)+1
      IF(l.eq.lold) THEN
      CALL dcopy(2*nspln,wns(1,iv-1,is),1,wns(1,iv,is),1)
      ELSE
      wns(:,iv,is)=0.0d0
      !CALL initzero(wns(1,1,iv,is),nspln)
        !IF(l.ne.skip(is)) THEN
          DO il=1,nspln
            xg=sqrt(ggng(il))*twopi/a_lattice
            !CALL bess(xg,l,meshw(is),rr(1,is),jl)
            DO ir=1,meshw(is)
              xgr=xg*rr(ir,is)
              jl(ir)=sph_bes(l-1,xgr)
              fnt(ir)=RW(IR,IS)*rr(IR,IS)*upf_gnl(ir,is,l)*jl(ir) !!CHECK RR AND RW
            ENDDO
            CALL simpsn(meshw(is),fnt,tmp)
            wns(il,iv,is)=fourpi*tmp
          ENDDO
           !CALL spline_inter(nspln,ggng,wns(1,1,iv,is),hg(1),wns(1,2,iv,is),hg(2)-hg(1),nspln)
      ENDIF
      lold=l
    ENDDO
  RETURN
  DEALLOCATE(GGNG)
  END SUBROUTINE upf_nlppgrid
  !---------------------------------------------------------------------------
  SUBROUTINE upf_nlpro
  USE system_data_types

  IMPLICIT NONE
  !rinforce.f90
  INTEGER   lpmax,is,iv,ig,l,lp,il
  REAL(KIND=DP)::    vol,xx,tw!,ylmb(ngrho_l,1,1)!todo
  REAL(KIND=DP)::    dasum
  external  dasum
 REAL(KIND=dp), DIMENSION(:,:,:),   POINTER ::ylmb
   REAL(KIND=dp), DIMENSION(:),   POINTER ::m_wns
  REAL(KIND=DP),DIMENSION(:),POINTER::  ggng
  ALLOCATE(m_wns(ngpw_l),ggng(nspln))
   m_wns=0.0d0
  
      DO il=1,nspln
        ggng(il)=(il-1)*(Gcutoff%pw)/DBLE(nspln-1)
     ENDDO

  IF(.NOT.ASSOCIATED(twnl))ALLOCATE(twnl(ngpw_l,upf_nghmax,spmax))
  !CALL initzero(ylmb,ngrho_l)
  !---------------------------------
  lpmax=0
  DO is=1,sp_t
    IF(upf_ngh(is).gt.0) THEN
      !if(l_upf(is))then
          IV=upf_NGH(IS)
          LPMAX=MAX(LPMAX,LPVAL(IV,IS))
      !ENDIF
    ENDIF
  ENDDO
  
  !---------------------------------
  allocate(ylmb(ngrho_l,lpmax,1))
  ylmb=0.0d0
  vol=1.d0/sqrt(cell_volume)
  DO lp=1,lpmax
  !   CALL ylmr2(lp,ngrho_l,hg,gvec,ylmb(1,lp,1))
      CALL compute_yL(lp,ngrho_l,hg,gvec,ylmb(1,lp,1))
  ENDDO
  DO is=1,sp_t
    DO iv=1,upf_ngh(is)
      LP=LPVAL(IV,IS)
      xx=dasum(nspln,wns(1,iv,is),1)
      CALL spline_inter(nspln,ggng,wns(1,iv,is),hg(1),m_wns(1),ngpw_l)
      IF(xx.gt.1.d-12) THEN
        DO ig=1,ngpw_l
          twnl(ig,iv,is)=ylmb(ig,lp,1)*vol*m_wns(ig)
        ENDDO
      ENDIF

      !IF(xx.gt.1.d-12) THEN
      !  DO ig=1,ngpw_l
      !    tw=wns(int(hg(ig)+1),2,iv,is)!curv2(hg(ig),nspln,ggng(1),wns(1,1,iv,is),wns(1,2,iv,is),0.0d0)
      !    twnl(ig,iv,is)=ylmb(ig,lp,1)*tw*vol
      !    
      !  ENDDO
      !ENDIF
    ENDDO
  ENDDO
  DEALLOCATE(ggng,m_wns,ylmb)
  END SUBROUTINE upf_nlpro

  SUBROUTINE upf_fnl(C0,nstat)
  USE max_parameter_pp
  USE system_data_types
  USE MYMPI
  USE constants
  USE math
  IMPLICIT NONE
  INTEGER   isa0,is,iv,ia,ig,isa,ikind,nat,ik,i,nstat,ii,NHXS
  REAL(KIND=DP)    er,ei,check,xg,xrg,flx,t,sum
  COMPLEX*16 ci,scr(ngpw_l,na_max),C0(NGPW_l,NSTAT)
  !REAL(KIND=DP), DIMENSION(:),ALLOCATABLE ::  SCRATCH  !
  !ALLOCATE(SCRATCH(NSTATE))
!  ALLOCATE(nl(1,atom_t,nghmax,nstat))
  ikind=1
  !CALL initzero(nl,nstat*atom_t)
 ALLOCATE(nl(1,atom_t,upf_nghmax,nstat))
  nl=0.0D0
  isa0=0
  DO is=1,sp_t
    DO iv=1,upf_ngh(is)
      ci=(0.0d0,-1.0d0)**UPF_nghtol(iv,is)
      IF(DABS(DREAL(ci)).gt.0.5d0) THEN
        DO ia=1,atom_p_sp(is)
          isa=isa0+ia
          DO ig=1,ngpw_l
            er=DREAL(eigr(ig,isa))
            ei=DIMAG(eigr(ig,isa))
            t=twnl(ig,iv,is)*(DREAL(ci))
            scr(ig,ia)=DCMPLX(t*er,t*ei)
      !      write(302,*)scr(ig,ia),ig,ia
          ENDDO
        ENDDO
      ELSE
        DO ia=1,atom_p_sp(is)
          isa=isa0+ia
          DO ig=1,ngpw_l
            er=DREAL(eigr(ig,isa))
            ei=DIMAG(eigr(ig,isa))
            t=twnl(ig,iv,is)*DIMAG(ci)
            scr(ig,ia)=DCMPLX(-t*ei,t*er)
       !     write(302,*)"*",scr(ig,ia),ig,ia
          ENDDO
        ENDDO
      ENDIF
      IF (G0_STAT)THEN
      DO ia=1,atom_p_sp(is)
        scr(1,ia)=0.5d0*scr(1,ia)
      ENDDO
      ENDIF
   !   print*,atom_p_sp(is),upf_nghmax,upf_ngh(is),size(nl),nstat
      IF(atom_p_sp(is).gt.1) THEN
        CALL dgemm('T','N',atom_p_sp(is),nstat,2*ngpw_l,2.d0,scr(1,1),&
      &          2*ngpw_l,c0(1,1),2*ngpw_l,0.0d0,&
      &          nl(1,isa0+1,iv,1),atom_t*upf_nghmax)
      ELSE
        CALL dgemv('t',2*ngpw_l,nstat,2.d0,c0(1,1),2*ngpw_l,scr(1,1),&
      &          1,0.0d0,nl(1,isa0+1,iv,1),atom_t*upf_nghmax)
      ENDIF
      !write(72,*)nl(1,isa0+1,iv,1),isa0
    ENDDO
    isa0=isa0+atom_p_sp(is)
   ENDDO
  CALL MPI_GlobSumR2(NL,nstat*atom_t*upf_NGHMAX)
!  print*,"NL",NL(1,1,1,1)
!  print*,"NL",NL(1,2,1,1)
!  print*,"NL",NL(1,1,2,1)
!  print*,"NL",NL(1,2,2,1)

!  print*,"NL",NL(1,1,3,1)
!  print*,"NL",NL(1,2,3,1)
!  print*,"NL",NL(1,1,4,1)
!  print*,"NL",NL(1,2,4,1)

!  print*,"NL",NL(1,1,5,1)
!  print*,"NL",NL(1,2,5,1)

  RETURN
  END SUBROUTINE upf_fnl
 SUBROUTINE nlppgrid(IS)
  USE max_parameter_pp
  USE system_data_types
  USE interpolgrid
  USE constants
  USE math
  USE bessel_func
  IMPLICIT NONE
  INTEGER   is,il,l,inghcom,j,iv,ir,lm,k,lold,ierr
  REAL(KIND=DP)::    gmax(splnmax),wsgtmp(splnmax)
  REAL(KIND=DP)::   wc(splnmax),wfc(splnmax),wtemp(splnmax)
  REAL(KIND=DP)::    ggng(1:nspln)      
  REAL(KIND=DP)::    jl(meshw(is))!,wns(nspln,2,1,sp_t)
  REAL(KIND=DP)::   xg,tmp,xgr
  REAL(KIND=DP)::   wfint(meshw(is)),fnt(splnmax),loggrid
  !ALLOCATE(wns(nspln,1,1,2))
fnt = 0.0d0
wfint = 0.0d0
gmax = 0.0d0
wsgtmp = 0.0d0
wc = 0.0d0
wfc = 0.0d0
wtemp = 0.0d0

  DO il=1,1!lread
    CALL intrp_grid(rr,meshw(is),rw(1,is),meshw(is),rps(1,is,il),splnmax)!,wc,wfc,wtemp)
  ENDDO

  CALL dcopy(meshv(is),vr(1,is,l_local(is)),1,gmax(1),1)

  IF(tkb(is)) THEN

    DO l=1,lmax(is)
      IF(l.eq.skip(is).or.l.eq.l_local(is)) THEN
        wsgtmp(l)=0.d0
      ELSE
        CALL dcopy(meshv(is),vr(1,is,l),1,gnl(1,is,l),1)
        DO ir=1,meshw(is)
          gnl(ir,is,l)=gnl(ir,is,l)-gmax(ir)
          IF(DABS(gnl(ir,is,l)).lt.1.d-10) gnl(ir,is,l) = 0.d0
        ENDDO

        DO ir=1,meshw(is)
          wfint(ir)=rps(ir,is,l)**2*gnl(ir,is,l)*rw(ir,is)
        ENDDO

        CALL simpsn(meshw(is),wfint,wsg(is,l))
        wsgtmp(l)=1.d0/((log(rr(meshw(is),is)/rr(1,is))/DBLE(meshw(is)-1))*wsg(is,l))
      ENDIF
    ENDDO
    k=0
    DO l=1,lmax(is)
      IF(l.ne.l_local(is).and.l.ne.skip(is)) THEN
        DO lm=1,2*l-1
          k=k+1
          wsg(is,k)=wsgtmp(l)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  DO il=1,nspln
            ggng(il)=(il-1)*(Gcutoff%pw)/DBLE(nspln-1)
  enddo
 IF(tkb(is)) THEN
    loggrid=(log(rr(meshw(is),is)/rr(1,is))/DBLE(meshw(is)-1))
    lold=-1
    DO iv=1,ngh(is)
      l=nghtol(iv,is)+1
      IF(l.eq.lold) THEN
      CALL dcopy(2*nspln,wns(1,iv-1,is),1,wns(1,iv,is),1)
      ELSE
      wns(:,iv,is)=0.0d0
      !CALL initzero(wns(1,1,iv,is),nspln)
        IF(l.ne.skip(is)) THEN
          DO il=1,nspln
       !     ggng(il)=(il-1)*(gcut_wfn)/DBLE(nspln-1)
            xg=sqrt(ggng(il))*twopi/a_lattice
            !CALL bess(xg,l,meshw(is),rw(1,is),jl)
           
            DO ir=1,meshw(is)
              xgr=xg*rw(ir,is)
              jl(ir)=sph_bes(l-1,xgr)
              fnt(ir)=rw(ir,is)**2*rps(ir,is,l)*gnl(ir,is,l)*jl(ir)
            ENDDO
           
            CALL simpsn(meshw(is),fnt,tmp)
            wns(il,iv,is)=loggrid*fourpi*tmp
          ENDDO
        ENDIF
      ENDIF
      lold=l
    ENDDO
  ENDIF
  RETURN

  END SUBROUTINE nlppgrid









END MODULE
