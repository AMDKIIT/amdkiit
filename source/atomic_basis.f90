MODULE atomic_basis  

  CONTAINS
  SUBROUTINE setbasis
  USE kinds
  USE system_data_types
  USE constants
  USE interpolgrid
  USE math
  USE bessel_func
  !       SLATER MINIMAL BASIS
  integer ity!e_config(4,7,99)
  !
  integer,dimension(:,:),allocatable ::    nqsto
  real(kind=dp), dimension(:,:), allocatable::stoexp
  real(kind=dp), dimension(:,:,:), allocatable::atwfr
  real(kind=dp), dimension(:), allocatable:: ggnh,ggng,ar

  real(kind=dp)  atrg(splnmax,spmax),an,DATOM(nspln,2),ARHO(splnmax),WRK(NSPLN)
  REAL(kind=dp)  XG,FAAC,TMP,OCCU,CURRCHARGE!,ar(NGRHO_l)
  complex*16 work(2*splnmax,5)
  real(kind=dp)  fint(2*splnmax),bsint(2*splnmax),gg(2*splnmax),temp(splnmax,3)
  real(kind=dp)  disc,xmax,rmin,gmin
  integer is,mmax,nwork,n2,n22,il,ishell,l,ir,ierr,ne(4,7),meshat,isa,isa0
  real(kind=dp)    valch,zc,dclog!,srules
  integer   iatom,ms,m,n,i,MUSED
  logical saved
  real(kind=dp)  fac(0:14)
  data    fac / 1.d0, 1.d0, 2.d0, 6.d0, 24.d0, 120.d0, 720.d0,&
  &        5040.d0, 40320.d0, 362880.d0, 36288.d2, 399168.d2,&
  &        4.790016d8, 6.2270208d9, 8.7178291d10/


  real(kind=dp) :: r(2),cc
  integer :: istate, ig,lxx,natst,ish,ll,iv,ly,lpp(5),iaorb,iat,ia,ixx
  real(kind=dp) :: norm
  complex*16::ci,tsfac
  data       lpp /0,1,4,9,16/

  ALLOCATE(ggnh(nspln),ggng(nspln),ar(ngrho_l))
  ALLOCATE(RHO_G(NGRHO_l,NLSD))  
  ALLOCATE(RHO(MAX_FFT,NLSD)) 
  DO il=1,nspln
    ggng(il)=(il-1)*(Gcutoff%pw)/DBLE(nspln-1)
    ggnh(il)=(il-1)*(Gcutoff%rho)/DBLE(nspln-1)
  ENDDO

  ! Total number of valence electron
  tn_velec=0
  DO i = 1, sp_t
     tn_velec = tn_velec + atom_p_sp( i ) * zv( i )
  ENDDO
  !nelec = nelec - tot_charge_

  !Total number of valence orbital
   nstate=0
   DO is=1,sp_t
    if(lopen_shell)then
      nstate=nstate+idint(dble(atom_p_sp(is)*(zv(is))))
    else
       nstate=nstate+idint(dble(atom_p_sp(is)*(zv(is)))/2.d0)
     IF(2.d0*dble(nstate).lt.dble(atom_p_sp(is)*(zv(is)))) THEN
        nstate=nstate+1
        PRINT *,"check NSTATE calculation"
        EXIT
     ENDIF
    endif
   ENDDO
  !Occupation number
  ALLOCATE(OCCUPATION(nstate))
   if(lopen_shell)then
     do i=1,nstate
       occupation(i)=1.0d0
     enddo
   else
     DO i=1,nstate
         OCCUPATION(i)=2.d0
     ENDDO
  endif
!-----------------------------------------
 L=0
      DO i=1,sp_t
        DO J=1,atom_p_sp(i)
          DO K=1,3
              L=L+1
              DTBYM(L)=5.D0*5.D0/(atwt(z(i))*1822.888485D0)!DT2BYM(IS)!*LSKCOR(K,IAT)
          ENDDO
        ENDDO
      ENDDO


!----------------------------------------
    nghmax=1
upf_nghmax=1
ngh=0

  if(l_upf) then
  do is=1,sp_t
  upf_nghmax=Max(upf_ngh(is),upf_nghmax)
  end do
  else
   DO is=1,sp_t
    iv=0
    DO il=1,lmax(is)
     IF(il.ne.l_local(is).and.il.ne.skip(is)) THEN
       l=il-1
       ngh(is)=ngh(is)+(2*l+1)
       !nhxs=Max(ngh(is),1)
       nghmax=Max(ngh(is),nghmax)
       IF(l.eq.0) THEN
          inghcom=0
       ELSE
          inghcom=l*l
       ENDIF
        DO j=1,2*l+1
         iv=iv+1
         nghtol(iv,is)=l
         nghcom(iv,is)=inghcom+j
        ENDDO
     ENDIF
    ENDDO
   ENDDO
  endif
!-----------------------------------------
      NL_ANGMOM=0
      DO IS=1,sp_t        
          !NL_ANGMOM=NL_ANGMOM+atom_p_sp(IS)*NGH(IS)  
      ENDDO

  bsint=0.0
  fint=0.0
  DO is=1,sp_t
    call atom_info
    !ibtype(is)=1
    iatom=z(is)
    valch=zv(is)
    zc=DBLE(iatom)-valch
    CALL icopy(28,e_config(1,1,iatom),1,ne,1)

    ms=0
    m=0
    do n=1,7
      do l=1,4
        m=m+ne(l,n)
        if(m.gt.nint(zc).and.ne(l,n).ne.0)ms=ms+1
      enddo
    enddo
    nshell(is)=ms
  ENDDO
  meshat=256
  dclog=1.049999881d0
  m1shl=0
  do is=1,sp_t
      if(m1shl.lt.nshell(is)) m1shl=nshell(is)
  enddo

  allocate(nqsto(m1shl,spmax),stoexp(m1shl,spmax),oc(m1shl,spmax),atwfr(splnmax,m1shl,spmax),cat(nspln,2,m1shl,spmax),NUMAOR(sp_t))

  DO is=1,sp_t
   ne=0
    CALL icopy(28,e_config(1,1,z(is)),1,ne,1)
      zc=DBLE(z(is))-zv(is)
        ms=0
        m=0
        do n=1,7
          do l=1,4
            m=m+ne(l,n)
            if(m.gt.nint(zc).and.ne(l,n).ne.0) then
              ms=ms+1
              nqsto(ms,is)=n
              lshell(ms,is)=l-1
              stoexp(ms,is)=slat_expo(z(is),ne,n,l-1)
              oc(ms,is)=ne(l,n)
            endif
          enddo
        enddo

    atrg(1,is)=.7142857d-03
    do ir=2,meshat
      atrg(ir,is)=dclog*atrg(ir-1,is)
    enddo
    do i=1,nshell(is)
       l=lshell(i,is)
      an=sqrt((2.d0*(stoexp(i,is)))**(2*(nqsto(i,is))+1)/fac(2*(nqsto(i,is))))
      do ir=1,meshat
        if((-(stoexp(i,is)*atrg(ir,is))).lt.(-400.D0)) then
          atwfr(ir,i,is) = 0.d0
        else
         atwfr(ir,i,is)= an*atrg(ir,is)**(nqsto(i,is))*exp(-(stoexp(i,is))*atrg(ir,is))
        endif
      enddo
    enddo
  ENDDO



     NATTOT=0
     DO IS=1,sp_t
        NUMAOR(IS)=0
        DO I=1,NSHELL(IS)
          L=LSHELL(I,IS)
          NATTOT=NATTOT+atom_p_sp(IS)*(2*L+1)
          NUMAOR(IS)=NUMAOR(IS)+(2*L+1)
        ENDDO
        NUMAORMAX=MAX(NUMAORMAX,NUMAOR(IS))
     ENDDO

      ity=0
      cat=0.0
      do is=1,sp_t
        n2=nint(log(1.d0*meshat)/log(2.d0)+0.499999d0)
        n22=2**n2
        rmin=log(atrg(1,is))

        if(ity.eq.1) then
          gmin=log(sqrt(Gcutoff%rho)*twopibya)-(meshat-1)*dlog(1.049999881D0)
        else
          gmin=log(sqrt(Gcutoff%pw)*twopibya)-(meshat-1)*dlog(1.049999881D0)
        endif
        do il=1,meshat
          gg(il)=(exp(gmin+(il-1)*dlog(1.049999881D0))/(twopibya))**2
        enddo
        do ishell=1,nshell(is)
          saved=.false.
          l=lshell(ishell,is)
          fint=0.0
          do ir=1,meshat
            fint(ir)=fourpi*atwfr(ir,ishell,is)/atrg(ir,is)
          enddo
!         fourier transformation
          call sbt(fint,bsint,l,rmin,gmin,dlog(1.049999881D0),n2,saved,&
          &work(1,1),work(1,2),work(1,3),work(1,5),2*splnmax,disc)
            call INTRP_grid(gg(1),meshat,ggng(1),nspln,bsint(1),splnmax)
            call dcopy(nspln,bsint(1),1,cat(1,1,ishell,is),1)
            if(l.gt.0.and.ggng(1).lt.1.d-12) cat(1,1,ishell,is)=0.0d0
             
        enddo
      enddo

!!!------------------------------------------------------------
    ISA0=0
    Rho_G=(0._dp,0._dp)
    DO IS=1,sp_t
      CURRCHARGE=0._dp
      ARHO=0._dp
        DO ISH=NSHELL(IS),1,-1
          IF(CURRCHARGE.GT.OC(ISH,IS)) THEN
            OCCU=0
            CURRCHARGE=CURRCHARGE-OC(ISH,IS)
          ELSE
            OCCU=OC(ISH,IS)-CURRCHARGE
            CURRCHARGE=0
          ENDIF
            DO IR=1,MESHAT
            ARHO(IR)=ARHO(IR)+OCCU*(ATWFR(IR,ISH,IS)/ATRG(IR,IS))**2
            ENDDO
        ENDDO
      XG=0.D0
        DO IR=MESHAT,1,-1
          XG=XG+ABS(ARHO(IR))
          IF(XG.GT.1.D-8) THEN
          MUSED = IR
          GOTO 100
          ENDIF
        ENDDO

      MUSED=MMAX
      100  CONTINUE
      DATOM=0._dp
      DO IL=1,NSPLN
        XG=SQRT(GGNH(IL))*twopibya
        CALL ov_sph_bes_f(ATRG(1,IS),DLOG(1.049999881D0),MUSED,ARHO,0,XG,ATRG(MUSED,IS),TMP)
        DATOM(IL,1)=TMP
      ENDDO
       !
        CALL ov_sph_bes_f(ATRG(1,IS),DLOG(1.049999881D0),MUSED,ARHO,0,0.0D0,ATRG(MUSED,IS),TMP)

      IF(TMP.LT.1.D-3) THEN
        FAAC=0.D0
      ELSE
        FAAC=(ZV(IS))/TMP!-ATCHG(IS))/TMP
      ENDIF

      CALL DSCAL(NSPLN,FAAC,DATOM(1,1),1)
      CALL spline_inter(NSPLN,GGNH(1),DATOM(1,1),HG(1),AR(1),NGRHO_l)

        DO IG=1,NGRHO_l
          TSFAC=DCMPLX(0.D0,0.D0)
          DO IA=1,atom_p_sp(IS)
            ISA=ISA0+IA
            TSFAC=TSFAC+eigr(ig,isa)
          ENDDO
          Rho_G(IG,NLSD)=Rho_G(IG,NLSD) + AR(IG)*(1/cell_volume)*TSFAC
        ENDDO
        ISA0=ISA0+atom_p_sp(IS)
    ENDDO !sp_t loop

DEALLOCATE(ggnh,ggng,nqsto,stoexp,ar,atwfr)
RETURN
!!!------------------------------------------------------------
   END SUBROUTINE setbasis

     function slat_expo(atmic_num,ne_shell,n,l)
!     slater exponents from the clementi raimondi table
      USE kinds
      implicit none
      real(kind=dp) slat_expo 
      integer atmic_num,ne_shell(4,7),n,l
      integer l_1,l_2,m,m_1,m_2,i
      real(kind=dp)  effec_ncharge,q_num(7),table_cle(2,3,18)
      data    table_cle  /1.0d0,5*0.0d0, 1.6875d0,5*0.0d0,&
          2.6906d0,0.0d0,0.6396d0,3*0.0d0,&
          3.6848d0,0.0d0,0.9560d0,3*0.0d0,&
          4.6795d0,0.0d0,1.2881d0,1.2107d0,2*0.0d0,&
          5.6727d0,0.0d0,1.6083d0,1.5679d0,2*0.0d0,&
          6.6651d0,0.0d0,1.9237d0,1.9170d0,2*0.0d0,&
          7.6579d0,0.0d0,2.2458d0,2.2266d0,2*0.0d0,&
          8.6501d0,0.0d0,2.5638d0,2.5500d0,2*0.0d0,&
          9.6421d0,0.0d0,2.8792d0,2.8792d0,2*0.0d0,&
          10.6259d0,0.d0,3.2857d0,3.4009d0,0.8359d0,0.d0,&
          11.6089d0,0.0d0,3.6960d0,3.9129d0,1.1025d0,0.0d0,&
          12.5910d0,0.0d0,4.1068d0,4.4817d0,1.3724d0,1.3552d0,&
          13.5724d0,0.0d0,5.5100d0,4.9725d0,1.6344d0,1.4284d0,&
          14.5578d0,0.0d0,4.9125d0,5.4806d0,1.8806d0,1.6288d0,&
          15.5409d0,0.0d0,5.3144d0,5.9885d0,2.1223d0,1.8273d0,&
          16.5239d0,0.0d0,5.7152d0,6.4966d0,2.3561d0,2.0387d0,&
          17.5075d0,0.0d0,6.1152d0,7.0041d0,2.5856d0,2.2547d0/
      data    q_num  /1.0d0,2.0d0,3.0d0,3.7d0,4.0d0,4.2d0,4.4d0/


      if(atmic_num.le.18.and.l.le.1) then
        slat_expo=table_cle(l+1,n,atmic_num)
        return
      endif
!     calculate the shielding
      effec_ncharge=0.0d0
!     the complete shell
      l_1=l+1
      if(l_1.eq.1) l_2=2
      if(l_1.eq.2) l_2=1
      if(l_1.eq.3) l_2=4
      if(l_1.eq.4) l_2=3

      if(n.eq.1) then
        m=ne_shell(1,1)
        effec_ncharge=effec_ncharge+0.3d0*(m-1)
      else
        m=ne_shell(l_1,n)+ne_shell(l_2,n)
        effec_ncharge=effec_ncharge+0.35d0*(m-1)
      endif
     !S = (0.35 * n_1^2) + (0.85 * n_2^2) + (1 * n_3^2)
      if(l_1+l_2.eq.3) then
        if(n.gt.1) then
          m_1=ne_shell(1,n-1)+ne_shell(2,n-1)+ne_shell(3,n-1)+ne_shell(4,n-1)
          m_2=0
          do i=1,n-2
            m_2=m_2+ne_shell(1,i)+ne_shell(2,i)+ne_shell(3,i)+ne_shell(4,i)
          enddo
          effec_ncharge=effec_ncharge+0.85d0*m_1+1.0d0*m_2
        endif
      else

        m=0
        do i=1,n-1
          m=m+ne_shell(1,i)+ne_shell(2,i)+ne_shell(3,i)+ne_shell(4,i)
        enddo
        effec_ncharge=effec_ncharge+1.0d0*m
      endif

      slat_expo = (dble(atmic_num) - effec_ncharge)/q_num(n)
!     ==--------------------------------------------------------------==
      return
      end
 END MODULE atomic_basis   
