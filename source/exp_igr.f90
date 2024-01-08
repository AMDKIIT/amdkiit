MOdule exp_igr
  CONTAINS

  SUBROUTINE exp_igrxfactor
  USE max_parameter_pp
  USE system_data_types
  USE kinds
  USE constants

  IMPLICIT NONE
  INTEGER    nh1,nh2,nh3,isa,is,ia,i,j,k,ig,natx
  COMPLEX*16, dimension(:,:),allocatable ::ei1,ei2,ei3
  COMPLEX*16 dp1,dp2,dp3,dm1,dm2,dm3,ei10,ei20,ei30,tm_dp,tm_dm,ei123
  REAL(KIND=dp)    sum1,sum2,sum3,sum_1,sum_2,sum_3
  INTEGER nr1s,nr2s,nr3s,iat
  REAL(kind=dp) b1(3),b2(3),b3(3)

  ALLOCATE(eigr(ngrho_l,atom_t),eigrxrhos(ngrho_l),eigrxvps(ngrho_l),eigr_pw(ngpw_l,atom_t))

  nr1s=nrgrids(1)
  nr2s=nrgrids(2)
  nr3s=nrgrids(3)

  nh1=nr1s/2+1
  nh2=nr2s/2+1
  nh3=nr3s/2+1
      
  natx=2*(atom_t/2) + 1

  ALLOCATE(ei1(natx,(2*nr1s-1)))
  ALLOCATE(ei2(natx,(2*nr2s-1)))
  ALLOCATE(ei3(natx,(2*nr3s-1)))
  ALLOCATE(iatpt(2,atom_t))

  ei1(:,:)=(0._dp,0._dp)
  ei2(:,:)=(0._dp,0._dp)
  ei3(:,:)=(0._dp,0._dp)

  eigrxrhos(:)=(0._dp,0._dp)
  eigrxvps(:)=(0._dp,0._dp)
  iat=0
  DO is=1,sp_t
    DO ia=1,atom_p_sp(is)
      iat=iat+1
      iatpt(1,iat)=ia
      iatpt(2,iat)=is
    ENDDO
  ENDDO

  DO isa =1,atom_t
    ia=iatpt(1,isa)
    is=iatpt(2,isa)
        sum1=b_cell(1,1)*atco(1,ia,is)+b_cell(2,1)*atco(2,ia,is)+&
     &       b_cell(3,1)*atco(3,ia,is)
        sum2=b_cell(1,2)*atco(1,ia,is)+b_cell(2,2)*atco(2,ia,is)+&
     &       b_cell(3,2)*atco(3,ia,is)
        sum3=b_cell(1,3)*atco(1,ia,is)+b_cell(2,3)*atco(2,ia,is)+&
     &       b_cell(3,3)*atco(3,ia,is)


    sum_1=twopibya*sum1
    sum_2=twopibya*sum2
    sum_3=twopibya*sum3

    ei1(isa,1)=DCMPLX(1.0d0,0.0d0)
    ei2(isa,1)=DCMPLX(1.0d0,0.0d0)
    ei3(isa,1)=DCMPLX(1.0d0,0.0d0)
    dp1=DCMPLX(DCOS(sum_1),-DSIN(sum_1))
    dp2=DCMPLX(DCOS(sum_2),-DSIN(sum_2))
    dp3=DCMPLX(DCOS(sum_3),-DSIN(sum_3))
    dm1=DCONJG(dp1)
    dm2=DCONJG(dp2)
    dm3=DCONJG(dp3)

    tm_dp=dp1
    tm_dm=dm1
    DO i=2,nr1s
      ei1(isa,i)=tm_dp
      tm_dp=tm_dp*dp1
      ei1(isa,nr1s+i-1)=tm_dm
      tm_dm=tm_dm*dm1
    ENDDO

    tm_dp=dp2
    tm_dm=dm2
    DO j=2,nr2s
      ei2(isa,j)=tm_dp
      tm_dp=tm_dp*dp2
      ei2(isa,nr2s+j-1)=tm_dm
      tm_dm=tm_dm*dm2
    ENDDO

    tm_dp=dp3
    tm_dm=dm3
    DO k=2,nr3s
      ei3(isa,k)=tm_dp
      tm_dp=tm_dp*dp3
      ei3(isa,nr3s+k-1)=tm_dm
      tm_dm=tm_dm*dm3
    ENDDO

    ei10=1.d0/ei1(isa,nh1)
    ei20=1.d0/ei2(isa,nh2)
    ei30=1.d0/ei3(isa,nh3)

    DO i=1,2*nr1s-1
      ei1(isa,i)=ei1(isa,i)*ei10
    ENDDO
    DO j=1,2*nr2s-1
      ei2(isa,j)=ei2(isa,j)*ei20
    ENDDO
    DO k=1,2*nr3s-1
      ei3(isa,k)=ei3(isa,k)*ei30
    ENDDO

  ENDDO

  DO isa=1,atom_t
    DO ig=1,ngrho_l
      eigr(ig,isa)=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*ei3(isa,inyh(3,ig))
    ENDDO
  ENDDO
   DO isa=1,atom_t
    DO ig=1,ngpw_l
      eigr_pw(ig,isa)=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*ei3(isa,inyh(3,ig))
    ENDDO
  ENDDO


  DO ig=1,ngrho_l
    DO isa=1,atom_t
      is=iatpt(2,isa)
      eigrxvps(ig)=eigrxvps(ig)+DCMPLX(DREAL(eigr(ig,isa))*vps(is,ig),DIMAG(eigr(ig,isa))*vps(is,ig))
      eigrxrhos(ig)=eigrxrhos(ig)+rhos(is,ig)*eigr(ig,isa)
   ENDDO
  ENDDO

  
  DEALLOCATE(ei1,ei2,ei3,iatpt)
END SUBROUTINE exp_igrxfactor

  SUBROUTINE form_factor_upf
  USE max_parameter_pp
  USE system_data_types
  USE interpolgrid
  USE constants
  USE math

  IMPLICIT NONE

  INTEGER   is,ir,il,ierr,ig,i,isa
  REAL(kind=dp)    fint(splnmax), dfint(splnmax)
  REAL(kind=dp)    xg,xrg,flx
  REAL(KIND=dp),   DIMENSION(:,:),       POINTER ::voo
  REAL(KIND=dp),   DIMENSION(:),       POINTER ::vo
  REAL(KIND=dp),   DIMENSION(:),       POINTER :: ggnh

  ALLOCATE(VPS(SPMAX,NGRHO_L),rhos(SPMAX,NGRHO_L))
  ALLOCATE(voo(1:nspln,1:sp_t))
  ALLOCATE(vo(ngrho_l),ggnh(nspln))


  DO il=1,nspln
      ggnh(il)=(il-1)*(Gcutoff%rho/DBLE(nspln-1))  !point grid from 0 to gcut
  ENDDO

  DO is=1,sp_t
    fint=0.0D0
    DO ir=1,meshv(is)
       fint(ir)=rr(ir,is)*rr(ir,is)*vr(ir,is,1)+zv(is)*derf(rr(ir,is)/raggio)*rr(ir,is)
    ENDDO

    voo=0.0D0 
    DO il=1,nspln
      !xg=sqrt(ggnh(il))*twopibya
      
      IF((sqrt(ggnh(il))*twopibya).gt.1.d-6) THEN
        DO ir=1,meshv(is)
          xrg = rr(ir,is)*sqrt(ggnh(il))*twopibya
          dfint (ir) = fint(ir)* SIN(xrg)/xrg*rw(ir,is)
        ENDDO
      ELSE
         DO IR=1,MESHV(IS)
            DFINT (IR) = FINT(IR) * RW(IR,IS)
         ENDDO
      ENDIF
      CALL simpsn (meshv(is), dfint(1), flx)
      VOO(IL,IS) = fourpi* FLX
    ENDDO
    CALL spline_inter(nspln,ggnh,voo(1,is),hg(1),vo(1),ngrho_l)
    DO ig=1,ngrho_l
      vps(is,ig)=(1.d0/cell_volume)*vo(ig)
    ENDDO

   ! CALL spline_inter(nspln,ggnh,voo(1,is),hg(1),voo(1,2,is),ngrho_l)
   ! DO ig=1,ngrho_l
   !   vps(is,ig)=(1.d0/cell_volume)*voo(int(hg(ig)+1),2,is)
   ! ENDDO
  ENDDO


   DO is=1,sp_t
    DO ig=1,ngrho_l
      rhos(is,ig)=-zv(is)*exp(-(0.25D0*(raggio*raggio)*hg(ig)*twopibya2))*(1.d0/cell_volume)
    ENDDO
   if(g0_stat) rhos(is,1)=-zv(is)*(1.d0/cell_volume)
  ENDDO
  RETURN
  END SUBROUTINE form_factor_upf

  SUBROUTINE form_factor
  USE max_parameter_pp
  USE system_data_types
  USE interpolgrid                                                                               
  USE constants                                                                             
  USE math

  IMPLICIT NONE
  INTEGER   is,ir,il,ierr,ig,i,isa
  REAL(kind=dp)    fint(splnmax), dfint(splnmax)
  REAL(kind=dp)    check,xg,xrg,flx
  REAL(KIND=dp),   DIMENSION(:,:),       POINTER ::voo
  REAL(KIND=dp),   DIMENSION(:),       POINTER ::vo
  REAL(KIND=dp),   DIMENSION(:),       POINTER :: ggnh
  LOGICAL  :: zer

  ALLOCATE(VPS(SPMAX,NGRHO_L),rhos(SPMAX,NGRHO_L))
  ALLOCATE(voo(1:nspln,1:sp_t))
  ALLOCATE(vo(ngrho_l),ggnh(nspln))

  DO is=1,sp_t
  fint=0.0D0
  !construct logarithmic grid from rr(1) to rr(mesh) with mesh points in rw(1:mesh).
  rw(1,is)=rr(1,is)
    DO ir=2,meshv(is)
      rw(ir,is)=rw((ir-1),is)*EXP(LOG(rr(meshv(is),is)/rr(1,is))/DBLE(meshv(is)-1))
    ENDDO

    DO il=1,lmax(is)
      CALL intrp_grid(rr(1,is),meshv(is),rw(1,is),meshv(is),vr(1,is,il),splnmax)
    ENDDO
    
    DO ir=1,meshv(is)
      fint(ir)=(rw(ir,is)*rw(ir,is)*rw(ir,is)*vr(ir,is,l_local(is))+&
      &zv(is)*derf(rw(ir,is)/raggio)*rw(ir,is)*rw(ir,is))
      !  the width of the ionic charge distribution ,raggio,default value 1.2
      zer = DABS(rw(ir,is)).lt.1.d-8
      IF(.not.zer) THEN
        check = vr(ir,is,l_local(is))+zv(is)*derf(rw(ir,is)/raggio)/rw(ir,is)
      ELSE
        check = vr(ir,is,l_local(is))+2.d0*zv(is)/(DSQRT(pi)*raggio)
      ENDIF
      IF(DABS(check).lt.1.d-8) fint(ir)=0.d0
    ENDDO
    
    voo=0.0D0
    DO il=1,nspln
      ggnh(il)=(il-1)*(Gcutoff%rho)/DBLE(nspln-1)
      xg=sqrt(ggnh(il))*twopibya


      IF(xg.gt.1.d-6) THEN
        DO ir=1,meshv(is)
          xrg = rw(ir,is)*xg 
          dfint (ir) = fint(ir)* SIN(xrg)/xrg
        ENDDO
      ELSE
        CALL dcopy(meshv(is),fint(1),1,dfint(1),1)
      ENDIF
      CALL simpsn (meshv(is), dfint(1), flx)
      voo(il,is) = (LOG(rr(meshv(is),is)/rr(1,is))/DBLE(meshv(is)-1))*fourpi* flx
    ENDDO
    CALL spline_inter(nspln,ggnh,voo(1,is),hg(1),vo(1),ngrho_l)
    DO ig=1,ngrho_l
      vps(is,ig)=(1.d0/cell_volume)*vo(ig)
    ENDDO
  ENDDO !loop over speciese

  
  DO is=1,sp_t
    DO ig=1,ngrho_l
      rhos(is,ig)=-zv(is)*exp(-(0.25D0*(raggio*raggio)*hg(ig)*twopibya2))*(1.d0/cell_volume)
    ENDDO
  IF(g0_stat) rhos(is,1)=-zv(is)*(1.d0/cell_volume)
  ENDDO
  DEALLOCATE(VOO,vo,ggnh)
  RETURN
  END SUBROUTINE form_factor

  END MOdule exp_igr
