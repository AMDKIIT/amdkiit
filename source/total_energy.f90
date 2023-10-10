MODULE total_energy

  CONTAINS
  SUBROUTINE e_total
  USE kinds
  USE system_data_types
  USE constants
  USE mympi

  IMPLICIT NONE
  INTEGER i,is,ig,isa
  COMPLEX(KIND=dp) :: tmp
   
  CALL e_hartree
  CALL e_sr
  estat= (DREAL(eh)*cell_volume) + esr-eself!/DBLE(ncpu)!nproc)
  CALL e_ke
  CALL e_xc   
! if(icpu==1)then
  etotal=enl+eloc+eke+exc+estat
! write(6,*)eke,etotal!,enl,dreal(eloc),eke,icpu

!if(IONODE.AND.DEBUG)then
!OPEN(7,FILE="energy.dat",STATUS="unknown",position="APPEND")
!  WRITE(7,'(1A,A30,A2,T38,F15.10,A5)'),"#",'TOTAL ENERGY',"=",etotal,"A.U."!ehatree
!  WRITE(7,'(1A,A30,A2,T38,F15.10,A5)')"#",'KINETIC ENERGY',"=",eke,"A.U."
!  WRITE(7,'(1A,A30,A2,T38,F15.10,A5)')"#",'ELECTROSTATIC ENERGY',"=",estat,"A.U."
!  WRITE(7,'(1A,A30,A2,T38,F15.10,A5)'),"#",'SELF ENERGY CONTRIBUTION',"=",dreal(eself),"A.U."
!  WRITE(7,'(1A,A30,A2,T38,F15.10,A5)'),"#",'OVERLAP ENERGY CONTRIBUTION',"=",esr,"A.U."
!  WRITE(7,'(1A,A30,A2,T38,F15.10,A5)')"#",'LOCAL PP ENERGY',"=",dreal(eloc),"A.U."
!  WRITE(7,'(1A,A30,A2,T38,F15.10,A5)')"#",'NON LOCAL PP ENERGY',"=",enl,"A.U."
!  WRITE(7,'(1A,A30,A2,T38,F15.10,A5)')"#",'EXCHANGE CORRELATION ENERGY',"=",exc,"A.U."
!  WRITE(7,*) "********"
! CLOSE(7)
!endif
  ENd SUBROUTINE

  SUBROUTINE e_ke
  USE kinds
  USE system_data_types
  USE constants
  use mympi
  IMPLICIT NONE
  INTEGER i,is,ig,isa
  REAL(kind=dp) :: ekin
  COMPLEX*16 :: tmp
 
  !write(6,*)"kinenergy",c_0(1,1)
  eke=0._dp
  DO i=1,nstate
     ekin=0._dp
    DO ig=1,ngpw_l
      tmp=c_0(ig,i)
      ekin=ekin+hg(ig)*DREAL(DCONJG(tmp)*tmp)
    END DO
    
    eke=eke+ekin*occupation(i)
  END DO
  eke=eke*twopibya2
  CALL MPI_GlobSumR2s(EKE)
  END SUBROUTINE


  SUBROUTINE e_hartree
  USE kinds
  USE system_data_types
  USE constants
  USE fft_interface
  USE MYMPI
  IMPLICIT NONE
  INTEGER is,ig,ig1
  COMPLEX*16 rhog
  COMPLEX(KIND=dp), DIMENSION(:), POINTER :: work

!Sudhir DBG 
    if(lopen_shell.and.nlsd.eq.2)then
       do ig=1,ngrho_l
          !rho_g(ig,1)=rho_g(ig,1)+rho_g(ig,2)
       enddo
     endif
!Sudhir DBG 

  if(g0_stat) then
    ig1=2
    rhog=rho_g(1,1)+eigrxrhos(1) !Sudhir DBG 
    eh=(0.0d0,0.0d0)!
    eion=(0.0d0,0.0d0)!
    ee=(0.0d0,0.0d0)! 
    elpp=0.5d0*eigrxvps(1)*DCONJG(rho_g(1,1)) !Sudhir DBG
  else
    ig1=1
    eh=(0.0d0,0.0d0)
    eion=(0.0d0,0.0d0)
    ee=(0.0d0,0.0d0)
    elpp=(0.0d0,0.0d0)
  endif
  !print*,(DCMPLX(fourpi/(twopibya2*hg(ig)),0.D0)),rhog(1,1)
  do ig=ig1,ngrho_l
    rhog=rho_g(ig,1)+eigrxrhos(ig)
    eh=eh+((DCMPLX(fourpi/(twopibya2*hg(ig)),0.D0))*rhog)*dconjg(rhog)
    eion=eion+(DCMPLX(fourpi/(twopibya2*hg(ig)),0.D0))*eigrxrhos(ig)*dconjg(eigrxrhos(ig))
    ee=ee+(DCMPLX(fourpi/(twopibya2*hg(ig)),0.D0))*rho_g(ig,1)*DCONJG(rho_g(ig,1)) 
    elpp=elpp+eigrxvps(ig)*DCONJG(rho_g(ig,1)) 
  enddo
  !

  CALL MPI_GlobSumC2s(elpp)
  CALL MPI_GlobSumC2s(eh)

  eself=0.d0
  do is=1,sp_t
    eself=eself+dble(atom_p_sp(is))*zv(is)*zv(is)/raggio
  enddo
    eself=eself/dsqrt(2.d0*pi)

 ! Correction for charged systems
  if(symmetry.NE.0) then
    do is=1,sp_t
      eself=eself+(charge*pi/cell_volume)*dble(atom_p_sp(is))*zv(is)*raggio**2
    enddo
  endif
!  endif
    !e_hat= (DREAL(eh)*cell_volume) + esr+eself/DBLE(ncpu)!nproc)
!Sudhir DBG 
    if(lopen_shell.and.nlsd.eq.2)then
       do ig=1,ngrho_l
          !rho_g(ig,1)=rho_g(ig,1)-rho_g(ig,2)
       enddo
     endif
!Sudhir DBG 
  END SUBROUTINE
  
  SUBROUTINE e_sr
  USE kinds
  USE system_data_types
  USE constants
  USE mympi
  IMPLICIT NONE
  REAL(KIND=dp)    ARGMAX          !DERFC(ARG) and DEXP(-ARG*ARG) underflow
  PARAMETER (ARGMAX=20.D0)
  REAL(KIND=dp)  RXLM1,RXLM2,RXLM3,ZV2,RCKJ,YLM,ZLM,ERRE2,RLM,ARG,ESRTZERO,ADDESR,ADDPRE,&
  &        REPAND,XLM,xlm_p,ylm_p,zlm_p
  REAL(KIND=dp)  TFION1,TFION2,TFION3,TFION4,TFION5,TFION6
  INTEGER ISUB,K,J,L,M,LAX,INF,IX,IY,IZ,ISHFT,IAT
  LOGICAL TZERO
  IF(.NOT.ASSOCIATED(force)) ALLOCATE(force(3,NA_MAX,sp_t))
  force=0._dp
  esr=0._dp
  iat=0
  DO k=1,sp_t
    DO j=k,sp_t
      zv2=zv(k)*zv(j)
      IF(ABS(zv2).lt.1.d-10) GOTO 2000
      rckj=DSQRT(raggio*raggio+raggio*raggio)
      lax=atom_p_sp(k)
        DO l=1,lax
         IF(IATPE(IAT+L).NE.icpu) GOTO 1000
          inf=1
          IF(k.eq.j)inf=l
            DO m=inf,atom_p_sp(j)

              IF(l.eq.m.and.k.eq.j) THEN
                xlm=0.d0
                ylm=0.d0
                zlm=0.d0
                tzero=.true.
              ELSE
                tzero=.false.
                xlm_p=atco(1,l,k)-atco(1,m,j)
                ylm_p=atco(2,l,k)-atco(2,m,j)
                zlm_p=atco(3,l,k)-atco(3,m,j)
                !print*,atco(1,l,k),atco(1,m,j)
                !CALL apply_pbc(xlm,ylm,zlm,xlm,ylm,zlm,1,a_lattice,symmetry)
                CALL apply_pbc(xlm_p,ylm_p,zlm_p,xlm,ylm,zlm,1,a_lattice,symmetry) 
              ENDIF
              IF(TFOR) THEN
                TFION1 = force(1,L,K)
                TFION2 = force(2,L,K)
                TFION3 = force(3,L,K)
                TFION4 = force(1,M,J)
                TFION5 = force(2,M,J)
                TFION6 = force(3,M,J)
              ENDIF
              DO ix=-nlp,nlp
                DO iy=-nlp,nlp
                  DO iz=-nlp,nlp
                    ishft=ix*ix+iy*iy+iz*iz
                    IF(.not.(tzero.and.ishft.eq.0)) THEN
                      rxlm1=xlm+ix*a_cell(1,1)+iy*a_cell(2,1)+iz*a_cell(3,1)
                      rxlm2=ylm+ix*a_cell(1,2)+iy*a_cell(2,2)+iz*a_cell(3,2)
                      rxlm3=zlm+ix*a_cell(1,3)+iy*a_cell(2,3)+iz*a_cell(3,3)
                      erre2=rxlm1*rxlm1+rxlm2*rxlm2+rxlm3*rxlm3
                      rlm=DSQRT(erre2)
                      arg=rlm/rckj
                     !write(50,*)"*",RXLM1,XLM,nlp,a_cell(1,1),a_cell(2,1),a_cell(3,1)
                     IF(arg.le.argmax) THEN !addesr,addpre /= 0
                       IF(tzero) THEN
                         esrtzero=0.5d0
                       ELSE
                         esrtzero=1.d0
                       ENDIF
                       addesr=zv2*derfc(arg)/rlm
                       esr=esr+addesr*esrtzero
                       IF(TFOR) THEN
                         ADDPRE=(2.D0*ZV2/DSQRT(PI))*DEXP(-ARG*ARG)/RCKJ
                         REPAND=ESRTZERO*(ADDESR+ADDPRE)/ERRE2
                         TFION1=TFION1+REPAND*RXLM1
                         TFION2=TFION2+REPAND*RXLM2
                         TFION3=TFION3+REPAND*RXLM3
                         TFION4=TFION4-REPAND*RXLM1
                         TFION5=TFION5-REPAND*RXLM2
                         TFION6=TFION6-REPAND*RXLM3
                      ENDIF
                     ENDIF
                     ENDIF
                  ENDDO  !ixc
                ENDDO    !iyc
              ENDDO      !izc
              IF(TFOR) THEN
                force(1,L,K) = TFION1
                force(2,L,K) = TFION2
                force(3,L,K) = TFION3
                force(1,M,J) = TFION4
                force(2,M,J) = TFION5
                force(3,M,J) = TFION6
              ENDIF
            ENDDO
         1000 CONTINUE
        ENDDO
       
      2000 CONTINUE
    ENDDO
    iat=iat+atom_p_sp(k)
  ENDDO
          !write(29,*)force(1,2,2),force(2,2,2),force(3,2,2)
          !write(29,*)"  **infw_esr  "

  CALL MPI_GlobSumR2s(Esr)
  END SUBROUTINE

  
  SUBROUTINE e_xc
  USE kinds
  USE system_data_types
  USE constants
  !USE functional
  USE potential
  USE XC
  USE gvectors
  use xc_f90_lib_m
  use mympi
  IMPLICIT NONE


  integer :: i!, vmajor, vminor, vmicro, func_id = 20

  INTEGER nr1,nr2,nr3
      

     nr1=nrgrids(1)
     nr2=nrgrids(2)
     nr3=nrgrids(3)
     exc = cell_volume/dble(nr1*nr2*nr3) * etxc
  END SUBROUTINE
 

END MODULE
