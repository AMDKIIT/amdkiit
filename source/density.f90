MODULE density

   CONTAINS

   SUBROUTINE eval_density  !!rhoofr.
   USE kinds
   USE system_data_types
   USE fft_interface
   USE constants
   IMPLICIT NONE
   COMPLEX(KIND=dp), DIMENSION(:), POINTER :: work
   INTEGER :: L,IG,IA,IB,IBB,istate,ID,IS1,IS2,NSTA,I
   integer :: ILSD1,ILSD2 !Sudhir DBG 
   REAL(kind=dp)     R1,R2,TT,COEF3,COEF4
   LOGICAL TFCAL
   ALLOCATE(WORK(MAX_FFT))!*NGROUP))
   rho=0.d0
   
      DO ID=1,NSTATE,2!*NGROUP
        NSTA=MIN((NSTATE-ID+2)/2,1)!NGROUP)
        DO I=ID,MIN(ID+2*NSTA-1,NSTATE)
          TFCAL=(OCCUPATION(I).NE.0.D0)
        ENDDO
        IF(TFCAL) THEN
          WORK=(0._dp,0._dp)
          DO IB=1,NSTA
            I=ID+2*(IB-1)
            IBB=(IB-1)*Ggridinfo(6,icpu)*nrlead(1)
            IS1=I
            IS2=I+1
            IF(IS2.GT.NSTATE) THEN
              DO IG=1,NGPW_l
                WORK(map_grid1d_p_pw(IG)+IBB)=C_0(IG,IS1)
                WORK(map_grid1d_m_pw(IG)+IBB)=DCONJG(C_0(IG,IS1))               
              ENDDO
            ELSE
              DO IG=1,NGPW_l
                WORK(map_grid1d_p_pw(IG)+IBB)=C_0(IG,IS1)+UIMAG*C_0(IG,IS2)
                WORK(map_grid1d_m_pw(IG)+IBB)=DCONJG(C_0(IG,IS1))+UIMAG*DCONJG(C_0(IG,IS2))
              ENDDO
            ENDIF
            
            IF(g0_stat.AND.(IS2.GT.NSTATE)) THEN
              WORK(map_grid1d_p_pw(1)+IBB)=C_0(1,IS1)
            ELSEIF(g0_stat) THEN
              WORK(map_grid1d_p_pw(1)+IBB)=C_0(1,IS1)+UIMAG*C_0(1,IS2)
            ENDIF
          ENDDO
            CALL fft_backward(WORK)

            IS1=ID
            IS2=ID+1
          COEF3=OCCUPATION(IS1)/cell_volume
          IF(IS2.GT.NSTATE) THEN
            COEF4=0.0D0
          ELSE
            COEF4=OCCUPATION(IS2)/cell_volume
          ENDIF
          IF(IS2.EQ.0) THEN
           COEF3=0.0D0
           COEF4=0.0D0
          ENDIF
!Sudhir DBG 
          IF(LOPEN_SHELL)THEN
             IF(IS1.LE.NEL_UP)ILSD1 = 1
             IF(IS2.LE.NEL_UP)ILSD2 = 1
             IF(IS1.GT.NEL_UP)ILSD1 = 2
             IF(IS2.GT.NEL_UP)ILSD2 = 2
!             print*,"**",ILSD1,ILSD2
             DO L=1,NNR1
                R1=DREAL(WORK(L))
                R2=DIMAG(WORK(L))
                !write(919,*)RHO(L,ILSD1),R1,RHO(L,ILSD2),R2
                RHO(L,ILSD1)=RHO(L,ILSD1)+COEF3*R1*R1
                RHO(L,ILSD2)=RHO(L,ILSD2)+COEF4*R2*R2
             
             ENDDO
          ELSE
            DO L=1,NNR1
              R1=DREAL(WORK(L))
              R2=DIMAG(WORK(L))
              TT=COEF3*R1*R1+COEF4*R2*R2
!              RHO(L)=RHO(L)+TT
              RHO(L,NLSD)=RHO(L,NLSD)+TT
              !write(6,*)L,R1,R2,COEF3,COEF4
            ENDDO 
          ENDIF !open shell  
!Sudhir DBG 
        ENDIF                   !endif TFCAL
      ENDDO                     !End loop over the electronic states

   DEALLOCATE(work)
   if(lOPEN_shell)then
   DO i=1,nnr1
   rho(i,1) = rho(i,1) + rho(i,2)
   enddo 
   endif
   END SUBROUTINE eval_density

   SUBROUTINE rhog2r
   USE kinds
   USE system_data_types
   USE fft_interface

   IMPLICIT NONE
   COMPLEX(KIND=dp), DIMENSION(:), POINTER :: work
   INTEGER IG,IR,ILSD

   ALLOCATE(WORK(MAX_FFT))

     !IF(.NOT.ASSOCIATED(RHO)) ALLOCATE(RHO(MAX_FFT,NLSD)) 
     DO ILSD=1,NLSD 
       WORK=(0.0d0,0.0d0)
       DO IG=1,NGRHO_l
         WORK(map_grid1d_m(IG))  = DCONJG(DCMPLX(RHO_G(IG,ILSD)))
         WORK(map_grid1d_P(IG))  = RHO_G(IG,ILSD)
       ENDDO
       IF(G0_stat) WORK(map_grid1d_P(1)) = RHO_G(1,ILSD)
       CALL fft_backward(WORK)
       DO IR=1,NNR1
         RHO(IR,ILSD)= DREAL(WORK(IR))
       ENDDO 
     ENDDO 
   DEALLOCATE(work)

   END SUBROUTINE rhog2r

   SUBROUTINE rhor2g
   USE kinds
   USE system_data_types
   USE fft_interface
   IMPLICIT NONE
   COMPLEX(KIND=dp), DIMENSION(:), POINTER :: work
   INTEGER IG,IR,ILSD 
   ALLOCATE(WORK(MAX_FFT))
!Sudhir DBG
   DO ILSD=1,NLSD  
     WORK=(0._dp,0._dp)
     DO IG=1,NNR1
       WORK(IG)  = DCMPLX(RHO(IG,ILSD),0.D0)
     ENDDO
     CALL fft_forward(WORK)
     DO IR=1,NGRHO_l
       RHO_G(IR,ILSD)= WORK(map_grid1d_p(IR))
     ENDDO
   ENDDO !lsd loop
   DEALLOCATE(work)
!Sudhir DBG 
   !if(icpu==1)write(*,*)"#1",rho_g(1),WORK(1),RHO(1),NNR1,NGRHO_l
   END SUBROUTINE rhor2g

!=-----------------------------------------------------------------------=!
   subroutine cubefile
   use system_data_types
   implicit none
   character*128  :: filename
   real*8 :: lowerleft(3),s(3),center(3) ! cube center i.e. geometric center  
   real*8,dimension(:,:),allocatable :: psi_tmp
   real*8,dimension(:,:,:),allocatable :: rho_tmp
   integer :: nris(3),is,ia,i,i_tot,ilowerleft(3),i1,i2,i3,ir,step,kr1s,kr2s,kr3s,j,k,ii,ii1,ii2
   integer, dimension(:),allocatable :: iatyp !atomic number of species 
!
   do i=1,nnr1
     RHO(i,1) = RHO(i,1)-RHO(i,2)!this sould go into PSI 
   enddo
!
   step=1 ! step=2 for halfmesh
   nris(1) = nrgrids(1) ! real space grid x,y,z 
   nris(2) = nrgrids(2)
   nris(3) = nrgrids(3)
!
   kr1s = nrlead(1) ! real space grid x+1,y+1,z+1
   kr2s = nrlead(2)
   kr3s = nrlead(3)
! 
   allocate(psi_tmp(kr2s,kr3s))
   allocate(rho_tmp(kr1s,kr2s,kr3s))
!
   ii=0
   do i=1,kr1s
     do j=1,kr1s
       do k=1,kr1s
         ii=ii+1
         rho_tmp(k,j,i) = rho(ii,1)
       enddo
     enddo
   enddo
  ii=0
   center(:) = 0.0d0
   i_tot = 0
   s(:) = 0.0d0
!
   do is=1,sp_t
     do ia=1,atom_p_sp(is)
        do i=1,3
           center(i) = center(i)+ATCO(i,IA,IS)
        enddo
        i_tot = i_tot + 1
     enddo
   enddo
!
   do i=1,3
     center(i)=center(i)/real(i_tot)
   enddo
!
  do i=1,3
     lowerleft(i) = center(i)-(0.5d0*(a_cell(1,i)+a_cell(2,i)+a_cell(3,i)))
     b_cell(1:3,i)=b_cell(1:3,i)/a_lattice
  enddo
!  s=matmul(transpose(b_cell),lowerleft)
  call DGEMV('T',3,3,1.0d0,b_cell,3,lowerleft,1,0d0,s,1)
  do i1=1,3
    ilowerleft(i1) = NINT(s(i1)*nris(i1))+1! generates 1..NR1S + n*NR1S
    lowerleft(i1) = 0d0
    do i2=1,3
      lowerleft(i1)=lowerleft(i1)+a_cell(i1,i2)*(ilowerleft(i2)-1)/DBLE(NRIS(i2))
     enddo
     ilowerleft(i1)=MOD((ilowerleft(i1)-1)+1000*NRIS(i1),NRIS(i1))+1
  enddo
  filename='rho_spin.cube'
   !open(200,file=filename)
   OPEN(200,FILE=filename,STATUS="unknown",position="APPEND")
   write(200,'(A)') ' CPMD CUBE FILE: '//trim(filename)
   write(200,'(A)') ' Total SCF Density'
   write(200,'(I5,3F12.6)')i_tot,(lowerleft(i),i=1,3)
   write(200,'(I5,3F12.6)')(nris(1)/step)+1,(dble(step)*a_cell(1,ir)/DBLE(nris(1)),ir=1,3)
   write(200,'(I5,3F12.6)')(nris(2)/step)+1,(dble(step)*a_cell(2,ir)/DBLE(nris(2)),ir=1,3)
   write(200,'(I5,3F12.6)')(nris(3)/step)+1,(dble(step)*a_cell(3,ir)/DBLE(nris(3)),ir=1,3)
   do is=1,sp_t
     do ia=1,atom_p_sp(is)
        write(200,'(I5,4F12.6)')z(is),dble(z(is)),(atco(ir,ia,is),ir=1,3)
      enddo
   enddo
   do i1=ilowerleft(1),nris(1)+ilowerleft(1)-1+step,step
      ii1=mod(i1-1,nris(1))+1
      call DCOPY(kr2s*kr3s,rho_tmp(ii1,1,1),kr1s, psi_tmp,1)
      do i2=ilowerleft(2),nris(2)+ilowerleft(2)-1+step,step
         ii2=mod(i2-1,nris(2))+1
         write(200,'(6E13.5)')(psi_tmp(ii2,mod(i3-1,nris(3))+1),i3=ilowerleft(3),nris(3)+ilowerleft(3)-1+step,step)
      enddo
   enddo
   close(200)
!
   do i=1,nnr1
     RHO(i,1) = RHO(i,1)+RHO(i,2) !back to alpha-density
   enddo
!
   deallocate(psi_tmp,rho_tmp)
   return
   end
END MODULE density
