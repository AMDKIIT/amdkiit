  module potential
  contains

  SUBROUTINE eval_pot
  USE kinds
  USE system_data_types
  USE constants
  USE fft_interface

  IMPLICIT NONE
  INTEGER is,ig,ig1,ILSD 
  COMPLEX*16 VTEMP(NGRHO_l,NLSD)  
  COMPLEX(KIND=dp), DIMENSION(:), POINTER :: work

  ALLOCATE(WORK(MAX_FFT))

  ALLOCATE(V_POT(MAX_FFT,NLSD))
  V_POT=(0.0d0,0.0d0)

  WORK=(0.0d0,0.0d0)
  if(g0_stat) then
    ig1=2
    vtemp(1,1)=0.0d0
    !WORK(map_grid1d_m(1))  = (0.0D0,0.0D0)
    !WORK(map_grid1d_P(1))  = (0.0D0,0.0D0)
  else
    ig1=1
  endif
 
    DO ig=ig1,ngrho_l
    vtemp(ig,1) =((DCMPLX(fourpi/(twopibya2*hg(ig)),0.D0))*(rho_g(ig,1)+eigrxrhos(ig)))+eigrxvps(ig)
    ENDDO

  DO ig=1,ngrho_l
    work(map_grid1d_m(ig))  = DCONJG(vtemp(ig,1))
    work(map_grid1d_p(ig))  = vtemp(ig,1)
  ENDDO
  if(g0_stat)work(map_grid1d_p(1))  = vtemp(1,1)
  CALL  fft_backward(work) ! fft routine is always one dimension

  do ig=1,nnr1
    v_pot(ig,1)=DREAL(work(ig))
    if(lopen_shell.and.nlsd.eq.2)v_pot(ig,2)=v_pot(ig,1)
  enddo
  DEALLOCATE(WORK)
  END SUBROUTINE
  END MODULE
