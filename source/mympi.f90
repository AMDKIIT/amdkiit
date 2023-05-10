MODULE mympi

CONTAINS

  SUBROUTINE MPI_Start()
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  INTEGER :: i_err
!#if defined (_MPI)
  CALL mpi_init(i_err)
!#endif
  END SUBROUTINE MPI_Start




  SUBROUTINE MPI_Stop()
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  INTEGER :: i_err
!#if defined (_MPI)
  call MPI_FINALIZE(i_err)
!#endif
  END SUBROUTINE MPI_Stop




  SUBROUTINE MPI_get_ncpu(ncpu)
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  INTEGER :: ncpu, i_err
  ncpu=1
!#if defined (_MPI)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,i_err)
!#endif
  END SUBROUTINE MPI_get_ncpu




  SUBROUTINE MPI_get_cpuid(icpu)
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  INTEGER :: icpu, i_err
  icpu=0
!#if defined (_MPI)
  CALL  MPI_COMM_RANK(MPI_COMM_WORLD,icpu,i_err)
!#endif
  icpu=icpu+1 !Number starts with 1
  END SUBROUTINE MPI_get_cpuid




  SUBROUTINE MPI_IBcast(myint,leng)
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  INTEGER :: leng, myint(*), i_err
!#if defined (_MPI)
  CALL MPI_BCAST(myint,leng,MPI_INTEGER,0,MPI_COMM_WORLD,i_err)
!#endif
  END SUBROUTINE MPI_IBcast




  SUBROUTINE MPI_RBcasts(myreal,srce)
  use kinds
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  REAL(kind=dp) :: myreal
  INTEGER :: leng, i_err,srce
  leng=1
!#if defined (_MPI)
  CALL MPI_BCAST(myreal,leng,MPI_DOUBLE_PRECISION,srce-1,MPI_COMM_WORLD,i_err)
!#endif
  END SUBROUTINE MPI_RBcasts




  SUBROUTINE MPI_RBcast(myreal,leng)
  use kinds
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  REAL(kind=dp) :: myreal(*)
  INTEGER :: leng, i_err
!#if defined (_MPI)
  CALL MPI_BCAST(myreal,leng,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,i_err)
!#endif
  END SUBROUTINE MPI_RBcast




  SUBROUTINE MPI_Sync_Procs
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  INTEGER i_err
!#if defined (_MPI)
  CALL MPI_Barrier(MPI_COMM_WORLD,i_err)
!#endif
  END SUBROUTINE MPI_Sync_Procs




  SUBROUTINE Set_IOnode(ionode)
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  LOGICAL :: ionode
  INTEGER :: icpu, i_err
  ionode=.false.
  icpu=0
!#if defined (_MPI)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD,icpu,i_err)
!#endif
  IF(icpu.EQ.0)ionode=.TRUE.
  END SUBROUTINE Set_IOnode



  
  SUBROUTINE MPI_GlobMaxR2s(myreal_in)
  USE kinds
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  REAL(kind=dp) :: myreal_in
  REAL(kind=dp) :: myreal_out
  INTEGER :: leng,i_err
  leng=1
!#if defined (_MPI)
  CALL MPI_Allreduce(myreal_in,myreal_out,leng,MPI_DOUBLE_PRECISION, &
                   MPI_MAX,MPI_COMM_WORLD,i_err)
!#endif
   myreal_in=myreal_out
  END SUBROUTINE MPI_GlobMaxR2s
  



  SUBROUTINE MPI_GlobSumR2(myreal_in,leng)
  USE kinds
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  REAL(kind=dp) :: myreal_in(*)
  INTEGER :: leng,i_err
  REAL(kind=dp), ALLOCATABLE :: myreal_out(:)
  ALLOCATE(myreal_out(leng))
!#if defined (_MPI)
  CALL MPI_Allreduce(myreal_in,myreal_out,leng,MPI_DOUBLE_PRECISION, &
                   MPI_SUM,MPI_COMM_WORLD,i_err)
!#endif
  myreal_in(1:leng)=myreal_out(1:leng)
  DEALLOCATE(myreal_out)
  END SUBROUTINE MPI_GlobSumR2
  



  SUBROUTINE MPI_GlobSumR2s(myreal_in)
  USE kinds
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  REAL(kind=dp) :: myreal_in
  REAL(kind=dp) :: myreal_out
  INTEGER :: leng,i_err
  leng=1
!#if defined (_MPI)
  CALL MPI_Allreduce(myreal_in,myreal_out,leng,MPI_DOUBLE_PRECISION, &
                   MPI_SUM,MPI_COMM_WORLD,i_err)
!#endif
   myreal_in=myreal_out
  END SUBROUTINE MPI_GlobSumR2s




  SUBROUTINE MPI_GlobSumR(myreal_in,myreal_out,leng)
  USE kinds
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  REAL(kind=dp) :: myreal_in(*), myreal_out(*)
  INTEGER :: leng,i_err
!#if defined (_MPI)
  CALL MPI_Allreduce(myreal_in,myreal_out,leng,MPI_DOUBLE_PRECISION, &
                   MPI_SUM,MPI_COMM_WORLD,i_err)
!#endif
  END SUBROUTINE MPI_GlobSumR




  SUBROUTINE MPI_GlobSumI(myint_in,myint_out,leng)
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  INTEGER :: myint_in(*), myint_out(*)
  INTEGER :: leng,i_err
!#if defined (_MPI)
  CALL MPI_Allreduce(myint_in,myint_out,leng,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,i_err)
!#endif
  END SUBROUTINE MPI_GlobSumI




  SUBROUTINE MPI_GlobSumI2(myint_in,leng)
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  INTEGER :: myint_in(*)
  INTEGER :: leng,i_err
  INTEGER, ALLOCATABLE::  myint_out(:)
!#if defined (_MPI)
  ALLOCATE(myint_out(leng))
  CALL MPI_Allreduce(myint_in,myint_out,leng,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,i_err)
  myint_in(1:leng)=myint_out(1:leng)
  DEALLOCATE(myint_out)
!#endif
  END SUBROUTINE MPI_GlobSumI2





  SUBROUTINE MPI_GlobSumI2s(myint_in)
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  INTEGER :: myint_in
  INTEGER :: leng,i_err
  INTEGER ::  myint_out
!#if defined (_MPI)
  leng=1
  CALL MPI_Allreduce(myint_in,myint_out,leng,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,i_err)
  myint_in=myint_out
!#endif
  END SUBROUTINE MPI_GlobSumI2s





  SUBROUTINE MPI_GlobSumC2(mycmplx_in,leng)
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  complex*16 :: mycmplx_in(*)
  INTEGER :: leng,i_err
  COMPLEX, ALLOCATABLE::  mycmplx_out(:)
!#if defined (_MPI)
  ALLOCATE(mycmplx_out(leng))
  CALL MPI_Allreduce(mycmplx_in,mycmplx_out,leng,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,i_err)
  mycmplx_in(1:leng)=mycmplx_out(1:leng)
  DEALLOCATE(mycmplx_out)
!#endif
  END SUBROUTINE MPI_GlobSumC2




  SUBROUTINE MPI_GlobSumC(mycmplx_in,mycmplx_out,leng)
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  COMPLEX*16 :: mycmplx_in(*), mycmplx_out(*)
  INTEGER :: leng,i_err
!#if defined (_MPI)
  CALL MPI_Allreduce(mycmplx_in,mycmplx_out,leng,MPI_DOUBLE_COMPLEX, &
                   MPI_SUM,MPI_COMM_WORLD,i_err)
!#endif
  END SUBROUTINE MPI_GlobSumC




  SUBROUTINE MPI_GlobSumC2s(mycmplx_in)
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  COMPLEX*16 :: mycmplx_in, mycmplx_out
  INTEGER :: leng,i_err
!#if defined (_MPI)
  leng=1
  CALL MPI_Allreduce(mycmplx_in,mycmplx_out,leng,MPI_DOUBLE_COMPLEX, &
                   MPI_SUM,MPI_COMM_WORLD,i_err)
  mycmplx_in=mycmplx_out
!#endif
  END SUBROUTINE MPI_GlobSumC2s




  SUBROUTINE MPI_SHIFT(MSEND,MRECV,MSGLEN,MEP,IP)
  use kinds
  IMPLICIT NONE
  INTEGER MSGLEN,MEP,IP
  REAL(kind=dp)  MSEND(*),MRECV(*)
  INCLUDE 'mpif.h'
  INTEGER status(MPI_STATUS_SIZE),IPSEND,IPRECV,&
          whoami,howmany, &
          ITYPE,IREQUEST,IERR,GID
  gid=MPI_COMM_WORLD

  call mpi_comm_rank(gid,whoami,ierr)
  call mpi_comm_size(gid,howmany,ierr)

  ipsend=mod(whoami+ip+howmany,howmany)
  iprecv=mod(whoami-ip+howmany,howmany)

  itype=1
  irequest=1
  CALL MPI_ISEND(msend,msglen,mpi_byte,ipsend,itype, &
               gid,irequest,ierr)
  CALL MPI_RECV(mrecv,msglen,mpi_byte,iprecv,itype, &
               gid,status,ierr)
  CALL MPI_WAIT(irequest,status,ierr)
  END SUBROUTINE MPI_Shift




  SUBROUTINE MPI_Concat(myint_in,myint_out,leng)
  IMPLICIT NONE
!#if defined (_MPI)
  INCLUDE 'mpif.h'
!#endif
  INTEGER :: myint_in(*)
  INTEGER :: leng,i_err
  INTEGER ::  myint_out(*)
!#if defined (_MPI)
  CALL MPI_AllGather(myint_in,leng,MPI_INTEGER,myint_out,leng,MPI_INTEGER,MPI_COMM_WORLD,i_err)
!#endif
  END SUBROUTINE MPI_Concat




!TODO: currently, taken from CPMD. Need to re-write
  SUBROUTINE MPI_trans(OUTMSG,INMSG,LENG)
  IMPLICIT NONE
!ifdef PARALLEL
  INCLUDE 'mpif.h'
!endif
! Arguments
  COMPLEX*16 OUTMSG(*),INMSG(*)
  INTEGER :: LENG
! Variables
  INTEGER IERR,GID

  GID=MPI_COMM_WORLD
  CALL MPI_ALLTOALL(OUTMSG,LENG,MPI_DOUBLE_COMPLEX,&
                        INMSG,LENG,MPI_DOUBLE_COMPLEX,GID,IERR)
  IF(IERR.NE.0) STOP 'MY_TRANS | IERR.NE.0'
  END SUBROUTINE MPI_trans

END MODULE mympi

