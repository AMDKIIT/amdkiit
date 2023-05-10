MODULE math
  CONTAINS

  SUBROUTINE inverse_matrix(a,b)
    !
    !Inverse of a matrix "a" is calculated and returned in "b"
    !Algorithm used: Gauss Jordan with pioviting
    
    USE system_data_types
    USE kinds
    IMPLICIT NONE
    REAL(KIND=dp), POINTER, DIMENSION(:,:) :: a,b
 
    REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE:: d
    INTEGER, DIMENSION(:), ALLOCATABLE :: ip
 
    REAL(KIND=dp) :: dmax, s
    INTEGER :: n, i, j, k, l, imax, ip_tmp
 
    n=SIZE(a,1)
 
    ALLOCATE(ip(n))
    ALLOCATE(d(n,n))
 
 
    b=0.0d0
    DO i=1,n
      b(i,i)=1.0D0
    END DO
    DO k=1,n
      ip(k)=k
    END DO
 
    DO k=1,n
      dmax=ABS(a(ip(k),ip(k)))
      imax=k
      DO l=k,n
        IF(ABS(a(ip(l),ip(k))) > dmax)THEN
          dmax=ABS(a(ip(l),ip(k)))
          imax=l
        END IF
      END DO
 
      ip_tmp = ip(imax)
      ip(imax)=ip(k)
      ip(k)=ip_tmp
      DO i=1,n
        IF(i/=k)THEN
!         Find the factor for elimination
          s=a(ip(i),k)/a(ip(k),k)
!         Correct the rest of the entries of "A" along column starting from k
          A(ip(i),k:n)=A(ip(i),k:n)-A(ip(k),k:n)*s
          !print*,A(ip(i),k:n)
!         Apply same correction to "b" matrix
          b(ip(i),1:n)=b(ip(i),1:n)-b(ip(k),1:n)*s
        END IF
      END DO
    END DO

    d=0._dp
    DO i=1,n
      d(i,1:n)=b(ip(i),1:n)/a(ip(i),i)
    END DO
    b(:,:)=d(:,:)
    DEALLOCATE(d,ip)
  END SUBROUTINE inverse_matrix

! Calculate Determinant of a square matrix "a" 
  RECURSIVE FUNCTION determinant(a) RESULT(d)
    USE kinds
    IMPLICIT NONE
    REAL(KIND=dp), POINTER, DIMENSION(:,:), INTENT(IN) :: a
    REAL(KIND=dp) :: d

    INTEGER :: n
    INTEGER :: i 
  
    REAL(KIND=dp), POINTER :: a_red(:,:)
    REAL(KIND=dp) :: minor

    n=SIZE(a,1)
    IF(n>1)THEN
      d=0._dp
      ALLOCATE(a_red(SIZE(a,1)-1,SIZE(a,2)-1))
      DO i=1,n
        a_red(1:i-1,1:n-1)=a(1:i-1,2:n)
        a_red(i:n-1,1:n-1)=a(i+1:n,2:n)
        minor=determinant(a_red)
        d=d+(1._dp)**(i+1)*a(i,1)*minor
      END DO
      DEALLOCATE(a_red)
    ELSE
      d=a(1,1)
    END IF
    
  END FUNCTION determinant
!

  SUBROUTINE sort_array(n,ra,ind)
  USE kinds
  USE constants

  implicit none
  !-input/output variables
  integer, intent(in) :: n
  integer, intent(inout) :: ind (*)
  real(KIND=dp), intent(inout) :: ra (*)
  !-local variables
  integer :: i, ir, j, l, iind
  real(DP) :: rra
  ! initialize index array
  if (ind (1) .eq.0) then
     do i = 1, n
        ind (i) = i
     enddo
  endif
  ! nothing to order
  if (n.lt.2) return
  ! initialize indices for hiring and retirement-promotion phase
  l = n / 2 + 1

  ir = n

  sorting: do

    ! still in hiring phase
    if ( l .gt. 1 ) then
       l    = l - 1
       rra  = ra (l)
       iind = ind (l)
       ! in retirement-promotion phase.
    else
       ! clear a space at the end of the array
       rra  = ra (ir)
       !
       iind = ind (ir)
       ! retire the top of the heap into it
       ra (ir) = ra (1)
       !
       ind (ir) = ind (1)
       ! decrease the size of the corporation
       ir = ir - 1
       ! done with the last promotion
       if ( ir .eq. 1 ) then
          ! the least competent worker at all !
          ra (1)  = rra
          !
          ind (1) = iind
          exit sorting
       endif
    endif
    ! wheter in hiring or promotion phase, we
    i = l
    ! set up to place rra in its proper level
    j = l + l
    !
    do while ( j .le. ir )
       if ( j .lt. ir ) then
          ! compare to better underling
          if ( abs(ra(j)-ra(j+1)).ge.eps8 ) then
             if (ra(j).lt.ra(j+1)) j = j + 1
          else
             ! this means ra(j) == ra(j+1) within tolerance
             if (ind (j) .lt.ind (j + 1) ) j = j + 1
          endif
       endif
       ! demote rra
       if ( abs(rra - ra(j)).ge.eps8 ) then
          if (rra.lt.ra(j)) then
             ra (i) = ra (j)
             ind (i) = ind (j)
             i = j
             j = j + j
          else
             ! set j to terminate do-while loop
             j = ir + 1
          end if
       else
          !this means rra == ra(j) within tolerance
          ! demote rra
          if (iind.lt.ind (j) ) then
             ra (i) = ra (j)
             ind (i) = ind (j)
             i = j
             j = j + j
          else
             ! set j to terminate do-while loop
             j = ir + 1
          endif
       end if
    enddo
    ra(i) = rra
    ind(i) = iind

  end do sorting
  !
  END SUBROUTINE sort_array


  SUBROUTINE sort_array2 ( n, arr, index )
        use kinds
        IMPLICIT NONE

        INTEGER n
        REAL(kind=dp)  arr(1:n)
        INTEGER index(1:n)

        INTEGER m, nstack
        PARAMETER ( m = 7, nstack = 50 )
        INTEGER i, l, ib, ir, jstack, k, istack(1:nstack), itemp, j
        REAL(kind=dp) a, temp

        DO i = 1, n
           index(i) = i
        END DO
        jstack = 0
        l = 1
        ir = n
1       IF (ir-l.lt.m) THEN
           DO j = l + 1, ir
              a = arr(j)
              ib = index(j)
              DO i = j - 1, 1, -1
                 IF (arr(i).le.a) GO TO 2
                 arr(i+1) = arr(i)
                 index(i+1) = index(i)
              END DO
              i = 0
2             arr(i+1) = a
              index(i+1) = ib
           END DO
           IF (jstack.eq.0) RETURN
           ir = istack(jstack)
           l = istack(jstack-1)
           jstack = jstack - 2
        ELSE
           k = (l+ir)/2
           temp = arr(k)
           arr(k) = arr(l+1)
           arr(l+1) = temp
           itemp = index(k)
           index(k) = index(l+1)
           index(l+1) = itemp
           IF (arr(l+1).gt.arr(ir)) THEN
              temp = arr(l+1)
              arr(l+1) = arr(ir)
              arr(ir) = temp
              itemp = index(l+1)
              index(l+1) = index(ir)
              index(ir) = itemp
           END IF
           IF (arr(l).gt.arr(ir)) THEN
              temp = arr(l)
              arr(l) = arr(ir)
              arr(ir) = temp
              itemp = index(l)
              index(l) = index(ir)
              index(ir) = itemp
           END IF
           IF (arr(l+1).gt.arr(l)) THEN
              temp = arr(l+1)
              arr(l+1) = arr(l)
              arr(l) = temp
              itemp = index(l+1)
              index(l+1) = index(l)
              index(l) = itemp
           END IF
           i = l + 1
           j = ir
           a = arr(l)
           ib = index(l)
3          CONTINUE
           i = i + 1
           IF (arr(i).lt.a) GO TO 3
4          CONTINUE
           j = j - 1
           IF (arr(j).gt.a) GO TO 4
           IF (j.lt.i) GO TO 5
           temp = arr(i)
           arr(i) = arr(j)
           arr(j) = temp
           itemp = index(i)
           index(i) = index(j)
           index(j) = itemp
           GO TO 3
5          arr(l) = arr(j)
           arr(j) = a
           index(l) = index(j)
           index(j) = ib
           jstack = jstack + 2
           IF (jstack.gt.nstack) STOP 'SORT2: Nstack too small'
           IF (ir-i+1.ge.j-l) THEN
              istack(jstack) = ir
              istack(jstack-1) = i
              ir = j - 1
           ELSE
              istack(jstack) = j - 1
              istack(jstack-1) = l
              l = i
           END IF
        END IF

        GO TO 1

      END SUBROUTINE sort_array2
      SUBROUTINE ICOPY(N,A,IA,B,IB)
      IMPLICIT NONE
      INTEGER N,IA,IB
      INTEGER A(IA*N),B(IB*N)
      INTEGER I,IAX,IBX
      DO I=1,N
        IAX=(I-1)*IA + 1
        IBX=(I-1)*IB + 1
        B(IBX) = A(IAX)
      ENDDO
      RETURN
      END SUBROUTINE icopy
     
      SUBROUTINE SIMPSN(N,INTE,SUM)
      USE KINDS
      IMPLICIT NONE
      INTEGER   N
      REAL(KIND=DP)    INTE(N),SUM
      REAL(KIND=DP)    C1,C2,C3,C4
      INTEGER   I
      PARAMETER (C1=109.D0/48.D0,C2=-5.D0/48.D0,&
     &     C3=63.D0/48.D0,C4=49.D0/48.D0)
      INTE(1)   = INTE(1)*C1
      INTE(2)   = INTE(2)*C2
      INTE(3)   = INTE(3)*C3
      INTE(4)   = INTE(4)*C4
      INTE(N-1) = INTE(N-1)*C1
      INTE(N-2) = INTE(N-2)*C2
      INTE(N-3) = INTE(N-3)*C3
      INTE(N-4) = INTE(N-4)*C4
      SUM=0.D0

      DO I=1,N-1
         SUM=SUM+INTE(I)
      ENDDO

      RETURN
      END



END MODULE math
