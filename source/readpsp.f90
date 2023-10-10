MODULE READpsp
  CONTAINS

  SUBROUTINE psp
  USE system_data_types
  IMPLICIT NONE
  ALLOCATE(slatr_ex(sp_t),tnum(sp_t))
  END SUBROUTINE psp


  SUBROUTINE read_psp(pspfile,isp)
  USE system_data_types
  USE readstring

  IMPLICIT NONE
  INTEGER :: ISP
  INTEGER   IUNIT,ios
  PARAMETER (IUNIT=20)
  CHARACTER, INTENT(in) :: pspfile*(*)
  INTEGER I,L

  CHARACTER LINE*80, PATHOFPP*42
  TNUM(ISP)=.FALSE.
     OPEN(UNIT=IUNIT,FILE=trim(PATHOFINPUT)//pspfile,IOSTAT=ios)
     !OPEN(UNIT=IUNIT,FILE=trim(PATHOFINPUT)//'PSEUDOPOTENTIAL/'//pspfile,IOSTAT=ios)
       IF(ios.NE.0)THEN
         IF(ionode)WRITE(*,*)"    Reading ",pspfile," is not successful"
         STOP
       ENDIF
       !ATOM section
       CALL search(iunit,"ATOM")
       CALL str2var(iunit,"Z ",line)
       READ(line,*)z(isp)
       CALL str2var(iunit,"ZV",line)
       READ(line,*)zv(isp)
       CALL str2var(iunit,"XC",line)
       CALL str2var(iunit,"TYPE",line)
       IF(INDEX(line,'NUMERIC').NE.0)tnum(isp)=.TRUE.
       CALL search(IUNIT,"END")
       !INFO section
       CALL search(IUNIT,"INFO")
       11 CONTINUE
       READ(iunit,END=20,ERR=20,FMT='(A)') line
       IF(INDEX(line,'&END').EQ.0) GOTO 11
       ! POTENTIAL section
       call search(iunit,"POTENTIAL")
       READ(iunit,END=20,ERR=20,FMT='(A)') line
       READ(line,*)meshv(isp)
         DO i=1,meshv(isp)
           READ(iunit,*) rr(i,isp),(vr(i,isp,l),l=1,lmax(isp))
         ENDDO
       call search(IUNIT,"END")

       ! WAVEFUNCTION section
       CALL search(iunit,"WAVEFUNCTION")
       READ(iunit,END=20,ERR=20,FMT='(A)') line
       READ(line,*)meshw(isp)
         DO i=1,meshw(isp)
           READ(iunit,*) rr(i,isp),(rps(i,isp,l),l=1,lmax(isp))
         ENDDO
       CALL search(IUNIT,"END")
      IF(meshw(isp).NE.meshv(isp))WRITE(6,*)"ERROR:INCOMPATIBLE MESHES....",meshw(isp),meshv(isp)
    RETURN
    
    20 WRITE(*,*) "Error occurred while reading file"
    STOP
    END SUBROUTINE read_psp



    SUBROUTINE search (id, string)
    implicit none
    integer :: id
    character (len=*) :: string  ! Label to be matched
    character (len=80) :: read_string ! String read from file
    integer :: ios

    ios = 0
    do while (ios.eq.0)
    read (id,iostat = ios, err = 99,FMT='(A)') read_string
    if (INDEX(read_string,"&"//string).NE.0 ) then
    return
    endif
    enddo
    99 print*, 'ERROR!!! No "',string,'" block found in psedopotential file, IOSTAT=', abs (ios)
    STOP
    end subroutine search
    
    SUBROUTINE str2var(id,string,var)
    USE system_data_types
    IMPLICIT NONE
    CHARACTER, INTENT(in) :: string*(*)
    CHARACTER, INTENT(out):: var*(*)
    CHARACTER (len=100) :: read_string
    INTEGER IOS,ID
     READ(id,iostat = ios,err=99,FMT='(A)')read_string
       IF(INDEX(read_string,string).NE.0 ) then
        var =trim(read_string(INDEX(read_string, '=')+1:len_trim(read_string)))
        RETURN
       ENDIF
     99 PRINT*, 'ERROR!!! No "',TRIM(string),'" block found in psedopotential file,IOSTAT=', abs (ios)
    STOP
    END SUBROUTINE str2var

END MODULE
