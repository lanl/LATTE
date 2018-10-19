!> Module to handle input output files.
!!
MODULE OPENFILES_MOD

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: GET_FILE_UNIT, OPEN_FILE, OPEN_FILE_TO_READ

CONTAINS

  !> Returns a unit number that is not in use.
  !! \param io_max Maximum units to search.
  !! \param get_file_unit Unit return to use for the file.
  !!
  INTEGER FUNCTION GET_FILE_UNIT(IO_MAX)
    IMPLICIT NONE
    INTEGER :: IO_MAX, IO, M, IOSTAT
    LOGICAL :: OPENED

    M = IO_MAX ; IF (M < 1) M = 97
    DO IO = M,1,-1
       INQUIRE (UNIT=IO, OPENED=OPENED, IOSTAT=IOSTAT)
       IF(IOSTAT.NE.0) CYCLE
       IF(.NOT.OPENED) EXIT
    END DO
    GET_FILE_UNIT = IO

  END FUNCTION GET_FILE_UNIT

  !> Opens a file to write.
  !! \param io Unit for the file.
  !! \param name Name of the file.
  !!
  SUBROUTINE OPEN_FILE(IO,NAME)
    IMPLICIT NONE
    CHARACTER(LEN=*) :: NAME
    CHARACTER(100) :: IO_NAME
    INTEGER :: IO

    IO=GET_FILE_UNIT(100)
    IO_NAME=TRIM(NAME)
    OPEN(IO,FILE=IO_NAME)

  END SUBROUTINE OPEN_FILE

  !> Opens a file to read.
  !! \param io Unit for the file.
  !! \param name Name of the file.
  !!
  SUBROUTINE OPEN_FILE_TO_READ(IO,NAME)
    IMPLICIT NONE
    CHARACTER(LEN=*) :: NAME
    CHARACTER(100) :: IO_NAME
    INTEGER :: IO
    LOGICAL :: EXISTS

    IO=GET_FILE_UNIT(100)
    IO_NAME=TRIM(NAME)

    INQUIRE(FILE=IO_NAME, EXIST=EXISTS)
    IF(.NOT.EXISTS)THEN
       WRITE(*,*)"FILE ",IO_NAME,"DOES NOT EXIST ..."
       STOP
    ENDIF
    OPEN(IO,FILE=IO_NAME,STATUS="OLD")

  END SUBROUTINE OPEN_FILE_TO_READ

END MODULE OPENFILES_MOD
