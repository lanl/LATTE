!> Module to handle input output files.
!! 
MODULE openfiles_mod

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: get_file_unit, open_file, open_file_to_read

CONTAINS

  !> Returns a unit number that is not in use.
  !! \param io_max Maximum units to search.
  !! \param get_file_unit Unit return to use for the file.
  !!
  INTEGER FUNCTION get_file_unit(io_max)
    IMPLICIT NONE
    INTEGER :: io_max, io, m, iostat
    LOGICAL :: opened

    m = io_max ; IF (m < 1) m = 97
    DO io = m,1,-1
       INQUIRE (unit=io, opened=opened, iostat=iostat)
       IF(iostat.NE.0) CYCLE
       IF(.NOT.opened) EXIT
    END DO
    get_file_unit = io

  END FUNCTION get_file_unit

  !> Opens a file to write.
  !! \param io Unit for the file.
  !! \param name Name of the file.
  !!
  SUBROUTINE open_file(io,name)
    IMPLICIT NONE
    CHARACTER(len=*) :: name
    CHARACTER(100) :: io_name
    INTEGER :: io

    io=get_file_unit(100)
    io_name=TRIM(name)
    OPEN(io,file=io_name)

  END SUBROUTINE open_file

  !> Opens a file to read.
  !! \param io Unit for the file.
  !! \param name Name of the file.
  !!
  SUBROUTINE open_file_to_read(io,name)
    IMPLICIT NONE
    CHARACTER(len=*) :: name
    CHARACTER(100) :: io_name
    INTEGER :: io
    LOGICAL :: exists

    io=get_file_unit(100)
    io_name=TRIM(name)

    INQUIRE(file=io_name, exist=exists) 
    IF(.NOT.exists)THEN 
       WRITE(*,*)"File ",io_name,"does not exist ..."
       STOP
    ENDIF
    OPEN(io,file=io_name,status="old")

  END SUBROUTINE open_file_to_read

END MODULE openfiles_mod
