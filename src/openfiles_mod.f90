!> Module to handle input output files.
!! 
module openfiles_mod

  implicit none

  private

  public :: get_file_unit, open_file, open_file_to_read

contains

  !> Returns a unit number that is not in use.
  !! \param io_max Maximum units to search.
  !! \param get_file_unit Unit return to use for the file.
  !!
  integer function get_file_unit(io_max)
    implicit none
    integer :: io_max, io, m, iostat
    logical :: opened
    
    m = io_max ; if (m < 1) m = 97
    do io = m,1,-1
      inquire (unit=io, opened=opened, iostat=iostat)
      if(iostat.ne.0) cycle
      if(.not.opened) exit
    end do
    get_file_unit = io
 
  end function get_file_unit

  !> Opens a file to write.
  !! \param io Unit for the file.
  !! \param name Name of the file.
  !!
  subroutine open_file(io,name)
    implicit none
    character(len=*) :: name
    character(100) :: io_name
    integer :: io

    io=get_file_unit(100)
    io_name=trim(name)
    open(io,file=io_name)

  end subroutine open_file

  !> Opens a file to read.
  !! \param io Unit for the file.
  !! \param name Name of the file.
  !!
  subroutine open_file_to_read(io,name)
    implicit none
    character(len=*) :: name
    character(100) :: io_name
    integer :: io
    logical :: exists
    
    io=get_file_unit(100)
    io_name=trim(name)
     
    inquire(file=io_name, exist=exists) 
    if(.not.exists)then 
      write(*,*)"File ",io_name,"does not exist ..."
      stop
    endif 
    open(io,file=io_name,status="old")

  end subroutine open_file_to_read

end module openfiles_mod
