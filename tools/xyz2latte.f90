! Small program to convert coordinates .xyz file to
! inputblock.dat which is the coordinate format used by latte
! Written by C. F. A. Negre. Nov. 2014. Los Alamos Nat. Lab.


program xyz2latte

  implicit none
  real(8) :: xmin,xmax,ymin,ymax,zmin,zmax
  integer :: N, i, io
  character :: dummy
  character(2), allocatable :: atom(:)
  real(8), allocatable :: x(:),y(:),z(:)

  write(*,*) 'Reading coordinates ...'

  open(1,file='coords.xyz') ! Default coordinates name.

  read(1,*) N ! Number of atoms.

  allocate(atom(N),x(N),y(N),z(N)) ! Allocate manes and coordinates.

  read(1,'(A1)',IOSTAT=io) dummy  ! Reads the dummy name of the xyz file.

  write(*,*)'Number of atoms', N

  do i=1,N ! Reading names and coordinates.
    read(1,*)atom(i),x(i),y(i),z(i)
    write(*,*)atom(i),x(i),y(i),z(i)
  enddo

  xmin=1000000 ! Initial guess for the boundaries.
  xmax=-100000
  ymin=100000
  ymax=-1000000
  zmin=1000000
  zmax=-10000000

  do i=1,N  ! Searching for the boundaries.

    if(xmin.GT.x(i))xmin=x(i)
    if(xmax.LT.x(i))xmax=x(i)

    if(ymin.GT.y(i))ymin=y(i)
    if(ymax.LT.y(i))ymax=y(i)

    if(zmin.GT.z(i))zmin=z(i)
    if(zmax.LT.z(i))zmax=z(i)

  enddo

  write(*,*)'xmin xmax', xmin,xmax
  write(*,*)'ymin ymax', ymin,ymax
  write(*,*)'zmin zmax', zmin,zmax

  xmin=xmin-5.0 ! Adding some space to the simulation box.
  xmax=xmax+5.0
  ymin=ymin-5.0
  ymax=ymax+5.0
  zmin=zmin-5.0
  zmax=zmax+5.0

  open(2,file='inputblock.dat')

  ! inputblock.dat format. See LATTE documentation files.
  !write(2,*)'NATS=',N
  write(2,*)N
  ! write(2,*)'1.0'
  write(2,"(3F20.5)")xmax-xmin,0.0,0.0
  write(2,"(3F20.5)")0.0,ymax-ymin,0.0
  write(2,"(3F20.5)")0.0,0.0,zmax-zmin
  do i=1,N
    write(2,'(A2,3F20.5)')atom(i),x(i),y(i),z(i)
  enddo

end
