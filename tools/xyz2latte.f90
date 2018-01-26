! Small program to convert coordinates .xyz file to
! inputblock.dat which is the coordinate format used by latte
! Written by C. F. A. Negre. Nov. 2014. Los Alamos Nat. Lab.


PROGRAM xyz2latte

  IMPLICIT NONE
  REAL(8) :: xmin,xmax,ymin,ymax,zmin,zmax
  INTEGER :: N, i, io
  CHARACTER :: dummy
  CHARACTER(2), ALLOCATABLE :: atom(:)
  REAL(8), ALLOCATABLE :: x(:),y(:),z(:)

  WRITE(*,*) 'Reading coordinates ...'

  OPEN(1,file='coords.xyz') ! Default coordinates name.

  READ(1,*) N ! Number of atoms.

  ALLOCATE(atom(N),x(N),y(N),z(N)) ! Allocate manes and coordinates.

  READ(1,'(A1)',IOSTAT=io) dummy  ! Reads the dummy name of the xyz file.

  WRITE(*,*)'Number of atoms', N

  DO i=1,N ! Reading names and coordinates.
     READ(1,*)atom(i),x(i),y(i),z(i)
     WRITE(*,*)atom(i),x(i),y(i),z(i)
  ENDDO

  xmin=1000000 ! Initial guess for the boundaries.
  xmax=-100000
  ymin=100000
  ymax=-1000000
  zmin=1000000
  zmax=-10000000

  DO i=1,N  ! Searching for the boundaries.

     IF(xmin.GT.x(i))xmin=x(i)
     IF(xmax.LT.x(i))xmax=x(i)

     IF(ymin.GT.y(i))ymin=y(i)
     IF(ymax.LT.y(i))ymax=y(i)

     IF(zmin.GT.z(i))zmin=z(i)
     IF(zmax.LT.z(i))zmax=z(i)

  ENDDO

  WRITE(*,*)'xmin xmax', xmin,xmax
  WRITE(*,*)'ymin ymax', ymin,ymax
  WRITE(*,*)'zmin zmax', zmin,zmax

  xmin=xmin-5.0 ! Adding some space to the simulation box.
  xmax=xmax+5.0
  ymin=ymin-5.0
  ymax=ymax+5.0
  zmin=zmin-5.0
  zmax=zmax+5.0

  OPEN(2,file='inputblock.dat')

  ! inputblock.dat format. See LATTE documentation files.
  WRITE(2,*)N
  WRITE(2,"(3F20.5)")xmax-xmin,0.0,0.0
  WRITE(2,"(3F20.5)")0.0,ymax-ymin,0.0
  WRITE(2,"(3F20.5)")0.0,0.0,zmax-zmin
  DO i=1,N
     WRITE(2,'(A2,3F20.5)')atom(i),x(i),y(i),z(i)
  ENDDO

END PROGRAM xyz2latte
