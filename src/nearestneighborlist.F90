subroutine nearestneighborlist(nrnnlist,nndist,nnRx,nnRy,nnRz,nnType,nnStruct,nrnnStruct, &
  R_X,R_Y,R_Z,LBox,Rcut,Nr_atoms,Max_Nr_Neigh)

use omp_lib
implicit none

integer,    parameter       :: PREC = 8
real(PREC), parameter       :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0
integer,    intent(in)      :: Nr_atoms, Max_Nr_Neigh
real(PREC), intent(in)      :: Rcut 
real(PREC), intent(in)      :: R_X(Nr_atoms), R_Y(Nr_atoms), R_Z(Nr_atoms), LBox(3)
real(PREC)                  :: Rx(10*Nr_atoms), Ry(10*Nr_atoms), Rz(10*Nr_atoms)
integer,    intent(out)     :: nrnnlist(Nr_atoms), nnType(Nr_atoms,Max_Nr_Neigh)
integer,    intent(out)     :: nnStruct(Nr_atoms,Nr_atoms), nrnnStruct(Nr_atoms)
real(PREC), intent(out)     :: nndist(Nr_atoms,Max_Nr_Neigh), nnRx(Nr_atoms,Max_Nr_Neigh)
real(PREC), intent(out)     :: nnRy(Nr_atoms,Max_Nr_Neigh), nnRz(Nr_atoms,Max_Nr_Neigh)

integer       :: i,j,k,l,m,t,nx,ny,nz,ss,cell,type(10*Nr_atoms),Nskin,N, NNtmp
integer       :: head(Nr_atoms*Nr_atoms), list(10*Nr_atoms), tmp(Nr_atoms), cnt, cnt2, cct
real(PREC)    :: Lx,Ly,Lz,Tx,Ty,Tz,RR(3),TT(3),dist,dLx,dLy,dLz

N = Nr_atoms
Lx = LBox(1) 
Ly = LBox(2) 
Lz = LBox(3) ! Dimensions of periodic BC
nx = floor(Lx/Rcut) 
ny = floor(Ly/Rcut) 
nz = floor(Lz/Rcut) ! Division into # cell boxes: nx, ny, nz
Rx(1:Nr_atoms) = R_X(1:Nr_atoms) !+Lx
Ry(1:Nr_atoms) = R_Y(1:Nr_atoms) !+Ly
Rz(1:Nr_atoms) = R_Z(1:Nr_atoms) !+Lz

nnStruct = ZERO
nrnnStruct = ZERO
if ((min(nx,ny,nz).lt.3).or.(Nr_atoms.lt.80)) then
!if (min(nx,ny,nz).lt.3) then

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,l,m,ss,RR,TT,Tx,Ty,Tz,dist,tmp,cnt,cnt2) 
 do i = 1,N
   cnt = 0
   tmp = ZERO
   RR(1) = Rx(i)
   RR(2) = Ry(i)
   RR(3) = Rz(i)
   do m = 1,N
     do j = -1,1
     do k = -1,1
     do l = -1,1
       Tx = Rx(m)+j*Lx  ! Search all neigbors within a single translation (multiple translations could be necessary for small systems!
       Ty = Ry(m)+k*Ly
       Tz = Rz(m)+l*Lz
       TT(1) = Tx
       TT(2) = Ty
       TT(3) = Tz
       dist = norm2(RR-TT)
       !if ((dist.lt.Rcut).and.(dist.gt.1e-12)) then ! Neighbors within Rcut inlcuidng translated atoms in the "skin"
       if ((dist < Rcut)) then
         cnt = cnt + 1
         nndist(i,cnt) = dist
         nnRx(i,cnt) = Tx
         nnRy(i,cnt) = Ty
         nnRz(i,cnt) = Tz
         nnType(i,cnt) = m  ! Neigbor is number of original ordering number m in the box that might have been stranslated to the skin
         tmp(m) = m
       endif
     enddo
     enddo
     enddo
   enddo
   nrnnlist(i) = cnt
   cnt2 = 0
   do ss = 1,N
     if (tmp(ss).gt.0) then  ! Includes only neighbors in the box within Rcut (without the skin)
       cnt2 = cnt2 + 1
       nnStruct(i,cnt2) = ss
     endif
   enddo
   nrnnStruct(i) = cnt2
 enddo
!$OMP END PARALLEL DO

else

 head = ZERO ! Linked list that keeps track of all atoms in the nx*ny*nz small boxes
 list = ZERO !
 do i = 1,N
  cell = 1 + floor(nx*Rx(i)/Lx) + floor(ny*Ry(i)/Ly)*nx + floor(nz*Rz(i)/Lz)*nx*ny
  list(i) = head(cell)
  head(cell) = i
  type(i) = i
 enddo

 !%%% And now add a skin or surface buffer to account for periodic BC, all 26 of them!
 cnt = 0
 do i = 1,nx*ny  ! All boxes in the first (z=0) layer
   t = head(i)
   do while (t.gt.0)    ! and all atoms of this first layer
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)
      Ry(N+cnt) = Ry(t)
      Rz(N+cnt) = Rz(t)+Lz ! Add skin atoms with coordinates translated by Lz
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = 1,nx*ny*nz,nx  ! All boxes in another (x=0) layer
   t = head(i)
   do while (t.gt.0)   ! and all atoms in that layer
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)+Lx ! Add skin atoms with coordinates translated by Lx
      Ry(N+cnt) = Ry(t)
      Rz(N+cnt) = Rz(t)
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = 1,nx*ny*nz,nx*ny  ! Continue ...
   do k = i,i+nx-1
     t = head(k)
     do while (t.gt.0) 
        cnt = cnt + 1
        Rx(N+cnt) = Rx(t)
        Ry(N+cnt) = Ry(t)+Ly
        Rz(N+cnt) = Rz(t)
        type(N+cnt) = t
        t = list(t)
     enddo
   enddo
 enddo
 cct = 0
 do i = nx*ny*(nz-1)+1,nx*ny*nz
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      cct = cct + 1
      Rx(N+cnt) = Rx(t)
      Ry(N+cnt) = Ry(t)
      Rz(N+cnt) = Rz(t)-Lz
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = nx,nx*ny*nz,nx
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)-Lx
      Ry(N+cnt) = Ry(t)
      Rz(N+cnt) = Rz(t)
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = nx*(ny-1)+1,nx*ny*nz,nx*ny
   do k = i,i+nx-1
     t = head(k)
     do while (t.gt.0)
        cnt = cnt + 1
        Rx(N+cnt) = Rx(t)
        Ry(N+cnt) = Ry(t)-Ly
        Rz(N+cnt) = Rz(t)
        type(N+cnt) = t
        t = list(t)
     enddo
   enddo
 enddo
 do i = 1,nx
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)
      Ry(N+cnt) = Ry(t)+Ly
      Rz(N+cnt) = Rz(t)+Lz
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = nx*(ny-1)+1,nx*(ny-1)+nx
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)
      Ry(N+cnt) = Ry(t)-Ly
      Rz(N+cnt) = Rz(t)+Lz
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = nx*ny*(nz-1)+1,nx*ny*(nz-1)+nx
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)
      Ry(N+cnt) = Ry(t)+Ly
      Rz(N+cnt) = Rz(t)-Lz
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = nx*ny*(nz-1)+nx*(ny-1)+1,nx*ny*(nz-1)+nx*(ny-1)+nx
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)
      Ry(N+cnt) = Ry(t)-Ly
      Rz(N+cnt) = Rz(t)-Lz
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = 1,nx*ny,nx
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)+Lx
      Ry(N+cnt) = Ry(t)
      Rz(N+cnt) = Rz(t)+Lz
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = nx*ny*(nz-1)+1,nx*ny*nz,nx
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)+Lx
      Ry(N+cnt) = Ry(t)
      Rz(N+cnt) = Rz(t)-Lz
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = nx*ny*(nz-1)+nx,nx*ny*nz,nx
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)-Lx
      Ry(N+cnt) = Ry(t)
      Rz(N+cnt) = Rz(t)-Lz
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = nx,nx*ny,nx
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)-Lx
      Ry(N+cnt) = Ry(t)
      Rz(N+cnt) = Rz(t)+Lz
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = 1,nx*ny*nz,nx*ny
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)+Lx
      Ry(N+cnt) = Ry(t)+Ly
      Rz(N+cnt) = Rz(t)
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = nx*(ny-1)+1,nx*ny*nz,nx*ny
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)+Lx
      Ry(N+cnt) = Ry(t)-Ly
      Rz(N+cnt) = Rz(t)
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = nx*(ny-1)+nx,nx*ny*nz,nx*ny
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)-Lx
      Ry(N+cnt) = Ry(t)-Ly
      Rz(N+cnt) = Rz(t)
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = nx,nx*ny*nz,nx*ny
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)-Lx
      Ry(N+cnt) = Ry(t)+Ly
      Rz(N+cnt) = Rz(t)
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 t = head(1)
 do while (t.gt.0)
    cnt = cnt + 1
    Rx(N+cnt) = Rx(t)+Lx
    Ry(N+cnt) = Ry(t)+Ly
    Rz(N+cnt) = Rz(t)+Lz
    type(N+cnt) = t
    t = list(t)
 enddo
 t = head(nx)
 do while (t.gt.0)
    cnt = cnt + 1
    Rx(N+cnt) = Rx(t)-Lx
    Ry(N+cnt) = Ry(t)+Ly
    Rz(N+cnt) = Rz(t)+Lz
    type(N+cnt) = t
    t = list(t)
 enddo
 t = head(nx*(ny-1)+1)
 do while (t.gt.0)
    cnt = cnt + 1
    Rx(N+cnt) = Rx(t)+Lx
    Ry(N+cnt) = Ry(t)-Ly
    Rz(N+cnt) = Rz(t)+Lz
    type(N+cnt) = t
    t = list(t)
 enddo
 t = head(nx*(ny-1)+nx)
 do while (t.gt.0)
    cnt = cnt + 1
    Rx(N+cnt) = Rx(t)-Lx
    Ry(N+cnt) = Ry(t)-Ly
    Rz(N+cnt) = Rz(t)+Lz
    type(N+cnt) = t
    t = list(t)
 enddo
 t = head(nx*ny*(nz-1)+1)
 do while (t.gt.0)
    cnt = cnt + 1
    Rx(N+cnt) = Rx(t)+Lx
    Ry(N+cnt) = Ry(t)+Ly
    Rz(N+cnt) = Rz(t)-Lz
    type(N+cnt) = t
    t = list(t)
 enddo
 t = head(nx*ny*(nz-1)+nx)
 do while (t.gt.0)
    cnt = cnt + 1
    Rx(N+cnt) = Rx(t)-Lx
    Ry(N+cnt) = Ry(t)+Ly
    Rz(N+cnt) = Rz(t)-Lz
    type(N+cnt) = t
    t = list(t)
 enddo
 t = head(nx*ny*(nz-1)+nx*(ny-1)+1)
 do while (t.gt.0)
    cnt = cnt + 1
    Rx(N+cnt) = Rx(t)+Lx
    Ry(N+cnt) = Ry(t)-Ly
    Rz(N+cnt) = Rz(t)-Lz
    type(N+cnt) = t
    t = list(t)
 enddo
 t = head(nx*ny*nz)
 do while (t.gt.0)
    cnt = cnt + 1
    Rx(N+cnt) = Rx(t)-Lx
    Ry(N+cnt) = Ry(t)-Ly
    Rz(N+cnt) = Rz(t)-Lz
    type(N+cnt) = t
    t = list(t)
 enddo
 Nskin = cnt

 !! And now create a list for everything including the Skin/buffer layer
! [max(Rz),min(Rz)]
 dLx = Lx/nx 
 dLy = Ly/ny 
 dLz = Lz/nz
 Rx = Rx+dLx 
 Ry = Ry+dLy 
 Rz = Rz+dLz ! Shift to avoid negative coordinates
 Lx = Lx + 2*dLx 
 Ly = Ly + 2*dLy 
 Lz = Lz + 2*dLz
 nx = nx+2
 ny = ny+2 
 nz = nz+2
 head = ZERO
 list = ZERO
 do i = 1,N+Nskin
  cell = 1 + floor(nx*Rx(i)/Lx) + floor(ny*Ry(i)/Ly)*nx + floor(nz*Rz(i)/Lz)*nx*ny
  list(i) = head(cell)
  head(cell) = i
 enddo

 do i = 1,N ! Go through all atoms
   cnt = 0
   cnt2 = 0
   do j = -1,1  ! Translate position to neighboring small boxes
   do k = -1,1
   do l = -1,1
     Tx = Rx(i)+j*dLx
     Ty = Ry(i)+k*dLy
     Tz = Rz(i)+l*dLz
     cell = 1 + floor(nx*Tx/Lx) + floor(ny*Ty/Ly)*nx + floor(nz*Tz/Lz)*nx*ny ! and extract all atoms in those small
     t = head(cell)                                                          ! neighbor boxes ...
     do while (t.gt.0)
       RR(1) = Rx(i) - Rx(t)
       RR(2) = Ry(i) - Ry(t)
       RR(3) = Rz(i) - Rz(t)
       dist = norm2(RR)
!       if (i == 70) then
!         write(*,*) ' dist = ', i,dist
!       endif
       if (dist.lt.Rcut) then  ! All atoms including the skin within Rcut
         cnt = cnt + 1
         nndist(i,cnt) = dist
         nnRx(i,cnt) = Rx(t)-dLx  ! Coordinates without the shift
         nnRy(i,cnt) = Ry(t)-dLy
         nnRz(i,cnt) = Rz(t)-dLz
         nnType(i,cnt) = type(t)
         if (t.le.N) then ! Keeps track of atoms within the original box without the skin
           cnt2 = cnt2 + 1
           nnStruct(i,cnt2) = type(t)
         endif
       endif
       t = list(t)
     enddo
   enddo
   enddo
   enddo
   nrnnlist(i) = cnt
   nrnnStruct(i) = cnt2
!   nrnnStruct(i) = cnt
 enddo
endif

end subroutine nearestneighborlist
