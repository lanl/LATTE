subroutine Get_dH(dHx,dHy,dHz,Nr_atoms,dx,HDIM,Max_Nr_Neigh,RX,RY,RZ,H_INDEX_START,H_INDEX_END, &
nrnnlist,nnRx,nnRy,nnRz,nnType,Element_Type,LBox)

use OMP_lib
implicit none
integer, parameter        :: PREC = 8
real(PREC), parameter     :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0, SIX = 6.D0
integer, intent(in)       :: Nr_atoms, HDIM, Max_Nr_Neigh
real(PREC), intent(in)    :: dx, RX(Nr_atoms), RY(Nr_atoms), RZ(Nr_atoms), LBox(3)
integer, intent(in)       :: H_INDEX_START(Nr_atoms), H_INDEX_END(Nr_atoms)
character(10), intent(in) :: Element_Type(Nr_atoms)
real(PREC), intent(out)   :: dHx(HDIM,HDIM), dHy(HDIM,HDIM), dHz(HDIM,HDIM)
character(10)             :: Type_pair(2)
real(PREC)                :: Rax_p(3), Rax_m(3), Ray_p(3), Ray_m(3), Raz_p(3), Raz_m(3)
real(PREC)                :: Ra(3), Rb(3), diagonal(2), h_0(4,4), Es, Ep, U
real(PREC)                :: dh0x_p(4,4), dh0x_m(4,4), dh0y_p(4,4), dh0y_m(4,4),dh0z_p(4,4), dh0z_m(4,4)
real(PREC)                :: fss_sigma(14),fsp_sigma(14),fps_sigma(14),fpp_sigma(14),fpp_pi(14)
real                      :: Start, Finish, t1, t2;

integer,    intent(in)     :: nrnnlist(Nr_atoms), nnType(Nr_atoms,Max_Nr_Neigh)
real(PREC), intent(in)     :: nnRx(Nr_atoms,Max_Nr_Neigh)
real(PREC), intent(in)     :: nnRy(Nr_atoms,Max_Nr_Neigh), nnRz(Nr_atoms,Max_Nr_Neigh)


integer :: IDim,JDim, IJ
integer :: I,J,K,II,JJ,II_H,JJ_H

dHx = ZERO
dHy = ZERO
dHz = ZERO
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,IJ,II,JJ,II_H,JJ_H,Ra,Rb,fss_sigma,fsp_sigma,fps_sigma,fpp_sigma,fpp_pi) &
!$OMP PRIVATE(Type_pair,diagonal,IDim,JDim,h_0,Rax_p,Rax_m,Ray_p,Ray_m,Raz_p,Raz_m,dh0x_p,dh0x_m) &
!$OMP PRIVATE(dh0y_p,dh0y_m,dh0z_p,dh0z_m)
do I = 1,Nr_atoms
  Type_pair(1) = Element_Type(I)
  Rax_p(1:3) = [RX(I)+dx,RY(I),RZ(I)] 
  Rax_m(1:3) = [RX(I)-dx,RY(I),RZ(I)]
  Ray_p(1:3) = [RX(I),RY(I)+dx,RZ(I)] 
  Ray_m(1:3) = [RX(I),RY(I)-dx,RZ(I)]
  Raz_p(1:3) = [RX(I),RY(I),RZ(I)+dx] 
  Raz_m(1:3) = [RX(I),RY(I),RZ(I)-dx]
  IDim = H_INDEX_END(I)-H_INDEX_START(I)+1
  do J = 1,nrnnlist(I)
    IJ = nnType(I,J)
    if (IJ.ne.I) then
      Type_pair(2) = Element_Type(IJ)
      Rb(1) = nnRx(I,J)
      Rb(2) = nnRy(I,J)
      Rb(3) = nnRz(I,J)
      JDim = H_INDEX_END(IJ)-H_INDEX_START(IJ)+1
      diagonal(1:2) = [1,1]
      call LoadBondIntegralParameters_H(fss_sigma,fsp_sigma,fps_sigma,fpp_sigma,fpp_pi,Es,Ep,U,Type_pair)

      call Slater_Koster_Pair(dh0x_p,Rax_p,Rb,LBox,Type_pair,fss_sigma,fsp_sigma,fps_sigma,fpp_sigma,fpp_pi,diagonal)
      call Slater_Koster_Pair(dh0x_m,Rax_m,Rb,LBox,Type_pair,fss_sigma,fsp_sigma,fps_sigma,fpp_sigma,fpp_pi,diagonal)

      call Slater_Koster_Pair(dh0y_p,Ray_p,Rb,LBox,Type_pair,fss_sigma,fsp_sigma,fps_sigma,fpp_sigma,fpp_pi,diagonal)
      call Slater_Koster_Pair(dh0y_m,Ray_m,Rb,LBox,Type_pair,fss_sigma,fsp_sigma,fps_sigma,fpp_sigma,fpp_pi,diagonal)

      call Slater_Koster_Pair(dh0z_p,Raz_p,Rb,LBox,Type_pair,fss_sigma,fsp_sigma,fps_sigma,fpp_sigma,fpp_pi,diagonal)
      call Slater_Koster_Pair(dh0z_m,Raz_m,Rb,LBox,Type_pair,fss_sigma,fsp_sigma,fps_sigma,fpp_sigma,fpp_pi,diagonal)

      do II = 1,IDim
        II_H = H_INDEX_START(I) + II - 1
        do JJ = 1,JDim
          JJ_H = H_INDEX_START(IJ) + JJ - 1
          dHx(II_H,JJ_H) = dHx(II_H,JJ_H) + (dh0x_p(II,JJ)-dh0x_m(II,JJ))/(2*dx)
          dHy(II_H,JJ_H) = dHy(II_H,JJ_H) + (dh0y_p(II,JJ)-dh0y_m(II,JJ))/(2*dx)
          dHz(II_H,JJ_H) = dHz(II_H,JJ_H) + (dh0z_p(II,JJ)-dh0z_m(II,JJ))/(2*dx)
        enddo
      enddo

    endif
  enddo
enddo
!$OMP END PARALLEL DO

end subroutine Get_dH
