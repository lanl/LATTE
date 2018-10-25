subroutine Build_H0(H0,Hubbard_U,Nr_atoms,HDIM,Max_Nr_Neigh,Element_Type,RX,RY,RZ,LBox,H_INDEX_START,H_INDEX_END, &
  nrnnlist,nnRx,nnRy,nnRz,nnType)

use omp_lib
implicit none

integer, parameter         :: PREC = 8
real(PREC), parameter      :: ZERO = 0, ONE = 1, TWO = 2
integer,    intent(in)     :: Nr_atoms, HDIM, Max_Nr_Neigh
integer                    :: IDim, JDim 
real(PREC), intent(in)     :: RX(Nr_atoms), RY(Nr_atoms), RZ(Nr_atoms), LBox(3)
integer,    intent(in)     :: H_INDEX_START(Nr_atoms), H_INDEX_END(Nr_atoms)
real(PREC), intent(out)    :: Hubbard_U(Nr_atoms)
real(PREC), intent(out)    :: H0(HDIM,HDIM)
real(PREC)                 :: Ra(3), Rb(3), diagonal(2), h_0(4,4)
real(PREC)                 :: fss_sigma(14),fsp_sigma(14),fps_sigma(14),fpp_sigma(14),fpp_pi(14)
real(PREC)                 :: Es, Ep, U

integer,    intent(in)     :: nrnnlist(Nr_atoms), nnType(Nr_atoms,Max_Nr_Neigh)
real(PREC), intent(in)     :: nnRx(Nr_atoms,Max_Nr_Neigh)
real(PREC), intent(in)     :: nnRy(Nr_atoms,Max_Nr_Neigh), nnRz(Nr_atoms,Max_Nr_Neigh)

character(10), intent(in)  :: Element_Type(Nr_atoms)
character(10)              :: Type_pair(2)

integer                    :: I,J,II,JJ,II_H,JJ_H,nthreads

nthreads = 8
call omp_set_num_threads(nthreads)

H0 = ZERO  ! Replace by sparse ELLPACK matrix
h_0 = ZERO

!!$OMP PARALLEL DO FIRSTPRIVATE(v_row_index,v_tmp,Row_index_Z_tmp) &
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,II,JJ,II_H,JJ_H,Ra,Rb,fss_sigma,fsp_sigma,fps_sigma,fpp_sigma,fpp_pi) & 
!$OMP PRIVATE(Es,Ep,U,Type_pair,diagonal,IDim,JDim,h_0)
do I = 1,Nr_atoms  ! OpenMP threading
  Type_pair(1) = Element_Type(I)
  Ra(1) = RX(I)
  Ra(2) = RY(I)
  Ra(3) = RZ(I)
  IDim = H_INDEX_END(I)-H_INDEX_START(I)+1

  do J = 1,nrnnlist(I)
    Type_pair(2) = Element_Type(nnType(I,J))
    Rb(1) = nnRx(I,J)
    Rb(2) = nnRy(I,J)
    Rb(3) = nnRz(I,J)
    JDim = H_INDEX_END(nnType(I,J))-H_INDEX_START(nnType(I,J))+1

! Hamiltonian block for a-b atom pair
    call LoadBondIntegralParameters_H(fss_sigma,fsp_sigma,fps_sigma,fpp_sigma,fpp_pi,Es,Ep,U,Type_pair) ! Used in BondIntegral(dR,fxx_xx)
    if (I.eq.nnType(I,J)) then
      Hubbard_U(I) = U
    endif
    diagonal(1) = Es
    diagonal(2) = Ep
    !call Slater_Koster_Block(h_0,Ra,Rb,LBox,Type_pair,fss_sigma,fsp_sigma,fps_sigma,fpp_sigma,fpp_pi,diagonal)
    call Slater_Koster_Pair(h_0,Ra,Rb,LBox,Type_pair,fss_sigma,fsp_sigma,fps_sigma,fpp_sigma,fpp_pi,diagonal)
    do II = 1,IDim
      II_H = H_INDEX_START(I) + II - 1
      do JJ = 1,JDim
        JJ_H = H_INDEX_START(nnType(I,J)) + JJ - 1
        H0(II_H,JJ_H) = h_0(II,JJ)  ! Replace by sparse ELLPACK matrix
        H0(JJ_H,II_H) = h_0(II,JJ)
      enddo
    enddo
  enddo

enddo
!$OMP END PARALLEL DO

end subroutine Build_H0
