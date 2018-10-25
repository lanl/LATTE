subroutine SlaterKosterForce(SKForce,Nr_atoms,HDIM,D,dHx,dHy,dHz,H_INDEX_START,H_INDEX_END)

use omp_lib
implicit none
integer,    parameter     :: PREC = 8
real(PREC), parameter     :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0, SIX = 6.D0
integer,    intent(in)    :: Nr_atoms, HDIM
real(PREC), intent(in)    :: D(HDIM,HDIM)
integer,    intent(in)    :: H_INDEX_START(Nr_atoms), H_INDEX_END(Nr_atoms)
real(PREC), intent(in)    :: dHx(HDIM,HDIM), dHy(HDIM,HDIM), dHz(HDIM,HDIM)
real(PREC), intent(out)   :: SKForce(3,Nr_atoms)
real(PREC)                :: Xtmp, Ytmp, Ztmp
integer                   :: I,J,K, I_A, I_B

SKForce = ZERO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,K,I_A,I_B,Xtmp,Ytmp,Ztmp)
do I = 1,Nr_atoms ! Slater-Koster Force SKForce from Tr[D*dH0/dR]
  I_A = H_INDEX_START(I)
  I_B = H_INDEX_END(I)
  Xtmp = ZERO
  Ytmp = ZERO
  Ztmp = ZERO
  do J = I_A,I_B
    do K = 1,HDIM
      Xtmp = Xtmp + dHx(J,K)*D(K,J);
      Ytmp = Ytmp + dHy(J,K)*D(K,J);
      Ztmp = Ztmp + dHz(J,K)*D(K,J);
    enddo
  enddo
  SKForce(1:3,I) = -2*[Xtmp,Ytmp,Ztmp]
enddo
!$OMP END PARALLEL DO

end subroutine SlaterKosterForce

