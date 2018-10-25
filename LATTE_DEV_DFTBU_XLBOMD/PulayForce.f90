subroutine PulayForce(FPUL,Nr_atoms,HDIM,Z,H,D,dSx,dSy,dSz,H_INDEX_START,H_INDEX_END) ! from 2 tr[Z*Z'*F*dS/dR]

use omp_lib
implicit none
integer,    parameter     :: PREC = 8
real(PREC), parameter     :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0, SIX = 6.D0
integer,    intent(in)    :: Nr_atoms, HDIM
real(PREC), intent(in)    :: D(HDIM,HDIM), H(HDIM,HDIM), Z(HDIM,HDIM)
integer,    intent(in)    :: H_INDEX_START(Nr_atoms), H_INDEX_END(Nr_atoms)
real(PREC), intent(in)    :: dSx(HDIM,HDIM), dSy(HDIM,HDIM), dSz(HDIM,HDIM)
real(PREC), intent(out)   :: FPUL(3,Nr_atoms)
real(PREC)                :: Xtmp, Ytmp, Ztmp, SIHD(HDIM,HDIM), SI(HDIM,HDIM), HD(HDIM,HDIM)

integer :: I,J,K, I_A, I_B

FPUL = ZERO
call MMult(ONE,Z,Z,ZERO,SI,'N','T',HDIM)
call MMult(ONE,H,D,ZERO,HD,'N','N',HDIM)
call MMult(TWO,SI,HD,ZERO,SIHD,'N','N',HDIM)
!SIHD = 2*Z*Z'*H*D;  % Pulay Force FPUL from 2Tr[ZZ'HD*dS/dR]
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,K,I_A,I_B,Xtmp,Ytmp,Ztmp)
do I = 1,Nr_atoms
  I_A = H_INDEX_START(I);
  I_B = H_INDEX_END(I);
  Xtmp = ZERO
  Ytmp = ZERO
  Ztmp = ZERO
  do J = I_A,I_B
    do K = 1,HDIM
      Xtmp = Xtmp + dSx(J,K)*SIHD(K,J)
      Ytmp = Ytmp + dSy(J,K)*SIHD(K,J)
      Ztmp = Ztmp + dSz(J,K)*SIHD(K,J)
    enddo
  enddo
  FPUL(1:3,I) = [Xtmp,Ytmp,Ztmp]
enddo
!$OMP END PARALLEL DO

end subroutine PulayForce
