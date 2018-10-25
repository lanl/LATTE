subroutine NonOrthCoulForce(FSCOUL,Nr_atoms,HDIM,dSx,dSy,dSz,H_INDEX_START,H_INDEX_END,D,q,Coulomb_Pot,Hubbard_U)

use omp_lib
implicit none
integer,    parameter     :: PREC = 8
real(PREC), parameter     :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0, SIX = 6.D0
integer,    intent(in)    :: Nr_atoms, HDIM
integer,    intent(in)    :: H_INDEX_START(Nr_atoms), H_INDEX_END(Nr_atoms)
real(PREC), intent(in)    :: D(HDIM,HDIM), q(Nr_atoms), Coulomb_Pot(Nr_atoms), Hubbard_U(Nr_atoms)
real(PREC), intent(in)    :: dSx(HDIM,HDIM), dSy(HDIM,HDIM), dSz(HDIM,HDIM)
real(PREC), intent(out)   :: FSCOUL(3,Nr_atoms)
real(PREC)                :: dXtmp, dYtmp, dZtmp
real(PREC)                :: dDSX(HDIM), dDSY(HDIM), dDSZ(HDIM), dQLxdR, dQLydR, dQLzdR

integer :: IDim,JDim
integer :: I,J,K
integer :: Ia, iq, Ia_1, Ia_A, Ia_B, j_a, j_b, threads, id

FSCOUL = ZERO
dDSX = ZERO
dDSY = ZERO
dDSZ = ZERO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(Ia,Ia_A,Ia_B,iq,dXtmp,dYtmp,dZtmp,I,dDSX,dDSY,dDsZ,J,j_a,j_b) &
!$OMP PRIVATE(dQLxdR,dQLydR,dQLzdR)
do Ia = 1,Nr_atoms ! Derivatives Ra
  Ia_A = H_INDEX_START(Ia)
  Ia_B = H_INDEX_END(Ia)
  do iq = 1,HDIM
    dXtmp = ZERO
    dYtmp = ZERO
    dZtmp = ZERO
    do I = Ia_A,Ia_B
      dXtmp = dXtmp + D(iq,I)*dSx(I,iq)
      dYtmp = dYtmp + D(iq,I)*dSy(I,iq)
      dZtmp = dZtmp + D(iq,I)*dSz(I,iq)
    enddo
    dDSX(iq) = dXtmp
    dDSY(iq) = dYtmp
    dDSZ(iq) = dZtmp
  enddo
  do iq = Ia_A,Ia_B
    dXtmp = ZERO
    dYtmp = ZERO
    dZtmp = ZERO
    do I = 1,HDIM
      dXtmp = dXtmp + D(iq,I)*dSx(iq,I)
      dYtmp = dYtmp + D(iq,I)*dSy(iq,I)
      dZtmp = dZtmp + D(iq,I)*dSz(iq,I)
    enddo
    dDSX(iq) = dDSX(iq) + dXtmp
    dDSY(iq) = dDSY(iq) + dYtmp
    dDSZ(iq) = dDSZ(iq) + DZtmp
  enddo
  do J = 1,Nr_atoms ! Get the Mulliken charges for all atoms
    j_a = H_INDEX_START(J)
    j_b = H_INDEX_END(J)
    dQLxdR = ZERO
    dQLydR = ZERO
    dQLzdR = ZERO
    do K = j_a,j_b 
      dQLxdR = dQLxdR + dDSX(K)  ! Derivative with respect to Ia of charge on atom J
      dQLydR = dQLydR + dDSY(K)  ! Derivative with respect to Ia of charge on atom J
      dQLzdR = dQLzdR + dDSZ(K)  ! Derivative with respect to Ia of charge on atom J
    enddo
    FSCOUL(1,Ia) = FSCOUL(1,Ia) - dQLxdR*(Hubbard_U(J)*q(J) + Coulomb_Pot(J))
    FSCOUL(2,Ia) = FSCOUL(2,Ia) - dQLydR*(Hubbard_U(J)*q(J) + Coulomb_Pot(J))
    FSCOUL(3,Ia) = FSCOUL(3,Ia) - dQLzdR*(Hubbard_U(J)*q(J) + Coulomb_Pot(J))
  enddo
enddo
!$OMP END PARALLEL DO

end subroutine NonOrthCoulForce
