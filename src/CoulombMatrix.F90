subroutine CoulombMatrix(CC,RX,RY,RZ,LBox,U,Element_Type,Nr_atoms,COULACC,TIMERATIO, &
                         nnRx,nnRy,nnRz,nrnnlist,nnType,HDIM,Max_Nr_Neigh)

implicit none

integer,    parameter     :: PREC = 8
real(PREC), parameter     :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0, SIX = 6.D0, THREE = 3.D0
real(PREC), parameter     :: FOURTYEIGHT = 48.D0, ELEVEN = 11.D0, SIXTEEN = 16.D0
real(PREC), parameter     :: pi = 3.14159265358979323846264D0
real(PREC), parameter     :: SQRTPI = 1.772453850905516D0
integer,    intent(in)    :: Nr_atoms, HDIM, Max_Nr_Neigh
real(PREC), intent(in)    :: COULACC, TIMERATIO
real(PREC)                :: TFACT, RELPERM, KECONST
real(PREC), intent(in)    :: RX(Nr_atoms), RY(Nr_atoms), RZ(Nr_atoms), LBox(3)
real(PREC), intent(in)    :: U(Nr_atoms)
real(PREC)                :: COULCUT, COULCUT2, Coulomb_Pot_Real(Nr_atoms),Coulomb_Force_real(3,Nr_atoms)
real(PREC)                :: Coulomb_Force_k(3,Nr_atoms),Coulomb_Force_Real_I(3)
character(10), intent(in) :: Element_Type(Nr_atoms)
integer,    intent(in)    :: nrnnlist(Nr_atoms), nnType(Nr_atoms,Max_Nr_Neigh)
real(PREC), intent(in)    :: nnRx(Nr_atoms,Max_Nr_Neigh)
real(PREC), intent(in)    :: nnRy(Nr_atoms,Max_Nr_Neigh), nnRz(Nr_atoms,Max_Nr_Neigh)
real(PREC), intent(out)   :: CC(Nr_atoms,Nr_atoms)
real(PREC)                :: COULOMBV, FCOUL(3), Coulomb_Pot(Nr_atoms), Coulomb_Pot_k(Nr_atoms)
real(PREC)                :: Ra(3), Rb(3), dR, Rab(3), Coulomb_Pot_Real_I
real(PREC)                :: COULVOL, SQRTX, CALPHA, DC(3), MAGR, MAGR2, Z
real(PREC)                :: CALPHA2, TI,TI2,TI3,TI4,TI6,SSA,SSB,SSC,SSD,SSE
real(PREC)                :: NUMREP_ERFC, TEST(10*Nr_atoms), CA, FORCE, EXPTI, EXPTJ
real(PREC)                :: TJ,TJ2,TJ3,TJ4,TJ6,TI2MTJ2A,SA,SB,SC,SD,SE,SF
real(PREC)                :: TI2MTJ2, TI2MTI2, TJ2MTI2, dq_J(Nr_atoms)

integer                   :: I,J,K, ccnt, nnI

CC = 0.D0
dq_J = 0.D0
Coulomb_Pot_k = 0.D0

  do J = 1,Nr_atoms
    dq_J(J) = 1.D0
!!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,Coulomb_Pot_Real_I,Coulomb_Force_Real_I)
    do I = 1,Nr_atoms
      call Ewald_Real_Space(Coulomb_Pot_Real_I,Coulomb_Force_Real_I,I,RX,RY,RZ,LBox, &
      dq_J,U,Element_Type,Nr_atoms,COULACC,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,HDIM,Max_Nr_Neigh)
      Coulomb_Pot_Real(I) = Coulomb_Pot_Real_I
    enddo
!!$OMP END PARALLEL DO
    call Ewald_k_Space(Coulomb_Pot_k,Coulomb_Force_k,RX,RY,RZ,LBox,dq_J,Nr_atoms,COULACC,TIMERATIO,Max_Nr_Neigh)
    Coulomb_Pot = Coulomb_Pot_Real+Coulomb_Pot_k
    CC(:,J) = Coulomb_Pot
    dq_J(J) = 0.D0
  enddo
  CC = 0.5d0*(CC+transpose(CC))

end subroutine CoulombMatrix
