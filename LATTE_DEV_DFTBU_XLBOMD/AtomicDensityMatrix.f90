subroutine AtomicDensityMatrix(D_atomic,Nr_atoms,H_INDEX_START,H_INDEX_END,HDIM,Znuc) 

implicit none

integer, parameter        :: PREC = 8
real(PREC), parameter     :: ZERO = 0.D0, ONE = 1.D0, TWO = 2.D0, THREE = 3.D0
integer,    intent(in)    :: Nr_atoms, HDIM
integer,    intent(in)    :: H_INDEX_START(Nr_Atoms), H_INDEX_END(Nr_atoms)
real(PREC), intent(out)   :: D_atomic(HDIM)
real(PREC), intent(in)    :: Znuc(Nr_atoms)
real(PREC)                :: OCC

integer                   :: I, INDEX, N_orb

D_atomic = ZERO
INDEX = ZERO
do I = 1,Nr_atoms
  N_orb = H_INDEX_END(I)-H_INDEX_START(I) + 1
  if (N_orb.eq.1) then
    INDEX = INDEX + 1
    D_atomic(INDEX) = Znuc(I)
  else
    if (Znuc(I).le.TWO) then
      INDEX = INDEX + 1
      D_atomic(INDEX) = Znuc(I)

      INDEX = INDEX + 1
      D_atomic(INDEX) = 0
      INDEX = INDEX + 1
      D_atomic(INDEX) = 0
      INDEX = INDEX + 1
      D_atomic(INDEX) = 0
    else
      INDEX = INDEX + 1
      D_atomic(INDEX) = TWO

      INDEX = INDEX + 1
      OCC = (Znuc(I)-TWO)/THREE
      D_atomic(INDEX) = OCC
      INDEX = INDEX + 1
      D_atomic(INDEX) = OCC
      INDEX = INDEX + 1
      D_atomic(INDEX) = OCC
    endif
  endif
enddo

end subroutine AtomicDensityMatrix

