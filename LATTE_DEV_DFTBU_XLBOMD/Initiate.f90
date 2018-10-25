subroutine Initiate(Coulomb_acc,TIMERATIO,LBox,RX,RY,RZ,Nr_atoms, &
  H_INDEX_START,H_INDEX_END,COULCUT,Element_Type,Znuc,Mnuc,NrOrb)

implicit none

integer,    parameter       :: PREC = 8
real(PREC), parameter       :: pi = 3.14159265358979323846264D0, ONE = 1.D0, SIX = 6.D0
real(PREC), parameter       :: sqrtpi = 1.772453850905516D0
integer,    intent(in)      :: Nr_atoms
real(PREC), intent(in)      :: Coulomb_acc, TIMERATIO, LBox(3)
real(PREC), intent(out)     :: RX(Nr_atoms), RY(Nr_atoms), RZ(Nr_atoms)
integer,    intent(out)     :: H_INDEX_START(Nr_atoms), H_INDEX_END(Nr_atoms)
real(PREC), intent(out)     :: Znuc(Nr_atoms), Mnuc(Nr_atoms), COULCUT
integer,    intent(out)     :: NrOrb(Nr_atoms)
character(10), intent(out)  :: Element_Type(Nr_atoms)

integer       :: i, CNT
real(PREC)    :: SQRTX, CALPHA, COULVOL

COULVOL = LBox(1)*LBox(2)*LBox(3)
SQRTX = sqrt(-log(Coulomb_acc))

CALPHA = sqrtpi*((TIMERATIO*Nr_atoms/(COULVOL*COULVOL))**(ONE/SIX))
COULCUT = SQRTX/CALPHA
if (COULCUT.GT.50.0D0) then
  COULCUT = 50.D0
endif

open(UNIT=22,STATUS="OLD",FILE="COORD.dat")
do i = 1,Nr_atoms
  read(22,*) Element_Type(i),RX(i),RY(i),RZ(i)
end do
close(12)

  do I = 1,Nr_atoms
    RZ(I) = RZ(I) + .2D0
    if (RX(I) > LBox(1)) then
      RX(I) = RX(I) - LBox(1)
    endif
    if (RY(I) > LBox(2)) then
      RY(I) = RY(I) - LBox(2)
    endif
    if (RZ(I) > LBox(3)) then
      RZ(I) = RZ(I) - LBox(3)
    endif
    if (RX(I) < 0.D0) then
      RX(I) = RX(I) + LBox(1)
    endif
    if (RY(I) < 0.D0) then
      RY(I) = RY(I) + LBox(2)
    endif
    if (RZ(I) < 0.D0) then
      RZ(I) = RZ(I) + LBox(3)
    endif
!    RX(I) = mod(RX(I)*1000000000.D0,1000000000.D0*LBox(1))/1000000000.D0
!    RY(I) = mod(RY(I)*1000000000.D0,1000000000.D0*LBox(2))/1000000000.D0
!    RZ(I) = mod(RZ(I)*1000000000.D0,1000000000.D0*LBox(3))/1000000000.D0
  enddo

CNT = 1
do i = 1,Nr_atoms
  if (Element_Type(i).eq.'H') then
    H_INDEX_START(i) = CNT
    NrOrb(i) = 1   ! a single s orbital
    CNT = CNT+1
    H_INDEX_END(i) = CNT-1
    Znuc(i) = 1.D0  ! For hydrogen
    Mnuc(i) = 1.0079D0
  else
    H_INDEX_START(i) = CNT;
    NrOrb(i) = 4.D0   ! one 1 + three p orbitals
    CNT = CNT+4 
    H_INDEX_END(i) = CNT-1 
    if (Element_Type(i).eq.'O') then
      Znuc(i) = 6.D0  ! For oxygen
      Mnuc(i) = 15.9994D0
    endif
    if (Element_Type(i).eq.'C') then
      Znuc(i) = 4.D0  ! For oxygen
      Mnuc(i) = 12.01D0 
    endif
    if (Element_Type(i).eq.'N') then
      Znuc(i) = 5.D0  ! For oxygen
      Mnuc(i) = 14.0067D0
    endif
  endif
enddo

end subroutine Initiate 
