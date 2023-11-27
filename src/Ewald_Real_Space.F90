subroutine  Ewald_Real_Space(COULOMBV,FCOUL,I,RX,RY,RZ,LBox, &
    DELTAQ,U,Element_Type,Nr_atoms,COULACC,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,HDIM,Max_Nr_Neigh)


implicit none

integer,    parameter     :: PREC = 8
real(PREC), parameter     :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0, SIX = 6.D0, THREE = 3.D0
real(PREC), parameter     :: FOURTYEIGHT = 48.D0, ELEVEN = 11.D0, SIXTEEN = 16.D0
real(PREC), parameter     :: pi = 3.14159265358979323846264D0
real(PREC), parameter     :: SQRTPI = 1.772453850905516D0
integer,    intent(in)    :: Nr_atoms, HDIM, Max_Nr_Neigh, I
real(PREC), intent(in)    :: COULACC, TIMERATIO
real(PREC)                :: TFACT, RELPERM, KECONST
real(PREC), intent(in)    :: RX(Nr_atoms), RY(Nr_atoms), RZ(Nr_atoms), LBox(3), DELTAQ(Nr_atoms)
real(PREC), intent(in)    :: U(Nr_atoms)
real(PREC)                :: COULCUT, COULCUT2
character(10), intent(in) :: Element_Type(Nr_atoms)
integer,    intent(in)    :: nrnnlist(Nr_atoms), nnType(Nr_atoms,Max_Nr_Neigh)
real(PREC), intent(in)    :: nnRx(Nr_atoms,Max_Nr_Neigh)
real(PREC), intent(in)    :: nnRy(Nr_atoms,Max_Nr_Neigh), nnRz(Nr_atoms,Max_Nr_Neigh)
real(PREC), intent(out)   :: COULOMBV, FCOUL(3)
real(PREC)                :: Ra(3), Rb(3), dR, Rab(3)
real(PREC)                :: COULVOL, SQRTX, CALPHA, DC(3), MAGR, MAGR2, Z
real(PREC)                :: CALPHA2, TI,TI2,TI3,TI4,TI6,SSA,SSB,SSC,SSD,SSE
real(PREC)                :: NUMREP_ERFC, TEST(10*Nr_atoms), CA, FORCE, EXPTI, EXPTJ
real(PREC)                :: TJ,TJ2,TJ3,TJ4,TJ6,TI2MTJ2A,SA,SB,SC,SD,SE,SF
real(PREC)                :: TI2MTJ2, TI2MTI2, TJ2MTI2

integer                   :: J,K, ccnt, nnI

COULVOL = LBox(1)*LBox(2)*LBox(3)
SQRTX = sqrt(-log(COULACC))

ccnt = 0
!COULCUT = 12.D0
!CALPHA = SQRTX/COULCUT
!COULCUT2 = COULCUT*COULCUT
!CALPHA2 = CALPHA*CALPHA

CALPHA = sqrtpi*((TIMERATIO*Nr_atoms/(COULVOL*COULVOL))**(ONE/SIX))
COULCUT = SQRTX/CALPHA
CALPHA2 =   CALPHA*CALPHA;
if (COULCUT.GT.50.0D0) then
  COULCUT = 50.D0
  CALPHA = SQRTX/COULCUT
endif
COULCUT2 = COULCUT*COULCUT
CALPHA2 = CALPHA*CALPHA

RELPERM = ONE
KECONST = 14.3996437701414D0*RELPERM
TFACT  = 16.0D0/(5.0D0*KECONST)


FCOUL = ZERO
COULOMBV = ZERO

FCOUL = ZERO
COULOMBV = 0.D0

TI = TFACT*U(I)
TI2 = TI*TI
TI3 = TI2*TI
TI4 = TI2*TI2
TI6 = TI4*TI2

SSA = TI
SSB = TI3/48.D0
SSC = 3.D0*TI2/16.D0
SSD = 11.D0*TI/16.D0
SSE = 1.D0

Ra(1) = RX(I)
Ra(2) = RY(I)
Ra(3) = RZ(I)

do nnI = 1,nrnnlist(I)
  Rb(1) = nnRx(I,nnI)
  Rb(2) = nnRy(I,nnI)
  Rb(3) = nnRz(I,nnI)
  J = nnType(I,nnI)

    Rab = Rb-Ra  ! OBS b - a !!!
    dR = norm2(Rab)
    MAGR = dR
    MAGR2 = dR*dR

    if ((dR <= COULCUT).and.(dR > 1e-12)) then

       TJ = TFACT*U(J)
       DC = Rab/dR

       ! Not Using Numerical Recipes ERFC
       Z = abs(CALPHA*MAGR)
       NUMREP_ERFC = erfc(Z)

       CA = NUMREP_ERFC/MAGR
       COULOMBV = COULOMBV + DELTAQ(J)*CA
       ccnt = ccnt + 1
       !TEST(ccnt) = DELTAQ(J)*CA
       CA = CA + TWO*CALPHA*exp( -CALPHA2*MAGR2 )/SQRTPI
       FORCE = -KECONST*DELTAQ(I)*DELTAQ(J)*CA/MAGR
       EXPTI = exp(-TI*MAGR )

       if (Element_Type(I).eq.Element_Type(J)) then
           COULOMBV = COULOMBV - DELTAQ(J)*EXPTI*(SSB*MAGR2 + SSC*MAGR + SSD + SSE/MAGR)
           ccnt = ccnt + 1
           !TEST(ccnt) = - DELTAQ(J)*EXPTI*(SSB*MAGR2 + SSC*MAGR + SSD + SSE/MAGR)
           FORCE = FORCE + (KECONST*DELTAQ(I)*DELTAQ(J)*EXPTI)*((SSE/MAGR2 - TWO*SSB*MAGR - SSC) &
             + SSA*(SSB*MAGR2 + SSC*MAGR + SSD + SSE/MAGR))
       else
           TJ2 = TJ*TJ
           TJ3 = TJ2*TJ
           TJ4 = TJ2*TJ2
           TJ6 = TJ4*TJ2
           EXPTJ = exp( -TJ*MAGR )
           TI2MTJ2 = TI2 - TJ2
           TJ2MTI2 = -TI2MTJ2
           SA = TI
           SB = TJ4*TI/(TWO*TI2MTJ2*TI2MTJ2)
           SC = (TJ6 - THREE*TJ4*TI2)/(TI2MTJ2*TI2MTJ2*TI2MTJ2)
           SD = TJ
           SE = TI4*TJ/(TWO * TJ2MTI2 * TJ2MTI2)
           SF = (TI6 - THREE*TI4*TJ2)/(TJ2MTI2*TJ2MTI2*TJ2MTI2)

           COULOMBV = COULOMBV - (DELTAQ(J)*(EXPTI*(SB - (SC/MAGR)) + EXPTJ*(SE - (SF/MAGR))))
           FORCE = FORCE + KECONST*DELTAQ(I)*DELTAQ(J)*((EXPTI*(SA*(SB - (SC/MAGR)) - (SC/MAGR2))) &
              + (EXPTJ*(SD*(SE - (SF/MAGR)) - (SF/MAGR2))))
       endif

       FCOUL(1) = FCOUL(1) + DC(1)*FORCE
       FCOUL(2) = FCOUL(2) + DC(2)*FORCE
       FCOUL(3) = FCOUL(3) + DC(3)*FORCE
    endif
enddo
COULOMBV = KECONST*COULOMBV

end subroutine Ewald_Real_Space
