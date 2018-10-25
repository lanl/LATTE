subroutine PairPotential_Forces(PairForces,ERep,RX,RY,RZ,LBox,Element_Type,Nr_atoms,PairPotType1,PairPotType2,PotCoef, &
     Max_Nr_Neigh,nrnnlist,nnRx,nnRy,nnRz,nnType)

use omp_lib
implicit none
integer,       parameter   :: PREC = 8
integer,       intent(in)  :: Nr_atoms, Max_Nr_Neigh
real(PREC),    parameter   :: ONE = 1, TWO = 2, ZERO = 0, SIX = 6
real(PREC),    intent(out) :: PairForces(3,Nr_atoms), ERep
real(PREC),    intent(in)  :: RX(Nr_atoms), RY(Nr_atoms), RZ(Nr_atoms), LBox(3), PotCoef(16,10)
character(10), intent(in)  :: Element_Type(Nr_atoms)
character(1),  intent(in)  :: PairPotType1(10), PairPotType2(10)

real(PREC)                :: UNIVPHI, CUTPHI, VIRUNIV, VIRCUT, FUNIV(3), FCUT(3), R1, MYR
real(PREC)                :: RXb, RYb, RZb, RCUT, RCUT2, Ra(3), Rb(3), Rab(3), dR, dR2
real(PREC)                :: POLYNOM, PHI, DPOLYNOM, DPHI(3), EXPTMP, DC(3), FTMP(3), FORCE
integer                   :: i, j, k, PPSEL, nr_shift_X, nr_shift_Y, nr_shift_Z, cnt, ij

integer,    intent(in)     :: nrnnlist(Nr_atoms), nnType(Nr_atoms,Max_Nr_Neigh)
real(PREC), intent(in)     :: nnRx(Nr_atoms,Max_Nr_Neigh)
real(PREC), intent(in)     :: nnRy(Nr_atoms,Max_Nr_Neigh), nnRz(Nr_atoms,Max_Nr_Neigh)
real                       :: t1, t2

UNIVPHI = ZERO
CUTPHI = ZERO

VIRUNIV = ZERO
VIRCUT =  ZERO
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,FUNIV,FCUT,Ra,j,ij,k,PPSEL,R1,RCUT,RCUT2,Rb,Rab,dR,dR2) &
!$OMP PRIVATE(DC,POLYNOM,PHI,DPOLYNOM,DPHI,EXPTMP,FTMP,MYR,FORCE) &
!$OMP REDUCTION(+:UNIVPHI,CUTPHI)
do i = 1,Nr_atoms
  FUNIV = ZERO
  FCUT = ZERO
  Ra(1) = RX(i) 
  Ra(2) = RY(i) 
  Ra(3) = RZ(i)
  do j = 1,nrnnlist(i)
  ij = nnType(i,j)
  if (i.ne.ij) then
    do k = 1,10
      if (((Element_Type(i) == PairPotType1(k)).and.(Element_Type(ij) == PairPotType2(k))).or. &
      (Element_Type(ij) == PairPotType1(k)).and.(Element_Type(i) == PairPotType2(k))) then
        PPSEL = k
        R1 = PotCoef(9,PPSEL)
        RCUT = PotCoef(10,PPSEL)
        RCUT2 = RCUT*RCUT
      endif
    enddo

    Rb(1) = nnRx(i,j) 
    Rb(2) = nnRy(i,j) 
    Rb(3) = nnRz(i,j) 

    Rab = Rb-Ra  ! OBS b - a !!!
    dR = norm2(Rab) 
    dR2 = dR*dR

    if (dR < RCUT) then
       DC = Rab/dR
       if (dR < R1) then
          POLYNOM = dR*(PotCoef(2,PPSEL) + dR*(PotCoef(3,PPSEL) + dR*(PotCoef(4,PPSEL) + dR*PotCoef(5,PPSEL))))
          PHI = PotCoef(1,PPSEL)*exp(POLYNOM)
          DPOLYNOM = PotCoef(2,PPSEL) + dR*(2*PotCoef(3,PPSEL) + dR*(3*PotCoef(4,PPSEL) + 4*PotCoef(5,PPSEL)*dR))
          DPHI = -DC*PHI*DPOLYNOM
          EXPTMP = PotCoef(6,PPSEL)*exp(PotCoef(7,PPSEL)*(dR - PotCoef(8,PPSEL)))
          UNIVPHI = UNIVPHI + PHI + EXPTMP
          FTMP = DC*PotCoef(7,PPSEL)*EXPTMP
          FUNIV = FUNIV - DPHI + FTMP
       else
          MYR = dR - R1
          CUTPHI =  CUTPHI + PotCoef(11,PPSEL) + MYR*(PotCoef(12,PPSEL) + MYR*(PotCoef(13,PPSEL) &
           + MYR*(PotCoef(14,PPSEL) + MYR*(PotCoef(15,PPSEL) + MYR*PotCoef(16,PPSEL)))))
          FORCE = PotCoef(12,PPSEL)  + MYR*(2*PotCoef(13,PPSEL) + MYR*(3*PotCoef(14,PPSEL) &
           + MYR*(4*PotCoef(15,PPSEL) + MYR*5*PotCoef(16,PPSEL))))
          FCUT = FCUT + DC*FORCE
       endif
    endif
  endif
  enddo
  PairForces(:,i) = FUNIV + FCUT
enddo
!$OMP END PARALLEL DO
ERep = 0.5*(UNIVPHI + CUTPHI)

end subroutine PairPotential_Forces

