subroutine Ewald_k_Space(COULOMBV,FCOUL,RX,RY,RZ,LBox,DELTAQ,Nr_atoms,COULACC,TIMERATIO,Max_Nr_Neigh)

implicit none

integer, parameter        :: PREC = 8
real(PREC), parameter     :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0, SIX = 6.D0, THREE = 3.D0, FOUR = 4.D0
real(PREC), parameter     :: FOURTYEIGHT = 48.D0, ELEVEN = 11.D0, SIXTEEN = 16.D0, EIGHT = 8.D0
real(PREC), parameter     :: pi = 3.14159265358979323846264D0
real(PREC), parameter     :: SQRTPI = 1.772453850905516D0
integer,    intent(in)    :: Nr_atoms, Max_Nr_Neigh
real(PREC), intent(in)    :: COULACC, TIMERATIO
real(PREC)                :: KECONST, TFACT, RELPERM
real(PREC), intent(in)    :: RX(Nr_atoms), RY(Nr_atoms), RZ(Nr_atoms), LBox(3), DELTAQ(Nr_atoms)
real(PREC)                :: COULCUT, COULCUT2
real(PREC), intent(out)   :: COULOMBV(Nr_atoms), FCOUL(3,Nr_atoms)
real(PREC)                :: Ra(3), Rb(3), dR, Rab(3)
real(PREC)                :: COULVOL, SQRTX, CALPHA, DC(3), MAGR, MAGR2, Z
real(PREC)                :: CALPHA2, TI,TI2,TI3,TI4,TI6,SSA,SSB,SSC,SSD,SSE
real(PREC)                :: CORRFACT,FOURCALPHA2, FORCE
real(PREC)                :: RECIPVECS(3,3),SINLIST(Nr_atoms),COSLIST(Nr_Atoms)
real(PREC)                :: K(3),L11,L12,L13,M21,M22,M23,K2,KCUTOFF,KCUTOFF2,PREFACTOR
real(PREC)                :: PREVIR, COSSUM,SINSUM,DOT, KEPREF, COSSUM2, SINSUM2

integer                   :: I,J,L,M,N, ccnt, nnI, LMAX,MMAX,NMAX,NMIN,MMIN

COULVOL = LBox(1)*LBox(2)*LBox(3)
SQRTX = sqrt(-log(COULACC))

ccnt = 0

CALPHA = sqrtpi*((TIMERATIO*Nr_atoms/(COULVOL*COULVOL))**(ONE/SIX))
COULCUT = SQRTX/CALPHA
CALPHA2 =   CALPHA*CALPHA;
if (COULCUT.GT.50.D0) then
  COULCUT = 50.D0
  CALPHA = SQRTX/COULCUT
endif

!COULCUT = 12.0D0
!CALPHA = SQRTX/COULCUT

COULCUT2 = COULCUT*COULCUT
KCUTOFF = TWO*CALPHA*SQRTX
KCUTOFF2 = KCUTOFF*KCUTOFF
CALPHA2 = CALPHA*CALPHA
FOURCALPHA2 = FOUR*CALPHA2

RECIPVECS = ZERO
RECIPVECS(1,1) = TWO*pi/LBox(1)
RECIPVECS(2,2) = TWO*pi/LBox(2)
RECIPVECS(3,3) = TWO*pi/LBox(3)
LMAX = floor(KCUTOFF / sqrt(RECIPVECS(1,1)*RECIPVECS(1,1)))
MMAX = floor(KCUTOFF / sqrt(RECIPVECS(2,2)*RECIPVECS(2,2)))
NMAX = floor(KCUTOFF / sqrt(RECIPVECS(3,3)*RECIPVECS(3,3)))
!write(*,*) ' Ewald k space Hej 2'

RELPERM = 1.D0
KECONST = 14.3996437701414D0*RELPERM

FCOUL = ZERO
COULOMBV = ZERO
SINLIST = ZERO
COSLIST = ZERO

do L = 0,LMAX

  if (L.eq.0) then
    MMIN = 0
  else
    MMIN = -MMAX
  endif

  L11 = L*RECIPVECS(1,1)
  L12 = L*RECIPVECS(1,2)
  L13 = L*RECIPVECS(1,3)

  do M = MMIN,MMAX

        NMIN = -NMAX
        if ((L==0).and.(M==0)) then
          NMIN = 1
        endif

        M21 = L11 + M*RECIPVECS(2,1)
        M22 = L12 + M*RECIPVECS(2,2)
        M23 = L13 + M*RECIPVECS(2,3)

        do N = NMIN,NMAX
           K(1) = M21 + N*RECIPVECS(3,1)
           K(2) = M22 + N*RECIPVECS(3,2)
           K(3) = M23 + N*RECIPVECS(3,3)
           K2 = K(1)*K(1) + K(2)*K(2) + K(3)*K(3)
           if (K2.le.KCUTOFF2) then
              PREFACTOR = EIGHT*pi*exp(-K2/(4.D0*CALPHA2))/(COULVOL*K2)
              PREVIR = (2.D0/K2) + (2.D0/(4.D0*CALPHA2));

              COSSUM = 0.D0
              SINSUM = 0.D0

              ! Doing the sin and cos sums
              do I = 1,Nr_atoms
                 DOT = K(1)*RX(I) + K(2)*RY(I) + K(3)*RZ(I)
                 ! We re-use these in the next loop...
                 SINLIST(I) = sin(DOT)
                 COSLIST(I) = cos(DOT)
                 COSSUM = COSSUM + DELTAQ(I)*COSLIST(I)
                 SINSUM = SINSUM + DELTAQ(I)*SINLIST(I)
              enddo
              COSSUM2 = COSSUM*COSSUM
              SINSUM2 = SINSUM*SINSUM

              ! Add up energy and force contributions

              KEPREF = KECONST*PREFACTOR
              do I = 1,Nr_atoms
                 COULOMBV(I) = COULOMBV(I) + KEPREF*(COSLIST(I)*COSSUM + SINLIST(I)*SINSUM)
                 FORCE = KEPREF*DELTAQ(I)*(SINLIST(I)*COSSUM - COSLIST(I)*SINSUM)
                 FCOUL(1,I) = FCOUL(1,I) + FORCE*K(1)
                 FCOUL(2,I) = FCOUL(2,I) + FORCE*K(2)
                 FCOUL(3,I) = FCOUL(3,I) + FORCE*K(3)
              enddo

              KEPREF = KEPREF*(COSSUM2 + SINSUM2)
           endif
        enddo
  enddo
enddo

! Point self energy
CORRFACT = 2.D0*KECONST*CALPHA/SQRTPI;
COULOMBV = COULOMBV - CORRFACT*DELTAQ;

end subroutine Ewald_k_Space
