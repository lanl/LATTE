!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    LATTE-XLBOMD Scaled down developers version of LATTE    !!!
!!!    Anders M.N. Niklasson et al. amn@lanl.gov               !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program XLBOMD

use omp_lib
implicit none
integer, parameter       :: PREC = 8, SEED = 78739

!integer, parameter       :: Nr_atoms = 30, HDIM = 60, Nocc = 40, Max_Nr_Neigh = 500   ! 10 water molecules
!real(PREC), parameter	  :: Lx = 5.0000000D0, Ly = 5.0000000D0, Lz = 5.0000000D0

!integer, parameter       :: Nr_atoms = 24, HDIM = 48, Nocc = 32, Max_Nr_Neigh = 500   ! 8 water molecules
!real(PREC), parameter	  :: Lx = 8.2670000D0, Ly = 8.2670000D0, Lz = 8.2670000D0

!integer, parameter       :: Nr_atoms = 300, HDIM = 600, Nocc = 400, Max_Nr_Neigh = 500   ! 100 water molecules
!real(PREC), parameter	  :: Lx = 14.404799, Ly = 14.404799, Lz = 14.404799

!integer, parameter       :: Nr_atoms = 55, HDIM = 220, Nocc = 110, Max_Nr_Neigh = 500   ! 55 Amorph C
!real(PREC), parameter    :: Lx = 10.0000000, Ly = 10.0000000, Lz = 10.0000000
!!real(PREC), parameter    :: T = 2000.00D0, dt = 1.00D0, C_mix = 0.16D0 ! Linear mix
!!!real(PREC), parameter    :: T = 1500.00D0, dt = 1.00D0, C_mix = 0.50D0, C_DFTB_U = 10.0D0
!real(PREC), parameter    :: T = 500.00D0, dt = 1.00D0, C_mix = 0.05D0, C_DFTB_U = 00.0D0
!integer, parameter       :: Full_Kernel_Initial = 0, Full_Kernel_MD = 0

!integer, parameter       :: Nr_atoms = 130, HDIM = 520, Nocc = 268, Max_Nr_Neigh = 500   ! non-PGM
!real(PREC), parameter    :: Lx = 9.6860000, Ly = 17.0480000, Lz = 20.0000000
!real(PREC), parameter    :: Lx = 20.0000, Ly = 17.0480000, Lz = 9.6860000

!integer, parameter       :: Nr_atoms = 3, HDIM = 6, Nocc = 4, Max_Nr_Neigh = 500   ! 1 Water
!real(PREC), parameter    :: Lx = 8.2670000D0, Ly = 8.2670000D0, Lz = 8.2670000D0
!real(PREC), parameter    :: T = 100.00D0, dt = 0.25D0, C_mix = 0.10D0, C_DFTB_U = 00.0D0
!integer, parameter       :: Full_Kernel_Initial = 0, Full_Kernel_MD = 0

!integer, parameter       :: Nr_atoms = 7, HDIM = 19, Nocc = 12, Max_Nr_Neigh = 500   ! 7 NM
!real(PREC), parameter    :: Lx = 12.1151720, Ly = 12.1151720, Lz = 12.1151720
!real(PREC), parameter    :: T = 100.00D0, dt = 0.25D0, C_mix = 0.15D0, C_DFTB_U = 2.0D0
!integer, parameter       :: Full_Kernel_Initial = 0, Full_Kernel_MD = 0

integer, parameter       :: Nr_atoms = 49, HDIM = 133, Nocc = 84, Max_Nr_Neigh = 500   ! 49 NM
real(PREC), parameter    :: Lx = 12.1151720, Ly = 12.1151720, Lz = 12.1151720
real(PREC), parameter    :: T = 100.00D0, dt = 0.25D0, C_mix = 0.45D0, C_DFTB_U = 0.0D0  ! 15 Iterations to convergence
integer, parameter       :: Full_Kernel_Initial = 0, Full_Kernel_MD = 0

!integer, parameter       :: Nr_atoms = 70, HDIM = 190, Nocc = 120, Max_Nr_Neigh = 500   ! 70 NM
!real(PREC), parameter    :: Lx = 9.6158180, Ly = 9.6158180, Lz = 9.6158180

!integer, parameter       :: Nr_atoms = 99, HDIM = 396, Nocc = 198, Max_Nr_Neigh = 500   ! 99 Graphene
!real(PREC), parameter    :: Lx = 14.216000, Ly = 21.274000, Lz = 12.361400
!real(PREC), parameter    :: T = 10000.00D0, dt = 0.5D0, C_mix = 0.05D0
!integer, parameter       :: Full_Kernel_Initial = 0, Full_Kernel_MD = 0

!integer, parameter       :: Nr_atoms = 240, HDIM = 600, Nocc = 300, Max_Nr_Neigh = 500   ! Benzene
!real(PREC), parameter    :: Lx = 11.078000, Ly = 11.078000, Lz = 11.078000

!!!!!!!!!!!!!!!!!!!!!!!!    
integer, parameter       :: m = 10, MD_Iter = 200, MaxIt = 20
real(PREC), parameter    :: Coulomb_acc = 1e-7, TIMERATIO = 10.D0, ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0, SIX = 6.D0
real(PREC), parameter    :: THREE = 3.D0

real(PREC), parameter    :: eps = 1e-7, SCF_EPS = 1e-6, dx = 0.0001D0  
real(PREC), parameter    :: OccErrLim = 1e-9

real(PREC), parameter    :: SEVEN = 7.D0
real(PREC), parameter    :: F2V = 0.01602176487D0/1.660548782D0
real(PREC), parameter    :: MVV2KE = 166.0538782D0/1.602176487D0
real(PREC), parameter    :: KE2T = 1.D0/0.000086173435D0
real(PREC)               :: RX(Nr_atoms), RY(Nr_atoms), RZ(Nr_atoms), q(Nr_atoms), LBox(3)
integer                  :: H_INDEX_START(Nr_atoms), H_INDEX_END(Nr_atoms)
real(PREC)               :: Znuc(Nr_atoms), Mnuc(Nr_atoms), Hubbard_U(Nr_atoms), COULCUT
real(PREC)               :: DFTB_U(HDIM), SU(HDIM,HDIM), HU(HDIM,HDIM), H_U(HDIM,HDIM), H_1(HDIM,HDIM)
real(PREC)               :: U(HDIM,HDIM)
real(PREC)               :: dEHub,EHub, EHub_0, EHub_p1, EHub_m1, EHub_m2, EHub_p2, EHubF
real(PREC)               :: D_atomic(HDIM), H0(HDIM,HDIM), H(HDIM,HDIM), H1(HDIM,HDIM)
real(PREC)               :: D0(HDIM,HDIM), D_0(HDIM,HDIM), D(HDIM,HDIM), EYE(HDIM,HDIM), II(Nr_atoms,Nr_atoms)
real(PREC)               :: Delta_DS(HDIM,HDIM), DS(HDIM,HDIM), dD_dP(HDIM,HDIM),dDO_dPO(HDIM,HDIM),nDS, ndDO_dU
real(PREC)               :: Delta_DO(HDIM,HDIM), DO(HDIM,HDIM), PO(HDIM,HDIM), PS(HDIM,HDIM)
real(PREC)               :: nDelta_DO(HDIM,HDIM), nDO
real(PREC)               :: DSPS(HDIM,HDIM), PSDS(HDIM,HDIM), PSPS(HDIM,HDIM), DSDS(HDIM,HDIM), dh
real(PREC)               :: dU(HDIM,HDIM),dDS_dU, qx(Nr_atoms), d2PO(HDIM,HDIM), dU_tmp(HDIM,HDIM)
real(PREC)               :: P(HDIM,HDIM), PO_0(HDIM,HDIM), PO_1(HDIM,HDIM), PO_2(HDIM,HDIM), PO_3(HDIM,HDIM)
real(PREC)               :: PO_4(HDIM,HDIM), PO_5(HDIM,HDIM), PO_6(HDIM,HDIM), PO_7(HDIM,HDIM)
real(PREC)               :: S(HDIM,HDIM), SI(HDIM,HDIM), Z(HDIM,HDIM), ZI(HDIM,HDIM), X(HDIM,HDIM), Y(HDIM,HDIM)
real(PREC)               :: Atmp(HDIM,HDIM), Btmp(HDIM,HDIM)
real(PREC)               :: S_p1(HDIM,HDIM), S_m1(HDIM,HDIM)
real(PREC)               :: alpha, beta, hh(HDIM), mu0, mu1, XI(HDIM,HDIM)
integer                  :: NrOrb(Nr_atoms), step
character(10)            :: Element_Type(Nr_atoms)
character(1)             :: TA, TB, PairPotType1(10), PairPotType2(10)
integer                  :: nrnnlist(Nr_atoms), nnType(Nr_atoms,Max_Nr_Neigh)
integer                  :: nnStruct(Nr_atoms,Nr_atoms), nrnnStruct(Nr_atoms)
real(PREC)               :: nndist(Nr_atoms,Max_Nr_Neigh), nnRx(Nr_atoms,Max_Nr_Neigh)
real(PREC)               :: nnRy(Nr_atoms,Max_Nr_Neigh), nnRz(Nr_atoms,Max_Nr_Neigh)
real(PREC)               :: HOrth(HDIM,HDIM)
real(PREC)               :: SCF_ERR
real(PREC)               :: Coulomb_Pot_Real(Nr_atoms), Coulomb_Force_Real(3,Nr_atoms)
real(PREC)               :: Coulomb_Force_k(3,Nr_atoms), Coulomb_Pot(Nr_atoms), Coulomb_Force(3,Nr_atoms)
real(PREC)               :: Coulomb_Force_Real_I(3), Coulomb_Pot_Real_I, Coulomb_Pot_k(Nr_atoms)
real(PREC)               :: KK(Nr_atoms,Nr_atoms), JJ(Nr_atoms,Nr_atoms)
real(PREC)               :: q_old(Nr_atoms), S_Ent, SS, Coulomb_Pot_a(Nr_atoms), Coulomb_Pot_b(Nr_atoms)
real(PREC)               :: ECoul_a, ECoul_b, q_q(Nr_atoms)
real(PREC)               :: PotCoef(16,10), PairForces(3,Nr_atoms), ERep, ECoul, ECoul0, EBand, EBand0, EEnt, EPOT
real(PREC)               :: dSx(HDIM,HDIM), dSy(HDIM,HDIM), dSz(HDIM,HDIM)
real(PREC)               :: dHx(HDIM,HDIM), dHy(HDIM,HDIM), dHz(HDIM,HDIM)
real(PREC)               :: SKForce(3,Nr_atoms), FPUL(3,Nr_atoms), FSCOUL(3,Nr_atoms), FTOT(3,Nr_atoms)
real(PREC)               :: HubForce(3,Nr_atoms)
real(PREC)               :: n(Nr_atoms),n_0(Nr_atoms),n_1(Nr_atoms),n_2(Nr_atoms),n_3(Nr_atoms),n_4(Nr_atoms)
real(PREC)               :: n_5(Nr_atoms), mu_x, mu_0, mu_1, mu_2, mu_3, mu_4, mu_5
real(PREC)               :: C0, C1, C2, C3, C4, C5, kappa, cc, Time
real(PREC)               :: VX(Nr_atoms), VY(Nr_atoms), VZ(Nr_atoms), EKIN, Energy
real(PREC)               :: dn2dt2(Nr_atoms), Temperature, W(Nr_atoms,Nr_atoms), v(Nr_atoms) 
!real(PREC)               :: u(Nr_atoms)
real(PREC)               :: JJV(Nr_atoms,Nr_atoms), dq_dv(Nr_atoms), Xtmp(Nr_atoms,Nr_atoms), X_Xtmp(Nr_atoms,Nr_atoms)
real(PREC)               :: Res(Nr_atoms), Res_old(Nr_atoms), QQ(Nr_atoms,Nr_atoms), RR(Nr_atoms,Nr_atoms), gq(Nr_atoms)
real(PREC)               :: DM_Res(Nr_atoms), DM_Res_old(Nr_atoms)
real(PREC)               :: QQ_tmp(Nr_atoms,Nr_atoms), RR_tmp(Nr_atoms,Nr_atoms)
real(PREC)               :: QxQ(HDIM,HDIM), ee(HDIM), Fe_vec(HDIM), Err, Rtmp, ET(3), FT(3),dEH(3),dFH(3), EH(3)
integer                  :: I,J,K,L, LL,SCF_IT, MD_step, It_v, It, MEM, MM
real(PREC)               :: norm_v, aa, bb
!real(PREC)               :: H_CHECK(HDIM,HDIM), Z_CHECK(HDIM,HDIM)


LBox(1) = Lx
LBox(2) = Ly
LBox(3) = Lz

EYE = ZERO
do I = 1, HDIM
  EYE(I,I) = ONE
enddo
II = ZERO
do I = 1, Nr_atoms
  II(I,I) = ONE
enddo

!Get the structure and index lists, unit conversion factors, mass and charge etc.
call Initiate(Coulomb_acc,TIMERATIO,LBox,RX,RY,RZ,Nr_atoms, &
  H_INDEX_START,H_INDEX_END,COULCUT,Element_Type,Znuc,Mnuc,NrOrb)

!Rtmp = RX(1)
!do MM = 1,3
!RX(1) = Rtmp + (2.D0-MM)*0.0001D0
!!!! Simple hack for Carbon only
U = ZERO
do I = 1, Nr_atoms   ! Orbital dependent U
   DFTB_U(H_INDEX_START(I)) = 0.D0  ! s-orbitals
   do J = H_INDEX_START(I)+1, H_INDEX_END(I)
      DFTB_U(J) = C_DFTB_U    ! C p-orbitals Hubbard U-J
   enddo
enddo
U = ZERO
do I = 1,HDIM
   U(I,I) = DFTB_U(I)
enddo

! Get the atomic diagonal density matrix
call AtomicDensityMatrix(D_atomic,Nr_atoms,H_INDEX_START,H_INDEX_END,HDIM,Znuc) 

call nearestneighborlist(nrnnlist,nndist,nnRx,nnRy,nnRz,nnType,nnStruct,nrnnStruct,RX,RY,RZ,LBox,4.0D0, &
 Nr_atoms,Max_Nr_Neigh)
call Build_H0(H0,Hubbard_U,Nr_atoms,HDIM,Max_Nr_Neigh,Element_Type,RX,RY,RZ,LBox,H_INDEX_START,H_INDEX_END, &
 nrnnlist,nnRx,nnRy,nnRz,nnType)
call Build_S(S,Nr_atoms,HDIM,Max_Nr_Neigh,Element_Type,RX,RY,RZ,LBox,H_INDEX_START,H_INDEX_END, &
 nrnnlist,nnRx,nnRy,nnRz,nnType)

call GetZ(Z,S,HDIM) !Z = S^(-1/2);  
call MMult(ONE,Z,S,ZERO,ZI,'T','N',HDIM)

mu0 = 0.D0  ! Initial guess that might need ot be adjusted
call Fermi_S(S_Ent,HOrth,D0,QxQ,ee,Fe_vec,mu0,H0,Z,Nocc,T,OccErrLim,MaxIt,HDIM) 
call MMult(ONE,Z,D0,ZERO,X,'N','N',HDIM)  
call MMult(ONE,X,ZI,ZERO,PS,'N','N',HDIM)
call MMult(TWO,X,Z,ZERO,D,'N','T',HDIM)

PO = D0

do I = 1, Nr_atoms
  q(I) = -Znuc(I)
  do K = H_INDEX_START(I), H_INDEX_END(I)
    q(I) = q(I) + dot_product(D(:,K),S(:,K))
  enddo
enddo

SCF_ERR = 1.D0
SCF_IT = 0

open(UNIT=22,STATUS="OLD",FILE="MD.xyz")
open(UNIT=23,STATUS="OLD",FILE="Energy.dat")

It_v = ZERO 
It = ZERO
Res = q
MEM = 500
RR = ZERO
QQ = ZERO
do while (SCF_ERR > SCF_EPS) 
  SCF_IT = SCF_IT + 1

  call nearestneighborlist(nrnnlist,nndist,nnRx,nnRy,nnRz,nnType,nnStruct,nrnnStruct,RX,RY,RZ,LBox,COULCUT, &
   Nr_atoms,Max_Nr_Neigh)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,Coulomb_Pot_Real_I,Coulomb_Force_Real_I)
  do I = 1,Nr_atoms
    call Ewald_Real_Space(Coulomb_Pot_Real_I,Coulomb_Force_Real_I,I,RX,RY,RZ,LBox, &
    q,Hubbard_U,Element_Type,Nr_atoms,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,HDIM,Max_Nr_Neigh)
    Coulomb_Pot_Real(I) = Coulomb_Pot_Real_I
    Coulomb_Force_Real(:,I) = Coulomb_Force_Real_I(:)
  enddo
!$OMP END PARALLEL DO
  call Ewald_k_Space(Coulomb_Pot_k,Coulomb_Force_k,RX,RY,RZ,LBox,q,Nr_atoms,Coulomb_acc,TIMERATIO,Max_Nr_Neigh)
  coulomb_Pot = Coulomb_Pot_Real+Coulomb_Pot_k

  h1 = ZERO
  do I = 1,Nr_atoms
     do J = H_INDEX_START(I),H_INDEX_END(I)
        H1(J,J) = Hubbard_U(I)*q(I) + Coulomb_Pot(I)
     enddo
  enddo
  call MMult(ONE/TWO,S,H1,ZERO,X,'N','N',HDIM)
  call MMult(ONE/TWO,H1,S,ONE,X,'N','T',HDIM)
  H = H0 + X

  !!!!!!!!!!! DFTB+U Term Begin
  do I = 1, HDIM
  do J = 1, HDIM
     SU(I,J) = S(I,J)*DFTB_U(J)
  enddo
  enddo
  call MMult(ONE,PS,SU,ZERO,X,'T','N',HDIM)
  call MMult(ONE,SU,PS,ZERO,Y,'N','N',HDIM)
  HU = 0.25D0*(SU-X-Y)
  do I = 1,HDIM
  do J = 1,HDIM
     H_U(I,J) = HU(I,J) + HU(J,I)
  enddo
  enddo
  H = H + H_U
  !!!!!!!!!!! DFTB+U Term End

  ! Diagonalize H and construct DM with Fermi occupation and HOrth etc
  call Fermi_S(S_Ent,HOrth,DO,QxQ,ee,Fe_vec,mu0,H,Z,Nocc,T,OccErrLim,MaxIt,HDIM) 
  write(*,*) SCF_IT, ' Gap = ', ee(Nocc+1)-ee(Nocc)

  Delta_DO = DO - PO
  nDO = 0.D0
  do I = 1, HDIM
    nDO = nDO + dot_product(Delta_DO(:,I),Delta_DO(:,I))
  enddo
  nDO =sqrt(nDO)
  nDelta_DO = Delta_DO/nDO

  call MMult(ONE,Z,Delta_DO,ZERO,X,'N','N',HDIM)
  call MMult(ONE,X,ZI,ZERO,Delta_DS,'N','N',HDIM)

  do I = 1, Nr_atoms
    Res(I) = 0.D0
    do K = H_INDEX_START(I), H_INDEX_END(I)
      Res(I) = Res(I) + 2*Delta_DS(K,K)
    enddo
  enddo

  call MMult(ONE,Delta_DS,SU,ZERO,X,'T','N',HDIM)
  call MMult(ONE,SU,Delta_DS,ZERO,Y,'N','N',HDIM)
  H1 = -0.25D0*(X+Y)/nDO
  do I = 1,HDIM
  do J = 1,HDIM
     H_1(I,J) = H1(I,J) + H1(J,I)  ! Response in H from the Hubbard energy term
  enddo
  enddo

  ! Canonical density matrix response with respect to normalized DM perturbation Delta_DS/nDO
  call V_DM_U_Kernel_Fermi_No_Diag(dDO_dPO,Delta_DS/nDO,HOrth,mu0,mu1,T,RX,RY,RZ,LBox,Hubbard_U,Element_Type,Nr_atoms, &
  MaxIt,eps,m,HDIM,Max_Nr_Neigh,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,H_INDEX_START, &
  H_INDEX_END,H,H_1,S,Z,Nocc,Znuc,QxQ,ee,Fe_vec)

  dU = dDO_dPO + ((1-C_mix)/C_mix)*nDelta_DO

  ndDO_dU = 0.D0
  do I = 1, HDIM
    ndDO_dU = ndDO_dU + dot_product(nDelta_DO(:,I),dU(:,I))
  enddo

  if (SCF_IT <= 1) then
    PO = PO + C_mix*Delta_DO   ! Regular linear mixing
  else    
    PO = PO + C_mix*Delta_DO + (C_mix*C_mix*nDO/(1.D0-C_mix*ndDO_dU))*dU ! Rank-1 updated kernel mixing
  endif

  call MMult(ONE,Z,PO,ZERO,X,'N','N',HDIM)
  call MMult(ONE,X,ZI,ZERO,PS,'N','N',HDIM)
  do I = 1, Nr_atoms
    q(I) = -Znuc(I)
    do K = H_INDEX_START(I), H_INDEX_END(I)
      q(I) = q(I) + 2*PS(K,K)  ! Updated charge for the next Coulomb potential
    enddo
  enddo

  It = It + 1
  SCF_ERR = norm2(Res)/sqrt(ONE*Nr_atoms)  ! RMS error in the residual charge
enddo
  write(*,*) '- RMS SCF_ERR = ', SCF_ERR, SCF_IT

call MMult(ONE,Z,DO,ZERO,X,'N','N',HDIM)  
call MMult(ONE,X,Z,ZERO,D,'N','T',HDIM)

call nearestneighborlist(nrnnlist,nndist,nnRx,nnRy,nnRz,nnType,nnStruct,nrnnStruct,RX,RY,RZ,LBox,4.0D0, &
 Nr_atoms,Max_Nr_Neigh)

call Build_S(S,Nr_atoms,HDIM,Max_Nr_Neigh,Element_Type,RX,RY,RZ,LBox,H_INDEX_START,H_INDEX_END, &
 nrnnlist,nnRx,nnRy,nnRz,nnType)

call Get_dS(dSx,dSy,dSz,Nr_atoms,dx,HDIM,Max_Nr_Neigh,RX,RY,RZ,H_INDEX_START,H_INDEX_END, &
            nrnnlist,nnRx,nnRy,nnRz,nnType,Element_Type,LBox)
P = D  ! For full regular Pulay force and energy without pseudofunctional 
call HubbardForce(HubForce,EHub,Nr_atoms,HDIM,U,D,P,S,dSx,dSy,dSz,H_INDEX_START,H_INDEX_END)

ECoul = ZERO
do I = 1,Nr_atoms
  ECoul = ECoul + q(I)*(Hubbard_U(I)*q(I) + Coulomb_Pot(I))
enddo

D = TWO*D
X = D
do I = 1, HDIM
  X(I,I) = X(I,I) - D_atomic(I)
enddo
EBand = ZERO
EBand0 = ZERO
!H = H-H_U
do I = 1, HDIM
  EBand = EBand + dot_product(X(:,I),H(:,I))    ! Regular band energy
  EBand0 = EBand0 + dot_product(X(:,I),H0(:,I)) ! = 2Tr[DH0] Band energy only for the charge independent H0 
enddo

EEnt = -2.D0*T*S_Ent

Coulomb_Force = Coulomb_Force_Real + Coulomb_Force_k

call Get_dH(dHx,dHy,dHz,Nr_atoms,dx,HDIM,Max_Nr_Neigh,RX,RY,RZ,H_INDEX_START,H_INDEX_END, &
            nrnnlist,nnRx,nnRy,nnRz,nnType,Element_Type,LBox)

call SlaterKosterForce(SKForce,Nr_atoms,HDIM,D,dHx,dHy,dHz,H_INDEX_START,H_INDEX_END) ! from Tr[D*dH0/dR]
call PulayForce(FPUL,Nr_atoms,HDIM,Z,H,D,dSx,dSy,dSz,H_INDEX_START,H_INDEX_END) ! from 2 tr[Z*Z'*F*dS/dR]
call NonOrthCoulForce(FSCOUL,Nr_atoms,HDIM,dSx,dSy,dSz,H_INDEX_START,H_INDEX_END,D,q,Coulomb_Pot,Hubbard_U)

call LoadPairPotParameters(PairPotType1,PairPotType2,PotCoef)
call nearestneighborlist(nrnnlist,nndist,nnRx,nnRy,nnRz,nnType,nnStruct,nrnnStruct,RX,RY,RZ,LBox,1.7D0, &
                         Nr_atoms,Max_Nr_Neigh)
call PairPotential_Forces(PairForces,ERep,RX,RY,RZ,LBox,Element_Type,Nr_atoms,PairPotType1,PairPotType2,PotCoef, &
                          Max_Nr_Neigh,nrnnlist,nnRx,nnRy,nnRz,nnType)

EPOT = EBand0 + 0.5D0*ECoul  + ERep + EEnt + 2.D0*EHub
FTOT = SKForce + PairForces + FPUL + Coulomb_Force + FSCOUL + 2.D0*HubForce  ! Total force

write(*,*) ' EPOT = ', EBand0 + 0.5D0*ECoul  + ERep + EEnt + 2.D0*EHub ! This total energy expression is used
write(*,*) ' EPOTX = ', EBand - 0.5D0*ECoul + ERep + EEnt  ! Not this one, which is equivalent for EHub = 0

! Initial BC for DM, assuming orthogonal representation
PO = DO
PO_0 = DO; PO_1 = DO; PO_2 = DO; PO_3 = DO; PO_4 = DO; PO_5 = DO;

! Coefficients for modified Verlet integration
C0 = -6.D0; C1 = 14.D0; C2 = -8.D0; C3 = -3.D0; C4 = 4.D0; C5 = -1.D0; 
kappa = 1.82D0  
alpha = 0.018D0

VX = 0*RX !Initialize velocities
VY = 0*RX 
VZ = 0*RX

step = 0
do MD_step = 1,MD_Iter  !! MAIN MD LOOP

  EKIN = ZERO
  do I = 1, Nr_atoms
    EKIN = EKIN + 0.5D0*MVV2KE*Mnuc(I)*(VX(I)**2+VY(I)**2+VZ(I)**2)
  enddo
  Temperature = (TWO/THREE)*KE2T*EKIN/Nr_atoms
  Energy = EKIN+EPOT
  Time = (MD_step-1)*dt
  write(*,*) ' Time = ', Time, ' Etotal = ', Energy, 'Temperature = ', Temperature , ' RMS = ',  norm2(q-n)/sqrt(ONE*Nr_atoms), &
  ' Gap = ', ee(Nocc+1)-ee(Nocc)
  write(23,*)  Time, Energy, Temperature, norm2(q-n)/sqrt(ONE*Nr_atoms),  ee(Nocc+1)-ee(Nocc)

  if (mod(MD_step,50)== 1) then
    step = step + 1
    write(22,*) Nr_atoms
    write(22,*) 'frame',step
    do I = 1,Nr_atoms
      write(22,*) Element_Type(I), RX(I),RY(I),RZ(I) ! Coordinate output forn animation etc
    enddo
  endif

  do I = 1, Nr_atoms
    VX(I) = VX(I) + 0.5D0*dt*F2V*FTOT(1,I)/Mnuc(I)       ! First 1/2 of Leapfrog step
    VY(I) = VY(I) + 0.5D0*dt*F2V*FTOT(2,I)/Mnuc(I)       ! First 1/2 of Leapfrog step
    VZ(I) = VZ(I) + 0.5D0*dt*F2V*FTOT(3,I)/Mnuc(I)       ! First 1/2 of Leapfrog step
  enddo

  RX = RX + dt*VX                               ! Update positions
  RY = RY + dt*VY
  RZ = RZ + dt*VZ

  do I = 1,Nr_atoms ! Take care of periodic boundary conditions (maybe not necessary)
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
  enddo

  call nearestneighborlist(nrnnlist,nndist,nnRx,nnRy,nnRz,nnType,nnStruct,nrnnStruct,RX,RY,RZ,LBox,4.0D0, &
                           Nr_atoms,Max_Nr_Neigh)
  call Build_H0(H0,Hubbard_U,Nr_atoms,HDIM,Max_Nr_Neigh,Element_Type,RX,RY,RZ,LBox,H_INDEX_START,H_INDEX_END, &
                nrnnlist,nnRx,nnRy,nnRz,nnType)
  call Build_S(S,Nr_atoms,HDIM,Max_Nr_Neigh,Element_Type,RX,RY,RZ,LBox,H_INDEX_START,H_INDEX_END, &
               nrnnlist,nnRx,nnRy,nnRz,nnType)

  if (Full_Kernel_MD == 0) then
    if ((MD_step <= 1000)) then  ! Standard scaled delta function approach
      ! Integrate extended electronic degrees of freedom
      ! using the modified Verlet scheme with some weak dissipation
      PO = 2.D0*PO_0 - PO_1 + C_mix*kappa*(DO-PO) + alpha*(C0*PO_0+C1*PO_1+C2*PO_2+C3*PO_3+C4*PO_4+C5*PO_5)
      PO_5 = PO_4; PO_4 = PO_3; PO_3 = PO_2; PO_2 = PO_1; PO_1 = PO_0; PO_0 = PO;

      call GetZ(Z,S,HDIM) !Z = S^(-1/2)
      call MMult(ONE,Z,S,ZERO,ZI,'T','N',HDIM) !ZI = S^(1/2)
      call MMult(ONE,Z,PO,ZERO,X,'N','N',HDIM)
      call MMult(ONE,X,ZI,ZERO,PS,'N','N',HDIM)
      call MMult(ONE,X,Z,ZERO,P,'N','T',HDIM)

      do I = 1, Nr_atoms
        qx(I) = -Znuc(I)
        do K = H_INDEX_START(I), H_INDEX_END(I)
          qx(I) = qx(I) + 2.D0*PS(K,K)
        enddo
      enddo
      n = qx
    else ! New approach to include rank-1 updated kernel response for the electronic equation of motion
      Delta_DO = DO - PO
      nDO = 0.D0
      do I = 1, HDIM
        nDO = nDO + dot_product(Delta_DO(:,I),Delta_DO(:,I))
      enddo
      nDO =sqrt(nDO)
      nDelta_DO = Delta_DO/nDO

      call MMult(ONE,Z,Delta_DO,ZERO,X,'N','N',HDIM)
      call MMult(ONE,X,ZI,ZERO,Delta_DS,'N','N',HDIM)

      do I = 1, Nr_atoms
        Res(I) = 0.D0
        do K = H_INDEX_START(I), H_INDEX_END(I)
          Res(I) = Res(I) + 2.D0*Delta_DS(K,K)
        enddo
      enddo

      call MMult(ONE,Delta_DS,SU,ZERO,X,'T','N',HDIM)
      call MMult(ONE,SU,Delta_DS,ZERO,Y,'N','N',HDIM)
      H1 = -0.25D0*(X+Y)/nDO
      do I = 1,HDIM
      do J = 1,HDIM
         H_1(I,J) = H1(I,J) + H1(J,I) ! Hubbard contribution
      enddo
      enddo

      ! Canonical density matrix response dDO_dPO with respect to normalized DM perturbation Delta_DS/nDO
      ! No diagonalization is needed since we can reuse eivenvectors from Fermi_S
      call V_DM_U_Kernel_Fermi_No_Diag(dDO_dPO,Delta_DS/nDO,HOrth,mu0,mu1,T,RX,RY,RZ,LBox,Hubbard_U,Element_Type,Nr_atoms, &
      MaxIt,eps,m,HDIM,Max_Nr_Neigh,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,H_INDEX_START, &
      H_INDEX_END,H,H_1,S,Z,Nocc,Znuc,QxQ,ee,Fe_vec)

      dU = dDO_dPO + ((1-C_mix)/C_mix)*nDelta_DO

      ndDO_dU = 0.D0
      do I = 1, HDIM
        ndDO_dU = ndDO_dU + dot_product(nDelta_DO(I,:),dU(:,I))
      enddo
  
      d2PO = C_mix*Delta_DO + (C_mix*C_mix*nDO/(1.D0-C_mix*ndDO_dU))*dU

      ! Integrate extended electronic degrees of freedom using rank-1 updated equation of motion
      ! using the modified Verlet scheme with some weak dissipation
      PO = 2.D0*PO_0 - PO_1 + kappa*d2PO + alpha*(C0*PO_0+C1*PO_1+C2*PO_2+C3*PO_3+C4*PO_4+C5*PO_5)
      PO_5 = PO_4; PO_4 = PO_3; PO_3 = PO_2; PO_2 = PO_1; PO_1 = PO_0; PO_0 = PO;

      call GetZ(Z,S,HDIM) !Z = S^(-1/2);
      call MMult(ONE,Z,S,ZERO,ZI,'T','N',HDIM)

      call MMult(ONE,Z,PO,ZERO,X,'N','N',HDIM)
      call MMult(ONE,X,ZI,ZERO,PS,'N','N',HDIM)
      call MMult(ONE,X,Z,ZERO,P,'N','T',HDIM)
      do I = 1, Nr_atoms
        qx(I) = -Znuc(I)
        do K = H_INDEX_START(I), H_INDEX_END(I)
          qx(I) = qx(I) + 2*PS(K,K)
        enddo
      enddo
      n = qx
    endif
  endif

  call nearestneighborlist(nrnnlist,nndist,nnRx,nnRy,nnRz,nnType,nnStruct,nrnnStruct,RX,RY,RZ,LBox,COULCUT, &
                           Nr_atoms,Max_Nr_Neigh)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,Coulomb_Pot_Real_I,Coulomb_Force_Real_I)
  do I = 1,Nr_atoms
    call Ewald_Real_Space(Coulomb_Pot_Real_I,Coulomb_Force_Real_I,I,RX,RY,RZ,LBox, &
    n,Hubbard_U,Element_Type,Nr_atoms,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,HDIM,Max_Nr_Neigh)
    Coulomb_Pot_Real(I) = Coulomb_Pot_Real_I
    Coulomb_Force_Real(:,I) = Coulomb_Force_Real_I(:)
  enddo
!$OMP END PARALLEL DO
  call Ewald_k_Space(Coulomb_Pot_k,Coulomb_Force_k,RX,RY,RZ,LBox,n,Nr_atoms,Coulomb_acc,TIMERATIO,Max_Nr_Neigh)
  Coulomb_Pot = Coulomb_Pot_Real+Coulomb_Pot_k
  Coulomb_Force = Coulomb_Force_Real + Coulomb_Force_k

  H1 = ZERO
  do I = 1,Nr_atoms
     do J = H_INDEX_START(I),H_INDEX_END(I)
         H1(J,J) = Hubbard_U(I)*n(I) + Coulomb_Pot(I)
     enddo
  enddo
  call MMult(ONE/TWO,S,H1,ZERO,X,'N','N',HDIM) ! Can be done faster since H1 is diagonal
  call MMult(ONE/TWO,H1,S,ONE,X,'N','T',HDIM)  ! and this is just the transpose
  H = H0 + X

!!!!!!!!!!! DFTB+U Term
  do I = 1, HDIM
  do J = 1, HDIM
     SU(I,J) = S(I,J)*DFTB_U(J)
  enddo
  enddo
  call MMult(ONE,PS,SU,ZERO,X,'T','N',HDIM)
  call MMult(ONE,SU,PS,ZERO,Y,'N','N',HDIM)
  HU = 0.25D0*(SU-X-Y)
  do I = 1,HDIM
  do J = 1,HDIM
     H_U(I,J) = HU(I,J) + HU(J,I)
  enddo
  enddo
  H = H + H_U
!!!!!!!!!!! DFTB+U Term

  call Fermi_S(S_Ent,HOrth,DO,QxQ,ee,Fe_vec,mu0,H,Z,Nocc,T,OccErrLim,MaxIt,HDIM)
  call MMult(ONE,Z,DO,ZERO,X,'N','N',HDIM)
  call MMult(ONE,X,Z,ZERO,D,'N','T',HDIM)
  call MMult(ONE,Z,DO,ZERO,X,'N','N',HDIM)
  call MMult(ONE,X,ZI,ZERO,DS,'N','N',HDIM)

! Net number of Mulliken electrons -Znuc
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K)
  do I = 1, Nr_atoms
    q(I) = -Znuc(I)
    do K = H_INDEX_START(I), H_INDEX_END(I)
      q(I) = q(I) + 2.D0*dot_product(D(:,K),S(:,K))
    enddo
  enddo
!$OMP END PARALLEL DO

  ECoul = ZERO
  ECoul0 = ZERO
  do I = 1,Nr_atoms
    ECoul = ECoul + n(I)*(Hubbard_U(I)*n(I) + Coulomb_Pot(I))
    ECoul0 = ECoul0 + (2.D0*q(I) - n(I))*(Hubbard_U(I)*n(I) + Coulomb_Pot(I))
  enddo

  call nearestneighborlist(nrnnlist,nndist,nnRx,nnRy,nnRz,nnType,nnStruct,nrnnStruct,RX,RY,RZ,LBox,1.7D0, &
                           Nr_atoms,Max_Nr_Neigh)
  call PairPotential_Forces(PairForces,ERep,RX,RY,RZ,LBox,Element_Type,Nr_atoms,PairPotType1,PairPotType2,PotCoef, &
                            Max_Nr_Neigh,nrnnlist,nnRx,nnRy,nnRz,nnType)

  D = TWO*D
  X = D
  do I = 1, HDIM
    X(I,I) = X(I,I) - D_atomic(I)
  enddo
  EBand = ZERO
  EBand0 = ZERO
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I) &
!$OMP REDUCTION(+:EBand0,EBand)
  do I = 1, HDIM
    EBand = EBand + dot_product(X(:,I),H(:,I))
    EBand0 = EBand0 + dot_product(X(:,I),H0(:,I))
  enddo
!$OMP END PARALLEL DO

  call nearestneighborlist(nrnnlist,nndist,nnRx,nnRy,nnRz,nnType,nnStruct,nrnnStruct,RX,RY,RZ,LBox,4.0D0, &
                           Nr_atoms,Max_Nr_Neigh)

  call Get_dS(dSx,dSy,dSz,Nr_atoms,dx,HDIM,Max_Nr_Neigh,RX,RY,RZ,H_INDEX_START,H_INDEX_END, &
              nrnnlist,nnRx,nnRy,nnRz,nnType,Element_Type,LBox)
  call Get_dH(dHx,dHy,dHz,Nr_atoms,dx,HDIM,Max_Nr_Neigh,RX,RY,RZ,H_INDEX_START,H_INDEX_END, &
              nrnnlist,nnRx,nnRy,nnRz,nnType,Element_Type,LBox)

  call SlaterKosterForce(SKForce,Nr_atoms,HDIM,D,dHx,dHy,dHz,H_INDEX_START,H_INDEX_END) ! from Tr[D*dH0/dR]
  call PulayForce(FPUL,Nr_atoms,HDIM,Z,H,D,dSx,dSy,dSz,H_INDEX_START,H_INDEX_END) ! from 2 tr[Z*Z'*F*dS/dR]
  call NonOrthCoulForce(FSCOUL,Nr_atoms,HDIM,dSx,dSy,dSz,H_INDEX_START,H_INDEX_END,D,n,Coulomb_Pot,Hubbard_U)
!!!!!!!!!!!!!!!!!!! TEST HUB FORCE !!!!!!!!!!!!
  call HubbardForce(HubForce,EHub,Nr_atoms,HDIM,U,0.5D0*D,P,S,dSx,dSy,dSz,H_INDEX_START,H_INDEX_END)
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  write(*,*) ' Hubbarde E and F = ', EHub, HubForce(1,1)

  EEnt = -2*T*S_Ent
  EPOT = EBand0 + 0.5D0*ECoul0  + ERep + EEnt + 2.D0*EHub
!  EPOT = EBand - 0.5D0*ECoul  + ERep + EEnt + 2.D0*EHub

  do I = 1, Nr_atoms
    Coulomb_Force(1,I) = (2.D0*q(I)-n(I))*Coulomb_Force(1,I)/n(I)  ! Adjust forces of linearized XLBOMD functional
    Coulomb_Force(2,I) = (2.D0*q(I)-n(I))*Coulomb_Force(2,I)/n(I) 
    Coulomb_Force(3,I) = (2.D0*q(I)-n(I))*Coulomb_Force(3,I)/n(I) 
  enddo

  FTOT = SKForce + PairForces + FPUL + Coulomb_Force + FSCOUL + 2.D0*HubForce  ! Total force

!write(*,*) ' EPOT = ', EPOT
!write(*,*) '-SKForce + PairForces + FPUL + Coulomb_Force + FSCOUL + HubForce  = ', &
!            SKForce(1,1),PairForces(1,1),FPUL(1,1),Coulomb_Force(1,1),FSCOUL(1,1),HubForce(1,1)

  do I = 1, Nr_atoms
    VX(I) = VX(I) + 0.5D0*dt*F2V*FTOT(1,I)/Mnuc(I)       ! second 1/2 of Leapfrog step
    VY(I) = VY(I) + 0.5D0*dt*F2V*FTOT(2,I)/Mnuc(I)       ! second 1/2 of Leapfrog step
    VZ(I) = VZ(I) + 0.5D0*dt*F2V*FTOT(3,I)/Mnuc(I)       ! second 1/2 of Leapfrog step
  enddo

enddo

close(22)
close(23)

end program XLBOMD
