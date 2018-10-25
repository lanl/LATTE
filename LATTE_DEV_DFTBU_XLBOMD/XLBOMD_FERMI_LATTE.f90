!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    LATTE-XLBOMD Scaled down developers version of LATTE    !!!
!!!    Anders M.N. Niklasson et al. amn@lanl.gov               !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program XLBOMD

use omp_lib
implicit none
integer, parameter       :: PREC = 8

!integer, parameter       :: Nr_atoms = 30, HDIM = 60, Nocc = 40, Max_Nr_Neigh = 500   ! 10 water molecules
!real(PREC), parameter	  :: Lx = 5.0000000D0, Ly = 5.0000000D0, Lz = 5.0000000D0
!real(PREC), parameter    :: T = 1000.00D0, dt = 0.50D0, C_mix = 0.50D0
!integer, parameter       :: Full_Kernel_Initial = 0, Full_Kernel_MD = 0

!integer, parameter       :: Nr_atoms = 24, HDIM = 48, Nocc = 32, Max_Nr_Neigh = 500   ! 8 water molecules
!real(PREC), parameter	  :: Lx = 8.2670000D0, Ly = 8.2670000D0, Lz = 8.2670000D0
!real(PREC), parameter    :: T = 300.00D0, dt = 0.4D0, C_mix = 0.15D0
!integer, parameter       :: Full_Kernel_Initial = 1, Full_Kernel_MD = 0

!integer, parameter       :: Nr_atoms = 300, HDIM = 600, Nocc = 400, Max_Nr_Neigh = 500   ! 100 water molecules
!real(PREC), parameter	  :: Lx = 14.404799, Ly = 14.404799, Lz = 14.404799

!integer, parameter       :: Nr_atoms = 55, HDIM = 220, Nocc = 110, Max_Nr_Neigh = 500   ! 55 Amorph C
!real(PREC), parameter    :: Lx = 10.0000000, Ly = 10.0000000, Lz = 10.0000000
!real(PREC), parameter    :: T = 500.00D0, dt = 1.00D0, C_mix = 0.05D0
!integer, parameter       :: Full_Kernel_Initial = 0, Full_Kernel_MD = 0

!integer, parameter       :: Nr_atoms = 130, HDIM = 520, Nocc = 268, Max_Nr_Neigh = 500   ! non-PGM
!real(PREC), parameter    :: Lx = 20.0000, Ly = 17.0480000, Lz = 9.6860000
!real(PREC), parameter    :: T = 5000.00D0, dt = 1.00D0, C_mix = 0.05D0
!integer, parameter       :: Full_Kernel_Initial = 0, Full_Kernel_MD = 2

!integer, parameter       :: Nr_atoms = 7, HDIM = 19, Nocc = 12, Max_Nr_Neigh = 500   ! 7 NM
!real(PREC), parameter    :: Lx = 12.1151720, Ly = 12.1151720, Lz = 12.1151720
!real(PREC), parameter    :: T = 100.00D0, dt = 0.00D0, C_mix = 0.15D0
!integer, parameter       :: Full_Kernel_Initial = 0, Full_Kernel_MD = 0

integer, parameter       :: Nr_atoms = 49, HDIM = 133, Nocc = 84, Max_Nr_Neigh = 500   ! 49 NM
real(PREC), parameter    :: Lx = 12.1151720, Ly = 12.1151720, Lz = 12.1151720
real(PREC), parameter    :: T = 100.00D0, dt = 0.25D0, C_mix = 0.10D0
integer, parameter       :: Full_Kernel_Initial = 0, Full_Kernel_MD = 0

!integer, parameter       :: Nr_atoms = 7, HDIM = 19, Nocc = 12, Max_Nr_Neigh = 500   ! 1 NM
!real(PREC), parameter    :: Lx = 12.1151720, Ly = 12.1151720, Lz = 12.1151720
!real(PREC), parameter    :: T = 10000.00D0, dt = 0.25D0, C_mix = 0.10D0
!integer, parameter       :: Full_Kernel_Initial = 0, Full_Kernel_MD = 0

!integer, parameter       :: Nr_atoms = 70, HDIM = 190, Nocc = 120, Max_Nr_Neigh = 500   ! 70 NM
!real(PREC), parameter    :: Lx = 9.6158180, Ly = 9.6158180, Lz = 9.6158180
!real(PREC), parameter    :: T = 1000.00D0, dt = 0.4D0, C_mix = 0.10D0
!integer, parameter       :: Full_Kernel_Initial = 0, Full_Kernel_MD = 0

!integer, parameter       :: Nr_atoms = 99, HDIM = 396, Nocc = 198, Max_Nr_Neigh = 500   ! 99 Graphene
!real(PREC), parameter    :: Lx = 14.216000, Ly = 21.274000, Lz = 12.361400
!real(PREC), parameter    :: T = 3000.00D0, dt = 1.0D0, C_mix = 0.05D0
!integer, parameter       :: Full_Kernel_Initial = 0, Full_Kernel_MD = 0

!integer, parameter       :: Nr_atoms = 240, HDIM = 600, Nocc = 300, Max_Nr_Neigh = 500   ! Benzene
!real(PREC), parameter    :: Lx = 11.078000, Ly = 11.078000, Lz = 11.078000

!!!!!!!!!!!!!!!!!!!!!!!!    
integer, parameter       :: m = 18, MD_Iter = 200, MaxIt = 20 
real(PREC), parameter    :: Coulomb_acc = 1e-7, TIMERATIO = 10.D0, ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0, SIX = 6.D0
real(PREC), parameter    :: THREE = 3.D0

real(PREC), parameter    :: eps = 1e-7, SCF_EPS = 1e-9, dx = 0.0001D0  
real(PREC), parameter    :: OccErrLim = 1e-9
!real(PREC), parameter    :: eps = 1e-7, SCF_EPS = 1e-3, dx = 0.0001D0  
!real(PREC), parameter    :: OccErrLim = 1e-3

real(PREC), parameter    :: SEVEN = 7.D0
real(PREC), parameter    :: F2V = 0.01602176487D0/1.660548782D0
real(PREC), parameter    :: MVV2KE = 166.0538782D0/1.602176487D0
real(PREC), parameter    :: KE2T = 1.D0/0.000086173435D0
real(PREC)               :: RX(Nr_atoms), RY(Nr_atoms), RZ(Nr_atoms), q(Nr_atoms), LBox(3)
integer                  :: H_INDEX_START(Nr_atoms), H_INDEX_END(Nr_atoms)
real(PREC)               :: Znuc(Nr_atoms), Mnuc(Nr_atoms), Hubbard_U(Nr_atoms), COULCUT
real(PREC)               :: D_atomic(HDIM), H0(HDIM,HDIM), H(HDIM,HDIM), H1(HDIM,HDIM)
real(PREC)               :: D0(HDIM,HDIM), D_0(HDIM,HDIM), D(HDIM,HDIM), EYE(HDIM,HDIM), II(Nr_atoms,Nr_atoms)
real(PREC)               :: S(HDIM,HDIM), SI(HDIM,HDIM), Z(HDIM,HDIM), X(HDIM,HDIM), Y(HDIM,HDIM)
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
real(PREC)               :: PotCoef(16,10), PairForces(3,Nr_atoms), ERep, ECoul, EBand, EEnt, EPOT
real(PREC)               :: dSx(HDIM,HDIM), dSy(HDIM,HDIM), dSz(HDIM,HDIM)
real(PREC)               :: dHx(HDIM,HDIM), dHy(HDIM,HDIM), dHz(HDIM,HDIM)
real(PREC)               :: SKForce(3,Nr_atoms), FPUL(3,Nr_atoms), FSCOUL(3,Nr_atoms), FTOT(3,Nr_atoms)
real(PREC)               :: n(Nr_atoms),n_0(Nr_atoms),n_1(Nr_atoms),n_2(Nr_atoms),n_3(Nr_atoms),n_4(Nr_atoms)
real(PREC)               :: n_5(Nr_atoms), n_6(Nr_atoms), mu_x, mu_0, mu_1, mu_2, mu_3, mu_4, mu_5, mu_6
real(PREC)               :: C0, C1, C2, C3, C4, C5, C6, kappa, cc, Time
real(PREC)               :: VX(Nr_atoms), VY(Nr_atoms), VZ(Nr_atoms), EKIN, Energy
real(PREC)               :: dn2dt2(Nr_atoms), Temperature, W(Nr_atoms,Nr_atoms), v(Nr_atoms) 
real(PREC)               :: u(Nr_atoms)
real(PREC)               :: JJV(Nr_atoms,Nr_atoms), dq_dv(Nr_atoms), Xtmp(Nr_atoms,Nr_atoms), X_Xtmp(Nr_atoms,Nr_atoms)
real(PREC)               :: Res(Nr_atoms), Res_old(Nr_atoms), QQ(Nr_atoms,Nr_atoms), RR(Nr_atoms,Nr_atoms), gq(Nr_atoms)
real(PREC)               :: QQ_tmp(Nr_atoms,Nr_atoms), RR_tmp(Nr_atoms,Nr_atoms)
real(PREC)               :: QxQ(HDIM,HDIM), ee(HDIM), Fe_vec(HDIM), Err
integer                  :: I,J,K,L, LL,SCF_IT, MD_step, It_v, It, MEM, MM
real(PREC)               :: norm_v, aa, bb

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

! Get the atomic diagonal density matrix
call AtomicDensityMatrix(D_atomic,Nr_atoms,H_INDEX_START,H_INDEX_END,HDIM,Znuc) 

call nearestneighborlist(nrnnlist,nndist,nnRx,nnRy,nnRz,nnType,nnStruct,nrnnStruct,RX,RY,RZ,LBox,4.0D0, &
 Nr_atoms,Max_Nr_Neigh)
call Build_H0(H0,Hubbard_U,Nr_atoms,HDIM,Max_Nr_Neigh,Element_Type,RX,RY,RZ,LBox,H_INDEX_START,H_INDEX_END, &
 nrnnlist,nnRx,nnRy,nnRz,nnType)
call Build_S(S,Nr_atoms,HDIM,Max_Nr_Neigh,Element_Type,RX,RY,RZ,LBox,H_INDEX_START,H_INDEX_END, &
 nrnnlist,nnRx,nnRy,nnRz,nnType)

call GetZ(Z,S,HDIM) !Z = S^(-1/2);  

mu0 = 0.D0  ! Initial guess that might need ot be adjusted
call Fermi_S(S_Ent,HOrth,D0,QxQ,ee,Fe_vec,mu0,H0,Z,Nocc,T,OccErrLim,MaxIt,HDIM) 

X = ZERO
!call RecONFermi(D0,X,HOrth,T,mu0,Nocc,m,eps,MaxIt,HDIM)
call MMult(TWO,Z,D0,ZERO,X,'N','N',HDIM)  
call MMult(ONE,X,Z,ZERO,D,'N','T',HDIM)

do I = 1, Nr_atoms
  q(I) = -Znuc(I)
  do K = H_INDEX_START(I), H_INDEX_END(I)
    q(I) = q(I) + dot_product(D(:,K),S(:,K))
  enddo
enddo
write(*,*) ' q0 = ', q(1:5)

call nearestneighborlist(nrnnlist,nndist,nnRx,nnRy,nnRz,nnType,nnStruct,nrnnStruct,RX,RY,RZ,LBox,COULCUT, &
 Nr_atoms,Max_Nr_Neigh)

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
  h = H0 + X

  call Fermi_S(S_Ent,HOrth,D0,QxQ,ee,Fe_vec,mu0,H,Z,Nocc,T,OccErrLim,MaxIt,HDIM)

  call MMult(TWO,Z,D0,ZERO,X,'N','N',HDIM)
  call MMult(ONE,X,Z,ZERO,D,'N','T',HDIM)

! Net number of Mulliken electrons -Znuc
  do I = 1, Nr_atoms
    gq(I) = -Znuc(I)
    do K = H_INDEX_START(I), H_INDEX_END(I)
      gq(I) = gq(I) + dot_product(D(:,K),S(:,K))
    enddo
  enddo

! UPDATE KERNEL FOR SCF OPTIMIZATION USING QUANTUM PERTURBATION THEORY
  Res = gq-q  ! Residual Res = q[n]-n

  if (Full_Kernel_Initial == 0) then
    if (SCF_IT <= 1) then
      mu1 = ZERO
      q = q + C_mix*Res    !! Linear mixing
    else
      norm_v = norm2(Res)
      v = Res/norm_v
      call V_Kernel_Fermi_No_Diag(D0,dq_dv,v,HOrth,mu0,mu1,T,RX,RY,RZ,LBox,Hubbard_U,Element_Type,Nr_atoms,MaxIt,eps,m, &
      HDIM,Max_Nr_Neigh,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,H_INDEX_START,H_INDEX_END,H,S,Z,Nocc, &
      Znuc,QxQ,ee,Fe_vec)
      u = dq_dv + ((1.D0-C_mix)/C_mix)*v
      q = q + C_mix*Res + C_mix**2*norm_v/(ONE-C_mix*dot_product(v,u))*u
!      write(*,*) ' 2nd q-Mix Coeff = ', C_mix**2*norm_v/(ONE-C_mix*dot_product(v,u))
    endif
  elseif (Full_Kernel_Initial == 1) then
    if (SCF_IT <= 2) then
      mu1 = ZERO
      q = q + C_mix*Res  !! Linear Mixing, do a few more initially if necessary
    else
      call Kernel_Fermi_NoDiag(KK,JJ,D0,mu0,mu1,T,RX,RY,RZ,LBox,Hubbard_U,Element_Type,Nr_atoms, &
      MaxIt,eps, m,HDIM, Max_Nr_Neigh,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,H_INDEX_START, &
      H_INDEX_END,H,S,Z,Nocc,Znuc,QxQ,ee,Fe_vec)
      q = q - MATMUL(KK,Res)
    endif
  else
     q = q + C_mix*Res  !! Linear mixing
  endif

  It = It + 1
  SCF_ERR = norm2(Res)/sqrt(ONE*Nr_atoms)
  write(*,*) ' RMS SCF_ERR = ', SCF_ERR, SCF_IT
enddo

ECoul = ZERO
do I = 1,Nr_atoms
  ECoul = ECoul + q(I)*(Hubbard_U(I)*q(I) + Coulomb_Pot(I))
enddo

X = D
do I = 1, HDIM
  X(I,I) = X(I,I) - D_atomic(I)
enddo
EBand = ZERO
do I = 1, HDIM
  EBand = EBand + dot_product(X(:,I),H(:,I))
enddo
write(*,*) ' EBand = ', EBand



EEnt = -2.D0*T*S_Ent

Coulomb_Force = Coulomb_Force_Real + Coulomb_Force_k

call nearestneighborlist(nrnnlist,nndist,nnRx,nnRy,nnRz,nnType,nnStruct,nrnnStruct,RX,RY,RZ,LBox,4.0D0, &
                         Nr_atoms,Max_Nr_Neigh)

call Get_dS(dSx,dSy,dSz,Nr_atoms,dx,HDIM,Max_Nr_Neigh,RX,RY,RZ,H_INDEX_START,H_INDEX_END, &
            nrnnlist,nnRx,nnRy,nnRz,nnType,Element_Type,LBox)
call Get_dH(dHx,dHy,dHz,Nr_atoms,dx,HDIM,Max_Nr_Neigh,RX,RY,RZ,H_INDEX_START,H_INDEX_END, &
            nrnnlist,nnRx,nnRy,nnRz,nnType,Element_Type,LBox)
!
call SlaterKosterForce(SKForce,Nr_atoms,HDIM,D,dHx,dHy,dHz,H_INDEX_START,H_INDEX_END) ! from Tr[D*dH0/dR]
call PulayForce(FPUL,Nr_atoms,HDIM,Z,H,D,dSx,dSy,dSz,H_INDEX_START,H_INDEX_END) ! from 2 tr[Z*Z'*F*dS/dR]
call NonOrthCoulForce(FSCOUL,Nr_atoms,HDIM,dSx,dSy,dSz,H_INDEX_START,H_INDEX_END,D,q,Coulomb_Pot,Hubbard_U)

call LoadPairPotParameters(PairPotType1,PairPotType2,PotCoef)
call nearestneighborlist(nrnnlist,nndist,nnRx,nnRy,nnRz,nnType,nnStruct,nrnnStruct,RX,RY,RZ,LBox,1.7D0, &
                         Nr_atoms,Max_Nr_Neigh)
call PairPotential_Forces(PairForces,ERep,RX,RY,RZ,LBox,Element_Type,Nr_atoms,PairPotType1,PairPotType2,PotCoef, &
                          Max_Nr_Neigh,nrnnlist,nnRx,nnRy,nnRz,nnType)

EPOT = EBand - 0.5D0*ECoul  + ERep + EEnt
write(*,*) ' EBand, ECoul, ERep, EEnt = ', EBand, ECoul, ERep, EEnt
write(*,*) ' EPOT = ', EPOT

FTOT = SKForce + PairForces + FPUL + Coulomb_Force + FSCOUL  ! Total force

write(*,*) ' SKForce + PairForces + FPUL + Coulomb_Force + FSCOUL = ', &
  SKForce(1,1), PairForces(1,1), FPUL(1,1), Coulomb_Force(1,1), FSCOUL(1,1)



! Initial BC for n
n = q; n_0 = q; n_1 = q; n_2 = q; n_3 = q; n_4 = q; n_5 = q; n_6 = q;

! Coefficients for modified Verlet integration
C0 = -14; C1 = 36; C2 = -27; C3 = -2; C4 = 12;  C5 = -6; C6 = 1;
kappa = 1.84D0
alpha = 0.00055D0

!C0 = -6.D0; C1 = 14.D0; C2 = -8.D0; C3 = -3.D0; C4 = 4.D0; C5 = -1.D0;  C6 = 0.D0;
!kappa = 1.82D0  
!alpha = 0.018D0

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
  write(23,*)  Time, Energy, Temperature, norm2(q-n)/sqrt(ONE*Nr_atoms), mu0 , ee(Nocc+1)-ee(Nocc)
  if (mod(MD_step,50)== 1) then
    step = step + 1
    write(22,*) Nr_atoms
    write(22,*) 'frame',step
    do I = 1,Nr_atoms
      write(22,*) Element_Type(I), RX(I),RY(I),RZ(I)
    enddo
  endif

  do I = 1, Nr_atoms
    VX(I) = VX(I) + 0.5D0*dt*F2V*FTOT(1,I)/Mnuc(I)       ! First 1/2 of Leapfrog step
    VY(I) = VY(I) + 0.5D0*dt*F2V*FTOT(2,I)/Mnuc(I)       ! First 1/2 of Leapfrog step
    VZ(I) = VZ(I) + 0.5D0*dt*F2V*FTOT(3,I)/Mnuc(I)       ! First 1/2 of Leapfrog step
  enddo

  ! UPDATE KERNEL FOR DIRECTION q-n 
  call nearestneighborlist(nrnnlist,nndist,nnRx,nnRy,nnRz,nnType,nnStruct,nrnnStruct,RX,RY,RZ,LBox,COULCUT, &
  Nr_atoms,Max_Nr_Neigh)

  if (Full_Kernel_MD == 0) then
    if (MD_step == 1) then
      dn2dt2 = ZERO*(q-n)
    else
      Res = q-n
      norm_v = norm2(Res)
      v = Res/norm_v
      call V_Kernel_Fermi_No_Diag(D0,dq_dv,v,HOrth,mu0,mu1,T,RX,RY,RZ,LBox,Hubbard_U,Element_Type,Nr_atoms,MaxIt,eps,m, &
      HDIM,Max_Nr_Neigh,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,H_INDEX_START,H_INDEX_END,H,S,Z,Nocc, &
      Znuc,QxQ,ee,Fe_vec)
      u = dq_dv + v*(1.D0-C_mix)/C_mix
      dn2dt2 =  C_mix*Res + (C_mix**2)*norm_v/(ONE-C_mix*dot_product(v,u))*u
    endif
  elseif (Full_Kernel_MD == 1) then
    Res = q-n
    if (mod(MD_step,500) == 1) then
     call Kernel_Fermi_NoDiag(KK,JJ,D0,mu0,mu1,T,RX,RY,RZ,LBox,Hubbard_U,Element_Type,Nr_atoms, &
     MaxIt,eps, m,HDIM, Max_Nr_Neigh,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,H_INDEX_START, &
     H_INDEX_END,H,S,Z,Nocc,Znuc,QxQ,ee,Fe_vec)
    endif
    dn2dt2 = -MATMUL(KK,(q-n))
  elseif (Full_Kernel_MD == 2) then
    Res = q-n
    if (mod(MD_step,Nr_atoms*1) == 1) then
      call Kernel_Fermi_NoDiag(KK,JJ,D0,mu0,mu1,T,RX,RY,RZ,LBox,Hubbard_U,Element_Type,Nr_atoms, &
      MaxIt,eps, m,HDIM, Max_Nr_Neigh,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,H_INDEX_START, &
      H_INDEX_END,H,S,Z,Nocc,Znuc,QxQ,ee,Fe_vec)
    endif
    dn2dt2 = -MATMUL(KK,(q-n))
  else
    Res = q-n
    dn2dt2 = C_mix*Res               !  Linear mixing as an alternative
  endif

  RX = RX + dt*VX                               ! Update positions
  RY = RY + dt*VY 
  RZ = RZ + dt*VZ 

  do I = 1,Nr_atoms
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

  call GetZ(Z,S,HDIM) !Z = S^(-1/2);

! Integrate equation of motion for charge density, d^2n/dt^2 = w^2(q-n)
  n = 2*n_0 - n_1 + kappa*dn2dt2 + alpha*(C0*n_0+C1*n_1+C2*n_2+C3*n_3+C4*n_4+C5*n_5+C6*n_6)
  n_6 = n_5; n_5 = n_4; n_4 = n_3; n_3 = n_2; n_2 = n_1; n_1 = n_0; n_0 = n

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
  call MMult(ONE/TWO,S,H1,ZERO,X,'N','N',HDIM)
  call MMult(ONE/TWO,H1,S,ONE,X,'N','T',HDIM)
  H = H0 + X

  call Fermi_S(S_Ent,HOrth,D0,QxQ,ee,Fe_vec,mu0,H,Z,Nocc,T,OccErrLim,MaxIt,HDIM)
  call MMult(TWO,Z,D0,ZERO,X,'N','N',HDIM)
  call MMult(ONE,X,Z,ZERO,D,'N','T',HDIM)

! Net number of Mulliken electrons -Znuc
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,K)
  do I = 1, Nr_atoms
    q(I) = -Znuc(I)
    do K = H_INDEX_START(I), H_INDEX_END(I)
      q(I) = q(I) + dot_product(D(:,K),S(:,K))
    enddo
  enddo
!$OMP END PARALLEL DO

  ECoul = ZERO
  do I = 1,Nr_atoms
    ECoul = ECoul + n(I)*(Hubbard_U(I)*n(I) + Coulomb_Pot(I))
  enddo

  call nearestneighborlist(nrnnlist,nndist,nnRx,nnRy,nnRz,nnType,nnStruct,nrnnStruct,RX,RY,RZ,LBox,1.7D0, &
                           Nr_atoms,Max_Nr_Neigh)
  call PairPotential_Forces(PairForces,ERep,RX,RY,RZ,LBox,Element_Type,Nr_atoms,PairPotType1,PairPotType2,PotCoef, &
                            Max_Nr_Neigh,nrnnlist,nnRx,nnRy,nnRz,nnType)

  X = D
  do I = 1, HDIM
    X(I,I) = X(I,I) - D_atomic(I)
  enddo
  EBand = ZERO
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I) &
!$OMP REDUCTION(+:EBand)
  do I = 1, HDIM
    EBand = EBand + dot_product(X(:,I),H(:,I))
  enddo
!$OMP END PARALLEL DO

  EEnt = -2*T*S_Ent
  EPOT = EBand - 0.5D0*ECoul  + ERep + EEnt

  call nearestneighborlist(nrnnlist,nndist,nnRx,nnRy,nnRz,nnType,nnStruct,nrnnStruct,RX,RY,RZ,LBox,4.0D0, &
                           Nr_atoms,Max_Nr_Neigh)

  call Get_dS(dSx,dSy,dSz,Nr_atoms,dx,HDIM,Max_Nr_Neigh,RX,RY,RZ,H_INDEX_START,H_INDEX_END, &
              nrnnlist,nnRx,nnRy,nnRz,nnType,Element_Type,LBox)
  call Get_dH(dHx,dHy,dHz,Nr_atoms,dx,HDIM,Max_Nr_Neigh,RX,RY,RZ,H_INDEX_START,H_INDEX_END, &
              nrnnlist,nnRx,nnRy,nnRz,nnType,Element_Type,LBox)

  call SlaterKosterForce(SKForce,Nr_atoms,HDIM,D,dHx,dHy,dHz,H_INDEX_START,H_INDEX_END) ! from Tr[D*dH0/dR]
  call PulayForce(FPUL,Nr_atoms,HDIM,Z,H,D,dSx,dSy,dSz,H_INDEX_START,H_INDEX_END) ! from 2 tr[Z*Z'*F*dS/dR]
  call NonOrthCoulForce(FSCOUL,Nr_atoms,HDIM,dSx,dSy,dSz,H_INDEX_START,H_INDEX_END,D,n,Coulomb_Pot,Hubbard_U)

  do I = 1, Nr_atoms
    Coulomb_Force(1,I) = (2.D0*q(I)-n(I))*Coulomb_Force(1,I)/n(I)  ! Adjust forces of linearized XLBOMD functional
    Coulomb_Force(2,I) = (2.D0*q(I)-n(I))*Coulomb_Force(2,I)/n(I) 
    Coulomb_Force(3,I) = (2.D0*q(I)-n(I))*Coulomb_Force(3,I)/n(I) 
  enddo

  FTOT = SKForce + PairForces + FPUL + Coulomb_Force + FSCOUL  ! Total force

  do I = 1, Nr_atoms
    VX(I) = VX(I) + 0.5D0*dt*F2V*FTOT(1,I)/Mnuc(I)       ! second 1/2 of Leapfrog step
    VY(I) = VY(I) + 0.5D0*dt*F2V*FTOT(2,I)/Mnuc(I)       ! second 1/2 of Leapfrog step
    VZ(I) = VZ(I) + 0.5D0*dt*F2V*FTOT(3,I)/Mnuc(I)       ! second 1/2 of Leapfrog step
  enddo

enddo

close(22)
close(23)

end program XLBOMD
