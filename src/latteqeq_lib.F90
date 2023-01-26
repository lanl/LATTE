!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!    LATTE-XLBOMD Scaled down developers version of LATTE   !!!
!!!    Anders M.N. Niklasson et al. amn@lanl.gov,             !!!
!!!            anders.niklasson@gmail.com                     !!!
!!!  See    Niklasson JCP, 152 104103 (2020) = Ref. [*]       !!!
!!!      Niklasson & Cawkwell JCP 147, 054103 (2017)          !!!
!!!      Niklasson & Cawkwell JCP 141, 164123 (2014)          !!!
!!!     Souvatzis & Niklasson JCP 240, 044117 (2014)          !!!
!!!      Cawkwell & Niklasson JCP 137, 134105 (2012)          !!!
!!!       Niklasson et al. JCP 130, 214109 (2009)             !!!
!!!          Niklasson PRL 100, 123004 (2008)                 !!!
!!!  Niklasson, Tymczak & Challacombe PRL 97, 123001 (2007)   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


MODULE QEQ_LIB


  IMPLICIT NONE

  PRIVATE
  
  
  INTEGER, PARAMETER :: LATTEQEQ_ABIVERSION = 20180622

  PUBLIC :: LATTEQEQ, LATTEQEQ_ABIVERSION

  CONTAINS



!! \param MAXITER Latte MAXITER keyword. If MAXITER = -1, only the Forces are computed.
!!        If MAXITER = 0, MAXITER is read from latte.in file.
!!        IF MAXITER > 0, MAXITER is passed trough the library call.
SUBROUTINE LATTEQEQ(NTYPES, TYPES, CR_IN, MASSES_IN, XLO, XHI, XY, XZ, YZ, FORCES, &
       MAXITER_IN, VENERG, VEL_IN, DT_IN, VIRIAL_INOUT, CURRENTSTEP, GRADX_IN, &
       NEWSYSTEM, EXISTERROR_INOUT) !FNAME)

use qeq_parameters
use omp_lib
implicit none

INTEGER,    INTENT(IN)   :: CURRENTSTEP
REAL(PREC), INTENT(IN)   :: CR_IN(:,:),VEL_IN(:,:), MASSES_IN(:),XLO(3),XHI(3)
REAL(PREC), INTENT(IN)   :: DT_IN, XY, XZ, YZ
REAL(PREC), INTENT(INOUT)  :: FORCES(:,:), VENERG
REAL(PREC), INTENT(IN):: GRADX_IN(:,:)
REAL(PREC), INTENT(OUT)  :: VIRIAL_INOUT(6)
INTEGER, INTENT(IN) ::  NTYPES, TYPES(:), MAXITER_IN
LOGICAL(1), INTENT(INOUT) :: EXISTERROR_INOUT
INTEGER, INTENT(INOUT) :: NEWSYSTEM
!CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: FNAME


integer,    parameter :: Max_Nr_Neigh = 500   ! 100 Water
integer,    parameter :: LL = 16  !! LL = max number of Rank-m
integer,    parameter :: m = 8, MaxIt = 20 
real(PREC), parameter :: T = 300.D0, C_mix = 0.1D0, FelTol = 0.01D0  ! FelTol is convergence required for Relative Err in KK0Res
real(PREC), parameter :: Coulomb_acc = 1e-9, TIMERATIO = 10.D0, ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0, SIX = 6.D0
real(PREC), parameter :: THREE = 3.D0
real(PREC), parameter :: eps = 1e-7, SCF_EPS = 1e-9, dx = 0.001D0  
real(PREC), parameter :: OccErrLim = 1e-9
real(PREC), parameter :: SEVEN = 7.D0

real(PREC) :: mls(8)
integer    :: timevector(8)

real(PREC) :: dt
real(PREC) :: t_start, t_end
integer    :: NATS
integer    :: MD_Iter = 40001

!!!!!!!!!!!!!!!!!!!!!!!!    
integer :: step
integer :: I,J,K,L,SCF_IT, MD_step, Cnt
CHARACTER*20 :: FNAME

FNAME = "latteqeq.in"
OPEN(UNIT=6, FILE=OUTFILE, FORM="formatted")

call cpu_time(t_start)
call date_and_time(values=timevector)
mls(1) = timevector(5)*60*60*1000.D0 + timevector(6)*60*1000 + &
         timevector(7)*1000 + timevector(8)

exact_solution = .False.
printcharges = .False.

INQUIRE(FILE=trim(FNAME), exist=LATTEINEXISTS )
IF (LATTEINEXISTS) THEN
   CALL PARSE_CONTROL(FNAME)
ELSE
   WRITE(6,*) 'latteqeq.in is not found!'
ENDIF

DO I = 1, 3
  LBox(I) = XHI(I) - XLO(I)
enddo

!open(UNIT=22,STATUS="OLD",FILE="MD.xyz")
open(UNIT=23,STATUS="UNKNOWN",FILE="Energy.dat", position="append")

if(printcharges) open(UNIT=24,STATUS="UNKNOWN",FILE="charges.dat", position="append")

NATS = SIZE(CR_IN,DIM=2)
dt = DT_IN * 1000.D0 ! ps to fs

!write(6,*) 'test:currentstep=', currentstep

if (.not.allocated(ATELE)) THEN
  ALLOCATE(ATELE(NATS))
  Allocate(H_INDEX_START(NATS), H_INDEX_END(NATS), NrOrb(NATS))
  ALLOCATE(RX(NATS), RY(NATS), RZ(NATS), q(NATS), qx(NATS), qqx(NATS))
  ALLOCATE(Znuc(NATS), Mnuc(NATS), Hubbard_U(NATS), Xi(NATS), bb(NATS+1))
  ALLOCATE(CC(NATS,NATS), AA(NATS+1,NATS+1), xx(NATS+1))
  ALLOCATE(AAI(NATS+1,NATS+1), KK(NATS,NATS))
  ALLOCATE(nrnnlist(NATS), nnType(NATS,Max_Nr_Neigh))
  ALLOCATE(nnStruct(NATS,NATS), nrnnStruct(NATS))
  ALLOCATE(nndist(NATS,Max_Nr_Neigh), nnRx(NATS,Max_Nr_Neigh))
  ALLOCATE(nnRy(NATS,Max_Nr_Neigh), nnRz(NATS,Max_Nr_Neigh))
  ALLOCATE(Coulomb_Pot_k(NATS))
  ALLOCATE(vi(NATS,LL), v(NATS), q_tmp(NATS), dq_dv(NATS))
  ALLOCATE(fi(NATS,LL), vK(NATS),vJ(NATS))
  ALLOCATE(dr(NATS), IdentRes(NATS), II(NATS,NATS))
  ALLOCATE(Coulomb_Pot_Real(NATS), Coulomb_Force_Real(3,NATS))
  ALLOCATE(Coulomb_Force_k(3,NATS), Coulomb_Pot(NATS), Coulomb_Force(3,NATS))
  ALLOCATE(JJ(NATS,NATS), KK0(NATS,NATS), DD(NATS,NATS))
  ALLOCATE(gq(NATS), KRes(NATS), S_Ent, VV(NATS), PairForces(3,NATS))
  ALLOCATE(SKForce(3,NATS), FPUL(3,NATS), FSCOUL(3,NATS), FTOT(3,NATS))
  ALLOCATE(n(NATS),n_0(NATS),n_1(NATS),n_2(NATS),n_3(NATS),n_4(NATS))
  ALLOCATE(n_5(NATS), n_6(NATS))
  ALLOCATE(VX(NATS), VY(NATS), VZ(NATS))
  ALLOCATE(dn2dt2(NATS), Res(NATS), K0Res(NATS))
ENDIF

!!Get the structure and index lists, unit conversion factors, mass and charge etc.
RX = CR_IN(1,:)
RY = CR_IN(2,:)
RZ = CR_IN(3,:)

!
CALL MASSES2SYMBOLS(TYPES,NTYPES,MASSES_IN,NATS,ATELE, Hubbard_U)


call Initiate_lmp(Coulomb_acc,TIMERATIO,LBox,RX,RY,RZ,NATS, &
  H_INDEX_START,H_INDEX_END,COULCUT,ATELE,Znuc,Mnuc,NrOrb)

HDIM = H_INDEX_END(NATS)
NOCC = 0
DO I = 1, NATS
  NOCC = NOCC + Znuc(I)
ENDDO
Nocc = Nocc/2

if (.not.allocated(D_atomic)) then
  ALLOCATE(D_atomic(HDIM), H0(HDIM,HDIM), H(HDIM,HDIM), H1(HDIM,HDIM))
  ALLOCATE(D0(HDIM,HDIM), D(HDIM,HDIM))
  ALLOCATE(S(HDIM,HDIM), Z(HDIM,HDIM), X(HDIM,HDIM), Y(HDIM,HDIM) )
  ALLOCATE(HOrth(HDIM,HDIM))
  ALLOCATE(dSx(HDIM,HDIM), dSy(HDIM,HDIM), dSz(HDIM,HDIM))
  ALLOCATE(dHx(HDIM,HDIM), dHy(HDIM,HDIM), dHz(HDIM,HDIM))
  ALLOCATE(QQ(HDIM,HDIM), ee(HDIM), Fe_vec(HDIM))  
endif


!! Get the atomic diagonal density matrix
call AtomicDensityMatrix(D_atomic,NATS,H_INDEX_START,H_INDEX_END,HDIM,Znuc) 

call get_qeq_params()

Do i = 1, NATS
  k = TYPES(i)
ENDDO

call date_and_time(values=timevector)
mls(2) = timevector(5)*60*60*1000.D0 + timevector(6)*60*1000 + &
         timevector(7)*1000 + timevector(8)

if (currentstep==0) then
!if (.true.) then
  call nearestneighborlist(nrnnlist,nndist,nnRx,nnRy,nnRz,nnType,nnStruct,nrnnStruct,RX,RY,RZ,LBox,4.0D0, &
   NATS,Max_Nr_Neigh)
  
  mu0 = 0.D0  ! Initial guess that might need ot be adjusted
  
  call nearestneighborlist(nrnnlist,nndist,nnRx,nnRy,nnRz,nnType,nnStruct,nrnnStruct,RX,RY,RZ,LBox,COULCUT, &
   NATS,Max_Nr_Neigh)
  
  ! qeq
  !write(6,*) 'CoulombMatrix calculation'
  call CoulombMatrix(CC,RX,RY,RZ,LBox,Hubbard_U,ATELE,NATS,Coulomb_acc,TIMERATIO, &
                     nnRx,nnRy,nnRz,nrnnlist,nnType,HDIM,Max_Nr_Neigh)
  
  bb = 0.d0
  DO I = 1, NATS
    ! use lammps value
    bb(I) = - GRADX_IN(I,1)
    ! use constant
    !bb(I) = - qeq_xi(nint(MASSES_IN(TYPES(I))))
  ENDDO
  
  AA(1:NATS,1:NATS) = CC
  II = 0.D0
  do I = 1,NATS
    AA(I,I) = AA(I,I) + Hubbard_U(I)
    II(I,I) = 1.D0
  enddo
  AA(1:NATS,NATS+1) = -1.D0
  AA(NATS+1,1:NATS) = 1.D0
  AA(NATS+1,NATS+1) = 0.D0
  call Invert(AA,AAI,NATS+1)
  xx = Matmul(AAI,bb)
  q(1:NATS) = xx(1:NATS)
  
  !! OLD
  !DD = -CC
  !do I = 1,NATS
  !  DD(I,:) = -CC(I,:)/Hubbard_U(I) 
  !enddo
  !do I = 1,NATS
  !  VV = sum(DD(:,I))/(1.D0*NATS)
  !  DD(:,I) = DD(:,I) - VV(:)
  !enddo
  !do I = 1,NATS
  !  DD(I,I) = DD(I,I) - 1.2D0 ! works better in most cases compared to -1.D0 
  !enddo
  !call Invert(DD,KK,NATS)
 
  ! NEW can be simplified
  do I = 1,NATS
  do J = 1,NATS
     DD(I,J) = -CC(I,J)/Hubbard_U(I)
  enddo
  enddo
  temp = 0.D0
  do I = 1,NATS
     temp = temp + 1.D0/Hubbard_U(I)
  enddo

  do I = 1,NATS
  do J = 1,NATS
    DD(I,J) = DD(I,J) + (-sum(DD(:,J))/temp)/Hubbard_U(I)
  enddo
  enddo
  !write(*,*) 'nTestar D = ', sum(DD(:,1))
  !write(*,*) 'nTestar D = ', sum(DD(:,5))
  do I = 1,NATS
    DD(I,I) = DD(I,I) - 1.0D0
  enddo
  call Invert(DD,KK,NATS)
  KK0 = KK

  do I = 1,NATS
    call Ewald_Real_Space(Coulomb_Pot_Real_I,Coulomb_Force_Real_I,I,RX,RY,RZ,LBox, &
    q,Hubbard_U,ATELE,NATS,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,HDIM,Max_Nr_Neigh)
    Coulomb_Pot_Real(I) = Coulomb_Pot_Real_I
    Coulomb_Force_Real(:,I) = Coulomb_Force_Real_I(:)
  enddo
  
  call Ewald_k_Space(Coulomb_Pot_k,Coulomb_Force_k,RX,RY,RZ,LBox,q,NATS,Coulomb_acc,TIMERATIO,Max_Nr_Neigh)

  coulomb_Pot = Coulomb_Pot_Real+Coulomb_Pot_k
  Coulomb_Force = Coulomb_Force_Real + Coulomb_Force_k
  VV = matmul(CC,q)
  
  ECoul = ZERO
  do I = 1,NATS
    ECoul = ECoul - q(I)*bb(I) + 0.5D0*q(I)*Hubbard_U(I)*q(I) + 0.5D0*q(I)*Coulomb_Pot(I)
  enddo
  
  call nearestneighborlist(nrnnlist,nndist,nnRx,nnRy,nnRz,nnType,nnStruct,nrnnStruct,RX,RY,RZ,LBox,4.0D0, &
                           NATS,Max_Nr_Neigh)
  
  EPOT = ECoul
  FTOT = Coulomb_Force 
  
  ! Initial BC for n, using the net Mulliken occupation per atom as extended electronic degrees of freedom
  if (currentstep==0) then
    n = q; n_0 = q; n_1 = q; n_2 = q; n_3 = q; n_4 = q; n_5 = q; n_6 = q;   ! Set all "old" n-vectors to the same at t = t0
    qx = q; qqx = q;
  
    !write(6,*) 'test-zy: current step=', currentstep
  endif

endif

call date_and_time(values=timevector)
mls(3) = timevector(5)*60*60*1000.D0 + timevector(6)*60*1000 + &
         timevector(7)*1000 + timevector(8)

write(6,*) 'Time of latteqeq initialization =', mls(2) - mls(1)
if (currentstep==0) write(6,*) 'Time of latteqeq step 0         =', mls(2) - mls(1)

!!! Coefficients for modified Verlet integration
C0 = -14; C1 = 36; C2 = -27; C3 = -2; C4 = 12;  C5 = -6; C6 = 1;
coeff(1) = -14; coeff(2) = 36; coeff(3) = -27; coeff(4)=-2;coeff(5)=12;coeff(6)=-6;coeff(7)=1;

kappa = 1.84D0; alpha = 0.0055D0;

!!!!! Coefficients for modified Verlet integration
!C0 = -6.D0; C1 = 14.D0; C2 = -8.D0; C3 = -3.D0; C4 = 4.D0; C5 = -1.D0;  C6 = 0.D0;
!kappa = 1.82D0; alpha = 0.018D0;


VX = 0*RX; VY = 0*RX; VZ = 0*RX; !Initialize velocities

VX = VEL_IN(1,1:NATS)/1000.0D0 ! convert from ang/ps to ang/fs
VY = VEL_IN(2,1:NATS)/1000.0D0 ! convert from ang/ps to ang/fs
VZ = VEL_IN(3,1:NATS)/1000.0D0 ! convert from ang/ps to ang/fs


step = 0
Cnt = 0
!! MAIN MD LOOP {dR2(0)/dt2: V(0)->V(1/2); dn2(0)/dt2: n(0)->n(1); V(1/2): R(0)->R(1); dR2(1)/dt2: V(1/2)->V(1)}


if (MAXITER_IN == -1) then
   if (currentstep==0) then
     MD_Iter = 1
     FORCES = FORCES +  FTOT   ! unit ev/ang
   else
     MD_Iter = 1
   endif
else if (MAXITER_IN > 0) THEN
  MD_Iter = MAXITER_IN
endif

IF (MD_Iter == 0) then
  EKIN = ZERO
  do I = 1, NATS
    EKIN = EKIN + 0.5D0*MVV2KE*Mnuc(I)*(VX(I)**2+VY(I)**2+VZ(I)**2)
  enddo
  Temperature = (TWO/THREE)*KE2T*EKIN/NATS
  Energy = EKIN+EPOT
  write(6,*) ' Etotal = ', Energy, 'Temperature = ', Temperature 
endif

do MD_step = 1,MD_Iter
  EKIN = ZERO
  do I = 1, NATS
    EKIN = EKIN + 0.5D0*MVV2KE*Mnuc(I)*(VX(I)**2+VY(I)**2+VZ(I)**2)
  enddo
  Temperature = (TWO/THREE)*KE2T*EKIN/NATS

  Time = (currentstep+MD_STEP-1)*dt

  !Time = (MD_step-1)*dt
  RMSE = norm2(qx-n)/sqrt(ONE*NATS)


  !if (mod(MD_step,1000).EQ.1) then
  if (mod(currentstep,500).EQ.0) then
     call nearestneighborlist(nrnnlist,nndist,nnRx,nnRy,nnRz,nnType,nnStruct,nrnnStruct,RX,RY,RZ,LBox,COULCUT, &
                            NATS,Max_Nr_Neigh)
     !write(6,*) 'CoulombMatrix calculation'
     call CoulombMatrix(CC,RX,RY,RZ,LBox,Hubbard_U,ATELE,NATS,Coulomb_acc,TIMERATIO, &
                      nnRx,nnRy,nnRz,nrnnlist,nnType,HDIM,Max_Nr_Neigh)

     !! old
     !DD = -CC
     !temp = 0.D0
     !do I = 1,NATS
     !  DD(I,:) = -CC(I,:)/Hubbard_U(I)
     !enddo
     !do I = 1,NATS
     !  VV = sum(DD(:,I))/(1.D0*NATS)
     !  DD(:,I) = DD(:,I) - VV(:)
     !enddo
     !do I = 1,NATS
     !  !DD(I,I) = DD(I,I) - 1.D0
     !  DD(I,I) = DD(I,I) - 1.2D0
     !enddo
     !call Invert(DD,KK,NATS)
     !KK0 = KK


     ! NEW
     DD = -CC
     do I = 1,NATS
     do J = 1,NATS
        DD(I,J) = -CC(I,J)/Hubbard_U(I)
     enddo
     enddo
     temp = 0.D0
     do I = 1,NATS
        temp = temp + 1.D0/Hubbard_U(I)
     enddo

     do I = 1,NATS
     do J = 1,NATS
       DD(I,J) = DD(I,J) + (-sum(DD(:,J))/temp)/Hubbard_U(I)
     enddo
     enddo
     do I = 1,NATS
       DD(I,I) = DD(I,I) - 1.0D0
     enddo
     call Invert(DD,KK,NATS)
     KK0 = KK

     !KK_IN = KK0
     write(6,*) 'preconditioner is generated'
  else
     !KK0 = KK_IN
     write(6,*) 'preconditioner is reused'
  endif

  call date_and_time(values=timevector)
  mls(4) = timevector(5)*60*60*1000.D0 + timevector(6)*60*1000 + &
           timevector(7)*1000 + timevector(8)
  write(6,*) 'time of getting kernel KK =', mls(4) - mls(3)

!!! Calculate KRes ~= -(K0*J)^(-1)*K0*(q[n]-n) as in Eq. (43) in Ref. [*] but without the w^2 as in KK0Res
!!! [*] Niklasson JCP, 152 104103 (2020)
  !if (MD_step > 1) then
  if (currentstep > 0) then
  Res = qx-n
  K0Res = MATMUL(KK0,Res)   !! EXTRA STEP FOR PRECONDITIONING
  dr = K0Res               !! Preconditioned residual
  I = 0                    !! Count number of rank updates
  Fel = 1.D0
  do while (Fel > FelTol)  !! Fel = "Error", Could potentially also use a highest number of allowed rank updates
      I = I + 1
      vi(:,I) = dr/norm2(dr)
      do J = 1,I-1
         vi(:,I) = vi(:,I) - dot_product(vi(:,I),vi(:,J))*vi(:,J)  !! Orthogonalized v_i as in Eq. (42) Ref. [*]
      enddo
      vi(:,I) = vi(:,I)/norm2(vi(:,I))
      v(:) = vi(:,I)  ! v_i

      call nearestneighborlist(nrnnlist,nndist,nnRx,nnRy,nnRz,nnType,nnStruct,nrnnStruct,RX,RY,RZ,LBox,COULCUT, &
                           NATS,Max_Nr_Neigh)
      do J = 1,NATS
        call Ewald_Real_Space(Coulomb_Pot_Real_I,Coulomb_Force_Real_I,J,RX,RY,RZ,LBox, &
        v,Hubbard_U,ATELE,NATS,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,HDIM,Max_Nr_Neigh)
        Coulomb_Pot_Real(J) = Coulomb_Pot_Real_I
      enddo
      call Ewald_k_Space(Coulomb_Pot_k,Coulomb_Force_k,RX,RY,RZ,LBox,v,NATS,Coulomb_acc,TIMERATIO,Max_Nr_Neigh)
      Coulomb_Pot = Coulomb_Pot_Real+Coulomb_Pot_k

      ! OLD
      !      temp = 0.D0
      !      do J = 1,NATS
      !         dq_dv(J) = -Coulomb_Pot(J)/Hubbard_U(J)
      !         temp = temp + Coulomb_Pot(J)/Hubbard_U(J)
      !      enddo
      !      dq_dv = dq_dv + temp/(1.D0*NATS)
      
      ! NEW
      temp = 0.D0
      temp2 = 0.D0
      do J = 1,NATS
         dq_dv(J) = -Coulomb_Pot(J)/Hubbard_U(J)
         temp = temp - dq_dv(J)
         temp2 = temp2 + 1.D0/Hubbard_U(J)
      enddo
      do J = 1,NATS
        dq_dv(J) = dq_dv(J) + (temp/temp2)/Hubbard_U(J)
      enddo
 
      dr = dq_dv - v       !! dr = df/dlambda, last row in Eq. (42) Ref[*]
      dr = MATMUL(KK0,dr)  !! dr = K0*(df/dlambda), last row in Eq. (42) Ref[*]
      fi(:,I) = dr  ! fv_i
!
      L = I
      allocate(OO(L,L), MM(L,L))
      do K = 1,L
      do J = 1,L
         OO(K,J) = dot_product(fi(:,K),fi(:,J))   ! O_KJ = < fv_i(K) | fv_i(J) >  see below Eq. (31)
      enddo
      enddo
      call Invert(OO,MM,L)                        ! M = O^(-1)
      IdentRes = 0.D0*K0Res
      KRes = 0.D0
      do K = 1,L
      do J = 1,L
         proj_tmp = MM(K,J)*dot_product(fi(:,J),K0Res(:))
         IdentRes = IdentRes + proj_tmp*fi(:,K)
         KRes = KRes + proj_tmp*vi(:,K)            !! KRes becomes the rank-L approximate of (K0*J)^(-1)*K0*(q[n]-n)
      enddo
      enddo
      !dn2dt2 = KRes
      Fel = norm2(IdentRes-K0Res)/norm2(IdentRes)  !! RELATIVE RESIDUAL ERROR ESTIMATE Eq. (48) Ref. [*]
      if (I.EQ.LL) then
         Fel = 0.D0
      endif
      write(6,*) '# I, L, Fel = ',I,L,Fel          !! Fel goes down with L
      deallocate(OO, MM)
   enddo
   endif
  
   call date_and_time(values=timevector)
   mls(5) = timevector(5)*60*60*1000.D0 + timevector(6)*60*1000 + &
            timevector(7)*1000 + timevector(8)
   write(6,*) 'time of getting KRES = ', mls(5) - mls(4)

   qqx = n - KRes ! Newton-Raphson 
   !write(23,'(f9.4,8f15.6)')  Time/1000, Energy, Temperature, norm2(qx-n)/sqrt(ONE*NATS),n(1),qx(1),qqx(1),q(1), sum(q)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! XL-BOMD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Integrate XL-BOMD equation of motion for charge density, kappa*d^2n/dt^2 = -dt^2*w^2*K*(q-n) !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dn2dt2 = matmul(KK0,(qx-n))
  dn2dt2 = KRes
!  if (MD_step>100) then
!     alpha = 0.D0
!  endif

  n = 2*n_0 - n_1 - kappa*dn2dt2 + alpha*(C0*n_0+C1*n_1+C2*n_2+C3*n_3+C4*n_4+C5*n_5+C6*n_6)
  n_6 = n_5; n_5 = n_4; n_4 = n_3; n_3 = n_2; n_2 = n_1; n_1 = n_0; n_0 = n

  !!! Exact solution
  if(exact_solution) then
  call nearestneighborlist(nrnnlist,nndist,nnRx,nnRy,nnRz,nnType,nnStruct,nrnnStruct,RX,RY,RZ,LBox,COULCUT, &
                            NATS,Max_Nr_Neigh)
  call CoulombMatrix(CC,RX,RY,RZ,LBox,Hubbard_U,ATELE,NATS,Coulomb_acc,TIMERATIO, &
                      nnRx,nnRy,nnRz,nrnnlist,nnType,HDIM,Max_Nr_Neigh)
  AA(1:NATS,1:NATS) = CC
  do I = 1,NATS
    AA(I,I) = AA(I,I) + Hubbard_U(I)
  enddo
  AA(NATS+1, NATS+1) = 0.D0

  call Invert(AA,AAI,NATS+1)
  xx = Matmul(AAI,bb)
  q(1:NATS) = xx(1:NATS)
  
  call date_and_time(values=timevector)
  mls(6) = timevector(5)*60*60*1000.D0 + timevector(6)*60*1000 + &
           timevector(7)*1000 + timevector(8)
  write(6,*) 'time of getting exact solution = ', mls(6) - mls(5)
  endif
  
  call nearestneighborlist(nrnnlist,nndist,nnRx,nnRy,nnRz,nnType,nnStruct,nrnnStruct,RX,RY,RZ,LBox,COULCUT, &
                           NATS,Max_Nr_Neigh)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,Coulomb_Pot_Real_I,Coulomb_Force_Real_I)
  do I = 1,NATS
    call Ewald_Real_Space(Coulomb_Pot_Real_I,Coulomb_Force_Real_I,I,RX,RY,RZ,LBox, &
    n,Hubbard_U,ATELE,NATS,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,HDIM,Max_Nr_Neigh)
    Coulomb_Pot_Real(I) = Coulomb_Pot_Real_I
    Coulomb_Force_Real(:,I) = Coulomb_Force_Real_I(:)
  enddo
!$OMP END PARALLEL DO
  call date_and_time(values=timevector)
  mls(7) = timevector(5)*60*60*1000.D0 + timevector(6)*60*1000 + &
           timevector(7)*1000 + timevector(8)
  write(6,*) 'time of getting ewald_real = ', mls(7) - mls(6)
  
  call Ewald_k_Space(Coulomb_Pot_k,Coulomb_Force_k,RX,RY,RZ,LBox,n,NATS,Coulomb_acc,TIMERATIO,Max_Nr_Neigh)
  Coulomb_Pot = Coulomb_Pot_Real+Coulomb_Pot_k
  Coulomb_Force = Coulomb_Force_Real + Coulomb_Force_k

  lambda = 0.D0
  temp = 0.D0
  do I = 1,NATS
     lambda = lambda + (-bb(I) + Coulomb_Pot(I))/Hubbard_U(I)
     temp = temp + 1./Hubbard_U(I)
  enddo
  lambda = lambda/temp

  !write(6,*) 'lambda =', lambda 
  do I = 1,NATS
     qx(I) = (bb(I)-Coulomb_Pot(I) + lambda)/Hubbard_U(I)
  enddo

  if (printcharges) then
     if (exact_solution) then
        write(24,'(f9.4,10000f15.6)')  Time/1000, qx, q, n
     else
        write(24,'(f9.4,10000f15.6)')  Time/1000, qx 
     endif
  endif

  ECoul = ZERO
  do I = 1,NATS
    ECoul = ECoul  - qx(I)*bb(I) + 0.5D0*qx(I)*Hubbard_U(I)*qx(I) + 0.5D0*(2.D0*qx(I)-n(I))* Coulomb_Pot(I)
  enddo

  EPOT = ECoul !+ energies from charge negativities

  do I = 1, NATS
    Coulomb_Force(1,I) = (2.D0*qx(I)-n(I))*Coulomb_Force(1,I)/n(I)  ! Adjust forces of linearized XLBOMD functional
    Coulomb_Force(2,I) = (2.D0*qx(I)-n(I))*Coulomb_Force(2,I)/n(I) 
    Coulomb_Force(3,I) = (2.D0*qx(I)-n(I))*Coulomb_Force(3,I)/n(I) 
  enddo

  FTOT = Coulomb_Force

  FORCES = FORCES + FTOT

  EKIN = ZERO
  do I = 1, NATS
    EKIN = EKIN + 0.5D0*MVV2KE*Mnuc(I)*(VX(I)**2+VY(I)**2+VZ(I)**2)
  enddo
  Temperature = (TWO/THREE)*KE2T*EKIN/NATS
  Energy = EKIN+EPOT
  VENERG = EPOT
   
  call date_and_time(values=timevector)
  mls(8) = timevector(5)*60*60*1000.D0 + timevector(6)*60*1000 + &
           timevector(7)*1000 + timevector(8)
  write(6,*) 'time of getting coulomb force = ', mls(8) - mls(7)
  write(23,'(f9.4,8f15.6)')  Time/1000, Energy, EPOT, Temperature, norm2(qx-n)/sqrt(ONE*NATS), sum(q)

  write(6,*) ' Time = ', Time, ' Etotal = ', Energy, 'Temperature = ', Temperature 
  write(6,*) ' ----- RMS = ',  norm2(qx-n)/sqrt(ONE*NATS), 'Nr Resp Cal = ', Cnt*1.D0/(MD_step*1.D0),  &
             ' Gap = ', ee(Nocc+1)-ee(Nocc)
enddo

!close(22)
close(23)
if(printcharges) close(24)

call cpu_time(t_end)
write(6,*) 'cpu time of latteqeq=', t_end - t_start

end subroutine LATTEQEQ

END MODULE QEQ_LIB
