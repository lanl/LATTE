subroutine Kernel_Fermi_NoDiag(KK,JJ,D0,mu0,mu1,T,RX,RY,RZ,LBox,Hubbard_U,Element_Type, &
Nr_atoms,MaxIt,eps, m,HDIM, Max_Nr_Neigh,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz, & 
nrnnlist,nnType,H_INDEX_START, H_INDEX_END,H,S,Z,Nocc,Znuc,QQ,ee,Fe_vec)

implicit none
integer,    parameter          :: PREC = 8
integer,    intent(in)         :: Nr_atoms, HDIM, Nocc, Max_Nr_Neigh
real(PREC), parameter          :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0
real(PREC), parameter          :: kB = 8.61739d-5 ! eV/K, kB = 6.33366256e-6 Ry/K, kB = 3.166811429e-6 Ha/K
real(PREC), intent(in)         :: Coulomb_acc, TIMERATIO
real(PREC)                     :: v(Nr_atoms)
real(PREC), intent(in)         :: RX(Nr_atoms), RY(Nr_atoms), RZ(Nr_atoms), LBox(3)
integer,    intent(in)         :: H_INDEX_START(Nr_atoms), H_INDEX_END(Nr_atoms)
real(PREC), intent(in)         :: Znuc(Nr_atoms), Hubbard_U(Nr_atoms)
real(PREC), intent(in)         :: H(HDIM,HDIM), S(HDIM,HDIM), Z(HDIM,HDIM)
real(PREC)                     :: H0(HDIM,HDIM), H1(HDIM,HDIM)
real(PREC), intent(inout)      :: D0(HDIM,HDIM)
real(PREC)                     :: X(HDIM,HDIM), XX(Nr_atoms,Nr_atoms)
real(PREC), intent(in)         :: T, eps
real(PREC), intent(inout)      :: mu0, mu1
integer,    intent(in)         :: m, MaxIt
character(10), intent(in)      :: Element_Type(Nr_atoms)
integer,    intent(in)         :: nrnnlist(Nr_atoms), nnType(Nr_atoms,Max_Nr_Neigh)
real(PREC), intent(in)         :: nnRx(Nr_atoms,Max_Nr_Neigh)
real(PREC), intent(in)         :: nnRy(Nr_atoms,Max_Nr_Neigh), nnRz(Nr_atoms,Max_Nr_Neigh)
real(PREC)                     :: Coulomb_Pot_Real(Nr_atoms), Coulomb_Pot(Nr_atoms)
real(PREC)                     :: Coulomb_Pot_Real_I, Coulomb_Pot_k(Nr_atoms)
real(PREC)                     :: Coulomb_Pot_Real_dq_v(Nr_atoms), Coulomb_Pot_dq_v(Nr_atoms)
real(PREC)                     :: Coulomb_Force_Real_I(3),Coulomb_Force_k(3,Nr_atoms)
real(PREC)                     :: D_dq_v(HDIM,HDIM), H_dq_v(HDIM,HDIM), dq_v(Nr_atoms)
real(PREC)                     :: dq_dv(Nr_atoms)
real(PREC), intent(out)        :: KK(Nr_atoms,Nr_atoms),JJ(Nr_atoms,Nr_atoms)
real(PREC), intent(in)         :: QQ(HDIM,HDIM), ee(HDIM), Fe_vec(HDIM)
integer                        :: I,J,K, ITER, mm
real                           :: t1, t2
!!!
  XX = ZERO
  X = ZERO
  H1 = ZERO
  Coulomb_Pot_dq_v = ZERO
  Coulomb_Pot_k = ZERO
  dq_v = ZERO
  H_dq_v = ZERO
  JJ = ZERO

  call MMult(ONE,Z,H,ZERO,X,'T','N',HDIM)  ! 
  call MMult(ONE,X,Z,ZERO,H0,'N','N',HDIM)   ! H0 = Z'*H*Z

!!!!!!!!! PROBLEMS WITH PERFORMANCE FOR OPENMP, IT MAKES IT MUCH SLOWER?

!!!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,K,dq_dv,Coulomb_Pot_Real_I) &
!!!$OMP PRIVATE(Coulomb_Pot_Real,Coulomb_Pot_k,Coulomb_Pot_dq_v,H_dq_v,X,H0,H1,D_dq_v) &
!!!$OMP PRIVATE(Coulomb_Force_k,Coulomb_Force_Real_I)

!!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(JJ,Hubbard_U,Element_Type,Nr_atoms,Coulomb_acc,TIMERATIO) &
!!$OMP SHARED(nnRx,nnRy,nnRz,nrnnlist,nnType,HDIM,Max_Nr_Neigh,S,Z,RX,RY,RZ,LBox) &
!!$OMP SHARED(Nocc,T,QQ,ee,Fe_vec,mu0,eps)
  do J = 1,Nr_atoms  ! TRIVIAL OPEN_MP OR MPI PARALLELISM
    dq_v(J) = ONE
    do I = 1,Nr_atoms
      call Ewald_Real_Space(Coulomb_Pot_Real_I,Coulomb_Force_Real_I,I,RX,RY,RZ,LBox, &
      dq_v,Hubbard_U,Element_Type,Nr_atoms,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist, &
      nnType,HDIM,Max_Nr_Neigh)
      Coulomb_Pot_Real(I) = Coulomb_Pot_Real_I
    enddo
    call Ewald_k_Space(Coulomb_Pot_k,Coulomb_Force_k,RX,RY,RZ,LBox,dq_v,Nr_atoms,Coulomb_acc, &
    TIMERATIO,HDIM,Max_Nr_Neigh)
    Coulomb_Pot_dq_v = Coulomb_Pot_Real+Coulomb_Pot_k

    H_dq_v = ZERO
    do I = 1,Nr_atoms
       do K = H_INDEX_START(I),H_INDEX_END(I)
           H_dq_v(K,K) = Hubbard_U(I)*dq_v(I) + Coulomb_Pot_dq_v(I)
       enddo
    enddo
    call MMult(ONE,S,H_dq_v,ZERO,X,'N','N',HDIM)
    call MMult(ONE,H_dq_v,S,ONE,X,'N','T',HDIM)
    H_dq_v =  (ONE/TWO)*X

    call MMult(ONE,Z,H_dq_v,ZERO,X,'T','N',HDIM)  ! 
    call MMult(ONE,X,Z,ZERO,H1,'N','N',HDIM)   ! H1 = Z'*H_dq_v*Z

    call get_deriv_finite_temp(D_dq_v,H0,H1,Nocc,T,QQ,ee,Fe_vec,mu0,eps,HDIM)

    call MMult(TWO,Z,D_dq_v,ZERO,X,'N','N',HDIM)
    call MMult(ONE,X,Z,ZERO,D_dq_v,'N','T',HDIM)

    call MMult(ONE,D_dq_v,S,ZERO,X,'N','N',HDIM)
    dq_dv = ZERO
    do I = 1, Nr_atoms
      do K = H_INDEX_START(I), H_INDEX_END(I)
        dq_dv(I) = dq_dv(I) + X(K,K)
      enddo
      JJ(I,J) = dq_dv(I)
    enddo
    dq_v = ZERO
  enddo
!!$OMP END PARALLEL DO
  XX = JJ
  do I = 1,Nr_atoms
     XX(I,I) = XX(I,I) - ONE
  enddo
  call Invert(XX,KK,Nr_atoms)
end subroutine Kernel_Fermi_NoDiag
