subroutine Kernel_Fermi(KK,D0,mu0,T,RX,RY,RZ,LBox,Hubbard_U,Element_Type,Nr_atoms,MaxIt,eps,m,HDIM, &
Max_Nr_Neigh,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,H_INDEX_START,H_INDEX_END,H,S,Z,Nocc,Znuc)

implicit none
integer, parameter             :: PREC = 8
integer, intent(in)            :: Nr_atoms, HDIM, Nocc,Max_Nr_Neigh
real(PREC), parameter          :: ONE = 1, TWO = 2, ZERO = 0
real(PREC), intent(in)         :: Coulomb_acc, TIMERATIO
real(PREC), intent(in)         :: RX(Nr_atoms), RY(Nr_atoms), RZ(Nr_atoms), LBox(3)
integer, intent(in)            :: H_INDEX_START(Nr_atoms), H_INDEX_END(Nr_atoms)
real(PREC), intent(in)         :: Znuc(Nr_atoms), Hubbard_U(Nr_atoms)
real(PREC), intent(in)         :: H(HDIM,HDIM), S(HDIM,HDIM), Z(HDIM,HDIM)
real(PREC)                     :: H0(HDIM,HDIM), H1(HDIM,HDIM)
real(PREC), intent(out)        :: D0(HDIM,HDIM)
real(PREC)                     :: D(HDIM,HDIM), XX(Nr_atoms,Nr_atoms)
real(PREC)                     :: X(HDIM,HDIM)
real(PREC), intent(in)         :: T, eps, mu0
integer, intent(in)            ::  m, MaxIt
character(10), intent(in)      :: Element_Type(Nr_atoms)
integer, intent(in)            :: nrnnlist(Nr_atoms), nnType(Nr_atoms,Max_Nr_Neigh)
real(PREC), intent(in)         :: nnRx(Nr_atoms,Max_Nr_Neigh)
real(PREC), intent(in)         :: nnRy(Nr_atoms,Max_Nr_Neigh), nnRz(Nr_atoms,Max_Nr_Neigh)
real(PREC)                     :: Coulomb_Pot_Real(Nr_atoms), Coulomb_Force_Real(3,Nr_atoms)
real(PREC)                     :: Coulomb_Force_k(3,Nr_atoms), Coulomb_Pot(Nr_atoms)
real(PREC)                     :: Coulomb_Force_Real_I(3), Coulomb_Pot_Real_I, Coulomb_Pot_k(Nr_atoms)
real(PREC)                     :: Coulomb_Pot_Real_dq_J(Nr_atoms), Coulomb_Pot_dq_J(Nr_atoms)
real(PREC), intent(out)        :: KK(Nr_atoms,Nr_atoms)
real(PREC)                     :: dq_J(Nr_atoms)
real(PREC)                     :: dqI_dqJ(Nr_atoms),mu1
real(PREC)                     :: D_dq_J(HDIM,HDIM), H_dq_J(HDIM,HDIM)
integer                        :: I,J,K, SCF_IT

  XX = ZERO
  X = ZERO
  H1 = ZERO
  Coulomb_Pot_dq_J = ZERO
  Coulomb_Pot_k = ZERO
  dq_J = ZERO
  H_dq_J = ZERO

  do J = 1,Nr_atoms  ! TRIVIAL OPEN_MP OR MPI PARALLELISM
    dq_J(J) = ONE
    do I = 1,Nr_atoms
      call Ewald_Real_Space(Coulomb_Pot_Real_I,Coulomb_Force_Real_I,I,RX,RY,RZ,LBox, &
      dq_J,Hubbard_U,Element_Type,Nr_atoms,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,HDIM,Max_Nr_Neigh)
      Coulomb_Pot_Real(I) = Coulomb_Pot_Real_I
    enddo
    call Ewald_k_Space(Coulomb_Pot_k,Coulomb_Force_k,RX,RY,RZ,LBox,dq_J,Nr_atoms,Coulomb_acc,TIMERATIO,HDIM,Max_Nr_Neigh)
    Coulomb_Pot_dq_J = Coulomb_Pot_Real+Coulomb_Pot_k

    H_dq_J = ZERO
    do I = 1,Nr_atoms
       do K = H_INDEX_START(I),H_INDEX_END(I)
           H_dq_J(K,K) = Hubbard_U(I)*dq_J(I) + Coulomb_Pot_dq_J(I)
       enddo
    enddo
    call MMult(ONE,S,H_dq_J,ZERO,X,'N','N',HDIM)
    call MMult(ONE,H_dq_J,S,ONE,X,'N','T',HDIM)
    H_dq_J =  (ONE/TWO)*X

    call MMult(ONE,Z,H,ZERO,X,'T','N',HDIM)  ! 
    call MMult(ONE,X,Z,ZERO,H0,'N','N',HDIM)   ! H0 = Z'*H*Z

    call MMult(ONE,Z,H_dq_J,ZERO,X,'T','N',HDIM)  ! 
    call MMult(ONE,X,Z,ZERO,H1,'N','N',HDIM)   ! H1 = Z'*H_dq_J*Z

    call Rec_Fermi_PRT1(D0,D_dq_J,H0,H1,T,mu0,mu1,Nocc,m,eps,MaxIt,HDIM)

    call MMult(TWO,Z,D_dq_J,ZERO,X,'N','N',HDIM)
    call MMult(ONE,X,Z,ZERO,D_dq_J,'N','T',HDIM)

    call MMult(ONE,D_dq_J,S,ZERO,X,'N','N',HDIM)
    dqI_dqJ = ZERO
    do I = 1, Nr_atoms
      do K = H_INDEX_START(I), H_INDEX_END(I)
        dqI_dqJ(I) = dqI_dqJ(I) + X(K,K)
      enddo
      XX(I,J) = dqI_dqJ(I)
    enddo
    dq_J = ZERO
  enddo
  do i = 1,Nr_atoms
     XX(i,i) = XX(i,i) - ONE
  enddo
  call Invert(XX,KK,Nr_atoms)

end subroutine Kernel_Fermi
