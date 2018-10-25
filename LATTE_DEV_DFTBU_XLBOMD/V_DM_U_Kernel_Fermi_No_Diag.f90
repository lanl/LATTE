subroutine V_DM_U_Kernel_Fermi_No_Diag(D_dq_v,DV,H0,mu0,mu1,T,RX,RY,RZ,LBox,Hubbard_U,Element_Type,Nr_atoms, &
MaxIt,eps,m,HDIM,Max_Nr_Neigh,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,H_INDEX_START, &
H_INDEX_END,H,H_1,S,Z,Nocc,Znuc,QQ,ee,Fe_vec)

implicit none
integer, parameter          :: PREC = 8
integer,    intent(in)      :: Nr_atoms, HDIM, Nocc,Max_Nr_Neigh
real(PREC), parameter       :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0
real(PREC), parameter       :: kB = 8.61739d-5 ! eV/K, kB = 6.33366256e-6 Ry/K, kB = 3.166811429e-6 Ha/K
real(PREC), intent(in)      :: Coulomb_acc, TIMERATIO
real(PREC), intent(in)      :: RX(Nr_atoms), RY(Nr_atoms), RZ(Nr_atoms), LBox(3)
integer,    intent(in)      :: H_INDEX_START(Nr_atoms), H_INDEX_END(Nr_atoms)
real(PREC), intent(in)      :: Znuc(Nr_atoms), Hubbard_U(Nr_atoms)
real(PREC), intent(in)      :: H(HDIM,HDIM), H_1(HDIM,HDIM), S(HDIM,HDIM), Z(HDIM,HDIM), H0(HDIM,HDIM)  ! H_1 in Non-Orth repre
real(PREC)                  :: H1(HDIM,HDIM), HX(HDIM,HDIM),v(Nr_atoms)
!real(PREC), intent(inout)   :: D0(HDIM,HDIM)
real(PREC)                  :: X(HDIM,HDIM)
real(PREC), intent(in)      :: T, eps, DV(HDIM,HDIM)
real(PREC), intent(inout)  :: mu0, mu1
integer,    intent(in)      :: m, MaxIt
character(10), intent(in)   :: Element_Type(Nr_atoms)
integer,    intent(in)      :: nrnnlist(Nr_atoms), nnType(Nr_atoms,Max_Nr_Neigh)
real(PREC), intent(in)      :: nnRx(Nr_atoms,Max_Nr_Neigh)
real(PREC), intent(in)      :: nnRy(Nr_atoms,Max_Nr_Neigh), nnRz(Nr_atoms,Max_Nr_Neigh)
real(PREC)                  :: Coulomb_Pot_Real(Nr_atoms), Coulomb_Pot(Nr_atoms)
real(PREC)                  :: Coulomb_Pot_Real_I, Coulomb_Pot_k(Nr_atoms)
real(PREC)                  :: Coulomb_Pot_Real_dq_v(Nr_atoms), Coulomb_Pot_dq_v(Nr_atoms)
real(PREC)                  :: Coulomb_Force_Real_I(3),Coulomb_Force_k(3,Nr_atoms)
real(PREC)                  :: H_dq_v(HDIM,HDIM), dq_v(Nr_atoms)
!real(PREC), intent(out)     :: dq_dv(Nr_atoms)
real(PREC), intent(out)     :: D_dq_v(HDIM,HDIM)
real(PREC), intent(in)      :: QQ(HDIM,HDIM), ee(HDIM), Fe_vec(HDIM)
integer                     :: I,J,K, ITER, mm
real                        :: t1, t2

    dq_v = 0.D0
    do I = 1,Nr_atoms
      do K = H_INDEX_START(I),H_INDEX_END(I)
        dq_v(I) = dq_v(I) + 2*DV(K,K) ! charge response from perturbation in DM
      enddo
    enddo

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,Coulomb_Pot_Real_I,Coulomb_Force_Real_I)
    do I = 1,Nr_atoms
      call Ewald_Real_Space(Coulomb_Pot_Real_I,Coulomb_Force_Real_I,I,RX,RY,RZ,LBox, &
      dq_v,Hubbard_U,Element_Type,Nr_atoms,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,HDIM,Max_Nr_Neigh)
      Coulomb_Pot_Real(I) = Coulomb_Pot_Real_I
    enddo
!$OMP END PARALLEL DO
    call Ewald_k_Space(Coulomb_Pot_k,Coulomb_Force_k,RX,RY,RZ,LBox,dq_v,Nr_atoms,Coulomb_acc,TIMERATIO,HDIM,Max_Nr_Neigh)
    Coulomb_Pot_dq_v = Coulomb_Pot_Real + Coulomb_Pot_k  ! Response in Coulomb potential

    H_dq_v = ZERO
    do I = 1,Nr_atoms
       do K = H_INDEX_START(I),H_INDEX_END(I)
           H_dq_v(K,K) = Hubbard_U(I)*dq_v(I) + Coulomb_Pot_dq_v(I) ! Response in Hamiltonian
       enddo
    enddo
    call MMult(ONE,S,H_dq_v,ZERO,X,'N','N',HDIM) ! Dense times diagonal can maybe be done faster
    call MMult(ONE,H_dq_v,S,ONE,X,'N','T',HDIM)
    H_dq_v =  (ONE/TWO)*X + H_1  ! Include also H_1, i.e. the direct response from the Hubbard U term (not via Coulomb)

    call MMult(ONE,Z,H_dq_v,ZERO,X,'T','N',HDIM)  ! 
    call MMult(ONE,X,Z,ZERO,H1,'N','N',HDIM)   ! H1 = Z'*H_dq_v*Z

    call get_deriv_finite_temp(D_dq_v,H0,H1,Nocc,T,QQ,ee,Fe_vec,mu0,eps,HDIM) ! Finite temperature quantum perturabtion theory
  
end subroutine V_DM_U_Kernel_Fermi_No_Diag
