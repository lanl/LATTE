!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CACULATES THE PRECONDITIONED MULTI-RANK APPROXIMATION OF KERNEL ACTING ON THE RESIDUAL !!
!!                  THEORY GIVEN IN Niklasson, JCP 152, 104103 (2020) [*]                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine KK0Res(KRes,KK0,Res,FelTol,L,H0,mu0,T,RX,RY,RZ,LBox,Hubbard_U,Element_Type,Nr_atoms,MaxIt,eps,nr_mom,&
       HDIM,Max_Nr_Neigh,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,H_INDEX_START,H_INDEX_END, &
       H,S,Z,Nocc,Znuc,QQ,ee,Fe_vec)

!! Res = q[n] - n
!! KK0 is preconditioner
!! KRes rank-L approximation of (K0*J)^(-1)*K0*(q[n]-n) with (K0*J)^(-1) as in Eq. (41) in Ref. [*]
!! QQ, ee Eigenvectors and eigen values of H0
!! Fe_vec Fermi occupation factors

implicit none
integer, parameter          :: PREC = 8
integer,    intent(in)      :: Nr_atoms, HDIM, Nocc,Max_Nr_Neigh
real(PREC), parameter       :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0
real(PREC), parameter       :: kB = 8.61739d-5 ! eV/K, kB = 6.33366256e-6 Ry/K, kB = 3.166811429e-6 Ha/K
real(PREC), intent(in)      :: Coulomb_acc, TIMERATIO, FelTol
real(PREC), intent(in)      :: RX(Nr_atoms), RY(Nr_atoms), RZ(Nr_atoms), LBox(3)
real(PREC), intent(in)      :: KK0(Nr_atoms,Nr_atoms), Res(Nr_atoms)
integer,    intent(in)      :: H_INDEX_START(Nr_atoms), H_INDEX_END(Nr_atoms)
real(PREC), intent(in)      :: Znuc(Nr_atoms), Hubbard_U(Nr_atoms)
real(PREC), intent(in)      :: H(HDIM,HDIM), S(HDIM,HDIM), Z(HDIM,HDIM), H0(HDIM,HDIM)
real(PREC)                  :: H1(HDIM,HDIM), HX(HDIM,HDIM), K0Res(Nr_atoms)
real(PREC)                  :: D0(HDIM,HDIM)
real(PREC)                  :: X(HDIM,HDIM)
real(PREC), intent(in)      :: T, eps
real(PREC), intent(inout)   :: mu0
integer,    intent(in)      :: nr_mom, MaxIt
character(10), intent(in)   :: Element_Type(Nr_atoms)
integer,    intent(in)      :: nrnnlist(Nr_atoms), nnType(Nr_atoms,Max_Nr_Neigh)
real(PREC), intent(in)      :: nnRx(Nr_atoms,Max_Nr_Neigh)
real(PREC), intent(in)      :: nnRy(Nr_atoms,Max_Nr_Neigh), nnRz(Nr_atoms,Max_Nr_Neigh)
real(PREC)                  :: Coulomb_Pot_Real(Nr_atoms), Coulomb_Pot(Nr_atoms)
real(PREC)                  :: Coulomb_Pot_Real_I, Coulomb_Pot_k(Nr_atoms),dq_dv(Nr_atoms)
real(PREC)                  :: Coulomb_Pot_Real_dq_v(Nr_atoms), Coulomb_Pot_dq_v(Nr_atoms)
real(PREC)                  :: Coulomb_Force_Real_I(3),Coulomb_Force_k(3,Nr_atoms)
real(PREC)                  :: D_dq_v(HDIM,HDIM), H_dq_v(HDIM,HDIM), dq_v(Nr_atoms), ZQ(HDIM,HDIM)
real(PREC), intent(out)     :: KRes(Nr_atoms)
real(PREC), intent(in)      :: QQ(HDIM,HDIM), ee(HDIM), Fe_vec(HDIM)
integer                     :: I,J,K
integer, intent(out)        :: L  !! Number of rank updates used to reach below threshold
real(PREC)                  :: Fel,mu1
real(PREC)                  :: vi(Nr_atoms,Nr_atoms), v(Nr_atoms), q_tmp(Nr_atoms)
real(PREC)                  :: fi(Nr_atoms,Nr_atoms)
real(PREC)                  :: dr(Nr_Atoms), proj_tmp, IdentRes(Nr_atoms)
real(PREC), allocatable     :: O(:,:),M(:,:)

     K0Res = MATMUL(KK0,Res)  !! EXTRA STEP FOR PRECONDITIONING
     dr = K0Res               !! Preconditioned residual
     I = 0                    !! Count number of rank updates
     Fel = 1.D0
     do while (Fel > FelTol)   !! Fel = "Error" in Swedish, Could potentially also use a highest number of allowed rank updates
        I = I + 1

        vi(:,I) = dr/norm2(dr)
        do J = 1,I-1
           vi(:,I) = vi(:,I) - dot_product(vi(:,I),vi(:,J))*vi(:,J)  !! Orthogonalized v_i as in Eq. (42) Ref. [*]
        enddo
        vi(:,I) = vi(:,I)/norm2(vi(:,I))
        v(:) = vi(:,I)  ! v_i

        !! Calculated dq_dv, which is the response in q(n) from change in directional input charge n = v
        call V_Kernel_Fermi_No_Diag(D0,dq_dv,v,H0,mu0,mu1,T,RX,RY,RZ,LBox,Hubbard_U,Element_Type,Nr_atoms, &
        MaxIt,eps,nr_mom,HDIM,Max_Nr_Neigh,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType, &
        H_INDEX_START,H_INDEX_END,H,S,Z,Nocc,Znuc,QQ,ee,Fe_vec)

        dr = dq_dv - v       !! dr = df/dlambda, last row in Eq. (42) Ref[*]
        dr = MATMUL(KK0,dr)  !! dr = K0*(df/dlambda), last row in Eq. (42) Ref[*]
        fi(:,I) = dr  ! fv_i

        L = I
        allocate(O(L,L), M(L,L))
        do K = 1,L
        do J = 1,L
           O(K,J) = dot_product(fi(:,K),fi(:,J))  ! O_KJ = < fv_i(K) | fv_i(J) >  see below Eq. (31)
        enddo
        enddo
        call Invert(O,M,L)                        ! M = O^(-1)
        IdentRes = 0.D0*K0Res
        KRes = 0.D0
        do K = 1,L
        do J = 1,L
           proj_tmp = M(K,J)*dot_product(fi(:,J),K0Res)
           IdentRes = IdentRes + proj_tmp*fi(:,K)
           KRes = KRes + proj_tmp*vi(:,K)            !! KRes becomes the rank-L approximate of (K0*J)^(-1)*K0*(q[n]-n) 
        enddo
        enddo
        Fel = norm2(IdentRes-K0Res)/norm2(IdentRes)  !! RELATIVE RESIDUAL ERROR ESTIMATE Eq. (48) Ref. [*]
        write(*,*) '# I, L, Fel = ',I,L,Fel, '  Norm(IdRes) = ', norm2(IdentRes)         !! Fel goes down with L
        deallocate(O, M)
     enddo

end subroutine KK0Res
