!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Simplified first-order Canonical Density Matrix Perturbation Theory, !!!
!!! assuming known mu0 and representation in H0's eigenbasis Q,          !!!
!!! where H0 and P0 become diagonal.                                     !!!
!!! (mu0, eigenvalues e and eigenvectors Q of H0 are assumed known)      !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!            Based on Niklasson, PRL 92, 193001 (2004), and            !!!
!!! Niklasson, Cawkwell, Rubensson, and Rudberg, PRE 92, 063301 (2015)   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine  Canon_DM_PRT(H1,beta,Q,e,mu0,m,HDIM)

  USE MYPRECISION
  USE SETUPARRAY

  implicit none
  integer, parameter      :: PREC = LATTEPREC

  !real(PREC), parameter   :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0
  integer, intent(in)     :: HDIM, m ! m = Number of recursion steps
  real(PREC), intent(in)  :: H1(HDIM,HDIM), Q(HDIM,HDIM), e(HDIM) ! Q and e are eigenvectors and eigenvalues of H0
  real(PREC), intent(in)  :: beta, mu0 ! Electronic temperature and chemical potential
!
  real(PREC), ALLOCATABLE :: P1(:,:)                  ! Density matrix response derivative with respect to perturbation H1 to H0
  real(PREC), ALLOCATABLE :: X(:,:), DX1(:,:), Y(:,:) ! Temporary matrices
  real(PREC), ALLOCATABLE :: h_0(:), p_0(:), dPdmu(:), p_02(:), iD0(:)
!
  real(PREC)              :: cnst, mu1
  real(PREC)              :: traceBO
  integer                 :: i, j, k, N 

  ALLOCATE(P1(HDIM,HDIM)) ! Density matrix response derivative with respect to perturbation H1 to H0
  ALLOCATE(X(HDIM,HDIM), DX1(HDIM,HDIM), Y(HDIM,HDIM)) ! Temporary matrices
  ALLOCATE(h_0(HDIM), p_0(HDIM), dPdmu(HDIM), p_02(HDIM), iD0(HDIM))

  h_0 = e                ! Diagonal Hamiltonian H0 respresented in the eigenbasis Q
  cnst = beta/(1.D0*2**(m+2)) ! Scaling constant

  p_0 = 0.5D0 + cnst*(h_0-mu0)  ! Initialization for P0 represented in eigenbasis Q

  N = HDIM
  CALL DGEMM('T', 'N', HDIM, HDIM, HDIM, ONE, &
          Q, HDIM, H1, HDIM, ZERO, X, HDIM)
  CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
          X, HDIM, Q, HDIM, ZERO, Y, HDIM)

!  call MMult(ONE,Q,H1,ZERO,X,'T','N',HDIM)      ! Main cost are the transformation 
!  call MMult(ONE,X,Q,ZERO,Y,'N','N',HDIM)       ! to the eigenbasis of H1
!  P1 = -cnst*H1    !(set mu1 = 0 for simplicity) ! Initialization of DM response in Q representation (not diagonal in Q)
  P1 = -cnst*Y    !(set mu1 = 0 for simplicity) ! Initialization of DM response (not diagonal in Q)

  do i = 1,m  ! Loop over m recursion steps
    p_02 = p_0*p_0
    do j = 1,HDIM
      do k = 1,HDIM
        DX1(k,j) = p_0(k)*P1(k,j) + P1(k,j)*p_0(j)
      enddo
    enddo
    iD0 = 1.D0/(2.D0*(p_02-p_0)+1.D0)
    p_0 = iD0*p_02
    do j = 1,HDIM
      do k = 1,HDIM
        P1(k,j) = iD0(k)*(DX1(k,j) + 2.D0*(P1(k,j)-DX1(k,j))*p_0(j))
      enddo
    enddo
  enddo
  dPdmu = beta*p_0*(1.D0-p_0)
  mu1 = 0.D0
  do i = 1,HDIM
    mu1 = mu1 + P1(i,i)
  enddo
  mu1 = -mu1/SUM(dPdmu)
  do i = 1,HDIM
    P1(i,i) = P1(i,i) + mu1*dPdmu(i)  ! Trace correction by adding (dP/dmu)*(dmu/dH1) to dP/dH1
  enddo

  CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
          Q, HDIM, P1, HDIM, ZERO, X, HDIM)
  CALL DGEMM('N', 'T', HDIM, HDIM, HDIM, ONE, &
          X, HDIM, Q, HDIM, ZERO, BO, HDIM)

  traceBO = 0.D0
  do i = 1, HDIM
    traceBO = traceBO + BO(I,I)
  enddo
  !write(6,*) 'trace of BO (DM) = ', traceBO

!  call MMult(ONE,Q,P1,ZERO,X,'N','N',HDIM) ! and back transformation of P1, with a total of
!  call MMult(ONE,X,Q,ZERO,P1,'N','T',HDIM) ! 4 matrix-matrix multiplication dominates the cost

  DEALLOCATE(P1, X, DX1, Y, H_0, P_0, DPDMU, P_02, ID0)

end subroutine Canon_DM_PRT


subroutine  Canon_DM_PRT_SPIN(H1UP,H1DOWN,beta,QUP,QDOWN,EUP,EDOWN,mu0,m,HDIM)

  USE MYPRECISION
  USE SETUPARRAY
  USE SPINARRAY, ONLY: RHOUP, RHODOWN

  implicit none
  integer, parameter      :: PREC = LATTEPREC

  integer, intent(in)     :: HDIM, m                                          ! m = Number of recursion steps
  real(PREC), intent(in)  :: H1UP(HDIM,HDIM), QUP(HDIM,HDIM), EUP(HDIM)       ! Q and e are eigenvectors and eigenvalues of H0
  real(PREC), intent(in)  :: H1DOWN(HDIM,HDIM), QDOWN(HDIM,HDIM), EDOWN(HDIM) ! Q and e are eigenvectors and eigenvalues of H0
  real(PREC), intent(in)  :: beta, mu0                                        ! Electronic temperature and chemical potential
  !
  ! Density matrix response derivative with respect to perturbation H1 to H0
  real(PREC), ALLOCATABLE :: P1UP(:,:), P1DOWN(:,:)
  real(PREC), ALLOCATABLE :: X(:,:), DX1(:,:), Y(:,:) ! Temporary matrices
  real(PREC), ALLOCATABLE :: h_0(:), p_0(:), p_02(:), iD0(:)
  REAL(PREC), ALLOCATABLE :: dPdmu(:), dpdmu_down(:) 
  !
  real(PREC)              :: cnst, mu1
  real(PREC)              :: traceBO
  integer                 :: i, j, k, N 

  N = HDIM

  ! Density matrix response derivative with respect to perturbation H1 to H0
  ALLOCATE(P1UP(HDIM,HDIM), P1DOWN(HDIM,HDIM)) 
  ALLOCATE(X(HDIM,HDIM), DX1(HDIM,HDIM), Y(HDIM,HDIM)) ! Temporary matrices
  ALLOCATE(h_0(HDIM), p_0(HDIM), dPdmu(HDIM), dPdmu_Down(HDIM), p_02(HDIM), iD0(HDIM))

  !1) SPIN UP
  h_0 = EUP                   ! Diagonal Hamiltonian H0 respresented in the eigenbasis Q
  cnst = beta/(1.D0*2**(m+2)) ! Scaling constant

  p_0 = 0.5D0 + cnst*(h_0-mu0)  ! Initialization for P0 represented in eigenbasis Q

  CALL DGEMM('T', 'N', HDIM, HDIM, HDIM, ONE, QUP, HDIM, H1UP, HDIM, ZERO, X, HDIM)
  CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, X,   HDIM, QUP,  HDIM, ZERO, Y, HDIM)

  P1UP = -cnst*Y    !(set mu1 = 0 for simplicity) ! Initialization of DM response (not diagonal in Q)

  do i = 1,m  ! Loop over m recursion steps
    p_02 = p_0*p_0
    do j = 1,HDIM
      do k = 1,HDIM
        DX1(k,j) = p_0(k)*P1UP(k,j) + P1UP(k,j)*p_0(j)
      enddo
    enddo
    iD0 = 1.D0/(2.D0*(p_02-p_0)+1.D0)
    p_0 = iD0*p_02
    do j = 1,HDIM
      do k = 1,HDIM
        P1UP(k,j) = iD0(k)*(DX1(k,j) + 2.D0*(P1UP(k,j)-DX1(k,j))*p_0(j))
      enddo
    enddo
  enddo

  dPdmu = beta*p_0*(1.D0-p_0)
  mu1 = 0.D0
  do i = 1,HDIM
    mu1 = mu1 + P1UP(i,i)
  enddo

  !2) SPIN DOWN
  h_0 = EDOWN                   ! Diagonal Hamiltonian H0 respresented in the eigenbasis Q
  p_0 = 0.5D0 + cnst*(h_0-mu0)  ! Initialization for P0 represented in eigenbasis Q

  CALL DGEMM('T', 'N', HDIM, HDIM, HDIM, ONE, QDOWN, HDIM, H1DOWN, HDIM, ZERO, X, HDIM)
  CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, X,     HDIM, QDOWN,  HDIM, ZERO, Y, HDIM)

  P1DOWN = -cnst*Y  

  do i = 1,m  ! Loop over m recursion steps
    p_02 = p_0*p_0
    do j = 1,HDIM
      do k = 1,HDIM
        DX1(k,j) = p_0(k)*P1DOWN(k,j) + P1DOWN(k,j)*p_0(j)
      enddo
    enddo
    iD0 = 1.D0/(2.D0*(p_02-p_0)+1.D0)
    p_0 = iD0*p_02

    do j = 1,HDIM
      do k = 1,HDIM
        P1DOWN(k,j) = iD0(k)*(DX1(k,j) + 2.D0*(P1DOWN(k,j)-DX1(k,j))*p_0(j))
      enddo
    enddo
  enddo


  dPdmu_DOWN = beta*p_0*(1.D0-p_0)

  do i = 1,HDIM
    mu1 = mu1 + P1DOWN(i,i)
  enddo

  mu1 = -mu1/(SUM(dPdmu) + SUM(dPdmu_DOWN))

  do i = 1,HDIM
    ! Trace correction by adding (dP/dmu)*(dmu/dH1) to dP/dH1
    P1UP(i,i)   = P1UP(i,i)   + mu1*dPdmu(i) 
    P1DOWN(i,i) = P1DOWN(i,i) + mu1*dPdmu_DOWN(i)  
  enddo

  CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, QUP, HDIM, P1UP, HDIM, ZERO, X,     HDIM)
  CALL DGEMM('N', 'T', HDIM, HDIM, HDIM, ONE, X,   HDIM, QUP,   HDIM, ZERO, RHOUP, HDIM)

  CALL DGEMM('N', 'N', HDIM,HDIM,HDIM,ONE,QDOWN,HDIM,P1DOWN, HDIM, ZERO, X,       HDIM)
  CALL DGEMM('N', 'T', HDIM,HDIM,HDIM,ONE,X,    HDIM,QDOWN,   HDIM, ZERO, RHODOWN, HDIM)

  DEALLOCATE(P1UP, P1DOWN, X, DX1, Y, H_0, P_0, DPDMU, DPDMU_DOWN, P_02, ID0)

end subroutine Canon_DM_PRT_SPIN
