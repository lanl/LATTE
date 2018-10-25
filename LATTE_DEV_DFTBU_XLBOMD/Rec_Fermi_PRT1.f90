subroutine  Rec_Fermi_PRT1(P0,P1,H0,H1,T,mu0,mu1,Ne,m,eps,MaxIt,HDIM)

integer, parameter      :: PREC = 8
integer, intent(in)     :: HDIM, Ne
real(PREC), parameter   :: ONE = 1.D0, ZERO = 0.D0, TWO = 2.D0
real(PREC), intent(in)  :: H0(HDIM,HDIM), H1(HDIM,HDIM)
real(PREC)              :: X(HDIM,HDIM), Y(HDIM,HDIM)
real(PREC)              :: P02(HDIM,HDIM), II(HDIM,HDIM), XI(HDIM,HDIM)
real(PREC)              :: DX1(HDIM,HDIM)
real(PREC)              :: ID0(HDIM,HDIM)
real(PREC), intent(out) :: P0(HDIM,HDIM), P1(HDIM,HDIM)
real(PREC), intent(in)  :: T, eps
real(PREC), intent(inout)  :: mu0, mu1
real(PREC)              ::  OccErr, TrdPdmu, TrP0, TrP1
real(PREC)              :: beta, cnst, kB
integer, intent(in)     :: m, MaxIt
integer                 :: N, Cnt, i, j, k,l

  mu1 = ZERO  ! Intial guess
  P0 = ZERO
  P1 = ZERO
  N = HDIM
  II = ZERO
  do j = 1,N
    II(j,j) = ONE
  enddo
  OccErr = ONE
  Cnt = ZERO
  kB = 8.61739d-5
  beta = 1.D0/(kB*T)     ! Temp in Kelvin
  cnst = beta/(2.D0**(2+m))
  do while (OccErr > eps)
    P0 = (ONE/TWO)*II - cnst*(H0-mu0*II)
    P1 = -cnst*(H1-mu1*II)
    do i = 1,m
      call MMult(ONE,P0,P0,ZERO,P02,'N','N',HDIM)
      call MMult(ONE,P0,P1,ZERO,X,'N','N',HDIM)
      do k = 1,N
      do l = 1,N
         DX1(k,l) = X(k,l) + X(l,k)
      enddo
      enddo
      X = TWO*(P02-P0) + II
      call Invert(X,ID0,HDIM)
      call MMult(ONE,ID0,P02,ZERO,P0,'N','N',HDIM)

      Y = TWO*(P1-DX1)
      call MMult(ONE,Y,P0,ZERO,X,'N','N',HDIM)
      Y = DX1+X
      call MMult(ONE,ID0,Y,ZERO,P1,'N','N',HDIM)
    enddo
    TrdPdmu = ZERO
    TrP1 = ZERO
    do j = 1,HDIM
      TrdPdmu = TrdPdmu + P0(j,j)
      TrP1 = TrP1 + P1(j,j)
    enddo
    TrP0 = TrdPdmu
    do i = 1,HDIM
    do j = 1,HDIM
      TrdPdmu = TrdPdmu - P0(i,j)**2
    enddo
    enddo
    TrdPdmu = beta*TrdPdmu
    if ((Ne - TrP0)/(TrdPdmu) > 1) then
      mu0 = mu0 + 1.D0
    elseif ((Ne - TrP0)/(TrdPdmu) < -1) then
      mu0 = mu0 - 1.D0
    else
      mu0 = mu0 + (Ne - TrP0)/(TrdPdmu)
    endif
    mu1 = mu1 + (0.D0 - TrP1)/(TrdPdmu)
    OccErr = abs(TrP0 + TrP1-Ne)
    Cnt = Cnt + 1
    if (Cnt > MaxIt) then
      OccErr = 0.D0
    endif
  enddo
  X = II-P0
  call MMult(ONE,P0,X,ZERO,Y,'N','N',HDIM)
  P0 = P0 + ((Ne - TrP0)/TrdPdmu)*beta*Y
  P1 = P1 - (TrP1/N)*II

end subroutine Rec_Fermi_PRT1
