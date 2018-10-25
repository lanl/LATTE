subroutine RecONFermi(P0,XI0,H0,T,mu0,Nocc,m,eps,MaxIt,HDIM)

implicit none

integer, parameter        :: PREC = 8
integer, intent(in)       :: HDIM, Nocc , MaxIt, m
real(PREC), parameter     :: ONE = 1.D0, ZERO = 0.D0, TWO = 2.D0
real(PREC), parameter     :: kB = 8.61739e-5 ! eV/K, kB = 6.33366256e-6 Ry/K, kB = 3.166811429e-6 Ha/K
real(PREC), intent(in)    :: H0(HDIM,HDIM)
real(PREC), intent(out)   :: P0(HDIM,HDIM)
real(PREC)                :: X(HDIM,HDIM), P02(HDIM,HDIM), II(HDIM,HDIM)
real(PREC), intent(in)    :: T, eps
real(PREC), intent(inout) :: mu0, XI0(HDIM,HDIM)
real(PREC)                :: XI(HDIM,HDIM), Y(HDIM,HDIM)
real(PREC)                :: OccErr, TrdPdmu, TrP0, Xtmp(HDIM,HDIM), III(HDIM,HDIM)
real(PREC)                :: beta, cnst, EBnd
integer                   :: N, Cnt, i, j, NrSchulz

  mu0 = 0.D0  
  N = HDIM
  II = ZERO 
  do j = 1, HDIM
    II(j,j) = ONE
  enddo
  OccErr = 1
  Cnt = 0
  do while (OccErr > eps)
    beta = 1/(kB*T)    ! Temp in Kelvin
    cnst = beta/(2**(2+m))
    P0 = 0.5D0*II - cnst*(H0-mu0*II)
    XI = XI0
    do j = 1,m
      call MMult(ONE,P0,P0,ZERO,P02,'N','T',HDIM)
      Y = 2*(P02-P0)+II
      if (abs(XI(1,1)) < eps) then  ! The first an only time do the full inverse or divide and conquer followed by Schulz It.
        call Invert(Y,XI,HDIM)
!        write(*,*) ' FULL INVERSION '
      else
        NrSchulz = 4  ! Could be reduced to 1 if XI is propagated in XL-BOMD
        do i = 1, NrSchulz
          Xtmp = XI
          call MMult(-ONE,Xtmp,Y,ZERO,X,'N','N',HDIM)
          call MMult(ONE,X,Xtmp,TWO,XI,'N','N',HDIM)
        enddo
     endif
!     call Inv(Y,XI,HDIM)  ! To check the exact inversion
      if (j == 1) then
        XI0 = XI
      endif
      call MMult(ONE,XI,P02,ZERO,P0,'N','N',HDIM)
    enddo
    TrdPdmu = ZERO
    do j = 1,HDIM
      TrdPdmu = TrdPdmu + P0(j,j)
    enddo
    TrP0 = TrdPdmu
    do i = 1,HDIM
    do j = 1,HDIM
      TrdPdmu = TrdPdmu - P0(i,j)**2
    enddo
    enddo
    TrdPdmu = beta*TrdPdmu
    mu0 = mu0 + (Nocc - TrP0)/(TrdPdmu)
    OccErr = abs(TrP0-Nocc)
    Cnt = Cnt + 1
    if (Cnt.ge.MaxIt) then
      OccErr = 0.D0
    endif
  enddo
  ! Adjust occupation
  X = II-P0
!  call MMult(ONE,P0,X,ZERO,XI,'N','N',HDIM)
!  P0 = P0 + ((Nocc - TrP0)/TrdPdmu)*beta*XI
  call MMult(ONE,P0,H0,ZERO,X,'N','N',HDIM)
  EBnd = 0.D0
  do I = 1,HDIM
    EBnd = EBnd + X(I,I)
  enddo
  write(*,*) ' RecONFermi EBnd = ', EBnd
  write(*,*) ' RecONFermi m = ', m
  

end subroutine RecONFermi
