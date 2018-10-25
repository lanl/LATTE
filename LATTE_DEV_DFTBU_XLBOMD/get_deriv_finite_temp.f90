subroutine  get_deriv_finite_temp(P1,H0,H1,Nocc,T,Q,ev,fe,mu0,eps,HDIM)

! call get_deriv_finite_temp(D_dq_J,H0,H1,Nocc,T,Q,E,mu0)
implicit none
integer, parameter      :: PREC = 8
integer, intent(in)       :: HDIM, Nocc
real(PREC), parameter   :: ONE = 1.D0, ZERO = 0.D0, TWO = 2.D0
real(PREC), intent(in)  :: H0(HDIM,HDIM), H1(HDIM,HDIM), Q(HDIM,HDIM), ev(HDIM), fe(HDIM)
real(PREC)              :: X(HDIM,HDIM), YY(HDIM,HDIM), fermifun
real(PREC)              :: P02(HDIM,HDIM), II(HDIM,HDIM), XI(HDIM,HDIM)
real(PREC)              :: DX1(HDIM,HDIM)
real(PREC)              :: ID0(HDIM,HDIM), T12(HDIM,HDIM), DDT(HDIM,HDIM)
real(PREC), intent(out) :: P1(HDIM,HDIM)
real(PREC), intent(in)  :: T, eps
real(PREC), intent(inout)  :: mu0
real(PREC)              ::  OccErr, TrdPdmu, TrP0, TrP1
real(PREC)              :: beta, cnst, kB, mu1, dh, dm, y, dy, trX, trDDT
integer                 :: N, Cnt, i, j, k,l, ind, ind2, fact

real                     :: t1, t2

  N = HDIM
  kB = 8.61739e-5
  beta = 1.D0/(kB*T)     ! Temp in Kelvin

  ! T12 = F0_ort_Vecs'*F1_ort*F0_ort_Vecs;
  call MMult(ONE,Q,H1,ZERO,X,'T','N',HDIM)  
  call MMult(ONE,X,Q,ZERO,T12,'N','N',HDIM)

  ! Divided differences matrix
  DDT = ZERO
  X = ZERO
  do i = 1,N
     do j = 1,N
        if (abs(ev(i)-ev(j)) < 1e-4) then
           ! divided difference from Taylor expansion
           y = (ev(i)+ev(j))/2.D0
           !DDT(i,j) = -beta*exp(beta*(y-mu0))/(exp(beta*(y-mu0))+1.D0)**2
           if (abs(beta*(y-mu0)) > 500.D0) then
              DDT(i,j) = 0.D0
           else
              DDT(i,j) = -beta*exp(beta*(y-mu0))/(exp(beta*(y-mu0))+1.D0)**2
           endif
        else
           DDT(i,j) = (fe(i)-fe(j))/(ev(i)-ev(j)) 
        endif
        X(i,j) = DDT(i,j)*T12(i,j) 
     end do
  end do
  mu1 = ZERO
  trX = ZERO
  trDDT = ZERO
  do i = 1,N
     trX = trX + X(i,i)
     trDDT = trDDT + DDT(i,i)
  end do
  if (abs(trDDT) > 1e-12) then
    mu1 = trX/trDDT
  else
    mu1 = 0.D0
  endif
  do i = 1,N
     X(i,i) = X(i,i)-DDT(i,i)*mu1
  end do
  call MMult(ONE,Q,X,ZERO,YY,'N','N',HDIM)
  call MMult(ONE,YY,Q,ZERO,P1,'N','T',HDIM)

end subroutine get_deriv_finite_temp
