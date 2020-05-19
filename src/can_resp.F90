subroutine  can_resp(H1,Nocc,beta,Q,ev,fe,mu0,eps,HDIM)

USE SETUPARRAY
implicit none

integer, parameter      :: PREC = 8
integer, intent(in)     :: HDIM, Nocc
!real(PREC), parameter   :: ONE = 1.D0, ZERO = 0.D0, TWO = 2.D0
real(PREC), intent(in)  :: H1(HDIM,HDIM), Q(HDIM,HDIM), ev(HDIM), fe(HDIM)
real(PREC)              :: X(HDIM,HDIM), YY(HDIM,HDIM), fermifun
real(PREC)              :: P02(HDIM,HDIM), II(HDIM,HDIM), XI(HDIM,HDIM)
real(PREC)              :: DX1(HDIM,HDIM)
real(PREC)              :: ID0(HDIM,HDIM), T12(HDIM,HDIM), DDT(HDIM,HDIM)
!real(PREC), intent(out) :: P1(HDIM,HDIM)
real(PREC), intent(in)  :: beta, eps
real(PREC), intent(inout)  :: mu0
real(PREC)              ::  OccErr, TrdPdmu, TrP0, TrP1
real(PREC)              :: trX, trDDT
real(PREC)              :: mu1, y
integer                 :: N, Cnt, i, j, k,l, ind, ind2, fact

  N = HDIM
  CALL DGEMM('T', 'N', HDIM, HDIM, HDIM, ONE, &
          Q, HDIM, H1, HDIM, ZERO, X, HDIM)
  CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
          X, HDIM, Q, HDIM, ZERO, T12, HDIM)
!  call MMult(ONE,Q,H1,ZERO,X,'T','N',HDIM)  
!  call MMult(ONE,X,Q,ZERO,T12,'N','N',HDIM)

  DDT = ZERO
  X = ZERO
  do i = 1,N
     do j = 1,N
        if (abs(ev(i)-ev(j)) < 1e-4) then
           y = (ev(i)+ev(j))/2.D0
           if (abs(beta*(y-mu0)) > 500.D0) then
              DDT(i,j) = 0.D0
           else
              DDT(i,j) = -beta*exp(beta*(y-mu0))/(exp(beta*(y-mu0))+1.D0)**2
           endif
        else
           DDT(i,j) = (fe(i)-fe(j))/(ev(i)-ev(j)) 
        endif
     end do
  end do
  do j = 1,N
    do i = 1,N
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
  CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
          Q, HDIM, X, HDIM, ZERO, YY, HDIM)
  CALL DGEMM('N', 'T', HDIM, HDIM, HDIM, ONE, &
          YY, HDIM, Q, HDIM, ZERO, T12, HDIM)
  BO = T12

!  call MMult(ONE,Q,X,ZERO,YY,'N','N',HDIM)
!  call MMult(ONE,YY,Q,ZERO,P1,'N','T',HDIM)

end subroutine can_resp
