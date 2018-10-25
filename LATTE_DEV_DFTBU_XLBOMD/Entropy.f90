subroutine Entropy(S,D0,HDIM)

implicit none
integer, parameter        :: PREC = 8
integer, intent(in)       :: HDIM
real(PREC), intent(in)    :: D0(HDIM,HDIM)
real(PREC), intent(out)   :: S
real(PREC), parameter     :: kB = 8.61739e-5, eps = 1e-14, ZERO = 0.D0
real(PREC)                :: p(HDIM), X(HDIM,HDIM), p_i
integer                   :: i, N

 call Eigx(D0,X,p,'N',HDIM)

 N = HDIM
 S = ZERO
 do i = 1,N
   p_i = p(i)
   if ((p_i > eps).and.((1.D0-p_i) > eps)) then
     S = S - kB*(p_i*log(p_i) + (1.D0-p_i)*log(1.D0-p_i))
   endif
 enddo

end subroutine Entropy
