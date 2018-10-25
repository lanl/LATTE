subroutine Gershgorin(e1,eN,A,HDIM)

implicit none
integer, parameter        :: PREC = 8
integer      , intent(in) :: HDIM
real(PREC), intent(in)    :: A(HDIM,HDIM)
real(PREC), intent(out)   :: e1,eN
integer                   :: I
real(PREC)                :: R, e_high, e_low

e1 = A(1,1)
eN = e1
do I = 1,HDIM
  R = sum(abs(A(I,:))) - abs(A(I,I))
  e_high = A(I,I) + R
  e_low = A(I,I) - R
  if (e_high > eN) then
    eN = e_high
  endif
  if (e_low < e1) then
    e1 = e_low
  endif
enddo

end subroutine Gershgorin
