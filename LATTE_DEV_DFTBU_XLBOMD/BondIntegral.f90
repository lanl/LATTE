subroutine BondIntegral(hxx,dR,f)

!%%% dR: distance between atoms
!%%% f: paramters/coefficeints for the bond integral

implicit none

integer, parameter :: PREC = 8
real(PREC)         :: f(14), dR, hxx
real(PREC)         :: RMOD, POLYNOM, RMINUSR1, X

  if (dR.le.f(7)) then
    RMOD = dR - f(6)
    POLYNOM = RMOD*(f(2) + RMOD*(f(3) + RMOD*(f(4) + f(5)*RMOD)))
    X = exp(POLYNOM)
  elseif ((dR.gt.f(7)).and.(dR.lt.f(8))) then
    RMINUSR1 = dR - f(7)
    X = f(9) + RMINUSR1*(f(10) + RMINUSR1*(f(11) + RMINUSR1*(f(12) + RMINUSR1*(f(13) + RMINUSR1*f(14)))))
  else
    X = 0
  endif
  hxx = f(1)*X

return
end

