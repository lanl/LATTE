subroutine Fermi_S(S,HO,D0,QQ,ee,Fe_vec,mu0,H,Z,Nocc,T,OccErrLim,MaxIt,HDIM)

!D0,dq_dv,v,mu0,mu1,T,RX,RY,RZ,LBox,Hubbard_U,Element_Type,Nr_atoms,MaxIt,eps,m,HDIM, &
!Max_Nr_Neigh,Coulomb_acc,TIMERATIO,nnRx,nnRy,nnRz,nrnnlist,nnType,H_INDEX_START,H_INDEX_END,H,S,Z,Nocc,Znuc)
!!!  In  Z,H,T,Nocc,OccErrLim,mu0, HDIM
!!!  Out Q,ee,Fe_vec,mu0

implicit none
integer,    parameter       :: PREC = 8
integer,    intent(in)      :: HDIM, Nocc
real(PREC), parameter       :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0
real(PREC), parameter       :: kB = 8.61739d-5 ! eV/K, kB = 6.33366256e-6 Ry/K, kB = 3.166811429e-6 Ha/K
real(PREC), intent(in)      :: H(HDIM,HDIM), Z(HDIM,HDIM)
real(PREC)                  :: Fe(HDIM,HDIM)
real(PREC), intent(out)     :: D0(HDIM,HDIM),QQ(HDIM,HDIM),ee(HDIM),Fe_vec(HDIM),S,HO(HDIM,HDIM)
real(PREC)                  :: X(HDIM,HDIM)
real(PREC), intent(in)      :: T, OccErrLim
real(PREC), intent(inout)   :: mu0
integer,    intent(in)      :: MaxIt
real(PREC)                  :: OccErr, Occ, dOcc, mu_0, beta, Occ_I, p_i, p(HDIM), EBnd
integer                     :: I,J,ITER

!!
  call MMult(ONE,Z,H,ZERO,X,'T','N',HDIM)  !
  call MMult(ONE,X,Z,ZERO,HO,'N','N',HDIM)   ! HO = Z'*H*Z
!  write(*,*) ' HO_Orth = ', HO(1:5,1)
  call Eigx(HO,QQ,ee,'V',HDIM)
!  write(*,*) ' FERMI_S                 Nocc =', Nocc,'     ee = ', ee(1:HDIM)

  OccErr = ONE
  Fe = ZERO
  beta = 1.D0/(kB*T)    ! Temp in Kelvin
  mu0 = 0.5D0*(ee(Nocc)+ee(Nocc+1)) ! If initial guess is hard ot find this is a good test, but better reuse from previous step
  ITER = 1
  do while (OccErr > OccErrLim)
    ITER = ITER + 1
    Occ = ZERO
    dOcc = ZERO
    do I = 1,HDIM
      Occ_I = ONE/(EXP(beta*(ee(I)-mu0))+ONE)
      Fe_vec(I) = Occ_I
      Occ = Occ + Occ_I
      dOcc = dOcc + beta*Occ_I*(ONE-Occ_I)
    enddo
    OccErr = abs(Nocc-Occ)
    if (abs(OccErr) > 1e-10) then
       mu0 = mu0 + (Nocc-Occ)/dOcc
    endif
    IF(ITER.GT.MaxIt) THEN
      OccErr = ZERO
      write(*,*) ' Warning occupation error in Fermi_S. Not converged within MaxIt '
    ENDIF
  enddo
  if (dOcc > 1e-12) then
    do I = 1,HDIM
      Fe_vec(I) = Fe_vec(I) + ((Nocc-Occ)/dOcc)*beta*Fe_vec(I)*(1.D0-Fe_vec(I))
      Fe(I,I) = Fe_vec(I)
    enddo
  else
    do I = 1,HDIM
      Fe_vec(I) = Fe_vec(I)
      Fe(I,I) = Fe_vec(I)
    enddo
  endif
  call MMult(ONE,QQ,Fe,ZERO,X,'N','N',HDIM)
  call MMult(ONE,X,QQ,ZERO,D0,'N','T',HDIM)

 S = ZERO
 do i = HDIM,1,-1
   p_i = Fe_vec(i)
   if ((p_i > 1e-14).and.((1.D0-p_i) > 1e-14)) then
     S = S - kB*(p_i*log(p_i) + (1.D0-p_i)*log(1.D0-p_i))
   endif
 enddo

!  call MMult(ONE,D0,H0,ZERO,X, 'N', 'N', HDIM)  
!  EBnd = 0.D0
!  do I = 1,HDIM
!    EBnd = EBnd + X(I,I)
!  enddo
!  write(*,*) ' Fermi_S EBnd = ', EBnd

!!! Round Off !!!
! do I = 1,HDIM
! do J = I,HDIM
!    D0(I,J) = floor(D0(I,J)*1000000.D0)/1000000.D0
!    D0(J,I) = D0(I,J)
! enddo
! enddo

end subroutine Fermi_S
