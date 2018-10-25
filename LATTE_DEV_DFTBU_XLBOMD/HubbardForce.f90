subroutine HubbardForce(HubForce,EHub,Nr_atoms,HDIM,U,D,P,S,dSx,dSy,dSz,H_INDEX_START,H_INDEX_END)

use omp_lib
implicit none
integer,    parameter     :: PREC = 8
real(PREC), parameter     :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0, SIX = 6.D0
integer,    intent(in)    :: Nr_atoms, HDIM
real(PREC), intent(in)    :: D(HDIM,HDIM), P(HDIM,HDIM), S(HDIM,HDIM), U(HDIM,HDIM)
integer,    intent(in)    :: H_INDEX_START(Nr_atoms), H_INDEX_END(Nr_atoms)
real(PREC), intent(in)    :: dSx(HDIM,HDIM), dSy(HDIM,HDIM), dSz(HDIM,HDIM)
real(PREC), intent(out)   :: HubForce(3,Nr_atoms), EHub
real(PREC)                :: Xtmp, Ytmp, Ztmp
real(PREC)                :: X(HDIM,HDIM),Y(HDIM,HDIM),Z(HDIM,HDIM),Mtmp(HDIM,HDIM)
real(PREC)                :: UD(HDIM,HDIM),UP(HDIM,HDIM),PS(HDIM,HDIM),DS(HDIM,HDIM)
real(PREC)                :: DSUP(HDIM,HDIM),UPSD(HDIM,HDIM),PSUD(HDIM,HDIM),UDSP(HDIM,HDIM)
real(PREC)                :: PSUP(HDIM,HDIM),UPSP(HDIM,HDIM), A(HDIM,HDIM),B(HDIM,HDIM)
real(PREC)                :: XP(HDIM,HDIM),YP(HDIM,HDIM),ZP(HDIM,HDIM)
integer                   :: I,J,K, I_A, I_B, L, II, JJ, J_A, J_B

!E1 = 0.5*trace(U*D*S - U*P*S*D*S - U*D*S*P*S + U*P*S*P*S)
!UD = U*D; UP = U*P; PS = P*S; DS = D*S;
!DSUP = DS*UP; UPSD = UP*DS'; PSUD = PS*UD; UDSP = UD*PS'; PSUP = PS*UP; UPSP = UP*PS';
!
!B = UD - DSUP - UPSD - PSUD - UDSP + PSUP + UPSP;
!A = B+B';
!F_1 = 0.5*trace(SR1*A)

HubForce = ZERO

call MMult(ONE,U,D,ZERO,UD,'N','N',HDIM)
call MMult(ONE,U,P,ZERO,UP,'N','N',HDIM)
call MMult(ONE,P,S,ZERO,PS,'N','N',HDIM)
call MMult(ONE,D,S,ZERO,DS,'N','N',HDIM)
call MMult(ONE,DS,UP,ZERO,DSUP,'N','N',HDIM)
call MMult(ONE,UP,DS,ZERO,UPSD,'N','T',HDIM)
call MMult(ONE,PS,UD,ZERO,PSUD,'N','N',HDIM)
call MMult(ONE,UD,PS,ZERO,UDSP,'N','T',HDIM)
call MMult(ONE,PS,UP,ZERO,PSUP,'N','N',HDIM)
call MMult(ONE,UP,PS,ZERO,UPSP,'N','T',HDIM)

!E1 = 0.5*trace(U*D*S - U*P*S*D*S - U*D*S*P*S + U*P*S*P*S)
A = UD - UPSD - UDSP + UPSP
call MMult(ONE,A,S,ZERO,B,'N','N',HDIM)
EHub = 0.D0
do I = 1,HDIM
  EHub = EHub + 0.5D0*B(I,I)
enddo

B = UD - DSUP - UPSD - PSUD - UDSP + PSUP + UPSP

do I = 1,HDIM
do J = 1,HDIM
  A(I,J) = B(I,J) + B(J,I)
enddo
enddo

!!! = (1/2)*Sum_I { U_I*( D(I,K)*SR(K,I) + SR(I,:)*D(:,I)*d_KI 
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,K,I_A,I_B,Xtmp,Ytmp,Ztmp)
do I = 1,Nr_atoms 
  Xtmp = ZERO
  Ytmp = ZERO
  Ztmp = ZERO
  I_A = H_INDEX_START(I)
  I_B = H_INDEX_END(I)
  do II = I_A,I_B
     do L = 1,HDIM
        Xtmp = Xtmp + dSx(II,L)*A(L,II)
        Ytmp = Ytmp + dSy(II,L)*A(L,II)
        Ztmp = Ztmp + dSz(II,L)*A(L,II)
     enddo
  enddo
  !HubForce(1:3,I) = -0.25D0*[Xtmp,Ytmp,Ztmp]
  HubForce(1:3,I) = -0.5D0*[Xtmp,Ytmp,Ztmp]
enddo
!$OMP END PARALLEL DO

end subroutine HubbardForce

