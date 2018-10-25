!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 2010.  Los Alamos National Security, LLC. This material was    !
! produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos !
! National Laboratory (LANL), which is operated by Los Alamos National     !
! Security, LLC for the U.S. Department of Energy. The U.S. Government has !
! rights to use, reproduce, and distribute this software.  NEITHER THE     !
! GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY,     !
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS         !
! SOFTWARE.  If software is modified to produce derivative works, such     !
! modified software should be clearly marked, so as not to confuse it      !
! with the version available from LANL.                                    !
!                                                                          !
! Additionally, this program is free software; you can redistribute it     !
! and/or modify it under the terms of the GNU General Public License as    !
! published by the Free Software Foundation; version 2.0 of the License.   !
! Accordingly, this program is distributed in the hope that it will be     !
! useful, but WITHOUT ANY WARRANTY; without even the implied warranty of   !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General !
! Public License for more details.                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE HUBBARDFORCE  ! AND ENERGY

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE COULOMBARRAY
  USE NONOARRAY
  USE KSPACEARRAY
  USE MYPRECISION
  USE DMARRAY

  IMPLICIT NONE

  INTEGER :: I, J, K
  INTEGER :: INDEX, NUMORB
  REAL(LATTEPREC) :: ES, EP, ED, EF
  REAL(LATTEPREC) :: HMOD, X(HDIM,HDIM), Y(HDIM,HDIM), HU(HDIM,HDIM), H_U(HDIM,HDIM), D(HDIM,HDIM)
  REAL(LATTEPREC) :: UD(HDIM,HDIM), UP(HDIM,HDIM), DS(HDIM,HDIM), P(HDIM,HDIM)
  REAL(LATTEPREC) :: DSUP(HDIM,HDIM), UPSD(HDIM,HDIM), PSUD(HDIM,HDIM), UDSP(HDIM,HDIM), PSUP(HDIM,HDIM)
  REAL(LATTEPREC) :: UPSP(HDIM,HDIM), A(HDIM,HDIM), B(HDIM,HDIM)
  COMPLEX(LATTEPREC) :: ZHMOD
  IF (EXISTERROR) RETURN

  D = 0.5D0*BO
  P = 0.5D0*BO
  call DGEMM('T','N',HDIM,HDIM,HDIM,ONE,BO,HDIM,SMAT,HDIM,ZERO,PS,HDIM)
  do I = 1, HDIM
  do J = 1, HDIM
     UD(I,J) = DFTB_U(I)*D(I,J)
     UP(I,J) = DFTB_U(I)*P(I,J)
  enddo
  enddo

!call MMult(ONE,U,D,ZERO,UD,'N','N',HDIM)
!call MMult(ONE,U,P,ZERO,UP,'N','N',HDIM)

!call MMult(ONE,P,S,ZERO,PS,'N','N',HDIM)
     call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,P,HDIM,SMAT,HDIM,ZERO,PS,HDIM)

!call MMult(ONE,D,S,ZERO,DS,'N','N',HDIM)
     call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,D,HDIM,SMAT,HDIM,ZERO,DS,HDIM)

!call MMult(ONE,DS,UP,ZERO,DSUP,'N','N',HDIM)
     call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,DS,HDIM,UP,HDIM,ZERO,DSUP,HDIM)

!call MMult(ONE,UP,DS,ZERO,UPSD,'N','T',HDIM)
     call DGEMM('N','T',HDIM,HDIM,HDIM,ONE,UP,HDIM,DS,HDIM,ZERO,UPSD,HDIM)

!call MMult(ONE,PS,UD,ZERO,PSUD,'N','N',HDIM)
     call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,PS,HDIM,UD,HDIM,ZERO,PSUD,HDIM)

!call MMult(ONE,UD,PS,ZERO,UDSP,'N','T',HDIM)
     call DGEMM('N','T',HDIM,HDIM,HDIM,ONE,UD,HDIM,PS,HDIM,ZERO,UDSP,HDIM)

!call MMult(ONE,PS,UP,ZERO,PSUP,'N','N',HDIM)
     call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,PS,HDIM,UP,HDIM,ZERO,PSUP,HDIM)

!call MMult(ONE,UP,PS,ZERO,UPSP,'N','T',HDIM)
     call DGEMM('N','T',HDIM,HDIM,HDIM,ONE,UP,HDIM,PS,HDIM,ZERO,UPSP,HDIM)

  A = UD - UPSD - UDSP + UPSP
!  call MMult(ONE,A,SMAT,ZERO,B,'N','N',HDIM)
     call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,A,HDIM,SMAT,HDIM,ZERO,B,HDIM)
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

!!!! = (1/2)*Sum_I { U_I*( D(I,K)*SR(K,I) + SR(I,:)*D(:,I)*d_KI
!!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,K,I_A,I_B,Xtmp,Ytmp,Ztmp)
!  do I = 1,Nr_atoms
!    Xtmp = ZERO
!    Ytmp = ZERO
!    Ztmp = ZERO
!    I_A = H_INDEX_START(I)
!    I_B = H_INDEX_END(I)
!    do II = I_A,I_B
!       do L = 1,HDIM
!          Xtmp = Xtmp + dSx(II,L)*A(L,II)
!          Ytmp = Ytmp + dSy(II,L)*A(L,II)
!          Ztmp = Ztmp + dSz(II,L)*A(L,II)
!       enddo
!    enddo
!    !HubForce(1:3,I) = -0.25D0*[Xtmp,Ytmp,Ztmp]
!    HubForce(1:3,I) = -0.5D0*[Xtmp,Ytmp,Ztmp]
!  enddo
!!$OMP END PARALLEL DO


  RETURN

END SUBROUTINE HUBBARDFORCE
