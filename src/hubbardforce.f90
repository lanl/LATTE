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
! subroutine HubbardForce(HubForce,EHub,Nr_atoms,HDIM,U,D,P,S,dSx,dSy,dSz,H_INDEX_START,H_INDEX_END)
! D*S*D = D, P = Dynamical variable density matrix in non-orthogonal atomic orbital representation
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
  REAL(LATTEPREC) :: HMOD, X(HDIM,HDIM), Y(HDIM,HDIM), HU(HDIM,HDIM), D(HDIM,HDIM)
  REAL(LATTEPREC) :: UD(HDIM,HDIM), UP(HDIM,HDIM), DS(HDIM,HDIM), P(HDIM,HDIM)
  REAL(LATTEPREC) :: DSUP(HDIM,HDIM), UPSD(HDIM,HDIM), PSUD(HDIM,HDIM), UDSP(HDIM,HDIM), PSUP(HDIM,HDIM)
  REAL(LATTEPREC) :: UPSP(HDIM,HDIM), A(HDIM,HDIM), B(HDIM,HDIM)
  REAL(LATTEPREC) :: Xtmp, Ytmp, Ztmp
  COMPLEX(LATTEPREC) :: ZHMOD
  INTEGER :: I_A, I_B, L, II, JJ, J_A, J_B  

  IF (EXISTERROR) RETURN

  D = 0.5D0*BO   ! D*S*D = D, BO*S*BO = 2*BO
  P = 0.5D0*PNO  ! D*S*D = D, SHOULD BE P = D for SCF CONVERGENCE

  call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,P,HDIM,SMAT,HDIM,ZERO,PS,HDIM)

  do I = 1, HDIM
  do J = 1, HDIM
     UD(I,J) = DFTB_U(I)*D(I,J)
     UP(I,J) = DFTB_U(I)*P(I,J)
  enddo
  enddo

  call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,P,HDIM,SMAT,HDIM,ZERO,PS,HDIM)
  call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,D,HDIM,SMAT,HDIM,ZERO,DS,HDIM)
  call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,DS,HDIM,UP,HDIM,ZERO,DSUP,HDIM)
  call DGEMM('N','T',HDIM,HDIM,HDIM,ONE,UP,HDIM,DS,HDIM,ZERO,UPSD,HDIM)
  call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,PS,HDIM,UD,HDIM,ZERO,PSUD,HDIM)
  call DGEMM('N','T',HDIM,HDIM,HDIM,ONE,UD,HDIM,PS,HDIM,ZERO,UDSP,HDIM)
  call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,PS,HDIM,UP,HDIM,ZERO,PSUP,HDIM)
  call DGEMM('N','T',HDIM,HDIM,HDIM,ONE,UP,HDIM,PS,HDIM,ZERO,UPSP,HDIM)
  A = UD - UPSD - UDSP + UPSP

  call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,A,HDIM,SMAT,HDIM,ZERO,B,HDIM)

  EHub = 0.D0 
  do I = 1,HDIM
    EHub = EHub + 0.5D0*B(I,I)  ! Hubbar Ebergy contribution
  enddo
  write(6,*)'EHUB=', EHUB

  B = UD - DSUP - UPSD - PSUD - UDSP + PSUP + UPSP

  do I = 1,HDIM
  do J = 1,HDIM
    A(I,J) = B(I,J) + B(J,I)
  enddo
  enddo
  AHub = A ! For input to tb_forces

  RETURN

END SUBROUTINE HUBBARDFORCE

