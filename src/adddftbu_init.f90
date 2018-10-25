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

SUBROUTINE ADDDFTBU_INIT

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
  REAL(LATTEPREC) :: HMOD, X(HDIM,HDIM), Y(HDIM,HDIM), HU(HDIM,HDIM), H_U(HDIM,HDIM)
  REAL(LATTEPREC) :: DS(HDIM,HDIM)
  COMPLEX(LATTEPREC) :: ZHMOD
  IF (EXISTERROR) RETURN

  !!! BRUTE FORCE INTRODUCTION OF HUBBARD U FOR THE DM SHIFT OPERATOR ONLY FOR SP MATERIALS
  C_DFTB_U = 0.D0
  do I = 1, NATS   ! Orbital dependent U
     DFTB_U(H_INDEX_START(I)) = 0.D0  ! s-orbitals
     do J = H_INDEX_START(I)+1, H_INDEX_END(I)
        DFTB_U(J) = C_DFTB_U    !  p-orbitals Hubbard U-J
     enddo
  enddo

  do I = 1, HDIM
  do J = 1, HDIM
     SU(I,J) = SMAT(I,J)*DFTB_U(J)
  enddo
  enddo
  call DGEMM('T','N',HDIM,HDIM,HDIM,ONE,BO,HDIM,SMAT,HDIM,ZERO,DS,HDIM)
  call DGEMM('T','N',HDIM,HDIM,HDIM,ONE,DS,HDIM,SU,HDIM,ZERO,X,HDIM)
  call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,SU,HDIM,DS,HDIM,ZERO,Y,HDIM)
  HU = 0.125D0*(SU-X-Y) ! since BO = 2*D, otherwise (1/4)*(...)
  do I = 1,HDIM
  do J = 1,HDIM
     H_U(I,J) = HU(I,J) + HU(J,I)
  enddo
  enddo
  H = H + H_U

  IF (KON .EQ. 0) THEN
     IF ( BASISTYPE .EQ. "ORTHO") THEN
          WRITE(*,*) ' DM HUBBARD U ONLY FOR NONORTHOGONAL REPRESENTATIONS '
     ENDIF
  ELSE
     WRITE(*,*) ' DM HUBBARD U ONLY FOR KON = 0 '
  ENDIF

  RETURN

END SUBROUTINE ADDDFTBU_INIT
