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

SUBROUTINE SP2PURE_SPARSE_PARALLEL_SIMPLE(MDITER)

  !
  ! PARALLEL SHARED MEMORY OPENMP VERSION
  ! This subroutine implements Niklasson's SP2 density matrix purification
  ! algorithm.
  !

  USE CONSTANTS_MOD
  USE TIMER_MOD
  USE SETUPARRAY
  USE PUREARRAY
  USE SPARSEARRAY
  USE NONOARRAY
  USE MYPRECISION
  USE MATRIXIO
  USE SPARSESP2
  USE HOMOLUMO

  IMPLICIT NONE

  INTEGER :: M, ITERZ, MDITER
  INTEGER :: NNZSTART, NNZEND
  INTEGER, ALLOCATABLE :: IIA(:), JJA(:,:)

  REAL(LATTEPREC), ALLOCATABLE :: AAN(:,:)
  IF (EXISTERROR) RETURN

  ! Calculate M - max number of non-zeroes per row
  M = NNZSTART(MSPARSE, HDIM)

  ! build the sparse H matrix and compute the bounds using
  ! Gershgorin circles
  CALL GERSHGORIN

  ! Allocate for parallelism
  ALLOCATE(IIA(HDIM))
  ALLOCATE(JJA(M,HDIM))
  ALLOCATE(AAN(M,HDIM))

  ! Put the normalized H into compressed sparse row format
  TX = START_TIMER(SP2ALL_TIMER)
  TX = START_TIMER(SP2SPARSE_TIMER)

  ! Convert dense H matrix to sparse
  TX = START_TIMER(DENSE2SPARSE_TIMER)
  CALL DENSE2SPARSE(H, HDIM, IIA, JJA, AAN)
  TX = STOP_TIMER(DENSE2SPARSE_TIMER)

  ! Calculate SP2 sequence of X^2 and 2X-X^2
  ! based on HOMO-LUMO gap
  CALL SP2SEQUENCE

  ! Do SP2 algorithm using calculated sequence
  TX = START_TIMER(DMBUILD_TIMER)
  CALL SP2SEQUENCELOOP(M, ITERZ, IIA, JJA, AAN)
  TX = STOP_TIMER(DMBUILD_TIMER)

  ! Calculate HOMO-LUMO gap estimate
  CALL HOMOLUMOGAP(ITERZ)

  TX = STOP_TIMER(SP2SPARSE_TIMER)

  ! Convert sparse to dense
  TX = START_TIMER(SPARSE2DENSE_TIMER)
  CALL SPARSE2DENSE(TWO, IIA, JJA, AAN, BO, HDIM)
  TX = STOP_TIMER(SPARSE2DENSE_TIMER)

  TX = STOP_TIMER(SP2ALL_TIMER)

  ! Reset number of non-zeroes
  MSPARSE = NNZEND(M, HDIM)

  ! Write out input and output matrix
  IF (DEBUGON .GE. 2) THEN
     CALL WRITEHMATRIX(HDIM, MSPARSE, H, NR_SP2_ITER, PP)
  ENDIF
  IF (DEBUGON .EQ. 3) THEN
     CALL WRITEDMATRIX(HDIM, BO)
  ENDIF

  ! Deallocate for sparse format
  DEALLOCATE(IIA)
  DEALLOCATE(JJA)
  DEALLOCATE(AAN)

  RETURN

END SUBROUTINE SP2PURE_SPARSE_PARALLEL_SIMPLE
