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

SUBROUTINE SP2PURE_SUBGRAPH_PARALLEL(MDITER)

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
  USE SUBGRAPH
  USE SPARSESP2
  USE SUBGRAPHSP2
  USE XBOARRAY
  USE HOMOLUMO
  USE MATRIXIO

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: MDITER
  INTEGER :: M, ITERZ
  INTEGER, ALLOCATABLE :: IIG(:), JJG(:,:), IIH(:), JJH(:,:)
  INTEGER :: NNZSTART, NNZEND

  REAL(LATTEPREC) :: TRACE
  REAL(LATTEPREC), ALLOCATABLE :: GGN(:,:), HHN(:,:)
  IF (EXISTERROR) RETURN

  ! First time through use current Hamiltonian and threshold
  IF (FIRST_STEP.EQ.1) THEN
     WRITE(*,*) '  '
     WRITE(*,*) '# WARNING PARTITIONED SUBGRAPH TECHNIQUE NOT YET EVAULATED, EXPERIMENTAL VERISON !!!!! '
     WRITE(*,*) '  '

     WRITE(*,*) '# SUBGRAPH FIRST '
     G = BO
     CALL THRESHOLDGRAPH
     TRACE = TRACEGRAPH(HDIM)
     !write(*,*) "graph trace = ", TRACE
     FIRST_STEP = 0
  ENDIF

  ! Calculate M
  M = NNZSTART(MSPARSE, HDIM)

  ! build the sparse H matrix and compute the bounds using
  ! Gershgorin circles
  CALL GERSHGORIN

  ! Allocate for parallelism
  ALLOCATE(IIG(HDIM),IIH(HDIM))
  ALLOCATE(JJG(M,HDIM), JJH(M,HDIM))
  ALLOCATE(GGN(M,HDIM), HHN(M,HDIM))

  ! Put the normalized H into compressed sparse row format
  TX = START_TIMER(SP2ALL_TIMER)
  TX = START_TIMER(SP2SPARSE_TIMER)

  ! Convert dense H matrix to sparse
  TX = START_TIMER(DENSE2SPARSE_TIMER)
  CALL DENSE2SPARSE(H, HDIM, IIH, JJH, HHN)
  CALL DENSE2SPARSEGRAPH(IIG, JJG, GGN)
  TX = STOP_TIMER(DENSE2SPARSE_TIMER)

  ! Start timer for density matrix build
  TX = START_TIMER(DMBUILD_TIMER)

  ! Calculate SP2 sequence of X^2 and 2X-X^2
  CALL SP2SEQUENCE

  ! Partition nodes into subgraphs
  CALL PARTITIONGRAPH

  ! Do progress loop
  CALL PROGRESSLOOP(IIH, JJH, HHN, IIG, JJG)

  ! Density matrix becomes new graph for next iteration
  ! Threshold graph
  G = BO
  CALL THRESHOLDGRAPH

  ! Calculate norm for each iteration
  CALL TRNORMGRAPH(ITERZ)

  ! End of density matrix build
  TX = STOP_TIMER(DMBUILD_TIMER)

  ! Calculate HOMO-LUMO gap estimate
  CALL HOMOLUMOGAP(ITERZ)

  TX = STOP_TIMER(SP2SPARSE_TIMER)
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
  DEALLOCATE(IIG,IIH)
  DEALLOCATE(JJG,JJH)
  DEALLOCATE(HHN,GGN)

  DEALLOCATE(NR_OF_NODES_IN_PART)
  DEALLOCATE(NODE_IN_PART)
  DEALLOCATE(VVX)

  RETURN

END SUBROUTINE SP2PURE_SUBGRAPH_PARALLEL
