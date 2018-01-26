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

MODULE SUBGRAPH

  USE PUREARRAY
  USE MYPRECISION

  IMPLICIT NONE
  SAVE

  INTEGER :: NR_NODES, NR_PART, FIRST_STEP
  INTEGER, ALLOCATABLE :: NODE_IN_PART(:,:), NR_OF_NODES_IN_PART(:)
  REAL(LATTEPREC), ALLOCATABLE :: G(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: VVX(:,:)
  REAL(LATTEPREC), PARAMETER :: GEPS = 1.0E-3

CONTAINS

  SUBROUTINE THRESHOLDGRAPH

    USE MYPRECISION
    USE CONSTANTS_MOD
    USE OMP_LIB

    IMPLICIT NONE
    INTEGER :: I, J

    !$OMP PARALLEL DO PRIVATE(I,J)
    DO J = 1,HDIM
       DO I = 1,HDIM
          IF (ABS(G(I,J)).LT.GEPS) G(I,J) = ZERO
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE THRESHOLDGRAPH


  SUBROUTINE DENSE2SPARSEGRAPH(IIG, JJG, GGN)

    USE MYPRECISION
    USE CONSTANTS_MOD

    IMPLICIT NONE
    INTEGER :: I, J, K, L
    INTEGER, INTENT(INOUT) :: IIG(:), JJG(:,:)

    REAL(LATTEPREC) :: XTMP
    REAL(LATTEPREC), INTENT(INOUT) :: GGN(:,:)

    !$OMP PARALLEL DO PRIVATE(I,J,K,L,XTMP) 
    DO I = 1,HDIM
       L = 0
       DO J = 1,HDIM
          XTMP = G(J,I)
          IF (J .EQ. I) THEN
             L = L + 1
             JJG(L,I) = J
             GGN(L,I) = XTMP
          ELSEIF (ABS(XTMP) .GE. GEPS) THEN
             L = L + 1
             JJG(L,I) = J
             GGN(L,I) = XTMP
          ENDIF
       ENDDO
       IIG(I) = L
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE DENSE2SPARSEGRAPH


  SUBROUTINE PARTITIONGRAPH

    USE MYPRECISION
    USE CONSTANTS_MOD
    USE SPARSEARRAY
    USE OMP_LIB

    IMPLICIT NONE
    INTEGER :: IT, I, J

    ! Partition graph
    ! Example liquid methane number of nodes = 8 = number of orbitals/molecule
    NR_NODES = 1 !!! Uniform partitioning for simplicity
    NR_PART = HDIM/NR_NODES
    !write(*,*) "number parts = ", NR_PART, " number nodes = ", NR_NODES

    ALLOCATE(NR_OF_NODES_IN_PART(NR_PART))
    ALLOCATE(NODE_IN_PART(NR_PART,NR_NODES))
    ALLOCATE(VVX(100,NR_PART))
    VVX = ZERO

!!! CREATE SIMPLE UNIFORM SUB-GRAPH PARTITIONING OF THE NODES
    IT = 0
    DO I = 1, NR_PART
       DO J = 1, NR_NODES
          IT = IT + 1
          NODE_IN_PART(I,J) = IT
       ENDDO
       NR_OF_NODES_IN_PART(I) = NR_NODES
    ENDDO

  END SUBROUTINE PARTITIONGRAPH


  SUBROUTINE TRNORMGRAPH(ITERZ)

    USE MYPRECISION
    USE CONSTANTS_MOD
    USE SPARSEARRAY

    IMPLICIT NONE
    INTEGER :: I, J
    INTEGER, INTENT(INOUT) :: ITERZ

    REAL(LATTEPREC) :: TRACE

    VV = ZERO

!!! STORE SQRT(Tr[(Tr(X-X^2)^2]) FOR HOMO-LUMO ESTIMATE
    DO J = 1, NR_SP2_ITER
       TRACE = ZERO
       DO I = 1, NR_PART
          TRACE = TRACE + VVX(J,I)
       ENDDO
       VV(J) = SQRT(TRACE)
    ENDDO
    ITERZ = NR_SP2_ITER

  END SUBROUTINE TRNORMGRAPH


  REAL(LATTEPREC) FUNCTION TRACEGRAPH(HDIM)

    USE MYPRECISION

    INTEGER :: I
    INTEGER, INTENT(IN) :: HDIM

    REAL(LATTEPREC) :: TRACE

    TRACE = 0
    DO I = 1,HDIM
       TRACE = TRACE + G(I,I)
    ENDDO

    TRACEGRAPH = TRACE

    RETURN

  END FUNCTION TRACEGRAPH

END MODULE SUBGRAPH
