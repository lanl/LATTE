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

MODULE SUBGRAPHSP2

  USE MYPRECISION

  IMPLICIT NONE
  SAVE

  INTEGER, ALLOCATABLE :: IX(:), JJN(:), JJP(:), LG(:)

  REAL(LATTEPREC), ALLOCATABLE :: XS(:,:)

CONTAINS

  SUBROUTINE PROGRESSLOOP(IIH, JJH, HHN, IIG, JJG)

    USE MYPRECISION
    USE CONSTANTS_MOD
    USE SPARSEARRAY
    USE SETUPARRAY
    USE SUBGRAPH
    USE OMP_LIB

    IMPLICIT NONE
    INTEGER :: I, L, LL
    INTEGER, INTENT(IN) :: IIH(:), JJH(:,:), IIG(:), JJG(:,:)

    REAL(LATTEPREC) :: OCC
    REAL(LATTEPREC), INTENT(IN) :: HHN(:,:)

    ALLOCATE(IX(HDIM))
    ALLOCATE(JJN(HDIM),JJP(HDIM))
    ALLOCATE(LG(HDIM))

    OCC = BNDFIL*FLOAT(HDIM)

!!! Progress loop
    IX = 0
    JJN = 0
    JJP = 0
    LG = 0

    !$OMP PARALLEL DO FIRSTPRIVATE(IX,JJP,JJN,LG) &
    !$OMP PRIVATE(I,L,LL,XS)
    DO I = 1, NR_PART

       LL = 0
       L = 0

!!! Extract small dense subgraph Hamiltonians
       CALL EXTRACTSUBGRAPH(I, IIH, JJH, HHN, IIG, JJG, IX, JJN, JJP, LG, L, LL)

       ! Normalize small dense subgraph
       ALLOCATE(XS(L,L))
       CALL NORMALIZESUBGRAPH(XS, JJN, L)

       ! Subgraph sp2 loop on each subgraph 
       CALL SUBGRAPHSP2LOOP(I, XS, JJP, L, LL)

       ! Collect each partial density matrix over subgraphs
       CALL SUBGRAPHSPARSE2DENSE(TWO, L, LL, JJN, JJP, BO, XS)
       DEALLOCATE(XS)

    ENDDO
    !$OMP END PARALLEL DO

    ! Deallocate arrays
    DEALLOCATE(IX)
    DEALLOCATE(JJN,JJP)
    DEALLOCATE(LG)

  END SUBROUTINE PROGRESSLOOP

  SUBROUTINE EXTRACTSUBGRAPH(I, IIH, JJH, HHN, IIG, JJG, IX, JJN, JJP, LG, L, LL)

    USE MYPRECISION
    USE CONSTANTS_MOD
    USE SUBGRAPH

    IMPLICIT NONE

    INTEGER :: J, II, JP, K, L0, LS
    INTEGER, INTENT(IN) :: I
    INTEGER, INTENT(INOUT) :: L, LL
    INTEGER, INTENT(IN) :: IIH(:), JJH(:,:), IIG(:), JJG(:,:)
    INTEGER, INTENT(INOUT) :: IX(:), JJN(:), JJP(:), LG(:)

    REAL(LATTEPREC) :: TRACE
    REAL(LATTEPREC), INTENT(IN) :: HHN(:,:)

    ! Extract elements from graph
    DO J = 1, NR_OF_NODES_IN_PART(I)
       II = NODE_IN_PART(I,J)

       ! NR_OF_NODES(:) = [2,3,3,1,4,...]
       ! NODE_PART(:,:) = [1,2; 3,4,5; 6,9,21; 12; 17,8,11,55; ...]

       ! From graph
       DO JP = 1, IIG(II)
          K = JJG(JP,II)
          IF (IX(K) .NE. I) THEN
             IX(K) = I
             L = L + 1
             JJN(L) = K   ! VECTOR FOR SUBGRAPH EXTRACTION IN a
             LG(K) = L
          ENDIF
          IF (K.EQ.II) THEN
             LL = LL + 1
             JJP(LL) = LG(K) !  ROW INDEX OF ESSENTIAL DM NODES IN SUBGRAPH PARTITIONING
          ENDIF
       ENDDO
    ENDDO
    !if (I .EQ. 1) write(*,*) "From graph 1: L LL = ", L, " ", LL
    !if (I .EQ. 2) write(*,*) "From graph 2: L LL = ", L, " ", LL

    ! Add possible new elements in H (Here the same as before, so no change!
    L0 = L
    DO J = 1, NR_OF_NODES_IN_PART(I)
       II = NODE_IN_PART(I,J)

       ! NR_OF_NODES(:) = [2,3,3,1,4,...]
       ! NODE_PART(:,:) = [1,2; 3,4,5; 6,9,21; 12; 17,8,11,55; ...]

       ! From Hamiltonian
       DO JP = 1, IIH(II)
          K = JJH(JP,II)
          IF (IX(K) .NE. I) THEN
             IX(K) = I
             L = L + 1
             JJN(L) = K   ! VECTOR FOR SUBGRAPH EXTRACTION IN a
          ENDIF
       ENDDO
    ENDDO
    !if(I .EQ. 1)write(*,*) "From H: L LL = ", L, " ", LL

    ! Perform a "double jump" for possible extra elements
    ! based on graph
    LS = L
    DO J = 1,LS
       II = JJN(J)
       DO JP = 1, IIG(II)
          K = JJG(JP,II)
          IF (IX(K) .NE. I) THEN
             IX(K) = I
             L = L + 1
             JJN(L) = K
          ENDIF
       ENDDO
    ENDDO
    !if(I.EQ.1)write(*,*) "From graph double jump: L LL = ", L, " ", LL

    ! Example of automatic partitioning 
    ! Includes 1st orbital or number of nodes plus halo
    IF (I.EQ.1) THEN
       WRITE(*,*) '# SUBGRAPH_1_SIZE: ', L, ' x ', L
    ENDIF
    IF (I.EQ.10) THEN
       WRITE(*,*) '# SUBGRAPH_10_SIZE: ', L, ' x ', L
    ENDIF

  END SUBROUTINE EXTRACTSUBGRAPH

  SUBROUTINE NORMALIZESUBGRAPH(XS, JJN, L)

    USE MYPRECISION
    USE CONSTANTS_MOD
    USE SETUPARRAY

    IMPLICIT NONE
    INTEGER :: J, JA, JB
    INTEGER, INTENT(IN) :: L
    INTEGER, INTENT(IN) ::JJN(:) 

    REAL(LATTEPREC), INTENT(INOUT) :: XS(:,:)

    ! Get values for subgraph from dense H
    DO JA = 1, L
       DO JB = 1, L
          XS(JA,JB) = H(JJN(JA),JJN(JB))  ! (*)
       ENDDO
    ENDDO

    XS = MINUSONE * XS
    DO J = 1,L
       XS(J,J) = XS(J,J) + MAXEVAL
    ENDDO
    XS = (ONE/(MAXEVAL-MINEVAL))*XS

  END SUBROUTINE NORMALIZESUBGRAPH

  SUBROUTINE SUBGRAPHSP2LOOP(I, XS, JJP, L, LL)

    USE MYPRECISION
    USE CONSTANTS_MOD
    USE PUREARRAY
    USE SPARSEARRAY
    USE SUBGRAPH

    IMPLICIT NONE
    INTEGER :: IT, J, JA, JB
    INTEGER, INTENT(IN) :: I, L, LL
    INTEGER, INTENT(IN) :: JJP(:)

    REAL(LATTEPREC) :: TRACE
    REAL(LATTEPREC), INTENT(INOUT) :: XS(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: XST(:,:)


    ! Subgraph sp2 loop on XS
    !
    ! Calculate X-X^2 for first iteration
    ! XST = -1.0 * XS^2 + 1.0 * XST
    ALLOCATE(XST(L,L))
    !!       XST = XS - MATMUL(XS,XS)
    XST = XS
    CALL DGEMM('N', 'N', L, L, L, MINUSONE, &
         XS, L, XS, L, ONE, XST, L)

    ! Calculate trace of (X-X^2)^2 for each iteration
    DO IT = 1,NR_SP2_ITER
       TRACE = ZERO
       DO J = 1,LL
          JA = JJP(J)
          DO JB = 1,L
             TRACE = TRACE + XST(JB,JA)*XST(JB,JA)   ! Tr[(X-X^2)^2]
          ENDDO
          !if (I.EQ.1) write(*,*)"IT trace = ", IT, " ", TRACE
       ENDDO
       VVX(IT,I) = VVX(IT,I) + TRACE  ! ADD IN SHARED MEMORY VVX

       ! Calculate X-X^2 for each iteration
       ! XS = XS +/- XST and XST = XS
       ! XST = -1.0 * XS^2 + 1.0 * XST
       XS = XS + (ONE-TWO*PP(IT))*XST
       !!          XST = XS - MATMUL(XS,XS)
       XST = XS
       CALL DGEMM('N', 'N', L, L, L, MINUSONE, &
            XS, L, XS, L, ONE, XST, L)
    ENDDO

    DEALLOCATE(XST)

  END SUBROUTINE SUBGRAPHSP2LOOP

  SUBROUTINE SUBGRAPHSPARSE2DENSE(SCALAR, L, LL, JJN, JJP, DARRAY, XS)

    USE MYPRECISION

    IMPLICIT NONE
    INTEGER :: JA, JB
    INTEGER, INTENT(IN) :: L, LL
    INTEGER, INTENT(IN) :: JJN(:), JJP(:)

    REAL(LATTEPREC), INTENT(IN) :: SCALAR
    REAL(LATTEPREC), INTENT(INOUT) :: XS(:,:)
    REAL(LATTEPREC), INTENT(INOUT) :: DARRAY(:,:)

    DO JA = 1,LL ! Collect density matrix over nodes for partial density matrix
       !          WRITE(*,*) ' JJP = ', JJP(JA)
       DO JB = 1, L
          DARRAY(JJN(JB),JJN(JJP(JA))) = SCALAR*XS(JB,JJP(JA))
       ENDDO
    ENDDO

  END SUBROUTINE SUBGRAPHSPARSE2DENSE

END MODULE SUBGRAPHSP2
