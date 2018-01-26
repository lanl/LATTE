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

MODULE SPARSEMATH

  USE MYPRECISION

  IMPLICIT NONE
  SAVE

  INTEGER, ALLOCATABLE :: IX(:), JJB(:,:)

  REAL(LATTEPREC), ALLOCATABLE :: X(:), Y(:)

CONTAINS

  !
  ! sparsex2 - Sparse matrix multiply  X^2
  !
  SUBROUTINE SPARSEX2(TTRX, TTRX2, II, JJ, VAL, II2, JJ2, VAL2)

    USE MYPRECISION
    USE CONSTANTS_MOD
    USE SPARSEARRAY
    USE OMP_LIB

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: II(:), JJ(:,:)
    INTEGER, INTENT(INOUT) :: II2(:), JJ2(:,:)
    INTEGER :: I, J, K, L, LL, JP, KP

    REAL(LATTEPREC) :: AX, XTMP
    REAL(LATTEPREC), INTENT(INOUT) :: TTRX, TTRX2
    REAL(LATTEPREC), INTENT(IN) :: VAL(:,:)
    REAL(LATTEPREC), INTENT(INOUT) :: VAL2(:,:)

    X = ZERO
    Y = ZERO
    IX = 0
    JJB = 0

    !$OMP PARALLEL DO FIRSTPRIVATE(IX,X) PRIVATE(I,J,L,LL,XTMP,AX,K,KP,JP) &
    !$OMP REDUCTION(+:TTRX,TTRX2)
    DO I = 1,HDIM
       L = 0
       DO JP = 1,II(I)
          AX = VAL(JP,I)
          J = JJ(JP,I)
          IF (J .EQ. I) THEN
             TTRX = TTRX + AX
          ENDIF
          DO KP = 1, II(J)
             K = JJ(KP,J)
             IF(IX(K) .EQ. 0) THEN
                L = L + 1
                JJ2(L,I) = K
                IX(K) = I
                X(K) = ZERO
             ENDIF
             X(K) = X(K) + AX*VAL(KP,J) ! TEMPORARY STORAGE VECTOR LENGTH FULL N
          ENDDO
       ENDDO

       ! Check if over number of non-zeroes limit
       IF (L > MSPARSE) THEN
          WRITE(6,*) "ERROR: Number of non-zeroes per row =",L,"> MSPARSE, Increase MSPARSE."
          CALL ERRORS("sparsemath","Number of non-zeroes per row > MSPARSE")
       ENDIF

       LL = 0
       DO J = 1, L
          JP = JJ2(J,I)
          XTMP = X(JP)
          IF (JP .EQ. I) THEN
             TTRX2 = TTRX2 + XTMP
             LL = LL + 1
             VAL2(LL,I) = XTMP
             JJ2(LL,I) = JP
          ELSEIF(ABS(XTMP) .GT. NUMTHRESH) THEN
             LL = LL + 1
             VAL2(LL,I) = XTMP
             JJ2(LL,I) = JP
          ENDIF
          IX(JP) = 0
          X(JP) = ZERO
       ENDDO
       II2(I) = LL
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE SPARSEX2

  !
  ! sparseadd - Sparse matrix add X = 2X - X^2
  !
  SUBROUTINE SPARSEADD(TTRNORM, II, JJ, VAL, II2, JJ2, VAL2)

    USE MYPRECISION
    USE CONSTANTS_MOD
    USE SPARSEARRAY
    USE OMP_LIB

    IMPLICIT NONE
    INTEGER :: I, J, K, L, LL, JP, KP
    INTEGER, INTENT(IN) :: II2(:), JJ2(:,:)
    INTEGER, INTENT(INOUT) :: II(:), JJ(:,:)

    REAL(LATTEPREC) :: XTMP
    REAL(LATTEPREC), INTENT(INOUT) :: TTRNORM
    REAL(LATTEPREC), INTENT(IN) :: VAL2(:,:)
    REAL(LATTEPREC), INTENT(INOUT) :: VAL(:,:)

    X = ZERO
    Y = ZERO
    IX = 0
    JJB = 0

    TTRNORM = ZERO

    !$OMP PARALLEL DO FIRSTPRIVATE(IX,X,Y)  PRIVATE(I,L,LL,XTMP,K,JP) &
    !$OMP REDUCTION(+:TTRNORM)

    DO I = 1,HDIM
       L = 0
       DO JP = 1, II(I)
          K = JJ(JP,I)
          IF (IX(K) .EQ. 0) THEN
             X(K) = ZERO
             IX(K) = I
             Y(K) = ZERO
             L = L + 1
             JJB(L,I) = K
          ENDIF
          X(K) = X(K) + TWO*VAL(JP,I)
          Y(K) = Y(K) + VAL(JP,I)
       ENDDO

       DO JP = 1, II2(I)
          K = JJ2(JP,I)
          IF (IX(K) .EQ. 0) THEN
             X(K) = ZERO
             IX(K) = I
             Y(K) = ZERO
             L = L + 1
             JJB(L,I) = K
          ENDIF
          X(K) = X(K) - VAL2(JP,I)
          Y(K) = Y(K) - VAL2(JP,I)
       ENDDO

       II(I) = L
       LL = 0
       DO JP = 1,L
          XTMP = X(JJB(JP,I))
          TTRNORM = TTRNORM + Y(JJB(JP,I))*Y(JJB(JP,I))
          IF (ABS(XTMP) .GT. NUMTHRESH) THEN ! THIS THRESHOLDING COULD BE IGNORED!?
             LL = LL + 1
             VAL(LL,I) = XTMP
             JJ(LL,I) = JJB(JP,I)
          ENDIF
          X(JJB(JP,I)) = ZERO
          IX(JJB(JP,I)) = 0
          Y(JJB(JP,I)) = ZERO
          JJB(JP,I) = 0
       ENDDO
       II(I) = LL
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE SPARSEADD

  !
  ! sparsesetx2 - Sparse matrix set to X^2
  !
  SUBROUTINE SPARSESETX2(TTRNORM, II, JJ, VAL, II2, JJ2, VAL2)

    USE MYPRECISION
    USE CONSTANTS_MOD
    USE SPARSEARRAY
    USE OMP_LIB

    IMPLICIT NONE
    INTEGER :: I, J, K, L, LL, JP, KP
    INTEGER, INTENT(INOUT) :: II(:), JJ(:,:)
    INTEGER, INTENT(IN) :: II2(:), JJ2(:,:)

    REAL(LATTEPREC) :: XTMP
    REAL(LATTEPREC), INTENT(INOUT) :: TTRNORM
    REAL(LATTEPREC), INTENT(IN) :: VAL2(:,:)
    REAL(LATTEPREC), INTENT(INOUT) :: VAL(:,:)

    X = ZERO
    Y = ZERO
    IX = 0
    JJB = 0

    TTRNORM = ZERO

    !$OMP PARALLEL DO FIRSTPRIVATE(IX,X,Y) PRIVATE(I,L,LL,K,XTMP,JP)&
    !$OMP REDUCTION(+:TTRNORM)
    DO I = 1,HDIM ! X = X^2

       L = 0
       DO JP = 1, II(I)
          K = JJ(JP,I)
          IF (IX(K) .EQ. 0) THEN
             X(K) = ZERO
             IX(K) = I
             Y(K) = ZERO
             L = L + 1
             JJB(L,I) = K
          ENDIF
          Y(K) = Y(K) + VAL(JP,I)
       ENDDO

       DO JP = 1, II2(I)
          K = JJ2(JP,I)
          IF (IX(K) .EQ. 0) THEN
             X(K) = ZERO
             IX(K) = I
             Y(K) = ZERO
             L = L + 1
             JJB(L,I) = K
          ENDIF
          X(K) = X(K) + VAL2(JP,I)
          Y(K) = Y(K) - VAL2(JP,I)
       ENDDO

       LL = 0
       DO JP = 1, L
          XTMP = X(JJB(JP,I))
          TTRNORM = TTRNORM + Y(JJB(JP,I))*Y(JJB(JP,I))
          IF (ABS(XTMP) .GT. NUMTHRESH) THEN ! THIS THRESHOLDING COULD BE IGNORED!?
             LL = LL + 1
             VAL(LL,I) = XTMP
             JJ(LL,I) = JJB(JP,I)
          ENDIF
          X(JJB(JP,I)) = ZERO
          IX(JJB(JP,I)) = 0
          Y(JJB(JP,I)) = ZERO
          JJB(JP,I) = 0
       ENDDO
       II(I) = LL
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE SPARSESETX2

END MODULE SPARSEMATH
