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

MODULE SPARSESP2

  IMPLICIT NONE

CONTAINS

  SUBROUTINE DENSE2SPARSE(HARRAY, HSIZE, II, JJ, VAL)

    USE MYPRECISION
    USE CONSTANTS_MOD
    USE SPARSEARRAY
    USE PUREARRAY
    USE NONOARRAY

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: HSIZE
    INTEGER, INTENT(INOUT) :: II(:), JJ(:,:)
    INTEGER :: I, J, L, K, BANDWIDTH

    REAL(LATTEPREC), INTENT(IN) :: HARRAY(:,:)
    REAL(LATTEPREC), INTENT(INOUT) :: VAL(:,:)
    REAL(LATTEPREC) :: XTMP

    VAL = ZERO
    JJ = 0
    II = 0

    BANDWIDTH = 0

    IF (BASISTYPE .EQ. "ORTHO") THEN

       ! New compressed row format

       !$OMP PARALLEL DO PRIVATE(I,J,K,L,XTMP) REDUCTION(MAX:BANDWIDTH)
       DO I = 1,HSIZE
          L = 0
          DO J = 1,HSIZE
             XTMP = HARRAY(J,I)
             IF (J .EQ. I) THEN
                L = L + 1
                JJ(L,I) = J
                VAL(L,I) = (MAXEVAL - XTMP)/MAXMINUSMIN
             ELSEIF (ABS(XTMP) .GE. HTHRESH) THEN
                L = L + 1
                JJ(L,I) = J
                VAL(L,I) = -XTMP/MAXMINUSMIN
                BANDWIDTH = MAX(BANDWIDTH, ABS(I-J))
             ENDIF
          ENDDO
          II(I) = L
       ENDDO
       !$OMP END PARALLEL DO

    ELSE

       !$OMP PARALLEL DO PRIVATE(I,J,L,XTMP) REDUCTION(MAX:BANDWIDTH)
       DO I = 1,HSIZE
          L = 0
          DO J = 1,HSIZE
             XTMP = ORTHOH(J,I)
             IF (J .EQ. I) THEN
                L = L + 1
                JJ(L,I) = J
                VAL(L,I) = (MAXEVAL - XTMP)/MAXMINUSMIN
             ELSEIF (ABS(XTMP) .GE. HTHRESH) THEN
                L = L + 1
                JJ(L,I) = J
                VAL(L,I) = -XTMP/MAXMINUSMIN
                BANDWIDTH = MAX(BANDWIDTH, ABS(I-J))
             ENDIF
          ENDDO
          II(I) = L
       ENDDO
       !$OMP END PARALLEL DO

    ENDIF

  END SUBROUTINE DENSE2SPARSE


  SUBROUTINE SPARSE2DENSE(SCALAR, II, JJ, VAL, DARRAY, HSIZE)

    USE MYPRECISION
    USE CONSTANTS_MOD

    IMPLICIT NONE
    INTEGER :: I,J
    INTEGER, INTENT(IN):: HSIZE
    INTEGER, INTENT(IN):: II(:), JJ(:,:)

    REAL(LATTEPREC), INTENT(IN) :: SCALAR
    REAL(LATTEPREC), INTENT(IN) :: VAL(:,:)
    REAL(LATTEPREC), INTENT(INOUT) :: DARRAY(:,:)

    DARRAY = ZERO

    DO I = 1, HSIZE
       DO J = 1,II(I)
          DARRAY(I,JJ(J,I)) = SCALAR*VAL(J,I)
       ENDDO
    ENDDO

  END SUBROUTINE SPARSE2DENSE

  SUBROUTINE SP2LOOP(MSIZE, ITER, II, JJ, VAL)

    USE MYPRECISION
    USE CONSTANTS_MOD
    USE SPARSEARRAY
    USE PUREARRAY
    USE SPARSEMATH
    USE OMP_LIB
    USE MATRIXIO

    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: MSIZE, ITER
    INTEGER, INTENT(INOUT) :: II(:), JJ(:,:)
    INTEGER :: BREAKLOOP, MAXM
    INTEGER, ALLOCATABLE :: IIC(:), JJC(:,:)

    REAL(LATTEPREC) :: TRX, TRNORM, OCC
    REAL(LATTEPREC) :: TRX2, TR2XX2, TRXOLD
    REAL(LATTEPREC) :: LIMDIFF, IDEMPERR, IDEMPERR1, IDEMPERR2
    REAL(LATTEPREC), INTENT(INOUT) :: VAL(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: CCN(:,:)
    REAL(LATTEPREC), PARAMETER :: IDEMTOL = 1.0E-14

    ! Allocate temporary arrays
    ALLOCATE(IX(HDIM))
    ALLOCATE(X(HDIM), Y(HDIM))
    ALLOCATE(IIC(HDIM))
    ALLOCATE(JJC(MSIZE,HDIM), JJB(MSIZE,HDIM))
    ALLOCATE(CCN(MSIZE,HDIM))

    OCC = BNDFIL*FLOAT(HDIM)

    IDEMPERR2 = ZERO
    IDEMPERR1 = ZERO
    IDEMPERR = ZERO

    MAXM = 0

    IX = 0
    X = 0
    Y = 0

    CCN = 0
    JJC = 0
    JJB = 0
    IIC = 0

    ITER = 0
    BREAKLOOP = 0

    DO WHILE ( BREAKLOOP .EQ. 0 .AND. ITER .LT. 100 )

       ITER = ITER + 1

       !
       ! Sparse BO * BO
       !

       TRX = ZERO
       TRX2 = ZERO

       ! Matrix multiply X^2
       CALL SPARSEX2(TRX, TRX2, II, JJ, VAL, IIC, JJC, CCN)

       MAXM = MAX(MAXM, 2*MAXVAL(IIC))

       TR2XX2 = TWO*TRX - TRX2
       TRXOLD = TRX
       LIMDIFF = ABS(TRX2 - OCC) - ABS(TR2XX2 - OCC)

       IF ( LIMDIFF .GT. IDEMTOL ) THEN ! X <= 2X-X^2

          TRX = TWO * TRX - TRX2

          PP(ITER) = 0

          ! Sparse matrix X = 2X - X^2
          CALL SPARSEADD(TRNORM, II, JJ, VAL, IIC, JJC, CCN)

       ELSEIF ( LIMDIFF .LT. -IDEMTOL ) THEN ! X <= X^2

          TRX = TRX2
          PP(ITER) = 1;

          ! Sparse matrix X = X^2
          CALL SPARSESETX2(TRNORM, II, JJ, VAL, IIC, JJC, CCN)

       ELSE

          TRX = TRXOLD
          BREAKLOOP = 1

       ENDIF

       VV(ITER) = SQRT(TRNORM)
       !WRITE(*,*) 'IT = ', ITER, ' VV(',ITER,',) = ', VV(ITER), PP(ITER)

       IDEMPERR2 = IDEMPERR1
       IDEMPERR1 = IDEMPERR
       IDEMPERR = ABS(TRX - TRXOLD)

       IF (SP2CONV .EQ. "REL" .AND.  ITER .GE. MINSP2ITER &
            .AND. (IDEMPERR .GE. IDEMPERR2) ) BREAKLOOP = 1

       IF (SP2CONV .EQ. "ABS" .AND. ABS(LIMDIFF) .LE. IDEMTOL) BREAKLOOP = 1

       ! NOT converging
       IF (ITER .EQ. 100) THEN
          CALL PANIC
          CALL ERRORS("sparsesp2","Sparse SP2 purification is not converging: STOP!")
       ENDIF

    ENDDO

    NR_SP2_ITER = ITER

    ! Reset number of non-zeroes per row - MSPARSE
    MSIZE = MAXM

    !  MSIZE = HDIM
    ! Deallocate temporary arrays
    DEALLOCATE(IX)
    DEALLOCATE(X, Y)
    DEALLOCATE(IIC)
    DEALLOCATE(JJC, JJB)
    DEALLOCATE(CCN)

  END SUBROUTINE SP2LOOP


  SUBROUTINE SP2SEQUENCELOOP(MSIZE, ITER, II, JJ, VAL)

    USE MYPRECISION
    USE CONSTANTS_MOD
    USE PUREARRAY
    USE SPARSEARRAY
    USE SPARSEMATH
    USE OMP_LIB

    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: MSIZE, ITER
    INTEGER :: BREAKLOOP, MAXM
    INTEGER, INTENT(INOUT) :: II(:), JJ(:,:)
    INTEGER, ALLOCATABLE :: IIC(:), JJC(:,:)

    REAL(LATTEPREC) :: TRX, TRNORM, OCC
    REAL(LATTEPREC) :: TRX2, TR2XX2, TRXOLD
    REAL(LATTEPREC) :: LIMDIFF, IDEMPERR, IDEMPERR1, IDEMPERR2
    REAL(LATTEPREC), INTENT(INOUT) :: VAL(:,:)
    REAL(LATTEPREC), ALLOCATABLE :: CCN(:,:)
    REAL(LATTEPREC), PARAMETER :: IDEMTOL = 1.0E-14

    ! Allocate temporary arrays
    ALLOCATE(IX(HDIM))
    ALLOCATE(X(HDIM), Y(HDIM))
    ALLOCATE(IIC(HDIM))
    ALLOCATE(JJC(MSIZE,HDIM), JJB(MSIZE,HDIM))
    ALLOCATE(CCN(MSIZE,HDIM))

    OCC = BNDFIL*FLOAT(HDIM)

    IDEMPERR2 = 0.0
    IDEMPERR1 = 0.0
    IDEMPERR = 0.0

    MAXM = 0

    IX = 0
    X = ZERO
    Y = ZERO

    CCN = ZERO
    JJC = 0
    IIC = 0

    ITER = 0
    BREAKLOOP = 0
    !  WRITE(*,*) ' NR_SP2_ITER = ',NR_SP2_ITER

    DO WHILE (ITER .LT. NR_SP2_ITER )

       ITER = ITER + 1

       !
       ! Sparse BO * BO
       !

       TRX = ZERO
       TRX2 = ZERO

       ! Matrix multiply X^2
       CALL SPARSEX2(TRX, TRX2, II, JJ, VAL, IIC, JJC, CCN)

       MAXM = MAX(MAXM, 2*MAXVAL(IIC))


       TR2XX2 = TWO*TRX - TRX2
       TRXOLD = TRX
       LIMDIFF = ABS(TRX2 - OCC) - ABS(TR2XX2 - OCC)

       IF (PP(ITER).EQ.0) THEN ! X <= 2X-X^2
          TRX = TWO * TRX - TRX2

          ! Sparse matrix X = 2X - X^2
          CALL SPARSEADD(TRNORM, II, JJ, VAL, IIC, JJC, CCN)

       ELSE ! X <= X^2

          TRX = TRX2

          ! Sparse matrix X = X^2
          CALL SPARSESETX2(TRNORM, II, JJ, VAL, IIC, JJC, CCN)

       ENDIF

       VV(ITER) = SQRT(TRNORM)
       !    WRITE(*,*) 'IT = ', ITER, ' VV(',ITER,',) = ', VV(ITER), PP(ITER)

       ! NOT converging
       IF (ITER .GE. 100) THEN
          CALL PANIC
          CALL ERRORS("sparsesp2","Sparse SP2 purification is not converging: STOP!")
       ENDIF

    ENDDO
    NR_SP2_ITER = ITER

    ! Reset number of non-zeroes per row - MSPARSE

    !  PRINT*, MAXM, HDIM, 2*MAXVAL(IIC)

    MSIZE = MAXM
    !  MSIZE = HDIM
    ! Deallocate temporary arrays
    DEALLOCATE(IX)
    DEALLOCATE(X, Y)
    DEALLOCATE(IIC)
    DEALLOCATE(JJC, JJB)
    DEALLOCATE(CCN)

  END SUBROUTINE SP2SEQUENCELOOP

END MODULE SPARSESP2
