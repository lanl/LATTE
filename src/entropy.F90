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

SUBROUTINE ENTROPY

  !
  ! This subroutine compute the electronic entropy if we're running
  ! with a finite electron temperature
  !

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE SPINARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J
  INTEGER :: LIWORK, LWORK, INFO
  INTEGER, ALLOCATABLE :: IWORK(:)
  REAL(LATTEPREC), ALLOCATABLE :: WORK(:), S_EVEC(:,:), S_EVAL(:)
  REAL(LATTEPREC), ALLOCATABLE :: MAT2(:,:), Y(:,:)
  REAL(LATTEPREC) :: OCC, S, NPOWER
  REAL(LATTEPREC) :: OCCLOGOCC_HOLES, OCCLOGOCC_ELECTRONS
  REAL(LATTEPREC) :: CLOSE2ZERO = 1.0D-16
  REAL(LATTEPREC) :: ELECTRON, HOLE, TMPELEC, TMPHOLE
  REAL(LATTEPREC) :: LN2, C4(2), C8(4)
  CHARACTER(LEN=1), PARAMETER :: JOBZ = "N",  UPLO = "U"

  IF (EXISTERROR) RETURN

  !
  ! Options:
  !
  ! ENTROPYKIND = 0: S=0 (for testing purposes)
  !
  ! Sanville's algorithm:
  !
  ! ENTROPYKIND = 1 : S = XlnX + (1-X)ln(1-X) (exact for Fermi-Dirac smearing)
  !
  ! THE NEXT 3 OPTIONS WERE INVENTED BY ANDERS NIKLASSON:
  !
  ! ENTROPYKIND = 2 : An option for use with running finite T SP2
  !
  ! N = 2^NORECS
  !
  ! S = 2N*[X(X^(1/N) -1)/(X^(1/N)+1) + (1-X)*((1-X)^(1/N)-1)/((1-X)^(1/N)+1)]
  !
  ! ENTROPYKIND = 3: 4th-order expansion: approximate but no diagonalization
  !
  ! Y = X(X-1)
  ! S = Tr(Y(C1 + C2*Y))
  !
  ! ENTROPYKIND = 4: As above, but with an 8th-order expansion
  !
  ! Y = X(X-1)
  ! S = Tr[C1*Y + Y^2*(C2*I + C3*Y + c4*Y^2)]
  !

  S = ZERO

  IF (ENTROPYKIND .EQ. 0) THEN

     S = ZERO

  ELSEIF (ENTROPYKIND .EQ. 1) THEN

     ALLOCATE(S_EVEC(HDIM, HDIM), S_EVAL(HDIM))

#ifdef XSYEV
     LWORK = 3*HDIM - 1
     ALLOCATE(WORK(LWORK))
#elif defined(XSYEVD)
     LWORK = 1 + 6*HDIM + 2*HDIM*HDIM
     LIWORK = 3 + 5*HDIM
     ALLOCATE(WORK(LWORK), IWORK(LIWORK))
#endif

     IF (SPINON .EQ. 0) THEN

        S_EVEC = HALF*BO

#ifdef XSYEV
        ! Pre-processing to select precision

#ifdef DOUBLEPREC
        CALL DSYEV(JOBZ, UPLO, HDIM, S_EVEC, HDIM, S_EVAL, WORK, LWORK,INFO)
#elif defined(SINGLEPREC)
        CALL SSYEV(JOBZ, UPLO, HDIM, S_EVEC, HDIM, S_EVAL, WORK, LWORK,INFO)
#endif

#elif defined(XSYEVD)

#ifdef DOUBLEPREC
        CALL DSYEVD(JOBZ, UPLO, HDIM, S_EVEC, HDIM, S_EVAL, WORK, LWORK, &
             IWORK, LIWORK, INFO)
#elif defined(SINGLEPREC)
        CALL SSYEVD(JOBZ, UPLO, HDIM, S_EVEC, HDIM, S_EVAL, WORK, LWORK, &
             IWORK, LIWORK, INFO)
#endif

#endif

        DO I = 1, HDIM

           OCC = S_EVAL(I)

           IF (OCC .LE. CLOSE2ZERO .OR. OCC .GE. (ONE - CLOSE2ZERO)) THEN

              OCCLOGOCC_ELECTRONS = ZERO
              OCCLOGOCC_HOLES = ZERO

           ELSE

              OCCLOGOCC_ELECTRONS = OCC * LOG(OCC)
              OCCLOGOCC_HOLES = (ONE - OCC) * LOG(ONE - OCC)

           ENDIF

           S = S + TWO*(OCCLOGOCC_ELECTRONS + OCCLOGOCC_HOLES)

        ENDDO

     ELSE

        !
        ! Do both density matrices separately - this is correct
        ! when we're doing spin-polarized calculations
        !

        S_EVEC = RHOUP

#ifdef XSYEV
        ! Pre-processing to select precision

#ifdef DOUBLEPREC
        CALL DSYEV(JOBZ, UPLO, HDIM, S_EVEC, HDIM, S_EVAL, WORK, LWORK,INFO)
#elif defined(SINGLEPREC)
        CALL SSYEV(JOBZ, UPLO, HDIM, S_EVEC, HDIM, S_EVAL, WORK, LWORK,INFO)
#endif

#elif defined(XSYEVD)

#ifdef DOUBLEPREC
        CALL DSYEVD(JOBZ, UPLO, HDIM, S_EVEC, HDIM, S_EVAL, WORK, LWORK, &
             IWORK, LIWORK, INFO)
#elif defined(SINGLEPREC)
        CALL SSYEVD(JOBZ, UPLO, HDIM, S_EVEC, HDIM, S_EVAL, WORK, LWORK, &
             IWORK, LIWORK, INFO)
#endif

#endif


        DO I = 1, HDIM

           OCC = S_EVAL(I)

           IF (OCC .LE. CLOSE2ZERO .OR. OCC .GE. (ONE - CLOSE2ZERO)) THEN

              OCCLOGOCC_ELECTRONS = ZERO
              OCCLOGOCC_HOLES = ZERO

           ELSE

              OCCLOGOCC_ELECTRONS = OCC * LOG(OCC)
              OCCLOGOCC_HOLES = (ONE - OCC) * LOG(ONE - OCC)

           ENDIF

           S = S + OCCLOGOCC_ELECTRONS + OCCLOGOCC_HOLES

        ENDDO

        S_EVEC = RHODOWN

#ifdef XSYEV
        ! Pre-processing to select precision

#ifdef DOUBLEPREC
        CALL DSYEV(JOBZ, UPLO, HDIM, S_EVEC, HDIM, S_EVAL, WORK, LWORK,INFO)
#elif defined(SINGLEPREC)
        CALL SSYEV(JOBZ, UPLO, HDIM, S_EVEC, HDIM, S_EVAL, WORK, LWORK,INFO)
#endif

#elif defined(XSYEVD)

#ifdef DOUBLEPREC
        CALL DSYEVD(JOBZ, UPLO, HDIM, S_EVEC, HDIM, S_EVAL, WORK, LWORK, &
             IWORK, LIWORK, INFO)
#elif defined(SINGLEPREC)
        CALL SSYEVD(JOBZ, UPLO, HDIM, S_EVEC, HDIM, S_EVAL, WORK, LWORK, &
             IWORK, LIWORK, INFO)
#endif

#endif


        DO I = 1, HDIM

           OCC = S_EVAL(I)

           IF (OCC .LE. CLOSE2ZERO .OR. OCC .GE. (ONE - CLOSE2ZERO)) THEN

              OCCLOGOCC_ELECTRONS = ZERO
              OCCLOGOCC_HOLES = ZERO

           ELSE

              OCCLOGOCC_ELECTRONS = OCC * LOG(OCC)
              OCCLOGOCC_HOLES = (ONE - OCC) * LOG(ONE - OCC)

           ENDIF

           S = S + OCCLOGOCC_ELECTRONS + OCCLOGOCC_HOLES

        ENDDO

     ENDIF

     DEALLOCATE(S_EVEC, S_EVAL)

#ifdef XSYEV
     DEALLOCATE(WORK)
#elif defined(XSYEVD)
     DEALLOCATE(WORK, IWORK)
#endif

  ELSEIF (ENTROPYKIND .EQ. 2) THEN

     IF (CONTROL .NE. 5) THEN
        CALL ERRORS("entropy","Only use ENTROPYKIND = 2 if you're using CONTROL = 5")
     ENDIF

     LWORK = 3*HDIM - 1

     ALLOCATE(S_EVEC(HDIM, HDIM), S_EVAL(HDIM), WORK(LWORK))

     NPOWER = TWO**NORECS

     IF (SPINON .EQ. 0) THEN

        S_EVEC = HALF*BO

#ifdef DOUBLEPREC
        CALL DSYEV(JOBZ, UPLO, HDIM, S_EVEC, HDIM, S_EVAL, WORK, LWORK,INFO)
#elif defined(SINGLEPREC)
        CALL SSYEV(JOBZ, UPLO, HDIM, S_EVEC, HDIM, S_EVAL, WORK, LWORK,INFO)
#endif

        DO I = 1, HDIM

           OCC = S_EVAL(I)

           IF (OCC .LE. CLOSE2ZERO) THEN
              ELECTRON = ZERO
           ENDIF

           IF (ONE - OCC .LE. CLOSE2ZERO) THEN
              HOLE = ZERO
           ENDIF

           IF (OCC .GT. CLOSE2ZERO .AND. &
                ONE - OCC .GT. CLOSE2ZERO) THEN

              TMPELEC = OCC**(ONE/NPOWER)
              TMPHOLE = (ONE - OCC)**(ONE/NPOWER)

              ELECTRON = OCC*(TMPELEC - ONE)/(TMPELEC + ONE)
              HOLE = (ONE - OCC)*(TMPHOLE - ONE)/(TMPHOLE + ONE)

           ENDIF

           S = S + FOUR*NPOWER*(ELECTRON + HOLE)

        ENDDO

     ELSE

        !
        ! Spin polarized
        !

        S_EVEC = RHOUP

#ifdef DOUBLEPREC
        CALL DSYEV(JOBZ, UPLO, HDIM, S_EVEC, HDIM, S_EVAL, WORK, LWORK,INFO)
#elif defined(SINGLEPREC)
        CALL SSYEV(JOBZ, UPLO, HDIM, S_EVEC, HDIM, S_EVAL, WORK, LWORK,INFO)
#endif

        DO I = 1, HDIM

           OCC = S_EVAL(I)

           IF (OCC .LE. CLOSE2ZERO) THEN
              ELECTRON = ZERO
           ENDIF

           IF (ONE - OCC .LE. CLOSE2ZERO) THEN
              HOLE = ZERO
           ENDIF

           IF (OCC .GT. CLOSE2ZERO .AND. &
                ONE - OCC .GT. CLOSE2ZERO) THEN

              TMPELEC = OCC**(ONE/NPOWER)
              TMPHOLE = (ONE - OCC)**(ONE/NPOWER)

              ELECTRON = OCC*(TMPELEC - ONE)/(TMPELEC + ONE)
              HOLE = (ONE - OCC)*(TMPHOLE - ONE)/(TMPHOLE + ONE)

           ENDIF

           S = S + TWO*NPOWER*(ELECTRON + HOLE)

        ENDDO

        S_EVEC = RHODOWN

#ifdef DOUBLEPREC
        CALL DSYEV(JOBZ, UPLO, HDIM, S_EVEC, HDIM, S_EVAL, WORK, LWORK,INFO)
#elif defined(SINGLEPREC)
        CALL SSYEV(JOBZ, UPLO, HDIM, S_EVEC, HDIM, S_EVAL, WORK, LWORK,INFO)
#endif

        DO I = 1, HDIM

           OCC = S_EVAL(I)

           IF (OCC .LE. CLOSE2ZERO) THEN
              ELECTRON = ZERO
           ENDIF

           IF (ONE - OCC .LE. CLOSE2ZERO) THEN
              HOLE = ZERO
           ENDIF

           IF (OCC .GT. CLOSE2ZERO .AND. &
                ONE - OCC .GT. CLOSE2ZERO) THEN

              TMPELEC = OCC**(ONE/NPOWER)
              TMPHOLE = (ONE - OCC)**(ONE/NPOWER)

              ELECTRON = OCC*(TMPELEC - ONE)/(TMPELEC + ONE)
              HOLE = (ONE - OCC)*(TMPHOLE - ONE)/(TMPHOLE + ONE)

           ENDIF

           S = S + TWO*NPOWER*(ELECTRON + HOLE)

        ENDDO

     ENDIF

     DEALLOCATE(S_EVEC, S_EVAL, WORK)

  ELSEIF (ENTROPYKIND .EQ. 3) THEN

     !
     ! 4th-order approximation
     !

     LN2 = LOG(TWO)

     C4(1) = EIGHT*LN2 - TWO
     C4(2) = SIXTEEN*LN2 - EIGHT

     ALLOCATE(MAT2(HDIM, HDIM))

     IF (SPINON .EQ. 0) THEN

        !
        ! Here we make use of *GEMM to do half*BO*half*BO
        !

#ifdef DOUBLEPREC
        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 0.25D0, &
             BO, HDIM, BO, HDIM, 0.0D0, MAT2, HDIM)
#elif defined(SINGLEPREC)
        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 0.25, &
             BO, HDIM, BO, HDIM, 0.0, MAT2, HDIM)
#endif

        DO I = 1, HDIM
           DO J = I, HDIM

              IF (I .EQ. J) THEN

                 S = S + (MAT2(I,I) - HALF*BO(I,I)) * &
                      (C4(1) + C4(2)*(MAT2(I,I) - HALF*BO(I,I)))

              ELSE

                 S = S + TWO*C4(2) * (MAT2(J,I) - HALF*BO(J,I)) * &
                      (MAT2(J,I) - HALF*BO(J,I))

              ENDIF

           ENDDO
        ENDDO

        S = TWO*S

     ELSE

        !
        ! First spin up:
        !

#ifdef DOUBLEPREC
        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
             RHOUP, HDIM, RHOUP, HDIM, 0.0D0, MAT2, HDIM)
#elif defined(SINGLEPREC)
        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
             RHOUP, HDIM, RHOUP, HDIM, 0.0, MAT2, HDIM)
#endif

        DO I = 1, HDIM
           DO J = I, HDIM

              IF (I .EQ. J) THEN

                 S = S + (MAT2(I,I) - RHOUP(I,I)) * &
                      (C4(1) + C4(2)*(MAT2(I,I) - RHOUP(I,I)))

              ELSE

                 S = S + TWO*C4(2) * (MAT2(J,I) - RHOUP(J,I)) * &
                      (MAT2(J,I) - RHOUP(J,I))

              ENDIF

           ENDDO
        ENDDO

        !
        ! Spin down
        !

#ifdef DOUBLEPREC
        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
             RHODOWN, HDIM, RHODOWN, HDIM, 0.0D0, MAT2, HDIM)
#elif defined(SINGLEPREC)
        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
             RHODOWN, HDIM, RHODOWN, HDIM, 0.0, MAT2, HDIM)
#endif

        DO I = 1, HDIM
           DO J = I, HDIM

              IF (I .EQ. J) THEN

                 S = S + (MAT2(I,I) - RHODOWN(I,I)) * &
                      (C4(1) + C4(2)*(MAT2(I,I) - RHODOWN(I,I)))

              ELSE

                 S = S + TWO*C4(2) * (MAT2(J,I) - RHODOWN(J,I)) * &
                      (MAT2(J,I) - RHODOWN(J,I))

              ENDIF

           ENDDO
        ENDDO

     ENDIF

     DEALLOCATE(MAT2)

  ELSEIF (ENTROPYKIND .EQ. 4) THEN

     !
     ! 8th order expansion
     !

     LN2 = LOG(TWO)

     C8(1) = 16.0D0*LN2 - (34.0D0/5.0D0)
     C8(2) = 96.0D0*LN2 - (844.0D0/15.0D0)
     C8(3) = 256.0D0*LN2 - (2336.0D0/15.0D0)
     C8(4) = 256.0D0*LN2 - (2368.0D0/15.0D0)

     ALLOCATE(MAT2(HDIM, HDIM), Y(HDIM, HDIM))

     IF (SPINON .EQ. 0) THEN

        !
        ! MAT2 = half*BO*half*BO, temporarily (note we take
        ! care of the 'halfs' in *GEMM
        !

#ifdef DOUBLEPREC
        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 0.25D0, &
             BO, HDIM, BO, HDIM, 0.0D0, MAT2, HDIM)
#elif defined(SINGLEPREC)
        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 0.25, &
             BO, HDIM, BO, HDIM, 0.0, MAT2, HDIM)
#endif

        Y = MAT2 - HALF*BO

        !
        ! MAT2 = Y*Y
        !

#ifdef DOUBLEPREC
        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
             Y, HDIM, Y, HDIM, 0.0D0, MAT2, HDIM)
#elif defined(SINGLEPREC)
        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
             Y, HDIM, Y, HDIM, 0.0, MAT2, HDIM)
#endif

        DO I = 1, HDIM
           DO J = I, HDIM

              IF (I .EQ. J) THEN

                 S = S + C8(1)*Y(I,I) + MAT2(I,I) * &
                      (C8(2) + C8(3)*Y(I,I) + C8(4)*MAT2(I,I))

              ELSE

                 S = S + C8(3)*MAT2(J,I)*Y(J,I) + &
                      C8(4)*MAT2(J,I)*MAT2(J,I)

              ENDIF

           ENDDO
        ENDDO

        S = TWO*S

     ELSE

        !
        ! For the up-spins:
        !

        ! MAT2 = X*X

#ifdef DOUBLEPREC
        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
             RHOUP, HDIM, RHOUP, HDIM, 0.0D0, MAT2, HDIM)
#elif defined(SINGLEPREC)
        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
             RHOUP, HDIM, RHOUP, HDIM, 0.0, MAT2, HDIM)
#endif

        Y = MAT2 - RHOUP

        ! MAT2 = Y*Y

#ifdef DOUBLEPREC
        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
             Y, HDIM, Y, HDIM, 0.0D0, MAT2, HDIM)
#elif defined(SINGLEPREC)
        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
             Y, HDIM, Y, HDIM, 0.0, MAT2, HDIM)
#endif

        DO I = 1, HDIM
           DO J = I, HDIM

              IF (I .EQ. J) THEN

                 S = S + C8(1)*Y(I,I) + MAT2(I,I) * &
                      (C8(2) + C8(3)*Y(I,I) + C8(4)*MAT2(I,I))

              ELSE

                 S = S + C8(3)*MAT2(J,I)*Y(J,I) + &
                      C8(4)*MAT2(J,I)*MAT2(J,I)

              ENDIF

           ENDDO
        ENDDO

        !
        ! Now for spin-down
        !

        ! MAT2 = X*X

#ifdef DOUBLEPREC
        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
             RHODOWN, HDIM, RHODOWN, HDIM, 0.0D0, MAT2, HDIM)
#elif defined(SINGLEPREC)
        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
             RHODOWN, HDIM, RHODOWN, HDIM, 0.0, MAT2, HDIM)
#endif

        Y = MAT2 - RHODOWN

        ! MAT2 = Y*Y

#ifdef DOUBLEPREC
        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
             Y, HDIM, Y, HDIM, 0.0D0, MAT2, HDIM)
#elif defined(SINGLEPREC)
        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
             Y, HDIM, Y, HDIM, 0.0, MAT2, HDIM)
#endif

        DO I = 1, HDIM
           DO J = I, HDIM

              IF (I .EQ. J) THEN

                 S = S + C8(1)*Y(I,I) + MAT2(I,I) * &
                      (C8(2) + C8(3)*Y(I,I) + C8(4)*MAT2(I,I))

              ELSE

                 S = S + C8(3)*MAT2(J,I)*Y(J,I) + &
                      C8(4)*MAT2(J,I)*MAT2(J,I)

              ENDIF

           ENDDO
        ENDDO

     ENDIF

     DEALLOCATE(MAT2, Y)

  ENDIF

  ENTE = MINUSONE*KBT*S

  RETURN

END SUBROUTINE ENTROPY
