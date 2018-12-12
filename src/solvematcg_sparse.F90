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

SUBROUTINE SOLVEMATCG

  USE CONSTANTS_MOD
  USE FERMICOMMON
  USE SETUPARRAY
  USE SPINARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, ITER, BREAKLOOP
  REAL(LATTEPREC) :: ERROR2
  REAL(LATTEPREC) :: R0VEC, P0VEC, R1VEC, XALPHA, XBETA

  IF (SPINON .EQ. 0) THEN

#ifdef DOUBLEPREC
     CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
          BO, HDIM, BO, HDIM, 0.0D0, X2, HDIM)
#elif defined(SINGLEPREC)
     CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
          BO, HDIM, BO, HDIM, 0.0, X2, HDIM)
#endif

     A = TWO*(X2 - BO)

     DO I = 1, HDIM
        A(I,I) = A(I,I) + ONE
     ENDDO

#ifdef DOUBLEPREC
     CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
          A, HDIM, BO, HDIM, 0.0D0, TMPMAT, HDIM)
#elif defined(SINGLEPREC)
     CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
          A, HDIM, BO, HDIM, 0.0, TMPMAT, HDIM)
#endif

     R0 = TMPMAT - X2
     P0 = MINUSONE*R0

     ITER = 0

     BREAKLOOP = 0

     DO WHILE (BREAKLOOP .EQ. 0)

        ITER = ITER + 1

        IF (ITER .EQ. 50) THEN
           CALL PANIC
           CALL ERRORS("solvematcg_sparse",'("SOLVEMATCG NOT CONVERGING")')
        ENDIF


#ifdef DOUBLEPREC
        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
             A, HDIM, P0, HDIM, 0.0D0, TMPMAT, HDIM)
#elif defined(SINGLEPREC)
        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
             A, HDIM, P0, HDIM, 0.0, TMPMAT, HDIM)
#endif

        ERROR2 = ZERO

        !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(GUIDED) &
        !$OMP SHARED (TMPMAT, R0, P0, BO,HDIM) &
        !$OMP PRIVATE(I, J, R0VEC, P0VEC, R1VEC, XALPHA, XBETA) &
        !$OMP REDUCTION(+ : ERROR2)

        DO I = 1, HDIM

           R0VEC = ZERO
           P0VEC = ZERO
           R1VEC = ZERO

           DO J = 1, HDIM

              P0VEC = P0VEC + P0(J,I)*TMPMAT(J,I)
              R0VEC = R0VEC + R0(J,I)*R0(J,I)

           ENDDO

           XALPHA = R0VEC/P0VEC

           DO J = 1, HDIM

              ! New density matrix

              BO(J,I) = BO(J,I) + P0(J,I)*XALPHA

              ! Calculating R1

              R0(J,I) = R0(J,I) + TMPMAT(J,I)*XALPHA

              R1VEC = R1VEC + R0(J,I)*R0(J,I)

           ENDDO

           ERROR2 = ERROR2 + R1VEC

           XBETA = R1VEC/R0VEC

           DO J = 1, HDIM

              P0(J,I) = P0(J,I)*XBETA - R0(J,I)

           ENDDO

        ENDDO

        !$OMP END PARALLEL DO

        IF (ERROR2 .LT. CGTOL2) THEN
           BREAKLOOP = 1
        ENDIF

        !        PRINT*, ITER, ERROR2

     ENDDO

  ELSE

     ! This is the spin-polarized version

#ifdef DOUBLEPREC
     CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
          RHOUP, HDIM, RHOUP, HDIM, 0.0D0, X2, HDIM)
#elif defined(SINGLEPREC)
     CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
          RHOUP, HDIM, RHOUP, HDIM, 0.0, X2, HDIM)
#endif

     A = TWO*(X2 - RHOUP)

     DO I = 1, HDIM
        A(I,I) = A(I,I) + ONE
     ENDDO

#ifdef DOUBLEPREC
     CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
          A, HDIM, RHOUP, HDIM, 0.0D0, TMPMAT, HDIM)
#elif defined(SINGLEPREC)
     CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
          A, HDIM, RHOUP, HDIM, 0.0, TMPMAT, HDIM)
#endif

     R0 = TMPMAT - X2
     P0 = MINUSONE*R0

     ITER = 0

     BREAKLOOP = 0

     DO WHILE (BREAKLOOP .EQ. 0)

        ITER = ITER + 1

        IF (ITER .EQ. 50) THEN
           CALL PANIC
           CALL ERRORS("solvematcg_sparse",'("SOLVEMATCG NOT CONVERGING: SPIN UP")')
        ENDIF


        !     PRINT*, ITER, ERROR2

#ifdef DOUBLEPREC
        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
             A, HDIM, P0, HDIM, 0.0D0, TMPMAT, HDIM)
#elif defined(SINGLEPREC)
        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
             A, HDIM, P0, HDIM, 0.0, TMPMAT, HDIM)
#endif

        ERROR2 = ZERO

        !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(GUIDED) &
        !$OMP SHARED (TMPMAT, R0, P0, RHOUP, HDIM) &
        !$OMP PRIVATE(I, J, R0VEC, P0VEC, R1VEC, XALPHA, XBETA) &
        !$OMP REDUCTION(+ : ERROR2)

        DO I = 1, HDIM

           R0VEC = ZERO
           P0VEC = ZERO
           R1VEC = ZERO

           DO J = 1, HDIM

              P0VEC = P0VEC + P0(J,I)*TMPMAT(J,I)
              R0VEC = R0VEC + R0(J,I)*R0(J,I)

           ENDDO

           XALPHA = R0VEC/P0VEC

           DO J = 1, HDIM

              ! New density matrix

              RHOUP(J,I) = RHOUP(J,I) + P0(J,I)*XALPHA

              ! Calculating R1

              R0(J,I) = R0(J,I) + TMPMAT(J,I)*XALPHA

              R1VEC = R1VEC + R0(J,I)*R0(J,I)

           ENDDO

           ERROR2 = ERROR2 + R1VEC

           XBETA = R1VEC/R0VEC

           DO J = 1, HDIM

              P0(J,I) = P0(J,I)*XBETA - R0(J,I)

           ENDDO

        ENDDO

        !$OMP END PARALLEL DO

        !        PRINT*, "UP ", ITER, ERROR2

        IF (ITER .GT. 3 .AND. ERROR2 .LT. CGTOL2) THEN
           BREAKLOOP = 1
        ENDIF

     ENDDO

#ifdef DOUBLEPREC
     CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
          RHODOWN, HDIM, RHODOWN, HDIM, 0.0D0, X2, HDIM)
#elif defined(SINGLEPREC)
     CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
          RHODOWN, HDIM, RHODOWN, HDIM, 0.0, X2, HDIM)
#endif

     A = TWO*(X2 - RHODOWN)

     DO I = 1, HDIM
        A(I,I) = A(I,I) + ONE
     ENDDO

#ifdef DOUBLEPREC
     CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
          A, HDIM, RHODOWN, HDIM, 0.0D0, TMPMAT, HDIM)
#elif defined(SINGLEPREC)
     CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
          A, HDIM, RHODOWN, HDIM, 0.0, TMPMAT, HDIM)
#endif

     R0 = TMPMAT - X2
     P0 = MINUSONE*R0

     ITER = 0

     BREAKLOOP = 0

     DO WHILE  (BREAKLOOP .EQ. 0)

        ITER = ITER + 1

        IF (ITER .EQ. 50) THEN
           CALL PANIC
           CALL ERRORS("solvematcg_sparse",'("SOLVEMATCG NOT CONVERGING: SPIN DOWN")')
        ENDIF

        !     PRINT*, ITER, ERROR2

#ifdef DOUBLEPREC
        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
             A, HDIM, P0, HDIM, 0.0D0, TMPMAT, HDIM)
#elif defined(SINGLEPREC)
        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
             A, HDIM, P0, HDIM, 0.0, TMPMAT, HDIM)
#endif

        ERROR2 = ZERO

        !$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(GUIDED) &
        !$OMP SHARED (TMPMAT, R0, P0, RHODOWN, HDIM) &
        !$OMP PRIVATE(I, J, R0VEC, P0VEC, R1VEC, XALPHA, XBETA) &
        !$OMP REDUCTION(+ : ERROR2)

        DO I = 1, HDIM

           R0VEC = ZERO
           P0VEC = ZERO
           R1VEC = ZERO

           DO J = 1, HDIM

              P0VEC = P0VEC + P0(J,I)*TMPMAT(J,I)
              R0VEC = R0VEC + R0(J,I)*R0(J,I)

           ENDDO

           XALPHA = R0VEC/P0VEC

           DO J = 1, HDIM

              ! New density matrix

              RHODOWN(J,I) = RHODOWN(J,I) + P0(J,I)*XALPHA

              ! Calculating R1

              R0(J,I) = R0(J,I) + TMPMAT(J,I)*XALPHA

              R1VEC = R1VEC + R0(J,I)*R0(J,I)

           ENDDO

           ERROR2 = ERROR2 + R1VEC

           XBETA = R1VEC/R0VEC

           DO J = 1, HDIM

              P0(J,I) = P0(J,I)*XBETA - R0(J,I)

           ENDDO

        ENDDO

        !$OMP END PARALLEL DO

        !        PRINT*, "DOWN ", ITER, ERROR2

        IF (ITER .GT. 3 .AND. ERROR2 .LT. CGTOL2) THEN
           BREAKLOOP = 1
        ENDIF

     ENDDO

  ENDIF

  RETURN

END SUBROUTINE SOLVEMATCG
