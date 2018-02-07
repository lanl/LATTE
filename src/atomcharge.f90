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

SUBROUTINE ATOMCHARGE(SWITCH)

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: SWITCH
  INTEGER :: I, J, MYINDEX, PREVMYINDEX, ALLOK, ITER, NEWBUILD, II
  REAL(LATTEPREC) :: ATCHG
  REAL(LATTEPREC) :: CORRCHG, MAXCHG
  REAL(LATTEPREC), PARAMETER :: CPROPMAX = 0.8D0, CPROPMIN = 0.01D0
  REAL(LATTEPREC) :: DELTAHII
  REAL(LATTEPREC) :: RIGHTFACT, WRONGFACT
  REAL(LATTEPREC), ALLOCATABLE :: PREVDQ(:), CPROP(:)
  IF (EXISTERROR) RETURN

  ALLOCATE(DELTAQ(NATS), PREVDQ(NATS), CPROP(NATS))

  CPROP = 0.25D0

  RIGHTFACT = 1.1D0
  WRONGFACT = 0.45D0

  MAXCHG = ZERO
  MYINDEX = 0
  PREVMYINDEX = 0
  ALLOK = 0

  !
  ! If we're applying full local charge neutrality, we need to run the
  ! whole caboodle until the DO WHILE loop is satisfied (ie, all atoms
  ! possess a number of electrons to within CHTOL
  !
  ! On the other hand, running one or two iterations of this should be
  ! good enough for an MD simulation when using XBO - we just need an
  ! infinitessimal improvement to ensure stability
  !

  SCFS_II = 0

  ! LCNON = 1: we're going full self-consistency

  IF (LCNON .EQ. 1 .OR. SWITCH .EQ. 0) THEN

     DO I = 1, NATS

        ATCHG = ZERO

        DO J = 1, NOELEM

           IF (ATELE(I) .EQ. ELE(J)) THEN
              CORRCHG = ATOCC(J)
              IF (BASIS(J) .EQ. "sp") THEN
                 MYINDEX = MYINDEX + 4
              ELSEIF (BASIS(J) .EQ. "s") THEN
                 MYINDEX = MYINDEX + 1
              ENDIF
           ENDIF

        ENDDO

        DO J = PREVMYINDEX+1, MYINDEX
           ATCHG = ATCHG + BO(J,J)
        ENDDO

        PREVMYINDEX = MYINDEX

        DELTAQ(I) = ATCHG - CORRCHG

        IF (ABS(DELTAQ(I)) .GT. CHTOL) THEN
           ALLOK = ALLOK + 1
        ENDIF

     ENDDO

     ITER = 0
     NEWBUILD = 0
     !     SCFS_II = 0

     DO WHILE (ALLOK .NE. 0)

        ITER = ITER + 1

        SCFS_II = SCFS_II + 1

        MYINDEX = 0
        PREVMYINDEX = 0

        DO I = 1, NATS

           IF (ITER .GT. 1) THEN

              IF (ABS(DELTAQ(I)) .LT. ABS(PREVDQ(I))) THEN

                 ! Going in the right direction, so go further
                 CPROP(I) = RIGHTFACT*CPROP(I)

              ELSEIF (ABS(DELTAQ(I)) .GT. ABS(PREVDQ(I))) THEN

                 ! Going in the wrong direction, so reverse
                 CPROP(I) = WRONGFACT*CPROP(I)

              ENDIF
           ENDIF

           CPROP(I) = MIN(CPROP(I), CPROPMAX)
           CPROP(I) = MAX(CPROP(I), CPROPMIN)

           DELTAHII = DELTAQ(I)*CPROP(I)

           DO J = 1, NOELEM
              IF (ATELE(I) .EQ. ELE(J)) THEN
                 IF (BASIS(J) .EQ. "sp") THEN
                    MYINDEX = MYINDEX + 4
                 ELSEIF (BASIS(J) .EQ. "s") THEN
                    MYINDEX = MYINDEX + 1
                 ENDIF
              ENDIF
           ENDDO

           DO J = PREVMYINDEX + 1, MYINDEX
              H(J,J) = H(J,J) + DELTAHII
           ENDDO

           PREVMYINDEX = MYINDEX

        ENDDO

        IF (ITER .EQ. 100 ) THEN

           DO I = 1, NATS
              WRITE(6,99) I, DELTAQ(I), ATELE(I), CPROP(I)
           ENDDO
99         FORMAT(I6,1X,F18.9,1X,A2,1X, F12.8)

           CALL ERRORS("atomcharge","LCN not converging: STOP!")

        ENDIF

        IF (CONTROL .EQ. 1) THEN

           CALL DIAGMYH()
           CALL BOEVECS()

        ELSEIF (CONTROL .EQ. 2) THEN

           CALL GERSHGORIN

           IF (SPARSEON .EQ. 0) THEN
              CALL SP2PURE
           ELSEIF (SPARSEON .EQ. 1) THEN

              CALL SP2PURE_SPARSE

           ENDIF

        ELSEIF (CONTROL .EQ. 3) THEN
           !           IF (SPARSEON .EQ. 0) THEN
           CALL FERMIEXPANS
           !           ELSEIF (SPARSEON .EQ. 1) THEN
           !              CALL ALLOCATEPURE
           !              CALL GERSHGORIN
           !              CALL INITSPARSESP2
           !              CALL DEALLOCATEPURE
           !              CALL FERMIEXPANSSPARSE
           !           ENDIF
        ENDIF

        ALLOK = 0

        MYINDEX = 0
        PREVMYINDEX = 0

        PREVDQ = DELTAQ

        DO I = 1, NATS

           ATCHG = ZERO

           DO J = 1, NOELEM

              IF (ATELE(I) .EQ. ELE(J)) THEN
                 CORRCHG = ATOCC(J)
                 IF (BASIS(J) .EQ. "sp") THEN
                    MYINDEX = MYINDEX + 4
                 ELSEIF (BASIS(J) .EQ. "s") THEN
                    MYINDEX = MYINDEX + 1
                 ENDIF
              ENDIF

           ENDDO

           DO J = PREVMYINDEX+1, MYINDEX
              ATCHG = ATCHG + BO(J,J)
           ENDDO

           PREVMYINDEX = MYINDEX

           DELTAQ(I) = ATCHG - CORRCHG

           IF (ABS(DELTAQ(I)) .GT. CHTOL) THEN
              ALLOK = ALLOK + 1
           ENDIF

        ENDDO

     ENDDO

     ! LCNON = 0: we're doing just some specified number of
     ! iterations (= LCNITER)

  ELSEIF (LCNON .EQ. 0 .AND. SWITCH .NE. 0) THEN

     DO II = 1, LCNITER

        MYINDEX = 0
        PREVMYINDEX = 0

        SCFS_II = SCFS_II + 1

        DO I = 1, NATS

           ATCHG = ZERO

           DO J = 1, NOELEM

              IF (ATELE(I) .EQ. ELE(J)) THEN
                 CORRCHG = ATOCC(J)
                 IF (BASIS(J) .EQ. "sp") THEN
                    MYINDEX = MYINDEX + 4
                 ELSEIF (BASIS(J) .EQ. "s") THEN
                    MYINDEX = MYINDEX + 1
                 ENDIF
              ENDIF

           ENDDO

           DO J = PREVMYINDEX+1, MYINDEX
              ATCHG = ATCHG + BO(J,J)
           ENDDO

           !
           ! Note that this factor of 0.25 is an empirical
           ! value that the user is free to adjust. I should probably
           ! move it to an input file
           !

           DELTAHII = QUARTER*(ATCHG - CORRCHG)

           DO J = PREVMYINDEX+1, MYINDEX
              H(J,J) = H(J,J) + DELTAHII
           ENDDO

           PREVMYINDEX = MYINDEX

        ENDDO

        IF ( II .LT. LCNITER ) THEN

           IF (CONTROL .EQ. 1) THEN

              CALL DIAGMYH()
              CALL BOEVECS()

           ELSEIF (CONTROL .EQ. 2) THEN

              CALL GERSHGORIN

              IF (SPARSEON .EQ. 0) THEN
                 CALL SP2PURE
              ELSEIF (SPARSEON .EQ. 1) THEN

                 CALL SP2PURE_SPARSE
              ENDIF

           ELSEIF (CONTROL .EQ. 3) THEN
              !              IF (SPARSEON .EQ. 0) THEN
              CALL FERMIEXPANS
              !              ELSEIF (SPARSEON .EQ. 1) THEN
              !                 CALL ALLOCATEPURE
              !                 CALL GERSHGORIN
              !                 CALL INITSPARSESP2
              !                 CALL DEALLOCATEPURE
              !                 CALL FERMIEXPANSSPARSE
              !              ENDIF
           ENDIF

        ENDIF

     ENDDO

  ENDIF

  DEALLOCATE(DELTAQ, PREVDQ, CPROP)

  RETURN

END SUBROUTINE ATOMCHARGE
