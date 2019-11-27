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

SUBROUTINE XBO(ITER)

  !
  ! This here subroutine implements Niklasson's eXtended
  ! Born-Oppenheimer, time-reversible, and most excellent MD thing.
  !
  ! Niklasson, PRL, vol. 100, p 123004, (2008)
  !

  USE CONSTANTS_MOD
  USE XBOARRAY
  USE SETUPARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, ITER
  IF (EXISTERROR) RETURN

  !
  ! If we have no damping:
  !
  ! Our time reversible guess for the next set of diagonal H-matrix
  ! elements depend on those from the last guess, the guess before that, and
  ! the H-matrix elements we calculated subject to contraints on Delta_q
  ! from the last iteration.
  !
  ! P(t) - last guess
  ! P(t - dt) - guess before that
  ! D(t) - the elements calculated from the last run
  !
  ! P(t + dt) = 2P(t) - P(t - dt) + 2[D(t) - P(t)]
  !
  ! With damping, it's the same general idea but we use guesses from even
  ! earlier iterations too. Preliminary testing shows better energy
  ! conservation with high 'K'.
  !

  ! Required for the Fast-QMMD scheme

  IF (QITER .EQ. 0) THEN
     KAPPA_SCALE = MDMIX
  ELSE
     KAPPA_SCALE = ONE
  ENDIF


  IF (ITER .EQ. 1) THEN

     IF (ELECTRO .EQ. 1) THEN

        IF (XBODISON .EQ. 0) THEN

           DO I = 1, NATS
              PNK(1,I) = DELTAQ(I)
              PNK(2,I) = PNK(1,I)
           ENDDO

        ELSEIF (XBODISON .EQ. 1) THEN

           DO I = 1, NATS
              PNK(1,I) = DELTAQ(I)

              DO J = 2, XBODISORDER + 1
                 PNK(J,I) = PNK(J-1,I)
              ENDDO

           ENDDO

        ENDIF

     ELSEIF (ELECTRO .EQ. 0) THEN

        IF (XBODISON .EQ. 0) THEN

           DO I = 1, NATS
              PNK(1,I) = LCNSHIFT(I)
              PNK(2,I) = PNK(1,I)
           ENDDO

        ELSEIF (XBODISON .EQ. 1) THEN

           DO I = 1, NATS
              PNK(1,I) = LCNSHIFT(I)

              DO J = 2, XBODISORDER + 1
                 PNK(J,I) = PNK(J-1,I)
              ENDDO

           ENDDO

        ENDIF

     ENDIF

  ELSEIF (ITER .GT. 1) THEN

     IF (XBODISON .EQ. 0) THEN

        IF (ELECTRO .EQ. 0) THEN

           DO I = 1, NATS

              LCNSHIFT(I) = TWO*PNK(1,I) - PNK(2,I) + &
                   TWO*(LCNSHIFT(I) - PNK(1,I))

              PNK(2,I) = PNK(1,I)
              PNK(1,I) = LCNSHIFT(I)

           ENDDO

        ELSEIF (ELECTRO .EQ. 1) THEN

           DO I = 1, NATS

              DELTAQ(I) = TWO*PNK(1,I) - PNK(2,I) + &
                   TWO*(DELTAQ(I) - PNK(1,I))

              PNK(2,I) = PNK(1,I)
              PNK(1,I) = DELTAQ(I)

           ENDDO

        ENDIF

     ELSEIF (XBODISON .EQ. 1) THEN

        IF (ELECTRO .EQ. 0) THEN

           DO I = 1, NATS

              LCNSHIFT(I) =  TWO*PNK(1,I) - PNK(2,I) + &
                   KAPPA_SCALE*KAPPA_XBO*(LCNSHIFT(I) - PNK(1,I))

              DO J = 1, XBODISORDER+1

                 LCNSHIFT(I) = LCNSHIFT(I) + ALPHA_XBO*CNK(J)*PNK(J,I)

              ENDDO

              DO J = 1, XBODISORDER

                 PNK(XBODISORDER+2 - J,I) = PNK(XBODISORDER+1 - J,I)

              ENDDO

              PNK(1,I) = LCNSHIFT(I)

           ENDDO

        ELSEIF (ELECTRO .EQ. 1) THEN

           DO I = 1, NATS

              DELTAQ(I) =  TWO*PNK(1,I) - PNK(2,I) + &
                   KAPPA_SCALE*KAPPA_XBO*(DELTAQ(I) - PNK(1,I))

              DO J = 1, XBODISORDER+1

                 DELTAQ(I) = DELTAQ(I) + ALPHA_XBO*CNK(J)*PNK(J,I)

              ENDDO

              DO J = 1, XBODISORDER

                 PNK(XBODISORDER+2 - J,I) = PNK(XBODISORDER+1 - J,I)

              ENDDO

              PNK(1,I) = DELTAQ(I)

           ENDDO

        ENDIF

     ENDIF

  ENDIF

  RETURN

END SUBROUTINE XBO
