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

SUBROUTINE PROPSPINS(ITER)

  USE CONSTANTS_MOD
  USE XBOARRAY
  USE SETUPARRAY
  USE SPINARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, ITER  
  IF (EXISTERROR) RETURN

  ! If iter = 1, then we just set up our arrays of previous guesses

  IF (ITER .EQ. 1) THEN

     ! Dissipation off

     IF (XBODISON .EQ. 0) THEN

        DO I = 1, DELTADIM
           SPIN_PNK(1,I) = DELTASPIN(I)
           SPIN_PNK(2,I) = SPIN_PNK(1,I)
        ENDDO

        ! Dissipation on

     ELSEIF (XBODISON .EQ. 1) THEN

        DO I = 1, DELTADIM

           SPIN_PNK(1,I) = DELTASPIN(I)

           DO J = 2, XBODISORDER + 1
              SPIN_PNK(J,I) = SPIN_PNK(J-1,I)
           ENDDO

        ENDDO

     ENDIF

  ELSE        ! ITER > 1

     IF (XBODISON .EQ. 0) THEN

        DO I = 1, DELTADIM

           DELTASPIN(I) = TWO*SPIN_PNK(1,I) - SPIN_PNK(2,I) + &
                TWO*(DELTASPIN(I) - SPIN_PNK(1,I))

           SPIN_PNK(2,I) = SPIN_PNK(1,I)
           SPIN_PNK(1,I) = DELTASPIN(I)

        ENDDO

     ELSE ! Dissipation on (XBODISON .EQ. 1)

        DO I = 1, DELTADIM

           DELTASPIN(I) = TWO*SPIN_PNK(1,I) - SPIN_PNK(2,I) + &
                KAPPA_SCALE*KAPPA_XBO*(DELTASPIN(I) - SPIN_PNK(1,I))

           DO J = 1, XBODISORDER + 1

              DELTASPIN(I) = DELTASPIN(I) + &
                   ALPHA_XBO*CNK(J)*SPIN_PNK(J,I)

           ENDDO

           DO J = 1, XBODISORDER 

              SPIN_PNK(XBODISORDER + 2 - J, I) = &
                   SPIN_PNK(XBODISORDER + 1 - J, I)

           ENDDO

           SPIN_PNK(1,I) = DELTASPIN(I)

        ENDDO

     ENDIF

  ENDIF

  RETURN

END SUBROUTINE PROPSPINS


