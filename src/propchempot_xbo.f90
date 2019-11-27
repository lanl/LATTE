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

SUBROUTINE PROPCHEMPOT(ITER)

  USE CONSTANTS_MOD
  USE XBOARRAY
  USE SETUPARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, ITER  
  IF (EXISTERROR) RETURN

  ! If iter = 1, then we just set up our arrays of previous guesses

  IF (ITER .EQ. 1) THEN

     ! Dissipation off

     IF (XBODISON .EQ. 0) THEN

        CHEMPOT_PNK(1) = CHEMPOT
        CHEMPOT_PNK(2) = CHEMPOT_PNK(1)

        ! Dissipation on

     ELSE

        CHEMPOT_PNK(1) = CHEMPOT

        DO I = 2, XBODISORDER + 1
           CHEMPOT_PNK(I) = CHEMPOT_PNK(I-1)
        ENDDO

     ENDIF

  ELSE   ! For ITER > 1

     IF (XBODISON .EQ. 0) THEN

        CHEMPOT = TWO*CHEMPOT_PNK(1) - CHEMPOT_PNK(2) + &
             TWO*(CHEMPOT - CHEMPOT_PNK(1))

        CHEMPOT_PNK(2) = CHEMPOT_PNK(1)
        CHEMPOT_PNK(1) = CHEMPOT

     ELSEIF (XBODISON .EQ. 1) THEN

        CHEMPOT = TWO*CHEMPOT_PNK(1) - CHEMPOT_PNK(2) + &
             KAPPA_SCALE*KAPPA_XBO*(CHEMPOT - CHEMPOT_PNK(1))

        DO I = 1, XBODISORDER + 1

           CHEMPOT = CHEMPOT + ALPHA_XBO*CNK(I)*CHEMPOT_PNK(I)

        ENDDO

        DO I = 1, XBODISORDER 

           CHEMPOT_PNK(XBODISORDER + 2 - I) = &
                CHEMPOT_PNK(XBODISORDER + 1 - I)

        ENDDO

        CHEMPOT_PNK(1) = CHEMPOT

     ENDIF

  ENDIF

  RETURN

END SUBROUTINE PROPCHEMPOT


