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

SUBROUTINE ALLOCATEXBO

  USE CONSTANTS_MOD
  USE XBOARRAY
  USE SPINARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: ARRAYDIM
  IF (EXISTERROR) RETURN

  CNK = ZERO

  !
  ! This depends on whether we're propagating the
  ! on-site H matrix elements (in the case of LCN calculations)
  ! or the atomic partial charges (full SCC-TB)
  !

  IF (ELECTRO .EQ. 0) THEN
     ARRAYDIM = NATS
     !     ARRAYDIM = HDIM
  ELSE
     ARRAYDIM = NATS
  ENDIF

  IF (XBODISON .EQ. 0) THEN

     ALLOCATE(PNK(2,ARRAYDIM))

     !
     ! If we need to propagate the chemical potential
     !

     IF (CONTROL .EQ. 1 .OR. CONTROL .EQ. 3 &
          .OR. CONTROL .EQ. 5) ALLOCATE( CHEMPOT_PNK(2) )


     !
     ! If we need to propagate the spin difference 
     ! densities too:
     !

     IF (SPINON .EQ. 1) ALLOCATE( SPIN_PNK(2, DELTADIM) )

  ELSE

     ALLOCATE( PNK(XBODISORDER + 1, ARRAYDIM) )

     !
     ! If we need to propagate the chemical potential
     !

     IF (CONTROL .EQ. 1 .OR. CONTROL .EQ. 3 .OR. &
          CONTROL .EQ. 5) ALLOCATE( CHEMPOT_PNK(XBODISORDER + 1) )

     !
     ! If we need to propagate the spin difference 
     ! densities too:
     !

     IF (SPINON .EQ. 1) ALLOCATE( SPIN_PNK(XBODISORDER + 1, DELTADIM) )

     IF (XBODISORDER .EQ. 3) THEN

        ALPHA_XBO = 150.0D-3
        KAPPA_XBO = 1.69D0
        CNK(1) = -2.0D0
        CNK(2) = 3.0D0
        CNK(3) = 0.0D0
        CNK(4) = -1.0D0

     ELSEIF (XBODISORDER .EQ. 4) THEN

        ALPHA_XBO = 57.0D-3
        KAPPA_XBO = 1.75D0
        CNK(1) = -3.0D0
        CNK(2) = 6.0D0
        CNK(3) = -2.0D0
        CNK(4) = -2.0D0
        CNK(5) = 1.0D0

     ELSEIF (XBODISORDER .EQ. 5) THEN

        ALPHA_XBO = 18.0D-3
        KAPPA_XBO = 1.82D0
        CNK(1) = -6.0D0
        CNK(2) = 14.0D0
        CNK(3) = -8.0D0
        CNK(4) = -3.0D0
        CNK(5) = 4.0D0
        CNK(6) = -1.0D0

     ELSEIF (XBODISORDER .EQ. 6) THEN

        ALPHA_XBO = 5.5D-3
        KAPPA_XBO = 1.84D0
        CNK(1) = -14.0D0
        CNK(2) = 36.0D0
        CNK(3) = -27.0D0
        CNK(4) = -2.0D0
        CNK(5) = 12.0D0
        CNK(6) = -6.0D0
        CNK(7) = 1.0D0

     ELSEIF (XBODISORDER .EQ. 7) THEN

        ALPHA_XBO = 1.6D-3
        KAPPA_XBO = 1.86D0
        CNK(1) = -36.0D0
        CNK(2) = 99.0D0
        CNK(3) = -88.0D0
        CNK(4) = 11.0D0
        CNK(5) = 32.0D0
        CNK(6) = -25.0D0
        CNK(7) = 8.0D0
        CNK(8) = -1.0D0

     ELSEIF (XBODISORDER .EQ. 8) THEN

        ALPHA_XBO = 0.44D-3
        KAPPA_XBO = 1.88D0
        CNK(1) = -99.0D0
        CNK(2) = 286.0D0
        CNK(3) = -286.0D0
        CNK(4) = 78.0D0
        CNK(5) = 78.0D0
        CNK(6) = -90.0D0
        CNK(7) = 42.0D0
        CNK(8) = -10.0D0
        CNK(9) = 1.0D0

     ELSEIF (XBODISORDER .EQ. 9) THEN

        ALPHA_XBO = 0.12D-3
        KAPPA_XBO = 1.89D0
        CNK(1) = -286.0D0
        CNK(2) = 858.0D0
        CNK(3) = -936.0D0
        CNK(4) = 364.0D0
        CNK(5) = 168.0D0
        CNK(6) = -300.0D0
        CNK(7) = 184.0D0
        CNK(8) = -63.0D0
        CNK(9) = 12.0D0
        CNK(10) = -1.0D0

     ENDIF

  ENDIF

  RETURN

END SUBROUTINE ALLOCATEXBO
