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

SUBROUTINE GETHDIM

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE SPINARRAY
  USE RESTARTARRAY
  USE KSPACEARRAY

  IMPLICIT NONE

  INTEGER :: I, J, NUMORB
  IF (EXISTERROR) RETURN

  DELTADIM = 0

  IF (RESTART .EQ. 0) THEN

     HDIM = 0

     DO I = 1, NATS

        SELECT CASE(BASIS(ELEMPOINTER(I)))

        CASE("s")
           NUMORB = 1
        CASE("p")
           NUMORB = 3
        CASE("d")
           NUMORB = 5
        CASE("f")
           NUMORB = 7
        CASE("sp")
           NUMORB = 4
        CASE("sd")
           NUMORB = 6
        CASE("sf")
           NUMORB = 8
        CASE("pd")
           NUMORB = 8
        CASE("pf")
           NUMORB = 10
        CASE("df")
           NUMORB = 12
        CASE("spd")
           NUMORB = 9
        CASE("spf")
           NUMORB = 11
        CASE("sdf")
           NUMORB = 13
        CASE("pdf")
           NUMORB = 15
        CASE("spdf") 
           NUMORB = 16
        END SELECT

        HDIM = HDIM + NUMORB

     ENDDO

  ELSEIF (RESTART .EQ. 1) THEN

     HDIM = TMPHDIM

  ENDIF

  ! For use when getting the partial charges

  ALLOCATE(QLIST(HDIM))

  ! Allocate the Hamiltonian matrix

  IF (KON .EQ. 0) THEN

     ! Real space

     ALLOCATE(H(HDIM,HDIM), HDIAG(HDIM))

  ELSE ! k-space

     ALLOCATE(HK(HDIM, HDIM, NKTOT), HKDIAG(HDIM, NKTOT))

  ENDIF

  ! If we're using a non-orthogonal basis:

  IF (BASISTYPE .EQ. "NONORTHO") THEN
     IF (KON .EQ. 0) THEN
        ALLOCATE(H0(HDIM, HDIM))
     ELSE
        ALLOCATE(HK0(HDIM, HDIM, NKTOT))
     ENDIF
  ENDIF

  IF (SPINON .EQ. 0) THEN

     ! No spins: allocate 1 double-occupied bond order matrix

     IF (KON .EQ. 0) THEN

        ALLOCATE(BO(HDIM,HDIM))

     ELSE

        ALLOCATE(KBO(HDIM, HDIM, NKTOT))

     ENDIF

  ELSEIF (SPINON .EQ. 1) THEN

     ! We're going to have two Hamiltonians because I can't
     ! figure out a more elegant way to do it just yet...

     ! With spins, we need spin-up and spin-down density matrices

     ALLOCATE(HUP(HDIM, HDIM), HDOWN(HDIM, HDIM))
     ALLOCATE(RHOUP(HDIM, HDIM), RHODOWN(HDIM, HDIM))
     ALLOCATE(H2VECT(HDIM))

     HUP = ZERO
     HDOWN = ZERO
     RHOUP = ZERO
     RHODOWN = ZERO

     ! And our spin-dependent H_(2) matrix

     DO I = 1, NATS

        SELECT CASE(BASIS(ELEMPOINTER(I)))

        CASE("s")
           NUMORB = 1
        CASE("p")
           NUMORB = 1
        CASE("d")
           NUMORB = 1
        CASE("f")

           NUMORB = 1

        CASE("sp")

           NUMORB = 2

        CASE("sd")

           NUMORB = 2

        CASE("sf")

           NUMORB = 2

        CASE("pd")
           NUMORB = 2
        CASE("pf")
           NUMORB = 2
        CASE("df")
           NUMORB = 2
        CASE("spd")
           NUMORB = 3
        CASE("spf")
           NUMORB = 3
        CASE("sdf")
           NUMORB = 3
        CASE("pdf")
           NUMORB = 3
        CASE("spdf") 
           NUMORB = 4
        END SELECT

        DELTADIM = DELTADIM + NUMORB

     ENDDO

     ! This array is required when we calculate Mulliken spin densities
     ! with a non-orthogonal basis

     IF (BASISTYPE .EQ. "NONORTHO") ALLOCATE(SPINLIST(HDIM))

     ALLOCATE(DELTASPIN(DELTADIM), OLDDELTASPIN(DELTADIM))

     DELTASPIN = ZERO
     OLDDELTASPIN = ZERO

  ENDIF

  RETURN

END SUBROUTINE GETHDIM
