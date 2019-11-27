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

SUBROUTINE GETMATINDLIST

  USE CONSTANTS_MOD
  USE SETUPARRAY

  IMPLICIT NONE

  INTEGER :: I, J, INDI
  IF (EXISTERROR) RETURN

  ! This arrays contains and index for the position of atoms
  ! in the Hamiltonian matrix. We only need to compute this once

  IF(.NOT. ALLOCATED(MATINDLIST)) ALLOCATE(MATINDLIST(NATS))

  DO I = 1, NATS
     INDI = 0
     DO J = 1, I - 1

        SELECT CASE(BASIS(ELEMPOINTER(J)))
        CASE("s")
           INDI = INDI + 1
        CASE("p")
           INDI = INDI + 3
        CASE("d")
           INDI = INDI + 5
        CASE("f")
           INDI = INDI + 7
        CASE("sp")
           INDI = INDI + 4
        CASE("sd")
           INDI = INDI + 6
        CASE("sf")
           INDI = INDI + 8
        CASE("pd")
           INDI = INDI + 8
        CASE("pf")
           INDI = INDI + 10
        CASE("df")
           INDI = INDI + 12
        CASE("spd")
           INDI = INDI + 9
        CASE("spf")
           INDI = INDI + 11
        CASE("sdf")
           INDI = INDI + 13
        CASE("pdf")
           INDI = INDI + 15
        CASE("spdf")
           INDI = INDI + 16
        END SELECT

     ENDDO

     MATINDLIST(I) = INDI

  ENDDO

  IF (SPINON .EQ. 1) THEN

     ALLOCATE(SPININDLIST(NATS))

     DO I = 1, NATS
        INDI = 0
        DO J = 1, I - 1

           ! The number of spins = number of orbitals

           SELECT CASE(BASIS(ELEMPOINTER(J)))
           CASE("s")
              INDI = INDI + 1
           CASE("p")
              INDI = INDI + 1
           CASE("d")
              INDI = INDI + 1
           CASE("f")
              INDI = INDI + 1
           CASE("sp")
              INDI = INDI + 2
           CASE("sd")
              INDI = INDI + 2
           CASE("sf")
              INDI = INDI + 2
           CASE("pd")
              INDI = INDI + 2
           CASE("pf")
              INDI = INDI + 2
           CASE("df")
              INDI = INDI + 2
           CASE("spd")
              INDI = INDI + 3
           CASE("spf")
              INDI = INDI + 3
           CASE("sdf")
              INDI = INDI + 3
           CASE("pdf")
              INDI = INDI + 3
           CASE("spdf")
              INDI = INDI + 4
           END SELECT

        ENDDO

        SPININDLIST(I) = INDI

     ENDDO

  ENDIF

  RETURN

END SUBROUTINE GETMATINDLIST
