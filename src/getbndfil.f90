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

SUBROUTINE GETBNDFIL()

  USE CONSTANTS_MOD
  USE SETUPARRAY

  IMPLICIT NONE

  INTEGER :: I, J
  INTEGER :: SUMBASIS, NUMORB
  IF (EXISTERROR) RETURN

  TOTNE = ZERO
  SUMBASIS = 0

  DO I = 1, NATS

     TOTNE = TOTNE + ATOCC(ELEMPOINTER(I))

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

     SUMBASIS = SUMBASIS + NUMORB

  ENDDO

  SUMBASIS = 2*SUMBASIS

  !
  ! TOTNE = total number of electrons = tr(rho_up) + tr(rho_down) = tr(bo)
  !

  TOTNE = TOTNE - REAL(CHARGE)

  BNDFIL = TOTNE/REAL(SUMBASIS)
  !  print*, bndfil
  RETURN

END SUBROUTINE GETBNDFIL
