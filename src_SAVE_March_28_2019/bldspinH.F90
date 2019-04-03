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

SUBROUTINE BLDSPINH

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE SPINARRAY
  USE NONOARRAY
  USE MYPRECISION

  IMPLICIT NONE

  !
  ! Here we build the spin-dependent Hamiltonian
  !
  ! Only the onsite blocks are non-zero, and there are only
  ! three unique elements in the case of sp-valent atoms, so we'll
  ! store this H matrix as a vector - suspended, to be inplemented later
  !

  INTEGER :: I, J, K
  INTEGER :: INDEX, HINDEX
  REAL(LATTEPREC) :: DSPINS, DSPINP
  REAL(LATTEPREC) :: HSS, HSP, HPP, HDD, HFF
  IF (EXISTERROR) RETURN

  INDEX = 0
  HINDEX = 0

  ! We start with the spin-indepentent H matrix and 
  ! either add or subtract the spin dependent one for
  ! the up-side and down-spin channels, respectively.


  DO I = 1, NATS

     SELECT CASE(BASIS(ELEMPOINTER(I)))

     CASE("s")

        HSS = WSS(ELEMPOINTER(I))*DELTASPIN(INDEX + 1)
        INDEX = INDEX + 1

        H2VECT(HINDEX + 1) = HSS
        HINDEX = HINDEX + 1

     CASE("p")

        HPP = WPP(ELEMPOINTER(I))*DELTASPIN(INDEX + 1)
        INDEX = INDEX + 1

        H2VECT(HINDEX + 1) = HPP
        H2VECT(HINDEX + 2) = HPP
        H2VECT(HINDEX + 3) = HPP
        HINDEX = HINDEX + 3

     CASE("d")

        HDD = WDD(ELEMPOINTER(I))*DELTASPIN(INDEX + 1)
        INDEX = INDEX + 1

        H2VECT(HINDEX + 1) = HDD
        H2VECT(HINDEX + 2) = HDD
        H2VECT(HINDEX + 3) = HDD
        H2VECT(HINDEX + 4) = HDD
        H2VECT(HINDEX + 5) = HDD
        HINDEX = HINDEX + 5

     CASE("f")

        HFF = WFF(ELEMPOINTER(I))*DELTASPIN(INDEX + 1)
        INDEX = INDEX + 1

        H2VECT(HINDEX + 1) = HFF
        H2VECT(HINDEX + 2) = HFF
        H2VECT(HINDEX + 3) = HFF
        H2VECT(HINDEX + 4) = HFF
        H2VECT(HINDEX + 5) = HFF
        H2VECT(HINDEX + 6) = HFF
        H2VECT(HINDEX + 7) = HFF
        HINDEX = HINDEX + 7

     CASE("sp")

        HSS = WSS(ELEMPOINTER(I))*DELTASPIN(INDEX + 1)
        HPP = WPP(ELEMPOINTER(I))*DELTASPIN(INDEX + 2)
        INDEX = INDEX + 2

        H2VECT(HINDEX + 1) = HSS
        H2VECT(HINDEX + 2) = HPP
        H2VECT(HINDEX + 3) = HPP
        H2VECT(HINDEX + 4) = HPP
        HINDEX = HINDEX + 4

     CASE("sd")

        HSS = WSS(ELEMPOINTER(I))*DELTASPIN(INDEX + 1)
        HDD = WDD(ELEMPOINTER(I))*DELTASPIN(INDEX + 2)
        INDEX = INDEX + 2

        H2VECT(HINDEX + 1) = HSS
        H2VECT(HINDEX + 2) = HDD
        H2VECT(HINDEX + 3) = HDD
        H2VECT(HINDEX + 4) = HDD
        H2VECT(HINDEX + 5) = HDD
        H2VECT(HINDEX + 6) = HDD
        HINDEX = HINDEX + 6

     CASE("sf")

        HSS = WSS(ELEMPOINTER(I))*DELTASPIN(INDEX + 1)
        HFF = WFF(ELEMPOINTER(I))*DELTASPIN(INDEX + 2)
        INDEX = INDEX + 2

        H2VECT(HINDEX + 1) = HSS
        H2VECT(HINDEX + 2) = HFF
        H2VECT(HINDEX + 3) = HFF
        H2VECT(HINDEX + 4) = HFF
        H2VECT(HINDEX + 5) = HFF
        H2VECT(HINDEX + 6) = HFF
        H2VECT(HINDEX + 7) = HFF
        H2VECT(HINDEX + 8) = HFF
        HINDEX = HINDEX + 8

     CASE("pd")

        HPP = WPP(ELEMPOINTER(I))*DELTASPIN(INDEX + 1)
        HDD = WDD(ELEMPOINTER(I))*DELTASPIN(INDEX + 2)
        INDEX = INDEX + 2

        H2VECT(HINDEX + 1) = HPP
        H2VECT(HINDEX + 2) = HPP
        H2VECT(HINDEX + 3) = HPP
        H2VECT(HINDEX + 4) = HDD
        H2VECT(HINDEX + 5) = HDD
        H2VECT(HINDEX + 6) = HDD
        H2VECT(HINDEX + 7) = HDD
        H2VECT(HINDEX + 8) = HDD
        HINDEX = HINDEX + 8

     CASE("pf")

        HPP = WPP(ELEMPOINTER(I))*DELTASPIN(INDEX + 1)
        HFF = WFF(ELEMPOINTER(I))*DELTASPIN(INDEX + 2)
        INDEX = INDEX + 2

        H2VECT(HINDEX + 1) = HPP
        H2VECT(HINDEX + 2) = HPP
        H2VECT(HINDEX + 3) = HPP
        H2VECT(HINDEX + 4) = HFF
        H2VECT(HINDEX + 5) = HFF
        H2VECT(HINDEX + 6) = HFF
        H2VECT(HINDEX + 7) = HFF
        H2VECT(HINDEX + 8) = HFF
        H2VECT(HINDEX + 9) = HFF
        H2VECT(HINDEX + 10) = HFF
        HINDEX = HINDEX + 10


     CASE("df")

        HDD = WDD(ELEMPOINTER(I))*DELTASPIN(INDEX + 1)
        HFF = WFF(ELEMPOINTER(I))*DELTASPIN(INDEX + 2)
        INDEX = INDEX + 2

        H2VECT(HINDEX + 1) = HDD
        H2VECT(HINDEX + 2) = HDD
        H2VECT(HINDEX + 3) = HDD
        H2VECT(HINDEX + 4) = HDD
        H2VECT(HINDEX + 5) = HDD
        H2VECT(HINDEX + 6) = HFF
        H2VECT(HINDEX + 7) = HFF
        H2VECT(HINDEX + 8) = HFF
        H2VECT(HINDEX + 9) = HFF
        H2VECT(HINDEX + 10) = HFF
        H2VECT(HINDEX + 11) = HFF
        H2VECT(HINDEX + 12) = HFF
        HINDEX = HINDEX + 12

     CASE("spd")

        HSS = WSS(ELEMPOINTER(I))*DELTASPIN(INDEX + 1)
        HPP = WPP(ELEMPOINTER(I))*DELTASPIN(INDEX + 2)
        HDD = WDD(ELEMPOINTER(I))*DELTASPIN(INDEX + 3)
        INDEX = INDEX + 3

        H2VECT(HINDEX + 1) = HSS
        H2VECT(HINDEX + 2) = HPP
        H2VECT(HINDEX + 3) = HPP
        H2VECT(HINDEX + 4) = HPP
        H2VECT(HINDEX + 5) = HDD
        H2VECT(HINDEX + 6) = HDD
        H2VECT(HINDEX + 7) = HDD
        H2VECT(HINDEX + 8) = HDD
        H2VECT(HINDEX + 9) = HDD
        HINDEX = HINDEX + 9

     CASE("spf")

        HSS = WSS(ELEMPOINTER(I))*DELTASPIN(INDEX + 1)
        HPP = WPP(ELEMPOINTER(I))*DELTASPIN(INDEX + 2)
        HFF = WFF(ELEMPOINTER(I))*DELTASPIN(INDEX + 3)
        INDEX = INDEX + 3

        H2VECT(HINDEX + 1) = HSS
        H2VECT(HINDEX + 2) = HPP
        H2VECT(HINDEX + 3) = HPP
        H2VECT(HINDEX + 4) = HPP
        H2VECT(HINDEX + 5) = HFF
        H2VECT(HINDEX + 6) = HFF
        H2VECT(HINDEX + 7) = HFF
        H2VECT(HINDEX + 8) = HFF
        H2VECT(HINDEX + 9) = HFF
        H2VECT(HINDEX + 10) = HFF
        H2VECT(HINDEX + 11) = HFF
        HINDEX = HINDEX + 11

     CASE("sdf")

        HSS = WSS(ELEMPOINTER(I))*DELTASPIN(INDEX + 1)
        HDD = WDD(ELEMPOINTER(I))*DELTASPIN(INDEX + 2)
        HFF = WFF(ELEMPOINTER(I))*DELTASPIN(INDEX + 3)
        INDEX = INDEX + 3

        H2VECT(HINDEX + 1) = HSS
        H2VECT(HINDEX + 2) = HDD
        H2VECT(HINDEX + 3) = HDD
        H2VECT(HINDEX + 4) = HDD
        H2VECT(HINDEX + 5) = HDD
        H2VECT(HINDEX + 6) = HDD
        H2VECT(HINDEX + 7) = HFF
        H2VECT(HINDEX + 8) = HFF
        H2VECT(HINDEX + 9) = HFF
        H2VECT(HINDEX + 10) = HFF
        H2VECT(HINDEX + 11) = HFF
        H2VECT(HINDEX + 12) = HFF
        H2VECT(HINDEX + 13) = HFF
        HINDEX = HINDEX + 13

     CASE("pdf")

        HPP = WPP(ELEMPOINTER(I))*DELTASPIN(INDEX + 1)
        HDD = WDD(ELEMPOINTER(I))*DELTASPIN(INDEX + 2)
        HFF = WFF(ELEMPOINTER(I))*DELTASPIN(INDEX + 3)
        INDEX = INDEX + 3

        H2VECT(HINDEX + 1) = HPP
        H2VECT(HINDEX + 2) = HPP
        H2VECT(HINDEX + 3) = HPP
        H2VECT(HINDEX + 4) = HDD
        H2VECT(HINDEX + 5) = HDD
        H2VECT(HINDEX + 6) = HDD
        H2VECT(HINDEX + 7) = HDD
        H2VECT(HINDEX + 8) = HDD
        H2VECT(HINDEX + 9) = HFF
        H2VECT(HINDEX + 10) = HFF
        H2VECT(HINDEX + 11) = HFF
        H2VECT(HINDEX + 12) = HFF
        H2VECT(HINDEX + 13) = HFF
        H2VECT(HINDEX + 14) = HFF
        H2VECT(HINDEX + 15) = HFF
        HINDEX = HINDEX + 15

     CASE("spdf")

        HSS = WSS(ELEMPOINTER(I))*DELTASPIN(INDEX + 1)
        HPP = WPP(ELEMPOINTER(I))*DELTASPIN(INDEX + 2)
        HDD = WDD(ELEMPOINTER(I))*DELTASPIN(INDEX + 3)
        HFF = WFF(ELEMPOINTER(I))*DELTASPIN(INDEX + 4)
        INDEX = INDEX + 4

        H2VECT(HINDEX + 1) = HSS
        H2VECT(HINDEX + 2) = HPP
        H2VECT(HINDEX + 3) = HPP
        H2VECT(HINDEX + 4) = HPP
        H2VECT(HINDEX + 5) = HDD
        H2VECT(HINDEX + 6) = HDD
        H2VECT(HINDEX + 7) = HDD
        H2VECT(HINDEX + 8) = HDD
        H2VECT(HINDEX + 9) = HDD
        H2VECT(HINDEX + 10) = HFF
        H2VECT(HINDEX + 11) = HFF
        H2VECT(HINDEX + 12) = HFF
        H2VECT(HINDEX + 13) = HFF
        H2VECT(HINDEX + 14) = HFF
        H2VECT(HINDEX + 15) = HFF
        H2VECT(HINDEX + 16) = HFF
        HINDEX = HINDEX + 16

     END SELECT

  ENDDO

  IF (BASISTYPE .EQ. "ORTHO") THEN

     HUP = H
     HDOWN = H

     DO I = 1, HDIM
        HUP(I,I) = HUP(I,I) + H2VECT(I)
        HDOWN(I,I) = HDOWN(I,I) - H2VECT(I)
     ENDDO

  ELSEIF (BASISTYPE .EQ. "NONORTHO") THEN

     ! Now we have H_2, we can form SH_2, and then symmetrize it

     DO I = 1, HDIM
        DO J = 1, HDIM

           SH2(J,I) = SMAT(J,I)*H2VECT(I)

        ENDDO
     ENDDO

     SH2 = (SH2 + TRANSPOSE(SH2))/TWO

     ! H_up = H + SH2 ; H_down = H = SH2

     HUP = H + SH2

     HDOWN = H - SH2

  ENDIF

  RETURN

END SUBROUTINE BLDSPINH




