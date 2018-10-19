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

SUBROUTINE GETDELTASPIN

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE SPINARRAY
  USE NONOARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, INDEX, DINDEX
  IF (EXISTERROR) RETURN

  INDEX = 0
  DINDEX = 0

  !  deltaspin = number of spin ups - number of spin downs

  IF (BASISTYPE .EQ. "ORTHO") THEN

     DO I = 1, NATS

        SELECT CASE(BASIS(ELEMPOINTER(I)))

        CASE("s")

           INDEX = INDEX + 1
           DINDEX = DINDEX + 1

           ! s

           DELTASPIN(DINDEX) = RHOUP(INDEX, INDEX) - RHODOWN(INDEX,INDEX)

        CASE("p")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) + &
                RHOUP(INDEX + 2, INDEX + 2) - RHODOWN(INDEX + 2, INDEX + 2) + &
                RHOUP(INDEX + 3, INDEX + 3) - RHODOWN(INDEX + 3, INDEX + 3)


           INDEX = INDEX + 3

        CASE("d")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) + &
                RHOUP(INDEX + 2, INDEX + 2) - RHODOWN(INDEX + 2, INDEX + 2) + &
                RHOUP(INDEX + 3, INDEX + 3) - RHODOWN(INDEX + 3, INDEX + 3) + &
                RHOUP(INDEX + 4, INDEX + 4) - RHODOWN(INDEX + 4, INDEX + 4) + &
                RHOUP(INDEX + 5, INDEX + 5) - RHODOWN(INDEX + 5, INDEX + 5)

           INDEX = INDEX + 5

        CASE("f")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) + &
                RHOUP(INDEX + 2, INDEX + 2) - RHODOWN(INDEX + 2, INDEX + 2) + &
                RHOUP(INDEX + 3, INDEX + 3) - RHODOWN(INDEX + 3, INDEX + 3) + &
                RHOUP(INDEX + 4, INDEX + 4) - RHODOWN(INDEX + 4, INDEX + 4) + &
                RHOUP(INDEX + 5, INDEX + 5) - RHODOWN(INDEX + 5, INDEX + 5) + &
                RHOUP(INDEX + 6, INDEX + 6) - RHODOWN(INDEX + 6, INDEX + 6) + &
                RHOUP(INDEX + 7, INDEX + 7) - RHODOWN(INDEX + 7, INDEX + 7)

           INDEX = INDEX + 7

        CASE("sp")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) 

           INDEX = INDEX + 1

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) + &
                RHOUP(INDEX + 2, INDEX + 2) - RHODOWN(INDEX + 2, INDEX + 2) + &
                RHOUP(INDEX + 3, INDEX + 3) - RHODOWN(INDEX + 3, INDEX + 3)

           INDEX = INDEX + 3 

        CASE("sd")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) 

           INDEX = INDEX + 1

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) + &
                RHOUP(INDEX + 2, INDEX + 2) - RHODOWN(INDEX + 2, INDEX + 2) + &
                RHOUP(INDEX + 3, INDEX + 3) - RHODOWN(INDEX + 3, INDEX + 3) + &
                RHOUP(INDEX + 4, INDEX + 4) - RHODOWN(INDEX + 4, INDEX + 4) + &
                RHOUP(INDEX + 5, INDEX + 5) - RHODOWN(INDEX + 5, INDEX + 5)

           INDEX = INDEX + 5 

        CASE("sf")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) 

           INDEX = INDEX + 1

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) + &
                RHOUP(INDEX + 2, INDEX + 2) - RHODOWN(INDEX + 2, INDEX + 2) + &
                RHOUP(INDEX + 3, INDEX + 3) - RHODOWN(INDEX + 3, INDEX + 3) + &
                RHOUP(INDEX + 4, INDEX + 4) - RHODOWN(INDEX + 4, INDEX + 4) + &
                RHOUP(INDEX + 5, INDEX + 5) - RHODOWN(INDEX + 5, INDEX + 5) + &
                RHOUP(INDEX + 6, INDEX + 6) - RHODOWN(INDEX + 6, INDEX + 6) + &
                RHOUP(INDEX + 7, INDEX + 7) - RHODOWN(INDEX + 7, INDEX + 7)


           INDEX = INDEX + 7

        CASE("pd")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) + &
                RHOUP(INDEX + 2, INDEX + 2) - RHODOWN(INDEX + 2, INDEX + 2) + &
                RHOUP(INDEX + 3, INDEX + 3) - RHODOWN(INDEX + 3, INDEX + 3)

           INDEX = INDEX + 3

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) + &
                RHOUP(INDEX + 2, INDEX + 2) - RHODOWN(INDEX + 2, INDEX + 2) + &
                RHOUP(INDEX + 3, INDEX + 3) - RHODOWN(INDEX + 3, INDEX + 3) + &
                RHOUP(INDEX + 4, INDEX + 4) - RHODOWN(INDEX + 4, INDEX + 4) + &
                RHOUP(INDEX + 5, INDEX + 5) - RHODOWN(INDEX + 5, INDEX + 5)

           INDEX = INDEX + 5

        CASE("pf")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) + &
                RHOUP(INDEX + 2, INDEX + 2) - RHODOWN(INDEX + 2, INDEX + 2) + &
                RHOUP(INDEX + 3, INDEX + 3) - RHODOWN(INDEX + 3, INDEX + 3)

           INDEX = INDEX + 3

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) + &
                RHOUP(INDEX + 2, INDEX + 2) - RHODOWN(INDEX + 2, INDEX + 2) + &
                RHOUP(INDEX + 3, INDEX + 3) - RHODOWN(INDEX + 3, INDEX + 3) + &
                RHOUP(INDEX + 4, INDEX + 4) - RHODOWN(INDEX + 4, INDEX + 4) + &
                RHOUP(INDEX + 5, INDEX + 5) - RHODOWN(INDEX + 5, INDEX + 5) + &
                RHOUP(INDEX + 6, INDEX + 6) - RHODOWN(INDEX + 6, INDEX + 6) + &
                RHOUP(INDEX + 7, INDEX + 7) - RHODOWN(INDEX + 7, INDEX + 7)

           INDEX = INDEX + 7

        CASE("df")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) + &
                RHOUP(INDEX + 2, INDEX + 2) - RHODOWN(INDEX + 2, INDEX + 2) + &
                RHOUP(INDEX + 3, INDEX + 3) - RHODOWN(INDEX + 3, INDEX + 3) + &
                RHOUP(INDEX + 4, INDEX + 4) - RHODOWN(INDEX + 4, INDEX + 4) + &
                RHOUP(INDEX + 5, INDEX + 5) - RHODOWN(INDEX + 5, INDEX + 5)

           INDEX = INDEX + 5

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) + &
                RHOUP(INDEX + 2, INDEX + 2) - RHODOWN(INDEX + 2, INDEX + 2) + &
                RHOUP(INDEX + 3, INDEX + 3) - RHODOWN(INDEX + 3, INDEX + 3) + &
                RHOUP(INDEX + 4, INDEX + 4) - RHODOWN(INDEX + 4, INDEX + 4) + &
                RHOUP(INDEX + 5, INDEX + 5) - RHODOWN(INDEX + 5, INDEX + 5) + &
                RHOUP(INDEX + 6, INDEX + 6) - RHODOWN(INDEX + 6, INDEX + 6) + &
                RHOUP(INDEX + 7, INDEX + 7) - RHODOWN(INDEX + 7, INDEX + 7)

           INDEX = INDEX + 7

        CASE("spd")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) 

           INDEX = INDEX + 1

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) + &
                RHOUP(INDEX + 2, INDEX + 2) - RHODOWN(INDEX + 2, INDEX + 2) + &
                RHOUP(INDEX + 3, INDEX + 3) - RHODOWN(INDEX + 3, INDEX + 3)

           INDEX = INDEX + 3

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) + &
                RHOUP(INDEX + 2, INDEX + 2) - RHODOWN(INDEX + 2, INDEX + 2) + &
                RHOUP(INDEX + 3, INDEX + 3) - RHODOWN(INDEX + 3, INDEX + 3) + &
                RHOUP(INDEX + 4, INDEX + 4) - RHODOWN(INDEX + 4, INDEX + 4) + &
                RHOUP(INDEX + 5, INDEX + 5) - RHODOWN(INDEX + 5, INDEX + 5)

           INDEX = INDEX + 5

        CASE("spf")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) 

           INDEX = INDEX + 1

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) + &
                RHOUP(INDEX + 2, INDEX + 2) - RHODOWN(INDEX + 2, INDEX + 2) + &
                RHOUP(INDEX + 3, INDEX + 3) - RHODOWN(INDEX + 3, INDEX + 3)

           INDEX = INDEX + 3

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) + &
                RHOUP(INDEX + 2, INDEX + 2) - RHODOWN(INDEX + 2, INDEX + 2) + &
                RHOUP(INDEX + 3, INDEX + 3) - RHODOWN(INDEX + 3, INDEX + 3) + &
                RHOUP(INDEX + 4, INDEX + 4) - RHODOWN(INDEX + 4, INDEX + 4) + &
                RHOUP(INDEX + 5, INDEX + 5) - RHODOWN(INDEX + 5, INDEX + 5) + &
                RHOUP(INDEX + 6, INDEX + 6) - RHODOWN(INDEX + 6, INDEX + 6) + &
                RHOUP(INDEX + 7, INDEX + 7) - RHODOWN(INDEX + 7, INDEX + 7)

           INDEX = INDEX + 7

        CASE("sdf")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) 

           INDEX = INDEX + 1

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) + &
                RHOUP(INDEX + 2, INDEX + 2) - RHODOWN(INDEX + 2, INDEX + 2) + &
                RHOUP(INDEX + 3, INDEX + 3) - RHODOWN(INDEX + 3, INDEX + 3) + &
                RHOUP(INDEX + 4, INDEX + 4) - RHODOWN(INDEX + 4, INDEX + 4) + &
                RHOUP(INDEX + 5, INDEX + 5) - RHODOWN(INDEX + 5, INDEX + 5)

           INDEX = INDEX + 5

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) + &
                RHOUP(INDEX + 2, INDEX + 2) - RHODOWN(INDEX + 2, INDEX + 2) + &
                RHOUP(INDEX + 3, INDEX + 3) - RHODOWN(INDEX + 3, INDEX + 3) + &
                RHOUP(INDEX + 4, INDEX + 4) - RHODOWN(INDEX + 4, INDEX + 4) + &
                RHOUP(INDEX + 5, INDEX + 5) - RHODOWN(INDEX + 5, INDEX + 5) + &
                RHOUP(INDEX + 6, INDEX + 6) - RHODOWN(INDEX + 6, INDEX + 6) + &
                RHOUP(INDEX + 7, INDEX + 7) - RHODOWN(INDEX + 7, INDEX + 7)

           INDEX = INDEX + 7


        CASE("pdf")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) + &
                RHOUP(INDEX + 2, INDEX + 2) - RHODOWN(INDEX + 2, INDEX + 2) + &
                RHOUP(INDEX + 3, INDEX + 3) - RHODOWN(INDEX + 3, INDEX + 3)

           INDEX = INDEX + 3

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) + &
                RHOUP(INDEX + 2, INDEX + 2) - RHODOWN(INDEX + 2, INDEX + 2) + &
                RHOUP(INDEX + 3, INDEX + 3) - RHODOWN(INDEX + 3, INDEX + 3) + &
                RHOUP(INDEX + 4, INDEX + 4) - RHODOWN(INDEX + 4, INDEX + 4) + &
                RHOUP(INDEX + 5, INDEX + 5) - RHODOWN(INDEX + 5, INDEX + 5)

           INDEX = INDEX + 5

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) + &
                RHOUP(INDEX + 2, INDEX + 2) - RHODOWN(INDEX + 2, INDEX + 2) + &
                RHOUP(INDEX + 3, INDEX + 3) - RHODOWN(INDEX + 3, INDEX + 3) + &
                RHOUP(INDEX + 4, INDEX + 4) - RHODOWN(INDEX + 4, INDEX + 4) + &
                RHOUP(INDEX + 5, INDEX + 5) - RHODOWN(INDEX + 5, INDEX + 5) + &
                RHOUP(INDEX + 6, INDEX + 6) - RHODOWN(INDEX + 6, INDEX + 6) + &
                RHOUP(INDEX + 7, INDEX + 7) - RHODOWN(INDEX + 7, INDEX + 7)

           INDEX = INDEX + 7

        CASE("spdf") 

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) 

           INDEX = INDEX + 1

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) + &
                RHOUP(INDEX + 2, INDEX + 2) - RHODOWN(INDEX + 2, INDEX + 2) + &
                RHOUP(INDEX + 3, INDEX + 3) - RHODOWN(INDEX + 3, INDEX + 3)

           INDEX = INDEX + 3

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) + &
                RHOUP(INDEX + 2, INDEX + 2) - RHODOWN(INDEX + 2, INDEX + 2) + &
                RHOUP(INDEX + 3, INDEX + 3) - RHODOWN(INDEX + 3, INDEX + 3) + &
                RHOUP(INDEX + 4, INDEX + 4) - RHODOWN(INDEX + 4, INDEX + 4) + &
                RHOUP(INDEX + 5, INDEX + 5) - RHODOWN(INDEX + 5, INDEX + 5)

           INDEX = INDEX + 5

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) =  &
                RHOUP(INDEX + 1, INDEX + 1) - RHODOWN(INDEX + 1, INDEX + 1) + &
                RHOUP(INDEX + 2, INDEX + 2) - RHODOWN(INDEX + 2, INDEX + 2) + &
                RHOUP(INDEX + 3, INDEX + 3) - RHODOWN(INDEX + 3, INDEX + 3) + &
                RHOUP(INDEX + 4, INDEX + 4) - RHODOWN(INDEX + 4, INDEX + 4) + &
                RHOUP(INDEX + 5, INDEX + 5) - RHODOWN(INDEX + 5, INDEX + 5) + &
                RHOUP(INDEX + 6, INDEX + 6) - RHODOWN(INDEX + 6, INDEX + 6) + &
                RHOUP(INDEX + 7, INDEX + 7) - RHODOWN(INDEX + 7, INDEX + 7)

           INDEX = INDEX + 7

        END SELECT

     ENDDO

  ELSEIF (BASISTYPE .EQ. "NONORTHO") THEN

     ! Mulliken spin densities in non-orthogonal basis:

     ! m = n_up - n_down

     ! n_sigma = partial_trace(rho_sigma S)

     SPINLIST = ZERO

     DO I = 1, HDIM
        DO J = 1, HDIM

           SPINLIST(I) = SPINLIST(I) + (RHOUP(J,I) - RHODOWN(J,I))*SMAT(J,I)

        ENDDO
     ENDDO

     DO I = 1, NATS

        SELECT CASE(BASIS(ELEMPOINTER(I)))

        CASE("s")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1)

           INDEX = INDEX + 1

        CASE("p")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1) + &
                SPINLIST(INDEX + 2) + SPINLIST(INDEX + 3)

           INDEX = INDEX + 3

        CASE("d")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1) + &
                SPINLIST(INDEX + 2) + SPINLIST(INDEX + 3) + &
                SPINLIST(INDEX + 4) + SPINLIST(INDEX + 5)

           INDEX = INDEX + 5

        CASE("f")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1) + &
                SPINLIST(INDEX + 2) + SPINLIST(INDEX + 3) + &
                SPINLIST(INDEX + 4) + SPINLIST(INDEX + 5)  + &
                SPINLIST(INDEX + 6) + SPINLIST(INDEX + 7)

           INDEX = INDEX + 7

        CASE("sp")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1)

           INDEX = INDEX + 1

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1) + &
                SPINLIST(INDEX + 2) + SPINLIST(INDEX + 3)

           INDEX = INDEX + 3

        CASE("sd")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1)

           INDEX = INDEX + 1

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1) + &
                SPINLIST(INDEX + 2) + SPINLIST(INDEX + 3) + &
                SPINLIST(INDEX + 4) + SPINLIST(INDEX + 5)

           INDEX = INDEX + 5

        CASE("sf")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1)

           INDEX = INDEX + 1

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1) + &
                SPINLIST(INDEX + 2) + SPINLIST(INDEX + 3) + &
                SPINLIST(INDEX + 4) + SPINLIST(INDEX + 5)  + &
                SPINLIST(INDEX + 6) + SPINLIST(INDEX + 7)

           INDEX = INDEX + 7           


        CASE("pd")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1) + &
                SPINLIST(INDEX + 2) + SPINLIST(INDEX + 3)

           INDEX = INDEX + 3           

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1) + &
                SPINLIST(INDEX + 2) + SPINLIST(INDEX + 3) + &
                SPINLIST(INDEX + 4) + SPINLIST(INDEX + 5)

           INDEX = INDEX + 5

        CASE("pf")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1) + &
                SPINLIST(INDEX + 2) + SPINLIST(INDEX + 3)

           INDEX = INDEX + 3 

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1) + &
                SPINLIST(INDEX + 2) + SPINLIST(INDEX + 3) + &
                SPINLIST(INDEX + 4) + SPINLIST(INDEX + 5)  + &
                SPINLIST(INDEX + 6) + SPINLIST(INDEX + 7)

           INDEX = INDEX + 7                  

        CASE("df")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1) + &
                SPINLIST(INDEX + 2) + SPINLIST(INDEX + 3) + &
                SPINLIST(INDEX + 4) + SPINLIST(INDEX + 5)

           INDEX = INDEX + 5

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1) + &
                SPINLIST(INDEX + 2) + SPINLIST(INDEX + 3) + &
                SPINLIST(INDEX + 4) + SPINLIST(INDEX + 5)  + &
                SPINLIST(INDEX + 6) + SPINLIST(INDEX + 7)

           INDEX = INDEX + 7     


        CASE("spd")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1)

           INDEX = INDEX + 1

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1) + &
                SPINLIST(INDEX + 2) + SPINLIST(INDEX + 3)

           INDEX = INDEX + 3

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1) + &
                SPINLIST(INDEX + 2) + SPINLIST(INDEX + 3) + &
                SPINLIST(INDEX + 4) + SPINLIST(INDEX + 5)

           INDEX = INDEX + 5

        CASE("spf")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1)

           INDEX = INDEX + 1

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1) + &
                SPINLIST(INDEX + 2) + SPINLIST(INDEX + 3)

           INDEX = INDEX + 3

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1) + &
                SPINLIST(INDEX + 2) + SPINLIST(INDEX + 3) + &
                SPINLIST(INDEX + 4) + SPINLIST(INDEX + 5)  + &
                SPINLIST(INDEX + 6) + SPINLIST(INDEX + 7)

           INDEX = INDEX + 7                

        CASE("sdf")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1)

           INDEX = INDEX + 1

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1) + &
                SPINLIST(INDEX + 2) + SPINLIST(INDEX + 3) + &
                SPINLIST(INDEX + 4) + SPINLIST(INDEX + 5)

           INDEX = INDEX + 5

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1) + &
                SPINLIST(INDEX + 2) + SPINLIST(INDEX + 3) + &
                SPINLIST(INDEX + 4) + SPINLIST(INDEX + 5)  + &
                SPINLIST(INDEX + 6) + SPINLIST(INDEX + 7)

           INDEX = INDEX + 7       

        CASE("pdf")

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1) + &
                SPINLIST(INDEX + 2) + SPINLIST(INDEX + 3)

           INDEX = INDEX + 3 

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1) + &
                SPINLIST(INDEX + 2) + SPINLIST(INDEX + 3) + &
                SPINLIST(INDEX + 4) + SPINLIST(INDEX + 5)

           INDEX = INDEX + 5

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1) + &
                SPINLIST(INDEX + 2) + SPINLIST(INDEX + 3) + &
                SPINLIST(INDEX + 4) + SPINLIST(INDEX + 5)  + &
                SPINLIST(INDEX + 6) + SPINLIST(INDEX + 7)

           INDEX = INDEX + 7    

        CASE("spdf") 

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1)

           INDEX = INDEX + 1

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1) + &
                SPINLIST(INDEX + 2) + SPINLIST(INDEX + 3)

           INDEX = INDEX + 3 

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1) + &
                SPINLIST(INDEX + 2) + SPINLIST(INDEX + 3) + &
                SPINLIST(INDEX + 4) + SPINLIST(INDEX + 5)

           INDEX = INDEX + 5

           DINDEX = DINDEX + 1

           DELTASPIN(DINDEX) = SPINLIST(INDEX + 1) + &
                SPINLIST(INDEX + 2) + SPINLIST(INDEX + 3) + &
                SPINLIST(INDEX + 4) + SPINLIST(INDEX + 5)  + &
                SPINLIST(INDEX + 6) + SPINLIST(INDEX + 7)

           INDEX = INDEX + 7   

        END SELECT

     ENDDO

  ENDIF

  RETURN

END SUBROUTINE GETDELTASPIN
