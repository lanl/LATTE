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

SUBROUTINE GENHONSITE

  USE CONSTANTS_MOD
  USE SETUPARRAY

  IMPLICIT NONE

  INTEGER :: I, J, K, INDEX, SUBI

  ALLOCATE(H_ONSITE(HDIM))
  
  INDEX = 0

  ! Build diagonal elements

  DO I = 1, NATS

     K = ELEMPOINTER(I)

     SELECT CASE(BASIS(K))

     CASE("s")

        INDEX = INDEX + 1
        H_ONSITE(INDEX) = HES(K)

     CASE("p")

        DO SUBI = 1, 3
           INDEX = INDEX + 1
           H_ONSITE(INDEX) = HEP(K)
        ENDDO

     CASE("d")

        DO SUBI = 1, 5
           INDEX = INDEX + 1
           H_ONSITE(INDEX) = HED(K)
        ENDDO

     CASE("f")

        DO SUBI = 1, 7
           INDEX = INDEX + 1
           H_ONSITE(INDEX) = HEF(K)
        ENDDO

     CASE("sp")

        DO SUBI = 1, 4

           INDEX = INDEX + 1
           IF (SUBI .EQ. 1) THEN
              H_ONSITE(INDEX) = HES(K)
           ELSE
              H_ONSITE(INDEX) = HEP(K)
           ENDIF

        ENDDO

     CASE("sd")

        DO SUBI = 1, 6

           INDEX = INDEX + 1
           IF (SUBI .EQ. 1) THEN
              H_ONSITE(INDEX) = HES(K)
           ELSE
              H_ONSITE(INDEX) = HED(K)
           ENDIF

        ENDDO

     CASE("sf")

        DO SUBI = 1, 8

           INDEX = INDEX + 1
           IF (SUBI .EQ. 1) THEN
              H_ONSITE(INDEX) = HES(K)
           ELSE
              H_ONSITE(INDEX) = HEF(K)
           ENDIF

        ENDDO

     CASE("pd")

        DO SUBI = 1, 8

           INDEX = INDEX + 1
           IF (SUBI .LE. 3) THEN
              H_ONSITE(INDEX) = HEP(K)
           ELSE
              H_ONSITE(INDEX) = HED(K)
           ENDIF

        ENDDO

     CASE("pf")

        DO SUBI = 1, 10

           INDEX = INDEX + 1
           IF (SUBI .LE. 3) THEN
              H_ONSITE(INDEX) = HEP(K)
           ELSE
              H_ONSITE(INDEX) = HEF(K)
           ENDIF

        ENDDO

     CASE("df")

        DO SUBI = 1, 12

           INDEX = INDEX + 1
           IF (SUBI .LE. 5) THEN
              H_ONSITE(INDEX) = HED(K)
           ELSE
              H_ONSITE(INDEX) = HEF(K)
           ENDIF

        ENDDO

     CASE("spd")

        DO SUBI = 1, 9

           INDEX = INDEX + 1
           IF (SUBI .EQ. 1) THEN
              H_ONSITE(INDEX) = HES(K)
           ELSEIF (SUBI .GT. 1 .AND. SUBI .LE. 4) THEN
              H_ONSITE(INDEX) = HEP(K)
           ELSE
              H_ONSITE(INDEX) = HED(K)
           ENDIF

        ENDDO

     CASE("spf")

        DO SUBI = 1, 11

           INDEX = INDEX + 1
           IF (SUBI .EQ. 1) THEN
              H_ONSITE(INDEX) = HES(K)
           ELSEIF (SUBI .GT. 1 .AND. SUBI .LE. 4) THEN
              H_ONSITE(INDEX) = HEP(K)
           ELSE
              H_ONSITE(INDEX) = HEF(K)
           ENDIF

        ENDDO

     CASE("sdf")

        DO SUBI = 1, 13

           INDEX = INDEX + 1
           IF (SUBI .EQ. 1) THEN
              H_ONSITE(INDEX) = HES(K)
           ELSEIF (SUBI .GT. 1 .AND. SUBI .LE. 6) THEN
              H_ONSITE(INDEX) = HED(K)
           ELSE
              H_ONSITE(INDEX) = HEF(K)
           ENDIF

        ENDDO


     CASE("pdf")

        DO SUBI = 1, 15

           INDEX = INDEX + 1
           IF (SUBI .LE. 3) THEN
              H_ONSITE(INDEX) = HEP(K)
           ELSEIF (SUBI .GT. 3 .AND. SUBI .LE. 8) THEN
              H_ONSITE(INDEX) = HED(K)
           ELSE
              H_ONSITE(INDEX) = HEF(K)
           ENDIF

        ENDDO

     CASE("spdf")

        DO SUBI = 1, 16

           INDEX = INDEX + 1
           IF (SUBI .EQ. 1) THEN
              H_ONSITE(INDEX) = HES(K)
           ELSEIF (SUBI .GT. 1 .AND. SUBI .LE. 4) THEN
              H_ONSITE(INDEX) = HEP(K)
           ELSEIF (SUBI .GT. 4 .AND. SUBI .LE. 9) THEN
              H_ONSITE(INDEX) = HED(K)
           ELSE
              H_ONSITE(INDEX) = HEF(K)
           ENDIF

        ENDDO

     END SELECT

  ENDDO
 
  RETURN

END SUBROUTINE GENHONSITE
  
