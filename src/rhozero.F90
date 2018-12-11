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

SUBROUTINE RHOZERO

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE SPINARRAY
  USE KSPACEARRAY
  USE RESTARTARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, K, SUBI, INDEX
  REAL(LATTEPREC) :: STOT, PTOT, NUMPERORB
  REAL(LATTEPREC) :: SUP, SDOWN, PUP, PDOWN, DUP, DDOWN, FUP, FDOWN
  REAL(LATTEPREC), PARAMETER :: SPINMAXS = ONE, SPINMAXP = THREE
  REAL(LATTEPREC), PARAMETER :: SPINMAXD = FIVE, SPINMAXF = SEVEN
  IF (EXISTERROR) RETURN

  ESPIN_ZERO = ZERO

  IF (SPINON .EQ. 0) THEN

     ! No spins

     ALLOCATE(BOZERO(HDIM))

     BOZERO = ZERO

     INDEX = 0

     DO I = 1, NATS

        SELECT CASE(BASIS(ELEMPOINTER(I)))

        CASE("s")

           BOZERO(INDEX + 1) = ATOCC(ELEMPOINTER(I))
           INDEX = INDEX + 1

        CASE("p")

           ! Distribute the electrons evenly between the orbitals

           NUMPERORB = ATOCC(ELEMPOINTER(I))/THREE
           BOZERO(INDEX + 1) = NUMPERORB
           BOZERO(INDEX + 2) = NUMPERORB
           BOZERO(INDEX + 3) = NUMPERORB

           INDEX = INDEX + 3


        CASE("d")

           NUMPERORB = ATOCC(ELEMPOINTER(I))/FIVE
           BOZERO(INDEX + 1) = NUMPERORB
           BOZERO(INDEX + 2) = NUMPERORB
           BOZERO(INDEX + 3) = NUMPERORB
           BOZERO(INDEX + 4) = NUMPERORB
           BOZERO(INDEX + 5) = NUMPERORB

           INDEX = INDEX + 5

        CASE("f")

           NUMPERORB = ATOCC(ELEMPOINTER(I))/SEVEN
           BOZERO(INDEX + 1) = NUMPERORB
           BOZERO(INDEX + 2) = NUMPERORB
           BOZERO(INDEX + 3) = NUMPERORB
           BOZERO(INDEX + 4) = NUMPERORB
           BOZERO(INDEX + 5) = NUMPERORB
           BOZERO(INDEX + 6) = NUMPERORB
           BOZERO(INDEX + 7) = NUMPERORB

           INDEX = INDEX + 7

        CASE("sp")

           ! First fill s, then p

           IF (ATOCC(ELEMPOINTER(I)) .LE. TWO) THEN

              BOZERO(INDEX + 1) = ATOCC(ELEMPOINTER(I))

              INDEX = INDEX + 1

              BOZERO(INDEX + 1) = ZERO
              BOZERO(INDEX + 2) = ZERO
              BOZERO(INDEX + 3) = ZERO

              INDEX = INDEX + 3

           ELSE

              BOZERO(INDEX + 1) = TWO

              INDEX = INDEX + 1

              NUMPERORB = (ATOCC(ELEMPOINTER(I)) - TWO)/THREE

              BOZERO(INDEX + 1) = NUMPERORB
              BOZERO(INDEX + 2) = NUMPERORB
              BOZERO(INDEX + 3) = NUMPERORB

              INDEX = INDEX + 3

           ENDIF

        CASE("sd")

           BOZERO(INDEX + 1) = TWO

           INDEX = INDEX + 1

           NUMPERORB = (ATOCC(ELEMPOINTER(I)) - TWO)/FIVE

           BOZERO(INDEX + 1) = NUMPERORB
           BOZERO(INDEX + 2) = NUMPERORB
           BOZERO(INDEX + 3) = NUMPERORB
           BOZERO(INDEX + 4) = NUMPERORB
           BOZERO(INDEX + 5) = NUMPERORB

           INDEX = INDEX + 5

        CASE("sf")

           ! Caution - f electron system!

           BOZERO(INDEX + 1) = TWO

           INDEX = INDEX + 1

           NUMPERORB = (ATOCC(ELEMPOINTER(I)) - TWO)/SEVEN

           BOZERO(INDEX + 1) = NUMPERORB
           BOZERO(INDEX + 2) = NUMPERORB
           BOZERO(INDEX + 3) = NUMPERORB
           BOZERO(INDEX + 4) = NUMPERORB
           BOZERO(INDEX + 5) = NUMPERORB
           BOZERO(INDEX + 6) = NUMPERORB
           BOZERO(INDEX + 7) = NUMPERORB

           INDEX = INDEX + 7

        CASE("pd")

           ! Fill d, then p

           NUMPERORB = (ATOCC(ELEMPOINTER(I)) - TEN)/THREE

           BOZERO(INDEX + 1) = NUMPERORB
           BOZERO(INDEX + 2) = NUMPERORB
           BOZERO(INDEX + 3) = NUMPERORB

           INDEX = INDEX + 3

           BOZERO(INDEX + 1) = TWO
           BOZERO(INDEX + 2) = TWO
           BOZERO(INDEX + 3) = TWO
           BOZERO(INDEX + 4) = TWO
           BOZERO(INDEX + 5) = TWO

           INDEX = INDEX + 5

        CASE("pf")

           ! Fill f, then p, I guess

           NUMPERORB = (ATOCC(ELEMPOINTER(I)) - FOURTEEN)/THREE

           BOZERO(INDEX + 1) = NUMPERORB
           BOZERO(INDEX + 2) = NUMPERORB
           BOZERO(INDEX + 3) = NUMPERORB

           INDEX = INDEX + 3

           BOZERO(INDEX + 1) = TWO
           BOZERO(INDEX + 2) = TWO
           BOZERO(INDEX + 3) = TWO
           BOZERO(INDEX + 4) = TWO
           BOZERO(INDEX + 5) = TWO
           BOZERO(INDEX + 6) = TWO
           BOZERO(INDEX + 7) = TWO

           INDEX = INDEX + 7


        CASE("df")

           ! For the light actinides we have 1 in the d, and the rest in the f
           ! (except Th, where we have 2d's and zero f's)

           IF (ATOCC(ELEMPOINTER(I)) .EQ. 2) THEN

              NUMPERORB = TWO/FIVE

              BOZERO(INDEX + 1) = NUMPERORB
              BOZERO(INDEX + 2) = NUMPERORB
              BOZERO(INDEX + 3) = NUMPERORB
              BOZERO(INDEX + 4) = NUMPERORB
              BOZERO(INDEX + 5) = NUMPERORB

              INDEX = INDEX + 5

              BOZERO(INDEX + 1) = ZERO
              BOZERO(INDEX + 2) = ZERO
              BOZERO(INDEX + 3) = ZERO
              BOZERO(INDEX + 4) = ZERO
              BOZERO(INDEX + 5) = ZERO
              BOZERO(INDEX + 6) = ZERO
              BOZERO(INDEX + 7) = ZERO

              INDEX = INDEX + 7

           ELSE

              NUMPERORB = ONE/FIVE

              BOZERO(INDEX + 1) = NUMPERORB
              BOZERO(INDEX + 2) = NUMPERORB
              BOZERO(INDEX + 3) = NUMPERORB
              BOZERO(INDEX + 4) = NUMPERORB
              BOZERO(INDEX + 5) = NUMPERORB

              INDEX = INDEX + 5

              NUMPERORB = (ATOCC(ELEMPOINTER(I)) - ONE)/SEVEN

              BOZERO(INDEX + 1) = NUMPERORB
              BOZERO(INDEX + 2) = NUMPERORB
              BOZERO(INDEX + 3) = NUMPERORB
              BOZERO(INDEX + 4) = NUMPERORB
              BOZERO(INDEX + 5) = NUMPERORB
              BOZERO(INDEX + 6) = NUMPERORB
              BOZERO(INDEX + 7) = NUMPERORB

              INDEX = INDEX + 7

           ENDIF

        CASE("spd")

           ! s, d, then p

           IF (ATOCC(ELEMPOINTER(I)) .LE. TWO) THEN

              BOZERO(INDEX + 1) = ATOCC(ELEMPOINTER(I))

              INDEX = INDEX + 1

              BOZERO(INDEX + 1) = ZERO
              BOZERO(INDEX + 2) = ZERO
              BOZERO(INDEX + 3) = ZERO

              INDEX = INDEX + 3

              BOZERO(INDEX + 1) = ZERO
              BOZERO(INDEX + 2) = ZERO
              BOZERO(INDEX + 3) = ZERO
              BOZERO(INDEX + 4) = ZERO
              BOZERO(INDEX + 5) = ZERO

              INDEX = INDEX + 5

           ELSE IF (ATOCC(ELEMPOINTER(I)) .GT. TWO .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWELVE) THEN

              BOZERO(INDEX + 1) = TWO

              INDEX = INDEX + 1

              NUMPERORB = (ATOCC(ELEMPOINTER(I)) - TWO)/FIVE

              BOZERO(INDEX + 1) = ZERO
              BOZERO(INDEX + 2) = ZERO
              BOZERO(INDEX + 3) = ZERO

              INDEX = INDEX + 3

              BOZERO(INDEX + 1) = NUMPERORB
              BOZERO(INDEX + 2) = NUMPERORB
              BOZERO(INDEX + 3) = NUMPERORB
              BOZERO(INDEX + 4) = NUMPERORB
              BOZERO(INDEX + 5) = NUMPERORB

              INDEX = INDEX + 5

           ELSE

              BOZERO(INDEX + 1) = TWO

              INDEX = INDEX + 1

              NUMPERORB = (ATOCC(ELEMPOINTER(I)) - TWELVE)/THREE

              BOZERO(INDEX + 1) = NUMPERORB
              BOZERO(INDEX + 2) = NUMPERORB
              BOZERO(INDEX + 3) = NUMPERORB

              INDEX = INDEX + 3

              BOZERO(INDEX + 1) = TWO
              BOZERO(INDEX + 2) = TWO
              BOZERO(INDEX + 3) = TWO
              BOZERO(INDEX + 4) = TWO
              BOZERO(INDEX + 5) = TWO

              INDEX = INDEX + 5

           ENDIF

        CASE("spf")

           ! s, f, then p

           IF (ATOCC(ELEMPOINTER(I)) .LE. TWO) THEN

              BOZERO(INDEX + 1) = ATOCC(ELEMPOINTER(I))

              INDEX = INDEX + 1

              BOZERO(INDEX + 1) = ZERO
              BOZERO(INDEX + 2) = ZERO
              BOZERO(INDEX + 3) = ZERO

              INDEX = INDEX + 3

              BOZERO(INDEX + 1) = ZERO
              BOZERO(INDEX + 2) = ZERO
              BOZERO(INDEX + 3) = ZERO
              BOZERO(INDEX + 4) = ZERO
              BOZERO(INDEX + 5) = ZERO
              BOZERO(INDEX + 6) = ZERO
              BOZERO(INDEX + 7) = ZERO

              INDEX = INDEX + 7

           ELSE IF (ATOCC(ELEMPOINTER(I)) .GT. TWO .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. SIXTEEN) THEN

              BOZERO(INDEX + 1) = TWO

              INDEX = INDEX + 1

              NUMPERORB = (ATOCC(ELEMPOINTER(I)) - TWO)/SEVEN

              BOZERO(INDEX + 1) = ZERO
              BOZERO(INDEX + 2) = ZERO
              BOZERO(INDEX + 3) = ZERO

              INDEX = INDEX + 3

              BOZERO(INDEX + 1) = NUMPERORB
              BOZERO(INDEX + 2) = NUMPERORB
              BOZERO(INDEX + 3) = NUMPERORB
              BOZERO(INDEX + 4) = NUMPERORB
              BOZERO(INDEX + 5) = NUMPERORB
              BOZERO(INDEX + 6) = NUMPERORB
              BOZERO(INDEX + 7) = NUMPERORB

              INDEX = INDEX + 7

           ELSE

              BOZERO(INDEX + 1) = TWO

              INDEX = INDEX + 1

              NUMPERORB = (ATOCC(ELEMPOINTER(I)) - SIXTEEN)/THREE

              BOZERO(INDEX + 1) = NUMPERORB
              BOZERO(INDEX + 2) = NUMPERORB
              BOZERO(INDEX + 3) = NUMPERORB

              INDEX = INDEX + 3

              BOZERO(INDEX + 1) = TWO
              BOZERO(INDEX + 2) = TWO
              BOZERO(INDEX + 3) = TWO
              BOZERO(INDEX + 4) = TWO
              BOZERO(INDEX + 5) = TWO
              BOZERO(INDEX + 6) = TWO
              BOZERO(INDEX + 7) = TWO

              INDEX = INDEX + 7

           ENDIF

        CASE("sdf")

           ! Let's build this for the light actinides

           ! 7s has 2 electrons, 6d2, 5f0 for Th, and 6d1, 5f? for the others

           IF (ATOCC(ELEMPOINTER(I)) .GE. THREE .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. FOUR) THEN

              BOZERO(INDEX + 1) = TWO

              INDEX = INDEX + 1

              NUMPERORB = (ATOCC(ELEMPOINTER(I)) - TWO)/FIVE

              BOZERO(INDEX + 1) = NUMPERORB
              BOZERO(INDEX + 2) = NUMPERORB
              BOZERO(INDEX + 3) = NUMPERORB
              BOZERO(INDEX + 4) = NUMPERORB
              BOZERO(INDEX + 5) = NUMPERORB

              INDEX = INDEX + 5

              BOZERO(INDEX + 1) = ZERO
              BOZERO(INDEX + 2) = ZERO
              BOZERO(INDEX + 3) = ZERO
              BOZERO(INDEX + 4) = ZERO
              BOZERO(INDEX + 5) = ZERO
              BOZERO(INDEX + 6) = ZERO
              BOZERO(INDEX + 7) = ZERO

              INDEX = INDEX + 7

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. FOUR) THEN

              BOZERO(INDEX + 1) = TWO

              INDEX = INDEX + 1

              !              NUMPERORB = ONE/FIVE
              NUMPERORB = ZERO

              BOZERO(INDEX + 1) = NUMPERORB
              BOZERO(INDEX + 2) = NUMPERORB
              BOZERO(INDEX + 3) = NUMPERORB
              BOZERO(INDEX + 4) = NUMPERORB
              BOZERO(INDEX + 5) = NUMPERORB

              INDEX = INDEX + 5

              !              NUMPERORB = (ATOCC(ELEMPOINTER(I)) - THREE)/SEVEN
              NUMPERORB = (ATOCC(ELEMPOINTER(I)) - TWO)/SEVEN

              BOZERO(INDEX + 1) = NUMPERORB
              BOZERO(INDEX + 2) = NUMPERORB
              BOZERO(INDEX + 3) = NUMPERORB
              BOZERO(INDEX + 4) = NUMPERORB
              BOZERO(INDEX + 5) = NUMPERORB
              BOZERO(INDEX + 6) = NUMPERORB
              BOZERO(INDEX + 7) = NUMPERORB

              INDEX = INDEX + 7

           ELSE

              CALL ERRORS("rhozero","Check the number of electrons &
                   & you're using with the sdf basis")

           ENDIF

        CASE("pdf")

           ! I suppose we fill the f, then the d, and finally the p

           IF (ATOCC(ELEMPOINTER(I)) .LE. FOURTEEN) THEN

              NUMPERORB = ATOCC(ELEMPOINTER(I))/SEVEN

              BOZERO(INDEX + 1) = NUMPERORB
              BOZERO(INDEX + 2) = NUMPERORB
              BOZERO(INDEX + 3) = NUMPERORB
              BOZERO(INDEX + 4) = NUMPERORB
              BOZERO(INDEX + 5) = NUMPERORB
              BOZERO(INDEX + 6) = NUMPERORB
              BOZERO(INDEX + 7) = NUMPERORB

              INDEX = INDEX + 7

              BOZERO(INDEX + 1) = ZERO
              BOZERO(INDEX + 2) = ZERO
              BOZERO(INDEX + 3) = ZERO
              BOZERO(INDEX + 4) = ZERO
              BOZERO(INDEX + 5) = ZERO

              INDEX = INDEX + 5

              BOZERO(INDEX + 1) = ZERO
              BOZERO(INDEX + 2) = ZERO
              BOZERO(INDEX + 3) = ZERO

              INDEX = INDEX + 3

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. FOURTEEN .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWENTYFOUR) THEN

              BOZERO(INDEX + 1) = TWO
              BOZERO(INDEX + 2) = TWO
              BOZERO(INDEX + 3) = TWO
              BOZERO(INDEX + 4) = TWO
              BOZERO(INDEX + 5) = TWO
              BOZERO(INDEX + 6) = TWO
              BOZERO(INDEX + 7) = TWO

              INDEX = INDEX + 7

              NUMPERORB = (ATOCC(ELEMPOINTER(I)) - FOURTEEN)/FIVE

              BOZERO(INDEX + 1) = NUMPERORB
              BOZERO(INDEX + 2) = NUMPERORB
              BOZERO(INDEX + 3) = NUMPERORB
              BOZERO(INDEX + 4) = NUMPERORB
              BOZERO(INDEX + 5) = NUMPERORB

              INDEX = INDEX + 5

              BOZERO(INDEX + 1) = ZERO
              BOZERO(INDEX + 2) = ZERO
              BOZERO(INDEX + 3) = ZERO

              INDEX = INDEX + 3

           ELSE

              BOZERO(INDEX + 1) = TWO
              BOZERO(INDEX + 2) = TWO
              BOZERO(INDEX + 3) = TWO
              BOZERO(INDEX + 4) = TWO
              BOZERO(INDEX + 5) = TWO
              BOZERO(INDEX + 6) = TWO
              BOZERO(INDEX + 7) = TWO

              INDEX = INDEX + 7

              BOZERO(INDEX + 1) = TWO
              BOZERO(INDEX + 2) = TWO
              BOZERO(INDEX + 3) = TWO
              BOZERO(INDEX + 4) = TWO
              BOZERO(INDEX + 5) = TWO

              INDEX = INDEX + 5

              NUMPERORB = (ATOCC(ELEMPOINTER(I)) - TWENTYFOUR)/THREE

              BOZERO(INDEX + 1) = NUMPERORB
              BOZERO(INDEX + 2) = NUMPERORB
              BOZERO(INDEX + 3) = NUMPERORB

              INDEX = INDEX + 3

           ENDIF

        CASE("spdf")

           ! s, then f, then d, then p...

           IF (ATOCC(ELEMPOINTER(I)) .LE. TWO) THEN

              NUMPERORB = ATOCC(ELEMPOINTER(I))

              BOZERO(INDEX + 1) = NUMPERORB

              BOZERO(INDEX + 1) = ZERO
              BOZERO(INDEX + 2) = ZERO
              BOZERO(INDEX + 3) = ZERO
              BOZERO(INDEX + 4) = ZERO
              BOZERO(INDEX + 5) = ZERO
              BOZERO(INDEX + 6) = ZERO
              BOZERO(INDEX + 7) = ZERO

              INDEX = INDEX + 7

              BOZERO(INDEX + 1) = ZERO
              BOZERO(INDEX + 2) = ZERO
              BOZERO(INDEX + 3) = ZERO
              BOZERO(INDEX + 4) = ZERO
              BOZERO(INDEX + 5) = ZERO

              INDEX = INDEX + 5

              BOZERO(INDEX + 1) = ZERO
              BOZERO(INDEX + 2) = ZERO
              BOZERO(INDEX + 3) = ZERO

              INDEX = INDEX + 3

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. TWO .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. SIXTEEN) THEN

              BOZERO(INDEX + 1) = TWO

              NUMPERORB = (ATOCC(ELEMPOINTER(I)) - TWO)/SEVEN

              BOZERO(INDEX + 1) = NUMPERORB
              BOZERO(INDEX + 2) = NUMPERORB
              BOZERO(INDEX + 3) = NUMPERORB
              BOZERO(INDEX + 4) = NUMPERORB
              BOZERO(INDEX + 5) = NUMPERORB
              BOZERO(INDEX + 6) = NUMPERORB
              BOZERO(INDEX + 7) = NUMPERORB

              INDEX = INDEX + 7

              BOZERO(INDEX + 1) = ZERO
              BOZERO(INDEX + 2) = ZERO
              BOZERO(INDEX + 3) = ZERO
              BOZERO(INDEX + 4) = ZERO
              BOZERO(INDEX + 5) = ZERO

              INDEX = INDEX + 5

              BOZERO(INDEX + 1) = ZERO
              BOZERO(INDEX + 2) = ZERO
              BOZERO(INDEX + 3) = ZERO

              INDEX = INDEX + 3

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. SIXTEEN .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWENTYSIX) THEN

              BOZERO(INDEX + 1) = TWO

              INDEX = INDEX + 1

              BOZERO(INDEX + 1) = TWO
              BOZERO(INDEX + 2) = TWO
              BOZERO(INDEX + 3) = TWO
              BOZERO(INDEX + 4) = TWO
              BOZERO(INDEX + 5) = TWO
              BOZERO(INDEX + 6) = TWO
              BOZERO(INDEX + 7) = TWO

              INDEX = INDEX + 7

              NUMPERORB = (ATOCC(ELEMPOINTER(I)) - SIXTEEN)/FIVE

              BOZERO(INDEX + 1) = NUMPERORB
              BOZERO(INDEX + 2) = NUMPERORB
              BOZERO(INDEX + 3) = NUMPERORB
              BOZERO(INDEX + 4) = NUMPERORB
              BOZERO(INDEX + 5) = NUMPERORB

              INDEX = INDEX + 5

              BOZERO(INDEX + 1) = ZERO
              BOZERO(INDEX + 2) = ZERO
              BOZERO(INDEX + 3) = ZERO

              INDEX = INDEX + 3

           ELSE

              BOZERO(INDEX + 1) = TWO

              INDEX = INDEX + 1

              BOZERO(INDEX + 1) = TWO
              BOZERO(INDEX + 2) = TWO
              BOZERO(INDEX + 3) = TWO
              BOZERO(INDEX + 4) = TWO
              BOZERO(INDEX + 5) = TWO
              BOZERO(INDEX + 6) = TWO
              BOZERO(INDEX + 7) = TWO

              INDEX = INDEX + 7

              BOZERO(INDEX + 1) = TWO
              BOZERO(INDEX + 2) = TWO
              BOZERO(INDEX + 3) = TWO
              BOZERO(INDEX + 4) = TWO
              BOZERO(INDEX + 5) = TWO

              INDEX = INDEX + 5

              NUMPERORB = (ATOCC(ELEMPOINTER(I)) - TWENTYSIX)/THREE

              BOZERO(INDEX + 1) = NUMPERORB
              BOZERO(INDEX + 2) = NUMPERORB
              BOZERO(INDEX + 3) = NUMPERORB

              INDEX = INDEX + 3

           ENDIF

        END SELECT

     ENDDO


     IF (KON .EQ. 0) THEN ! Real-space

        IF (RESTART .EQ. 0) THEN

           !
           ! Initialize the diagonal elements of BO to those for free atoms
           !

           DO I = 1, HDIM
              BO(I,I) = BOZERO(I)
           ENDDO

        ELSEIF (RESTART .EQ. 1) THEN

           DO I = 1, HDIM
              BO(I,I) = TMPBODIAG(I)
           ENDDO

           DEALLOCATE(TMPBODIAG)

        ENDIF

     ELSE ! k-space

        IF (RESTART .EQ. 0) THEN

           !
           ! Initialize the diagonal elements of BO to those for free atoms
           !

           DO K = 1, NKTOT
              DO I = 1, HDIM
                 KBO(I,I,K) = CMPLX(BOZERO(I))
              ENDDO
           ENDDO


        ELSEIF (RESTART .EQ. 1) THEN

           DO K = 1, NKTOT
              DO I = 1, HDIM
                 KBO(I,I,K) = CMPLX(TMPBODIAG(I))
              ENDDO
           ENDDO

           DEALLOCATE(TMPBODIAG)

        ENDIF

     ENDIF


  ELSEIF (SPINON .EQ. 1) THEN

     ! With spins

     ALLOCATE(RHOUPZERO(HDIM), RHODOWNZERO(HDIM))

     RHOUPZERO = ZERO
     RHODOWNZERO = ZERO

     INDEX = 0

     DO I = 1, NATS

        SELECT CASE(BASIS(ELEMPOINTER(I)))

        CASE("s")

           IF (ATOCC(ELEMPOINTER(I)) .LE. SPINMAXS) THEN

              RHOUPZERO(INDEX+1) = ATOCC(ELEMPOINTER(I))
              RHODOWNZERO(INDEX+1) = ZERO

              INDEX = INDEX + 1

           ELSE

              RHOUPZERO(INDEX+1) = ONE
              RHODOWNZERO(INDEX+1) = ATOCC(ELEMPOINTER(I)) - SPINMAXS

              INDEX = INDEX + 1

           ENDIF

        CASE("p")

           IF (ATOCC(ELEMPOINTER(I)) .LE. SPINMAXP) THEN

              PUP = ATOCC(ELEMPOINTER(I))/SPINMAXP
              PDOWN = ZERO

           ELSE

              PUP = ONE
              PDOWN = (ATOCC(ELEMPOINTER(I)) - SPINMAXP)/SPINMAXP

           ENDIF

           RHOUPZERO(INDEX + 1) = PUP
           RHOUPZERO(INDEX + 2) = PUP
           RHOUPZERO(INDEX + 3) = PUP

           RHODOWNZERO(INDEX + 1) = PDOWN
           RHODOWNZERO(INDEX + 2) = PDOWN
           RHODOWNZERO(INDEX + 3) = PDOWN

           INDEX = INDEX + 3

        CASE("d")

           IF (ATOCC(ELEMPOINTER(I)) .LE. SPINMAXD) THEN

              DUP = ATOCC(ELEMPOINTER(I))/SPINMAXD
              DDOWN = ZERO

           ELSE

              DUP = ONE
              DDOWN = (ATOCC(ELEMPOINTER(I)) - SPINMAXD)/SPINMAXD

           ENDIF

           RHOUPZERO(INDEX + 1) = DUP
           RHOUPZERO(INDEX + 2) = DUP
           RHOUPZERO(INDEX + 3) = DUP
           RHOUPZERO(INDEX + 4) = DUP
           RHOUPZERO(INDEX + 5) = DUP


           RHODOWNZERO(INDEX + 1) = DDOWN
           RHODOWNZERO(INDEX + 2) = DDOWN
           RHODOWNZERO(INDEX + 3) = DDOWN
           RHODOWNZERO(INDEX + 4) = DDOWN
           RHODOWNZERO(INDEX + 5) = DDOWN

           INDEX = INDEX + 5

        CASE("f")

           IF (ATOCC(ELEMPOINTER(I)) .LE. SPINMAXF) THEN

              FUP = ATOCC(ELEMPOINTER(I))/SPINMAXF
              FDOWN = ZERO

           ELSE

              FUP = ONE
              FDOWN = (ATOCC(ELEMPOINTER(I)) - SPINMAXF)/SPINMAXF

           ENDIF


           RHOUPZERO(INDEX + 1) = FUP
           RHOUPZERO(INDEX + 2) = FUP
           RHOUPZERO(INDEX + 3) = FUP
           RHOUPZERO(INDEX + 4) = FUP
           RHOUPZERO(INDEX + 5) = FUP
           RHOUPZERO(INDEX + 6) = FUP
           RHOUPZERO(INDEX + 7) = FUP

           RHODOWNZERO(INDEX + 1) = FDOWN
           RHODOWNZERO(INDEX + 2) = FDOWN
           RHODOWNZERO(INDEX + 3) = FDOWN
           RHODOWNZERO(INDEX + 4) = FDOWN
           RHODOWNZERO(INDEX + 5) = FDOWN
           RHODOWNZERO(INDEX + 6) = FDOWN
           RHODOWNZERO(INDEX + 7) = FDOWN

           INDEX = INDEX + 7

        CASE("sp")

           IF (ATOCC(ELEMPOINTER(I)) .LE. SPINMAXS) THEN

              SUP = ATOCC(ELEMPOINTER(I))
              SDOWN = ZERO
              PUP = ZERO
              PDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. SPINMAXS .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWO*SPINMAXS) THEN

              SUP = ONE
              SDOWN = ATOCC(ELEMPOINTER(I)) - SPINMAXS
              PUP = ZERO
              PDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. TWO*SPINMAXS .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWO*SPINMAXS + SPINMAXP) THEN


              SUP = ONE
              SDOWN = ONE
              PUP = (ATOCC(ELEMPOINTER(I)) - TWO*SPINMAXS)/SPINMAXP
              PDOWN = ZERO

           ELSE

              SUP = ONE
              SDOWN = ONE
              PUP = ONE
              PDOWN = (ATOCC(ELEMPOINTER(I)) - &
                   (TWO*SPINMAXS + SPINMAXP))/SPINMAXP

           ENDIF

           RHOUPZERO(INDEX + 1) = SUP
           RHOUPZERO(INDEX + 2) = PUP
           RHOUPZERO(INDEX + 3) = PUP
           RHOUPZERO(INDEX + 4) = PUP

           RHODOWNZERO(INDEX + 1) = SDOWN
           RHODOWNZERO(INDEX + 2) = PDOWN
           RHODOWNZERO(INDEX + 3) = PDOWN
           RHODOWNZERO(INDEX + 4) = PDOWN

           INDEX = INDEX + 4

        CASE("sd")

           IF (ATOCC(ELEMPOINTER(I)) .LE. SPINMAXS) THEN

              SUP = ATOCC(ELEMPOINTER(I))
              SDOWN = ZERO
              DUP = ZERO
              DDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. SPINMAXS .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWO*SPINMAXS) THEN

              SUP = ONE
              SDOWN = ATOCC(ELEMPOINTER(I)) - SPINMAXS
              DUP = ZERO
              DDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. TWO*SPINMAXS .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWO*SPINMAXS + SPINMAXD) THEN


              SUP = ONE
              SDOWN = ONE
              DUP = (ATOCC(ELEMPOINTER(I)) - TWO*SPINMAXS)/SPINMAXD
              DDOWN = ZERO

           ELSE

              SUP = ONE
              SDOWN = ONE
              DUP = ONE
              DDOWN = (ATOCC(ELEMPOINTER(I)) - &
                   (TWO*SPINMAXS + SPINMAXD))/SPINMAXD

           ENDIF

           RHOUPZERO(INDEX + 1) = SUP
           RHOUPZERO(INDEX + 2) = DUP
           RHOUPZERO(INDEX + 3) = DUP
           RHOUPZERO(INDEX + 4) = DUP
           RHOUPZERO(INDEX + 5) = DUP
           RHOUPZERO(INDEX + 6) = DUP


           RHODOWNZERO(INDEX + 1) = SDOWN
           RHODOWNZERO(INDEX + 2) = DDOWN
           RHODOWNZERO(INDEX + 3) = DDOWN
           RHODOWNZERO(INDEX + 4) = DDOWN
           RHODOWNZERO(INDEX + 5) = DDOWN
           RHODOWNZERO(INDEX + 6) = DDOWN

           INDEX = INDEX + 6

        CASE("sf")

           IF (ATOCC(ELEMPOINTER(I)) .LE. SPINMAXS) THEN

              SUP = ATOCC(ELEMPOINTER(I))
              SDOWN = ZERO
              FUP = ZERO
              FDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. SPINMAXS .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWO*SPINMAXS) THEN

              SUP = ONE
              SDOWN = ATOCC(ELEMPOINTER(I)) - ONE
              FUP = ZERO
              FDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. TWO*SPINMAXS .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWO*SPINMAXS + SPINMAXF) THEN


              SUP = ONE
              SDOWN = ONE
              FUP = (ATOCC(ELEMPOINTER(I)) - TWO*SPINMAXS)/SPINMAXF
              FDOWN = ZERO

           ELSE

              SUP = ONE
              SDOWN = ONE
              FUP = ONE
              FDOWN = (ATOCC(ELEMPOINTER(I)) - &
                   (TWO*SPINMAXS + SPINMAXF))/SPINMAXF

           ENDIF

           RHOUPZERO(INDEX + 1) = SUP
           RHOUPZERO(INDEX + 2) = FUP
           RHOUPZERO(INDEX + 3) = FUP
           RHOUPZERO(INDEX + 4) = FUP
           RHOUPZERO(INDEX + 5) = FUP
           RHOUPZERO(INDEX + 6) = FUP
           RHOUPZERO(INDEX + 7) = FUP
           RHOUPZERO(INDEX + 8) = FUP

           RHODOWNZERO(INDEX + 1) = SDOWN
           RHODOWNZERO(INDEX + 2) = FDOWN
           RHODOWNZERO(INDEX + 3) = FDOWN
           RHODOWNZERO(INDEX + 4) = FDOWN
           RHODOWNZERO(INDEX + 5) = FDOWN
           RHODOWNZERO(INDEX + 6) = FDOWN
           RHODOWNZERO(INDEX + 7) = FDOWN
           RHODOWNZERO(INDEX + 8) = FDOWN

           INDEX = INDEX + 8

        CASE("pd")

           ! We'll fill the d, then the p

           IF (ATOCC(ELEMPOINTER(I)) .LE. SPINMAXD) THEN

              DUP = ATOCC(ELEMPOINTER(I))/SPINMAXD
              DDOWN = ZERO
              PUP = ZERO
              PDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. SPINMAXD .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWO*SPINMAXD) THEN

              DUP = ONE
              DDOWN = (ATOCC(ELEMPOINTER(I)) - SPINMAXD)/SPINMAXD
              PUP = ZERO
              PDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. TWO*SPINMAXD .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWO*SPINMAXD + SPINMAXP) THEN


              DUP = ONE
              DDOWN = ONE
              PUP = (ATOCC(ELEMPOINTER(I)) - TWO*SPINMAXD)/SPINMAXP
              PDOWN = ZERO

           ELSE

              DUP = ONE
              DDOWN = ONE
              PUP = ONE
              PDOWN = (ATOCC(ELEMPOINTER(I)) - &
                   (TWO*SPINMAXD + SPINMAXP))/SPINMAXP

           ENDIF

           RHOUPZERO(INDEX + 1) = PUP
           RHOUPZERO(INDEX + 2) = PUP
           RHOUPZERO(INDEX + 3) = PUP
           RHOUPZERO(INDEX + 4) = DUP
           RHOUPZERO(INDEX + 5) = DUP
           RHOUPZERO(INDEX + 6) = DUP
           RHOUPZERO(INDEX + 7) = DUP
           RHOUPZERO(INDEX + 8) = DUP

           RHODOWNZERO(INDEX + 1) = PDOWN
           RHODOWNZERO(INDEX + 2) = PDOWN
           RHODOWNZERO(INDEX + 3) = PDOWN
           RHODOWNZERO(INDEX + 4) = DDOWN
           RHODOWNZERO(INDEX + 5) = DDOWN
           RHODOWNZERO(INDEX + 6) = DDOWN
           RHODOWNZERO(INDEX + 7) = DDOWN
           RHODOWNZERO(INDEX + 8) = DDOWN

           INDEX = INDEX + 8

        CASE("pf")

           ! First the f, then the p

           IF (ATOCC(ELEMPOINTER(I)) .LE. SPINMAXF) THEN

              FUP = ATOCC(ELEMPOINTER(I))/SPINMAXF
              FDOWN = ZERO
              PUP = ZERO
              PDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. SPINMAXF .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWO*SPINMAXF) THEN

              FUP = ONE
              FDOWN = (ATOCC(ELEMPOINTER(I)) - SPINMAXF)/SPINMAXF
              PUP = ZERO
              PDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. TWO*SPINMAXF .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWO*SPINMAXF + SPINMAXP) THEN


              FUP = ONE
              FDOWN = ONE
              PUP = (ATOCC(ELEMPOINTER(I)) - TWO*SPINMAXF)/SPINMAXP
              PDOWN = ZERO

           ELSE

              FUP = ONE
              FDOWN = ONE
              PUP = ONE
              PDOWN = (ATOCC(ELEMPOINTER(I)) - &
                   (TWO*SPINMAXF + SPINMAXP))/SPINMAXP

           ENDIF

           RHOUPZERO(INDEX + 1) = PUP
           RHOUPZERO(INDEX + 2) = PUP
           RHOUPZERO(INDEX + 3) = PUP
           RHOUPZERO(INDEX + 4) = FUP
           RHOUPZERO(INDEX + 5) = FUP
           RHOUPZERO(INDEX + 6) = FUP
           RHOUPZERO(INDEX + 7) = FUP
           RHOUPZERO(INDEX + 8) = FUP
           RHOUPZERO(INDEX + 9) = FUP
           RHOUPZERO(INDEX + 10) = FUP

           RHODOWNZERO(INDEX + 1) = PDOWN
           RHODOWNZERO(INDEX + 2) = PDOWN
           RHODOWNZERO(INDEX + 3) = PDOWN
           RHODOWNZERO(INDEX + 4) = FDOWN
           RHODOWNZERO(INDEX + 5) = FDOWN
           RHODOWNZERO(INDEX + 6) = FDOWN
           RHODOWNZERO(INDEX + 7) = FDOWN
           RHODOWNZERO(INDEX + 8) = FDOWN
           RHODOWNZERO(INDEX + 9) = FDOWN
           RHODOWNZERO(INDEX + 10) = FDOWN

           INDEX = INDEX + 10


        CASE("df")

           ! Let's imagine these are the light actinides

           ! Thorium

           IF (ATOCC(ELEMPOINTER(I)) .GT. FOUR) THEN
              CALL ERRORS("rhozero","The df basis is limited to 4 electrons (U)&
                   & at the moment. Modify the source to go beyond U")
           ENDIF

           IF (ATOCC(ELEMPOINTER(I)) .LE. TWO) THEN

              DUP = ATOCC(ELEMPOINTER(I))/SPINMAXD
              DDOWN = ZERO
              FUP = ZERO
              FDOWN = ZERO

           ELSE

              DUP = ONE/SPINMAXD
              DDOWN = ZERO
              FUP = (ATOCC(ELEMPOINTER(I)) - ONE)/SPINMAXF
              FDOWN = ZERO

           ENDIF

           RHOUPZERO(INDEX + 1) = DUP
           RHOUPZERO(INDEX + 2) = DUP
           RHOUPZERO(INDEX + 3) = DUP
           RHOUPZERO(INDEX + 4) = DUP
           RHOUPZERO(INDEX + 5) = DUP
           RHOUPZERO(INDEX + 6) = FUP
           RHOUPZERO(INDEX + 7) = FUP
           RHOUPZERO(INDEX + 8) = FUP
           RHOUPZERO(INDEX + 9) = FUP
           RHOUPZERO(INDEX + 10) = FUP
           RHOUPZERO(INDEX + 11) = FUP
           RHOUPZERO(INDEX + 12) = FUP

           RHODOWNZERO(INDEX + 1) = DDOWN
           RHODOWNZERO(INDEX + 2) = DDOWN
           RHODOWNZERO(INDEX + 3) = DDOWN
           RHODOWNZERO(INDEX + 4) = DDOWN
           RHODOWNZERO(INDEX + 5) = DDOWN
           RHODOWNZERO(INDEX + 6) = FDOWN
           RHODOWNZERO(INDEX + 7) = FDOWN
           RHODOWNZERO(INDEX + 8) = FDOWN
           RHODOWNZERO(INDEX + 9) = FDOWN
           RHODOWNZERO(INDEX + 10) = FDOWN
           RHODOWNZERO(INDEX + 11) = FDOWN
           RHODOWNZERO(INDEX + 12) = FDOWN

           INDEX = INDEX + 12

        CASE("spd")

           ! s, then d, then p...

           IF (ATOCC(ELEMPOINTER(I)) .LE. SPINMAXS) THEN

              SUP = ATOCC(ELEMPOINTER(I))
              SDOWN = ZERO
              DUP = ZERO
              DDOWN = ZERO
              PUP = ZERO
              PDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. SPINMAXS .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWO*SPINMAXS) THEN

              SUP = ONE
              SDOWN = ATOCC(ELEMPOINTER(I)) - SPINMAXS
              DUP = ZERO
              DDOWN = ZERO
              PUP = ZERO
              PDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. TWO*SPINMAXS .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWO*SPINMAXS + SPINMAXD) THEN


              SUP = ONE
              SDOWN = ONE
              DUP = (ATOCC(ELEMPOINTER(I)) - TWO*SPINMAXS)/SPINMAXD
              DDOWN = ZERO
              PUP = ZERO
              PDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. TWO*SPINMAXS + SPINMAXD .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWO*SPINMAXS + TWO*SPINMAXD) THEN

              SUP = ONE
              SDOWN = ONE
              DUP = ONE
              DDOWN = (ATOCC(ELEMPOINTER(I)) - &
                   (TWO*SPINMAXS + SPINMAXD))/SPINMAXD
              PUP = ZERO
              PDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. TWO*SPINMAXS + TWO*SPINMAXD &
                .AND. ATOCC(ELEMPOINTER(I)) .LE. &
                TWO*SPINMAXS + TWO*SPINMAXD + SPINMAXP ) THEN

              SUP = ONE
              SDOWN = ONE
              DUP = ONE
              DDOWN = ONE
              PUP = (ATOCC(ELEMPOINTER(I)) - &
                   (TWO*SPINMAXS + TWO*SPINMAXD))/SPINMAXP
              PDOWN = ZERO

           ELSE

              SUP = ONE
              SDOWN = ONE
              DUP = ONE
              DDOWN = ONE
              PUP = ONE
              PDOWN = (ATOCC(ELEMPOINTER(I)) - &
                   (TWO*SPINMAXS + TWO*SPINMAXD + SPINMAXP))/SPINMAXP

           ENDIF

           RHOUPZERO(INDEX + 1) = SUP
           RHOUPZERO(INDEX + 2) = PUP
           RHOUPZERO(INDEX + 3) = PUP
           RHOUPZERO(INDEX + 4) = PUP
           RHOUPZERO(INDEX + 5) = DUP
           RHOUPZERO(INDEX + 6) = DUP
           RHOUPZERO(INDEX + 7) = DUP
           RHOUPZERO(INDEX + 8) = DUP
           RHOUPZERO(INDEX + 9) = DUP

           RHODOWNZERO(INDEX + 1) = SDOWN
           RHODOWNZERO(INDEX + 2) = PDOWN
           RHODOWNZERO(INDEX + 3) = PDOWN
           RHODOWNZERO(INDEX + 4) = PDOWN
           RHODOWNZERO(INDEX + 5) = DDOWN
           RHODOWNZERO(INDEX + 6) = DDOWN
           RHODOWNZERO(INDEX + 7) = DDOWN
           RHODOWNZERO(INDEX + 8) = DDOWN
           RHODOWNZERO(INDEX + 9) = DDOWN

           INDEX = INDEX + 9


        CASE("spf")

           ! s, then f, then p

           IF (ATOCC(ELEMPOINTER(I)) .LE. SPINMAXS) THEN

              SUP = ATOCC(ELEMPOINTER(I))
              SDOWN = ZERO
              FUP = ZERO
              FDOWN = ZERO
              PUP = ZERO
              PDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. SPINMAXS .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWO*SPINMAXS) THEN

              SUP = ONE
              SDOWN = ATOCC(ELEMPOINTER(I)) - SPINMAXS
              FUP = ZERO
              FDOWN = ZERO
              PUP = ZERO
              PDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. TWO*SPINMAXS .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWO*SPINMAXS + SPINMAXF) THEN


              SUP = ONE
              SDOWN = ONE
              FUP = (ATOCC(ELEMPOINTER(I)) - TWO*SPINMAXS)/SPINMAXF
              FDOWN = ZERO
              PUP = ZERO
              PDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. TWO*SPINMAXS + SPINMAXF .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWO*SPINMAXS + TWO*SPINMAXF) THEN

              SUP = ONE
              SDOWN = ONE
              FUP = ONE
              FDOWN = (ATOCC(ELEMPOINTER(I)) - &
                   (TWO*SPINMAXS + SPINMAXF))/SPINMAXF
              PUP = ZERO
              PDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. TWO*SPINMAXS + TWO*SPINMAXF &
                .AND. ATOCC(ELEMPOINTER(I)) .LE. &
                TWO*SPINMAXS + TWO*SPINMAXF + SPINMAXP ) THEN

              SUP = ONE
              SDOWN = ONE
              FUP = ONE
              FDOWN = ONE
              PUP = (ATOCC(ELEMPOINTER(I)) - &
                   (TWO*SPINMAXS + TWO*SPINMAXF))/SPINMAXP
              PDOWN = ZERO

           ELSE

              SUP = ONE
              SDOWN = ONE
              FUP = ONE
              FDOWN = ONE
              PUP = ONE
              PDOWN = (ATOCC(ELEMPOINTER(I)) - &
                   (TWO*SPINMAXS + TWO*SPINMAXF + SPINMAXP))/SPINMAXP

           ENDIF

           RHOUPZERO(INDEX + 1) = SUP
           RHOUPZERO(INDEX + 2) = PUP
           RHOUPZERO(INDEX + 3) = PUP
           RHOUPZERO(INDEX + 4) = PUP
           RHOUPZERO(INDEX + 5) = FUP
           RHOUPZERO(INDEX + 6) = FUP
           RHOUPZERO(INDEX + 7) = FUP
           RHOUPZERO(INDEX + 8) = FUP
           RHOUPZERO(INDEX + 9) = FUP
           RHOUPZERO(INDEX + 10) = FUP
           RHOUPZERO(INDEX + 11) = FUP

           RHODOWNZERO(INDEX + 1) = SDOWN
           RHODOWNZERO(INDEX + 2) = PDOWN
           RHODOWNZERO(INDEX + 3) = PDOWN
           RHODOWNZERO(INDEX + 4) = PDOWN
           RHODOWNZERO(INDEX + 5) = FDOWN
           RHODOWNZERO(INDEX + 6) = FDOWN
           RHODOWNZERO(INDEX + 7) = FDOWN
           RHODOWNZERO(INDEX + 8) = FDOWN
           RHODOWNZERO(INDEX + 9) = FDOWN
           RHODOWNZERO(INDEX + 10) = FDOWN
           RHODOWNZERO(INDEX + 11) = FDOWN

           INDEX = INDEX + 11

        CASE("sdf")

           ! Let's do the light actinides explicitly

           ! Thorium

           IF (ATOCC(ELEMPOINTER(I)) .LT. SPINMAXS) THEN

              SUP = ATOCC(ELEMPOINTER(I))
              SDOWN = ZERO
              FUP = ZERO
              FDOWN = ZERO
              DUP = ZERO
              DDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. SPINMAXS .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWO*SPINMAXS) THEN

              SUP = ONE
              SDOWN = ATOCC(ELEMPOINTER(I)) - SPINMAXS
              FUP = ZERO
              FDOWN = ZERO
              DUP = ZERO
              DDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. TWO*SPINMAXS .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWO*SPINMAXS + TWO) THEN

              SUP = ONE
              SDOWN = ONE
              FUP = ZERO
              FDOWN = ZERO
              DUP = (ATOCC(ELEMPOINTER(I)) - TWO*SPINMAXS)/SPINMAXD
              DDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. TWO*SPINMAXS + TWO) THEN

              ! Having one electron in the 6d isn't working...
              ! Let's go with 4 in the 5 f instead

              !              SUP = ONE
              !              SDOWN = ONE
              !              FUP = (ATOCC(ELEMPOINTER(I)) - ONE)/SPINMAXF
              !              FDOWN = ZERO
              !              DUP = ONE/SPINMAXD
              !              DDOWN = ZERO

              SUP = ONE
              SDOWN = ONE
              FUP = (ATOCC(ELEMPOINTER(I)) - TWO)/SPINMAXF
              FDOWN = ZERO
              DUP = ZERO
              DDOWN = ZERO

           ENDIF

           RHOUPZERO(INDEX + 1) = SUP
           RHOUPZERO(INDEX + 2) = DUP
           RHOUPZERO(INDEX + 3) = DUP
           RHOUPZERO(INDEX + 4) = DUP
           RHOUPZERO(INDEX + 5) = DUP
           RHOUPZERO(INDEX + 6) = DUP
           RHOUPZERO(INDEX + 7) = FUP
           RHOUPZERO(INDEX + 8) = FUP
           RHOUPZERO(INDEX + 9) = FUP
           RHOUPZERO(INDEX + 10) = FUP
           RHOUPZERO(INDEX + 11) = FUP
           RHOUPZERO(INDEX + 12) = FUP
           RHOUPZERO(INDEX + 13) = FUP

           RHODOWNZERO(INDEX + 1) = SDOWN
           RHODOWNZERO(INDEX + 2) = DDOWN
           RHODOWNZERO(INDEX + 3) = DDOWN
           RHODOWNZERO(INDEX + 4) = DDOWN
           RHODOWNZERO(INDEX + 5) = DDOWN
           RHODOWNZERO(INDEX + 6) = DDOWN
           RHODOWNZERO(INDEX + 7) = FDOWN
           RHODOWNZERO(INDEX + 8) = FDOWN
           RHODOWNZERO(INDEX + 9) = FDOWN
           RHODOWNZERO(INDEX + 10) = FDOWN
           RHODOWNZERO(INDEX + 11) = FDOWN
           RHODOWNZERO(INDEX + 12) = FDOWN
           RHODOWNZERO(INDEX + 13) = FDOWN

           INDEX = INDEX + 13

        CASE("pdf")

           ! f, then d, then p

           IF (ATOCC(ELEMPOINTER(I)) .LT. SPINMAXF) THEN

              FUP = ATOCC(ELEMPOINTER(I))/SPINMAXF
              FDOWN = ZERO
              DUP = ZERO
              DDOWN = ZERO
              PUP = ZERO
              PDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. SPINMAXF .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWO*SPINMAXF) THEN

              FUP = ONE
              FDOWN = (ATOCC(ELEMPOINTER(I)) - SPINMAXF)/SPINMAXF
              DUP = ZERO
              DDOWN = ZERO
              PUP = ZERO
              PDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. TWO*SPINMAXF .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWO*SPINMAXF + SPINMAXD) THEN

              FUP = ONE
              FDOWN = ONE
              DUP = (ATOCC(ELEMPOINTER(I)) - TWO*SPINMAXF)/SPINMAXD
              DDOWN = ZERO
              PUP = ZERO
              PDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. TWO*SPINMAXF + SPINMAXD .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWO*SPINMAXF + TWO*SPINMAXD) THEN

              FUP = ONE
              FDOWN = ONE
              DUP = ONE
              DDOWN = (ATOCC(ELEMPOINTER(I)) - &
                   (TWO*SPINMAXF + SPINMAXD))/SPINMAXD
              PUP = ZERO
              PDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. &
                TWO*SPINMAXF + TWO*SPINMAXD .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. &
                TWO*SPINMAXF + TWO*SPINMAXD + SPINMAXP) THEN

              FUP = ONE
              FDOWN = ONE
              DUP = ONE
              DDOWN = ONE
              PUP = (ATOCC(ELEMPOINTER(I)) - &
                   (TWO*SPINMAXF + TWO*SPINMAXD))/SPINMAXP
              PDOWN = ZERO

           ELSE

              FUP = ONE
              FDOWN = ONE
              DUP = ONE
              DDOWN = ONE
              PUP = ONE
              PDOWN = (ATOCC(ELEMPOINTER(I)) - &
                   (TWO*SPINMAXF + TWO*SPINMAXD + SPINMAXP))/SPINMAXP

           ENDIF

           RHOUPZERO(INDEX + 1) = PUP
           RHOUPZERO(INDEX + 2) = PUP
           RHOUPZERO(INDEX + 3) = PUP
           RHOUPZERO(INDEX + 4) = DUP
           RHOUPZERO(INDEX + 5) = DUP
           RHOUPZERO(INDEX + 6) = DUP
           RHOUPZERO(INDEX + 7) = DUP
           RHOUPZERO(INDEX + 8) = DUP
           RHOUPZERO(INDEX + 9) = FUP
           RHOUPZERO(INDEX + 10) = FUP
           RHOUPZERO(INDEX + 11) = FUP
           RHOUPZERO(INDEX + 12) = FUP
           RHOUPZERO(INDEX + 13) = FUP
           RHOUPZERO(INDEX + 14) = FUP
           RHOUPZERO(INDEX + 15) = FUP

           RHODOWNZERO(INDEX + 1) = PDOWN
           RHODOWNZERO(INDEX + 2) = PDOWN
           RHODOWNZERO(INDEX + 3) = PDOWN
           RHODOWNZERO(INDEX + 4) = DDOWN
           RHODOWNZERO(INDEX + 5) = DDOWN
           RHODOWNZERO(INDEX + 6) = DDOWN
           RHODOWNZERO(INDEX + 7) = DDOWN
           RHODOWNZERO(INDEX + 8) = DDOWN
           RHODOWNZERO(INDEX + 9) = FDOWN
           RHODOWNZERO(INDEX + 10) = FDOWN
           RHODOWNZERO(INDEX + 11) = FDOWN
           RHODOWNZERO(INDEX + 12) = FDOWN
           RHODOWNZERO(INDEX + 13) = FDOWN
           RHODOWNZERO(INDEX + 14) = FDOWN
           RHODOWNZERO(INDEX + 15) = FDOWN

           INDEX = INDEX + 15

        CASE("spdf")

           ! s, f, d, then p

           IF (ATOCC(ELEMPOINTER(I)) .LE. SPINMAXS) THEN

              SUP = ATOCC(ELEMPOINTER(I))
              SDOWN = ZERO
              FUP = ZERO
              FDOWN = ZERO
              DUP = ZERO
              DDOWN = ZERO
              PUP = ZERO
              PDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. SPINMAXS .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWO*SPINMAXS) THEN

              SUP = ONE
              SDOWN = (ATOCC(ELEMPOINTER(I)) - SPINMAXF)
              FUP = ZERO
              FDOWN = ZERO
              DUP = ZERO
              DDOWN = ZERO
              PUP = ZERO
              PDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. TWO*SPINMAXS .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWO*SPINMAXS + SPINMAXF) THEN


              SUP = ONE
              SDOWN = ONE
              FUP = (ATOCC(ELEMPOINTER(I)) - TWO*SPINMAXS)/SPINMAXF
              FDOWN = ZERO
              DUP = ZERO
              DDOWN = ZERO
              PUP = ZERO
              PDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. TWO*SPINMAXS + SPINMAXF .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. TWO*SPINMAXS + TWO*SPINMAXF) THEN

              SUP = ONE
              SDOWN = ONE
              FUP = ONE
              FDOWN = (ATOCC(ELEMPOINTER(I)) - &
                   (TWO*SPINMAXS + SPINMAXF))/SPINMAXF
              DUP = ZERO
              DDOWN = ZERO
              PUP = ZERO
              PDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. &
                TWO*SPINMAXS + TWO*SPINMAXF .AND. &
                ATOCC(ELEMPOINTER(I)) .LE. &
                TWO*SPINMAXS + TWO*SPINMAXF + SPINMAXD) THEN

              SUP = ONE
              SDOWN = ONE
              FUP = ONE
              FDOWN = ONE
              DUP = (ATOCC(ELEMPOINTER(I)) - &
                   (TWO*SPINMAXS + TWO*SPINMAXF))/SPINMAXD
              DDOWN = ZERO
              PUP = ZERO
              PDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. &
                TWO*SPINMAXS + TWO*SPINMAXF + SPINMAXD.AND. &
                ATOCC(ELEMPOINTER(I)) .LE.  &
                TWO*SPINMAXS + TWO*SPINMAXF + TWO*SPINMAXD) THEN

              SUP = ONE
              SDOWN = ONE
              FUP = ONE
              FDOWN = ONE
              DUP = ONE
              DDOWN = (ATOCC(ELEMPOINTER(I)) - &
                   (TWO*SPINMAXS + TWO*SPINMAXF + SPINMAXD))/SPINMAXD
              PUP = ZERO
              PDOWN = ZERO

           ELSEIF (ATOCC(ELEMPOINTER(I)) .GT. &
                TWO*SPINMAXS + TWO*SPINMAXF + TWO*SPINMAXD .AND. &
                ATOCC(ELEMPOINTER(I)) .LE.  &
                TWO*SPINMAXS + TWO*SPINMAXF + TWO*SPINMAXD + SPINMAXP) THEN

              SUP = ONE
              SDOWN = ONE
              FUP = ONE
              FDOWN = ONE
              DUP = ONE
              DDOWN = ONE
              PUP = (ATOCC(ELEMPOINTER(I)) - &
                   (TWO*SPINMAXS + TWO*SPINMAXF + TWO*SPINMAXD))/SPINMAXP
              PDOWN = ZERO

           ELSE

              SUP = ONE
              SDOWN = ONE
              FUP = ONE
              FDOWN = ONE
              DUP = ONE
              DDOWN = ONE
              PUP = ONE
              PDOWN = (ATOCC(ELEMPOINTER(I)) - &
                   (TWO*SPINMAXS + TWO*SPINMAXF + &
                   TWO*SPINMAXD + SPINMAXP))/SPINMAXP
           ENDIF

           RHOUPZERO(INDEX + 1) = SUP
           RHOUPZERO(INDEX + 2) = PUP
           RHOUPZERO(INDEX + 3) = PUP
           RHOUPZERO(INDEX + 4) = PUP
           RHOUPZERO(INDEX + 5) = DUP
           RHOUPZERO(INDEX + 6) = DUP
           RHOUPZERO(INDEX + 7) = DUP
           RHOUPZERO(INDEX + 8) = DUP
           RHOUPZERO(INDEX + 9) = DUP
           RHOUPZERO(INDEX + 10) = FUP
           RHOUPZERO(INDEX + 11) = FUP
           RHOUPZERO(INDEX + 12) = FUP
           RHOUPZERO(INDEX + 13) = FUP
           RHOUPZERO(INDEX + 14) = FUP
           RHOUPZERO(INDEX + 15) = FUP
           RHOUPZERO(INDEX + 16) = FUP

           RHODOWNZERO(INDEX + 1) = SDOWN
           RHODOWNZERO(INDEX + 2) = PDOWN
           RHODOWNZERO(INDEX + 3) = PDOWN
           RHODOWNZERO(INDEX + 4) = PDOWN
           RHODOWNZERO(INDEX + 5) = DDOWN
           RHODOWNZERO(INDEX + 6) = DDOWN
           RHODOWNZERO(INDEX + 7) = DDOWN
           RHODOWNZERO(INDEX + 8) = DDOWN
           RHODOWNZERO(INDEX + 9) = DDOWN
           RHODOWNZERO(INDEX + 10) = FDOWN
           RHODOWNZERO(INDEX + 11) = FDOWN
           RHODOWNZERO(INDEX + 12) = FDOWN
           RHODOWNZERO(INDEX + 13) = FDOWN
           RHODOWNZERO(INDEX + 14) = FDOWN
           RHODOWNZERO(INDEX + 15) = FDOWN
           RHODOWNZERO(INDEX + 16) = FDOWN

           INDEX = INDEX + 16

        END SELECT

     ENDDO

     IF (RESTART .EQ. 0) THEN

        !
        ! Let's initialize our self-consistent spin densities as those for
        ! free atoms
        !

        DO I = 1, HDIM

           RHOUP(I,I) = RHOUPZERO(I)
           RHODOWN(I,I) = RHODOWNZERO(I)

        ENDDO

     ELSEIF (RESTART .EQ. 1) THEN

        DO I = 1, HDIM

           RHOUP(I,I) = TMPRHOUP(I)
           RHODOWN(I,I) = TMPRHODOWN(I)

        ENDDO

        DEALLOCATE(TMPRHOUP, TMPRHODOWN)

     ENDIF

     !     print*,"a"
     !     CALL GETDELTASPIN
     !     print*, "b"
     !     CALL GETSPINE
     !     print*, "c"
     !     ESPIN_ZERO = ESPIN

  ENDIF

  RETURN

END SUBROUTINE RHOZERO
