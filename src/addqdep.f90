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

SUBROUTINE ADDQDEP

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE COULOMBARRAY
  USE NONOARRAY
  USE KSPACEARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, K
  INTEGER :: INDEX, NUMORB
  REAL(LATTEPREC) :: ES, EP, ED, EF
  REAL(LATTEPREC) :: HMOD
  COMPLEX(LATTEPREC) :: ZHMOD
  IF (EXISTERROR) RETURN

  INDEX = 0

  IF (KON .EQ. 0) THEN

     IF ( BASISTYPE .EQ. "ORTHO") THEN

        DO I = 1, NATS

           HMOD = HUBBARDU(ELEMPOINTER(I))*DELTAQ(I) + COULOMBV(I)

           SELECT CASE(BASIS(ELEMPOINTER(I)))

           CASE("s")

              H(INDEX + 1, INDEX + 1) = HDIAG(INDEX + 1) + HMOD
              INDEX = INDEX + 1

           CASE("p")

              H(INDEX + 1, INDEX + 1) = HDIAG(INDEX + 1) + HMOD
              H(INDEX + 2, INDEX + 2) = HDIAG(INDEX + 2) + HMOD
              H(INDEX + 3, INDEX + 3) = HDIAG(INDEX + 3) + HMOD
              INDEX = INDEX + 3

           CASE("d")

              H(INDEX + 1, INDEX + 1) = HDIAG(INDEX + 1) + HMOD
              H(INDEX + 2, INDEX + 2) = HDIAG(INDEX + 2) + HMOD
              H(INDEX + 3, INDEX + 3) = HDIAG(INDEX + 3) + HMOD
              H(INDEX + 4, INDEX + 4) = HDIAG(INDEX + 4) + HMOD
              H(INDEX + 5, INDEX + 5) = HDIAG(INDEX + 5) + HMOD
              INDEX = INDEX + 5

           CASE("f")

              H(INDEX + 1, INDEX + 1) = HDIAG(INDEX + 1) + HMOD
              H(INDEX + 2, INDEX + 2) = HDIAG(INDEX + 2) + HMOD
              H(INDEX + 3, INDEX + 3) = HDIAG(INDEX + 3) + HMOD
              H(INDEX + 4, INDEX + 4) = HDIAG(INDEX + 4) + HMOD
              H(INDEX + 5, INDEX + 5) = HDIAG(INDEX + 5) + HMOD
              H(INDEX + 6, INDEX + 6) = HDIAG(INDEX + 6) + HMOD
              H(INDEX + 7, INDEX + 7) = HDIAG(INDEX + 7) + HMOD
              INDEX = INDEX + 7

           CASE("sp")

              H(INDEX + 1, INDEX + 1) = HDIAG(INDEX + 1) + HMOD
              H(INDEX + 2, INDEX + 2) = HDIAG(INDEX + 2) + HMOD
              H(INDEX + 3, INDEX + 3) = HDIAG(INDEX + 3) + HMOD
              H(INDEX + 4, INDEX + 4) = HDIAG(INDEX + 4) + HMOD
              INDEX = INDEX + 4

           CASE("sd")

              H(INDEX + 1, INDEX + 1) = HDIAG(INDEX + 1) + HMOD
              H(INDEX + 2, INDEX + 2) = HDIAG(INDEX + 2) + HMOD
              H(INDEX + 3, INDEX + 3) = HDIAG(INDEX + 3) + HMOD
              H(INDEX + 4, INDEX + 4) = HDIAG(INDEX + 4) + HMOD
              H(INDEX + 5, INDEX + 5) = HDIAG(INDEX + 5) + HMOD
              H(INDEX + 6, INDEX + 6) = HDIAG(INDEX + 6) + HMOD
              INDEX = INDEX + 6

           CASE("sf")

              H(INDEX + 1, INDEX + 1) = HDIAG(INDEX + 1) + HMOD
              H(INDEX + 2, INDEX + 2) = HDIAG(INDEX + 2) + HMOD
              H(INDEX + 3, INDEX + 3) = HDIAG(INDEX + 3) + HMOD
              H(INDEX + 4, INDEX + 4) = HDIAG(INDEX + 4) + HMOD
              H(INDEX + 5, INDEX + 5) = HDIAG(INDEX + 5) + HMOD
              H(INDEX + 6, INDEX + 6) = HDIAG(INDEX + 6) + HMOD
              H(INDEX + 7, INDEX + 7) = HDIAG(INDEX + 7) + HMOD
              H(INDEX + 8, INDEX + 8) = HDIAG(INDEX + 8) + HMOD
              INDEX = INDEX + 8

           CASE("pd")

              H(INDEX + 1, INDEX + 1) = HDIAG(INDEX + 1) + HMOD
              H(INDEX + 2, INDEX + 2) = HDIAG(INDEX + 2) + HMOD
              H(INDEX + 3, INDEX + 3) = HDIAG(INDEX + 3) + HMOD
              H(INDEX + 4, INDEX + 4) = HDIAG(INDEX + 4) + HMOD
              H(INDEX + 5, INDEX + 5) = HDIAG(INDEX + 5) + HMOD
              H(INDEX + 6, INDEX + 6) = HDIAG(INDEX + 6) + HMOD
              H(INDEX + 7, INDEX + 7) = HDIAG(INDEX + 7) + HMOD
              H(INDEX + 8, INDEX + 8) = HDIAG(INDEX + 8) + HMOD
              INDEX = INDEX + 8

           CASE("pf")

              H(INDEX + 1, INDEX + 1) = HDIAG(INDEX + 1) + HMOD
              H(INDEX + 2, INDEX + 2) = HDIAG(INDEX + 2) + HMOD
              H(INDEX + 3, INDEX + 3) = HDIAG(INDEX + 3) + HMOD
              H(INDEX + 4, INDEX + 4) = HDIAG(INDEX + 4) + HMOD
              H(INDEX + 5, INDEX + 5) = HDIAG(INDEX + 5) + HMOD
              H(INDEX + 6, INDEX + 6) = HDIAG(INDEX + 6) + HMOD
              H(INDEX + 7, INDEX + 7) = HDIAG(INDEX + 7) + HMOD
              H(INDEX + 8, INDEX + 8) = HDIAG(INDEX + 8) + HMOD
              H(INDEX + 9, INDEX + 9) = HDIAG(INDEX + 9) + HMOD
              H(INDEX + 10, INDEX + 10) = HDIAG(INDEX + 10) + HMOD
              INDEX = INDEX + 10

           CASE("df")

              H(INDEX + 1, INDEX + 1) = HDIAG(INDEX + 1) + HMOD
              H(INDEX + 2, INDEX + 2) = HDIAG(INDEX + 2) + HMOD
              H(INDEX + 3, INDEX + 3) = HDIAG(INDEX + 3) + HMOD
              H(INDEX + 4, INDEX + 4) = HDIAG(INDEX + 4) + HMOD
              H(INDEX + 5, INDEX + 5) = HDIAG(INDEX + 5) + HMOD
              H(INDEX + 6, INDEX + 6) = HDIAG(INDEX + 6) + HMOD
              H(INDEX + 7, INDEX + 7) = HDIAG(INDEX + 7) + HMOD
              H(INDEX + 8, INDEX + 8) = HDIAG(INDEX + 8) + HMOD
              H(INDEX + 9, INDEX + 9) = HDIAG(INDEX + 9) + HMOD
              H(INDEX + 10, INDEX + 10) = HDIAG(INDEX + 10) + HMOD
              H(INDEX + 11, INDEX + 11) = HDIAG(INDEX + 11) + HMOD
              H(INDEX + 12, INDEX + 12) = HDIAG(INDEX + 12) + HMOD
              INDEX = INDEX + 12

           CASE("spd")

              H(INDEX + 1, INDEX + 1) = HDIAG(INDEX + 1) + HMOD
              H(INDEX + 2, INDEX + 2) = HDIAG(INDEX + 2) + HMOD
              H(INDEX + 3, INDEX + 3) = HDIAG(INDEX + 3) + HMOD
              H(INDEX + 4, INDEX + 4) = HDIAG(INDEX + 4) + HMOD
              H(INDEX + 5, INDEX + 5) = HDIAG(INDEX + 5) + HMOD
              H(INDEX + 6, INDEX + 6) = HDIAG(INDEX + 6) + HMOD
              H(INDEX + 7, INDEX + 7) = HDIAG(INDEX + 7) + HMOD
              H(INDEX + 8, INDEX + 8) = HDIAG(INDEX + 8) + HMOD
              H(INDEX + 9, INDEX + 9) = HDIAG(INDEX + 9) + HMOD
              INDEX = INDEX + 9

           CASE("spf")

              H(INDEX + 1, INDEX + 1) = HDIAG(INDEX + 1) + HMOD
              H(INDEX + 2, INDEX + 2) = HDIAG(INDEX + 2) + HMOD
              H(INDEX + 3, INDEX + 3) = HDIAG(INDEX + 3) + HMOD
              H(INDEX + 4, INDEX + 4) = HDIAG(INDEX + 4) + HMOD
              H(INDEX + 5, INDEX + 5) = HDIAG(INDEX + 5) + HMOD
              H(INDEX + 6, INDEX + 6) = HDIAG(INDEX + 6) + HMOD
              H(INDEX + 7, INDEX + 7) = HDIAG(INDEX + 7) + HMOD
              H(INDEX + 8, INDEX + 8) = HDIAG(INDEX + 8) + HMOD
              H(INDEX + 9, INDEX + 9) = HDIAG(INDEX + 9) + HMOD
              H(INDEX + 10, INDEX + 10) = HDIAG(INDEX + 10) + HMOD
              H(INDEX + 11, INDEX + 11) = HDIAG(INDEX + 11) + HMOD
              INDEX = INDEX + 11

           CASE("sdf")

              H(INDEX + 1, INDEX + 1) = HDIAG(INDEX + 1) + HMOD
              H(INDEX + 2, INDEX + 2) = HDIAG(INDEX + 2) + HMOD
              H(INDEX + 3, INDEX + 3) = HDIAG(INDEX + 3) + HMOD
              H(INDEX + 4, INDEX + 4) = HDIAG(INDEX + 4) + HMOD
              H(INDEX + 5, INDEX + 5) = HDIAG(INDEX + 5) + HMOD
              H(INDEX + 6, INDEX + 6) = HDIAG(INDEX + 6) + HMOD
              H(INDEX + 7, INDEX + 7) = HDIAG(INDEX + 7) + HMOD
              H(INDEX + 8, INDEX + 8) = HDIAG(INDEX + 8) + HMOD
              H(INDEX + 9, INDEX + 9) = HDIAG(INDEX + 9) + HMOD
              H(INDEX + 10, INDEX + 10) = HDIAG(INDEX + 10) + HMOD
              H(INDEX + 11, INDEX + 11) = HDIAG(INDEX + 11) + HMOD
              H(INDEX + 12, INDEX + 12) = HDIAG(INDEX + 12) + HMOD
              H(INDEX + 13, INDEX + 13) = HDIAG(INDEX + 13) + HMOD
              INDEX = INDEX + 13

           CASE("pdf")

              H(INDEX + 1, INDEX + 1) = HDIAG(INDEX + 1) + HMOD
              H(INDEX + 2, INDEX + 2) = HDIAG(INDEX + 2) + HMOD
              H(INDEX + 3, INDEX + 3) = HDIAG(INDEX + 3) + HMOD
              H(INDEX + 4, INDEX + 4) = HDIAG(INDEX + 4) + HMOD
              H(INDEX + 5, INDEX + 5) = HDIAG(INDEX + 5) + HMOD
              H(INDEX + 6, INDEX + 6) = HDIAG(INDEX + 6) + HMOD
              H(INDEX + 7, INDEX + 7) = HDIAG(INDEX + 7) + HMOD
              H(INDEX + 8, INDEX + 8) = HDIAG(INDEX + 8) + HMOD
              H(INDEX + 9, INDEX + 9) = HDIAG(INDEX + 9) + HMOD
              H(INDEX + 10, INDEX + 10) = HDIAG(INDEX + 10) + HMOD
              H(INDEX + 11, INDEX + 11) = HDIAG(INDEX + 11) + HMOD
              H(INDEX + 12, INDEX + 12) = HDIAG(INDEX + 12) + HMOD
              H(INDEX + 13, INDEX + 13) = HDIAG(INDEX + 13) + HMOD
              H(INDEX + 14, INDEX + 14) = HDIAG(INDEX + 14) + HMOD
              H(INDEX + 15, INDEX + 15) = HDIAG(INDEX + 15) + HMOD
              INDEX = INDEX + 15

           CASE("spdf") 

              H(INDEX + 1, INDEX + 1) = HDIAG(INDEX + 1) + HMOD
              H(INDEX + 2, INDEX + 2) = HDIAG(INDEX + 2) + HMOD
              H(INDEX + 3, INDEX + 3) = HDIAG(INDEX + 3) + HMOD
              H(INDEX + 4, INDEX + 4) = HDIAG(INDEX + 4) + HMOD
              H(INDEX + 5, INDEX + 5) = HDIAG(INDEX + 5) + HMOD
              H(INDEX + 6, INDEX + 6) = HDIAG(INDEX + 6) + HMOD
              H(INDEX + 7, INDEX + 7) = HDIAG(INDEX + 7) + HMOD
              H(INDEX + 8, INDEX + 8) = HDIAG(INDEX + 8) + HMOD
              H(INDEX + 9, INDEX + 9) = HDIAG(INDEX + 9) + HMOD
              H(INDEX + 10, INDEX + 10) = HDIAG(INDEX + 10) + HMOD
              H(INDEX + 11, INDEX + 11) = HDIAG(INDEX + 11) + HMOD
              H(INDEX + 12, INDEX + 12) = HDIAG(INDEX + 12) + HMOD
              H(INDEX + 13, INDEX + 13) = HDIAG(INDEX + 13) + HMOD
              H(INDEX + 14, INDEX + 14) = HDIAG(INDEX + 14) + HMOD
              H(INDEX + 15, INDEX + 15) = HDIAG(INDEX + 15) + HMOD
              H(INDEX + 16, INDEX + 16) = HDIAG(INDEX + 16) + HMOD
              INDEX = INDEX + 16

           END SELECT

        ENDDO

     ELSEIF ( BASISTYPE .EQ. "NONORTHO" ) THEN

        !
        ! When we have a non-orthogonal basis, the electrostatic
        ! potential enters as SH_1 
        !

        DO I = 1, NATS

           HMOD = HUBBARDU(ELEMPOINTER(I))*DELTAQ(I) + COULOMBV(I)

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

           DO J = 1, NUMORB

              INDEX = INDEX + 1
              HJJ(INDEX) = HMOD

           ENDDO

        ENDDO

        ! H = H_0 + S*H_1

        DO I = 1, HDIM
           DO J = 1, HDIM

              H(J,I) = H0(J,I) + SMAT(J,I)*(HJJ(I) + HJJ(J))/TWO

           ENDDO
        ENDDO

     ENDIF

  ELSE
     !  IF (KON .EQ. 1) THEN 
     ! k-space - we have to add the potential to all NKTOT Hamiltonians

     ! Orthogonal basis only at the moment

     IF ( BASISTYPE .EQ. "ORTHO") THEN

        DO K = 1, NKTOT

           INDEX = 0

           DO I = 1, NATS

              ZHMOD = CMPLX(HUBBARDU(ELEMPOINTER(I))*DELTAQ(I) &
                   + COULOMBV(I))

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

              DO J = 1, NUMORB

                 INDEX = INDEX + 1
                 HK(INDEX, INDEX, K) = HKDIAG(INDEX, K) + ZHMOD

              ENDDO

           ENDDO

        ENDDO

     ELSEIF (BASISTYPE .EQ. "NONORTHO") THEN ! k-space and non-orthogonal basis

        DO I = 1, NATS

           HMOD = HUBBARDU(ELEMPOINTER(I))*DELTAQ(I) + COULOMBV(I)

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

           DO J = 1, NUMORB

              INDEX = INDEX + 1
              ZHJJ(INDEX) = CMPLX(HMOD)

           ENDDO

        ENDDO

        ! H = H_0 + S*H_1                                                             

        DO K = 1, NKTOT
           DO I = 1, HDIM
              DO J = 1, HDIM

                 HK(J,I,K) = HK0(J,I,K) + SK(J,I,K)*(ZHJJ(I) + ZHJJ(J))/TWO

              ENDDO
           ENDDO
        ENDDO

     ENDIF



  ENDIF

  IF (DEBUGON .EQ. 1) THEN

     OPEN(UNIT=31, STATUS="UNKNOWN", FILE="myH.dat")

     PRINT*, "Caution - the H0+H1 matrix is being written to file"

     IF (KON .EQ. 0) THEN

        DO I = 1, HDIM
           WRITE(31,10) (H(I,J), J = 1, HDIM)
        ENDDO

     ELSE

        DO K = 1, NKTOT
           WRITE(31,*) K
           DO I = 1, HDIM
              WRITE(31,12) (HK(I,J,K), J = 1, HDIM)
           ENDDO
        ENDDO

     ENDIF

     CLOSE(31)

  ENDIF

10 FORMAT(100G18.9)  
11 FORMAT(I5,100G18.9)  
12 FORMAT(100F8.3)  
  RETURN

END SUBROUTINE ADDQDEP


