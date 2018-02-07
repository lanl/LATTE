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

SUBROUTINE GETDELTAQ

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE SPINARRAY
  USE NONOARRAY
  USE COULOMBARRAY
  USE MYPRECISION
  USE KSPACEARRAY

  IMPLICIT NONE

  INTEGER :: I, J, K, INDEX, NUMORB
  COMPLEX(LATTEPREC) :: QTMP
  IF (EXISTERROR) RETURN
  !  COMPLEX(LATTEPREC), ALLOCATABLE :: TMPQ(:)

  !  ALLOCATE(TMPQ(HDIM))

  INDEX = 0

  QLIST = ZERO
  MYCHARGE = ZERO

  IF (KON .EQ. 0 ) THEN

     IF (BASISTYPE .EQ. "ORTHO") THEN

        IF (SPINON .EQ. 0) THEN

           DO I = 1, HDIM
              QLIST(I) = BO(I,I)
           ENDDO

        ELSE

           DO I = 1, HDIM
              QLIST(I) = RHOUP(I,I) + RHODOWN(I,I)
           ENDDO

        ENDIF

     ELSEIF (BASISTYPE .EQ. "NONORTHO") THEN

        IF (SPINON .EQ. 0) THEN

           DO I = 1, HDIM
              DO J = 1, HDIM
                 QLIST(I) = QLIST(I) + BO(J,I)*SMAT(J,I)
              ENDDO
           ENDDO

        ELSE

           DO I = 1, HDIM
              DO J = 1, HDIM

                 QLIST(I) = QLIST(I) + (RHOUP(J,I)+RHODOWN(J,I))*SMAT(J,I)

              ENDDO
           ENDDO

        ENDIF

     ENDIF

  ELSE

     ! Loop over k-points

     IF (BASISTYPE .EQ. "ORTHO") THEN

        DO K = 1, NKTOT
           DO I = 1, HDIM
              QLIST(I) = QLIST(I) + REAL(KBO(I,I,K))
           ENDDO
        ENDDO

        QLIST = QLIST/REAL(NKTOT)

     ELSE

        !        TMPQ = (ZERO, ZERO)
        DO K = 1, NKTOT
           DO I = 1, HDIM
              DO J = 1, HDIM
                 !                 TMPQ(I) = TMPQ(I) + KBO(J,I,K)*CONJG(SK(I,J,K))
                 QTMP = HALF*((KBO(J,I,K))*SK(I,J,K) + &
                      (KBO(I,J,K))*SK(J,I,K))
                 !                 QTMP = HALF*((KBO(J,I,K))*SK(I,J,K) + &
                 !                      (KBO(J,I,K))*SK(I,J,K))
                 QLIST(I) = QLIST(I) + REAL(QTMP)
                 !                 QLIST(I) = QLIST(I) + REAL(KBO(J,I,K)*CONJG(SK(J,I,K)))
              ENDDO
           ENDDO
        ENDDO
        !        TMPQ = TMPQ/REAL(NKTOT)

        QLIST = QLIST/REAL(NKTOT)
        !         DO I =1, HDIM
        !           PRINT*, I, TMPQ(I), QLIST(I)
        !        ENDDO


        !        PRINT*, "K-space not yet implemented with non-orthogonal basis"
        !        STOP

     ENDIF

  ENDIF


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

     !     MYCHARGE = ZERO
     DO J = 1, NUMORB

        INDEX = INDEX + 1
        MYCHARGE(I) = MYCHARGE(I) + QLIST(INDEX)

     ENDDO

     DELTAQ(I) = MYCHARGE(I) - ATOCC(ELEMPOINTER(I))

  ENDDO
  !  deallocate(tmpq)

  RETURN

END SUBROUTINE GETDELTAQ
