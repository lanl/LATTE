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

SUBROUTINE GETRESPF

  ! Use our already computed eigenvalues and eigenvectors to compute
  ! the response function while enforcing LCN

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE DIAGARRAY
  USE KSPACEARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, K, KK, INDEX, NUMORB
  REAL(LATTEPREC), ALLOCATABLE :: RESPF(:)
  IF (EXISTERROR) RETURN

  ALLOCATE(RESPF(HDIM))

  RESPF = ZERO

  IF (KON .EQ. 0) THEN

     DO I = 1, HDIM

        IF (EVALS(I) .LE. CHEMPOT) THEN

           DO J = I, HDIM

              IF (EVALS(J) .GT. CHEMPOT .AND. ABS(EVALS(I)-EVALS(J)) .GT. 1.0D-4) THEN

                 DO K = 1, HDIM

                    RESPF(K) = RESPF(K) + &
                         TWO*EVECS(K,I)*EVECS(K,I)*EVECS(K,J)*EVECS(K,J)/ &
                         (EVALS(I) - EVALS(J))

                 ENDDO

              ENDIF
           ENDDO

        ENDIF

     ENDDO

  ELSE ! K-space integration now

     DO KK = 1, NKTOT

        DO I = 1, HDIM

           IF (KEVALS(I,KK) .LE. CHEMPOT) THEN

              DO J = I, HDIM

                 IF (KEVALS(J,KK) .GT. CHEMPOT) THEN

                    DO K = 1, HDIM

                       RESPF(K) = RESPF(K) - &
                            (FOUR/(KEVALS(I,KK) - KEVALS(J,KK)))*&
                            REAL(CONJG(KEVECS(K,I,KK))*KEVECS(K,I,KK) &
                            *CONJG(KEVECS(K,J,KK))*KEVECS(K,J,KK))


                    ENDDO

                 ENDIF
              ENDDO

           ENDIF

        ENDDO

     ENDDO

     RESPF = RESPF/REAL(NKTOT)

  ENDIF

  RESPCHI = ZERO
  INDEX = 0

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

     DO J = 1, NUMORB
        INDEX = INDEX + 1
        RESPCHI(I) = RESPCHI(I) + RESPF(INDEX)
     ENDDO

     RESPCHI(I) = RESPCHI(I)/REAL(NUMORB)

  ENDDO

  !  DO I = 1, NATS
  !     PRINT*, RESPCHI(I)
  !  ENDDO


  DEALLOCATE(RESPF)

  RETURN

END SUBROUTINE GETRESPF

