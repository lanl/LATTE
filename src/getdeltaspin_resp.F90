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

SUBROUTINE GETDELTASPIN_RESP

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE SPINARRAY
  USE KSPACEARRAY
  USE NONOARRAY
  USE MYPRECISION

  IMPLICIT NONE

  REAL(LATTEPREC), ALLOCATABLE :: SPINDEN(:,:)
  INTEGER :: I, J, K, INDEX, DINDEX
  IF (EXISTERROR) RETURN

  INDEX = 0
  DINDEX = 0

  ALLOCATE(SPINDEN(HDIM,NSPIN))
  !  deltaspin = number of spin ups - number of spin downs

  DELTA_QS = ZERO
  IF (BASISTYPE .EQ. "ORTHO") THEN

     SPINDEN = ZERO
     DO I = 1, HDIM
        SPINDEN(I,1) = RHOUP(I,I)
        SPINDEN(I,2) = RHODOWN(I,I)
     ENDDO

  ELSEIF (BASISTYPE .EQ. "NONORTHO") THEN

     !
     ! Mulliken spin densities in non-orthogonal basis:

     ! m = n_up - n_down

     ! n_sigma = partial_trace(rho_sigma S)

     IF (KON .EQ. 0) THEN

        SPINLIST = ZERO
        SPINDEN = ZERO

        DO I = 1, HDIM
           DO J = 1, HDIM

              SPINLIST(I) = SPINLIST(I) + (RHOUP(J,I) - RHODOWN(J,I))*SMAT(J,I)
              SPINDEN(I,1) = SPINDEN(I,1) + RHOUP(J,I) * SMAT(J,I)
              SPINDEN(I,2) = SPINDEN(I,2) + RHODOWN(J,I) * SMAT(J,I)

           ENDDO
        ENDDO

     ELSE ! If we're doing magnetic TB in k-space
        
        ZSPINLIST = (ZERO, ZERO)

        DO K = 1, NKTOT
           DO I = 1, HDIM
              DO J = 1, HDIM

                 ZSPINLIST(I) = ZSPINLIST(I) + ((KRHOUP(J,I,K) - KRHODOWN(J,I,K))*SK(I,J,K) + &
                      (KRHOUP(I,J,K) - KRHODOWN(I,J,K))*SK(J,I,K))/TWO
                 
              ENDDO
           ENDDO
        ENDDO

        ZSPINLIST = ZSPINLIST/REAL(NKTOT)

        SPINLIST = REAL(ZSPINLIST)

     ENDIF

     DO I = 1, NATS

        SELECT CASE(BASIS(ELEMPOINTER(I)))

        CASE("s")

           DINDEX = DINDEX + 1
           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J)
           ENDDO

           INDEX = INDEX + 1

        CASE("p")

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1, J) + &
                SPINDEN(INDEX + 2, J) + SPINDEN(INDEX + 3, J)
           ENDDO

           INDEX = INDEX + 3

        CASE("d")

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J) + &
                SPINDEN(INDEX + 2,J) + SPINDEN(INDEX + 3,J) + &
                SPINDEN(INDEX + 4,J) + SPINDEN(INDEX + 5,J)
           ENDDO

           INDEX = INDEX + 5

        CASE("f")

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J) + &
                SPINDEN(INDEX + 2,J) + SPINDEN(INDEX + 3,J) + &
                SPINDEN(INDEX + 4,J) + SPINDEN(INDEX + 5,J)  + &
                SPINDEN(INDEX + 6,J) + SPINDEN(INDEX + 7,J)
           ENDDO

           INDEX = INDEX + 7

        CASE("sp")

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J)
           ENDDO

           INDEX = INDEX + 1

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J) + &
                SPINDEN(INDEX + 2,J) + SPINDEN(INDEX + 3,J)
           ENDDO

           INDEX = INDEX + 3

        CASE("sd")

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J)
           ENDDO

           INDEX = INDEX + 1

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J) + &
                SPINDEN(INDEX + 2,J) + SPINDEN(INDEX + 3,J) + &
                SPINDEN(INDEX + 4,J) + SPINDEN(INDEX + 5,J)
           ENDDO

           INDEX = INDEX + 5

        CASE("sf")

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J)
           ENDDO

           INDEX = INDEX + 1

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J) + &
                SPINDEN(INDEX + 2,J) + SPINDEN(INDEX + 3,J) + &
                SPINDEN(INDEX + 4,J) + SPINDEN(INDEX + 5,J)  + &
                SPINDEN(INDEX + 6,J) + SPINDEN(INDEX + 7,J)
           ENDDO

           INDEX = INDEX + 7           


        CASE("pd")

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J) + &
                SPINDEN(INDEX + 2,J) + SPINDEN(INDEX + 3,J)
           ENDDO

           INDEX = INDEX + 3           

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J) + &
                SPINDEN(INDEX + 2,J) + SPINDEN(INDEX + 3,J) + &
                SPINDEN(INDEX + 4,J) + SPINDEN(INDEX + 5,J)
           ENDDO

           INDEX = INDEX + 5

        CASE("pf")

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J) + &
                SPINDEN(INDEX + 2,J) + SPINDEN(INDEX + 3,J)
           ENDDO

           INDEX = INDEX + 3 

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J) + &
                SPINDEN(INDEX + 2,J) + SPINDEN(INDEX + 3,J) + &
                SPINDEN(INDEX + 4,J) + SPINDEN(INDEX + 5,J)  + &
                SPINDEN(INDEX + 6,J) + SPINDEN(INDEX + 7,J)
           ENDDO

           INDEX = INDEX + 7                  

        CASE("df")

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J) + &
                SPINDEN(INDEX + 2,J) + SPINDEN(INDEX + 3,J) + &
                SPINDEN(INDEX + 4,J) + SPINDEN(INDEX + 5,J)
           ENDDO

           INDEX = INDEX + 5

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J) + &
                SPINDEN(INDEX + 2,J) + SPINDEN(INDEX + 3,J) + &
                SPINDEN(INDEX + 4,J) + SPINDEN(INDEX + 5,J)  + &
                SPINDEN(INDEX + 6,J) + SPINDEN(INDEX + 7,J)
           ENDDO

           INDEX = INDEX + 7     


        CASE("spd")

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J)
           ENDDO

           INDEX = INDEX + 1

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J) + &
                SPINDEN(INDEX + 2,J) + SPINDEN(INDEX + 3,J)
           ENDDO

           INDEX = INDEX + 3

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J) + &
                SPINDEN(INDEX + 2,J) + SPINDEN(INDEX + 3,J) + &
                SPINDEN(INDEX + 4,J) + SPINDEN(INDEX + 5,J)
           ENDDO

           INDEX = INDEX + 5

        CASE("spf")

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J)
           ENDDO

           INDEX = INDEX + 1

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J) + &
                SPINDEN(INDEX + 2,J) + SPINDEN(INDEX + 3,J)
           ENDDO

           INDEX = INDEX + 3

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J) + &
                SPINDEN(INDEX + 2,J) + SPINDEN(INDEX + 3,J) + &
                SPINDEN(INDEX + 4,J) + SPINDEN(INDEX + 5,J)  + &
                SPINDEN(INDEX + 6,J) + SPINDEN(INDEX + 7,J)
           ENDDO

           INDEX = INDEX + 7                

        CASE("sdf")

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J)
           ENDDO

           INDEX = INDEX + 1

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J) + &
                SPINDEN(INDEX + 2,J) + SPINDEN(INDEX + 3,J) + &
                SPINDEN(INDEX + 4,J) + SPINDEN(INDEX + 5,J)
           ENDDO

           INDEX = INDEX + 5

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J) + &
                SPINDEN(INDEX + 2,J) + SPINDEN(INDEX + 3,J) + &
                SPINDEN(INDEX + 4,J) + SPINDEN(INDEX + 5,J)  + &
                SPINDEN(INDEX + 6,J) + SPINDEN(INDEX + 7,J)
           ENDDO

           INDEX = INDEX + 7

        CASE("pdf")

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J) + &
                SPINDEN(INDEX + 2,J) + SPINDEN(INDEX + 3,J)
           ENDDO

           INDEX = INDEX + 3 

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J) + &
                SPINDEN(INDEX + 2,J) + SPINDEN(INDEX + 3,J) + &
                SPINDEN(INDEX + 4,J) + SPINDEN(INDEX + 5,J)
           ENDDO

           INDEX = INDEX + 5

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J) + &
                SPINDEN(INDEX + 2,J) + SPINDEN(INDEX + 3,J) + &
                SPINDEN(INDEX + 4,J) + SPINDEN(INDEX + 5,J)  + &
                SPINDEN(INDEX + 6,J) + SPINDEN(INDEX + 7,J)
           ENDDO

           INDEX = INDEX + 7    

        CASE("spdf") 

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J)
           ENDDO

           INDEX = INDEX + 1

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J) + &
                SPINDEN(INDEX + 2,J) + SPINDEN(INDEX + 3,J)
           ENDDO

           INDEX = INDEX + 3 

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J) + &
                SPINDEN(INDEX + 2,J) + SPINDEN(INDEX + 3,J) + &
                SPINDEN(INDEX + 4,J) + SPINDEN(INDEX + 5,J)
           ENDDO

           INDEX = INDEX + 5

           DINDEX = DINDEX + 1

           DO J = 1, 2
             DELTA_QS(DINDEX,J) = SPINDEN(INDEX + 1,J) + &
                SPINDEN(INDEX + 2,J) + SPINDEN(INDEX + 3,J) + &
                SPINDEN(INDEX + 4,J) + SPINDEN(INDEX + 5,J)  + &
                SPINDEN(INDEX + 6,J) + SPINDEN(INDEX + 7,J)
           ENDDO

           INDEX = INDEX + 7                  
        END SELECT

     ENDDO

  ENDIF


  DEALLOCATE(SPINDEN)

  RETURN

END SUBROUTINE GETDELTASPIN_RESP


SUBROUTINE REDUCE_DELTASPIN(NATS, DELTADIM, DELTASPIN, DELTAQ,OPT)
USE MYPRECISION
USE SETUPARRAY, ONLY: BASIS, ELEMPOINTER, ATOCC
IMPLICIT NONE

INTEGER,         INTENT(IN)  :: NATS, DELTADIM
REAL(LATTEPREC), INTENT(IN)  :: DELTASPIN(DELTADIM)
REAL(LATTEPREC), INTENT(OUT) :: DELTAQ(NATS)
INTEGER,         INTENT(IN)  :: OPT !WEHTER - ATOCC
!
INTEGER :: I, J, INDEX2, NUMORB2

DELTAQ = 0.0_LATTEPREC
INDEX2 = 0
DO I = 1, NATS

   SELECT CASE(BASIS(ELEMPOINTER(I)))

   CASE("s")

      NUMORB2 = 1

   CASE("p")

      NUMORB2 = 1

   CASE("d")

      NUMORB2 = 1

   CASE("f")

      NUMORB2 = 1

   CASE("sp")

      NUMORB2 = 2

   CASE("sd")

      NUMORB2 = 2

   CASE("sf")

      NUMORB2 = 2

   CASE("pd")

      NUMORB2 = 2

   CASE("pf")

      NUMORB2 = 2

   CASE("df")

      NUMORB2 = 2

   CASE("spd")

      NUMORB2 = 3

   CASE("spf")

      NUMORB2 = 3

   CASE("sdf")

      NUMORB2 = 3

   CASE("pdf")

      NUMORB2 = 3

   CASE("spdf") 

      NUMORB2 = 4

   END SELECT

   DO J = 1, NUMORB2
     INDEX2 = INDEX2 + 1
     DELTAQ(I) = DELTAQ(I) + DELTASPIN(INDEX2)
   ENDDO

   IF (OPT==2) DELTAQ(I) = DELTAQ(I) - ATOCC(ELEMPOINTER(I))
ENDDO

END SUBROUTINE REDUCE_DELTASPIN
