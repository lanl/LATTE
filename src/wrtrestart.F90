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

SUBROUTINE WRTRESTART(ITER)

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE KSPACEARRAY
  USE MDARRAY
  USE SPINARRAY
#ifdef MPI_ON
  USE MPI
#endif

  IMPLICIT NONE

  INTEGER :: I, K, ITER
  INTEGER :: MYID, IERR
  COMPLEX(LATTEPREC) :: KSUM
  CHARACTER(LEN=100) :: FLNM
  IF (EXISTERROR) RETURN

  IF( VERBOSE < 0 ) RETURN

  IF (MDON .EQ. 1) THEN

     IF (PARREP .EQ. 0) THEN

        IF (ITER .LT. 10) THEN
           WRITE(FLNM,'("Restarts/restartMD.", I1,".dat")') ITER
        ELSEIF (ITER .GE. 10 .AND. ITER .LT. 100) THEN
           WRITE(FLNM,'("Restarts/restartMD.", I2,".dat")') ITER
        ELSEIF (ITER .GE. 100 .AND. ITER .LT. 1000) THEN
           WRITE(FLNM,'("Restarts/restartMD.", I3,".dat")') ITER
        ELSEIF (ITER .GE. 1000 .AND. ITER .LT. 10000) THEN
           WRITE(FLNM,'("Restarts/restartMD.", I4,".dat")') ITER
        ELSEIF (ITER .GE. 10000 .AND. ITER .LT. 100000) THEN
           WRITE(FLNM,'("Restarts/restartMD.", I5,".dat")') ITER
        ELSEIF (ITER .GE. 100000 .AND. ITER .LT. 1000000) THEN
           WRITE(FLNM,'("Restarts/restartMD.", I6,".dat")') ITER
        ELSEIF (ITER .GE. 1000000 .AND. ITER .LT. 10000000) THEN
           WRITE(FLNM,'("Restarts/restartMD.", I7,".dat")') ITER
        ENDIF

        OPEN (UNIT = 21, STATUS="UNKNOWN", FILE=FLNM)
        OPEN (UNIT = 19, STATUS="UNKNOWN", &
             FILE="Restarts/restartMD.last.dat")

     ELSE

#ifdef MPI_ON
        CALL MPI_COMM_RANK( MPI_COMM_WORLD, MYID,  IERR )
#endif

        IF (MYID .LT. 10) THEN

           IF (ITER .LT. 10) THEN
              WRITE(FLNM,'(I1,"/Restarts/restartMD.", I1,".dat")') MYID,ITER
           ELSEIF (ITER .GE. 10 .AND. ITER .LT. 100) THEN
              WRITE(FLNM,'(I1,"/Restarts/restartMD.", I2,".dat")') MYID,ITER
           ELSEIF (ITER .GE. 100 .AND. ITER .LT. 1000) THEN
              WRITE(FLNM,'(I1,"/Restarts/restartMD.", I3,".dat")') MYID,ITER
           ELSEIF (ITER .GE. 1000 .AND. ITER .LT. 10000) THEN
              WRITE(FLNM,'(I1,"/Restarts/restartMD.", I4,".dat")') MYID,ITER
           ELSEIF (ITER .GE. 10000 .AND. ITER .LT. 100000) THEN
              WRITE(FLNM,'(I1,"/Restarts/restartMD.", I5,".dat")') MYID,ITER
           ELSEIF (ITER .GE. 100000 .AND. ITER .LT. 1000000) THEN
              WRITE(FLNM,'(I1,"/Restarts/restartMD.", I6,".dat")') MYID,ITER
           ELSEIF (ITER .GE. 1000000 .AND. ITER .LT. 10000000) THEN
              WRITE(FLNM,'(I1,"/Restarts/restartMD.", I7,".dat")') MYID,ITER
           ENDIF

        ELSEIF ( MYID .GE. 10 .AND. MYID .LT. 100) THEN

           IF (ITER .LT. 10) THEN
              WRITE(FLNM,'(I2,"/Restarts/restartMD.", I1,".dat")') MYID,ITER
           ELSEIF (ITER .GE. 10 .AND. ITER .LT. 100) THEN
              WRITE(FLNM,'(I2,"/Restarts/restartMD.", I2,".dat")') MYID,ITER
           ELSEIF (ITER .GE. 100 .AND. ITER .LT. 1000) THEN
              WRITE(FLNM,'(I2,"/Restarts/restartMD.", I3,".dat")') MYID,ITER
           ELSEIF (ITER .GE. 1000 .AND. ITER .LT. 10000) THEN
              WRITE(FLNM,'(I2,"/Restarts/restartMD.", I4,".dat")') MYID,ITER
           ELSEIF (ITER .GE. 10000 .AND. ITER .LT. 100000) THEN
              WRITE(FLNM,'(I2,"/Restarts/restartMD.", I5,".dat")') MYID,ITER
           ELSEIF (ITER .GE. 100000 .AND. ITER .LT. 1000000) THEN
              WRITE(FLNM,'(I2,"/Restarts/restartMD.", I6,".dat")') MYID,ITER
           ELSEIF (ITER .GE. 1000000 .AND. ITER .LT. 10000000) THEN
              WRITE(FLNM,'(I2,"/Restarts/restartMD.", I7,".dat")') MYID,ITER
           ENDIF

        ELSEIF ( MYID .GE. 100 .AND. MYID .LT. 1000) THEN

           IF (ITER .LT. 10) THEN
              WRITE(FLNM,'(I3,"/Restarts/restartMD.", I1,".dat")') MYID,ITER
           ELSEIF (ITER .GE. 10 .AND. ITER .LT. 100) THEN
              WRITE(FLNM,'(I3,"/Restarts/restartMD.", I2,".dat")') MYID,ITER
           ELSEIF (ITER .GE. 100 .AND. ITER .LT. 1000) THEN
              WRITE(FLNM,'(I3,"/Restarts/restartMD.", I3,".dat")') MYID,ITER
           ELSEIF (ITER .GE. 1000 .AND. ITER .LT. 10000) THEN
              WRITE(FLNM,'(I3,"/Restarts/restartMD.", I4,".dat")') MYID,ITER
           ELSEIF (ITER .GE. 10000 .AND. ITER .LT. 100000) THEN
              WRITE(FLNM,'(I3,"/Restarts/restartMD.", I5,".dat")') MYID,ITER
           ELSEIF (ITER .GE. 100000 .AND. ITER .LT. 1000000) THEN
              WRITE(FLNM,'(I3,"/Restarts/restartMD.", I6,".dat")') MYID,ITER
           ELSEIF (ITER .GE. 1000000 .AND. ITER .LT. 10000000) THEN
              WRITE(FLNM,'(I3,"/Restarts/restartMD.", I7,".dat")') MYID,ITER
           ENDIF

        ENDIF

        OPEN (UNIT = 21, STATUS="UNKNOWN", FILE=FLNM)

        IF (MYID .LT. 10) THEN
           WRITE(FLNM,'(I1,"/Restarts/restartMD.last.dat")') MYID
        ELSEIF (MYID .GE. 10 .AND. MYID .LT. 100) THEN
           WRITE(FLNM,'(I2,"/Restarts/restartMD.last.dat")') MYID
        ELSEIF (MYID .GE. 100 .AND. MYID .LT. 1000) THEN
           WRITE(FLNM,'(I3,"/Restarts/restartMD.last.dat")') MYID
        ENDIF

        OPEN (UNIT = 19, STATUS="UNKNOWN", FILE=FLNM)

     ENDIF

  ELSEIF (RELAXME .EQ. 1) THEN

     OPEN (UNIT = 21, STATUS="UNKNOWN", FILE="restartREL.dat")

  ELSEIF (RELAXME .EQ. 0 .AND. MDON .EQ. 0) THEN

     OPEN (UNIT = 21, STATUS="UNKNOWN", FILE="restart_singlepoint.dat")

  ENDIF

  IF (MDON .EQ. 1)  WRITE(21,16) "Iter= ", ITER
  WRITE(21,*) NATS
  WRITE(21,*) BOX(1,1), BOX(1,2), BOX(1,3)
  WRITE(21,*) BOX(2,1), BOX(2,2), BOX(2,3)
  WRITE(21,*) BOX(3,1), BOX(3,2), BOX(3,3)


  !  WRITE(21,10) "Nats= ", NATS
  !  WRITE(21,11) "1.0"
  !  WRITE(21,12) BOX(1,1), BOX(2,1), BOX(1,2), BOX(2,2), &
  !       BOX(1,3), BOX(2,3)

  DO I = 1, NATS
     WRITE(21,14) ATELE(I), CR(1,I), CR(2,I), CR(3,I)
  ENDDO

  WRITE(21,17) CHEMPOT

  IF (SPINON .EQ. 0) THEN

     WRITE(21,18) HDIM

     IF (KON .EQ. 0) THEN

        DO I = 1, HDIM
           WRITE(21,19) BO(I,I)/TWO, BO(I,I)/TWO
        ENDDO

     ELSE

        DO I = 1, HDIM
           KSUM = CMPLX(ZERO)

           DO K = 1, NKTOT
              KSUM = KSUM + KBO(I,I,K)
           ENDDO

           WRITE(21,19) REAL(KSUM)/(TWO*REAL(NKTOT)), &
                REAL(KSUM)/(TWO*REAL(NKTOT))
        ENDDO

     ENDIF

  ELSEIF (SPINON .EQ. 1) THEN

     WRITE(21,18) HDIM
     DO I = 1, HDIM
        WRITE(21,19) RHOUP(I,I), RHODOWN(I,I)
     ENDDO

  ENDIF

  IF (MDON .EQ. 1) THEN
     DO I = 1, NATS
        WRITE(21,15) V(1,I), V(2,I), V(3,I)
     ENDDO
  ENDIF

  CLOSE(21)

  IF (MDON .EQ. 1) THEN

     WRITE(19,16) "Iter= ", ITER
     WRITE(19,*) NATS
     WRITE(19,*) BOX(1,1), BOX(1,2), BOX(1,3)
     WRITE(19,*) BOX(2,1), BOX(2,2), BOX(2,3)
     WRITE(19,*) BOX(3,1), BOX(3,2), BOX(3,3)

     DO I = 1, NATS
        WRITE(19,14) ATELE(I), CR(1,I), CR(2,I), CR(3,I)
     ENDDO

     WRITE(19,17) CHEMPOT

     IF (SPINON .EQ. 0) THEN

        WRITE(19,18) HDIM

        IF (KON .EQ. 0) THEN

           DO I = 1, HDIM
              WRITE(19,19) BO(I,I)/TWO, BO(I,I)/TWO
           ENDDO

        ELSE

           DO I = 1, HDIM
              KSUM = CMPLX(ZERO)

              DO K = 1, NKTOT
                 KSUM = KSUM + KBO(I,I,K)
              ENDDO

              WRITE(19,19) REAL(KSUM)/(TWO*REAL(NKTOT)), &
                   REAL(KSUM)/(TWO*REAL(NKTOT))
           ENDDO

        ENDIF

     ELSEIF (SPINON .EQ. 1) THEN

        WRITE(19,18) HDIM
        DO I = 1, HDIM
           WRITE(19,19) RHOUP(I,I), RHODOWN(I,I)
        ENDDO

     ENDIF

     DO I = 1, NATS
        WRITE(19,15) V(1,I), V(2,I), V(3,I)
     ENDDO

     CLOSE(19)

  ENDIF

10 FORMAT(A6,1X,I7)
11 FORMAT(A3)
12 FORMAT(6(E24.16,1X))
  !13 FORMAT(3(E24.16,1X),A2,1X,3(I4,1X))
14 FORMAT(A2, 1X, 3G24.14)
15 FORMAT(3(E24.16,1X))
16 FORMAT(A6, 1X, I10)
17 FORMAT(E24.16)
18 FORMAT(I12)
19 FORMAT(2E24.16)

  RETURN

END SUBROUTINE WRTRESTART
