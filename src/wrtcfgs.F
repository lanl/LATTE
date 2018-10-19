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

SUBROUTINE WRTCFGS(ITER)

  !
  ! This subroutine will write out all quantities during calculated
  ! during a MD run to .cfg files so they can be read into
  ! Ju Li's most awesome (and most free) Atomeye visualization code directly.
  !
  ! http://mt.seas.upenn.edu/Archive/Graphics/A/
  !
  ! I suppose that coordinates, velocity, forces, atomic temperature,
  ! partials charges, and atomic spin densities would come handy
  !

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MDARRAY
  USE SPINARRAY
  USE MYPRECISION
#ifdef MPI_ON
  USE MPI
#endif

  IMPLICIT NONE

  INTEGER :: I, J, L
  INTEGER :: ITER, NOMYELEMENTS, MYINDEX
  INTEGER :: MYID, IERR, IPIV(3), INFO
  REAL(LATTEPREC), ALLOCATABLE :: T(:), MAGF(:), ATSPIN(:)
  REAL(LATTEPREC) :: V2, MYMASS, MASSMYELEMENTS(4)
  REAL(LATTEPREC) :: WORK(3), BOXINV(3,3), S(3)
  !  REAL(LATTEPREC), PARAMETER :: CONV = 1.66053D6/(3.0D0*1.38062D0)
  CHARACTER(LEN=50) :: FLNM
  CHARACTER(LEN=2) :: MYELEMENTS(92)
  IF (EXISTERROR) RETURN

  ALLOCATE(T(NATS), MAGF(NATS))

  !
  ! Let's first get the temperature and |F| for each atom
  !


  IF (MDON .EQ. 1) THEN

     IF (.NOT. ALLOCATED(V)) THEN
        ALLOCATE(V(3,NATS)); V = 0.0d0
     ENDIF

     DO I = 1, NATS

        V2 = V(1,I)*V(1,I) + V(2,I)*V(2,I) + V(3,I)*V(3,I)

        T(I) = MASS(ELEMPOINTER(I))*V2

        MAGF(I) = SQRT(FTOT(1,I)*FTOT(1,I) + FTOT(2,I)*FTOT(2,I) + &
             FTOT(3,I)*FTOT(3,I))

     ENDDO

  ELSE

     T = ZERO

     DO I = 1, NATS

        MAGF(I) = SQRT(FTOT(1,I)*FTOT(1,I) + FTOT(2,I)*FTOT(2,I) + &
             FTOT(3,I)*FTOT(3,I))

     ENDDO

  ENDIF

  T = T*MVV2T/THREE

  !
  ! Now let's figure out how many elements we've got and the cfg
  ! formats is most efficient when they're all written out in the same
  ! block
  !

  MYELEMENTS = "XX"
  NOMYELEMENTS = 1
  MYELEMENTS(NOMYELEMENTS) = ATELE(1)
  MASSMYELEMENTS(NOMYELEMENTS) = MASS(ELEMPOINTER(1))

  DO J = 2, NATS

     IF (ATELE(J) .NE. MYELEMENTS(1) &
          .AND.  ATELE(J) .NE. MYELEMENTS(2) &
          .AND.  ATELE(J) .NE. MYELEMENTS(3) &
          .AND.  ATELE(J) .NE. MYELEMENTS(4) ) THEN

        NOMYELEMENTS = NOMYELEMENTS + 1
        MYELEMENTS(NOMYELEMENTS) = ATELE(J)
        MASSMYELEMENTS(NOMYELEMENTS) = MASS(ELEMPOINTER(J))

     ENDIF

  ENDDO

  IF (PBCON .EQ. 1) THEN

     ! Write CFG files if periodic boundary conditions are on

     IF (ITER .LT. 0) THEN
        WRITE(FLNM,'("animate/dump2atomeye.PANIC.cfg")')
        IF (ITER .EQ. -999) WRITE(FLNM,'("lastsystem.cfg")')
     ENDIF

     IF (PARREP .EQ. 0) THEN

        IF (ITER .GE. 0 .AND. ITER .LT. 10) THEN
           WRITE(FLNM,'("animate/dump2atomeye.",I1,".cfg")') ITER
        ELSEIF (ITER .GE. 10 .AND. ITER .LT. 100) THEN
           WRITE(FLNM,'("animate/dump2atomeye.",I2,".cfg")') ITER
        ELSEIF (ITER .GE. 100 .AND. ITER .LT. 1000) THEN
           WRITE(FLNM,'("animate/dump2atomeye.",I3,".cfg")') ITER
        ELSEIF (ITER .GE. 1000 .AND. ITER .LT. 10000) THEN
           WRITE(FLNM,'("animate/dump2atomeye.",I4,".cfg")') ITER
        ELSEIF (ITER .GE. 10000 .AND. ITER .LT. 100000) THEN
           WRITE(FLNM,'("animate/dump2atomeye.",I5,".cfg")') ITER
        ELSEIF (ITER .GE. 100000 .AND. ITER .LT. 1000000) THEN
           WRITE(FLNM,'("animate/dump2atomeye.",I6,".cfg")') ITER
        ELSEIF (ITER .GE. 1000000 .AND. ITER .LT. 10000000) THEN
           WRITE(FLNM,'("animate/dump2atomeye.",I7,".cfg")') ITER
        ELSEIF (ITER .GE. 10000000 .AND. ITER .LT. 100000000) THEN
           WRITE(FLNM,'("animate/dump2atomeye.",I8,".cfg")') ITER
        ENDIF

        OPEN(UNIT = 23, STATUS="UNKNOWN", FILE=FLNM)

     ELSE

#ifdef MPI_ON

        CALL MPI_COMM_RANK( MPI_COMM_WORLD, MYID, IERR )

#endif

        IF (MYID .LT. 10) THEN

           IF (ITER .GE. 0 .AND. ITER .LT. 10) THEN
              WRITE(FLNM,'(I1,"/animate/dump2atomeye.",I1,".cfg")') MYID, ITER
           ELSEIF (ITER .GE. 10 .AND. ITER .LT. 100) THEN
              WRITE(FLNM,'(I1,"/animate/dump2atomeye.",I2,".cfg")') MYID,ITER
           ELSEIF (ITER .GE.100 .AND. ITER .LT. 1000) THEN
              WRITE(FLNM,'(I1,"/animate/dump2atomeye.",I3,".cfg")') MYID,ITER
           ELSEIF (ITER .GE. 1000 .AND. ITER .LT. 10000) THEN
              WRITE(FLNM,'(I1,"/animate/dump2atomeye.",I4,".cfg")') MYID,ITER
           ELSEIF (ITER .GE. 10000 .AND. ITER .LT. 100000) THEN
              WRITE(FLNM,'(I1,"/animate/dump2atomeye.",I5,".cfg")') MYID,ITER
           ELSEIF (ITER .GE. 100000 .AND. ITER .LT. 1000000) THEN
              WRITE(FLNM,'(I1,"/animate/dump2atomeye.",I6,".cfg")') MYID,ITER
           ELSEIF (ITER .GE. 1000000 .AND. ITER .LT. 10000000) THEN
              WRITE(FLNM,'(I1,"/animate/dump2atomeye.",I7,".cfg")') MYID,ITER
           ELSEIF (ITER .GE. 10000000 .AND. ITER .LT. 100000000) THEN
              WRITE(FLNM,'(I1,"/animate/dump2atomeye.",I8,".cfg")') MYID,ITER
           ENDIF

        ELSEIF (MYID .GE. 10 .AND. MYID .LT. 100) THEN

           IF (ITER .GE. 0 .AND. ITER .LT. 10) THEN
              WRITE(FLNM,'(I2,"/animate/dump2atomeye.",I1,".cfg")') MYID, ITER
           ELSEIF (ITER .GE. 10 .AND. ITER .LT. 100) THEN
              WRITE(FLNM,'(I2,"/animate/dump2atomeye.",I2,".cfg")') MYID,ITER
           ELSEIF (ITER .GE.100 .AND. ITER .LT. 1000) THEN
              WRITE(FLNM,'(I2,"/animate/dump2atomeye.",I3,".cfg")') MYID,ITER
           ELSEIF (ITER .GE. 1000 .AND. ITER .LT. 10000) THEN
              WRITE(FLNM,'(I2,"/animate/dump2atomeye.",I4,".cfg")') MYID,ITER
           ELSEIF (ITER .GE. 10000 .AND. ITER .LT. 100000) THEN
              WRITE(FLNM,'(I2,"/animate/dump2atomeye.",I5,".cfg")') MYID,ITER
           ELSEIF (ITER .GE. 100000 .AND. ITER .LT. 1000000) THEN
              WRITE(FLNM,'(I2,"/animate/dump2atomeye.",I6,".cfg")') MYID,ITER
           ELSEIF (ITER .GE. 1000000 .AND. ITER .LT. 10000000) THEN
              WRITE(FLNM,'(I2,"/animate/dump2atomeye.",I7,".cfg")') MYID,ITER
           ELSEIF (ITER .GE. 10000000 .AND. ITER .LT. 100000000) THEN
              WRITE(FLNM,'(I2,"/animate/dump2atomeye.",I8,".cfg")') MYID,ITER
           ENDIF

        ELSEIF (MYID .GE. 100 .AND. MYID .LT. 1000) THEN

           IF (ITER .GE. 0 .AND. ITER .LT. 10) THEN
              WRITE(FLNM,'(I3,"/animate/dump2atomeye.",I1,".cfg")') MYID, ITER
           ELSEIF (ITER .GE. 10 .AND. ITER .LT. 100) THEN
              WRITE(FLNM,'(I3,"/animate/dump2atomeye.",I2,".cfg")') MYID,ITER
           ELSEIF (ITER .GE.100 .AND. ITER .LT. 1000) THEN
              WRITE(FLNM,'(I3,"/animate/dump2atomeye.",I3,".cfg")') MYID,ITER
           ELSEIF (ITER .GE. 1000 .AND. ITER .LT. 10000) THEN
              WRITE(FLNM,'(I3,"/animate/dump2atomeye.",I4,".cfg")') MYID,ITER
           ELSEIF (ITER .GE. 10000 .AND. ITER .LT. 100000) THEN
              WRITE(FLNM,'(I3,"/animate/dump2atomeye.",I5,".cfg")') MYID,ITER
           ELSEIF (ITER .GE. 100000 .AND. ITER .LT. 1000000) THEN
              WRITE(FLNM,'(I3,"/animate/dump2atomeye.",I6,".cfg")') MYID,ITER
           ELSEIF (ITER .GE. 1000000 .AND. ITER .LT. 10000000) THEN
              WRITE(FLNM,'(I3,"/animate/dump2atomeye.",I7,".cfg")') MYID,ITER
           ELSEIF (ITER .GE. 10000000 .AND. ITER .LT. 100000000) THEN
              WRITE(FLNM,'(I3,"/animate/dump2atomeye.",I8,".cfg")') MYID,ITER
           ENDIF

        ENDIF

        OPEN(UNIT = 23, STATUS="UNKNOWN", FILE=FLNM)

     ENDIF


     WRITE(23,10) "Number of particles = ", NATS
     WRITE(23,'("A = 1.0 Angstrom")')
     WRITE(23,11) "H0(1,1) = ", BOX(1,1), " A"
     WRITE(23,11) "H0(1,2) = ", BOX(1,2), " A"
     WRITE(23,11) "H0(1,3) = ", BOX(1,3), " A"
     WRITE(23,11) "H0(2,1) = ", BOX(2,1), " A"
     WRITE(23,11) "H0(2,2) = ", BOX(2,2), " A"
     WRITE(23,11) "H0(2,3) = ", BOX(2,3), " A"
     WRITE(23,11) "H0(3,1) = ", BOX(3,1), " A"
     WRITE(23,11) "H0(3,2) = ", BOX(3,2), " A"
     WRITE(23,11) "H0(3,3) = ", BOX(3,3), " A"

     ! Box inverse

     BOXINV = BOX

     CALL DGETRF(3, 3, BOXINV, 3, IPIV, INFO)

     CALL DGETRI(3, BOXINV, 3, IPIV, WORK, 3, INFO)


     !
     ! No partial charges, no spins
     !
     ! Write: x, y, z, vx, vy, vz, fx, fy, fz, |F|, temperature
     !

     IF (MDON .EQ. 1) THEN

        IF(SPINON .EQ. 0) THEN

           !
           ! Same as above, but now we add partial charges too
           !

           WRITE(23,'("entry_count = 12")')
           WRITE(23,13) "auxiliary[0] = fx "
           WRITE(23,13) "auxiliary[1] = fy "
           WRITE(23,13) "auxiliary[2] = fz  "
           WRITE(23,13) "auxiliary[3] = |F|"
           WRITE(23,13) "auxiliary[4] = T/K"
           WRITE(23,13) "auxiliary[5] = q  "

           DO I = 1, NOMYELEMENTS

              WRITE(23,14) MASSMYELEMENTS(I)
              WRITE(23,15) MYELEMENTS(I)

              DO J = 1, NATS

                 IF (ATELE(J) .EQ. MYELEMENTS(I)) THEN

                    CALL DGEMV('T', 3, 3, ONE, BOXINV, 3, CR(1,J), 1, ZERO, S, 1)
                    WRITE(23,17) S(1), S(2), S(3), V(1,J), V(2,J), V(3,J), &
                         FTOT(1,J), FTOT(2,J), FTOT(3,J), &
                         MAGF(J), T(J), DELTAQ(J)

                 ENDIF

              ENDDO

           ENDDO

        ELSEIF (SPINON .EQ. 1) THEN

           !
           ! Now we add the net spin per atom too
           !

           ALLOCATE(ATSPIN(NATS))

           MYINDEX = 0

           DO I = 1, NATS

              IF (BASIS(ELEMPOINTER(I)) .EQ. "s") THEN

                 MYINDEX = MYINDEX + 1
                 ATSPIN(I) = DELTASPIN(MYINDEX)

              ELSEIF (BASIS(ELEMPOINTER(I)) .EQ. "sp") THEN

                 ATSPIN(I) = DELTASPIN(MYINDEX + 1) + DELTASPIN(MYINDEX + 2)

                 MYINDEX = MYINDEX + 2

              ENDIF

           ENDDO

           WRITE(23,'("entry_count = 10")')
           WRITE(23,13) "auxiliary[0] = |F|"
           WRITE(23,13) "auxiliary[1] = T/K"
           WRITE(23,13) "auxiliary[2] = q  "
           WRITE(23,13) "auxiliary[3] = spn"

           DO I = 1, NOMYELEMENTS

              WRITE(23,14) MASSMYELEMENTS(I)
              WRITE(23,15) MYELEMENTS(I)

              DO J = 1, NATS

                 IF (ATELE(J) .EQ. MYELEMENTS(I)) THEN

	            CALL DGEMV('T', 3, 3, ONE, BOXINV, 3, CR(1,J), 1, ZERO, S, 1)

                    WRITE(23,18) S(1), S(2), S(3), V(1,J), V(2,J), V(3,J), &
                         MAGF(J), T(J), DELTAQ(J), ATSPIN(J)

                 ENDIF

              ENDDO

           ENDDO

           DEALLOCATE(ATSPIN)

        ENDIF

     ELSE

        ! No MD so no velocities

        IF (SPINON .EQ. 0) THEN

           WRITE(23,'(".NO_VELOCITY.")')
           WRITE(23,'("entry_count = 8")')
           WRITE(23,13) "auxiliary[0] = fx "
           WRITE(23,13) "auxiliary[1] = fy "
           WRITE(23,13) "auxiliary[2] = fz "
           WRITE(23,13) "auxiliary[3] = |F|"
           WRITE(23,13) "auxiliary[4] = q  "

           DO I = 1, NOMYELEMENTS

              WRITE(23,14) MASSMYELEMENTS(I)
              WRITE(23,15) MYELEMENTS(I)

              DO J = 1, NATS

                 IF (ATELE(J) .EQ. MYELEMENTS(I)) THEN

		    CALL DGEMV('T', 3, 3, ONE, BOXINV, 3, CR(1,J), 1, ZERO, S, 1)

                    WRITE(23,22) S(1), S(2), S(3), &
                         FTOT(1,J), FTOT(2,J), FTOT(3,J), &
                         MAGF(J), DELTAQ(J)

                 ENDIF

              ENDDO

           ENDDO

        ELSEIF (SPINON .EQ. 1) THEN

           WRITE(23,'(".NO_VELOCITY.")')
           WRITE(23,'("entry_count = 9")')
           WRITE(23,13) "auxiliary[0] = fx "
           WRITE(23,13) "auxiliary[1] = fy "
           WRITE(23,13) "auxiliary[2] = fz "
           WRITE(23,13) "auxiliary[3] = |F|"
           WRITE(23,13) "auxiliary[4] = q  "
           WRITE(23,13) "auxiliary[5] = spn"

           ALLOCATE(ATSPIN(NATS))

           MYINDEX = 0

           DO I = 1, NATS

              IF (BASIS(ELEMPOINTER(I)) .EQ. "s") THEN

                 MYINDEX = MYINDEX + 1
                 ATSPIN(I) = DELTASPIN(MYINDEX)

              ELSEIF (BASIS(ELEMPOINTER(I)) .EQ. "sp") THEN

                 ATSPIN(I) = DELTASPIN(MYINDEX + 1) + DELTASPIN(MYINDEX + 2)

                 MYINDEX = MYINDEX + 2

              ENDIF

           ENDDO

           DO I = 1, NOMYELEMENTS

              WRITE(23,14) MASSMYELEMENTS(I)
              WRITE(23,15) MYELEMENTS(I)

              DO J = 1, NATS

                 IF (ATELE(J) .EQ. MYELEMENTS(I)) THEN

		    CALL DGEMV('T', 3, 3, ONE, BOXINV, 3, CR(1,J), 1, ZERO, S, 1)

                    WRITE(23,23) S(1), S(2), S(3), &
                         FTOT(1,J), FTOT(2,J), FTOT(3,J), &
                         MAGF(J), DELTAQ(J), ATSPIN(J)

                 ENDIF

              ENDDO

           ENDDO

           DEALLOCATE(ATSPIN)

        ENDIF

     ENDIF

     !     DEALLOCATE(T, MAGF)

     CLOSE(23)

  ELSE

     ! No PBCs? write an .xyz file instead

     IF (ITER .EQ. 0) THEN

        WRITE(FLNM,'("animate/myXYZfile.xyz")')

        OPEN(UNIT=23, STATUS="UNKNOWN", FILE=FLNM)

     ELSEIF (ITER .LT. 0) THEN

        WRITE(FLNM,'("animate/myXYZfile.PANIC.xyz")')
        IF (ITER .EQ. -999) WRITE(FLNM,'("lastsystem.xyz")')

        OPEN(UNIT=23, STATUS="UNKNOWN", FILE=FLNM)

     ENDIF



     WRITE(23,20) NATS
     IF (ITER .EQ. -999) THEN
        WRITE(23,'("Last atomic configuration")')
     ELSE
        WRITE(23,'("LATTE MD: dt = ", F8.3, " Time step = ", I9)') DT, ITER
     ENDIF

     DO I = 1, NATS
        WRITE(23,21) ATELE(I), CR(1,I), CR(2,I), CR(3,I)
     ENDDO

     IF (ITER .EQ. -999) CLOSE(23)

     !     IF (ITER .LT. 0) THEN

     !        WRITE(FLNM,'("animate/myXYZfile.PANIC.xyz")')
     !        IF (ITER .EQ. -999) WRITE(FLNM,'("lastsystem.xyz")')

     !        OPEN(UNIT=24, STATUS="UNKNOWN", FILE=FLNM)

     !        WRITE(24,20) NATS
     !        WRITE(24,'("LATTE MD: dt = ", F8.3, " Time step = ", I9)') DT, ITER
     !        DO I = 1, NATS
     !           WRITE(24,21) ATELE(I), CR(1,I), CR(2,I), CR(3,I)
     !        ENDDO

     !        CLOSE (24)

     !     ENDIF

  END IF

  IF (MDON .EQ. 1) THEN
     DEALLOCATE(T, MAGF)
  ELSE
     DEALLOCATE(MAGF)
  ENDIF

10 FORMAT(A22, 1X, I6)
11 FORMAT(A10,1X,F18.9, A2)
12 FORMAT(A15)
13 FORMAT(A18)
14 FORMAT(F5.1)
15 FORMAT(A2)
16 FORMAT(3(F8.6,1X), 3(F12.8,1X), F12.6, 1X, F12.2)
17 FORMAT(3(F8.6,1X), 6(F15.8,1X), F12.6, 1X, F12.2, 1X, F10.4)
18 FORMAT(3(F8.6,1X), 3(F12.8,1X), F12.6, 1X, F12.2, 1X, 2(F10.4, 1X))
20 FORMAT(I6)
21 FORMAT(A2, 1X, 3F12.2)

22 FORMAT(3(G14.6,1X), 3(G14.6,1X), G14.6, 1X, 1(G14.6, 1X))
23 FORMAT(3(G14.6,1X), 3(G14.6,1X), G14.6, 1X, 2(G14.6, 1X))


  RETURN

END SUBROUTINE WRTCFGS
