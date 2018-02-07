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

SUBROUTINE KDIAGMYH

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE DIAGARRAY
  USE SPINARRAY
  USE NONOARRAY
  USE KSPACEARRAY
  USE MYPRECISION
#ifdef MPI_ON
  USE MPI
#endif

  IMPLICIT NONE

  INTEGER :: INFO
  INTEGER :: I, J, K
#ifdef MPI_ON
  INTEGER :: MYID, IERR, NUMPROCS, MYKPOINT, STATUS(MPI_STATUS_SIZE)
  REAL(LATTEPREC), ALLOCATABLE :: TMPVALS(:)
  COMPLEX(LATTEPREC), ALLOCATABLE :: TMPVECS(:,:)
  IF (EXISTERROR) RETURN

  ALLOCATE(TMPVECS(HDIM, HDIM), TMPVALS(HDIM)) 

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NUMPROCS, IERR)

#endif


  ! Now we're doing k-space. We have to get the eigenvalues and 
  ! eigenvectors of all NKTOT of our H matrices               

  IF (BASISTYPE .EQ. "ORTHO") THEN
     KEVECS = HK
  ELSEIF (BASISTYPE .EQ. "NONORTHO") THEN
     KEVECS = KORTHOH
  ENDIF

#ifdef MPI_OFF

  DO I = 1, NKTOT

     CALL ZHEEV('V', 'U', HDIM, KEVECS(:,:,I), HDIM, KEVALS(:,I), &
          DIAG_ZWORK, DIAG_LZWORK, DIAG_RWORK, INFO)

  ENDDO


#elif defined(MPI_ON)

  DO I = 1, NKTOT

     IF (MOD(I, NUMPROCS) .EQ. MYID) THEN

        CALL ZHEEV('V', 'U', HDIM, KEVECS(:,:,I), HDIM, KEVALS(:,I), &
             DIAG_ZWORK, DIAG_LZWORK, DIAG_RWORK, INFO)

     ENDIF

  ENDDO

  !  CALL MPI_Barrier(MPI_COMM_WORLD, IERR)  

  ! Loop through all the k-points. Every processor is going to 
  ! send all the sub-matrices it computed to all the others (coz
  ! they don't have them yet)
  ! We will label the sends with the k-point ID in MPI_TAG

  IF (MYID .NE. 0) THEN ! Try to collect on the master

     DO I = 1, NKTOT

        ! Check whether the rank owns the data...

        IF (MOD(I,NUMPROCS) .EQ. MYID) THEN 

           ! If so, sent it to the others


           CALL MPI_SEND(KEVALS(1,I), HDIM, MPI_DOUBLE_PRECISION, &
                0, I, MPI_COMM_WORLD, IERR)

        ENDIF

     ENDDO

  ELSE

     DO I = 1, NKTOT - (NKTOT/NUMPROCS)

        CALL MPI_RECV(TMPVALS, HDIM, MPI_DOUBLE_PRECISION, &
             MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, STATUS, IERR)

        MYKPOINT = STATUS(MPI_TAG)
        !        PRINT*, MYID, MYKPOINT
        DO K = 1, HDIM
           KEVALS(K,MYKPOINT) = TMPVALS(K)
        ENDDO

     ENDDO

  ENDIF

  CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

  CALL MPI_BCAST(KEVALS, HDIM*NKTOT, MPI_DOUBLE_PRECISION, 0, &
       MPI_COMM_WORLD, IERR)

  IF (MYID .NE. 0) THEN ! Try to collect on the master

     DO I = 1, NKTOT

        ! Check whether the rank owns the data...

        IF (MOD(I,NUMPROCS) .EQ. MYID) THEN 

           ! If so, sent it to the myid = 0


           CALL MPI_SEND(KEVECS(1,1,I), HDIM*HDIM, MPI_DOUBLE_COMPLEX, &
                0, I, MPI_COMM_WORLD, IERR)

        ENDIF

     ENDDO

  ELSE

     ! myid = 0 receives all the data

     DO I = 1, NKTOT - (NKTOT/NUMPROCS)

        CALL MPI_RECV(TMPVECS, HDIM*HDIM, MPI_DOUBLE_COMPLEX, &
             MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, STATUS, IERR)

        MYKPOINT = STATUS(MPI_TAG)
        !        PRINT*, MYID, MYKPOINT
        DO J = 1, HDIM
           DO K = 1, HDIM
              KEVECS(K,J,MYKPOINT) = TMPVECS(K,J)
           ENDDO
        ENDDO

     ENDDO

  ENDIF

  CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

  CALL MPI_BCAST(KEVECS, HDIM*HDIM*NKTOT, MPI_DOUBLE_COMPLEX, 0, &
       MPI_COMM_WORLD, IERR)

  DEALLOCATE(TMPVECS, TMPVALS) 

#endif

  RETURN

END SUBROUTINE KDIAGMYH
