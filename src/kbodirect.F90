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

SUBROUTINE KBOEVECS

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE DIAGARRAY
  USE MYPRECISION
  USE KSPACEARRAY
  USE NEBLISTARRAY
  USE MDARRAY
#ifdef MPI_ON
  USE MPI
#endif

  IMPLICIT NONE

  INTEGER :: I, J, K
  INTEGER :: ITER, BREAKLOOP, LOOPTARGET, COUNT, MYKPOINT
  INTEGER :: CPLOC(1)
  REAL(LATTEPREC) :: OCCTARGET, OCC, FDIRAC, DFDIRAC
  REAL(LATTEPREC) :: OCCERROR, SHIFTCP, FDIRACARG, EXPARG
  REAL(LATTEPREC), PARAMETER :: MAXSHIFT = ONE
  REAL(LATTEPREC) :: S, OCCLOGOCC_ELECTRONS, OCCLOGOCC_HOLES
  REAL(LATTEPREC) :: TRACE, EBAND
  COMPLEX(LATTEPREC) :: KEBAND, ZTRACE, ZFDIRAC, ZONE
#ifdef MPI_ON
  INTEGER :: MYID, IERR, NUMPROCS, STATUS(MPI_STATUS_SIZE)
  COMPLEX(LATTEPREC), ALLOCATABLE :: TMPKBO(:,:)
  IF (EXISTERROR) RETURN

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NUMPROCS, IERR)

  ALLOCATE(TMPKBO(HDIM, HDIM))
#endif

  ZONE = CMPLX(ONE)
  KBO = CMPLX(ZERO) ! Initialize the density matrix

  OCCTARGET = BNDFIL*REAL(HDIM*NKTOT)

  !  PRINT*, TOTNE, OCCTARGET

  ITER = 0

  BREAKLOOP = 0

  OCCERROR = 1000000000.0

  !
  ! The do-while loop uses a Newton-Raphson optimization of the chemical
  ! potential to obtain the correct occupation
  !

  IF (KBT .GT. 0.001) THEN  ! This bit is for a finite electronic temperature

     DO WHILE (ABS(OCCERROR) .GT. BREAKTOL .AND. ITER .LT. 100)

        ITER = ITER + 1
        OCC = ZERO
        DFDIRAC = ZERO


        ! Loop over all k-points too

        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP SHARED(NKTOT, HDIM, KEVALS, CHEMPOT, KBT) &
        !$OMP PRIVATE(K, I, FDIRACARG, EXPARG, FDIRAC) &
        !$OMP REDUCTION(+: OCC, DFDIRAC)

        DO K = 1, NKTOT
           DO I = 1, HDIM

              FDIRACARG = (KEVALS(I, K) - CHEMPOT)/KBT

              FDIRACARG = MAX(FDIRACARG, -EXPTOL)
              FDIRACARG = MIN(FDIRACARG, EXPTOL)

              EXPARG = EXP(FDIRACARG)
              FDIRAC = ONE/(ONE + EXPARG)
              OCC = OCC + FDIRAC
              DFDIRAC = DFDIRAC + EXPARG*FDIRAC*FDIRAC

           ENDDO
        ENDDO

        !$OMP END PARALLEL DO

        DFDIRAC = DFDIRAC/REAL(NKTOT)

        DFDIRAC = DFDIRAC/KBT

        OCCERROR = (OCCTARGET - OCC)/REAL(NKTOT)

        IF (ABS(DFDIRAC) .LT. NUMLIMIT) DFDIRAC = SIGN(NUMLIMIT, DFDIRAC)

        SHIFTCP = OCCERROR/DFDIRAC

        IF (ABS(SHIFTCP) .GT. MAXSHIFT) SHIFTCP = SIGN(MAXSHIFT, SHIFTCP)

        CHEMPOT = CHEMPOT + SHIFTCP

     ENDDO

     IF (ITER .EQ. 100) THEN
        CALL ERRORS("kbodirect","Newton-Raphson scheme to find the Chemical potential does not converge")
     ENDIF

     ! Now we have the chemical potential we can build the density matrix

     S = ZERO

     IF (MDON .EQ. 0 .OR. &
          (MDON .EQ. 1 .AND. MOD(ENTROPYITER, WRTFREQ) .EQ. 0 )) THEN

        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP SHARED(NKTOT, HDIM, KEVALS, CHEMPOT, KBT) &
        !$OMP PRIVATE(K, I, FDIRACARG, FDIRAC, OCCLOGOCC_HOLES, OCCLOGOCC_ELECTRONS) &
        !$OMP REDUCTION(+: S)

        DO K = 1, NKTOT
           DO I = 1, HDIM

              FDIRACARG = (KEVALS(I, K) - CHEMPOT)/KBT

              FDIRACARG = MAX(FDIRACARG, -EXPTOL)
              FDIRACARG = MIN(FDIRACARG, EXPTOL)

              FDIRAC = ONE/(ONE + EXP(FDIRACARG))

              OCCLOGOCC_ELECTRONS = FDIRAC * LOG(FDIRAC)
              OCCLOGOCC_HOLES = (ONE - FDIRAC) * LOG(ONE - FDIRAC)

              S = S + TWO*(OCCLOGOCC_ELECTRONS + OCCLOGOCC_HOLES)

           ENDDO
        ENDDO

        !$OMP END PARALLEL DO

        S = S/REAL(NKTOT)

        ! Compute the gap only when we have to...

        ! If we have an even number of electrons

        !        IF (MOD(INT(TOTNE),2) .EQ. 0) THEN
        !           EGAP = EVALS(INT(OCCTARGET) + 1) - EVALS(INT(OCCTARGET))
        !        ELSE
        !           EGAP = ZERO
        !        ENDIF

     ENDIF

     ENTE = -KBT*S

#ifdef MPI_OFF

     DO K = 1, NKTOT
        DO I = 1, HDIM

           FDIRACARG = (KEVALS(I,K) - CHEMPOT)/KBT

           FDIRACARG = MAX(FDIRACARG, -EXPTOL)
           FDIRACARG = MIN(FDIRACARG, EXPTOL)

           ZFDIRAC = CMPLX(ONE/(ONE + EXP(FDIRACARG)))

           CALL ZGERC(HDIM, HDIM, ZFDIRAC, KEVECS(:,I,K), 1, KEVECS(:,I,K), 1, KBO(:,:,K), HDIM)

        ENDDO
     ENDDO

#elif defined(MPI_ON)

     DO K = 1, NKTOT

        IF (MOD(K,NUMPROCS) .EQ. MYID) THEN

           DO I = 1, HDIM

              FDIRACARG = (KEVALS(I,K) - CHEMPOT)/KBT

              FDIRACARG = MAX(FDIRACARG, -EXPTOL)
              FDIRACARG = MIN(FDIRACARG, EXPTOL)

              ZFDIRAC = CMPLX(ONE/(ONE + EXP(FDIRACARG)))

              CALL ZGERC(HDIM, HDIM, ZFDIRAC, KEVECS(:,I,K), 1, KEVECS(:,I,K), 1, KBO(:,:,K), HDIM)

           ENDDO
        ENDIF
     ENDDO

     !     CALL MPI_Barrier(MPI_COMM_WORLD, IERR)

     ! Collect KBO on rank 0 then broadcast

     IF (MYID .NE. 0) THEN ! Try to collect on the master

        DO I = 1, NKTOT

           ! Check whether the rank owns the data...

           IF (MOD(I,NUMPROCS) .EQ. MYID) THEN

              ! If so, sent it to the others

              CALL MPI_SEND(KBO(1,1,I), HDIM*HDIM, MPI_DOUBLE_COMPLEX, &
                   0, I, MPI_COMM_WORLD, IERR)

           ENDIF

        ENDDO

     ELSE

        ! The master receives everything and puts them in the right place
        ! using the MPI_TAG

        DO I = 1, NKTOT - (NKTOT/NUMPROCS)

           CALL MPI_RECV(TMPKBO, HDIM*HDIM, MPI_DOUBLE_COMPLEX, &
                MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, STATUS, IERR)

           MYKPOINT = STATUS(MPI_TAG)
           !        PRINT*, MYID, MYKPOINT
           DO J = 1, HDIM
              DO K = 1, HDIM
                 KBO(K,J,MYKPOINT) = TMPKBO(K,J)
              ENDDO
           ENDDO

        ENDDO

     ENDIF

     CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

     CALL MPI_BCAST(KBO, HDIM*HDIM*NKTOT, MPI_DOUBLE_COMPLEX, 0, &
          MPI_COMM_WORLD, IERR)

#endif

  ELSE ! This bit is for zero electronic temperature

     IF (MOD(INT(TOTNE),2) .NE. 0) THEN
        CALL ERRORS("kbodirect","Odd number of electrons - run a &
             & spin-polarized calculation or use a finite electron temperature")
     ENDIF

     !
     ! This definition of the chemical potential is a little arbitrary
     !

     LOOPTARGET = INT(OCCTARGET)

     ! Find the chemical potential - we need the looptarget'th lowest
     ! eigenvalue

     COUNT = 0
     DO I = 1, NKTOT
        DO J = 1, HDIM
           COUNT = COUNT + 1
           CPLIST(COUNT) = KEVALS(J,I)
        ENDDO
     ENDDO

     DO I = 1, LOOPTARGET
        CHEMPOT = MINVAL(CPLIST)
        CPLOC = MINLOC(CPLIST)
        CPLIST(CPLOC(1)) = 1.0D12 ! We do this so we don't get this eigenvalue again on the next loop
     ENDDO

     EGAP = MINVAL(CPLIST) - CHEMPOT

#ifdef MPI_OFF

     COUNT = 0
     DO K = 1, NKTOT
        DO I = 1, HDIM

           IF (KEVALS(I,K) .LE. CHEMPOT .AND. COUNT .LT. LOOPTARGET) THEN
              COUNT = COUNT + 1
              CALL ZGERC(HDIM, HDIM, ZONE, KEVECS(:,I,K), 1, &
                   KEVECS(:,I,K), 1, KBO(:,:,K), HDIM)

	   ENDIF

        ENDDO
     ENDDO

#elif defined(MPI_ON)

     DO K = 1, NKTOT

        IF (MOD(K,NUMPROCS) .EQ. MYID) THEN

           DO I = 1, HDIM

              IF (KEVALS(I,K) .LE. CHEMPOT) THEN

                 CALL ZGERC(HDIM, HDIM, ZONE, KEVECS(:,I,K), 1, &
                      KEVECS(:,I,K), 1, KBO(:,:,K), HDIM)

              ENDIF
           ENDDO
        ENDIF
     ENDDO

     !     CALL MPI_Barrier(MPI_COMM_WORLD, IERR)

     ! Collect KBO on rank 0 then broadcast

     IF (MYID .NE. 0) THEN ! Try to collect on the master

        DO I = 1, NKTOT

           ! Check whether the rank owns the data...

           IF (MOD(I,NUMPROCS) .EQ. MYID) THEN

              ! If so, sent it to the others

              CALL MPI_SEND(KBO(1,1,I), HDIM*HDIM, MPI_DOUBLE_COMPLEX, &
                   0, I, MPI_COMM_WORLD, IERR)

           ENDIF

        ENDDO

     ELSE

        ! The master receives everything and puts them in the right place
        ! using the MPI_TAG

        DO I = 1, NKTOT - (NKTOT/NUMPROCS)

           CALL MPI_RECV(TMPKBO, HDIM*HDIM, MPI_DOUBLE_COMPLEX, &
                MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, STATUS, IERR)

           MYKPOINT = STATUS(MPI_TAG)
           !        PRINT*, MYID, MYKPOINT
           DO J = 1, HDIM
              DO K = 1, HDIM
                 KBO(K,J,MYKPOINT) = TMPKBO(K,J)
              ENDDO
           ENDDO

        ENDDO

     ENDIF

     CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)

     CALL MPI_BCAST(KBO, HDIM*HDIM*NKTOT, MPI_DOUBLE_COMPLEX, 0, &
          MPI_COMM_WORLD, IERR)

#endif


  ENDIF

  !  PRINT*, COUNT, OCCTARGET

  KBO = KBO*CMPLX(TWO)

  IF (DEBUGON .EQ. 1) THEN

     OPEN(UNIT=31, STATUS="UNKNOWN", FILE="mykBO.dat")

     DO K = 1, NKTOT
        WRITE(31,*) K
        DO I = 1, HDIM
           WRITE(31,12) (KBO(I,J,K), J = 1, HDIM)
        ENDDO
     ENDDO

     CLOSE(31)

  ENDIF

12 FORMAT(100F8.3)

#ifdef MPI_ON
  DEALLOCATE(TMPKBO)
#endif

  RETURN

END SUBROUTINE KBOEVECS
