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

PROGRAM LATTE

  USE CONSTANTS_MOD
  USE TIMER_MOD
  USE SETUPARRAY
  USE PPOTARRAY
  USE PUREARRAY
  USE COULOMBARRAY
  USE SPINARRAY
  USE SPARSEARRAY
  USE MYPRECISION
  USE VIRIALARRAY
  USE DIAGARRAY
  USE KSPACEARRAY
  USE LATTEPARSER_LATTE_MOD
#ifdef MPI_ON
  USE MPI
#endif
#ifdef DBCSR_ON
  USE DBCSR_VAR_MOD
#endif

#ifdef PROGRESSON
  USE PRG_SYSTEM_MOD ! FROM PROGRESS
  USE PRG_PULAYMIXER_MOD
  USE MIXER_MOD
  USE BML
#endif


  IMPLICIT NONE

  INTEGER :: I, J
  INTEGER :: START_CLOCK, STOP_CLOCK, CLOCK_RATE, CLOCK_MAX
  INTEGER :: MYID = 0
  REAL :: TARRAY(2), RESULT, SYSTDIAG, SYSTPURE
  CHARACTER(LEN=50) :: FLNM

#ifdef MPI_ON
  INTEGER :: IERR, STATUS(MPI_STATUS_SIZE), NUMPROCS

  CALL MPI_INIT( IERR )
  CALL MPI_COMM_RANK( MPI_COMM_WORLD, MYID, IERR )
  CALL MPI_COMM_SIZE( MPI_COMM_WORLD, NUMPROCS, IERR )
#endif

  NUMSCF = 0
  CHEMPOT = ZERO

  ! Start timers
  TX = INIT_TIMER()
  TX = START_TIMER(LATTE_TIMER)
  IF (VERBOSE >= 1) CALL TIMEDATE_TAG("# LATTE started at : ")

  INQUIRE( FILE="latte.in", exist=LATTEINEXISTS )
  IF (LATTEINEXISTS) THEN
     IF(.NOT. LIBINIT) CALL PARSE_CONTROL("latte.in")

#ifdef PROGRESSON
     CALL PRG_PARSE_MIXER(MX,"latte.in")
#endif

  ELSE
     CALL READCONTROLS
  ENDIF

  CALL READTB

  IF (RESTART .EQ. 0) THEN
     CALL READCR
  ELSE
     CALL READRESTART
  ENDIF

  IF (PPOTON .EQ. 1) CALL READPPOT
  IF (PPOTON .EQ. 2) CALL READPPOTTAB
  IF (PPOTON .EQ. 3) CALL READPPOTSPLINE

  IF (DEBUGON .EQ. 1) THEN
     CALL PLOTUNIV
     IF (PPOTON .EQ. 1) CALL PLOTPPOT
  ENDIF

  CALL GETHDIM

  CALL GETMATINDLIST

  CALL RHOZERO

  CALL GETBNDFIL()

#ifdef GPUON

  CALL initialize( NGPU )

#endif

#ifdef DBCSR_ON


  IF (CONTROL .EQ. 2 .AND. SPARSEON .EQ. 1) CALL INIT_DBCSR

#endif

  IF (KBT .LT. 0.0000001 .OR. CONTROL .EQ. 2) ENTE = ZERO

  IF (MDON .EQ. 0 .AND. RELAXME .EQ. 0 .AND. DOSFITON .EQ. 0 &
       .AND. PPFITON .EQ. 0 .AND. ALLFITON .EQ. 0) THEN

     !
     ! Start the timers
     !

     CALL SYSTEM_CLOCK(START_CLOCK, CLOCK_RATE, CLOCK_MAX)
     CALL DTIME(TARRAY, RESULT)

     ! Set up neighbor lists for building the H and pair potentials

     CALL ALLOCATENEBARRAYS

     IF (ELECTRO .EQ. 1) THEN

        CALL ALLOCATECOULOMB

        CALL INITCOULOMB

     ENDIF

     IF (BASISTYPE .EQ. "NONORTHO") CALL ALLOCATENONO

     CALL NEBLISTS(0)

     ! Build the charge independent H matrix

     IF (KON .EQ. 0) THEN

        IF (SPONLY .EQ. 0) THEN
           CALL BLDNEWHS_SP
        ELSE
           CALL BLDNEWHS
        ENDIF

     ELSE

        CALL KBLDNEWH

     ENDIF


     !
     ! If we're starting from a restart file, we need to modify H such
     ! that it agrees with the density matrix elements read from file
     !

     IF (RESTART .EQ. 1) CALL IFRESTART


     !
     ! See whether we need spin-dependence too
     !

     IF (SPINON .EQ. 1) THEN
        CALL GETDELTASPIN
        CALL BLDSPINH
     ENDIF

     IF (CONTROL .EQ. 1) THEN
        CALL ALLOCATEDIAG
     ELSEIF (CONTROL .EQ. 2 .OR. CONTROL .EQ. 4 .OR. CONTROL .EQ. 5) THEN
        CALL ALLOCATEPURE
     ELSEIF (CONTROL .EQ. 3) THEN
        CALL FERMIALLOCATE
     ENDIF

     IF (CONTROL .EQ. 5) THEN

        CALL GERSHGORIN
        CALL SP2FERMIINIT

     ENDIF

     IF (ELECTRO .EQ. 0) CALL QNEUTRAL(0,1) ! Local charge neutrality

     IF (ELECTRO .EQ. 1) CALL QCONSISTENCY(0,1) ! Self-consistent charges

     ! We have to build our NKTOT complex H matrices and compute the
     ! self consistent density matrix

     ! Tr[rho dH/dR], Pulay force, and Tr[rho H] need to de-orthogonalized rho

     IF (KON .EQ. 1) CALL KGETDOS

     IF (DEBUGON .EQ. 1 .AND. SPINON .EQ. 0 .AND. KON .EQ. 0) THEN

        PRINT*, "Caution - you're writing to file the density matrix!"

        OPEN(UNIT=31, STATUS="UNKNOWN", FILE="myrho.dat")

        DO I = 1, HDIM
           WRITE(31,10) (BO(I,J), J = 1, HDIM)
        ENDDO

        CLOSE(31)

10      FORMAT(100G18.8)

     ENDIF

     IF (COMPFORCE .EQ. 1) CALL GETFORCE

     EREP = ZERO
     IF (PPOTON .EQ. 1) THEN
        CALL PAIRPOT
     ENDIF

     IF (PPOTON .EQ. 2) THEN
        CALL PAIRPOTTAB
     ENDIF

     IF (PPOTON .EQ. 3) THEN
        CALL PAIRPOTSPLINE
     ENDIF

     CALL TOTENG
     
     ECOUL = ZERO
     IF (ELECTRO .EQ. 1) CALL GETCOULE

     ESPIN = ZERO
     IF (SPINON .EQ. 1) CALL GETSPINE

     IF (CONTROL .NE. 1 .AND. CONTROL .NE. 2 .AND. KBT .GT. 0.000001 ) THEN

 	! We get the entropy automatically when using diagonalization.
        ! This is only required when employing the recursive expansion
        ! of the Fermi-operator at finite electronic temperature

        CALL ENTROPY

     ENDIF

     CALL WRTRESTART(0)

     IF (CONTROL .EQ. 1) THEN
        !        CALL DEALLOCATEDIAG
     ELSEIF (CONTROL .EQ. 2 .OR. CONTROL .EQ. 4 .OR. CONTROL .EQ. 5) THEN
        CALL DEALLOCATEPURE
     ELSEIF (CONTROL .EQ. 3) THEN
        CALL FERMIDEALLOCATE
     ENDIF

     !
     ! Stop the clocks
     !

     TX = STOP_TIMER(LATTE_TIMER)
     CALL DTIME(TARRAY, RESULT)
     CALL SYSTEM_CLOCK(STOP_CLOCK, CLOCK_RATE, CLOCK_MAX)

     CALL GETPRESSURE

!     WRITE(6,*) F(1,1), F(2,1), F(3,1)
!     WRITE(6,*) FPP(1,1), FPP(2,1), FPP(3,1)
!     WRITE(6,*) FSLCN(1,1), FSLCN(2,1), FSLCN(3,1)
!     WRITE(6,*) FPUL(1,1), FPUL(2,1), FPUL(3,1)


     !     WRITE(6,*) "Force ", FPP(1,1), FPP(2,1), FPP(3,1)
     !     PRINT*, "PCHECK ", (1.0/3.0)*(VIRBOND(1)+VIRBOND(2) + VIRBOND(3)), &
     !          (1.0/3.0)*(VIRCOUL(1)+VIRCOUL(2) + VIRCOUL(3)), &
     !          (1.0/3.0)*(VIRPAIR(1)+VIRPAIR(2) + VIRPAIR(3)), &
     !          (1.0/3.0)*(VIRPUL(1)+VIRPUL(2) + VIRPUL(3)), &
     !          (1.0/3.0)*(VIRSCOUL(1)+VIRSCOUL(2) + VIRSCOUL(3))

#ifdef DBCSR_ON

     IF (CONTROL .EQ. 2 .AND. SPARSEON .EQ. 1 .AND.  MYNODE .EQ. 0) THEN

#endif

        IF (MYID .EQ. 0) THEN
           CALL FITTINGOUTPUT(0) ! This has to come first (MJC)
           CALL SUMMARY

           !     IF (SPINON .EQ. 0) CALL NORMS

           PRINT*, "# System time  = ", TARRAY(1)
           PRINT*, "# Wall time = ", FLOAT(STOP_CLOCK - START_CLOCK)/FLOAT(CLOCK_RATE)
           PRINT*, "# Wall time per SCF =", &
                FLOAT(STOP_CLOCK - START_CLOCK)/(FLOAT(CLOCK_RATE)*FLOAT(NUMSCF))
           !     PRINT*, HDIM, FLOAT(STOP_CLOCK - START_CLOCK)/FLOAT(CLOCK_RATE)
           TX = TIMER_RESULTS()
           PRINT*, "# NUMSCF = ", NUMSCF

        ENDIF
#ifdef DBCSR_ON

     ENDIF

#endif


     !     CALL ASSESSOCC

     IF (ELECTRO .EQ. 1) CALL DEALLOCATECOULOMB

     IF (BASISTYPE .EQ. "NONORTHO") CALL DEALLOCATENONO

     CALL DEALLOCATENEBARRAYS

  ELSEIF (MDON .EQ. 1 .AND. RELAXME .EQ. 0) THEN

     IF (BASISTYPE .EQ. "NONORTHO") CALL ALLOCATENONO

     IF (XBOON .EQ. 1) CALL ALLOCATEXBO

     IF (ELECTRO .EQ. 1) THEN
        CALL ALLOCATECOULOMB
        CALL INITCOULOMB
     ENDIF

     ! Start the timers

     CALL SYSTEM_CLOCK(START_CLOCK, CLOCK_RATE, CLOCK_MAX)
     CALL DTIME(TARRAY, RESULT)

     !
     ! Call TBMD
     !

     CALL TBMD


#ifdef MPI_ON
     IF (PARREP .EQ. 1) CALL MPI_BARRIER (MPI_COMM_WORLD, IERR )
#endif


     ! Stop the timers

     CALL DTIME(TARRAY, RESULT)
     CALL SYSTEM_CLOCK(STOP_CLOCK, CLOCK_RATE, CLOCK_MAX)

     CALL SUMMARY

     IF (PBCON .EQ. 0) CLOSE(23)

     IF (BASISTYPE .EQ. "NONORTHO") CALL DEALLOCATENONO

     IF (XBOON .EQ. 1) CALL DEALLOCATEXBO

     IF (ELECTRO .EQ. 1) CALL DEALLOCATECOULOMB

     !     SYSTPURE = TARRAY(1)
     !     WRITE(6,'("# System time for MD run = ", F12.2, " s")') SYSTPURE
     WRITE(6,'("# Wall time for MD run = ", F12.2, " s")') &
          FLOAT(STOP_CLOCK - START_CLOCK)/FLOAT(CLOCK_RATE)

  ELSEIF (MDON .EQ. 0 .AND. RELAXME .EQ. 1) THEN

     CALL MSRELAX

  ELSEIF (MDON .EQ. 0 .AND. RELAXME .EQ. 0 .AND. DOSFITON .EQ. 1) THEN

     IF (SCLTYPE .EQ. "EXP") THEN
        CALL DOSFIT
     ELSE
        CALL DOSFITTAB
     ENDIF

     IF (KON .EQ. 1) CALL KGETDOS

  ELSEIF  (MDON .EQ. 0 .AND. RELAXME .EQ. 0 .AND. DOSFITON .EQ. 2) THEN


     CALL MOFIT

  ELSEIF (MDON .EQ. 0 .AND. RELAXME .EQ. 0 .AND. DOSFITON .EQ. 3) THEN

     CALL MOFITPLATO

  ELSEIF (MDON .EQ. 0 .AND. RELAXME .EQ. 0 .AND. PPFITON .EQ. 1) THEN

     CALL PPFIT

  ELSEIF (MDON .EQ. 0 .AND. RELAXME .EQ. 0 .AND. ALLFITON .EQ. 1) THEN

     CALL ALLFIT

  ELSE

     CALL ERRORS("latte","You can't have RELAXME = 1 and MDON = 1")

  ENDIF

#ifdef GPUON

  CALL shutdown()

#endif


  CALL DEALLOCATEALL

#ifdef DBCSR_ON

  !ends mpi

  IF (CONTROL .EQ. 2 .AND. SPARSEON .EQ. 1) CALL SHUTDOWN_DBCSR

#endif

  ! Done with timers
  TX = STOP_TIMER(LATTE_TIMER)
  TX = TIMER_RESULTS()
  TX = SHUTDOWN_TIMER()

#ifdef MPI_ON
  CALL MPI_FINALIZE( IERR )
#endif

END PROGRAM LATTE
