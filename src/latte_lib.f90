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

MODULE LATTE_LIB

  USE CONSTANTS_MOD
  USE TIMER_MOD
  USE SETUPARRAY
  USE PPOTARRAY
  USE PUREARRAY
  USE COULOMBARRAY
  USE SPINARRAY
  USE SPARSEARRAY
  USE MDARRAY
  USE MYPRECISION
  USE VIRIALARRAY
  USE DIAGARRAY
  USE KSPACEARRAY
  USE LATTEPARSER_LATTE_MOD
  USE NEBLISTARRAY
  USE NONOARRAY

#ifdef PROGRESSON
  USE SYSTEM_MOD ! FROM PROGRESS
  USE BML
#endif

#ifdef MPI_ON
  USE MPI
#endif
#ifdef DBCSR_ON
  USE DBCSR_VAR_MOD
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: LATTE

CONTAINS

  !! Latte subroutine
  !! \param FLAGS Different control flags that can be passed to LATTE (not in use yet)
  !! \param NATS Number of atoms
  !! \param COORDS Coordinates. Example: y-coordinate of atom 1 = COORDS(2,1)
  !! \param TYPES An index for all the different atoms in the system.
  !! \param NTYPES Number of different elements in the system
  !! \param MASSES Element masses for every different element of the system.
  !! \param XLO Lowest dimensions of the box
  !! \param XHI Highest dimensions of the box
  !! \param FORCES Forces for every atom as output.
  !! \param MAXITER Latte MAXITER keyword. If MAXITER = -1, only the Forces are computed.
  !!        If MAXITER = 0, MAXITER is read from latte.in file.
  !!        IF MAXITER > 0, MAXITER is passed trough the library call.
  !! \param VENERG This is the potential Energy that is given back from latte to the hosting code.
  !! \param VEL Velocities passed to latte.
  !! \param DT integration step passed to latte.
  !!
  !! \brief This routine will be used load call latte_lib from a C/C++ program:
  !!
  !! \brief Note: To get the mass of atom 3 we do:
  !! \verbatim
  !!      MASS(TYPES(3))
  !! \endverbatim
  !!
  !! \brief Note: To get the lattice vectors as formated in LATTE we do:
  !! \verbatim
  !!      BOX(1,1) = XHI(1) - XLO(1); ...
  !! \endverbatim
  !!
  !! \brief Note: All units are LATTE units by default. See https://github.com/losalamos/LATTE/blob/master/Manual/LATTE_manual.pdf
  !!
  SUBROUTINE LATTE(NTYPES,TYPES,CR_IN,MASSES_IN,XLO,XHI,FTOT_OUT, &
                   MAXITER_IN, VENERG, VEL_IN, DT_IN)

    IMPLICIT NONE

    INTEGER :: I, J
    INTEGER :: START_CLOCK, STOP_CLOCK, CLOCK_RATE, CLOCK_MAX
    INTEGER :: MYID = 0
    REAL :: TARRAY(2), RESULT, SYSTDIAG, SYSTPURE
    CHARACTER(LEN=50) :: FLNM
    LOGICAL :: EXISTS

    REAL(LATTEPREC), INTENT(IN) :: CR_IN(:,:),VEL_IN(:,:), MASSES_IN(:),XLO(3),XHI(3)
    REAL(LATTEPREC), INTENT(IN) :: DT_IN
    REAL(LATTEPREC), INTENT(OUT) :: FTOT_OUT(:,:), VENERG
    INTEGER, INTENT(IN) ::  NTYPES, TYPES(:), MAXITER_IN

#ifdef PROGRESSON
    TYPE(SYSTEM_TYPE) :: SY
#endif

#ifdef MPI_ON
    INTEGER :: IERR, STATUS(MPI_STATUS_SIZE), NUMPROCS

    CALL MPI_INIT( IERR )
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, MYID, IERR )
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, NUMPROCS, IERR )
#endif

    IF(.NOT. INITIALIZED)THEN
      LIBCALLS = 0 ; MAXITER = -10
    ELSE
      LIBCALLS = LIBCALLS + 1
    ENDIF

    OPEN(UNIT=6, FILE="log.latte", FORM="formatted")

    IF(.NOT. INITIALIZED)THEN

      WRITE(6,*)"The log file for latte_lib"
      WRITE(6,*)""

      NUMSCF = 0
      CHEMPOT = ZERO

      ! Start timers
      TX = INIT_TIMER()
      TX = START_TIMER(LATTE_TIMER)

      INQUIRE( FILE="latte.in", exist=EXISTS )
      IF (EXISTS) THEN
        IF(.NOT. INITIALIZED) CALL PARSE_CONTROL("latte.in")
      ELSE
        IF(.NOT. INITIALIZED) CALL READCONTROLS
      ENDIF

      CALL READTB

      IF (RESTART .EQ. 0) THEN

        BOX = 0.0d0
        BOX(1,1) = xhi(1) - xlo(1)
        BOX(2,2) = xhi(2) - xlo(2)
        BOX(3,3) = xhi(3) - xlo(3)

        NATS = SIZE(CR_IN,DIM=2)

        IF(.NOT.ALLOCATED(CR)) ALLOCATE(CR(3,NATS))
        CR = CR_IN

        !     IF(.NOT.INITIALIZED)then
        ALLOCATE(ATELE(NATS))
        CALL MASSES2SYMBOLS(TYPES,NTYPES,MASSES_IN,NATS,ATELE)
        !     ENDIF

        CALL READCR

      ELSE
        CALL READRESTART
      ENDIF

      IF (PPOTON .EQ. 1) CALL READPPOT
      IF (PPOTON .EQ. 2) CALL READPPOTTAB


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

    ELSE

      BOX = 0.0d0
      BOX(1,1) = xhi(1) - xlo(1)
      BOX(2,2) = xhi(2) - xlo(2)
      BOX(3,3) = xhi(3) - xlo(3)

      NATS = SIZE(CR_IN,DIM=2)

      IF(.NOT.ALLOCATED(CR)) ALLOCATE(CR(3,NATS))
      CR = CR_IN

    ENDIF

    !END OF INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef PROGRESSON
    ! SY%NATS = NATS
    ! IF(.NOT. ALLOCATED(SY%COORDINATE))ALLOCATE(SY%COORDINATE(3,NATS))
    ! SY%COORDINATE = CR
    ! SY%SYMBOL = ATELE
    ! SY%lattice_vector = BOX
    ! call write_trajectory(sy,LIBCALLS,1,0.01d0,"traj","pdb")
    ! call write_system(sy,"sy","pdb")
#endif


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

      ! OPEN(UNIT=31, STATUS="UNKNOWN", FILE="myrho.dat")

      !         DO I = 1, HDIM
      !           DO J = 1,HDIM

      !              IF (ABS(BO(J,I)) .GT. 1.0D-5) WRITE(31,99) I, J

      !                  ENDDO
      !                  ENDDO

      !99 FORMAT(2I9)

      !     CLose(31)
      IF (DEBUGON .EQ. 1 .AND. SPINON .EQ. 0 .AND. KON .EQ. 0) THEN

        PRINT*, "Caution - you're writing to file the density matrix!"

        OPEN(UNIT=31, STATUS="UNKNOWN", FILE="myrho.dat")

        DO I = 1, HDIM
          WRITE(31,10) (BO(I,J), J = 1, HDIM)
        ENDDO

        CLOSE(31)

        10      FORMAT(100G18.8)

      ENDIF

      FTOT = ZERO

      IF (COMPFORCE .EQ. 1) THEN

        IF (KON .EQ. 0) THEN

          IF (SPONLY .EQ. 0) THEN
            CALL GRADHSP
          ELSE
            CALL GRADH
          ENDIF

        ELSE
          CALL KGRADH
        ENDIF

        FTOT = TWO * F

      ENDIF

      EREP = ZERO
      IF (PPOTON .EQ. 1) THEN
        CALL PAIRPOT
        FTOT = FTOT + FPP
      ENDIF

      IF (PPOTON .EQ. 2) THEN
        CALL PAIRPOTTAB
        FTOT = FTOT + FPP
      ENDIF

      IF (ELECTRO .EQ. 1) FTOT = FTOT + FCOUL

      IF (BASISTYPE .EQ. "NONORTHO") THEN

        IF (SPONLY .EQ. 0) THEN
          ! s/sp orbitals only so we can use the analytic code
          CALL FCOULNONO_SP
          CALL PULAY_SP
          IF (SPINON .EQ. 1) CALL FSPINNONO_SP
        ELSE
          ! Otherwise use the complex but general expansions Josh
          ! Coe implemented
          CALL FCOULNONO
          CALL PULAY
          IF (SPINON .EQ. 1) CALL FSPINNONO
        ENDIF

        FTOT = FTOT - TWO*FPUL + FSCOUL

        IF (SPINON .EQ. 1) FTOT = FTOT + FSSPIN

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
          CALL SUMMARY
          CALL FITTINGOUTPUT(0)

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

    ELSEIF (MDON .EQ. 1 .AND. RELAXME .EQ. 0 .AND. MAXITER_IN < 0 ) THEN

      DT = DT_IN ! Get the integration step from the hosting code.

      IF(LIBCALLS == 0)THEN

        IF (BASISTYPE .EQ. "NONORTHO") CALL ALLOCATENONO

        IF (XBOON .EQ. 1) CALL ALLOCATEXBO

        IF (ELECTRO .EQ. 1) THEN
          CALL ALLOCATECOULOMB
          CALL INITCOULOMB
        ENDIF

        ! Start the timers

        CALL SYSTEM_CLOCK(START_CLOCK, CLOCK_RATE, CLOCK_MAX)
        CALL DTIME(TARRAY, RESULT)

        Call SETUPTBMD

      ENDIF

      IF (MOD(LIBCALLS, UDNEIGH) .EQ. 0) CALL NEBLISTS(1)

      IF (QITER .NE. 0) THEN
        ECOUL = ZERO
        IF (ELECTRO .EQ. 1) CALL GETCOULE
      ENDIF

      write(*,*)"LIBCALLS",LIBCALLS

      IF(LIBCALLS > 0) CALL GETMDF(1, LIBCALLS)

      CALL TOTENG

      ! For the 0 SCF MD the coulomb energy is calculated in GETMDF

      IF (PPOTON .EQ. 1) THEN
        CALL PAIRPOT
      ELSEIF (PPOTON .EQ. 2) THEN
        CALL PAIRPOTTAB
      ENDIF

      IF (QITER .NE. 0) THEN
        ECOUL = ZERO
        IF (ELECTRO .EQ. 1) CALL GETCOULE
      ENDIF

      ESPIN = ZERO
      IF (SPINON .EQ. 1) CALL GETSPINE

      !CALL GETPRESSURE

      IF (CONTROL .NE. 1 .AND. CONTROL .NE. 2 .AND. KBT .GT. 0.000001 ) THEN

        ! Only required when using the recursive expansion of the Fermi operator

        ! 2/26/13
        ! The entropy is now calculated when we get the density
        ! matrix in the spin polarized case with diagonalization,
        ! as it should be...

        CALL ENTROPY

      ENDIF

      VENERG = TRRHOH + EREP - ENTE - ECOUL + ESPIN

      write(6,*)"Energy Components (TRRHOH, EREP, ENTE, ECOUL)",TRRHOH, EREP, ENTE, ECOUL
      write(6,*)"Epot", VENERG

      FTOT_OUT = FTOT


      INITIALIZED = .true.

      RETURN


  ELSEIF (MDON .EQ. 1 .AND. RELAXME .EQ. 0 .AND. MAXITER_IN >= 0) THEN


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

      CALL SYSTEM_CLOCK(START_CLOCK, CLOCK_RATE, CLOCK_MAX)

      CALL DOSFIT

      CALL SYSTEM_CLOCK(STOP_CLOCK, CLOCK_RATE, CLOCK_MAX)

      WRITE(6,'("# Wall time = ", F12.2, " s")') &
        FLOAT(STOP_CLOCK - START_CLOCK)/FLOAT(CLOCK_RATE)

    ELSEIF  (MDON .EQ. 0 .AND. RELAXME .EQ. 0 .AND. DOSFITON .EQ. 2) THEN

      CALL MOFIT

    ELSEIF (MDON .EQ. 0 .AND. RELAXME .EQ. 0 .AND. DOSFITON .EQ. 3) THEN

      CALL MOFITPLATO

    ELSEIF (MDON .EQ. 0 .AND. RELAXME .EQ. 0 .AND. PPFITON .EQ. 1) THEN

      CALL PPFIT

    ELSEIF (MDON .EQ. 0 .AND. RELAXME .EQ. 0 .AND. ALLFITON .EQ. 1) THEN

      CALL ALLFIT

  ELSE

     WRITE(6,*) "You can't have RELAXME = 1 and MDON = 1"
     STOP

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
  TX = SHUTDOWN_TIMER()

#ifdef MPI_ON
  CALL MPI_FINALIZE( IERR )
#endif

    INITIALIZED = .true.

  END SUBROUTINE LATTE

END MODULE LATTE_LIB
