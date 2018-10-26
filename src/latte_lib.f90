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
  USE CONSTRAINTS_MOD

#ifdef PROGRESSON
  USE PRG_SYSTEM_MOD ! FROM PROGRESS
  USE PRG_PULAYMIXER_MOD
  USE PRG_EXTRAS_MOD
  USE MIXER_MOD
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

  ! Defines the version of the binary interface.
  ! Adjust if non-backward compatible changes are made to the interface.
  ! Use in codes using the library interface to make certain a compatible
  ! version of the LATTE library is used.
  INTEGER, PARAMETER :: LATTE_ABIVERSION = 20180622

  PUBLIC :: LATTE, LATTE_ABIVERSION

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
  !! \param XY Tilt factor.
  !! \param XZ Tilt factor.
  !! \param YZ Tilt factor. The lattice vectors are constructed as:
  !! a = (xhi-xlo,0,0); b = (xy,yhi-ylo,0); c = (xz,yz,zhi-zlo).
  !! \param FORCES Forces for every atom as output.
  !! \param MAXITER Latte MAXITER keyword. If MAXITER = -1, only the Forces are computed.
  !!        If MAXITER = 0, MAXITER is read from latte.in file.
  !!        IF MAXITER > 0, MAXITER is passed trough the library call.
  !! \param VENERG This is the potential Energy that is given back from latte to the hosting code.
  !! \param VEL Velocities passed to latte.
  !! \param DT integration step passed to latte.
  !! \param VIRIAL_INOUT Components of the second virial coefficient.
  !! \param NEWSYSTEM Tells LATTE if a new system is passed.
  !! \param EXISTERROR_INOUT Returns an error flag (.true.) to the hosting code.
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
  SUBROUTINE LATTE(NTYPES, TYPES, CR_IN, MASSES_IN, XLO, XHI, XY, XZ, YZ, FTOT_OUT, &
       MAXITER_IN, VENERG, VEL_IN, DT_IN, VIRIAL_INOUT, NEWSYSTEM, EXISTERROR_INOUT)

    USE CONSTANTS_MOD, ONLY: EXISTERROR

    IMPLICIT NONE

    INTEGER :: I, J
    INTEGER :: START_CLOCK, STOP_CLOCK, CLOCK_RATE, CLOCK_MAX
    INTEGER :: MYID = 0
    REAL :: TARRAY(2), RESULT, SYSTDIAG, SYSTPURE
    REAL(LATTEPREC) :: DBOX
    CHARACTER(LEN=50) :: FLNM

    REAL(LATTEPREC), INTENT(IN)  :: CR_IN(:,:),VEL_IN(:,:), MASSES_IN(:),XLO(3),XHI(3)
    REAL(LATTEPREC), INTENT(IN)  :: DT_IN, XY, XZ, YZ
    REAL(LATTEPREC), INTENT(OUT) :: FTOT_OUT(:,:), VENERG
    REAL(LATTEPREC), INTENT(OUT) :: VIRIAL_INOUT(6)
    INTEGER, INTENT(IN) ::  NTYPES, TYPES(:), MAXITER_IN
    LOGICAL(1), INTENT(INOUT) :: EXISTERROR_INOUT
    REAL(LATTEPREC) :: MLSI, LUMO, HOMO
    INTEGER, INTENT(INOUT) :: NEWSYSTEM

#ifdef PROGRESSON
    TYPE(SYSTEM_TYPE) :: SY
#endif

#ifdef MPI_ON
    INTEGER :: IERR, STATUS(MPI_STATUS_SIZE), NUMPROCS
    CALL MPI_INIT( IERR )
    CALL MPI_COMM_RANK( MPI_COMM_WORLD, MYID, IERR )
    CALL MPI_COMM_SIZE( MPI_COMM_WORLD, NUMPROCS, IERR )
#endif

    EXISTERROR = .FALSE. !We assume we start the lib call without errors

    IF(.NOT. LIBINIT .OR. NEWSYSTEM == 1)THEN

       CALL DEALLOCATEALL()

       LIBRUN = .TRUE.

       LIBCALLS = 0 ; MAXITER = -10

      ! Only LATTE main code will create the animate folder
      !  INQUIRE( FILE="animate/.", exist=LATTEINEXISTS)
      !  IF (.NOT. LATTEINEXISTS) CALL SYSTEM("mkdir animate")

       NUMSCF = 0
       CHEMPOT = ZERO

       ! Start timers
       TX = INIT_TIMER()
       TX = START_TIMER(LATTE_TIMER)

       INQUIRE( FILE="latte.in", exist=LATTEINEXISTS )

       IF (LATTEINEXISTS) THEN
          CALL PARSE_CONTROL("latte.in")

#ifdef PROGRESSON
          CALL PRG_PARSE_MIXER(MX,"latte.in")
#endif

       ELSE
          CALL READCONTROLS
       ENDIF

       IF(VERBOSE <= 0)THEN
          OPEN(UNIT=6, FILE="/dev/null", FORM="formatted")
       ELSE
          OPEN(UNIT=6, FILE="log.latte", FORM="formatted")
       ENDIF

       IF(VERBOSE >= 1)THEN
          WRITE(*,*)"# The log file for latte_lib"
          WRITE(*,*)""
          CALL TIMEDATE_TAG("LATTE started at : ")
       ENDIF

       CALL READTB

       IF (RESTART .EQ. 0) THEN

          BOX = 0.0d0
          BOX(1,1) = xhi(1) - xlo(1)
          BOX(2,1) = XY
          BOX(2,2) = xhi(2) - xlo(2)
          BOX(3,1) = XZ
          BOX(3,2) = YZ
          BOX(3,3) = xhi(3) - xlo(3)

          IF(VERBOSE >= 1)THEN
             WRITE(*,*)"Lattice vectors:"
             WRITE(*,*)"a=",BOX(1,1),BOX(1,2),BOX(1,3)
             WRITE(*,*)"b=",BOX(2,1),BOX(2,2),BOX(2,3)
             WRITE(*,*)"c=",BOX(3,1),BOX(3,2),BOX(3,3)
             WRITE(*,*)""
          ENDIF

          BOX_OLD = BOX

          NATS = SIZE(CR_IN,DIM=2)

          IF (.NOT.ALLOCATED(CR)) ALLOCATE(CR(3,NATS))
          CR = CR_IN

          IF(VERBOSE >= 1)WRITE(*,*)"Converting masses to symbols ..."
          IF(.NOT. ALLOCATED(ATELE)) ALLOCATE(ATELE(NATS))
          CALL MASSES2SYMBOLS(TYPES,NTYPES,MASSES_IN,NATS,ATELE)

          !Forces, charges and element pointers are allocated in readcr
          CALL READCR

          CALL FLUSH(6)

       ELSE

          IF(VERBOSE >= 1)WRITE(*,*)"Restarting calculation from file ..."
          CALL READRESTART

       ENDIF

       IF (VERBOSE >= 1) WRITE(*,*)"Reading ppots from file (if PPOTON >= 1) ..."
       IF (PPOTON .EQ. 1) CALL READPPOT
       IF (PPOTON .EQ. 2) CALL READPPOTTAB
       IF (PPOTON .EQ. 3) CALL READPPOTSPLINE

       IF (DEBUGON .EQ. 1) THEN
          CALL PLOTUNIV
          IF (PPOTON .EQ. 1) CALL PLOTPPOT
       ENDIF

       CALL GETHDIM

       CALL GETMATINDLIST

       IF (VERBOSE >= 1) WRITE(*,*)"Getting rho0 ..."
       CALL RHOZERO

       CALL GETBNDFIL()

       CALL FLUSH(6)

#ifdef GPUON

       CALL INITIALIZE( NGPU )

#endif

#ifdef DBCSR_ON

       IF (CONTROL .EQ. 2 .AND. SPARSEON .EQ. 1) CALL INIT_DBCSR

#endif

       IF (KBT .LT. 0.0000001 .OR. CONTROL .EQ. 2) ENTE = ZERO

       IF (.NOT. ALLOCATED(V)) ALLOCATE(V(3,NATS))
       IF(VERBOSE >= 1)WRITE(*,*)"End of INITIALIZATION"

    ELSE

       BOX = 0.0d0
       BOX(1,1) = xhi(1) - xlo(1)
       BOX(2,1) = XY
       BOX(2,2) = xhi(2) - xlo(2)
       BOX(3,1) = XZ
       BOX(3,2) = YZ
       BOX(3,3) = xhi(3) - xlo(3)

       IF(VERBOSE >= 1)THEN
          WRITE(*,*)"Lattice vectors:"
          WRITE(*,*)"a=",BOX(1,1),BOX(1,2),BOX(1,3)
          WRITE(*,*)"b=",BOX(2,1),BOX(2,2),BOX(2,3)
          WRITE(*,*)"c=",BOX(3,1),BOX(3,2),BOX(3,3)
          WRITE(*,*)""
       ENDIF

       LIBCALLS = LIBCALLS + 1

       NATS = SIZE(CR_IN,DIM=2)

       IF(.NOT.ALLOCATED(CR)) ALLOCATE(CR(3,NATS))
       CR = CR_IN

       CALL FLUSH(6)

    ENDIF
    !End of initialization


    IF (MDON .EQ. 0 .AND. RELAXME .EQ. 0 .AND. DOSFITON .EQ. 0 &
         .AND. PPFITON .EQ. 0 .AND. ALLFITON .EQ. 0) THEN


       ! IF (.NOT. LIBINIT) THEN

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

       !  ELSE
       !     CALL ERRORS("latte_lib","Attemting to perform multiple single point calculations. &
       !          & MDON .EQ. 0 .AND. RELAXME .EQ. 0 can be done &
       !          &for only one geometry (A single call to the library). Please reduce &
       !          &the number of steps (md of relaxation) at the host code")
       !     EXISTERROR_INOUT = EXISTERROR
       !     RETURN
       !  ENDIF

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

       ! DO I = 1, HDIM
       !   DO J = 1,HDIM

       !     IF (ABS(BO(J,I)) .GT. 1.0D-5) WRITE(31,99) I, J

       !   ENDDO
       ! ENDDO

       !99 FORMAT(2I9)
       !CLOSE(31)

       IF (DEBUGON .EQ. 1 .AND. SPINON .EQ. 0 .AND. KON .EQ. 0) THEN

          PRINT*, "Caution - you're writing to file the density matrix!"

          OPEN(UNIT=31, STATUS="UNKNOWN", FILE="myrho.dat")

          DO I = 1, HDIM
             WRITE(31,10) (BO(I,J), J = 1, HDIM)
          ENDDO

          CLOSE(31)

10        FORMAT(100G18.8)

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

       IF (PPOTON .EQ. 3) THEN
          CALL PAIRPOTSPLINE
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
          !  CALL DEALLOCATEDIAG
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

       !  WRITE(*,*) "Force", FPP(1,1), FPP(2,1), FPP(3,1)
       !  PRINT*, "PCHECK ", (1.0/3.0)*(VIRBOND(1)+VIRBOND(2) + VIRBOND(3)), &
       !       (1.0/3.0)*(VIRCOUL(1)+VIRCOUL(2) + VIRCOUL(3)), &
       !       (1.0/3.0)*(VIRPAIR(1)+VIRPAIR(2) + VIRPAIR(3)), &
       !       (1.0/3.0)*(VIRPUL(1)+VIRPUL(2) + VIRPUL(3)), &
       !       (1.0/3.0)*(VIRSCOUL(1)+VIRSCOUL(2) + VIRSCOUL(3))

#ifdef DBCSR_ON

       IF (CONTROL .EQ. 2 .AND. SPARSEON .EQ. 1 .AND.  MYNODE .EQ. 0) THEN

#endif

          IF (MYID .EQ. 0) THEN
             CALL FITTINGOUTPUT(0)
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

       CALL DEALLOCATEALL

       LIBINIT = .FALSE.

       RETURN

    ELSEIF (MDON .EQ. 1 .AND. RELAXME .EQ. 0 .AND. MAXITER_IN < 0 ) THEN

       IF(VERBOSE >= 1)WRITE(*,*)"Insie MDON= 1 and RELAXME= 0 ..."

       DT = DT_IN ! Get the integration step from the hosting code.

       V = VEL_IN/1000.0d0  !Convert from Ang/ps to Ang/fs

       !Control for implicit geometry optimization.
       !This will need to be replaced by a proper flag.
       IF (DT_IN == 0) THEN
          IF (VERBOSE >= 1) WRITE(*,*)"NOTE: DT = 0 => FULLQCONV = 1"
          IF (VERBOSE >= 1) WRITE(*,*)"NOTE: DT = 0 => MDMIX = QMIX"
          FULLQCONV = 1
          MDMIX = QMIX
          CALL FLUSH(6)
       ENDIF

       IF (LIBCALLS == 0) THEN

          IF (VERBOSE >= 1)WRITE(*,*)"Allocating nonorthogonal arrays ..."
          IF (BASISTYPE .EQ. "NONORTHO") CALL ALLOCATENONO

          IF (VERBOSE >= 1)WRITE(*,*)"Allocating XLBO arrays ..."
          IF (XBOON .EQ. 1) CALL ALLOCATEXBO

          IF (VERBOSE >= 1)WRITE(*,*)"Allocating COULOMB arrays ..."
          IF (ELECTRO .EQ. 1) THEN
             CALL ALLOCATECOULOMB
             CALL INITCOULOMB
          ENDIF

          ! Start the timers
          IF (VERBOSE >= 1)WRITE(*,*)"Starting timers ..."
          CALL SYSTEM_CLOCK(START_CLOCK, CLOCK_RATE, CLOCK_MAX)
          CALL DTIME(TARRAY, RESULT)

          IF (VERBOSE >= 1)WRITE(*,*)"Setting up TBMD ..."
          CALL SETUPTBMD(NEWSYSTEM)

          CALL FLUSH(6)

       ELSEIF (LIBCALLS > 0 .AND. RESTARTLIB == 0) THEN

          DBOX = SQRT((BOX(1,1)-BOX_OLD(1,1))**2 + &
               & (BOX(2,2)-BOX_OLD(2,2))**2 + (BOX(3,3)-BOX_OLD(3,3))**2)

          ! Reinitializing Coulombic contribution if Kpoints are used
          IF ((ELECTRO .EQ. 1 .AND. ELECMETH .EQ. 0) .OR. DBOX .GT. 0.0d0) THEN
             CALL INITCOULOMB
          ENDIF

          IF (MOD(LIBCALLS, UDNEIGH) .EQ. 0) THEN
             !If box is changing
             IF (DBOX .GT. 0.0d0) THEN
                CALL NEBLISTS(0)
                CALL INITCOULOMB !If the box is changing we need to recompute the kspace list
             ELSE
                CALL NEBLISTS(1)
             ENDIF
          ENDIF

          BOX_OLD = BOX

       ENDIF

       IF (QITER .NE. 0) THEN
          ECOUL = ZERO
          IF (ELECTRO .EQ. 1) CALL GETCOULE
       ENDIF

       IF(VERBOSE >= 1) WRITE(*,*)"LIBCALLS",LIBCALLS

       IF(VERBOSE >= 1) MLSI = TIME_MLS()
       IF(LIBCALLS > 0) CALL GETMDF(1, LIBCALLS)
       IF(VERBOSE >= 1) WRITE(*,*)"Time for GETMDF =", TIME_MLS()-MLSI

       CALL TOTENG

       ! For the 0 SCF MD the coulomb energy is calculated in GETMDF

       IF (PPOTON .EQ. 1) THEN
          CALL PAIRPOT
       ELSEIF (PPOTON .EQ. 2) THEN
          CALL PAIRPOTTAB
       ELSEIF (PPOTON .EQ. 3) THEN
          CALL PAIRPOTSPLINE
       ENDIF

       IF (QITER .NE. 0) THEN
          ECOUL = ZERO
          IF (ELECTRO .EQ. 1) CALL GETCOULE
       ENDIF

       ESPIN = ZERO
       IF (SPINON .EQ. 1) CALL GETSPINE

       IF (CONTROL .NE. 1 .AND. CONTROL .NE. 2 .AND. KBT .GT. 0.000001 ) THEN

          ! Only required when using the recursive expansion of the Fermi operator

          ! 2/26/13
          ! The entropy is now calculated when we get the density
          ! matrix in the spin polarized case with diagonalization,
          ! as it should be...

          CALL ENTROPY

       ENDIF

       IF(VERBOSE >= 1) WRITE(*,*)"Energy Components (TRRHOH, EREP, ENTE, ECOUL)",TRRHOH, EREP, ENTE, ECOUL

       IF (FREEZE .EQ. 1) CALL FREEZE_ATOMS(FTOT,V)

       IF(MAXVAL(FTOT_OUT) .NE. 0.0d0)THEN
          IF(VERBOSE >= 1) WRITE(*,*)"Adding force components and energies from applicacion code ..."
          IF(VERBOSE >= 1) WRITE(*,*)"APPCODE,LATTE",VENERG,TRRHOH + EREP - ENTE - ECOUL + ESPIN
          VENERG = TRRHOH + EREP - ENTE - ECOUL + ESPIN
          FTOT_OUT = FTOT_OUT +  FTOT
       ELSE
          VENERG = TRRHOH + EREP - ENTE - ECOUL + ESPIN
          FTOT_OUT = FTOT
       ENDIF


       ! Get the seccond virial coefficient to pass it to the application program
       IF (ELECTRO .EQ. 0) VIRCOUL = ZERO
       VIRIAL = VIRBOND + VIRPAIR + VIRCOUL

       IF (SPINON .EQ. 1) VIRIAL = VIRIAL + VIRSSPIN

       IF (BASISTYPE .EQ. "NONORTHO") THEN
          VIRIAL = VIRIAL - VIRPUL + VIRSCOUL
       ENDIF

       !       CALL GETPRESSURE

       VIRIAL_INOUT = -VIRIAL

       LIBINIT = .TRUE.
       NEWSYSTEM = 0 !Setting newsystem back to 0.

#ifdef PROGRESSON
       IF(MOD(LIBCALLS,WRTFREQ) == 0)THEN
          IF(VERBOSE >= 1) THEN
             WRITE(*,*)"Writing trajectory into trajectory.pdb ..."
             SY%NATS = NATS
             IF(.NOT. ALLOCATED(SY%COORDINATE))ALLOCATE(SY%COORDINATE(3,NATS))
             SY%COORDINATE = CR
             SY%SYMBOL = ATELE
             SY%LATTICE_VECTOR = BOX
             SY%NET_CHARGE = DELTAQ
             CALL PRG_WRITE_TRAJECTORY(SY,LIBCALLS,WRTFREQ,DT_IN,"trajectory","pdb")

             WRITE(*,*)"Writing trajectory into trajectory.xyz ..."
             IF(LIBCALLS .EQ. 0)THEN
                OPEN(UNIT=20,FILE="trajectory.xyz",STATUS='unknown')
             ELSE
                OPEN(UNIT=20,FILE="trajectory.xyz",ACCESS='append',STATUS='old')
             ENDIF
             !Extended xyz file.
             WRITE(20,*)NATS
             WRITE(20,*) 'Lattice="',BOX(1,1),BOX(1,2),BOX(1,3),&
                  &BOX(2,1),BOX(2,2),BOX(2,3),BOX(3,1),BOX(3,2),BOX(3,3),'"',&
                  &"Properties=species:S:1:pos:R:3:vel:R:3:for:R:3:cha:R:1  Time=",LIBCALLS*DT_IN
             DO I=1,NATS
                WRITE(20,*)ATELE(I),CR(1,I),CR(2,I),CR(3,I),V(1,I),V(2,I),V(3,I),&
                     &FTOT(1,I),FTOT(2,I),FTOT(3,I),-DELTAQ(I)
             ENDDO
             CLOSE(20)
          ENDIF
       ENDIF
#else
       IF(MOD(LIBCALLS,WRTFREQ) == 0)THEN
          IF(VERBOSE >= 1) THEN
             WRITE(*,*)"Writing trajectory into trajectory.xyz ..."
             IF(LIBCALLS .EQ. 0)THEN
                OPEN(UNIT=20,FILE="trajectory.xyz",STATUS='unknown')
             ELSE
                OPEN(UNIT=20,FILE="trajectory.xyz",ACCESS='append',STATUS='old')
             ENDIF
             !Extended xyz file.
             WRITE(20,*)NATS
             WRITE(20,*) 'Lattice="',BOX(1,1),BOX(1,2),BOX(1,3),&
                  &BOX(2,1),BOX(2,2),BOX(2,3),BOX(3,1),BOX(3,2),BOX(3,3),'"',&
                  &"Properties=species:S:1:pos:R:3:vel:R:3:for:R:3:cha:R:1  Time=",LIBCALLS*DT_IN
             DO I=1,NATS
                WRITE(20,*)ATELE(I),CR(1,I),CR(2,I),CR(3,I),V(1,I),V(2,I),V(3,I),&
                     &FTOT(1,I),FTOT(2,I),FTOT(3,I),-DELTAQ(I)
             ENDDO
             CLOSE(20)
          ENDIF
       ENDIF
#endif

       IF(VERBOSE >= 1  .AND. CONTROL == 1)THEN
          IF(SPINON == 0) THEN
             HOMO = EVALS(FLOOR(BNDFIL*FLOAT(HDIM)))
             LUMO = EVALS(FLOOR(BNDFIL*FLOAT(HDIM))+1)
             WRITE(*,*)"HOMO=",HOMO, "LUMO=",LUMO
             WRITE(*,*)"EGAP=",LUMO - HOMO
          ELSE
             HOMO = MAX(DOWNEVALS(FLOOR(BNDFIL*FLOAT(HDIM))),UPEVALS(FLOOR(BNDFIL*FLOAT(HDIM))))
             LUMO = MIN(DOWNEVALS(FLOOR(BNDFIL*FLOAT(HDIM))+1),UPEVALS(FLOOR(BNDFIL*FLOAT(HDIM))+1))
             WRITE(*,*)"HOMO=",HOMO, "LUMO=",LUMO
             WRITE(*,*)"EGAP=",LUMO - HOMO
          ENDIF
       ENDIF

       IF (MOD(LIBCALLS, RSFREQ) .EQ. 0)THEN
          IF(VERBOSE >= 0) CALL WRTRESTARTLIB(LIBCALLS)
       ENDIF

       IF(RESTARTLIB == 1 .AND. LIBCALLS == 0)THEN
          CALL READRESTARTLIB(LIBCALLS)
       ENDIF

#ifdef PROGRESSON
       !  IF(VERBOSE >= 1)THEN
       !    CALL GETDIPOLE(DIPOLEMAG,DIPOLEVECOUT=DIPOLEVECOUT)
       !    WRITE(*,*)"Dipole Magnitude = DIPOLEMAG"
       !    WRITE(*,*)"Dipole Vector:"
       !    WRITE(*,*)DIPOLEVECOUT(1),DIPOLEVECOUT(2),DIPOLEVECOUT(3)
       !  ENDIF
#endif

       CALL FLUSH(6) !To force writing to file at every call

       EXISTERROR_INOUT = EXISTERROR

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

       CALL ERRORS("latte_lib","This option was not tested for the &
            &library version of LATTE: RELAXME= 1")
       RETURN

       CALL MSRELAX

    ELSEIF (MDON .EQ. 0 .AND. RELAXME .EQ. 0 .AND. DOSFITON .EQ. 1) THEN

       CALL ERRORS("latte_lib","This option was not tested for the &
            &library version of LATTE: DOSFITON= 1")
       RETURN

       CALL SYSTEM_CLOCK(START_CLOCK, CLOCK_RATE, CLOCK_MAX)

       CALL DOSFIT

       CALL SYSTEM_CLOCK(STOP_CLOCK, CLOCK_RATE, CLOCK_MAX)

       WRITE(6,'("# Wall time = ", F12.2, " s")') &
            FLOAT(STOP_CLOCK - START_CLOCK)/FLOAT(CLOCK_RATE)

    ELSEIF  (MDON .EQ. 0 .AND. RELAXME .EQ. 0 .AND. DOSFITON .EQ. 2) THEN

       CALL ERRORS("latte_lib","This option was not tested for the &
            &library version of LATTE: DOSFITON= 3")
       RETURN

       CALL MOFIT

    ELSEIF (MDON .EQ. 0 .AND. RELAXME .EQ. 0 .AND. DOSFITON .EQ. 3) THEN

       CALL ERRORS("latte_lib","This option was not tested for the &
            &library version of LATTE: DOSFITON= 3")
       RETURN

       CALL MOFITPLATO

    ELSEIF (MDON .EQ. 0 .AND. RELAXME .EQ. 0 .AND. PPFITON .EQ. 1) THEN

       CALL ERRORS("latte_lib","This option was not tested for the &
            &library version of LATTE: PPFITON= 1")
       RETURN

       CALL PPFIT

    ELSEIF (MDON .EQ. 0 .AND. RELAXME .EQ. 0 .AND. ALLFITON .EQ. 1) THEN

       CALL ERRORS("latte_lib","This option was not tested for the &
            &library version of LATTE: ALLFITON= 1")
       RETURN

       CALL ALLFIT

    ELSE

       CALL ERRORS("latte_lib","You can't have RELAXME = 1 and MDON = 1")

    ENDIF

#ifdef GPUON

    CALL SHUTDOWN()

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

    LIBINIT = .TRUE.
    NEWSYSTEM = 0 !Setting newsystem back to 0.
    EXISTERROR_INOUT = EXISTERROR

  END SUBROUTINE LATTE

END MODULE LATTE_LIB
