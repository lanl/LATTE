!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

SUBROUTINE TBMD

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE PPOTARRAY
  USE MDARRAY
  USE NEBLISTARRAY
  USE COULOMBARRAY
  USE SPINARRAY
  USE VIRIALARRAY
  USE NONOARRAY
  USE MYPRECISION
  USE LATTEPARSER_LATTE_MOD

  IMPLICIT NONE

  INTEGER :: I
  INTEGER :: ITER
  INTEGER :: CURRITER, TOTSCF
  INTEGER :: START_CLOCK, STOP_CLOCK, CLOCK_RATE, CLOCK_MAX
  REAL(LATTEPREC) :: THETIME, NEWESPIN, NEWECOUL
  REAL(LATTEPREC) :: RN, MYVOL
  INTEGER :: FLAGAND

  IF (EXISTERROR) RETURN

  !
  ! Read MDcontroller to determine what kind of MD simulation to do
  !
  IF (LATTEINEXISTS) THEN
     CALL PARSE_MD("latte.in")
  ELSE
     CALL READMDCONTROLLER
  ENDIF
  !
  ! Allocate stuff for building the neighbor lists, then build them
  !

  CALL ALLOCATENEBARRAYS

  CALL NEBLISTS(0)

  !
  ! Allocate things depending on which method we're using
  ! to get the bond-order
  !

  IF (CONTROL .EQ. 1) THEN
     CALL ALLOCATEDIAG
  ELSEIF (CONTROL .EQ. 2 .OR. CONTROL .EQ. 4 .OR. CONTROL .EQ. 5) THEN
     CALL ALLOCATEPURE
  ELSEIF (CONTROL .EQ. 3) THEN
     CALL FERMIALLOCATE
  ENDIF

  !
  ! Set the seed for the randon number generator
  ! Applied to initializing the velocities and the
  ! Langevin and Andersen thermostats.
  !

  IF (SEEDINIT .EQ. "RANDOM") CALL INITRNG


  IF (RESTART .EQ. 0) THEN

     ALLOCATE (V(3,NATS))

     !
     ! Initialize velocities if TOINITTEMP = 1
     !

     IF (TOINITTEMP .EQ. 1) THEN
        CALL INITIALV
     ELSE
        WRITE(6,*) "Caution: you haven't initialized velocities"
     ENDIF

  ENDIF

  !
  ! If we're going to run with the Hugoniostat, some things need
  ! to be initialized
  !

  IF (SHOCKON .EQ. 1) CALL INITSHOCKCOMP


  !
  ! Get forces - we need these at this point only if we're running
  ! NVE MD with the velocity verlet algorithm
  !

  CURRITER = 0

  IF (RESTART .EQ. 0) THEN

     CALL GETMDF(0,1)

  ELSEIF (RESTART .EQ. 1) THEN

     !
     ! If we've read from a restart file then we don't need to run
     ! qconsistency to full self-consistency at the first time step of
     ! the new run - we're already there...
     !

     ! Yes we do until this is fixed

     CALL GETMDF(0,1)

  ENDIF

  IF (RESTART .EQ. 0) THEN
     ITER = -1
  ELSE
     ITER = CONTITER - 1
  ENDIF

  IF (GETHUG .EQ. 1 .AND. NVTON .NE. 0) THEN
     CALL ERRORS("tbmd","Please don't have GETHUG = 1 and NVTON .NE. 0")
  ENDIF

  IF (GETHUG .EQ. 1 .AND. NPTON .NE. 0) THEN
     CALL ERRORS("tbmd","Please don't have GETHUG = 1 and NPTON .NE. 0")
  ENDIF

  IF (NVTON .EQ. 1 .OR. NPTON .EQ. 1) THEN

     ALLOCATE(THIST(AVEPER))

     THIST = TTARGET

  ENDIF

  IF (NPTON .EQ. 1) THEN

     CALL GETPRESSURE

     IF (NPTTYPE .EQ. "ISO") THEN

        ! Change box dims assuming a hydrostatic pressure

        ALLOCATE( PHIST(AVEPER) )
        PHIST = PRESSURE

     ELSE ! Allow each box vector to change independently

        ALLOCATE (PHISTX(AVEPER), PHISTY(AVEPER), PHISTZ(AVEPER))

        PHISTX = PRESSURE
        PHISTY = PRESSURE
        PHISTZ = PRESSURE

     ENDIF

     CALL GETDENSITY

  ENDIF

  IF (GETHUG .EQ. 1) THEN

     ! Allocate the arrays for the averages

     CALL GETPRESSURE
     CALL GETKE

     MYVOL = ABS(BOX(1,1)*(BOX(2,2)*BOX(3,3) - BOX(3,2)*BOX(2,3)) + &
          BOX(1,2)*(BOX(2,1)*BOX(3,3) - BOX(3,1)*BOX(2,3)) + &
          BOX(1,3)*(BOX(2,1)*BOX(3,2) - BOX(3,1)*BOX(2,2)))

     ALLOCATE(PHIST(AVEPER/WRTFREQ), THIST(AVEPER/WRTFREQ), EHIST(AVEPER/WRTFREQ))
     ALLOCATE(VHIST(AVEPER/WRTFREQ))

     PHIST = PRESSURE
     THIST = TEMPERATURE
     EHIST = E0
     VHIST = MYVOL

  ENDIF

  IF (PARREP .EQ. 1) CALL PARAFILEOPEN


  CUMDT = ZERO

  TOTSCF = 0

  !  OPEN(UNIT=30, STATUS="UNKNOWN", FILE="AMD_gap_Ef.dat")

  WRITE(6,17) "#","Time (ps)", "Free energy (eV)", "T (K)", "Pressure (GPa)"

17 FORMAT(A1, 2X, A10, 6X, A16, 2X, A5, 3X, A14)

  !  CALL SYSTEM_CLOCK(START_CLOCK, CLOCK_RATE, CLOCK_MAX)

  CURRITER = 0

  DO WHILE (CURRITER .LE. MAXITER)

     TOTSCF = TOTSCF + SCFS_II

     CURRITER = CURRITER + 1

     ITER = ITER + 1
     ENTROPYITER = ITER

     IF (SHOCKON .EQ. 1 .AND. ITER .GE. SHOCKSTART .AND. &
          ITER .LT. SHOCKSTOP) THEN

        CALL SHOCKCOMP

        !
        ! Since we're chaning the dimensions of the box, we
        ! should probably also adjust the reciprocal lattice vectors
        ! used in the Ewald sum. There's no harm in doing this
        ! every time step while the box size is changing
        !
        ! We may as well reinitialize the Coulomb stuff since
        ! it also optimizes the real-space cut-off based on the volume
        ! of the cell
        !

        CALL INITCOULOMB

     ENDIF

     IF ((NVTON .EQ. 0 .AND. NPTON .EQ. 0 .AND. GETHUG .EQ. 0) .OR. CURRITER .GT. THERMRUN) THEN

        CALL VELVERLET(CURRITER)

     ELSEIF (NVTON .EQ. 1 .AND. CURRITER .LE. THERMRUN) THEN

        ! Velocity rescaling thermostat

        !
        ! To smooth things out, let's average the
        ! temperature over the previous AVEPER time steps
        !

        CALL AVETEMP

        IF (MOD(ITER, THERMPER) .EQ. 0 .AND. CURRITER .GT. 1) THEN

           CALL NVTRESCALE

        ELSE

           CALL VELVERLET(CURRITER)

        ENDIF

     ELSEIF (NVTON .EQ. 2 .AND. CURRITER .LE. THERMRUN) THEN

        ! Langevin thermostat

        CALL NVTLANGEVIN(ITER)

     ELSEIF (NVTON .EQ. 3 .AND. CURRITER .LE. THERMRUN) THEN

        ! Andersen thermostat


        CUMDT = CUMDT + DT

        CALL VELVERLET(CURRITER)

        CALL NVTANDERSEN

     ELSEIF (NVTON .EQ. 4 .AND. CURRITER .LE. THERMRUN) THEN

        ! NOSE thermostat

        CALL NVTNH

     ELSEIF (NPTON .EQ. 1 .AND. CURRITER .LE. THERMRUN) THEN

        ! Velocity rescaling thermostat

        !
        ! To smooth things out, let's average the
        ! temperature over the previous AVEPER time steps
        !

        CALL NVTLANGEVIN(ITER)

        CALL AVEPRESS

        IF (MOD(ITER, THERMPER) .EQ. 0 .AND. CURRITER .GT. 1) THEN

           CALL NPTRESCALE
           CALL GETDENSITY
           CALL INITCOULOMB

        ENDIF


     ELSEIF (GETHUG .EQ. 1 .AND. CURRITER .LE. THERMRUN) THEN

        CALL NVTLANGEVIN(ITER)

        IF (MOD(ITER, THERMPER) .EQ. 0 .AND. CURRITER .GT. 1) THEN

           CALL HUGRESCALE
           CALL GETDENSITY
           CALL INITCOULOMB

        ENDIF

     ENDIF

     IF (MOD(ITER,WRTFREQ) .EQ. 0) THEN

        ! The non-orthogonal density matrix is obtained from
        ! GETMDF (when the Pulay force is computed)

        IF (PPOTON .EQ. 1) THEN
           CALL PAIRPOT
        ELSEIF (PPOTON .EQ. 2) THEN
           CALL PAIRPOTTAB
        ELSEIF (PPOTON .EQ. 3) THEN
           CALL PAIRPOTSPLINE
        ENDIF

        CALL GETKE

        CALL TOTENG

        ! For the 0 SCF MD the coulomb energy is calculated in GETMDF

        IF (QITER .NE. 0) THEN
           ECOUL = ZERO
           IF (ELECTRO .EQ. 1) CALL GETCOULE
        ENDIF

        ESPIN = ZERO
        IF (SPINON .EQ. 1) CALL GETSPINE

        CALL GETPRESSURE

        IF (CONTROL .NE. 1 .AND. CONTROL .NE. 2 .AND. KBT .GT. 0.000001 ) THEN

           ! Only required when using the recursive expansion of the Fermi operator

           ! 2/26/13
           ! The entropy is now calculated when we get the density
           ! matrix in the spin polarized case with diagonalization,
           ! as it should be...

           CALL ENTROPY

        ENDIF

        TOTE = TRRHOH + EREP + KEE - ENTE - ECOUL + ESPIN

        !         write(*,*)"Ekin", KEE
        !         write(*,*)"Epot", TRRHOH + EREP - ENTE - ECOUL + ESPIN
        !         write(*,*)"components",TRRHOH, EREP, ENTE, ECOUL

        IF (GETHUG .EQ. 1) CALL AVESFORHUG(PRESSURE, TOTE, TEMPERATURE, SYSVOL)

        IF (ABS(TOTE) .GT. 10000000000.0) THEN
           CALL ERRORS("tbmd","The calculation has diverged - check SCF parameters")
        ENDIF

        !        write(70,*) ITER, CR(1,1)

        THETIME = REAL(ITER)*DT/THOUSAND

        IF (PARREP .EQ. 0) THEN

           IF (NPTON .EQ. 0 .AND. GETHUG .EQ. 0 .AND. NVTON .NE. 0) THEN

              IF (NVTON .NE. 4) THEN

                 WRITE(6,99)"Data", THETIME, TOTE, TEMPERATURE, PRESSURE, EGAP, &
                      CHEMPOT !, TRRHOH, EREP, KEE, ECOUL, REAL(NUMSCF)

              ELSE

                 ! Special case for Nose Hoover
                 WRITE(6,99)"Data", THETIME, TOTE, TOTE+CONSMOT, TEMPERATURE, PRESSURE, &
                      EGAP, CHEMPOT !, STRTEN(1), STRTEN(2), STRTEN(3), STRTEN(4), &
                 !STRTEN(5), STRTEN(6)

              ENDIF

           ENDIF

           IF (NVTON .EQ. 0 .AND. NPTON .EQ. 0 .AND. GETHUG .EQ. 0) THEN

              WRITE(6,99)"Data", THETIME, TOTE, TEMPERATURE, PRESSURE, EGAP, &
                   CHEMPOT

           ENDIF

           IF (NPTON .NE. 0 .AND. NVTON .EQ. 0 .AND. GETHUG .EQ. 0) THEN

              WRITE(6,99)"Data", THETIME, TOTE, TEMPERATURE, PRESSURE, &
                   EGAP, MASSDEN, BOX(1,1), BOX(2,2), BOX(3,3), &
                   SYSVOL

           ENDIF

           IF (GETHUG .EQ. 1 .AND. NVTON .EQ. 0 .AND. NPTON .EQ. 0) THEN

              WRITE(6,99)"Data", THETIME, TOTE, TEMPERATURE, PRESSURE, &
                   EGAP, MASSDEN, HG, TTARGET, BOX(1,1), BOX(2,2), BOX(3,3), &
                   SYSVOL

           ENDIF

           FLUSH(6)

        ELSE

           CALL PARAWRITE(THETIME)

        ENDIF


99      FORMAT(A4,20G18.9)

16      FORMAT(F12.5, F20.8, 1X, F9.1, 1X, F12.3, 1X, G18.9, 1X, G18.9)

     ENDIF

     IF (MOD(ITER, DUMPFREQ) .EQ. 0) CALL WRTCFGS(ITER)

     IF (MOD(ITER, RSFREQ) .EQ. 0) CALL WRTRESTART(ITER)

     IF (MOD(ITER, UDNEIGH) .EQ. 0) CALL NEBLISTS(1)


  ENDDO

  !  CALL SYSTEM_CLOCK(STOP_CLOCK, CLOCK_RATE, CLOCK_MAX)
  !  PRINT*, HDIM, FLOAT(STOP_CLOCK - START_CLOCK)/FLOAT(MAXITER*CLOCK_RATE)

  IF (CONTROL .EQ. 1) THEN
     !     CALL DEALLOCATEDIAG
  ELSEIF (CONTROL .EQ. 2 .OR. CONTROL .EQ. 4 .OR. CONTROL .EQ. 5) THEN
     CALL DEALLOCATEPURE
  ELSEIF (CONTROL .EQ. 3) THEN
     CALL FERMIDEALLOCATE
  ENDIF

  IF (NVTON .EQ. 1) DEALLOCATE(THIST)
  IF (NPTON .EQ. 1) THEN
     IF (NPTTYPE .EQ. "ISO") THEN
        DEALLOCATE(THIST, PHIST)
     ELSE
        DEALLOCATE(THIST, PHISTX, PHISTY, PHISTZ)
     ENDIF
  ENDIF

  RETURN

END SUBROUTINE TBMD
