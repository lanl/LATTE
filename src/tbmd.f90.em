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

  IMPLICIT NONE

  INTEGER :: I
  INTEGER :: ITER
  INTEGER :: CURRITER, TOTSCF
  INTEGER :: START_CLOCK, STOP_CLOCK, CLOCK_RATE, CLOCK_MAX
  REAL(LATTEPREC) :: THETIME, NEWESPIN, NEWECOUL
  REAL(LATTEPREC) :: RN
  INTEGER :: FLAGAND

  !
  ! Read MDcontroller to determine what kind of MD simulation to do
  !

  INQUIRE( FILE="latte.in", exist=EXISTS )
    IF (EXISTS) THEN  
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
  
!!$  IF (SEEDINIT .EQ. "RANDOM") CALL INITRNG
  CALL INITRNG

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
  CONSMOT = ZERO

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

  IF (NVTON .GE. 1 .OR. NPTON .EQ. 1) THEN
     
     ALLOCATE(THIST(AVEPER))

     THIST = TTARGET

  ENDIF
  IF (NVTON .EQ. 2) THEN
     
     ALLOCATE(FRANPREV(3,NATS))
     FRANPREV = ZERO

  ENDIF
  IF (NVTON .EQ. 3) THEN

     CUMDT = ZERO
     QITERAND = 10
     QITERIN = QITER
     FLAGAND = ZERO

  ENDIF
  IF (NVTON .EQ. 4) THEN

     GAMMA = ZERO
     DGAMMA = ZERO
     
  ENDIF
  IF (NPTON .EQ. 1) THEN 
     ALLOCATE( PHIST(AVEPER) )
     PHIST = ZERO

     CALL GETDENSITY
  ENDIF

     
  TOTSCF = 0

  WRITE(6,17) "#","Time (ps)", "Free energy (eV)", "T (K)", "Pressure (GPa)"

17 FORMAT(A1, 2X, A10, 6X, A16, 2X, A5, 3X, A14) 

!  CALL SYSTEM_CLOCK(START_CLOCK, CLOCK_RATE, CLOCK_MAX)

  CURRITER = 0
  
  DO WHILE (ITER .LE. MAXITER)

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

     IF ((NVTON .EQ. 0 .AND. NPTON .EQ. 0) .OR. ITER .GT. THERMRUN) THEN

        CALL VELVERLET(CURRITER)

     ELSEIF (NVTON .EQ. 1 .AND. ITER .LE. THERMRUN) THEN
        
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

     ELSEIF (NVTON .EQ. 2 .AND. ITER .LE. THERMRUN) THEN
        
        ! Langevin thermostat

        !
        ! To smooth things out, let's average the
        ! temperature over the previous AVEPER time steps
        !

!!$        CALL AVETEMP
        
!!$        IF (CURRITER .GT. 1) THEN
        IF (MOD(ITER, THERMPER) .EQ. 0 .AND. CURRITER .GT. 0) THEN
          
           CALL NVTLANGEVIN(CURRITER)

        ELSE

           CALL VELVERLET(CURRITER)

        ENDIF

     ELSEIF (NVTON .EQ. 3 .AND. ITER .LE. THERMRUN) THEN
        
        ! Andersen thermostat

        !
        ! To smooth things out, let's average the
        ! temperature over the previous AVEPER time steps
        !

!!$        CALL AVETEMP

        CUMDT = CUMDT + DT
        IF(FLAGAND .EQ. ZERO) THEN
           CALL RANDOM_NUMBER(RN)
           TAU = -LOG(RN)/FRICTION
        ENDIF
        
!!$        IF (CURRITER .GT. 1) THEN
        IF (TAU .LE. CUMDT .AND. CURRITER .GT. 0) THEN
          
           CALL NVTANDERSEN(ITER)
           !CALL VELVERLET(CURRITER)
           CUMDT = ZERO
           FLAGAND = ZERO
        ELSE

           CALL VELVERLET(CURRITER)
           FLAGAND = ONE

        ENDIF

     ELSEIF (NVTON .EQ. 4 .AND. ITER .LE. THERMRUN) THEN
        
        ! NOSE thermostat

!!$        !
!!$        ! To smooth things out, let's average the
!!$        ! temperature over the previous AVEPER time steps
!!$        !
!!$
!!$        CALL AVETEMP

        IF (MOD(ITER, THERMPER) .EQ. 0 .AND. CURRITER .GT. 0) THEN

           CALL NVTNH

        ELSE

           CALL VELVERLET(CURRITER)

        ENDIF


     ELSEIF (NPTON .EQ. 1 .AND. ITER .LE. THERMRUN) THEN
        
        ! Velocity rescaling thermostat

        !
        ! To smooth things out, let's average the
        ! temperature over the previous AVEPER time steps
        !

        CALL AVETEMP
        CALL AVEPRESS

        IF (MOD(ITER, THERMPER) .EQ. 0 .AND. CURRITER .GT. 1) THEN

           CALL NPTRESCALE
           CALL NEBLISTS(1)
           CALL GETDENSITY

        ELSE

           CALL VELVERLET(CURRITER)
           
        ENDIF

     ENDIF

     IF (MOD(ITER,WRTFREQ) .EQ. 0) THEN

        ! The non-orthogonal density matrix is obtained from
        ! GETMDF (when the Pulay force is computed)

        CALL PAIRPOT

        CALL GETKE

        CALL TOTENG

        ECOUL = ZERO
        IF (ELECTRO .EQ. 1) CALL GETCOULE

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

        IF (ABS(TOTE) .GT. 10000000000.0) THEN
           WRITE(6,*) "# The calculation has diverged - check SCF parameters"
           STOP
        ENDIF

        THETIME = FLOAT(ITER)*DT/THOUSAND

        IF (NPTON .EQ. 0) THEN
           WRITE(6,99) THETIME, TOTE, TEMPERATURE, PRESSURE, TOTE + CONSMOT, EGAP, CHEMPOT, &
                STRTEN(1), STRTEN(2), STRTEN(3), STRTEN(4), STRTEN(5), STRTEN(6)
        ELSE
           WRITE(6,99) THETIME, TOTE, TEMPERATURE, PRESSURE, &
                MASSDEN, EGAP, TOTE + CONSMOT
        ENDIF

99      FORMAT(20G18.9)

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

  IF (NVTON .GE. 1) DEALLOCATE(THIST)
  IF (NPTON .EQ. 1) DEALLOCATE(THIST, PHIST)
  IF (NVTON .EQ. 2) DEALLOCATE(FRANPREV)

  RETURN

END SUBROUTINE TBMD
