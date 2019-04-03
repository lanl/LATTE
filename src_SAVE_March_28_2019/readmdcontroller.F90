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

SUBROUTINE READMDCONTROLLER

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MDARRAY
  USE NEBLISTARRAY

  IMPLICIT NONE

  CHARACTER(LEN=20) :: HD
  IF (EXISTERROR) RETURN

  OPEN (UNIT=15, STATUS="OLD", FILE="MDcontroller")

  !
  ! MAXITER = run this many MD time steps
  !

  READ(15,*) HD, MAXITER

  !
  ! UDNEIGH = update the neighbor lists every UDNEIGH time steps
  !

  READ(15,*) HD, UDNEIGH

  !
  ! DT = size of the time step in fs
  !

  READ(15,*) HD, DT

  !
  ! TTARGET = temperature in K were initialize and aim for during NVT MD
  ! RNDIST = Type of distribution of random numbers used to initialize T
  !        = GAUSSIAN or UNIFORM
  ! SEEDINIT = Type of seed used in the generation of random numbers
  !          = RANDOM - seed changes every time
  !          = DEFAULT - use the same, default seed every time
  !

  READ(15,*) HD, TTARGET, HD, RNDIST, HD, SEEDINIT

  !
  ! DUMPFREQ: Write a dump file every DUMPFREQ time steps
  !

  READ(15,*) HD, DUMPFREQ

  !
  ! RSFREQ: Write a restart file every RSFREQ time steps
  !

  READ(15,*) HD, RSFREQ

  !
  ! WRTFREQ: Output energy and temperature every WRTFREQ time steps
  !

  READ(15,*) HD, WRTFREQ

  IF ( WRTFREQ <= 0 ) THEN 
    CALL ERRORS("latteparser_latte_mod","You cannot have WRTFREQ <= 0.& 
                &Set this variable to a very high value to avoid frequent printing")
  ENDIF


  !
  ! TOINITTEMP: Whether or not we are going to initialize velocities
  ! using a random number generator (sometimes during a restart we
  ! may not want to reinitialize the temperature
  !

  READ(15,*) HD, TOINITTEMP

  !
  ! THERMPER: If we're running NVT, rescale velocities every THERMPER
  ! time steps.
  !

  READ(15,*) HD, THERMPER

  !
  ! THERMRUN: Thermalize over this many time steps when NVT is on
  !

  READ(15,*) HD, THERMRUN

  !
  ! NVTON: 0 = running NVE MD, 1 = running NVT MD
  ! AVEPER: Average the temperature over AVEPER time steps when determining
  ! how to rescale velocities
  !

  READ(15,*) HD, NVTON, HD, NPTON, HD, AVEPER, HD, FRICTION, &
       HD, SEEDTH

  IF (NVTON .EQ. 1 .AND. NPTON .EQ. 1) THEN
     CALL ERRORS("readmdcontroller","You can't have NVTON = 1 and NPTON = 1")
  ENDIF

  ! PTARGET = Target pressure (in GPa) when running NPT
  ! NPTTYPE = ISO or ANISO

  READ(15,*) HD, PTARGET, HD, NPTTYPE

  !
  ! The following are for the Hugoniostat
  !

  ! On (1) or off (0)?

  READ(15,*) HD, SHOCKON

  !
  ! SHOCKSTART = the MD iteration where we will start to compress
  ! the iteration when we stop depends on the size of the block and Us
  !

  READ(15,*) HD, SHOCKSTART

  !
  ! SHOCKDIR is the cartensian direction (1 = X, 2 = Y, 3 = Z),
  ! parallel to which we're going to compress uniaxially
  !

  READ(15,*) HD, SHOCKDIR

  !
  ! And finally, the particle and shock velocities
  ! IN UNITS OF METRES PER SECOND
  !

  READ(15,*) HD, UPARTICLE, HD, USHOCK, HD, C0

  ! Adapt SCF on the fly?

  READ(15,*) HD, MDADAPT

  ! Calculating Hugoniot points?

  READ(15,*) HD, GETHUG, HD, E0, HD, V0, HD, P0


  CLOSE(15)

  RETURN

END SUBROUTINE READMDCONTROLLER
