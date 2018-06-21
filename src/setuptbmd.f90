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

SUBROUTINE SETUPTBMD(NEWSYSTEM)

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
  INTEGER :: FLAGAND, NEWSYSTEM

  IF (EXISTERROR) RETURN

  !
  ! Read MDcontroller to determine what kind of MD simulation to do
  !

  IF(NEWSYSTEM == 1 .OR. (.NOT.LIBINIT))THEN

     IF (LATTEINEXISTS) THEN
        CALL PARSE_MD("latte.in")
     ELSE
        CALL READMDCONTROLLER
     ENDIF

     !
     ! Allocate stuff for building the neighbor lists, then build them
     !
     CALL ALLOCATENEBARRAYS

     CALL FLUSH(6)

  ENDIF

  CALL NEBLISTS(0)

  !
  ! Allocate things depending on which method we're using
  ! to get the bond-order
  !

  IF(NEWSYSTEM == 1 .OR. .NOT.LIBINIT)THEN
     IF (CONTROL .EQ. 1) THEN
        CALL ALLOCATEDIAG
     ELSEIF (CONTROL .EQ. 2 .OR. CONTROL .EQ. 4 .OR. CONTROL .EQ. 5) THEN
        CALL ALLOCATEPURE
     ELSEIF (CONTROL .EQ. 3) THEN
        CALL FERMIALLOCATE
     ENDIF
  ENDIF

  !
  IF(VERBOSE >= 1)WRITE(*,*)"Getting MD forces ..."
  IF (RESTART .EQ. 0) CALL GETMDF(0,1)

  CUMDT = ZERO

  TOTSCF = 0

  RETURN

END SUBROUTINE SETUPTBMD
