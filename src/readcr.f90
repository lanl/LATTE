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

SUBROUTINE READCR

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MDARRAY
  USE NEBLISTARRAY
  USE KSPACEARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, MYINDEX
  REAL(LATTEPREC) :: LN(6)
  CHARACTER(LEN=20) :: HEADER
  IF (EXISTERROR) RETURN

  IF(.NOT.ALLOCATED(CR))THEN

     OPEN(UNIT=12, STATUS="OLD", FILE=COORDSFILE)

     READ(12,*) NATS

     ALLOCATE(CR(3,NATS), ATELE(NATS), F(3,NATS), FPP(3,NATS), FTOT(3,NATS))
     ALLOCATE(DELTAQ(NATS), MYCHARGE(NATS))
     ALLOCATE(ELEMPOINTER(NATS))

     IF (ELECTRO .EQ. 0) THEN
        ALLOCATE(LCNSHIFT(NATS))
        LCNSHIFT = ZERO
     ENDIF

     IF (KON .EQ. 1) ALLOCATE(KF(3,NATS))

     IF (BASISTYPE .EQ. "NONORTHO") THEN
        ALLOCATE(FPUL(3,NATS))
        IF (ELECTRO .EQ. 0) ALLOCATE(FSLCN(3,NATS))
        IF (ELECTRO .EQ. 1) ALLOCATE(FSCOUL(3,NATS))
        IF (SPINON .EQ. 1) ALLOCATE(FSSPIN(3,NATS))
     ENDIF

     READ(12,*) BOX(1,1), BOX(1,2), BOX(1,3)
     READ(12,*) BOX(2,1), BOX(2,2), BOX(2,3)
     READ(12,*) BOX(3,1), BOX(3,2), BOX(3,3)

     DO I = 1, NATS
        READ(12,*) ATELE(I), CR(1,I), CR(2,I), CR(3,I)
     ENDDO

     CLOSE(12)

  ELSE

     ALLOCATE(F(3,NATS), FPP(3,NATS), FTOT(3,NATS))
     ALLOCATE(DELTAQ(NATS), MYCHARGE(NATS))
     ALLOCATE(ELEMPOINTER(NATS))

     IF (ELECTRO .EQ. 0) THEN
        ALLOCATE(LCNSHIFT(NATS))
        LCNSHIFT = ZERO
     ENDIF

     IF (KON .EQ. 1) ALLOCATE(KF(3,NATS))

     IF (BASISTYPE .EQ. "NONORTHO") THEN
        IF (SPINON .EQ. 0) THEN
           ALLOCATE(FPUL(3,NATS), FSCOUL(3,NATS))
        ELSE
           ALLOCATE(FPUL(3,NATS), FSCOUL(3,NATS), FSSPIN(3,NATS))
        ENDIF
     ENDIF

  ENDIF

  ! Set up pointer to the data in TBparam/electrons.dat
  ELEMPOINTER = 0

  DO I = 1, NATS
     DO J = 1, NOELEM
        IF (ATELE(I) .EQ. ELE(J)) ELEMPOINTER(I) = J
     ENDDO
  ENDDO

  ! Let's compute the total mass - this comes in handy when computing the
  ! density later

  SUMMASS = ZERO

  DO I = 1, NATS
     SUMMASS = SUMMASS + MASS(ELEMPOINTER(I))
  ENDDO

  ! Let's check whether we have only sp elements. If so, we can
  ! use a much faster verison of gradH

  ! SPONLY = 0: use GRADHSP
  ! SPONLY = 1: use Josh Coe's implementation of the automatic H build

  SPONLY = 0
  DO I = 1, NATS
     IF (BASIS(ELEMPOINTER(I)) .NE. "s" .AND. &
          BASIS(ELEMPOINTER(I)) .NE. "sp") SPONLY = 1
  ENDDO

  ! At the moment we must run Josh's implementation if we do tabulated integrals

  IF (SCLTYPE .EQ. "TABLE") SPONLY = 1

  ! Print a warning

  IF (SPONLY .EQ. 1) THEN
     PRINT*, "#FYI: d or f orbitals detected so we're using the"
     PRINT*, "#slower, general SK expansions"
  ENDIF

  !  SPONLY = 1

  ! If we're enforcing LCN and we're using diagonalization

  IF (CONTROL .EQ. 1 .AND. ELECTRO .EQ. 0) ALLOCATE(RESPCHI(NATS))

  RETURN

END SUBROUTINE READCR
