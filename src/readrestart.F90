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

SUBROUTINE READRESTART

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MDARRAY
  USE RESTARTARRAY
  USE KSPACEARRAY
#ifdef MPI_ON
  USE MPI
#endif
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, MYINDEX
  INTEGER :: MYID, IERR
  REAL(LATTEPREC) :: LN(6), JUNK, TMP1, TMP2
  CHARACTER(LEN=20) :: HEADER
  CHARACTER(LEN=50) :: FLNM
  IF (EXISTERROR) RETURN


  IF (MDON .EQ. 1 .AND. PARREP .EQ. 1) THEN

#ifdef MPI_ON
     CALL MPI_COMM_RANK( MPI_COMM_WORLD, MYID, IERR )
#endif


     IF (MYID .LT. 10) THEN
        WRITE(FLNM,'(I1,"/bl/restart.dat")') MYID
     ELSEIF (MYID .GE. 10 .AND. MYID .LT. 100) THEN
        WRITE(FLNM,'(I2,"/bl/restart.dat")') MYID
     ELSEIF (MYID .GE. 100 .AND. MYID .LT. 1000) THEN
        WRITE(FLNM,'(I3,"/bl/restart.dat")') MYID
     ENDIF

     PRINT*, MYID, FLNM

     OPEN(UNIT=12, STATUS="OLD", FILE=FLNM)

  ELSE

     ! Regular restart

     OPEN(UNIT=12, STATUS="OLD", FILE="bl/restart.dat")

  ENDIF

  IF (MDON .EQ. 1)  READ(12,*) HEADER, CONTITER
  READ(12,*) NATS

  ALLOCATE(CR(3,NATS), ATELE(NATS), F(3,NATS), FPP(3,NATS), FTOT(3,NATS)) 
  ALLOCATE(DELTAQ(NATS), MYCHARGE(NATS))
  ALLOCATE(ELEMPOINTER(NATS))

  IF (ELECTRO .EQ. 0) THEN
     ALLOCATE(LCNSHIFT(NATS))
     LCNSHIFT = ZERO
  ENDIF

  IF (BASISTYPE .EQ. "NONORTHO") THEN
     IF (SPINON .EQ. 0) THEN
        ALLOCATE(FPUL(3,NATS), FSCOUL(3,NATS))
     ELSE
        ALLOCATE(FPUL(3,NATS), FSCOUL(3,NATS), FSSPIN(3,NATS))
     ENDIF
  ENDIF

  IF (KON .EQ. 1) ALLOCATE(KF(3,NATS))

  IF (MDON .EQ. 1) ALLOCATE(V(3,NATS))

  READ(12,*) BOX(1,1), BOX(1,2), BOX(1,3)
  READ(12,*) BOX(2,1), BOX(2,2), BOX(2,3)
  READ(12,*) BOX(3,1), BOX(3,2), BOX(3,3)


  !  READ(12,*) (ATELE(I), CR(1,I), CR(2,I), CR(3,I), I = 1, NATS)

  DO I = 1, NATS
     READ(12,*) ATELE(I), CR(1,I), CR(2,I), CR(3,I)
  ENDDO


  ! Set up pointer to the data in TBparam/electrons.dat

  DO I = 1, NATS
     DO J = 1, NOELEM
        IF (ATELE(I) .EQ. ELE(J)) ELEMPOINTER(I) = J
     ENDDO
  ENDDO

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

  IF (SCLTYPE .EQ. "TABLE") SPONLY = 1

  READ(12,*) CHEMPOT
  !  PRINT*, CHEMPOT

  IF (SPINON .EQ. 0) THEN

     READ(12,*) TMPHDIM
     ALLOCATE(TMPBODIAG(TMPHDIM))

     DO I = 1, TMPHDIM
        READ(12,*) TMP1, TMP2
        TMPBODIAG(I) = TMP1 + TMP2
     ENDDO

  ELSEIF (SPINON .EQ. 1) THEN

     READ(12,*) TMPHDIM
     ALLOCATE(TMPRHOUP(TMPHDIM), TMPRHODOWN(TMPHDIM))
     READ(12,*) (TMPRHOUP(I), TMPRHODOWN(I), I = 1, TMPHDIM)

  ENDIF

  IF (MDON .EQ. 1) THEN
     DO I = 1, NATS
        READ(12,*) V(1,I), V(2,I), V(3,I)
     ENDDO
  ENDIF

  CLOSE(12)

  !  CR = ALAT * CR

  IF (PBCON .EQ. 0) THEN

     IF (CONTITER .LT. 10) THEN
        WRITE(FLNM,'("animate/myXYZfile_restart.",I1,".xyz")') CONTITER
     ELSEIF (CONTITER .GE. 10 .AND. CONTITER .LT. 100) THEN
        WRITE(FLNM,'("animate/myXYZfile_restart.",I2,".xyz")') CONTITER
     ELSEIF (CONTITER .GE. 100 .AND. CONTITER .LT. 1000) THEN
        WRITE(FLNM,'("animate/myXYZfile_restart.",I3,".xyz")') CONTITER
     ELSEIF (CONTITER .GE. 1000 .AND. CONTITER .LT. 10000) THEN
        WRITE(FLNM,'("animate/myXYZfile_restart.",I4,".xyz")') CONTITER
     ELSEIF (CONTITER .GE. 10000 .AND. CONTITER .LT. 100000) THEN
        WRITE(FLNM,'("animate/myXYZfile_restart.",I5,".xyz")') CONTITER
     ELSEIF (CONTITER .GE. 100000 .AND. CONTITER .LT. 1000000) THEN
        WRITE(FLNM,'("animate/myXYZfile_restart.",I6,".xyz")') CONTITER
     ELSEIF (CONTITER .GE. 1000000 .AND. CONTITER .LT. 10000000) THEN
        WRITE(FLNM,'("animate/myXYZfile_restart.",I7,".xyz")') CONTITER 
     ENDIF

     OPEN(UNIT=23, STATUS="UNKNOWN", FILE=FLNM)

  ENDIF


  RETURN

END SUBROUTINE READRESTART
