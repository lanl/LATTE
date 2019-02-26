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

SUBROUTINE READPPOTPLUSD()

  USE CONSTANTS_MOD
  USE PPOTARRAY

  IMPLICIT NONE

  INTEGER :: I, J, K
  CHARACTER(LEN=20) :: HD
  LOGICAL :: FILEEXISTS

  IF (EXISTERROR) RETURN

  IF (BASISTYPE .EQ. "ORTHO") THEN
     INQUIRE( FILE=TRIM(PARAMPATH)//"/ppotsplusD.ortho", exist=FILEEXISTS)
     IF (.NOT. FILEEXISTS) THEN
        CALL ERRORS("readppot","ppot.ortho file does not exist. &
             & Please either set PPOTON= 0 or add a file for the pair potentials.")
     ELSE
        OPEN(UNIT=14,STATUS="OLD", FILE=TRIM(PARAMPATH)//"/ppotsplusD.ortho")
     END IF
  ELSEIF (BASISTYPE .EQ. "NONORTHO") THEN
     INQUIRE( FILE=TRIM(PARAMPATH)//"/ppotsplusD.nonortho", exist=FILEEXISTS)
     IF (.NOT. FILEEXISTS) THEN
        CALL ERRORS("readppot","ppot.ortho file does not exist. &
             & Please either set PPOTON= 0 or add a file for the pair potentials.")
     ELSE
        OPEN(UNIT=14, STATUS="OLD", FILE=TRIM(PARAMPATH)//"/ppotsplusD.nonortho")
     END IF
  END IF

  READ(14,*) HD, NOPPD

  ALLOCATE(PPELE(NOPPD), RZERO(NOPPD), C6(NOPPD))

  READ(14,*) HD, PLUSDS6
  READ(14,*) HD, PLUSDGAMMA
  READ(14,*) HD, PLUSDCUT
  READ(14,*) HD, HD, HD
  DO I = 1, NOPPD
     READ(14,*) PPELE(I), RZERO(I), C6(I)
  ENDDO

  RETURN

END SUBROUTINE READPPOTPLUSD
