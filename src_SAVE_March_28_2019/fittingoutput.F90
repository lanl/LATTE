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

SUBROUTINE FITTINGOUTPUT(BECLEAN)

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE PPOTARRAY
  USE PUREARRAY
  USE SPARSEARRAY
  USE COULOMBARRAY
  USE SPINARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I
  INTEGER, INTENT(IN) :: BECLEAN
  REAL(LATTEPREC) :: MYDIPOLE
  REAL(LATTEPREC) :: BIGNO = 9.9E9
  IF (EXISTERROR) RETURN

  OPEN(UNIT=66, STATUS="UNKNOWN", FILE="fittingoutput.dat")

  IF (BECLEAN .EQ. 0) THEN

     IF (SPINON .EQ. 0) THEN
        TOTE = TRRHOH + EREP - ECOUL - ENTE
     ELSE
        TOTE = TRRHOH + EREP - ECOUL - ENTE + ESPIN
     ENDIF

     WRITE(66,10) TOTE
     DO I = 1, NATS
        WRITE(66,11) FTOT(1,I), FTOT(2,I), FTOT(3,I)
     ENDDO

     CALL GETDIPOLE(MYDIPOLE)

     WRITE(66,12) MYDIPOLE

     WRITE(66,12) EGAP

  ELSEIF (BECLEAN .EQ. 1) THEN ! There was a problem with the calculation

     WRITE(66,10) BIGNO
     DO I = 1, NATS
        WRITE(66,11) BIGNO, BIGNO, BIGNO
     ENDDO

     WRITE(66,12) BIGNO

  ENDIF

10 FORMAT(G18.9)
11 FORMAT(3G18.9)
12 FORMAT(G18.9)

  CLOSE(66)

END SUBROUTINE FITTINGOUTPUT
