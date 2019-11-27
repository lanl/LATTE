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

SUBROUTINE GETFORCE

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MYPRECISION

  IMPLICIT NONE
  IF (EXISTERROR) RETURN

  FTOT = ZERO

  IF (KON .EQ. 0) THEN

     IF (BASISTYPE .EQ. "ORTHO") THEN
        
        CALL GRADH
        
        FTOT = TWO * F
        
     ELSEIF (BASISTYPE .EQ. "NONORTHO") THEN
        
        CALL TBFORCES
        
        FTOT = F + FPUL

     ENDIF

     
  ELSEIF (KON .EQ. 1) THEN

     IF (BASISTYPE .EQ. "ORTHO") THEN
        
        CALL KGRADH
        
        FTOT = TWO*F

     ELSEIF (BASISTYPE .EQ. "NONORTHO") THEN

        CALL KTBFORCES
        
        FTOT = F + FPUL

     ENDIF

  ENDIF
  
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

  IF (PLUSDON .EQ. 1) THEN
     CALL PAIRPOTPLUSD
     FTOT = FTOT + FPLUSD
  ENDIF

  IF (ELECTRO .EQ. 1) FTOT = FTOT + FCOUL

  RETURN

END SUBROUTINE GETFORCE
