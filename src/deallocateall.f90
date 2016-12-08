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

SUBROUTINE DEALLOCATEALL

  USE SETUPARRAY
  USE UNIVARRAY
  USE PPOTARRAY
  USE NEBLISTARRAY
  USE MDARRAY
  USE KSPACEARRAY
  USE CONSTANTS_MOD
  USE SPINARRAY

  IMPLICIT NONE

  IF (MDON .EQ. 1) DEALLOCATE(V)

  DEALLOCATE(CR, BTYPE, HES, HEP, HED, HEF, ATOCC, &
       ELE, ELE1, ELE2, BASIS, ATELE, MASS, HUBBARDU, DELTAQ, QLIST, MYCHARGE)
  DEALLOCATE(ELEMPOINTER)

  IF (SPINON .EQ. 0) THEN
     DEALLOCATE(BOZERO)
     IF (KON .EQ. 0) DEALLOCATE(BO)
  ELSEIF (SPINON .EQ. 1) THEN
     DEALLOCATE(RHOUP, RHODOWN, HUP, HDOWN, DELTASPIN, OLDDELTASPIN, &
          WSS, WPP, WDD, WFF, RHOUPZERO, RHODOWNZERO)
  ENDIF

  DEALLOCATE(BOND)
  IF (PPOTON .EQ. 1) DEALLOCATE(POTCOEF, PPELE1, PPELE2)
  DEALLOCATE(F, FPP, FTOT)

  IF (KON .EQ. 0) THEN

     DEALLOCATE(H, HDIAG)

  ELSE

     DEALLOCATE(HK, HKDIAG, KBO)

  ENDIF

  IF (BASISTYPE .EQ. "NONORTHO") THEN
     
     DEALLOCATE(OVERL)
     IF (SPINON .EQ. 0) THEN
        DEALLOCATE(FPUL, FSCOUL)
     ELSE
        DEALLOCATE(FPUL, FSCOUL, FSSPIN)
     ENDIF
  ENDIF
   
  DEALLOCATE(NEBTB)
  IF (PPOTON .EQ. 1) DEALLOCATE(NEBPP)
  IF (ELECTRO .EQ. 1) DEALLOCATE(NEBCOUL)
  IF (ELECTRO .EQ. 0) DEALLOCATE(LCNSHIFT)
  

  RETURN

END SUBROUTINE DEALLOCATEALL
