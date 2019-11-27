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

SUBROUTINE DEORTHOMYRHOPRG

#ifdef  PROGRESSON

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE NONOARRAY
  USE SPINARRAY
  USE MYPRECISION
  USE LATTEPARSER_LATTE_MOD

  USE BML
  USE PRG_NONORTHO_MOD
  USE GENXPROGRESS

  IMPLICIT NONE

  TYPE(BML_MATRIX_T) :: ORTHOBO_BML, BO_BML

  IF (EXISTERROR) RETURN

  !
  ! RHO = X ORTHORHO X^dag
  !

  IF (SPINON .EQ. 0) THEN

     !! Convert Density matrix to bml format
     CALL BML_ZERO_MATRIX(LT%BML_TYPE, BML_ELEMENT_REAL, &
          LATTEPREC, HDIM, LT%MDIM, ORTHOBO_BML)

     CALL BML_IMPORT_FROM_DENSE(LT%BML_TYPE, &
          BO, BO_BML, ZERO, LT%MDIM)

     CALL PRG_DEORTHOGONALIZE(BO_BML,ZMAT_BML,ORTHOBO_BML,&
          LT%THRESHOLD,LT%BML_TYPE,LT%VERBOSE)

     CALL BML_EXPORT_TO_DENSE(ORTHOBO_BML, BO)

     CALL BML_DEALLOCATE(BO_BML)
     CALL BML_DEALLOCATE(ORTHOBO_BML)


  ELSE

     !! For rhoup
     CALL BML_ZERO_MATRIX(LT%BML_TYPE, BML_ELEMENT_REAL, &
          LATTEPREC, HDIM, LT%MDIM, ORTHOBO_BML)

     CALL BML_IMPORT_FROM_DENSE(LT%BML_TYPE, &
          RHOUP, BO_BML, ZERO, LT%MDIM)

     CALL PRG_DEORTHOGONALIZE(BO_BML,ZMAT_BML,ORTHOBO_BML,&
          LT%THRESHOLD,LT%BML_TYPE,LT%VERBOSE)

     CALL BML_EXPORT_TO_DENSE(ORTHOBO_BML, RHOUP)

     CALL BML_DEALLOCATE(BO_BML)
     CALL BML_DEALLOCATE(ORTHOBO_BML)

     !! For rhodown
     CALL BML_ZERO_MATRIX(LT%BML_TYPE, BML_ELEMENT_REAL, &
          LATTEPREC, HDIM, LT%MDIM, ORTHOBO_BML)

     CALL BML_IMPORT_FROM_DENSE(LT%BML_TYPE, &
          RHODOWN, BO_BML, ZERO, LT%MDIM)

     CALL PRG_DEORTHOGONALIZE(BO_BML,ZMAT_BML,ORTHOBO_BML,&
          LT%THRESHOLD,LT%BML_TYPE,LT%VERBOSE)

     CALL BML_EXPORT_TO_DENSE(ORTHOBO_BML, RHODOWN)

     CALL BML_DEALLOCATE(BO_BML)
     CALL BML_DEALLOCATE(ORTHOBO_BML)

  ENDIF

  RETURN

#endif

END SUBROUTINE DEORTHOMYRHOPRG
