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

SUBROUTINE ORTHOMYHPRG

#ifdef  PROGRESSON

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE NONOARRAY
  USE SPINARRAY
  USE MYPRECISION

  USE BML
  USE PRG_NONORTHO_MOD
  USE PRG_EXTRAS_MOD
  USE GENXPROGRESS
  USE LATTEPARSER_LATTE_MOD

  IMPLICIT NONE

  INTEGER I, J
  REAL(LATTEPREC) :: MLSI

  TYPE(BML_MATRIX_T) :: ORTHOH_BML, HAM_BML
  IF (EXISTERROR) RETURN

  !
  ! ORTHOH = X^dag H X
  !
  MLSI = MLS()
  IF (SPINON .EQ. 0) THEN
     !! Convert Hamiltonian to bml format
     IF(LT%MDIM == -1)LT%MDIM = HDIM

     CALL BML_ZERO_MATRIX(LT%BML_TYPE, BML_ELEMENT_REAL, &
          LATTEPREC, HDIM, LT%MDIM, ORTHOH_BML)

     CALL BML_IMPORT_FROM_DENSE(LT%BML_TYPE, &
          H, HAM_BML, LT%THRESHOLD, LT%MDIM)

     CALL PRG_ORTHOGONALIZE(HAM_BML,ZMAT_BML,ORTHOH_BML,&
          LT%THRESHOLD,LT%BML_TYPE,LT%VERBOSE)

     CALL BML_EXPORT_TO_DENSE(ORTHOH_BML, ORTHOH)

     CALL BML_DEALLOCATE(HAM_BML)
     CALL BML_DEALLOCATE(ORTHOH_BML)

     IF (DEBUGON .EQ. 1) THEN

        OPEN (UNIT=33, STATUS="UNKNOWN", FILE="myXHX.dat")

        DO I = 1, HDIM

           WRITE(33,10) (ORTHOH(I,J), J = 1, HDIM)

        ENDDO

        CLOSE(33)

10      FORMAT(100G18.9)

     ENDIF

  ELSE

     CALL BML_ZERO_MATRIX(LT%BML_TYPE, BML_ELEMENT_REAL, &
          LATTEPREC, HDIM, LT%MDIM, ORTHOH_BML)

     CALL BML_IMPORT_FROM_DENSE(LT%BML_TYPE, &
          HUP, HAM_BML, ZERO, LT%MDIM)

     CALL PRG_ORTHOGONALIZE(HAM_BML,ZMAT_BML,ORTHOH_BML,&
          LT%THRESHOLD,LT%BML_TYPE,LT%VERBOSE)

     CALL BML_EXPORT_TO_DENSE(ORTHOH_BML, ORTHOHUP)

     CALL BML_IMPORT_FROM_DENSE(LT%BML_TYPE, &
          HDOWN, HAM_BML, ZERO, LT%MDIM)

     CALL PRG_ORTHOGONALIZE(HAM_BML,ZMAT_BML,ORTHOH_BML,&
          LT%THRESHOLD,LT%BML_TYPE,LT%VERBOSE)

     CALL BML_EXPORT_TO_DENSE(ORTHOH_BML, ORTHOHDOWN)

     CALL BML_DEALLOCATE(HAM_BML)
     CALL BML_DEALLOCATE(ORTHOH_BML)

  ENDIF

  IF (VERBOSE >= 1) WRITE(*,*)"Time for ORTHOMYH =", MLS()-MLSI

  RETURN

#endif

END SUBROUTINE ORTHOMYHPRG
