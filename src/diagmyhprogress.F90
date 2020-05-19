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

SUBROUTINE DIAGMYHPRG()

#ifdef PROGRESSON

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE DIAGARRAY
  USE SPINARRAY
  USE NONOARRAY
  USE KSPACEARRAY
  USE MYPRECISION
  USE BML

  IMPLICIT NONE

  INTEGER :: INFO
  INTEGER :: I, J, K, M
  CHARACTER(LEN=1), PARAMETER :: JOBZ = 'V',  UPLO = 'U'
  TYPE(BML_MATRIX_T) :: ORTHOH_BML, EVECS_BML
  TYPE(BML_MATRIX_T) :: ORTHOHUP_BML, ORTHOHDOWN_BML 
  TYPE(BML_MATRIX_T) :: UPEVECS_BML, DOWNEVECS_BML 

  IF (EXISTERROR) RETURN

  IF (VERBOSE >= 1) WRITE(*,*)"In DIAGMYHPRG ..."

  IF (SPINON .EQ. 0) THEN

     CALL BML_ZERO_MATRIX(BML_MATRIX_DENSE, BML_ELEMENT_REAL, &
          LATTEPREC, HDIM, HDIM, ORTHOH_BML)
     CALL BML_ZERO_MATRIX(BML_MATRIX_DENSE, BML_ELEMENT_REAL, &
          LATTEPREC, HDIM, HDIM, EVECS_BML)

     IF (BASISTYPE .EQ. "ORTHO") THEN
        CALL BML_IMPORT_FROM_DENSE(BML_MATRIX_DENSE, &
             H, ORTHOH_BML, ZERO, HDIM)
     ELSE
        CALL BML_IMPORT_FROM_DENSE(BML_MATRIX_DENSE, &
             ORTHOH, ORTHOH_BML, ZERO, HDIM)
     ENDIF

     CALL BML_IMPORT_FROM_DENSE(BML_MATRIX_DENSE, &
          ORTHOH, ORTHOH_BML, ZERO, HDIM)

     CALL BML_DIAGONALIZE(ORTHOH_BML,EVALS,EVECS_BML)

     CALL BML_EXPORT_TO_DENSE(EVECS_BML, EVECS)
     IF (DOKERNEL .EQV. .TRUE.) CALL BML_EXPORT_TO_DENSE(ORTHOH_BML, ORTHOH)
     IF (DFTBU .EQV. .TRUE.) CALL BML_EXPORT_TO_DENSE(ORTHOH_BML, ORTHOH)

     CALL BML_DEALLOCATE(ORTHOH_BML)
     CALL BML_DEALLOCATE(EVECS_BML)


  ELSE

     !Up evecs
     CALL BML_ZERO_MATRIX(BML_MATRIX_DENSE, BML_ELEMENT_REAL, &
          LATTEPREC, HDIM, HDIM, ORTHOHUP_BML)
     CALL BML_ZERO_MATRIX(BML_MATRIX_DENSE, BML_ELEMENT_REAL, &
          LATTEPREC, HDIM, HDIM, UPEVECS_BML)
 
     IF (BASISTYPE .EQ. "ORTHO") THEN
        CALL BML_IMPORT_FROM_DENSE(BML_MATRIX_DENSE, &
             HUP, ORTHOHUP_BML, ZERO, HDIM)
     ELSE
        CALL BML_IMPORT_FROM_DENSE(BML_MATRIX_DENSE, &
             ORTHOHUP, ORTHOHUP_BML, ZERO, HDIM)
     ENDIF

     CALL BML_DIAGONALIZE(ORTHOHUP_BML,UPEVALS,UPEVECS_BML)

     CALL BML_EXPORT_TO_DENSE(UPEVECS_BML, UPEVECS)

     CALL BML_DEALLOCATE(ORTHOHUP_BML)
     CALL BML_DEALLOCATE(UPEVECS_BML)


     !Down evecs
     CALL BML_ZERO_MATRIX(BML_MATRIX_DENSE, BML_ELEMENT_REAL, &
          LATTEPREC, HDIM, HDIM, ORTHOHDOWN_BML)
     CALL BML_ZERO_MATRIX(BML_MATRIX_DENSE, BML_ELEMENT_REAL, &
          LATTEPREC, HDIM, HDIM, DOWNEVECS_BML)
 
     IF (BASISTYPE .EQ. "ORTHO") THEN
        CALL BML_IMPORT_FROM_DENSE(BML_MATRIX_DENSE, &
             HUP, ORTHOHDOWN_BML, ZERO, HDIM)
     ELSE
        CALL BML_IMPORT_FROM_DENSE(BML_MATRIX_DENSE, &
             ORTHOHDOWN, ORTHOHDOWN_BML, ZERO, HDIM)
     ENDIF

     CALL BML_DIAGONALIZE(ORTHOHDOWN_BML,DOWNEVALS,DOWNEVECS_BML)

     CALL BML_EXPORT_TO_DENSE(DOWNEVECS_BML, DOWNEVECS)

     CALL BML_DEALLOCATE(ORTHOHDOWN_BML)
     CALL BML_DEALLOCATE(DOWNEVECS_BML)

  ENDIF

  RETURN

#endif

END SUBROUTINE DIAGMYHPRG
