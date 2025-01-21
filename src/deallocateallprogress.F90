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

SUBROUTINE DEALLOCATEALLPROGRESS

  USE SETUPARRAY
  USE nonoarrayprogress
  USE CONSTANTS_MOD
  USE GENXPROGRESS
  USE BML 

  IMPLICIT NONE
  IF (EXISTERROR) RETURN

  IF (BML_ALLOCATED (HAM_BML)) CALL BML_DEALLOCATE(HAM_BML)
  IF (BML_ALLOCATED (ZMAT_BML)) CALL BML_DEALLOCATE(ZMAT_BML)
  IF (BML_ALLOCATED (OVER_BML)) CALL BML_DEALLOCATE(OVER_BML)
  IF (BML_ALLOCATED (ORTHOH_BML)) CALL BML_DEALLOCATE(ORTHOH_BML)
  IF (BML_ALLOCATED (ORTHOBO_BML)) CALL BML_DEALLOCATE(ORTHOBO_BML)
  IF (BML_ALLOCATED (BO_BML)) CALL BML_DEALLOCATE(BO_BML)

  RETURN

END SUBROUTINE DEALLOCATEALLPROGRESS
