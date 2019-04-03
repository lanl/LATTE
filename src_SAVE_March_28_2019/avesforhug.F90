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

SUBROUTINE AVESFORHUG(MYPRESSURE, MYENERGY, MYTEMPERATURE, MYVOL)

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MDARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I
  REAL(LATTEPREC), INTENT(IN) :: MYPRESSURE, MYENERGY, MYTEMPERATURE, MYVOL
  IF (EXISTERROR) RETURN

  ! Record the pressure, volume, and energy

  DO I = 1, (AVEPER/WRTFREQ) - 1

     PHIST(I) = PHIST(I + 1)
     VHIST(I) = VHIST(I + 1)
     EHIST(I) = EHIST(I + 1)
     THIST(I) = THIST(I + 1)

  ENDDO

  PHIST(AVEPER/WRTFREQ) = MYPRESSURE
  VHIST(AVEPER/WRTFREQ) = MYVOL
  EHIST(AVEPER/WRTFREQ) = MYENERGY
  THIST(AVEPER/WRTFREQ) = MYTEMPERATURE

  RETURN

END SUBROUTINE AVESFORHUG
