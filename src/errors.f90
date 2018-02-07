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

!> Soubroutine to handle errors.
!! \brief This will decide what to do if the program has to stop for some
!! reason. It will either stop the program or send an error flag up.
!! \param SUB Subroutine where the program would stop.
!! \param TAG Error message to be passed.
!!
SUBROUTINE ERRORS(SUB,TAG)
  USE CONSTANTS_MOD, ONLY: LIBRUN, EXISTERROR, LIBCALLS
  USE TIMER_MOD, ONLY: TIMEDATE_TAG
  USE MYPRECISION, ONLY: LATTEPREC
  USE LATTE_LIB, ONLY: LATTE
  CHARACTER(*), INTENT(IN) :: SUB,TAG

  IF(LIBRUN)THEN
     WRITE(*,*) ""
     WRITE(*,*) "# *************ERROR**************"
     WRITE(*,*) "# LATTE stopped due to the following error at subroutine ",SUB,":"
     WRITE(*,*) "# ",TAG
     CALL TIMEDATE_TAG("The error occurred at: ")
     WRITE(*,*) "# The error will be reported back to the host code"
     WRITE(*,*) "# ********************************"
     WRITE(*,*) ""
     CALL FLUSH(6)
     EXISTERROR = .TRUE.
  ELSE
     WRITE(*,*) ""
     WRITE(*,*) "# *************ERROR**************"
     WRITE(*,*) "# LATTE stopped due to the following error at subroutine ",SUB,":"
     WRITE(*,*) "# ",TAG
     CALL TIMEDATE_TAG("The error occurred at: ")
     WRITE(*,*) "# ********************************"
     WRITE(*,*) ""
     STOP
  ENDIF

END SUBROUTINE ERRORS
