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

SUBROUTINE ORTHOMYRHO

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE NONOARRAY
  USE SPINARRAY
  USE MYPRECISION

  IMPLICIT NONE
  IF (EXISTERROR) RETURN

  !
  ! ORTHORHO = X^dag RHO X
  !

  IF (SPINON .EQ. 0) THEN

#ifdef DOUBLEPREC

     CALL DGEMM('T', 'N', HDIM, HDIM, HDIM, ONE, &
          XMAT, HDIM, BO, HDIM, ZERO, NONOTMP, HDIM)

     CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
          NONOTMP, HDIM, XMAT, HDIM, ZERO, BO, HDIM)

#elif defined(SINGLEPREC)

     CALL SGEMM('T', 'N', HDIM, HDIM, HDIM, ONE, &
          XMAT, HDIM, BO, HDIM, ZERO, NONOTMP, HDIM)

     CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
          NONOTMP, HDIM, XMAT, HDIM, ZERO, BO, HDIM)

#endif

  ELSE

#ifdef DOUBLEPREC

     CALL DGEMM('T', 'N', HDIM, HDIM, HDIM, ONE, &
          XMAT, HDIM, RHOUP, HDIM, ZERO, NONOTMP, HDIM)

     CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
          NONOTMP, HDIM, XMAT, HDIM, ZERO, RHOUP, HDIM)

     CALL DGEMM('T', 'N', HDIM, HDIM, HDIM, ONE, &
          XMAT, HDIM, RHODOWN, HDIM, ZERO, NONOTMP, HDIM)

     CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
          NONOTMP, HDIM, XMAT, HDIM, ZERO, RHODOWN, HDIM)

#elif defined(SINGLEPREC)

     CALL SGEMM('T', 'N', HDIM, HDIM, HDIM, ONE, &
          XMAT, HDIM, RHOUP, HDIM, ZERO, NONOTMP, HDIM)

     CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
          NONOTMP, HDIM, XMAT, HDIM, ZERO, RHOUP, HDIM)

     CALL SGEMM('T', 'N', HDIM, HDIM, HDIM, ONE, &
          XMAT, HDIM, RHODOWN, HDIM, ZERO, NONOTMP, HDIM)

     CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
          NONOTMP, HDIM, XMAT, HDIM, ZERO, RHODOWN, HDIM)     

#endif

  ENDIF

  RETURN

END SUBROUTINE ORTHOMYRHO

