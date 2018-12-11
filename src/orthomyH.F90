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

SUBROUTINE ORTHOMYH

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE NONOARRAY
  USE SPINARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER I, J
  IF (EXISTERROR) RETURN

  !
  ! ORTHOH = X^dag H X
  !

  !
  ! Don't overwrite H - we update this with new charges and build it 
  ! from scratch only after moving the atoms
  !

  IF (SPINON .EQ. 0) THEN 

#ifdef DOUBLEPREC

     CALL DGEMM('T', 'N', HDIM, HDIM, HDIM, ONE, &
          XMAT, HDIM, H, HDIM, ZERO, NONOTMP, HDIM)
     CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
          NONOTMP, HDIM, XMAT, HDIM, ZERO, ORTHOH, HDIM)

#elif defined(SINGLEPREC)

     CALL SGEMM('T', 'N', HDIM, HDIM, HDIM, ONE, &
          XMAT, HDIM, H, HDIM, ZERO, NONOTMP, HDIM)
     CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
          NONOTMP, HDIM, XMAT, HDIM, ZERO, ORTHOH, HDIM)

#endif  

     IF (DEBUGON .EQ. 1) THEN

        OPEN (UNIT=33, STATUS="UNKNOWN", FILE="myXHX.dat")

        DO I = 1, HDIM

           WRITE(33,10) (ORTHOH(I,J), J = 1, HDIM)

        ENDDO

        CLOSE(33)

10      FORMAT(100G18.9)

     ENDIF

  ELSE

#ifdef DOUBLEPREC

     CALL DGEMM('T', 'N', HDIM, HDIM, HDIM, ONE, &
          XMAT, HDIM, HUP, HDIM, ZERO, NONOTMP, HDIM)
     CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
          NONOTMP, HDIM, XMAT, HDIM, ZERO, ORTHOHUP, HDIM)

     CALL DGEMM('T', 'N', HDIM, HDIM, HDIM, ONE, &
          XMAT, HDIM, HDOWN, HDIM, ZERO, NONOTMP, HDIM)
     CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
          NONOTMP, HDIM, XMAT, HDIM, ZERO, ORTHOHDOWN, HDIM)

#elif defined(SINGLEPREC)

     CALL SGEMM('T', 'N', HDIM, HDIM, HDIM, ONE, &
          XMAT, HDIM, HUP, HDIM, ZERO, NONOTMP, HDIM)
     CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
          NONOTMP, HDIM, XMAT, HDIM, ZERO, ORTHOHUP, HDIM)

     CALL SGEMM('T', 'N', HDIM, HDIM, HDIM, ONE, &
          XMAT, HDIM, HDOWN, HDIM, ZERO, NONOTMP, HDIM)
     CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
          NONOTMP, HDIM, XMAT, HDIM, ZERO, ORTHOHDOWN, HDIM)

#endif

  ENDIF

  RETURN

END SUBROUTINE ORTHOMYH

