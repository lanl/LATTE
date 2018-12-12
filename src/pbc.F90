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

SUBROUTINE PBC

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, IPIV(3), INFO
  REAL(LATTEPREC) :: WORK(3), BOXINV(3,3), S(3)
  IF (EXISTERROR) RETURN

  ! reduced coordinates s = box^-1 X real coordinates

  BOXINV = BOX

#ifdef DOUBLEPREC

  CALL DGETRF(3, 3, BOXINV, 3, IPIV, INFO)

  CALL DGETRI(3, BOXINV, 3, IPIV, WORK, 3, INFO)

#elif defined(SINGLEPREC)

  CALL SGETRF(3, 3, BOXINV, 3, IPIV, INFO)

  CALL SGETRI(3, BOXINV, 3, IPIV, WORK, 3, INFO)

#endif

  DO I = 1, NATS

#ifdef DOUBLEPREC     

     CALL DGEMV('T', 3, 3, ONE, BOXINV, 3, CR(1,I), 1, ZERO, S, 1)

#elif defined(SINGLEPREC)

     CALL SGEMV('T', 3, 3, ONE, BOXINV, 3, CR(1,I), 1, ZERO, S, 1)

#endif

     ! If we're outside the box, add or subtract the corresponding vector accordingly

     IF (S(1) .GT. ONE) CR(:,I) = CR(:,I) - BOX(1,:) 
     IF (S(1) .LT. ZERO) CR(:,I) = CR(:,I) + BOX(1,:)

     IF (S(2) .GT. ONE) CR(:,I) = CR(:,I) - BOX(2,:)
     IF (S(2) .LT. ZERO) CR(:,I) = CR(:,I) + BOX(2,:)

     IF (S(3) .GT. ONE) CR(:,I) = CR(:,I) - BOX(3,:)
     IF (S(3) .LT. ZERO) CR(:,I) = CR(:,I) + BOX(3,:)

  ENDDO


  RETURN

END SUBROUTINE PBC
