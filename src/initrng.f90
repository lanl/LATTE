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

SUBROUTINE INITRNG

  USE CONSTANTS_MOD

#ifdef MPI_ON
  USE MPI
#endif

  IMPLICIT NONE

  INTEGER :: I, N, COUNT, CLOCK
  INTEGER, ALLOCATABLE :: SEED(:)
  INTEGER :: MYID, IERR
  IF (EXISTERROR) RETURN

  MYID = 0

#ifdef MPI_ON
  CALL MPI_COMM_RANK( MPI_COMM_WORLD, MYID, IERR )
#endif


  CALL RANDOM_SEED(SIZE = N)

  ALLOCATE(SEED(N))

  CALL SYSTEM_CLOCK(COUNT = CLOCK)

  ! Adding MPI rank to the seed should make each rank use a different seed...

  IF (PARREP .EQ. 0) MYID = 0

  DO I = 1, N
     SEED(I) = CLOCK + 37*(I - 1) + MYID
  ENDDO

  CALL RANDOM_SEED(PUT = SEED)

  DEALLOCATE(SEED)

  RETURN

END SUBROUTINE INITRNG
