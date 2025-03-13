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

!> Subroutine for computing the weight (contributions) from the cores to the 
!! eigenvectors of the full system 
!! \brief Given a subsystem, a set of eigenvectors for the subsystem, the number
!! of orbitals for a core region belonging to the subsystem. This routine would 
!! give back an array containing information on the contribution from the core
!! region to every eigenvector. The size of the output is N (total number of 
!! orbitals of the system).  
!! \param SYEVECS 2D array containing the eigenvectors of the system. 
!! \param NCORES Number of orbitals in the core region. 
!!
SUBROUTINE GETDVALS(NCORES)
#ifdef PROGRESSON
  USE BML
#endif
  USE DIAGARRAY
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NCORES
  INTEGER :: NORBS
  INTEGER :: I, J 

  NORBS = SIZE(EVECS, DIM=1)
  IF (ALLOCATED(DVALS)) DEALLOCATE(DVALS)
  ALLOCATE(DVALS(NORBS))
  DVALS = 0.D0

  ! DVALS = SUM(EVECS[1:NCORES, 1:NORBS] ** 2)
  DO I = 1, NORBS 
    DO J = 1, NCORES
        DVALS(I) = DVALS(I) + (EVECS(J, I) ** 2)
    ENDDO
  ENDDO 

  RETURN 

END SUBROUTINE GETDVALS
