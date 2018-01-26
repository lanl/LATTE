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

SUBROUTINE ASSESSOCC

  USE CONSTANTS_MOD
  USE SETUPARRAY  
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, LWORK, INFO
  REAL(LATTEPREC), ALLOCATABLE :: WORK(:), EVECS(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: HEVALS(:), RHOEVALS(:)
  CHARACTER(LEN=1), PARAMETER :: JOBZ = "V",  UPLO = "U" 

  LWORK = 3*HDIM - 1

  ALLOCATE(WORK(LWORK), EVECS(HDIM, HDIM), HEVALS(HDIM), RHOEVALS(HDIM))

  !
  ! First get the eigenvalues of the Hamiltonian
  !

  EVECS = H

#ifdef DOUBLEPREC
  CALL DSYEV(JOBZ, UPLO, HDIM, EVECS, HDIM, HEVALS, WORK, LWORK,INFO)
#elif defined(SINGLEPREC)
  CALL SSYEV(JOBZ, UPLO, HDIM, EVECS, HDIM, HEVALS, WORK, LWORK,INFO)
#endif


  !
  ! Now the eigenvalues of the density matrix
  !

  EVECS = HALF*BO

#ifdef DOUBLEPREC
  CALL DSYEV(JOBZ, UPLO, HDIM, EVECS, HDIM, RHOEVALS, WORK, LWORK,INFO)
#elif defined(SINGLEPREC)
  CALL SSYEV(JOBZ, UPLO, HDIM, EVECS, HDIM, RHOEVALS, WORK, LWORK,INFO)
#endif

  !
  ! ... and print them out
  !

  !
  ! Blas lists the eigenvalues of H in assending order to we'll have 
  ! to flip this:
  !

  OPEN(UNIT=50, STATUS="UNKNOWN", FILE="checkoccupancy.dat")

  DO I = 1, HDIM
     WRITE(50,10) HEVALS(I), ONE-RHOEVALS(I)
  ENDDO

10 FORMAT(2F18.9)

  CLOSE(50)

  DEALLOCATE(WORK, EVECS, HEVALS, RHOEVALS)

  RETURN

END SUBROUTINE ASSESSOCC

