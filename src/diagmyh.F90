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

SUBROUTINE DIAGMYH()

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE DIAGARRAY
  USE SPINARRAY
  USE NONOARRAY
  USE KSPACEARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: INFO
  INTEGER :: I, J, K, M
  !  REAL(8), ALLOCATABLE :: TMPHMAT(:,:), TMPEVALS(:), DWORK(:)
  CHARACTER(LEN=1), PARAMETER :: JOBZ = "V",  UPLO = "U"
  IF (EXISTERROR) RETURN

  IF (SPINON .EQ. 0) THEN

     IF (BASISTYPE .EQ. "ORTHO") THEN
        EVECS = H
     ELSE
        EVECS = ORTHOH
     ENDIF

#ifdef XSYEV

#ifdef DOUBLEPREC


     CALL DSYEV(JOBZ, UPLO, HDIM, EVECS, HDIM, EVALS, &
          DIAG_WORK, DIAG_LWORK,INFO)

#elif defined(SINGLEPREC)

     !     ALLOCATE(TMPHMAT(HDIM, HDIM), TMPEVALS(HDIM), DWORK(3*HDIM - 1))
     !     TMPHMAT = DBLE(H)

     CALL SSYEV(JOBZ, UPLO, HDIM, EVECS, HDIM, EVALS, &
          DIAG_WORK, DIAG_LWORK,INFO)
     !     CALL DSYEV(JOBZ, UPLO, HDIM, TMPHMAT, HDIM, TMPEVALS, &
     !          DWORK, DIAG_LWORK,INFO)

     !     EVECS = REAL(TMPHMAT)
     !     EVALS = REAL(TMPEVALS)
     !     DEALLOCATE(TMPHMAT, TMPEVALS, DWORK)


#endif

#elif defined(XSYEVD)

#ifdef DOUBLEPREC
     CALL DSYEVD(JOBZ, UPLO, HDIM, EVECS, HDIM, EVALS, &
          DIAG_WORK, DIAG_LWORK, DIAG_IWORK, DIAG_LIWORK, INFO)

          IF (INFO .NE. 0) THEN

        INFO = 0

        CALL DSYEV(JOBZ, UPLO, HDIM, EVECS, HDIM, EVALS, &
          DIAG_WORK, DIAG_LWORK, INFO)

     ENDIF


#elif defined(SINGLEPREC)
     CALL SSYEVD(JOBZ, UPLO, HDIM, EVECS, HDIM, EVALS, &
          DIAG_WORK, DIAG_LWORK, DIAG_IWORK, DIAG_LIWORK, INFO)
#endif

#endif


  ELSE

     IF (BASISTYPE .EQ. "ORTHO") THEN
        UPEVECS = HUP
        DOWNEVECS = HDOWN
     ELSE
        UPEVECS = ORTHOHUP
        DOWNEVECS = ORTHOHDOWN
     ENDIF

#ifdef XSYEV

#ifdef DOUBLEPREC
     CALL DSYEV(JOBZ, UPLO, HDIM, UPEVECS, HDIM, UPEVALS, &
          DIAG_WORK, DIAG_LWORK,INFO)
#elif defined(SINGLEPREC)
     CALL SSYEV(JOBZ, UPLO, HDIM, UPEVECS, HDIM, UPEVALS, &
          DIAG_WORK, DIAG_LWORK,INFO)
#endif

#ifdef DOUBLEPREC
     CALL DSYEV(JOBZ, UPLO, HDIM, DOWNEVECS, HDIM, DOWNEVALS, &
          DIAG_WORK, DIAG_LWORK,INFO)
#elif defined(SINGLEPREC)
     CALL SSYEV(JOBZ, UPLO, HDIM, DOWNEVECS, HDIM, DOWNEVALS, &
          DIAG_WORK, DIAG_LWORK,INFO)
#endif

#elif defined(XSYEVD)

#ifdef DOUBLEPREC
     CALL DSYEVD(JOBZ, UPLO, HDIM, UPEVECS, HDIM, UPEVALS, &
          DIAG_WORK, DIAG_LWORK,DIAG_IWORK, DIAG_LIWORK, INFO)
     
     IF (INFO .NE. 0) THEN

        INFO = 0

        CALL DSYEV(JOBZ, UPLO, HDIM, UPEVECS, HDIM, UPEVALS, &
             DIAG_WORK, DIAG_LWORK,INFO)

     ENDIF


#elif defined(SINGLEPREC)
     CALL SSYEVD(JOBZ, UPLO, HDIM, UPEVECS, HDIM, UPEVALS, &
          DIAG_WORK, DIAG_LWORK, DIAG_IWORK, DIAG_LIWORK,INFO)
#endif

#ifdef DOUBLEPREC
     CALL DSYEVD(JOBZ, UPLO, HDIM, DOWNEVECS, HDIM, DOWNEVALS, &
          DIAG_WORK, DIAG_LWORK, DIAG_IWORK, DIAG_LIWORK,INFO)

     IF (INFO .NE. 0) THEN

        INFO = 0

        CALL DSYEV(JOBZ, UPLO, HDIM, DOWNEVECS, HDIM, DOWNEVALS, &
             DIAG_WORK, DIAG_LWORK,INFO)

     ENDIF


#elif defined(SINGLEPREC)
     CALL SSYEVD(JOBZ, UPLO, HDIM, DOWNEVECS, HDIM, DOWNEVALS, &
          DIAG_WORK, DIAG_LWORK, DIAG_IWORK, DIAG_LIWORK, INFO)
#endif

#endif

  ENDIF


  RETURN

END SUBROUTINE DIAGMYH
