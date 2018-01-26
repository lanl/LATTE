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

MODULE MATRIXIO

  IMPLICIT NONE

CONTAINS

  !
  ! writehmatrix - Write out hamiltonian matrix to a file.
  !
  SUBROUTINE WRITEHMATRIX(HSIZE, MSIZE, HARRAY, NITER, PVEC)

    USE MYPRECISION

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: HSIZE, MSIZE, NITER
    INTEGER, INTENT(IN) :: PVEC(:)
    INTEGER :: I, J
    REAL(LATTEPREC), INTENT(IN) :: HARRAY(:,:)
    CHARACTER(LEN=100) :: FLNM

    WRITE(FLNM,'("hmatrix.in.dat")')

    OPEN (UNIT = 10, STATUS="UNKNOWN", FILE=FLNM)

    WRITE(10,*) HSIZE, MSIZE
    WRITE(10,*) NITER, (PVEC(I), I = 1,NITER)
    DO I = 1,HSIZE
       WRITE(10,*) (HARRAY(J,I), J = 1,HSIZE)
       !!WRITE(10,10) (HARRAY(J,I), J = 1,HSIZE)
    ENDDO
    !!10  FORMAT(100(E15.5,3X))

    CLOSE(10)

  END SUBROUTINE WRITEHMATRIX

  !
  ! writedmatrix - Write out density matrix to a file.
  !
  SUBROUTINE WRITEDMATRIX(HSIZE, DARRAY)

    USE MYPRECISION
    USE SETUPARRAY

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: HSIZE
    INTEGER :: I, J
    REAL(LATTEPREC), INTENT(IN) :: DARRAY(:,:)
    CHARACTER(LEN=100) :: FLNM

    WRITE(FLNM,'("dmatrix.out.dat")')

    OPEN (UNIT = 10, STATUS="UNKNOWN", FILE=FLNM)
    WRITE(10,*) HSIZE
    DO I = 1,HSIZE
       WRITE(10,10) (DARRAY(J,I), J = 1,HSIZE)
    ENDDO
10  FORMAT(100(E15.5,3X))
    CLOSE(10)

  END SUBROUTINE WRITEDMATRIX

  SUBROUTINE WRITEMTX(ITER, HSIZE, II, JJ, VAL)

    USE MYPRECISION

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: HSIZE, ITER
    INTEGER, INTENT(IN) :: II(:), JJ(:,:)
    INTEGER :: I, J, MSUM
    REAL(LATTEPREC), INTENT(IN) :: VAL(:,:)
    CHARACTER(LEN=100) :: FLNM
    CHARACTER(LEN=20) :: FFMT

    FFMT = '(A9,I2.2,A4)'

    WRITE(FLNM, FFMT) "spmatrix_", ITER, ".mtx"

    OPEN (UNIT = 10, STATUS="UNKNOWN", FILE=FLNM)
    WRITE(10,*) "%%MatrixMarket sparse coordinate real general"
    MSUM = 0
    DO I = 1, HSIZE
       MSUM = MSUM + II(I)
    ENDDO
    WRITE(10,*) HSIZE, HSIZE, MSUM
    DO I = 1,HSIZE
       DO J = 1,II(I)
          WRITE(10,*) I, JJ(J,I), VAL(J,I) 
       ENDDO
    ENDDO
    !10  FORMAT(100(E15.5,3X))
    CLOSE(10)

  END SUBROUTINE WRITEMTX

END MODULE MATRIXIO
