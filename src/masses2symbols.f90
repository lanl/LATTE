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

!> Soubroutine to get the Symbols out of the masses of the elements.
!! \param TYPES atom type index.
!! \param NTYPES Number of types.
!! \param MASSES_IN Masses for every type.
!! \param NATSIN Number of total atoms.
!! \param SYMBOLS Symbols for every atom.
!!
SUBROUTINE MASSES2SYMBOLS(TYPES,NTYPES,MASSES_IN,NATSIN,SYMBOLS)

  USE SETUPARRAY
  USE MYPRECISION
  USE CONSTANTS_MOD

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: NTYPES
  INTEGER, INTENT(IN) :: TYPES(NTYPES)
  INTEGER, INTENT(IN) :: NATSIN
  REAL(LATTEPREC), INTENT(IN) :: MASSES_IN(NTYPES)
  CHARACTER(LEN=2), INTENT(INOUT) :: SYMBOLS(NATSIN)
  INTEGER :: NZ = 103
  INTEGER, PARAMETER :: DP = LATTEPREC
  CHARACTER(2), ALLOCATABLE :: ELEMENT_SYMBOL(:), TYPE_SYMBOLS(:)
  CHARACTER(20) :: DUMMYC1, DUMMYC2, DUMMYC3, DUMMYC4
  REAL(DP), ALLOCATABLE :: ELEMENT_MASS(:)
  REAL(DP) :: DUMMYR(20)
  CHARACTER(20) :: DUMMYC(20)
  REAL(DP) :: DUMMYR1, DUMMYR2
  INTEGER :: I, J

  !> Read ptable
  !!
  OPEN(UNIT=14,STATUS="OLD", FILE=TRIM(PARAMPATH)//"/electrons.dat")

  ALLOCATE(ELEMENT_SYMBOL(NZ))
  ALLOCATE(ELEMENT_MASS(NZ))

  READ(14,*) DUMMYC1, NZ
  READ(14,*) (DUMMYC(J),J=1,13)
  DO I = 1,NZ
     READ(14,*) ELEMENT_SYMBOL(I), DUMMYC1, (DUMMYR(J),J=1,5), ELEMENT_MASS(I), (DUMMYR(J),J=6,10)
  ENDDO

  CLOSE(14)

  ALLOCATE(TYPE_SYMBOLS(NTYPES))
  TYPE_SYMBOLS=""

  DO I=1,NTYPES
     DO J=1,NZ
        IF(ABS(MASSES_IN(I) - ELEMENT_MASS(J)) < 0.01) THEN
           TYPE_SYMBOLS(I) = ELEMENT_SYMBOL(J)
           IF (EXISTERROR) RETURN
           EXIT
        ENDIF
     ENDDO
     IF(TYPE_SYMBOLS(I) == "")THEN
        WRITE(*,*)"ERROR: Mass of element",I,"cannot be identified."
        CALL ERRORS("masses2symbols","The mass of an element in the &
             & coordinates file cannot be identified. Please verify that &
             & masses in electrons.dat coincide with the masses in the input file")
     ENDIF
  ENDDO

  DO I=1,NATSIN
     SYMBOLS(I) = TYPE_SYMBOLS(TYPES(I))
  ENDDO

  DEALLOCATE(TYPE_SYMBOLS)
  DEALLOCATE(ELEMENT_SYMBOL)
  DEALLOCATE(ELEMENT_MASS)

END SUBROUTINE MASSES2SYMBOLS
