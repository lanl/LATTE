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
  INTEGER :: NATSIN

  REAL(LATTEPREC), INTENT(IN) :: MASSES_IN(NTYPES)
  CHARACTER(LEN=2), INTENT(INOUT) :: SYMBOLS(NATSIN)

  INTEGER, PARAMETER :: NZ = 103
  INTEGER, PARAMETER :: DP = LATTEPREC
  CHARACTER(2) :: ELEMENT_SYMBOL(NZ), DUMMY
  REAL(DP) :: ELEMENT_MASS(NZ)
  INTEGER :: I, J
  CHARACTER(LEN=2), ALLOCATABLE :: TYPE_SYMBOLS(:)

  !> Read ptable
  !!
  OPEN(UNIT=14,STATUS="OLD", FILE=trim(PARAMPATH)//"/ptable.dat")

  READ(14,*) DUMMY
  DO I = 1,NZ
    READ(14,*) ELEMENT_SYMBOL(I), ELEMENT_MASS(I)
  ENDDO

  CLOSE(14)

  ALLOCATE(TYPE_SYMBOLS(NTYPES))

  DO I=1,NTYPES
    DO J=1,NZ
      IF(ABS(MASSES_IN(I) - ELEMENT_MASS(J)) < 0.001) THEN
        TYPE_SYMBOLS(I) = ELEMENT_SYMBOL(J)
        EXIT
      ENDIF
    ENDDO
  ENDDO

  DO I=1,NATSIN
    SYMBOLS(I) = TYPE_SYMBOLS(TYPES(I))
  ENDDO

  DEALLOCATE(TYPE_SYMBOLS)

END SUBROUTINE MASSES2SYMBOLS
