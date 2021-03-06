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

SUBROUTINE GETPRESSURE

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE VIRIALARRAY
  USE MDARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J
  REAL(LATTEPREC) :: MYMASS, TPRESS
  REAL(LATTEPREC) :: MYP(5)
  IF (EXISTERROR) RETURN
  KETEN = ZERO

  !  SYSVOL = BOXDIMS(1)*BOXDIMS(2)*BOXDIMS(3)

  SYSVOL = ABS(BOX(1,1)*(BOX(2,2)*BOX(3,3) - BOX(3,2)*BOX(2,3)) + &
       BOX(1,2)*(BOX(2,1)*BOX(3,3) - BOX(3,1)*BOX(2,3)) + &
       BOX(1,3)*(BOX(2,1)*BOX(3,2) - BOX(3,1)*BOX(2,2)))


  VIRIAL = VIRBOND

  IF (PPOTON .NE. 0) VIRIAL = VIRIAL + VIRPAIR

  IF (ELECTRO .EQ. 1) VIRIAL = VIRIAL + VIRCOUL

  IF (BASISTYPE .EQ. "NONORTHO") THEN

     IF (SPONLY .EQ. 0) THEN

        VIRIAL = VIRIAL - VIRPUL
        
        IF (ELECTRO .EQ. 0) VIRIAL = VIRIAL + VIRSLCN
        IF (ELECTRO .EQ. 1) VIRIAL = VIRIAL + VIRSCOUL
        IF (SPINON .EQ. 1) VIRIAL = VIRIAL + VIRSSPIN

     ELSE 

        VIRIAL = VIRIAL + VIRPUL

     ENDIF

  ENDIF

  IF (PLUSDON .EQ. 1) VIRIAL = VIRIAL + VIRPLUSD

  IF ( MDON .EQ. 1 ) THEN

     DO I = 1, NATS

        MYMASS = MASS(ELEMPOINTER(I))

        KETEN(1) = KETEN(1) + MYMASS*V(1,I)*V(1,I)
        KETEN(2) = KETEN(2) + MYMASS*V(2,I)*V(2,I)
        KETEN(3) = KETEN(3) + MYMASS*V(3,I)*V(3,I)
        KETEN(4) = KETEN(4) + MYMASS*V(1,I)*V(2,I)
        KETEN(5) = KETEN(5) + MYMASS*V(2,I)*V(3,I)
        KETEN(6) = KETEN(6) + MYMASS*V(3,I)*V(1,I)

     ENDDO

     CALL GETKE 

     TPRESS = (REAL(NATS)/SYSVOL)*TEMPERATURE/KE2T

  ENDIF

  ! Minus the virial sum because we have Rij = Rj - Ri. See Allen 
  ! and Tildesley.

  STRTEN(1) = ( -VIRIAL(1) + KETEN(1)/F2V ) / SYSVOL
  STRTEN(2) = ( -VIRIAL(2) + KETEN(2)/F2V ) / SYSVOL
  STRTEN(3) = ( -VIRIAL(3) + KETEN(3)/F2V ) / SYSVOL
  STRTEN(4) = ( -VIRIAL(4) + KETEN(4)/F2V ) / SYSVOL
  STRTEN(5) = ( -VIRIAL(5) + KETEN(5)/F2V ) / SYSVOL
  STRTEN(6) = ( -VIRIAL(6) + KETEN(6)/F2V ) / SYSVOL


  !
  ! That's right folkS: we're going to be reporting pressure in 
  ! GPa since that's what some of us materials science/shock physics
  ! people like to deal with. STRTEN is still in eV/A^3
  !

  STRTEN = STRTEN * TOGPA

  PRESSURE = (STRTEN(1) + STRTEN(2) + STRTEN(3))/THREE

  RETURN

END SUBROUTINE GETPRESSURE
