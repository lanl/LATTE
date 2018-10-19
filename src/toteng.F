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

SUBROUTINE TOTENG

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE SPINARRAY
  USE NONOARRAY
  USE KSPACEARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, K
  COMPLEX(LATTEPREC) :: ZTRRHOH
  IF (EXISTERROR) RETURN

  TRRHOH = ZERO
  ZTRRHOH = CMPLX(ZERO,ZERO)

  !  IF (ELECTRO .EQ. 1) THEN

  IF (SPINON .EQ. 0) THEN

     IF (KON .EQ. 0) THEN

        DO I = 1, HDIM
           DO J = 1, HDIM

              TRRHOH = TRRHOH + BO(J,I)*H(J,I)

           ENDDO

           TRRHOH = TRRHOH - BOZERO(I)*H(I,I)

        ENDDO

     ELSE

        DO K = 1, NKTOT

           DO I = 1, HDIM
              DO J = 1, HDIM

                 ZTRRHOH = ZTRRHOH + KBO(J,I,K)*(HK(I,J,K))

              ENDDO

              ZTRRHOH = ZTRRHOH - CMPLX(BOZERO(I))*HK(I,I,K)

           ENDDO

        ENDDO

        TRRHOH = REAL(ZTRRHOH)/REAL(NKTOT)

     ENDIF

  ELSE

     !
     ! This is everything: tr(rhoup - rhoupzero)*Hup + 
     ! tr(rhodown - rhodownzero)*Hdown
     !
     ! Hup = H(Slater-Koster) + H(electrostatic) + H(spin)
     ! Hdown = H(Slater-Koster) + H(electrostatic) - H(spin)
     !
     ! Thus, we're calculating Covalent (from SK) + electrostatic (from H1) +
     ! spin (from H2) and we just need to add on the entropy and pairwise 
     ! bits to get the total energy
     !

     DO I = 1, HDIM
        DO J = 1, HDIM

           TRRHOH = TRRHOH + (RHOUP(J,I) + RHODOWN(J,I))*H(J,I)

        ENDDO

        TRRHOH = TRRHOH - (RHOUPZERO(I) + RHODOWNZERO(I))*H(I,I)

     ENDDO

  ENDIF

  !  ELSE

  !     IF (KON .EQ. 0) THEN

  !        DO I = 1, HDIM
  !           DO J = 1, HDIM!
  !
  !              TRRHOH = TRRHOH + BO(J,I)*H(J,I)
  !
  !           ENDDO
  !        ENDDO

  !     ELSE

  !        DO K = 1, NKTOT

  !           DO I = 1, HDIM
  !              DO J = 1, HDIM

  !                 ZTRRHOH = ZTRRHOH + KBO(J,I,K)*HK(I,J,K)

  !              ENDDO
  !           ENDDO

  !        ENDDO

  !        TRRHOH = REAL(ZTRRHOH)/REAL(NKTOT)
  !     
  !     ENDIF

  !  ENDIF

  ! Check for something bad happening


  RETURN

END SUBROUTINE TOTENG
