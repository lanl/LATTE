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

SUBROUTINE GERSHGORIN

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE SPINARRAY
  USE NONOARRAY
  USE SPARSEARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J
  REAL(LATTEPREC) :: RADIUS, DOWNRAD, UPRAD, ABSHAM
  IF (EXISTERROR) RETURN

  MINEVAL = 100000000000.0
  MAXEVAL = -100000000000.0

  NNZ = 0

  IF (SPINON .EQ. 0) THEN

     IF (BASISTYPE .EQ. "ORTHO") THEN

        ! No spin, orthogonal

        DO I = 1, HDIM

           RADIUS = ZERO

           DO J = 1, HDIM

              ABSHAM = ABS( H(J,I) )
              RADIUS = RADIUS + ABSHAM

              IF ( ABSHAM .GT. HTHRESH ) NNZ = NNZ + 1

           ENDDO

           RADIUS = RADIUS - ABS(H(I,I))

           MAXEVAL = MAX(MAXEVAL, H(I,I) + RADIUS)
           MINEVAL = MIN(MINEVAL, H(I,I) - RADIUS)

        ENDDO

     ELSE

        ! No spin, non-orthogonal

        DO I = 1, HDIM

           RADIUS = ZERO

           DO J = 1, HDIM

              ABSHAM = ABS( ORTHOH(J,I) )
              RADIUS = RADIUS + ABSHAM

              IF ( ABSHAM .GT. HTHRESH ) NNZ = NNZ + 1

           ENDDO

           RADIUS = RADIUS - ABS(ORTHOH(I,I))

           MAXEVAL = MAX(MAXEVAL, ORTHOH(I,I) + RADIUS)
           MINEVAL = MIN(MINEVAL, ORTHOH(I,I) - RADIUS)

        ENDDO

     ENDIF

  ELSE

     IF (BASISTYPE .EQ. "ORTHO") THEN

        ! Spin polarized, orthogonal

        DO I = 1, HDIM

           UPRAD = ZERO
           DOWNRAD = ZERO

           DO J = 1, HDIM

              UPRAD = UPRAD + ABS( HUP(J,I) )
              DOWNRAD = DOWNRAD + ABS( HDOWN(J,I) )

           ENDDO

           UPRAD = UPRAD - ABS(HUP(I,I))
           DOWNRAD = DOWNRAD - ABS(HDOWN(I,I))

           MAXEVAL = MAX(MAXEVAL, HUP(I,I) + UPRAD, &
                HDOWN(I,I) + DOWNRAD)
           MINEVAL = MIN(MINEVAL, HUP(I,I) - UPRAD, &
                HDOWN(I,I) - DOWNRAD)

        ENDDO

     ELSE

        ! Spin-polarized, non-orthogonal

        DO I = 1, HDIM

           UPRAD = ZERO
           DOWNRAD = ZERO

           DO J = 1, HDIM

              UPRAD = UPRAD + ABS( ORTHOHUP(J,I) )
              DOWNRAD = DOWNRAD + ABS( ORTHOHDOWN(J,I) )

           ENDDO

           UPRAD = UPRAD - ABS(ORTHOHUP(I,I))
           DOWNRAD = DOWNRAD - ABS(ORTHOHDOWN(I,I))

           MAXEVAL = MAX(MAXEVAL, ORTHOHUP(I,I) + UPRAD, &
                ORTHOHDOWN(I,I) + DOWNRAD)
           MINEVAL = MIN(MINEVAL, ORTHOHUP(I,I) - UPRAD, &
                ORTHOHDOWN(I,I) - DOWNRAD)

        ENDDO

     ENDIF

  ENDIF

  MAXMINUSMIN = MAXEVAL - MINEVAL

  RETURN

END SUBROUTINE GERSHGORIN
