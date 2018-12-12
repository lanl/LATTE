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

MODULE HOMOLUMO

  IMPLICIT NONE

CONTAINS

  SUBROUTINE HOMOLUMOGAP(ITERZ)

    USE MYPRECISION
    USE CONSTANTS_MOD
    USE PUREARRAY
    USE SPARSEARRAY

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ITERZ
    INTEGER :: I, J
    REAL(8) :: PRECOMP, X_A, X_B, Y_A, Y_B, HGAMMA

!!!!!!!!!!!!!!!!!! CALCULATE HOMO-LUMO ESTIMATES !!!!!!!!!!
    X_A = 0.0D0
    X_B = 1.0D0
    Y_A = 0.0D0
    Y_B = 0.0D0

    I = ITERZ
    HGAMMA = SIX - FOUR*SQRT(TWO)
    HGAMMA = HGAMMA*( ONE - HGAMMA )

    !write(*,*) "First I = ", ITERZ, " HGAMMA = ", HGAMMA, " VV = ", VV(I)

    DO WHILE (DBLE(VV(I)) .LT. HGAMMA)

       !write(*,*) "ITERZ = ", ITERZ, "I = ", I, "VV = ", VV(I)

       PRECOMP = SQRT( ONE - FOUR*VV(I))
       Y_A = HALF * ( ONE + PRECOMP )
       Y_B = HALF * ( ONE - PRECOMP )

       DO J = I-1,1,-1                !!! SHIFT -1 IN I
          IF (PP(J) .GT. 0) THEN
             Y_A = SQRT(Y_A)
             Y_B = SQRT(Y_B)
          ELSE
             Y_A = 1.0D0 - SQRT(ONE - Y_A)
             Y_B = 1.0D0 - SQRT(ONE - Y_B)
          ENDIF
       ENDDO

       X_A = MAX(X_A,Y_A)
       X_B = MIN(X_B,Y_B)

       I = I - 1
       IF (I .LT. 1) THEN
          WRITE(*,*) "HomoLumo i = ", I
       ENDIF
    ENDDO

    EHOMO = MAXEVAL - REAL(X_A)*(MAXEVAL - MINEVAL)
    ELUMO = MAXEVAL - REAL(X_B)*(MAXEVAL - MINEVAL)

    EGAP = ELUMO - EHOMO

    !  PRINT*, "EGAP = ", EGAP

    !WRITE(*,*) '### EHOMO = ', EHOMO, '   ELUMO = ', ELUMO
    !WRITE(*,*) '### MAXEVAL = ', MAXEVAL, '   MINEVAL = ', MINEVAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  END SUBROUTINE HOMOLUMOGAP

  SUBROUTINE SP2SEQUENCE()

    USE MYPRECISION
    USE CONSTANTS_MOD
    USE PUREARRAY
    USE SPARSEARRAY

    IMPLICIT NONE

    INTEGER :: IT
    REAL(8) :: EH, EL, ERR, SGM
    !REAL(LATTEPREC), PARAMETER :: ERRLIMIT=1E-20
    REAL(8), PARAMETER :: ERRLIMIT=1D-16
    REAL(8) :: DMAXEVAL, DMINEVAL, DEHOMO, DELUMO

    DMAXEVAL = DBLE(MAXEVAL)
    DMINEVAL = DBLE(MINEVAL)
    DEHOMO = DBLE(EHOMO)
    DELUMO = DBLE(ELUMO)

!!! GET SEQUENCE OF X^2 AND 2X-X^2 -> PP(1:NIT) = 1 FOR X^2 AND 0 FOR 2X-X^2

    !  WRITE(*,*) ' FORE NR_SP2_ITER = ',NR_SP2_ITER

    EH = (DMAXEVAL - DEHOMO)/(DMAXEVAL - DMINEVAL)
    EL = (DMAXEVAL - DELUMO)/(DMAXEVAL - DMINEVAL)

    ERR = ONE 

    IT = 0

    DO WHILE (ERR .GT. ERRLIMIT)

       IT = IT + 1

       IF ( (ABS(ONE - EH*EH) + ABS(EL*EL)) .LT. &
            (ABS(ONE - (TWO*EH - EH*EH) + ABS(TWO*EL - EL*EL)))) THEN

          PP(IT) = 1

          EH = EH*EH
          EL = EL*EL

       ELSE

          PP(IT) = 0

          EH = TWO*EH - EH*EH
          EL = TWO*EL - EL*EL

       ENDIF

       ERR = ABS(ONE - EH) + ABS(EL)

       !write(*,*) "SGM = ", SGM, " EH = ", EH, " EL = ", EL, " ERR = ", ERR

       IF (IT.GE.100) THEN
          ERR = ZERO
          WRITE(*,*) 'SP2SEQUENCE WARNING NOT CONVERGING IN SP2'
       ENDIF

    ENDDO

    NR_SP2_ITER = IT
    !  WRITE(*,*) ' #SIMPLE NR_SP2_ITER = ', NR_SP2_ITER

  END SUBROUTINE SP2SEQUENCE

END MODULE HOMOLUMO
