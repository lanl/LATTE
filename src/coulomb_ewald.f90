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

SUBROUTINE COULOMBEWALD

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE COULOMBARRAY
  USE VIRIALARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, L, M, N
  INTEGER :: LMIN, MMIN, NMIN
  REAL(LATTEPREC) :: PREFACTOR, KEPREF, K2, K(3), DOT
  REAL(LATTEPREC) :: COSSUM, SINSUM, COSSUM2, SINSUM2
  REAL(LATTEPREC) :: FLL, FLM, FLN
  REAL(LATTEPREC) :: FORCE
  REAL(LATTEPREC) :: PREVIR, CORRFACT
  REAL(LATTEPREC) :: L11, L12, L13, M21, M22, M23
  IF (EXISTERROR) RETURN

  LMIN = 0

  ! LMAX, MMAX, NMAX computed during initialization

  DO L = LMIN, LMAX

     IF (L .EQ. 0) THEN
        MMIN = 0
     ELSE
        MMIN = -MMAX
     ENDIF

     L11 = REAL(L)*RECIPVECS(1,1)
     L12 = REAL(L)*RECIPVECS(1,2)
     L13 = REAL(L)*RECIPVECS(1,3)

     DO M = MMIN, MMAX

        NMIN = -NMAX

        IF (L .EQ. 0 .AND. M .EQ. 0) NMIN = 1

        M21 = L11 + REAL(M)*RECIPVECS(2,1)
        M22 = L12 + REAL(M)*RECIPVECS(2,2)
        M23 = L13 + REAL(M)*RECIPVECS(2,3)

        DO N = NMIN, NMAX

           K(1) = M21 + REAL(N)*RECIPVECS(3,1)
           K(2) = M22 + REAL(N)*RECIPVECS(3,2)
           K(3) = M23 + REAL(N)*RECIPVECS(3,3)

           K2 = K(1)*K(1) + K(2)*K(2) + K(3)*K(3)

           IF (K2 .LE. KCUTOFF2) THEN

              PREFACTOR = EIGHTPI*EXP( -K2/(FOURCALPHA2) ) / &
                   (COULVOL*K2)

              PREVIR = (TWO/K2) + (TWO/(FOURCALPHA2))

              COSSUM = ZERO
              SINSUM = ZERO

              ! Doing the sin and cos sums

              DO I = 1, NATS

                 DOT = K(1)*CR(1,I) + K(2)*CR(2,I) + K(3)*CR(3,I)

                 ! We re-use these in the next loop...

                 SINLIST(I) = SIN(DOT)
                 COSLIST(I) = COS(DOT)

                 COSSUM = COSSUM + DELTAQ(I)*COSLIST(I)
                 SINSUM = SINSUM + DELTAQ(I)*SINLIST(I)

              ENDDO

              COSSUM2 = COSSUM*COSSUM
              SINSUM2 = SINSUM*SINSUM

              ! Add up energy and force contributions

              KEPREF = KECONST*PREFACTOR

              DO I = 1, NATS

                 COULOMBV(I) = COULOMBV(I) +  &
                      KEPREF*(COSLIST(I)*COSSUM + SINLIST(I)*SINSUM)

                 FORCE = KEPREF * DELTAQ(I) * &
                      (SINLIST(I)*COSSUM - COSLIST(I)*SINSUM)

                 FCOUL(1,I) = FCOUL(1,I) + FORCE*K(1)
                 FCOUL(2,I) = FCOUL(2,I) + FORCE*K(2)
                 FCOUL(3,I) = FCOUL(3,I) + FORCE*K(3)

              ENDDO

              KEPREF = -KEPREF * (COSSUM2 + SINSUM2)/TWO

              VIRCOUL(1) = VIRCOUL(1) + KEPREF * &
                   (ONE - PREVIR*K(1)*K(1))
              VIRCOUL(2) = VIRCOUL(2) + KEPREF * &
                   (ONE - PREVIR*K(2)*K(2))
              VIRCOUL(3) = VIRCOUL(3) + KEPREF * &
                   (ONE - PREVIR*K(3)*K(3))
              VIRCOUL(4) = VIRCOUL(4) - KEPREF*PREVIR* &
                   K(1)*K(2)
              VIRCOUL(5) = VIRCOUL(5) - KEPREF*PREVIR* &
                   K(2)*K(3)
              VIRCOUL(6) = VIRCOUL(6) - KEPREF*PREVIR* &
                   K(3)*K(1)

           ENDIF

        ENDDO
     ENDDO
  ENDDO

  ! Point self energy

  CORRFACT = TWO*KECONST*CALPHA/SQRTPI

  COULOMBV = COULOMBV - CORRFACT*DELTAQ

  RETURN

END SUBROUTINE COULOMBEWALD

