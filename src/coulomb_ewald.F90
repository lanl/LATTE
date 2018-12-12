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

  INTEGER :: I, J, L, M, N, IK
  INTEGER :: LMIN, MMIN, NMIN
  REAL(LATTEPREC) :: PREFACTOR, KEPREF, K2, K(3), DOT
  REAL(LATTEPREC) :: COSSUM, SINSUM, COSSUM2, SINSUM2
  REAL(LATTEPREC) :: FLL, FLM, FLN
  REAL(LATTEPREC) :: FORCE
  REAL(LATTEPREC) :: PREVIR, CORRFACT
  REAL(LATTEPREC) :: L11, L12, L13, M21, M22, M23

  IF (EXISTERROR) RETURN

  ! $OMP PARALLEL DO DEFAULT(NONE) &
  ! $OMP SHARED(NATS,DELTAQ,K1_LIST,K2_LIST,K3_LIST,KSQ_LIST) &
  ! $OMP SHARED(KCUTOFF2,EIGHTPI,FOURCALPHA2,COULVOL,KECONST,CR) &
  ! $OMP PRIVATE(I,IK,DOT,FORCE,SINLIST,COSLIST,COSSUM,SINSUM) &
  ! $OMP PRIVATE(K2,PREFACTOR,PREVIR,COSSUM2,SINSUM2,KEPREF) &
  ! $OMP REDUCTION(+:FCOUL,VIRCOUL,COULOMBV)
  DO IK = 1, NK

     K2 = KSQ_LIST(IK)

     IF (K2 .LE. KCUTOFF2) THEN

        PREFACTOR = EIGHTPI*EXP( -K2/(FOURCALPHA2) ) / &
             (COULVOL*K2)

        PREVIR = (TWO/K2) + (TWO/(FOURCALPHA2))

        COSSUM = ZERO
        SINSUM = ZERO

        ! Doing the sin and cos sums

        DO I = 1, NATS

           DOT = K1_LIST(IK)*CR(1,I) + K2_LIST(IK)*CR(2,I) + K3_LIST(IK)*CR(3,I)

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

           FCOUL(1,I) = FCOUL(1,I) + FORCE*K1_LIST(IK)
           FCOUL(2,I) = FCOUL(2,I) + FORCE*K2_LIST(IK)
           FCOUL(3,I) = FCOUL(3,I) + FORCE*K3_LIST(IK)

        ENDDO

        KEPREF = -KEPREF * (COSSUM2 + SINSUM2)/TWO

        VIRCOUL(1) = VIRCOUL(1) + KEPREF * &
             (ONE - PREVIR*K1_LIST(IK)*K1_LIST(IK))
        VIRCOUL(2) = VIRCOUL(2) + KEPREF * &
             (ONE - PREVIR*K2_LIST(IK)*K2_LIST(IK))
        VIRCOUL(3) = VIRCOUL(3) + KEPREF * &
             (ONE - PREVIR*K3_LIST(IK)*K3_LIST(IK))
        VIRCOUL(4) = VIRCOUL(4) - KEPREF*PREVIR* &
             K1_LIST(IK)*K2_LIST(IK)
        VIRCOUL(5) = VIRCOUL(5) - KEPREF*PREVIR* &
             K2_LIST(IK)*K3_LIST(IK)
        VIRCOUL(6) = VIRCOUL(6) - KEPREF*PREVIR* &
             K3_LIST(IK)*K1_LIST(IK)

     ENDIF

  END DO
  ! $OMP END PARALLEL DO

  ! Point self energy

  CORRFACT = TWO*KECONST*CALPHA/SQRTPI

  COULOMBV = COULOMBV - CORRFACT*DELTAQ

  RETURN

END SUBROUTINE COULOMBEWALD
