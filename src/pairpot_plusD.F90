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

SUBROUTINE PAIRPOTPLUSD

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE PPOTARRAY
  USE NEBLISTARRAY
  USE VIRIALARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, NEWJ, J, K, GOTI, GOTJ, PPIDI, PPIDJ
  INTEGER :: PBCI, PBCJ, PBCK
  REAL(LATTEPREC) :: JR1, JRCUT, R1, RCUT2, RCUT
  REAL(LATTEPREC) :: FORCE, DC(3), RIJ(3)
  REAL(LATTEPREC) :: MAGR2, MAGR
  REAL(LATTEPREC) :: R6, FTMP(3)
  REAL(LATTEPREC) :: MYR0, MYC6, FDAMP

  IF (EXISTERROR) RETURN

  VIRPLUSD = ZERO
  EPLUSD  = ZERO
  FPLUSD = ZERO
  
  DO I = 1, NATS
     
     FTMP = ZERO

     GOTI = 0
     DO K = 1, NOPPD
        IF (ATELE(I) .EQ. PPELE(K)) THEN
           GOTI = 1
           PPIDI = K
        ENDIF
     ENDDO
        
     
     ! Loop over all neighbors of I
     
     DO NEWJ = 1, TOTNEBPP(I)
        
        J = NEBPP(1, NEWJ, I)
        
        PBCI = NEBPP(2, NEWJ, I)
        PBCJ = NEBPP(3, NEWJ, I)
        PBCK = NEBPP(4, NEWJ, I)
        
        GOTJ = 0
        DO K = 1, NOPPD
           IF (ATELE(J) .EQ. PPELE(K)) THEN
              PPIDJ = K
              GOTJ = 1
           ENDIF
        ENDDO
        
        RIJ(1) = CR(1,J) + REAL(PBCI)*BOX(1,1) + REAL(PBCJ)*BOX(2,1) + &
             REAL(PBCK)*BOX(3,1) - CR(1,I)
        
        RIJ(2) = CR(2,J) + REAL(PBCI)*BOX(1,2) + REAL(PBCJ)*BOX(2,2) + &
             REAL(PBCK)*BOX(3,2) - CR(2,I)
        
        RIJ(3) = CR(3,J) + REAL(PBCI)*BOX(1,3) + REAL(PBCJ)*BOX(2,3) + &
             REAL(PBCK)*BOX(3,3) - CR(3,I)
        
        MAGR2 = RIJ(1)*RIJ(1) + RIJ(2)*RIJ(2) + RIJ(3)*RIJ(3)
        
        IF (MAGR2 .LE. PLUSDCUT*PLUSDCUT) THEN
           
           MAGR = SQRT(MAGR2) 
           
           ! Direction cosines
           
           DC = RIJ/MAGR
           
           MYC6 = SQRT(C6(PPIDI) * C6(PPIDJ))
           MYR0 = RZERO(PPIDI) + RZERO(PPIDJ)
!           MYR0 = HALF*( RZERO(PPIDI) + RZERO(PPIDJ))
           
           MYC6 = MYC6*PLUSDS6

           FDAMP = ( ONE + EXP(-PLUSDGAMMA * ((MAGR/MYR0) - ONE)))
           
           R6 = MAGR2*MAGR2*MAGR2
           
           EPLUSD = EPLUSD - (MYC6/R6)/FDAMP * REAL(GOTI*GOTJ)
           
           FTMP = DC*(-(MYC6/R6)*&
                (PLUSDGAMMA/MYR0)*EXP(-PLUSDGAMMA*((MAGR/MYR0) - ONE))/(FDAMP*FDAMP) & 
                 + (6.0d0*MYC6/(R6*MAGR))/FDAMP)
           
           FTMP = FTMP * REAL(GOTI*GOTJ)
           
           VIRPLUSD(1) = VIRPLUSD(1) + RIJ(1) * FTMP(1)
           VIRPLUSD(2) = VIRPLUSD(2) + RIJ(2) * FTMP(2)
           VIRPLUSD(3) = VIRPLUSD(3) + RIJ(3) * FTMP(3)
           VIRPLUSD(4) = VIRPLUSD(4) + RIJ(1) * FTMP(2)
           VIRPLUSD(5) = VIRPLUSD(5) + RIJ(2) * FTMP(3)
           VIRPLUSD(6) = VIRPLUSD(6) + RIJ(3) * FTMP(1)
           
           FPLUSD(1,I) = FPLUSD(1,I) + FTMP(1)
           FPLUSD(2,I) = FPLUSD(2,I) + FTMP(2)
           FPLUSD(3,I) = FPLUSD(3,I) + FTMP(3)
           
           
        ENDIF
     ENDDO
  ENDDO
    
  EPLUSD = HALF*EPLUSD
  VIRPLUSD = HALF*VIRPLUSD

  RETURN

END SUBROUTINE PAIRPOTPLUSD

