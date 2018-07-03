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

SUBROUTINE PLOTUNIV

  USE CONSTANTS_MOD
  USE UNIVARRAY
  USE PPOTARRAY
  USE SETUPARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: IGS, I, J
  INTEGER, PARAMETER :: NOSAMPLES = 1001
  REAL(LATTEPREC), ALLOCATABLE :: INTVALUE(:), DINTVALUE(:)
  REAL(LATTEPREC) :: MAGR, TMP, DTMP, POLYNOM, DPOLYNOM
  REAL(LATTEPREC) :: RMOD, RMINUSR1, MAXCUT
  IF (EXISTERROR) RETURN

  ALLOCATE( INTVALUE(NOINT), DINTVALUE(NOINT) )

  OPEN(UNIT=40, STATUS="UNKNOWN", FILE="H_scaling.dat")
  OPEN(UNIT=41, STATUS="UNKNOWN", FILE="dH_scaling.dat")

  MAXCUT = ZERO
  DO I = 1, NOINT
     IF (HCUT(I) .GT. MAXCUT) MAXCUT = HCUT(I)
  ENDDO

  IF (SCLTYPE .EQ. "EXP") THEN

     DO I = 1, NOSAMPLES
        
        MAGR = HALF + (MAXCUT - HALF)*REAL(I-1)/REAL(NOSAMPLES-1)
        
        DO J = 1, NOINT
           
           IF (MAGR .LE. BOND(7,J)) THEN
              
              RMOD = MAGR - BOND(6,J)
              
              POLYNOM = RMOD*(BOND(2,J) + RMOD*(BOND(3,J) + &
                   RMOD*(BOND(4,J) + BOND(5,J)*RMOD)))
              
              DPOLYNOM =  BOND(2,J) + RMOD*(TWO*BOND(3,J) + &
                   RMOD*(THREE*BOND(4,J) + FOUR*BOND(5,J)*RMOD))
              
              TMP = EXP(POLYNOM)
              
              DTMP = DPOLYNOM*TMP
              
           ELSEIF (MAGR .GT. BOND(7,J) .AND. MAGR .LE. BOND(8,J)) THEN
              
              RMINUSR1 = MAGR - BOND(7,J)
              
              TMP = BOND(9,J) + RMINUSR1*(BOND(10,J) + &
                   RMINUSR1*(BOND(11,J)+ RMINUSR1*(BOND(12,J) + &
                   RMINUSR1*(BOND(13,J) + RMINUSR1*BOND(14,J)))))
              
              DTMP = BOND(10,J) + RMINUSR1*(TWO*BOND(11,J) + &
                   RMINUSR1*(THREE*BOND(12,J) + RMINUSR1*(FOUR*BOND(13,J) + &
                   RMINUSR1*FIVE*BOND(14,J))))
              
           ELSE
              
              TMP = ZERO
              DTMP = ZERO
              
           ENDIF
           
           INTVALUE(J) = BOND(1,J)*TMP
           DINTVALUE(J) = -BOND(1,J)*DTMP
        
        ENDDO

        WRITE(40,10) MAGR, (INTVALUE(J), J = 1, NOINT)
        WRITE(41,10) MAGR, (DINTVALUE(J), J = 1, NOINT)
        
     ENDDO
     
     CLOSE(40)
     CLOSE(41)

     IF (BASISTYPE .EQ. "NONORTHO") THEN
        
        MAXCUT = ZERO
        DO I = 1, NOINT
           IF (OVERL(8,I) .GT.MAXCUT) MAXCUT = OVERL(8,I)
        ENDDO
        
        OPEN(UNIT=40, STATUS="UNKNOWN", FILE="S_scaling.dat")
        
        DO I = 1, NOSAMPLES
           
           MAGR = HALF + (MAXCUT - HALF)*REAL(I-1)/REAL(NOSAMPLES-1)
           
           DO J = 1, NOINT
              
              IF (MAGR .LE. OVERL(7,J)) THEN
                 
                 RMOD = MAGR - OVERL(6,J)
                 POLYNOM = RMOD*(OVERL(2,J) + RMOD*(OVERL(3,J) + &
                      RMOD*(OVERL(4,J) + OVERL(5,J)*RMOD)))
                 
                 TMP = EXP(POLYNOM)
                 
              ELSEIF (MAGR .GT. OVERL(7,J) .AND. MAGR .LE. OVERL(8,J)) THEN
                 
                 RMINUSR1 = MAGR - OVERL(7,J)
                 
                 TMP = OVERL(9,J) + RMINUSR1*(OVERL(10,J) + &
                      RMINUSR1*(OVERL(11,J)+ RMINUSR1*(OVERL(12,J) + &
                      RMINUSR1*(OVERL(13,J) + RMINUSR1*OVERL(14,J)))))
                 
              ELSE 
                 
                 TMP = ZERO
                 
              ENDIF
              
              INTVALUE(J) = OVERL(1,J)*TMP
              
           ENDDO
           
           WRITE(40,10) MAGR, (INTVALUE(J), J = 1, NOINT)
           
        ENDDO
        
        CLOSE(40)
        
     ENDIF
     
10   FORMAT(F12.6, 1X, 500(G12.6, 1X))
     DEALLOCATE(INTVALUE, DINTVALUE)

  ELSEIF (SCLTYPE .EQ. "TABLE") THEN
     
     OPEN(UNIT=40, STATUS="UNKNOWN", FILE="H_scaling.dat")
     OPEN(UNIT=41, STATUS="UNKNOWN", FILE="S_scaling.dat")
     
     DO I = 1, LENTABINT(1)
        
        WRITE(40,20) TABR(I,1), (TABH(I,J), J = 1, NOINT)
        WRITE(41,20) TABR(I,1), (TABS(I,J), J = 1, NOINT)
        
     ENDDO
     
     CLOSE(40)
     CLOSE(41)
     
  ENDIF
  
20 FORMAT(100F12.6)

  RETURN

END SUBROUTINE PLOTUNIV
