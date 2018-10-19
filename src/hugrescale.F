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


SUBROUTINE HUGRESCALE

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MDARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, K
  REAL(LATTEPREC) :: PREFACTOR, PREPREFACT
  REAL(LATTEPREC) :: CHI, TMPTEMP, MYTEMP
  REAL(LATTEPREC) :: MAXCHI = 1.05D0, MAXDBOX = 0.001D0
  REAL(LATTEPREC) :: PSCALE = 100.0D0
  REAL(LATTEPREC) :: TMPPRESS, MYPRESS, DBOX, DTEMP
  REAL(LATTEPREC) :: TMPPRESSX, TMPPRESSY, TMPPRESSZ, TMPENERGY, AVEENERGY
  REAL(LATTEPREC) :: DBOXX, DBOXY, DBOXZ, MYPRESSX, MYPRESSY, MYPRESSZ
  REAL(LATTEPREC) :: AVEVOL, TMPVOL
  REAL(LATTEPREC), PARAMETER :: MAXDT = 10.0
  IF (EXISTERROR) RETURN

  TMPTEMP = ZERO
  ! Isotropic
  TMPPRESS = ZERO
  ! non-isotropic
  TMPPRESSX = ZERO
  TMPPRESSY = ZERO
  TMPPRESSZ = ZERO

  ! Energy and volume too

  TMPENERGY = ZERO
  TMPVOL = ZERO


  DO I = 1, AVEPER/WRTFREQ

     TMPENERGY = TMPENERGY + EHIST(I)
     TMPVOL = TMPVOL + VHIST(I)
     TMPTEMP = TMPTEMP + THIST(I)
     IF (NPTTYPE .EQ. "ISO") THEN
        TMPPRESS = TMPPRESS + PHIST(I)
     ELSE
        TMPPRESSX = TMPPRESSX + PHISTX(I)
        TMPPRESSY = TMPPRESSY + PHISTY(I)
        TMPPRESSZ = TMPPRESSZ + PHISTZ(I)
     ENDIF

  ENDDO

  MYTEMP = TMPTEMP/REAL(AVEPER/WRTFREQ)

  MYPRESS = TMPPRESS/REAL(AVEPER/WRTFREQ)

  MYPRESSX = TMPPRESSX/REAL(AVEPER/WRTFREQ)
  MYPRESSY = TMPPRESSY/REAL(AVEPER/WRTFREQ)
  MYPRESSZ = TMPPRESSZ/REAL(AVEPER/WRTFREQ)

  AVEENERGY = TMPENERGY/REAL(AVEPER/WRTFREQ)

  AVEVOL = TMPVOL/REAL(AVEPER/WRTFREQ)



  ! Change box dimensions depending on the average pressure

  IF (NPTTYPE .EQ. "ISO") THEN

     DBOX = MIN(ABS(MYPRESS - PTARGET)*MAXDBOX, MAXDBOX)

     DBOX = SIGN(DBOX, MYPRESS - PTARGET)

     BOX = BOX * (ONE + DBOX)

     CR = CR * (ONE + DBOX)

  ELSE ! Allow the three vectors to change length independently

     DBOXX = MIN(ABS(MYPRESSX - PTARGET)*MAXDBOX, MAXDBOX)
     DBOXY = MIN(ABS(MYPRESSY - PTARGET)*MAXDBOX, MAXDBOX)
     DBOXZ = MIN(ABS(MYPRESSZ - PTARGET)*MAXDBOX, MAXDBOX)

     DBOXX = SIGN(DBOXX, MYPRESSX - PTARGET)
     DBOXY = SIGN(DBOXY, MYPRESSY - PTARGET)
     DBOXZ = SIGN(DBOXZ, MYPRESSZ - PTARGET)

     BOX(1,1) = BOX(1,1) * (ONE + DBOXX) 
     BOX(2,2) = BOX(2,2) * (ONE + DBOXY) 
     BOX(3,3) = BOX(3,3) * (ONE + DBOXZ) 

     DO I = 1, NATS

        CR(1,I) = CR(1,I) * (ONE + DBOXX)
        CR(2,I) = CR(2,I) * (ONE + DBOXY)
        CR(3,I) = CR(3,I) * (ONE + DBOXZ)

     ENDDO

  ENDIF

  ! SOLVING HG = (E - E0) - 0.5*(P + P0)*(V0 - V)

  HG = (AVEENERGY - E0) - HALF*(MYPRESS + P0)*(V0 - AVEVOL)/TOGPA

  DTEMP = 0.025D0*TTARGET

  DTEMP = MIN(DTEMP, MAXDT)

  IF (HG .LT. ZERO) THEN
     TTARGET = TTARGET + DTEMP
  ELSE
     TTARGET = TTARGET - DTEMP
  ENDIF

  CHI = SQRT(TTARGET/MYTEMP)

  V = V * CHI

  RETURN

END SUBROUTINE HUGRESCALE

