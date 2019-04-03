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

SUBROUTINE SP2GAP

  !
  ! This subroutine the HOMO-LUMO gap-based version of Niklasson's SP2
  ! method
  !


  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE PUREARRAY
  USE SPINARRAY
  USE NONOARRAY
  USE HOMOLUMO
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, ITER
  REAL(LATTEPREC) :: TRX, TRX2, TRXT, GERSHFACT
  IF (EXISTERROR) RETURN

  ! Estimate the largest and smallest eigenvalues

  CALL GERSHGORIN

  CALL SP2SEQUENCE

  IF (BASISTYPE .EQ. "ORTHO") THEN
     BO = -H/MAXMINUSMIN	 
  ELSE
     BO = -ORTHOH/MAXMINUSMIN
  ENDIF

  GERSHFACT =  MAXEVAL/MAXMINUSMIN

  TRX = ZERO  
  DO I = 1, HDIM     
     BO(I,I) = GERSHFACT + BO(I,I)
     TRX = TRX + BO(I,I)
  ENDDO

  ITER = 0

  DO WHILE (ITER .LT. NR_SP2_ITER )

     ITER = ITER + 1

#ifdef DOUBLEPREC     
     CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
          BO, HDIM, BO, HDIM, ZERO, X2, HDIM)
#elif defined(SINGLEPREC)
     CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, ONE, &
          BO, HDIM, BO, HDIM, ZERO, X2, HDIM)
#endif    

     TRX2 = ZERO
     DO I = 1, HDIM
        TRX2 = TRX2 + X2(I,I)
     ENDDO

     TRXT = ZERO

     DO I = 1, HDIM
        DO J = 1, HDIM

           TRXT = TRXT + (BO(J,I) - X2(J,I))*(BO(J,I) - X2(J,I))

        ENDDO
     ENDDO

     IF (PP(ITER) .EQ. 0) THEN

        TRX = TWO*TRX - TRX2

        BO = TWO*BO - X2

     ELSE

        TRX = TRX2

        BO = X2

     ENDIF

     VV(ITER) = SQRT(TRXT)

  ENDDO

  BO = TWO*BO

  CALL HOMOLUMOGAP(ITER)

  RETURN

END SUBROUTINE SP2GAP
