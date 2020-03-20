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

!SUBROUTINE ADDDFTBU_INIT  ! ANDERS CHANGE 
!
!  USE CONSTANTS_MOD
!  USE SETUPARRAY
!  USE COULOMBARRAY
!  USE NONOARRAY
!  USE KSPACEARRAY
!  USE MYPRECISION
!  USE DMARRAY
!
!  IMPLICIT NONE
!
!  INTEGER :: I, J, K
!  INTEGER :: INDEX, NUMORB
!  REAL(LATTEPREC) :: ES, EP, ED, EF
!  REAL(LATTEPREC) :: HMOD, X(HDIM,HDIM), Y(HDIM,HDIM), HU(HDIM,HDIM), H_U(HDIM,HDIM)
!  REAL(LATTEPREC) :: DS(HDIM,HDIM)
!  IF (EXISTERROR) RETURN
!
!
!  do J = 1, HDIM
!  do I = 1, HDIM
!     SU(I,J) = SMAT(I,J)*DFTB_U(J)
!  enddo
!  enddo
!  call DGEMM('T','N',HDIM,HDIM,HDIM,ONE,BO,HDIM,SMAT,HDIM,ZERO,DS,HDIM)  !! Should be PNO instead of BO for propagation in XL-BOMD
!  call DGEMM('T','N',HDIM,HDIM,HDIM,ONE,DS,HDIM,SU,HDIM,ZERO,X,HDIM)
!
!  call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,SU,HDIM,DS,HDIM,ZERO,Y,HDIM)
!  HU = 0.125D0*(2.D0*SU-X-Y) ! since BO = 2*D, otherwise (1/4)*(...)
!  do J = 1,HDIM
!  do I = 1,HDIM
!     H_U(I,J) = HU(I,J) + HU(J,I)
!     !write(6,*) J,I,H_U(I,J),H(I,J)
!  enddo
!  enddo
!  H = H + H_U
!
!  IF (KON .EQ. 0) THEN
!     IF ( BASISTYPE .EQ. "ORTHO") THEN
!          WRITE(*,*) ' DM HUBBARD U ONLY FOR NONORTHOGONAL REPRESENTATIONS '
!     ENDIF
!  ELSE
!     WRITE(*,*) ' DM HUBBARD U ONLY FOR KON = 0 '
!  ENDIF
!
!  !WRITE(6,*) 'ENTERING ADDFTBU_INIT,done'
!  RETURN
!
!END SUBROUTINE ADDDFTBU_INIT


SUBROUTINE ADDDFTBU(FIRSTCALL)

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE COULOMBARRAY
  USE NONOARRAY
  USE KSPACEARRAY, ONLY : NKTOT, SK, HK, KBO
  USE MYPRECISION
  USE DMARRAY

  IMPLICIT NONE

  LOGICAL, INTENT(IN) :: FIRSTCALL

  INTEGER :: I, J, K
  INTEGER :: INDEX, NUMORB
  REAL(LATTEPREC) :: ES, EP, ED, EF
  REAL(LATTEPREC) :: HMOD
  REAL(LATTEPREC) :: X(HDIM,HDIM), Y(HDIM,HDIM), HU(HDIM,HDIM)
  COMPLEX(LATTEPREC) :: CONE, CZERO

  IF (EXISTERROR) RETURN

  IF (KON .EQ. 0) THEN

    IF ( BASISTYPE .EQ. "ORTHO") THEN
         WRITE(*,*) ' DM HUBBARD U ONLY FOR NONORTHOGONAL REPRESENTATIONS '
         STOP
    ENDIF

    DO J = 1, HDIM
    DO I = 1, HDIM
       SU(I,J) = SMAT(I,J)*DFTB_U(J)
    ENDDO
    ENDDO
  
    IF (FIRSTCALL) THEN
      call DGEMM('T','N',HDIM,HDIM,HDIM,ONE,BO,HDIM,SMAT,HDIM,ZERO,PS,HDIM)  !! Should be PNO instead of BO for propagation in XL-BOMD
    ELSE
      !!! call DGEMM('T','N',HDIM,HDIM,HDIM,ONE,BO,HDIM,SMAT,HDIM,ZERO,PS,HDIM) !! SHOUDL BE USED FOR SCF CONVERGENCE
      call DGEMM('T','N',HDIM,HDIM,HDIM,ONE,PNO,HDIM,SMAT,HDIM,ZERO,PS,HDIM)    !! PNO = propagated density matrix dynamical variable
    ENDIF
    call DGEMM('T','N',HDIM,HDIM,HDIM,ONE,PS,HDIM,SU,HDIM,ZERO,X,HDIM)
    call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,SU,HDIM,PS,HDIM,ZERO,Y,HDIM)
    HU = 0.125D0*(2.D0*SU-X-Y)
    DO J = 1,HDIM
    DO I = 1,HDIM
       H_U(I,J) = HU(I,J) + HU(J,I)
    ENDDO
    ENDDO
    H = H + H_U

  ELSE
    CZERO = CMPLX(ZERO)
    CONE = CMPLX(ONE)

    ! ZY: this part is moved to initiatedm.f90
    DO K = 1, NKTOT
      DO J = 1, HDIM
      DO I = 1, HDIM
          SUK(I,J,K) = SK(I,J,K)*DFTB_U(J)
      ENDDO
      ENDDO
    ENDDO

    DO K = 1, NKTOT
       IF (FIRSTCALL) THEN
           CALL ZGEMM('C','N',HDIM,HDIM,HDIM,CONE,KBO(:,:,K),HDIM,SK(:,:,K),HDIM,CZERO,PSK,HDIM)  !! Should be PNO instead of BO for propagation in XL-BOMD
       ELSE
           WRITE(6,*) 'DFBTU FOR K SPACE IS NOT READY YET!'
           STOP
           !CALL DGEMM('T','N',HDIM,HDIM,HDIM,CONE,PNO,HDIM,SMAT,HDIM,CZERO,PSK,HDIM)    !! PNO = propagated density matrix dynamical variable
       ENDIF

       CALL ZGEMM('C','N',HDIM,HDIM,HDIM,CONE,PSK,HDIM,SUK(1,1,K),HDIM,CZERO,XK,HDIM)
       CALL ZGEMM('N','N',HDIM,HDIM,HDIM,CONE,SUK(1,1,K),HDIM,PSK,HDIM,CZERO,YK,HDIM)

       HK_U = 0.125D0*(2.D0*SUK(:,:,K)-XK-YK)

       HK_U = HK_U + TRANSPOSE(CONJG(HK_U))

       HK(:,:,K) = HK(:,:,K) + HK_U

    ENDDO

  ENDIF

  RETURN

END SUBROUTINE ADDDFTBU

