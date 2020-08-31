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

SUBROUTINE ADDDFTBUPRG(FIRSTCALL)

#ifdef PROGRESSON

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE COULOMBARRAY
  USE NONOARRAY
  USE KSPACEARRAY, ONLY : NKTOT, SK, HK, KBO
  USE MYPRECISION
  USE DMARRAY
  USE LATTEPARSER
  USE GENXPROGRESS
  USE SPARSEARRAY
  USE SPINARRAY
  USE BML

  IMPLICIT NONE

  LOGICAL, INTENT(IN) :: FIRSTCALL

  INTEGER, PARAMETER :: dp = LATTEPREC
  INTEGER :: I, J, K
  INTEGER :: INDEX, NUMORB
  REAL(LATTEPREC) :: ES, EP, ED, EF
  REAL(LATTEPREC) :: HMOD
  REAL(LATTEPREC),allocatable :: X(:,:), Y(:,:), HU(:,:)
  COMPLEX(LATTEPREC) :: CONE, CZERO
  TYPE(BML_MATRIX_T) :: BO_BML, PS_BML, AUX_BML, SU_BML, X_BML

  IF (EXISTERROR) RETURN

  ALLOCATE(X(HDIM,HDIM), Y(HDIM,HDIM), HU(HDIM,HDIM))

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

    !IF(bml_allocated(OVER_BML)) THEN
    !   write(6,*) 'test bml1'
    !ENDIF
    IF(bml_allocated(OVER_BML)) CALL BML_DEALLOCATE(OVER_BML)

    CALL BML_ZERO_MATRIX(LT%BML_TYPE,BML_ELEMENT_REAL,LATTEPREC,HDIM,LT%MDIM,OVER_BML)
    CALL BML_IMPORT_FROM_DENSE(LT%BML_TYPE, SMAT, OVER_BML, ZERO, LT%MDIM)

    CALL BML_ZERO_MATRIX(LT%BML_TYPE,BML_ELEMENT_REAL,LATTEPREC,HDIM,LT%MDIM,BO_BML)
    CALL BML_ZERO_MATRIX(LT%BML_TYPE,BML_ELEMENT_REAL,LATTEPREC,HDIM,LT%MDIM,PS_BML)
    CALL BML_ZERO_MATRIX(LT%BML_TYPE,BML_ELEMENT_REAL,LATTEPREC,HDIM,LT%MDIM,AUX_BML)
    CALL BML_ZERO_MATRIX(LT%BML_TYPE,BML_ELEMENT_REAL,LATTEPREC,HDIM,LT%MDIM,SU_BML)
    CALL BML_ZERO_MATRIX(LT%BML_TYPE,BML_ELEMENT_REAL,LATTEPREC,HDIM,LT%MDIM,X_BML)
  
    CALL BML_IMPORT_FROM_DENSE(LT%BML_TYPE, SU, SU_BML, ZERO, LT%MDIM)

    IF (FIRSTCALL) THEN
      !! Should be PNO instead of BO for propagation in XL-BOMD
      !call DGEMM('T','N',HDIM,HDIM,HDIM,ONE,BO,HDIM,SMAT,HDIM,ZERO,PS,HDIM)

      CALL BML_IMPORT_FROM_DENSE(LT%BML_TYPE, BO, BO_BML, ZERO, LT%MDIM)
      CALL BML_TRANSPOSE(BO_BML, AUX_BML)
      CALL BML_MULTIPLY(AUX_BML,OVER_BML,PS_BML,1.0_dp,0.0_dp,NUMTHRESH)
    ELSE
      !! PNO = propagated density matrix dynamical variable
      !call DGEMM('T','N',HDIM,HDIM,HDIM,ONE,PNO,HDIM,SMAT,HDIM,ZERO,PS,HDIM)

      CALL BML_IMPORT_FROM_DENSE(LT%BML_TYPE, PNO, BO_BML, ZERO, LT%MDIM)
      CALL BML_TRANSPOSE(BO_BML, AUX_BML)
      CALL BML_MULTIPLY(AUX_BML,OVER_BML,PS_BML,1.0_dp,0.0_dp,NUMTHRESH)
    ENDIF

    !call DGEMM('T','N',HDIM,HDIM,HDIM,ONE,PS,HDIM,SU,HDIM,ZERO,X,HDIM)
    !call DGEMM('N','N',HDIM,HDIM,HDIM,ONE,SU,HDIM,PS,HDIM,ZERO,Y,HDIM)

    !CALL BML_IMPORT_FROM_DENSE(LT%BML_TYPE, PS, PS_BML, ZERO, LT%MDIM)

    CALL BML_TRANSPOSE(PS_BML, AUX_BML)
    CALL BML_MULTIPLY(AUX_BML,SU_BML,X_BML,1.0_dp,0.0_dp,NUMTHRESH)
    CALL BML_EXPORT_TO_DENSE(X_BML, X)

    CALL BML_MULTIPLY(SU_BML,PS_BML,X_BML,1.0_dp,0.0_dp,NUMTHRESH)
    CALL BML_EXPORT_TO_DENSE(X_BML, Y)

    !CALL BML_IMPORT_FROM_DENSE(LT%BML_TYPE, X, X_BML, ZERO, LT%MDIM)
    !CALL BML_EXPORT_TO_DENSE(X_BML, X)
    !CALL BML_IMPORT_FROM_DENSE(LT%BML_TYPE, Y, X_BML, ZERO, LT%MDIM)
    !CALL BML_EXPORT_TO_DENSE(X_BML, Y)

    CALL BML_DEALLOCATE(BO_BML)
    CALL BML_DEALLOCATE(PS_BML)
    CALL BML_DEALLOCATE(SU_BML)
    CALL BML_DEALLOCATE(X_BML)
    CALL BML_DEALLOCATE(AUX_BML)

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

  DEALLOCATE(X,Y,HU)
  RETURN
#endif

END SUBROUTINE ADDDFTBUPRG
