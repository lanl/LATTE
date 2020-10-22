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
! subroutine HubbardForce(HubForce,EHub,Nr_atoms,HDIM,U,D,P,S,dSx,dSy,dSz,H_INDEX_START,H_INDEX_END)
! D*S*D = D, P = Dynamical variable density matrix in non-orthogonal atomic orbital representation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE HUBBARDFORCEPRG  ! AND ENERGY

#ifdef PROGRESSON 

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE COULOMBARRAY
  USE NONOARRAY
  USE KSPACEARRAY
  USE MYPRECISION
  USE SPARSEARRAY
  USE DMARRAY
  USE BML
  USE LATTEPARSER
  USE GENXPROGRESS
  USE NONOARRAYPROGRESS

  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = LATTEPREC
  INTEGER :: I, J, K
  INTEGER :: INDEX, NUMORB
  REAL(LATTEPREC) :: ES, EP, ED, EF
  !REAL(LATTEPREC) :: HMOD, D(HDIM,HDIM), P(HDIM,HDIM)
  REAL(LATTEPREC) :: UD(HDIM,HDIM), UP(HDIM,HDIM)
  !REAL(LATTEPREC) :: DSUP(HDIM,HDIM), UPSD(HDIM,HDIM), PSUD(HDIM,HDIM), UDSP(HDIM,HDIM), PSUP(HDIM,HDIM)
  REAL(LATTEPREC), ALLOCATABLE :: B(:,:)
  REAL(LATTEPREC) :: UPSP(HDIM,HDIM), A(HDIM,HDIM)
  REAL(LATTEPREC) :: Xtmp, Ytmp, Ztmp
  COMPLEX(LATTEPREC) :: ZHMOD
  COMPLEX(LATTEPREC) :: CONE, CZERO,CTMP
  INTEGER :: I_A, I_B, L, II, JJ, J_A, J_B  
  TYPE(BML_MATRIX_T) :: UD_BML, UP_BML, AUX1_BML,AUX2_BML

  IF (EXISTERROR) RETURN

  ALLOCATE(B(HDIM,HDIM))

  IF (KON==0) then
    !D = 0.5D0*BO   ! D*S*D = D, BO*S*BO = 2*BO
    !P = 0.5D0*PNO  ! D*S*D = D, SHOULD BE P = D for SCF CONVERGENCE
  
    do I = 1, HDIM
    do J = 1, HDIM
       UD(I,J) = 0.5_DP*DFTB_U(I)*BO(I,J)
       UP(I,J) = 0.5_DP*DFTB_U(I)*PNO(I,J)
    enddo
    enddo

    !CALL BML_ZERO_MATRIX(LT%BML_TYPE, BML_ELEMENT_REAL,LATTEPREC, HDIM, HDIM, D_BML)
    !CALL BML_ZERO_MATRIX(LT%BML_TYPE, BML_ELEMENT_REAL,LATTEPREC, HDIM, HDIM, P_BML)
    CALL BML_ZERO_MATRIX(LT%BML_TYPE, BML_ELEMENT_REAL,LATTEPREC, HDIM, HDIM, UD_BML)
    CALL BML_ZERO_MATRIX(LT%BML_TYPE, BML_ELEMENT_REAL,LATTEPREC, HDIM, HDIM, UP_BML)
    CALL BML_ZERO_MATRIX(LT%BML_TYPE, BML_ELEMENT_REAL,LATTEPREC, HDIM, HDIM, AUX1_BML)
    CALL BML_ZERO_MATRIX(LT%BML_TYPE, BML_ELEMENT_REAL,LATTEPREC, HDIM, HDIM, AUX2_BML)
 
    !CALL BML_IMPORT_FROM_DENSE(ZSP%BML_TYPE, D, D_BML, ZERO, ZSP%MDIM)
    CALL BML_IMPORT_FROM_DENSE(ZSP%BML_TYPE, PNO, PNO_BML, ZERO, ZSP%MDIM)
    CALL BML_IMPORT_FROM_DENSE(ZSP%BML_TYPE, UD, UD_BML, ZERO, ZSP%MDIM)
    CALL BML_IMPORT_FROM_DENSE(ZSP%BML_TYPE, UP, UP_BML, ZERO, ZSP%MDIM)
 
    !D*s,  UP * DS^T; PS*UD
    !CALL BML_MULTIPLY(D_BML,OVER_BML,AUX2_BML,1.0_dp,0.0_dp,NUMTHRESH)
    ! 0.5 * D * S
    CALL BML_MULTIPLY(BO_BML,OVER_BML,AUX2_BML,0.5_dp,0.0_dp,NUMTHRESH)
    CALL BML_TRANSPOSE(AUX2_BML, AUX1_BML)
    CALL BML_MULTIPLY(UP_BML,AUX1_BML,AUX2_BML,1.0_dp,0.0_dp,NUMTHRESH)

    ! UD - UPSD
    CALL BML_ADD(UD_BML,AUX2_BML,1.0_dp,-1.0_dp,NUMTHRESH)
 
    !( UD - UPSD) * SMat
    CALL BML_MULTIPLY(UD_BML,OVER_BML,AUX1_BML,1.0_dp,0.0_dp,NUMTHRESH)
  
    CALL BML_EXPORT_TO_DENSE(AUX1_BML, B)
    EHub = 0.D0 
    do I = 1,HDIM
      EHub = EHub + 0.5D0*B(I,I)  ! Hubbar Ebergy contribution
    enddo
    write(6,*)'EHUB=', EHUB
  
    ! PSUD stored in AUX2_BML
    !0.5*PNO, so there is a 0.5 factor
    CALL BML_MULTIPLY(PNO_BML,OVER_BML,AUX1_BML,0.5_dp,0.0_dp,NUMTHRESH)
    CALL BML_MULTIPLY(AUX1_BML,UD_BML,AUX2_BML,1.0_dp,0.0_dp,NUMTHRESH)

    ! B=(UD - UPSD - PSUD)
    ! UD already contains -USPD, so, here only sub PSUD
    CALL BML_ADD(UD_BML,AUX2_BML,1.0_dp,-1.0_dp,NUMTHRESH)

    ! B + B^T -> Ahub
    CALL BML_TRANSPOSE(UD_BML, AUX1_BML)
    CALL BML_ADD(UD_BML,AUX1_BML,1.0_dp, 1.0_dp,NUMTHRESH)

    CALL BML_EXPORT_TO_DENSE(UD_BML, Ahub)

    !CALL BML_DEALLOCATE(D_BML)
    !CALL BML_DEALLOCATE(P_BML)
    CALL BML_DEALLOCATE(UD_BML)
    CALL BML_DEALLOCATE(UP_BML)
    CALL BML_DEALLOCATE(AUX1_BML)
    CALL BML_DEALLOCATE(AUX2_BML)
    
  ELSE
    CZERO = CMPLX(ZERO)
    CONE = CMPLX(ONE)

    CTMP = CZERO
    DO K = 1, NKTOT
      DO I = 1, HDIM
      DO J = 1, HDIM
        UDK(I,J) = DFTB_U(I) * KBO(I,J,K) * 0.5D0
      ENDDO
      ENDDO

      CALL ZGEMM('N','N',HDIM,HDIM,HDIM,CONE*0.5D0,KBO(:,:,K),HDIM,SK(:,:,K),HDIM,CZERO,PSK,HDIM)

      CALL ZGEMM('N','C',HDIM,HDIM,HDIM,CONE,UDK,HDIM,PSK,HDIM,CZERO,UPSDK,HDIM)
      CALL ZGEMM('N','N',HDIM,HDIM,HDIM,CONE,PSK,HDIM,UDK,HDIM,CZERO,PSUDK,HDIM)

      XK = UDK - UPSDK

      CALL ZGEMM('N','N',HDIM,HDIM,HDIM,CONE,XK,HDIM,SK(:,:,K),HDIM,CZERO,YK,HDIM)
      DO I = 1,HDIM
        CTMP = CTMP + 0.5D0*YK(I,I)  ! Hubbar Ebergy contribution
        !EHub = EHub + 0.5D0*YK(I,I)  ! Hubbar Ebergy contribution
      ENDDO

      YK = UDK - UPSDK - PSUDK

      DO I = 1, HDIM
      DO J = 1, HDIM
        AHUBK(I,J,K) = YK(I,J) + dconjg(YK(J,I))
      ENDDO
      ENDDO
      !
    ENDDO ! k point loop

    EHUB = REAL(CTMP)/REAL(NKTOT)

  ENDIF  

  DEALLOCATE(B)

  RETURN

#endif

END SUBROUTINE HUBBARDFORCEPRG
