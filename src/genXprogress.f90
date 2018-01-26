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

!> To produce a matrix \f$Z\f$ which is needed to orthogonalize \f$H\f$.
!!
!! \f$ H_{orth} = Z^{\dagger}HZ \f$
!! See Negre 2016 \cite Negre2016
!!
MODULE GENXPROGRESS

#ifdef PROGRESSON

  USE BML
  USE PRG_GENZ_MOD

  USE CONSTANTS_MOD
  USE NONOARRAY
  USE MYPRECISION

  PRIVATE

  PUBLIC :: GENXBML

  TYPE(BML_MATRIX_T), PUBLIC                :: OVER_BML !FOR STORING S IN BML FORMAT
  TYPE(BML_MATRIX_T), PUBLIC                :: ZK1_BML, ZK2_BML, ZK3_BML
  TYPE(BML_MATRIX_T), PUBLIC                :: ZK4_BML, ZK5_BML, ZK6_BML, ZMAT_BML
  INTEGER, PUBLIC                           :: IGENX = 0!COUTER TO KEEP TRACK OF THE TIMES ZMAT IS COMPUTED.
  TYPE(GENZSPINP), PUBLIC                   :: ZSP

CONTAINS

  SUBROUTINE GENXBML

    IMPLICIT NONE

    !> Parsing Z sparse propagation and build.
    IF(IGENX == 0)THEN
       CALL PRG_PARSE_ZSP(ZSP,"latte.in")
       IF(ZSP%BML_TYPE == BML_MATRIX_DENSE .AND. SPARSEON == 1) STOP 'If CONTROL{ SPARSEON= 1 } then ZSP{ BMLType= Ellpack }'
       IF(ZSP%BML_TYPE .EQ. BML_MATRIX_ELLPACK)THEN
          IF(ZSP%ZSP .EQV. .FALSE.)STOP 'If ZSP{ ZSP= F } then ZSP{ BMLType= Dense }'
          IF(SPARSEON == 0) STOP 'If CONTROL{ SPARSEON= 0 } then ZSP{ BMLType= Dense }'
       ENDIF
    ENDIF

    IF(bml_allocated(OVER_BML)) CALL BML_DEALLOCATE(OVER_BML)

    IF(ZSP%VERBOSE.GE.1)WRITE(*,*)"Inside genx ..."

    IF(ZSP%MDIM < 0) ZSP%MDIM = HDIM


    IGENX = IGENX + 1 !Counter to keep track of the iterations (md and optimization)

    IF(VERBOSE >= 1) WRITE(*,*)"IGENX =",IGENX

    CALL BML_ZERO_MATRIX(ZSP%BML_TYPE,BML_ELEMENT_REAL,LATTEPREC,HDIM,ZSP%MDIM,OVER_BML)

    CALL BML_IMPORT_FROM_DENSE(ZSP%BML_TYPE, &
         SMAT, OVER_BML, ZERO, ZSP%MDIM)

    !    CALL BML_PRINT_MATRIX("OVER_BML",OVER_BML,1,10,1,10)

    IF(ZSP%ZSP)THEN !Congruence transformation.

       IF(IGENX == 0) THEN
          CALL PRG_INIT_ZSPMAT(IGENX,ZK1_BML,ZK2_BML,ZK3_BML&
               ,ZK4_BML,ZK5_BML,ZK6_BML,ZSP%MDIM,ZSP%BML_TYPE)
       ENDIF

       CALL PRG_BUILDZSPARSE(OVER_BML,ZMAT_BML,IGENX,ZSP%MDIM,&
            ZSP%BML_TYPE, ZK1_BML,ZK2_BML,ZK3_BML&
            ,ZK4_BML,ZK5_BML,ZK6_BML,ZSP%NFIRST,ZSP%NREFI,ZSP%NREFF,&
            ZSP%NUMTHRESI,ZSP%NUMTHRESF,ZSP%INTEGRATION,ZSP%VERBOSE)

    ELSE

       !Build z matrix using diagonalization (usual method).
       CALL PRG_BUILDZDIAG(OVER_BML,ZMAT_BML,ZSP%NUMTHRESF,ZSP%MDIM,ZSP%BML_TYPE)

    ENDIF

    CALL FLUSH(6)

    CALL BML_EXPORT_TO_DENSE(ZMAT_BML, XMAT)

    IF(ZSP%VERBOSE.GE.1)WRITE(*,*)"Out of genx ..."

  END SUBROUTINE GENXBML

#endif

END MODULE GENXPROGRESS
