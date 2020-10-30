!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

MODULE DMARRAY

  USE MYPRECISION
#ifdef PROGRESSON
  USE BML 
#endif

  IMPLICIT NONE
  SAVE

  REAL(LATTEPREC), ALLOCATABLE :: PO_0(:,:), PO_1(:,:), PO_2(:,:), PO_3(:,:)  ! For DM integration
  REAL(LATTEPREC), ALLOCATABLE :: PO_4(:,:), PO_5(:,:), PO_6(:,:), PO_7(:,:)  ! For DM integration
  REAL(LATTEPREC), ALLOCATABLE :: POrth(:,:), DOrth(:,:), DELTAQDM(:) !, C_DFTB_U ! C_DFTB_U is the chosen Hubbard U for the DM term
  REAL(LATTEPREC), ALLOCATABLE :: PNO(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: DFTB_U(:), IXMAT(:,:), PS(:,:),SU(:,:), Delta_DS(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: HubForce(:,:), DOrth_old(:,:), AHub(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: H_U(:,:)
  REAL(LATTEPREC) :: EHub 
  INTEGER, ALLOCATABLE :: H_INDEX_START(:), H_INDEX_END(:) ! Index list of start and end position of matrix elements of atom I

  ! kspace array
  COMPLEX(LATTEPREC), ALLOCATABLE :: SUK(:,:,:)
  COMPLEX(LATTEPREC), ALLOCATABLE :: HK_U(:,:)
  COMPLEX(LATTEPREC), ALLOCATABLE :: XK(:,:), YK(:,:), PSK(:,:), HUK(:,:)
  COMPLEX(LATTEPREC), ALLOCATABLE :: UDK(:,:),UPSDK(:,:),PSUDK(:,:),AHUBK(:,:,:)
  COMPLEX(LATTEPREC), ALLOCATABLE :: DORK(:,:,:),DORK_OLD(:,:,:) !Orthogonal DM for k points

#ifdef PROGRESSON
  TYPE(BML_MATRIX_T), PUBLIC :: PO0_BML, PO1_BML, PO2_BML, PO3_BML
  TYPE(BML_MATRIX_T), PUBLIC :: PO4_BML, PO5_BML, PO6_BML
  TYPE(BML_MATRIX_T), PUBLIC :: PO_BML
#endif
  
END MODULE DMARRAY
