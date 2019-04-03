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

  IMPLICIT NONE
  SAVE

  REAL(LATTEPREC), ALLOCATABLE :: PO_0(:,:), PO_1(:,:), PO_2(:,:), PO_3(:,:)  ! For DM integration
  REAL(LATTEPREC), ALLOCATABLE :: PO_4(:,:), PO_5(:,:), PO_6(:,:), PO_7(:,:)  ! For DM integration
  REAL(LATTEPREC), ALLOCATABLE :: POrth(:,:), DOrth(:,:), DELTAQDM(:), C_DFTB_U ! C_DFTB_U is the chosen Hubbard U for the DM term
  REAL(LATTEPREC), ALLOCATABLE :: DFTB_U(:), IXMAT(:,:), PS(:,:),SU(:,:), Delta_DS(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: HubForce(:,:), DOrth_old(:,:), AHub(:,:)
  REAL(LATTEPREC) :: EHub 
  INTEGER, ALLOCATABLE :: H_INDEX_START(:), H_INDEX_END(:) ! Index list of start and end position of matrix elements of atom I

END MODULE DMARRAY
