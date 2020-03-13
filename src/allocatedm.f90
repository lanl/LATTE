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

SUBROUTINE ALLOCATEDM

  USE CONSTANTS_MOD
!  USE XBOARRAY
!  USE SPINARRAY
  USE MYPRECISION
  USE DMARRAY

  IMPLICIT NONE

  IF (EXISTERROR) RETURN

  !
  ! This depends on whether we're propagating the
  ! on-site H matrix elements (in the case of LCN calculations)
  ! or the atomic partial charges (full SCC-TB)


  ALLOCATE(PNO(HDIM,HDIM), PO_0(HDIM,HDIM), PO_1(HDIM,HDIM), PO_2(HDIM,HDIM), PO_3(HDIM,HDIM))
  ALLOCATE(PO_4(HDIM,HDIM), PO_5(HDIM,HDIM), PO_6(HDIM,HDIM), PO_7(HDIM,HDIM))
  ALLOCATE(POrth(HDIM,HDIM), DOrth(HDIM,HDIM), DOrth_old(HDIM,HDIM), IXMAT(HDIM,HDIM))
  ALLOCATE(H_INDEX_START(NATS), H_INDEX_END(NATS), DELTAQDM(NATS), DFTB_U(HDIM))
  ALLOCATE(PS(HDIM,HDIM), SU(HDIM,HDIM), Delta_DS(HDIM,HDIM))
  ALLOCATE(HubForce(3,NATS),AHub(HDIM,HDIM))

  RETURN

END SUBROUTINE ALLOCATEDM

