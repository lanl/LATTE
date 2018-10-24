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

!> To implement mixing schemes from the progress library
!!
MODULE MIXER_MOD

#ifdef PROGRESSON

  USE MYPRECISION
  USE COULOMBARRAY
  USE SETUPARRAY
  USE PRG_PULAYMIXER_MOD

  PRIVATE

  PUBLIC :: QMIXPRG

  !For mixing scheme
  LOGICAL, PUBLIC                      ::  MIXINIT = .FALSE.
  REAL(LATTEPREC), ALLOCATABLE, PUBLIC  ::  DQIN(:,:), DQOUT(:,:)
  REAL(LATTEPREC), PUBLIC              ::  SCFERROR
  TYPE(MX_TYPE), PUBLIC                ::  MX

CONTAINS

  SUBROUTINE QMIXPRG(PITER)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: PITER

    IF(MX%MIXERTYPE == "Linear")THEN
       CALL PRG_LINEARMIXER(DELTAQ,OLDDELTAQS,SCFERROR,MX%MIXCOEFF,MX%VERBOSE)
    ELSEIF(MX%MIXERTYPE == "Pulay")THEN
       CALL PRG_QMIXER(DELTAQ,OLDDELTAQS,DQIN,DQOUT,SCFERROR,PITER,MX%MIXCOEFF,MX%MPULAY,MX%VERBOSE)
    ELSEIF(MX%MIXERTYPE == "PulayQlist")THEN
       CALL PRG_QMIXER(QLIST,OLDDELTAQS,DQIN,DQOUT,SCFERROR,PITER,MX%MIXCOEFF,MX%MPULAY,MX%VERBOSE)
    ELSE
       CALL ERRORS("mixer_mod:qmixprg","Mixing scheme not implemented. &
            & Check MixerType keyword in the input file")
    ENDIF

  END SUBROUTINE QMIXPRG

#endif

END MODULE MIXER_MOD
