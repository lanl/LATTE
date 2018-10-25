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
    INTEGER :: I,J,NUMORB, INDEX
    CHARACTER(20) :: MYMIXERTYPE

    MYMIXERTYPE = MX%MIXERTYPE

    IF(VERBOSE >= 1) WRITE(*,*)"MixerType=", MYMIXERTYPE

    IF(MYMIXERTYPE == "PulayLinear" .AND. PITER >= 10) MYMIXERTYPE = "Linear"

    IF(MYMIXERTYPE == "Linear")THEN

       CALL PRG_LINEARMIXER(DELTAQ,OLDDELTAQS,SCFERROR,MX%MIXCOEFF,MX%VERBOSE)

    ELSEIF(MYMIXERTYPE == "Pulay")THEN

       CALL PRG_QMIXER(DELTAQ,OLDDELTAQS,DQIN,DQOUT,SCFERROR,PITER,MX%MIXCOEFF,MX%MPULAY,MX%VERBOSE)

    ELSEIF(MYMIXERTYPE == "PulayLinear")THEN

        CALL PRG_QMIXER(DELTAQ,OLDDELTAQS,DQIN,DQOUT,SCFERROR,PITER,MX%MIXCOEFF,MX%MPULAY,MX%VERBOSE)

    ELSEIF(MYMIXERTYPE == "PulayQlist")THEN

       IF(PITER == 1) OLDQLIST = QLIST
       CALL PRG_QMIXER(QLIST,OLDQLIST,DQIN,DQOUT,SCFERROR,PITER,MX%MIXCOEFF,MX%MPULAY,MX%VERBOSE)
       IF(.NOT. ALLOCATED(MYCHARGE)) ALLOCATE(MYCHARGE(NATS))
       INDEX = 0
       MYCHARGE = 0.0d0

       DO I = 1, NATS

          SELECT CASE(BASIS(ELEMPOINTER(I)))

          CASE("s")

             NUMORB = 1

          CASE("p")

             NUMORB = 3

          CASE("d")

             NUMORB = 5

          CASE("f")

             NUMORB = 7

          CASE("sp")

             NUMORB = 4

          CASE("sd")

             NUMORB = 6

          CASE("sf")

             NUMORB = 8

          CASE("pd")

             NUMORB = 8

          CASE("pf")

             NUMORB = 10

          CASE("df")

             NUMORB = 12

          CASE("spd")

             NUMORB = 9

          CASE("spf")

             NUMORB = 11

          CASE("sdf")

             NUMORB = 13

          CASE("pdf")

             NUMORB = 15

          CASE("spdf")

             NUMORB = 16

          END SELECT

          !     MYCHARGE = ZERO
          DO J = 1, NUMORB

             INDEX = INDEX + 1
             MYCHARGE(I) = MYCHARGE(I) + QLIST(INDEX)

          ENDDO

          DELTAQ(I) = MYCHARGE(I) - ATOCC(ELEMPOINTER(I))

       ENDDO

       OLDDELTAQS = DELTAQ

    ELSE
       CALL ERRORS("mixer_mod:qmixprg","Mixing scheme not implemented. &
            & Check MixerType keyword in the input file")
    ENDIF


  END SUBROUTINE QMIXPRG

#endif

END MODULE MIXER_MOD
