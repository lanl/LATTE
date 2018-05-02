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

MODULE TIMER_MOD

  IMPLICIT NONE
  SAVE

  INTEGER :: TSTART_CLOCK, TSTOP_CLOCK, TCLOCK_RATE, TCLOCK_MAX
  INTEGER :: TX, NUM_TIMERS
  INTEGER :: LATTE_TIMER, DENSE2SPARSE_TIMER, DMBUILD_TIMER, SPARSE2DENSE_TIMER
  INTEGER :: SP2ALL_TIMER, SP2SPARSE_TIMER
  INTEGER, ALLOCATABLE :: TSTART(:), TTOTAL(:), TCOUNT(:)
  REAL, ALLOCATABLE :: TAVG(:), TSUM(:), TPERCENT(:)
  CHARACTER(LEN=20), ALLOCATABLE :: TNAME(:)

CONTAINS

  ! Initialize timers
  !
  FUNCTION INIT_TIMER()

    INTEGER :: I, INIT_TIMER

    NUM_TIMERS = 6

    IF(.NOT.ALLOCATED(TSTART)) ALLOCATE(TSTART(NUM_TIMERS), TTOTAL(NUM_TIMERS), TCOUNT(NUM_TIMERS))
    IF(.NOT.ALLOCATED(TNAME)) ALLOCATE(TNAME(NUM_TIMERS))
    IF(.NOT.ALLOCATED(TAVG))THEN
       ALLOCATE(TAVG(NUM_TIMERS),TSUM(NUM_TIMERS),TPERCENT(NUM_TIMERS))
    ENDIF

    ! Timer handles, names, and counters
    LATTE_TIMER = 1
    SP2ALL_TIMER = 2
    SP2SPARSE_TIMER = 3
    DENSE2SPARSE_TIMER = 4
    DMBUILD_TIMER = 5
    SPARSE2DENSE_TIMER = 6

    TNAME(LATTE_TIMER) = "LATTE"
    TNAME(SP2ALL_TIMER) = "Sp2All"
    TNAME(SP2SPARSE_TIMER) = "  Sp2Sparse"
    TNAME(DENSE2SPARSE_TIMER) = "    Dense2Sparse"
    TNAME(DMBUILD_TIMER) = "    DMBuild"
    TNAME(SPARSE2DENSE_TIMER) = "  Sparse2Dense"

    TTOTAL = 0
    TCOUNT = 0

    INIT_TIMER = NUM_TIMERS

  END FUNCTION INIT_TIMER

  ! Done with timers
  !
  FUNCTION SHUTDOWN_TIMER()

    INTEGER :: SHUTDOWN_TIMER

    DEALLOCATE(TSTART, TTOTAL, TCOUNT)
    DEALLOCATE(TNAME)

    SHUTDOWN_TIMER = NUM_TIMERS

  END FUNCTION SHUTDOWN_TIMER

  ! Start Timing
  !
  FUNCTION START_TIMER(ITIMER)

    INTEGER :: ITIMER, START_TIMER

    CALL SYSTEM_CLOCK(TSTART_CLOCK, TCLOCK_RATE, TCLOCK_MAX)
    TSTART(ITIMER) = TSTART_CLOCK

    START_TIMER = TSTART_CLOCK

  END FUNCTION START_TIMER

  ! Stop timing
  !
  FUNCTION STOP_TIMER(ITIMER)

    INTEGER :: ITIMER, TDELTA, STOP_TIMER

    CALL SYSTEM_CLOCK(TSTOP_CLOCK, TCLOCK_RATE, TCLOCK_MAX)
    TDELTA = TSTOP_CLOCK - TSTART(ITIMER)
    TTOTAL(ITIMER) = TTOTAL(ITIMER) + TDELTA
    TCOUNT(ITIMER) = TCOUNT(ITIMER) + 1

    STOP_TIMER = TSTOP_CLOCK

  END FUNCTION STOP_TIMER

  ! Print performance results
  !
  FUNCTION TIMER_RESULTS()

    INTEGER :: I, TIMER_RESULTS

    PRINT *, ""
    WRITE(6,*) "Timer                 # Calls  Avg/Call (s)     Total (s)       % Time"
    PRINT *, ""

    DO I = 1, NUM_TIMERS

       IF (TCOUNT(I) .GT. 0) THEN

          TAVG(I) = (FLOAT(TTOTAL(I))/FLOAT(TCLOCK_RATE))/FLOAT(TCOUNT(I))
          TSUM(I) = FLOAT(TTOTAL(I))/FLOAT(TCLOCK_RATE)
          TPERCENT(I) = (TSUM(I) / TSUM(1)) * 100.0

          WRITE(6,10) TNAME(I), TCOUNT(I), TAVG(I), TSUM(I), TPERCENT(I)
10        FORMAT(A25, I4, 3G16.6)
       ENDIF

    ENDDO

    TIMER_RESULTS = NUM_TIMERS

  END FUNCTION TIMER_RESULTS

  ! Print a tag time and date
  !
  SUBROUTINE TIMEDATE_TAG(TAG)
    CHARACTER(LEN=*) :: TAG
    INTEGER :: VALUES(8)

    CALL DATE_AND_TIME(VALUES=VALUES)

    WRITE(*,'(A2,1X,A,1X,I2,A1,I2,2X,A2,1X,I2,A1,I2,A1,I4)') " #",TRIM(TAG),&
         & VALUES(5),":",VALUES(6),&
         & "on",VALUES(2),"/",VALUES(3),"/",VALUES(1)

  END SUBROUTINE TIMEDATE_TAG

  ! Get the actual time in mls
  !
  FUNCTION TIME_MLS()
    REAL(8) :: TIME_MLS
    INTEGER :: TIMEVECTOR(8)

    TIME_MLS = 0.0d0
    CALL DATE_AND_TIME(VALUES=TIMEVECTOR)
    TIME_MLS=TIMEVECTOR(5)*60.0d0*60.0d0*1000.0d0 + TIMEVECTOR(6)*60.0d0*1000.0d0 &
       & + TIMEVECTOR(7)*1000.0d0 + TIMEVECTOR(8)

  END FUNCTION TIME_MLS


END MODULE TIMER_MOD
