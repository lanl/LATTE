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

SUBROUTINE INITSHOCKCOMP

  USE CONSTANTS_MOD
  USE MDARRAY

  IMPLICIT NONE
  IF (EXISTERROR) RETURN

  ! First convert units of Up and Us from m/s to A/fs

  UPARTICLE = UPARTICLE*1.0D-5
  USHOCK = USHOCK*1.0D-5

  ! If USHOCK is positive, then we'll use it. If it's negative then
  ! we'll compute it using the sound velocity (input) and the Universal 
  ! liquid Hugoniot

  IF (USHOCK .LT. 0.0) THEN

     ! Converting sound velocity from m/s to A/fs

     C0 = C0 * 1.0D-5

     USHOCK = C0 * (1.37D0 - 0.37D0*EXP( -TWO*UPARTICLE / C0 )) &
          + 1.62D0*UPARTICLE

  ENDIF

  !  PRINT*, "Up = ", UPARTICLE, "Us = ", USHOCK
  !
  ! The duration of the shock, i.e., the time taken for the 
  ! shock front to traverse the simulation cell is t_dur = l_0/Us
  !
  ! So... we'll start the Hugoniostat at timestep SHOCKSTART and 
  ! turn it off INT(L_O/(U_s *dt)) time steps later
  !

  ! If the box has 90 degree angles, then its length in the three
  ! directions is BOX(1,1), BOX(2,2), AND BOX(3,3)...

  SHOCKSTOP = SHOCKSTART + &
       INT(BOX(SHOCKDIR,SHOCKDIR)/(USHOCK * DT))
  !  SHOCKSTOP = SHOCKSTART + &
  !       INT((BOX(2,SHOCKDIR) - BOX(1,SHOCKDIR))/(USHOCK * DT))

  !  PRINT*, "start, stop = ", SHOCKSTART, SHOCKSTOP

END SUBROUTINE INITSHOCKCOMP
  
