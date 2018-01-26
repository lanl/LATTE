SUBROUTINE INITSHOCKCOMP

  USE CONSTANTS_MOD
  USE MDARRAY

  IMPLICIT NONE

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
  
