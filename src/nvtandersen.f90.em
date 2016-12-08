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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\

!!$The algorithm developed by Gronbech-Jensen and Farago is implemented
!!$Molecular Physics 111 (2013) 983-991
!!$
!!$Contribution from Enrique Martinez

SUBROUTINE NVTANDERSEN

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MDARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, K, N
  INTEGER :: ITER
  REAL(LATTEPREC) :: RESCALE
  REAL(LATTEPREC) :: CHI, TMPTEMP, MYTEMP, MOMENTUM(3)
  REAL(LATTEPREC) :: LIM
  INTEGER, ALLOCATABLE :: SEED(:)
  REAL(LATTEPREC) :: RN
  REAL(LATTEPREC), EXTERNAL :: GAUSSAND
  REAL(LATTEPREC) :: MEAN, STDDEV
  INTEGER :: OPTION, TRANSFREE
  REAL(LATTEPREC) :: VELFACTOR
  REAL(LATTEPREC) :: BOLTZ

  OPTION = 3
  TRANSFREE = 1
  SETTH = 0

  if(OPTION.eq.1) then
     
!!$     TMPTEMP = ZERO
!!$     
!!$     DO I = 1, AVEPER
!!$        TMPTEMP = TMPTEMP + THIST(I)
!!$     ENDDO
!!$     
!!$     MYTEMP = TMPTEMP/FLOAT(AVEPER)
     
     SETTH = 0
     
!!$     CALL RANDOM_SEED(SIZE = N)
!!$     
!!$     ALLOCATE(SEED(N))
!!$     
!!$     DO I = 1, N
!!$        SEED(I) = (SEEDTH + ITER) + 3*(I - 1)
!!$     ENDDO
!!$     
!!$     CALL RANDOM_SEED(PUT = SEED)
!!$     DEALLOCATE(SEED)
     
     MEAN = ZERO
     STDDEV = ONE
     MOMENTUM = ZERO
     LIM = FRICTION * CUMDT
     DO I = 1, NATS
        CALL RANDOM_NUMBER(RN)
        IF(RN .LE. LIM) THEN
           VELFACTOR = SQRT(ONE/MASS(ELEMPOINTER(I)))
           V(1,I) = GAUSSAND(MEAN,STDDEV) * VELFACTOR
           V(2,I) = GAUSSAND(MEAN,STDDEV) * VELFACTOR
           V(3,I) = GAUSSAND(MEAN,STDDEV) * VELFACTOR
!!$        write (*,*) V(1,I)
!!$        write (*,*) V(2,I)
!!$        write (*,*) V(3,I)
        ENDIF
        MOMENTUM(1) = MOMENTUM(1) + V(1,I)*MASS(ELEMPOINTER(I))
        MOMENTUM(2) = MOMENTUM(2) + V(2,I)*MASS(ELEMPOINTER(I))
        MOMENTUM(3) = MOMENTUM(3) + V(3,I)*MASS(ELEMPOINTER(I))
     ENDDO
     MOMENTUM = MOMENTUM/SUMMASS
     
     DO I = 1, NATS
        V(1,I) = V(1,I) - MOMENTUM(1)
        V(2,I) = V(2,I) - MOMENTUM(2)
        V(3,I) = V(3,I) - MOMENTUM(3)
     ENDDO

     CALL GETKE
     
     RESCALE = SQRT(TTARGET/TEMPERATURE)
     
     V = RESCALE * V

  ELSEIF(OPTION.EQ.2) THEN

!!$     TMPTEMP = ZERO
!!$     
!!$     DO I = 1, AVEPER
!!$        TMPTEMP = TMPTEMP + THIST(I)
!!$     ENDDO
!!$     
!!$     MYTEMP = TMPTEMP/FLOAT(AVEPER)
!!$     
!!$     SETTH = 0
!!$     
!!$     CALL RANDOM_SEED(SIZE = N)
!!$     
!!$     ALLOCATE(SEED(N))
!!$     
!!$     DO I = 1, N
!!$        SEED(I) = (SEEDTH + ITER) + 3*(I - 1)
!!$     ENDDO
!!$     
!!$     CALL RANDOM_SEED(PUT = SEED)
!!$     DEALLOCATE(SEED)
     
     MEAN = ZERO
     STDDEV = ONE
     MOMENTUM = ZERO
     LIM = FRICTION * CUMDT
     DO I = 1, NATS
        
        VELFACTOR = SQRT(ONE/MASS(ELEMPOINTER(I)))
        V(1,I) = GAUSSAND(MEAN,STDDEV) * VELFACTOR
        V(2,I) = GAUSSAND(MEAN,STDDEV) * VELFACTOR
        V(3,I) = GAUSSAND(MEAN,STDDEV) * VELFACTOR
!!$        write (*,*) V(1,I)
!!$        write (*,*) V(2,I)
!!$        write (*,*) V(3,I)
        
        MOMENTUM(1) = MOMENTUM(1) + V(1,I)*MASS(ELEMPOINTER(I))
        MOMENTUM(2) = MOMENTUM(2) + V(2,I)*MASS(ELEMPOINTER(I))
        MOMENTUM(3) = MOMENTUM(3) + V(3,I)*MASS(ELEMPOINTER(I))
     ENDDO
     MOMENTUM = MOMENTUM/SUMMASS
     
     DO I = 1, NATS
        V(1,I) = V(1,I) - MOMENTUM(1)
        V(2,I) = V(2,I) - MOMENTUM(2)
        V(3,I) = V(3,I) - MOMENTUM(3)
     ENDDO

     CALL GETKE
     
     RESCALE = SQRT(TTARGET/TEMPERATURE)
     
     V = RESCALE * V

  ELSEIF(OPTION.EQ.3) THEN
     
     MEAN = ZERO
     MOMENTUM = ZERO
     BOLTZ = 1.0/KE2T
	write(*,*) 'TTARGET= ',TTARGET
     DO I = 1, NATS
        
        STDDEV = SQRT(BOLTZ*TTARGET/MASS(ELEMPOINTER(I))/MVV2KE)
        !STDDEV = SQRT(BOLTZ*TTARGET/MASS(ELEMPOINTER(I)))

        V(1,I) = GAUSSAND(MEAN,STDDEV)
        V(2,I) = GAUSSAND(MEAN,STDDEV)
        V(3,I) = GAUSSAND(MEAN,STDDEV)
        
        MOMENTUM(1) = MOMENTUM(1) + V(1,I)*MASS(ELEMPOINTER(I))
        MOMENTUM(2) = MOMENTUM(2) + V(2,I)*MASS(ELEMPOINTER(I))
        MOMENTUM(3) = MOMENTUM(3) + V(3,I)*MASS(ELEMPOINTER(I))

     ENDDO
     MOMENTUM = MOMENTUM/SUMMASS
     
     IF (TRANSFREE .EQ. 1) THEN
	open(303, file='velocities.out', status='REPLACE')
        DO I = 1, NATS
           V(1,I) = V(1,I) - MOMENTUM(1)
           V(2,I) = V(2,I) - MOMENTUM(2)
           V(3,I) = V(3,I) - MOMENTUM(3)
	   write(303,*) MASS(ELEMPOINTER(I)),V(1,I),V(2,I),V(3,I)
        ENDDO
	close(303)
     ENDIF

     CALL GETKE
     
  ENDIF

  QITER = QITERAND
  CALL GETMDF(1, 100)
  QITER = QITERIN

  RETURN

END SUBROUTINE NVTANDERSEN
        
FUNCTION GAUSSAND(MEAN,STDDEV)

  USE MDARRAY
  USE MYPRECISION
  
  IMPLICIT NONE
  
  !
  ! Based on GASDEV from Numerical Recipes
  !
  ! We generate 2 random numbers at a time so we have to be
  ! a little careful. 
  !

  REAL(LATTEPREC) :: MEAN, STDDEV
  REAL(LATTEPREC) :: GAUSSAND
  REAL(LATTEPREC) :: V1, V2, RN(2), R, FAC, Y1, Y2
  REAL(LATTEPREC), SAVE :: G2

  IF (SETTH .EQ. 0) THEN

     R = TWO

     DO WHILE (R .GE. ONE .OR. R .EQ. ZERO) 

        CALL RANDOM_NUMBER(RN)
        
        V1 = TWO*RN(1) - ONE
        V2 = TWO*RN(2) - ONE
        
        R = V1*V1 + V2*V2

     ENDDO

     FAC = SQRT( -TWO * LOG(R) / R )
     Y1 = V1 * FAC
     Y2 = V2 * FAC

     G2 = (Y1 * STDDEV) + MEAN
     GAUSSAND = (Y2 * STDDEV) + MEAN

     SETTH = 0

  ELSE 

     GAUSSAND = G2
     
     SETTH = 0

  ENDIF

  RETURN

END FUNCTION GAUSSAND
