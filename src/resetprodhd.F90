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

SUBROUTINE RESETPRODHD

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE KSPACEARRAY
  USE MYPRECISION

  IMPLICIT NONE

  !
  ! Let's try keeping Sum_i H_ii D_ii constant = 0
  !

  INTEGER :: I, K, J, NUMORB, INDEX
  REAL(LATTEPREC) :: MYSUM
  COMPLEX(LATTEPREC) :: ZMYSUM
  IF (EXISTERROR) RETURN

  MYSUM = ZERO
  ZMYSUM = CMPLX(ZERO)

  IF (KON .EQ. 0) THEN ! Real-space

     DO I = 1, HDIM
        MYSUM = MYSUM + H(I,I)*BO(I,I)
     ENDDO

     MYSUM = MYSUM/FLOAT(HDIM)

     DO I = 1, HDIM
        H(I,I) = H(I,I) - MYSUM
     ENDDO

  ELSE ! k-space

     DO K = 1, NKTOT
        INDEX = 0
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

           DO J = 1, NUMORB
              INDEX = INDEX + 1
              ZMYSUM = ZMYSUM + LCNSHIFT(I)*KBO(INDEX,INDEX,K)
           ENDDO
        ENDDO
     ENDDO

     ZMYSUM = ZMYSUM/REAL(HDIM*NKTOT)

     !     PRINT*, REAL(ZMYSUM)

     DO I = 1, NATS

        LCNSHIFT(I) = LCNSHIFT(I) - REAL(ZMYSUM)

     ENDDO

  ENDIF

  RETURN

END SUBROUTINE RESETPRODHD

