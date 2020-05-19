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

SUBROUTINE GENORBITALLIST

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J

  ALLOCATE(ORBITAL_LIST(5,NATS))
  
  ORBITAL_LIST = 0

  DO I = 1, NATS

     SELECT CASE(BASIS(ELEMPOINTER(I)))

     CASE("s")
        ORBITAL_LIST(1,I) = 0
        ORBITAL_LIST(2,I) = -1
     CASE("p")
        ORBITAL_LIST(1,I) = 1
        ORBITAL_LIST(2,I) = -1
     CASE("d")
        ORBITAL_LIST(1,I) = 2
        ORBITAL_LIST(2,I) = -1
     CASE("f")
        ORBITAL_LIST(1,I) = 3
        ORBITAL_LIST(2,I) = -1
     CASE("sp") 
        ORBITAL_LIST(1,I) = 0
        ORBITAL_LIST(2,I) = 1
        ORBITAL_LIST(3,I) = -1
     CASE("sd") 
        ORBITAL_LIST(1,I) = 0
        ORBITAL_LIST(2,I) = 2
        ORBITAL_LIST(3,I) = -1
     CASE("sf") 
        ORBITAL_LIST(1,I) = 0
        ORBITAL_LIST(2,I) = 3
        ORBITAL_LIST(3,I) = -1
     CASE("pd") 
        ORBITAL_LIST(1,I) = 1
        ORBITAL_LIST(2,I) = 2
        ORBITAL_LIST(3,I) = -1
     CASE("pf") 
        ORBITAL_LIST(1,I) = 1
        ORBITAL_LIST(2,I) = 3
        ORBITAL_LIST(3,I) = -1
     CASE("df") 
        ORBITAL_LIST(1,I) = 2
        ORBITAL_LIST(2,I) = 3
        ORBITAL_LIST(3,I) = -1
     CASE("spd") 
        ORBITAL_LIST(1,I) = 0
        ORBITAL_LIST(2,I) = 1
        ORBITAL_LIST(3,I) = 2
        ORBITAL_LIST(4,I) = -1
     CASE("spf") 
        ORBITAL_LIST(1,I) = 0
        ORBITAL_LIST(2,I) = 1
        ORBITAL_LIST(3,I) = 3
        ORBITAL_LIST(4,I) = -1
     CASE("sdf") 
        ORBITAL_LIST(1,I) = 0
        ORBITAL_LIST(2,I) = 2
        ORBITAL_LIST(3,I) = 3
        ORBITAL_LIST(4,I) = -1
     CASE("pdf") 
        ORBITAL_LIST(1,I) = 1
        ORBITAL_LIST(2,I) = 2
        ORBITAL_LIST(3,I) = 3
        ORBITAL_LIST(4,I) = -1
     CASE("spdf") 
        ORBITAL_LIST(1,I) = 0
        ORBITAL_LIST(2,I) = 1
        ORBITAL_LIST(3,I) = 2
        ORBITAL_LIST(4,I) = 3
        ORBITAL_LIST(5,I) = -1
     END SELECT
     
  ENDDO

  RETURN
  
END SUBROUTINE GENORBITALLIST
