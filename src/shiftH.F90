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

SUBROUTINE SHIFTH(CHI)

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE COULOMBARRAY
  USE NONOARRAY
  USE KSPACEARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, K
  INTEGER :: INDEX, NUMORB
  REAL(LATTEPREC) :: ES, EP, ED, EF
  REAL(LATTEPREC) :: HMOD, CHI
  COMPLEX(LATTEPREC) :: ZHMOD
  

  INDEX = 0

  IF (KON .EQ. 0) THEN

     IF ( BASISTYPE .EQ. "ORTHO") THEN
        
        DO I = 1, NATS
           
           ! Using a constant
           
           LCNSHIFT(I) = LCNSHIFT(I) + CHI * DELTAQ(I)
           
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
              H(INDEX, INDEX) = HDIAG(INDEX) + LCNSHIFT(I)
           ENDDO
 
        ENDDO
        
     ELSEIF ( BASISTYPE .EQ. "NONORTHO" ) THEN
        
        DO I = 1, NATS
           
           LCNSHIFT(I) = LCNSHIFT(I) + CHI * DELTAQ(I)
           
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
              HJJ(INDEX) = LCNSHIFT(I)
!              H(INDEX, INDEX) = HDIAG(INDEX) + LCNSHIFT(I)
                           
           ENDDO
           
        ENDDO
        
        DO I = 1, HDIM
           DO J = 1, HDIM

              H(J,I) = H0(J,I) + SMAT(J,I)*(HJJ(I) + HJJ(J))/TWO

           ENDDO
        ENDDO

        
     ENDIF
     
  ELSE

     !  IF (KON .EQ. 1) THEN 
     ! k-space - we have to add the potential to all NKTOT Hamiltonians
     
     ! Orthogonal basis only at the moment
     
     IF ( BASISTYPE .EQ. "ORTHO") THEN
        
        DO I = 1, NATS           

!           IF (CONTROL .EQ. 1) THEN

              !Using the response function calculated in GETRESPF               

!              LCNSHIFT(I) = LCNSHIFT(I) + RESPCHI(I)*DELTAQ(I)

!           ELSE

              ! Using a constant                                                

              LCNSHIFT(I) = LCNSHIFT(I) + CHI * DELTAQ(I)

!           ENDIf
           
        ENDDO

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
                 HK(INDEX, INDEX, K) = HKDIAG(INDEX, K) + CMPLX(LCNSHIFT(I))
              ENDDO
           
           ENDDO

        ENDDO

     ELSEIF (BASISTYPE .EQ. "NONORTHO") THEN

         DO I = 1, NATS

           LCNSHIFT(I) = LCNSHIFT(I) + CHI * DELTAQ(I)
 
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
              ZHJJ(INDEX) = CMPLX(LCNSHIFT(I))

           ENDDO

        ENDDO
        
        ! H = H_0 + S*H_1                                                     \
                                                                               

        DO K = 1, NKTOT
           DO I = 1, HDIM
              DO J = 1, HDIM

                 HK(J,I,K) = HK0(J,I,K) + SK(J,I,K)*(ZHJJ(I) + ZHJJ(J))/TWO

              ENDDO
           ENDDO
        ENDDO
        
     ENDIF

        
  ENDIF

  IF (DEBUGON .EQ. 1) THEN
     
     OPEN(UNIT=31, STATUS="UNKNOWN", FILE="myH.dat")
     
     PRINT*, "Caution - the H0+H1 matrix is being written to file"
     
     IF (KON .EQ. 0) THEN
        
        DO I = 1, HDIM
           WRITE(31,10) (H(I,J), J = 1, HDIM)
        ENDDO
        
     ELSE

        DO K = 1, NKTOT
           WRITE(31,*) K
           DO I = 1, HDIM
              WRITE(31,10) (HK(I,J,K), J = 1, HDIM)
           ENDDO
        ENDDO

     ENDIF

     CLOSE(31)
     
  ENDIF
  
!  PRINT*, LCNSHIFT(1)


10 FORMAT(100G18.9)  
11 FORMAT(I5,100G18.9)  
  
  RETURN

END SUBROUTINE SHIFTH

     
