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

SUBROUTINE NEBLISTS(AMIALLO)

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE UNIVARRAY
  USE PPOTARRAY
  USE NEBLISTARRAY
  USE COULOMBARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, K, L, M
  INTEGER :: II, JJ, KK
  INTEGER, INTENT(IN) :: AMIALLO
  INTEGER :: XRANGE, YRANGE, ZRANGE
  INTEGER :: MAXNEBTB, MAXNEBPP, MAXNEBCOUL
  INTEGER, ALLOCATABLE :: DIMTB(:), DIMPP(:), DIMCOUL(:)
  REAL(LATTEPREC) :: RIJ(3), MAGR2
  REAL(LATTEPREC) :: MAGA, MAGB, MAGC
  REAL(LATTEPREC) :: RCUTTB, RCUTCOUL, PPMAX
  REAL(LATTEPREC), PARAMETER :: MINR = 0.01

  IF (PBCON .EQ. 1) CALL PBC

  TOTNEBTB = 0
  IF (PPOTON .GT. 0) TOTNEBPP = 0
  IF (ELECTRO .EQ. 1) TOTNEBCOUL = 0


  IF (AMIALLO .NE. 0) THEN

     !
     ! Allow the arrays to grow. It's pointless to let them shrink.
     !
     
     IF (MAXDIMTB .GT. PREVDIMTB) THEN
        DEALLOCATE(NEBTB)
        ALLOCATE(NEBTB( 4, MAXDIMTB, NATS ))
     ENDIF

     IF (PPOTON .GT. 0) THEN
        IF (MAXDIMPP .GT. PREVDIMPP) THEN
           DEALLOCATE(NEBPP)
           ALLOCATE(NEBPP( 4, MAXDIMPP, NATS ))
        ENDIF
     ENDIF

     IF (ELECTRO .EQ. 1) THEN
        IF (MAXDIMCOUL .GT. PREVDIMCOUL) THEN
           DEALLOCATE(NEBCOUL)
           ALLOCATE(NEBCOUL( 4, MAXDIMCOUL, NATS))
        ENDIF
     ENDIF

  ELSE

     ! This bit is only done on the first neighbor list build

     ! Let's get the cut-offs for our interactions

     PREVNEBTB = 0
     PREVNEBPP = 0
     PREVNEBCOUL = 0

     RCUTTB = ZERO
     PPMAX = ZERO

     DO I = 1, NATS
        DO J = I, NATS

           DO K = 1, NOINT

              IF ( (ATELE(I) .EQ. ELE1(K) .AND. &
                   ATELE(J) .EQ. ELE2(K)) .OR. & 
                   (ATELE(J) .EQ. ELE1(K) .AND. &
                   ATELE(I) .EQ. ELE2(K) )) THEN
                 
                 IF (BOND(8,K) .GT. RCUTTB ) RCUTTB = BOND(8,K)

                 IF (BASISTYPE .EQ. "NONORTHO") THEN
                    IF (OVERL(8,K) .GT. RCUTTB ) RCUTTB = OVERL(8,K)
                 ENDIF
                 
              ENDIF
              
           ENDDO
           
           IF (PPOTON .GT. 0) THEN

              DO K = 1, NOPPS
                 
                 IF ( (ATELE(I) .EQ. PPELE1(K) .AND. &
                      ATELE(J) .EQ. PPELE2(K)) .OR. & 
                      (ATELE(J) .EQ. PPELE1(K) .AND. &
                      ATELE(I) .EQ. PPELE2(K)) ) THEN
                    
                    IF (PPOTON .EQ. 1 .AND. POTCOEF(10,K) .GT. PPMAX ) PPMAX = POTCOEF(10,K)
                    
                    IF (PPOTON .EQ. 2 .AND. PPR(PPTABLENGTH(K), K) .GT. PPMAX) &
                         PPMAX = PPR(PPTABLENGTH(K), K)


                 ENDIF

              ENDDO
           
           ENDIF
        ENDDO
     ENDDO

!     print*, "RCUTTB = ", RCUTTB
!     print*, "RCUTPP = ", PPMAX

     RCUTTB = RCUTTB + SKIN
     RCUTTB2 = RCUTTB*RCUTTB
     
     IF (PPOTON .GT. 0) THEN
        PPMAX = PPMAX + SKIN
        PPMAX2 = PPMAX*PPMAX
     ELSE
        PPMAX = ZERO
        PPMAX2 = ZERO
     ENDIF
        
     RCUTCOUL = COULCUT + SKIN
     RCUTCOUL2 = RCUTCOUL * RCUTCOUL

     IF (ELECTRO .EQ. 0) RCUTCOUL = ZERO

     MAXCUT = MAX(RCUTTB, PPMAX, RCUTCOUL)

!     PRINT*, RCUTTB, PPMAX, RCUTCOUL

     ! Now let's estimate the size of the arrays we need for to 
     ! store the neighbor lists, plus some

     

     IF (PBCON .EQ. 1) THEN

        MAGA = SQRT(BOX(1,1)*BOX(1,1) + BOX(1,2)*BOX(1,2) + BOX(1,3)*BOX(1,3))
        XRANGE = INT( MAXCUT/MAGA ) + 1
        MAGB = SQRT(BOX(2,1)*BOX(2,1) + BOX(2,2)*BOX(2,2) + BOX(2,3)*BOX(2,3))
        YRANGE = INT( MAXCUT/MAGB ) + 1
        MAGC = SQRT(BOX(3,1)*BOX(3,1) + BOX(3,2)*BOX(3,2) + BOX(3,3)*BOX(3,3))
        ZRANGE = INT( MAXCUT/MAGC ) + 1

     ELSE
        XRANGE = 0
        YRANGE = 0
        ZRANGE = 0
     ENDIF

     ! FOR ATOM 1
     
     MAXDIMTB = 0
     MAXDIMPP = 0
     MAXDIMCOUL = 0

     ALLOCATE(DIMTB(NATS), DIMPP(NATS), DIMCOUL(NATS))
     DIMTB = 0
     DIMPP = 0
     DIMCOUL = 0
     
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(NATS, XRANGE, YRANGE, ZRANGE, BOX, CR, RCUTTB2, PPMAX2, RCUTCOUL2) &
!$OMP SHARED(DIMTB, DIMPP, DIMCOUL, PPOTON, ELECTRO) &
!$OMP PRIVATE(I, J, II, JJ, KK,  RIJ, MAGR2)     

     DO I = 1, NATS        
        DO J = 1, NATS
           DO II = -XRANGE, XRANGE
              DO JJ = -YRANGE, YRANGE
                 DO KK = -ZRANGE, ZRANGE
                    
                    RIJ(1) = CR(1,J) + REAL(II)*BOX(1,1) + &
                         REAL(JJ)*BOX(2,1) + REAL(KK)*BOX(3,1) - CR(1,I)

                    RIJ(2) = CR(2,J) + REAL(II)*BOX(1,2) + &
                         REAL(JJ)*BOX(2,2) + REAL(KK)*BOX(3,2) - CR(2,I)

                    RIJ(3) = CR(3,J) + REAL(II)*BOX(1,3) + &
                         REAL(JJ)*BOX(2,3) + REAL(KK)*BOX(3,3) - CR(3,I)
                    
                    MAGR2 = RIJ(1)*RIJ(1) + RIJ(2)*RIJ(2) + RIJ(3)*RIJ(3)
                    
                    
                    IF (MAGR2 .GT. MINR) THEN

                       IF (MAGR2 .LT. RCUTTB2) DIMTB(I) = DIMTB(I) + 1
                       IF (PPOTON .GT. 0 .AND. &
                            MAGR2 .LT. PPMAX2) DIMPP(I) = DIMPP(I) + 1
                       IF (ELECTRO .EQ. 1 .AND. &
                            MAGR2 .LT. RCUTCOUL2) DIMCOUL(I) = DIMCOUL(I) + 1
                       
                    ENDIF
                 
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO

!$OMP END PARALLEL DO
     
     ! Let's be very conservative the with allocation here

!     PRINT*, MAXDIMTB, MAXDIMPP, MAXDIMCOUL

     MAXDIMTB = INT(REAL(MAXVAL(DIMTB))*1.5)
     MAXDIMPP = INT(REAL(MAXVAL(DIMPP))*1.5)
     MAXDIMCOUL = INT(REAL(MAXVAL(DIMCOUL))*1.5)

!     PRINT*, MAXDIMTB, MAXDIMPP, MAXDIMCOUL

     DEALLOCATE(DIMTB, DIMPP, DIMCOUL)

     ALLOCATE ( NEBTB( 4, MAXDIMTB, NATS ) )
     IF (PPOTON .GT. 0) ALLOCATE(NEBPP( 4, MAXDIMPP, NATS ) ) 
     IF (ELECTRO .EQ. 1) ALLOCATE(NEBCOUL( 4, MAXDIMCOUL, NATS ))
     
  ENDIF

  IF (PBCON .EQ. 1) THEN
     
     MAGA = SQRT(BOX(1,1)*BOX(1,1) + BOX(1,2)*BOX(1,2) + BOX(1,3)*BOX(1,3))
     XRANGE = INT( MAXCUT/MAGA ) + 1
     MAGB = SQRT(BOX(2,1)*BOX(2,1) + BOX(2,2)*BOX(2,2) + BOX(2,3)*BOX(2,3))
     YRANGE = INT( MAXCUT/MAGB ) + 1
     MAGC = SQRT(BOX(3,1)*BOX(3,1) + BOX(3,2)*BOX(3,2) + BOX(3,3)*BOX(3,3))
     ZRANGE = INT( MAXCUT/MAGC ) + 1

  ELSE
     XRANGE = 0
     YRANGE = 0
     ZRANGE = 0
  ENDIF

!  PRINT*, "RANGE = ", XRANGE, YRANGE, ZRANGE

!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP SHARED(NATS, XRANGE, YRANGE, ZRANGE, BOX, CR, RCUTTB2, PPMAX2, RCUTCOUL2) &
!$OMP SHARED(TOTNEBTB, TOTNEBPP, TOTNEBCOUL, NEBTB, NEBPP, NEBCOUL, PPOTON, ELECTRO, MAXDIMTB, MAXDIMPP, MAXDIMCOUL) &
!$OMP PRIVATE(I, J, II, JJ, KK,  RIJ, MAGR2)

  DO I = 1, NATS
     DO J = 1, NATS

        ! 
        ! Do the H-matrix neighbor list
        !
        
        ! Looping over the neighboring boxes

        DO II = -XRANGE, XRANGE

           DO JJ = -YRANGE, YRANGE

              DO KK = -ZRANGE, ZRANGE

                 RIJ(1) = CR(1,J) + REAL(II)*BOX(1,1) + &
                      REAL(JJ)*BOX(2,1) + REAL(KK)*BOX(3,1) - CR(1,I)
                 
                 RIJ(2) = CR(2,J) + REAL(II)*BOX(1,2) + &
                      REAL(JJ)*BOX(2,2) + REAL(KK)*BOX(3,2) - CR(2,I)
                 
                 RIJ(3) = CR(3,J) + REAL(II)*BOX(1,3) + &
                      REAL(JJ)*BOX(2,3) + REAL(KK)*BOX(3,3) - CR(3,I)
                 
                 
                 MAGR2 = RIJ(1)*RIJ(1) + RIJ(2)*RIJ(2) + RIJ(3)*RIJ(3)
                 
                 IF (MAGR2 .GT. MINR) THEN
                    
                    ! List for the Hamiltonian build

                    IF (MAGR2 .LT. RCUTTB2) THEN
                    
                       TOTNEBTB(I) = TOTNEBTB(I) + 1
                       
                       IF (TOTNEBTB(I) .GT. MAXDIMTB) THEN
                          PRINT*, "NUMBER OF NEIGHBORS EXCEEDS ARRAY DIMENSION (TB)"
                          STOP
                       ENDIF
                          
                       NEBTB( 1, TOTNEBTB(I), I ) = J
                       NEBTB( 2, TOTNEBTB(I), I ) = II
                       NEBTB( 3, TOTNEBTB(I), I ) = JJ
                       NEBTB( 4, TOTNEBTB(I), I ) = KK
                    
                    ENDIF
                    
                    !
                    ! Now the neighbor list for the pair potential
                    !

                    IF (PPOTON .GT. 0 ) THEN
                       
                       IF (MAGR2 .LT. PPMAX2) THEN
                             
                          TOTNEBPP(I) = TOTNEBPP(I) + 1
                          
                          IF (TOTNEBPP(I) .GT. MAXDIMPP) THEN
                             PRINT*, "NUMBER OF NEIGHBORS EXCEEDS ARRAY DIMENSION (PP)"
                             STOP
                          ENDIF
                                                    
                          NEBPP( 1, TOTNEBPP(I), I ) = J 
                          NEBPP( 2, TOTNEBPP(I), I ) = II 
                          NEBPP( 3, TOTNEBPP(I), I ) = JJ
                          NEBPP( 4, TOTNEBPP(I), I ) = KK
                          
                       ENDIF
                       
                    ENDIF
                       
                 !   ENDIF

                 

                    ! List for the coulomb sum

                    IF (ELECTRO .EQ. 1) THEN
                       
                       IF (MAGR2 .LE. RCUTCOUL2) THEN
                          
                          TOTNEBCOUL(I) = TOTNEBCOUL(I) + 1
                          
                          IF (TOTNEBCOUL(I) .GT. MAXDIMCOUL) THEN
                             PRINT*, "NUMBER OF NEIGHBORS EXCEEDS ARRAY DIMENSION (COUL)"
                             STOP
                          ENDIF

                          NEBCOUL( 1, TOTNEBCOUL(I), I ) = J
                          NEBCOUL( 2, TOTNEBCOUL(I), I ) = II
                          NEBCOUL( 3, TOTNEBCOUL(I), I ) = JJ
                          NEBCOUL( 4, TOTNEBCOUL(I), I ) = KK
                          
                       ENDIF
                       
                    ENDIF

                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
!$OMP END PARALLEL DO


  ! Let's get the dimensions of the arrays about right for the next
  ! loop through here

  MAXNEBTB = MAXVAL(TOTNEBTB)
  IF (PPOTON .GT. 0) MAXNEBPP = MAXVAL(TOTNEBPP) 
  IF (ELECTRO .EQ. 1) MAXNEBCOUL = MAXVAL(TOTNEBCOUL)
  
  PREVDIMTB = MAXDIMTB
  PREVDIMPP = MAXDIMPP
  PREVDIMCOUL = MAXDIMCOUL

  ! If we have more neighbors this time around increase the allocation

  ! Allocate more storage to be safe

  IF (MAXNEBTB .GT. PREVNEBTB) THEN
     PREVNEBTB = MAXNEBTB
     MAXDIMTB = INT(REAL(MAXNEBTB)*1.5)
  ENDIF

  IF (PPOTON .GT. 0 .AND. MAXNEBPP .GT. PREVNEBPP) THEN
     PREVNEBPP = MAXNEBPP
     MAXDIMPP = INT(REAL(MAXNEBPP)*1.5)
  ENDIF

  IF (ELECTRO .EQ. 1 .AND. MAXNEBCOUL .GT. PREVNEBCOUL) THEN
     PREVNEBCOUL = MAXNEBCOUL
     MAXDIMCOUL = INT(REAL(MAXNEBCOUL)*1.5) 
  ENDIF

  RETURN

END SUBROUTINE NEBLISTS
