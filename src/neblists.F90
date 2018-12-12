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
  INTEGER :: A, B, C, MYATOMI, MYATOMJ
  INTEGER :: PBCX, PBCY, PBCZ
  INTEGER :: BOXX, BOXY, BOXZ
  INTEGER :: II, JJ, KK
  INTEGER, INTENT(IN) :: AMIALLO
  INTEGER :: XRANGE, YRANGE, ZRANGE
  INTEGER :: MAXNEBTB, MAXNEBPP, MAXNEBCOUL
  INTEGER :: NCELL(3), NUMCELL
  INTEGER :: IPIV(3), INFO
  INTEGER :: MYCELL, BOXID(3), COUNT
  !  INTEGER, SAVE :: ALLOCEST
  INTEGER, ALLOCATABLE :: TOTINCELL(:), CELLLIST(:,:)
  !  INTEGER, ALLOCATABLE :: DIMTB(:), DIMPP(:), DIMCOUL(:)
  REAL(LATTEPREC) :: RIJ(3), MAGR2
  REAL(LATTEPREC) :: MAGA(3)
  REAL(LATTEPREC) :: RCUTTB, RCUTCOUL, PPMAX, MAXCUT2
  REAL(LATTEPREC) :: WORK(3), BOXINV(3,3), S(3)
  REAL(LATTEPREC), PARAMETER :: MINR = 0.01
  IF (EXISTERROR) RETURN


  IF (PBCON .EQ. 1) CALL PBC

  TOTNEBTB = 0
  IF (PPOTON .GT. 0) TOTNEBPP = 0
  IF (ELECTRO .EQ. 1) TOTNEBCOUL = 0

  IF (AMIALLO .NE. 0) THEN

     ! Reallocate the neighbor lists based on their size the last time
     DEALLOCATE(NEBTB)
     ALLOCATE(NEBTB( 4, MAXDIMTB, NATS ))

     IF (PPOTON .NE. 0) THEN
        DEALLOCATE(NEBPP)
        ALLOCATE(NEBPP( 4, MAXDIMPP, NATS ))
     ENDIF

     IF (ELECTRO .EQ. 1) THEN
        DEALLOCATE(NEBCOUL)
        ALLOCATE(NEBCOUL( 4, MAXDIMCOUL, NATS))
     ENDIF

  ELSE

     ! This bit is only done on the first neighbor list build

     ! Let's get the cut-offs for our interactions

     RCUTTB = ZERO
     PPMAX = ZERO

     ! Find the maximum cut off

     DO K = 1, NOINT

        IF (HCUT(K) .GT. RCUTTB ) RCUTTB = HCUT(K)

        IF (BASISTYPE .EQ. "NONORTHO") THEN
           IF (SCUT(K) .GT. RCUTTB ) RCUTTB = SCUT(K)
        ENDIF

     ENDDO

     IF (PPOTON .GT. 0) THEN

        DO K = 1, NOPPS

           IF (PPOTON .EQ. 1 .AND. POTCOEF(10,K) .GT. PPMAX ) PPMAX = POTCOEF(10,K)

           IF (PPOTON .EQ. 2 .AND. PPR(PPTABLENGTH(K), K) .GT. PPMAX) &
                PPMAX = PPR(PPTABLENGTH(K), K)

           IF (PPOTON .EQ. 3 .AND. PPRK(1,K) .GT. PPMAX) PPMAX = PPRK(1,K)

        ENDDO

     ENDIF

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

     MAXCUT2 = MAXCUT*MAXCUT

     ! Now let's estimate the size of the arrays we need for to
     ! store the neighbor lists, plus some

     IF (PBCON .EQ. 1) THEN

        XRANGE = INT(MAXCUT/BOX(1,1)) + 1
        YRANGE = INT(MAXCUT/BOX(2,2)) + 1
        ZRANGE = INT(MAXCUT/BOX(3,3)) + 1
        !        print*, maxcut, xrange, yrange, zrange

        ! Here we're hoping atom 1 is in a typical environment

        COUNT = 0
        DO J = 1, NATS
           DO II = -XRANGE, XRANGE
              DO JJ = -YRANGE, YRANGE
                 DO KK = -ZRANGE, ZRANGE

                    RIJ(1) = CR(1,J) + REAL(II)*BOX(1,1) + &
                         REAL(JJ)*BOX(2,1) + REAL(KK)*BOX(3,1) - CR(1,1)

                    RIJ(2) = CR(2,J) + REAL(II)*BOX(1,2) + &
                         REAL(JJ)*BOX(2,2) + REAL(KK)*BOX(3,2) - CR(2,1)

                    RIJ(3) = CR(3,J) + REAL(II)*BOX(1,3) + &
                         REAL(JJ)*BOX(2,3) + REAL(KK)*BOX(3,3) - CR(3,1)

                    MAGR2 = RIJ(1)*RIJ(1) + RIJ(2)*RIJ(2) + RIJ(3)*RIJ(3)

                    IF (MAGR2 .LE. MAXCUT2) COUNT = COUNT + 1

                 ENDDO
              ENDDO
           ENDDO
        ENDDO

        MAXDIMTB = 2*COUNT
        MAXDIMPP = 2*COUNT
        MAXDIMCOUL = 2*COUNT

     ELSEIF (PBCON .EQ. 0) THEN

        DIMLIST = 0
        DO I = 1, NATS
           COUNT = 0
           DO J = 1, NATS

              RIJ(1) = CR(1,J) - CR(1,I)
              RIJ(2) = CR(2,J) - CR(2,I)
              RIJ(3) = CR(3,J) - CR(3,I)

              MAGR2 = RIJ(1)*RIJ(1) + RIJ(2)*RIJ(2) + RIJ(3)*RIJ(3)

              IF (MAGR2 .LE. MAXCUT2) COUNT = COUNT + 1

           ENDDO

           IF (COUNT .GT. DIMLIST) DIMLIST = COUNT

        ENDDO

        MAXDIMTB = DIMLIST
        MAXDIMPP = DIMLIST
        MAXDIMCOUL = DIMLIST

     ENDIF

     IF( ALLOCATED(NEBTB) ) DEALLOCATE(NEBTB)
     ALLOCATE ( NEBTB( 4, MAXDIMTB, NATS ) )

     IF (PPOTON .GT. 0) THEN
        IF( ALLOCATED(NEBPP) ) DEALLOCATE(NEBPP)
        ALLOCATE(NEBPP( 4, MAXDIMPP, NATS ) )
     ENDIF

     IF (ELECTRO .EQ. 1) THEN
        IF( ALLOCATED(NEBCOUL) ) DEALLOCATE(NEBCOUL)
        ALLOCATE(NEBCOUL( 4, MAXDIMCOUL, NATS ))
     ENDIF

  ENDIF

  ! Now build the neighbor list

  ! With periodic boundaries first:

  IF (PBCON .EQ. 1) THEN

     MAGA(1) = SQRT(BOX(1,1)*BOX(1,1) + BOX(1,2)*BOX(1,2) + BOX(1,3)*BOX(1,3))
     MAGA(2) = SQRT(BOX(2,1)*BOX(2,1) + BOX(2,2)*BOX(2,2) + BOX(2,3)*BOX(2,3))
     MAGA(3) = SQRT(BOX(3,1)*BOX(3,1) + BOX(3,2)*BOX(3,2) + BOX(3,3)*BOX(3,3))

     XRANGE = INT(MAXCUT/MAGA(1)) + 1
     YRANGE = INT(MAXCUT/MAGA(2)) + 1
     ZRANGE = INT(MAXCUT/MAGA(3)) + 1

     ! This gives the number of sub-cells along each lattice vector

     NCELL(1) = MAX(INT(MAGA(1)/MAXCUT),1)
     NCELL(2) = MAX(INT(MAGA(2)/MAXCUT),1)
     NCELL(3) = MAX(INT(MAGA(3)/MAXCUT),1)

     NUMCELL = NCELL(1)*NCELL(2)*NCELL(3)

     !     PRINT*, NCELL(1), NCELL(2), NCELL(3), NUMCELL

     IF (AMIALLO .EQ. 0) ALLOCEST = 2*NATS/NUMCELL

     ALLOCATE(TOTINCELL(NUMCELL), CELLLIST(ALLOCEST, NUMCELL))

     TOTINCELL = 0

     BOXINV = BOX

     CALL DGETRF(3, 3, BOXINV, 3, IPIV, INFO)

     CALL DGETRI(3, BOXINV, 3, IPIV, WORK, 3, INFO)

     ! Put the atoms into the sub-cells

     DO I = 1, NATS

        CALL DGEMV('T', 3, 3, ONE, BOXINV, 3, CR(1,I), 1, ZERO, S, 1)

        ! Dangerous condition caught below (MJC)

        IF (S(1) .GE. ONE) S(1) = ZERO
        IF (S(2) .GE. ONE) S(2) = ZERO
        IF (S(3) .GE. ONE) S(3) = ZERO

        BOXID(1) = INT(S(1)*NCELL(1))
        BOXID(2) = INT(S(2)*NCELL(2))
        BOXID(3) = INT(S(3)*NCELL(3))

        MYCELL = BOXID(1) + NCELL(1)*BOXID(2) + NCELL(1)*NCELL(2)*BOXID(3) + 1

        TOTINCELL(MYCELL) = TOTINCELL(MYCELL) + 1
        CELLLIST(TOTINCELL(MYCELL), MYCELL) = I

     ENDDO

     ALLOCEST = 2*MAXVAL(TOTINCELL)

     ! Loop over the subcells and build the lists

     DO II = 1, NUMCELL

        ! Indices of subcell II

        BOXZ = (II - 1)/(NCELL(1)*NCELL(2))
        BOXY = (II - 1 - BOXZ*NCELL(1)*NCELL(2))/NCELL(1)
        BOXX = (II - 1 - NCELL(1)*BOXY - NCELL(1)*NCELL(2)*BOXZ)
        !        print*, boxx, boxy, boxz, totincell(ii), celllist(1,ii), NUMCELL

        ! Loop over atoms in cell II

        DO I = 1, TOTINCELL(II)

           MYATOMI = CELLLIST(I,II)

           ! Loop over the neighboring subcells

           DO A = BOXZ - ZRANGE, BOXZ + ZRANGE
              DO B = BOXY - YRANGE, BOXY + YRANGE
                 DO C = BOXX - XRANGE, BOXX + XRANGE

                    PBCX = 0
                    PBCY = 0
                    PBCZ = 0

                    BOXID(1) = C
                    BOXID(2) = B
                    BOXID(3) = A

                    IF (A .LT. 0) THEN
                       PBCZ = A
                       BOXID(3) = NCELL(3) - 1
                    ELSEIF (A .GE. NCELL(3)) THEN
                       PBCZ = A - NCELL(3) + 1
                       BOXID(3) = 0
                    ENDIF

                    IF (B .LT. 0) THEN
                       PBCY = B
                       BOXID(2) = NCELL(2) - 1
                    ELSEIF (B .GE. NCELL(2)) THEN
                       PBCY = B - NCELL(2) + 1
                       BOXID(2) = 0
                    ENDIF

                    IF (C .LT. 0) THEN
                       PBCX = C
                       BOXID(1) = NCELL(1) - 1
                    ELSEIF (C .GE. NCELL(1)) THEN
                       PBCX = C - NCELL(1) + 1
                       BOXID(1) = 0
                    ENDIF

                    MYCELL = BOXID(1) + NCELL(1)*BOXID(2) + NCELL(1)*NCELL(2)*BOXID(3) + 1

                    ! Loop over the atoms in the neighboring cell

                    DO J = 1, TOTINCELL(MYCELL)

                       MYATOMJ = CELLLIST(J, MYCELL)

                       RIJ(1) = CR(1,MYATOMJ) + REAL(PBCX)*BOX(1,1) + &
                            REAL(PBCY)*BOX(2,1) + REAL(PBCZ)*BOX(3,1) - CR(1,MYATOMI)

                       RIJ(2) = CR(2,MYATOMJ) + REAL(PBCX)*BOX(1,2) + &
                            REAL(PBCY)*BOX(2,2) + REAL(PBCZ)*BOX(3,2) - CR(2,MYATOMI)

                       RIJ(3) = CR(3,MYATOMJ) + REAL(PBCX)*BOX(1,3) + &
                            REAL(PBCY)*BOX(2,3) + REAL(PBCZ)*BOX(3,3) - CR(3,MYATOMI)

                       MAGR2 = RIJ(1)*RIJ(1) + RIJ(2)*RIJ(2) + RIJ(3)*RIJ(3)

                       IF (MAGR2 .GT. MINR .AND. MAGR2 .LT. RCUTTB2) THEN

                          TOTNEBTB(MYATOMI) = TOTNEBTB(MYATOMI) + 1

                          IF (TOTNEBTB(MYATOMI) .GT. MAXDIMTB) THEN
                             CALL ERRORS("neblists","NUMBER OF NEIGHBORS EXCEEDS ARRAY DIMENSION (TB)")
                             RETURN
                          ENDIF

                          NEBTB( 1, TOTNEBTB(MYATOMI), MYATOMI ) = MYATOMJ
                          NEBTB( 2, TOTNEBTB(MYATOMI), MYATOMI ) = PBCX
                          NEBTB( 3, TOTNEBTB(MYATOMI), MYATOMI ) = PBCY
                          NEBTB( 4, TOTNEBTB(MYATOMI), MYATOMI ) = PBCZ

                       ENDIF

                       IF (PPOTON .NE. 0) THEN

                          IF (MAGR2 .GT. MINR .AND. MAGR2 .LT. PPMAX2) THEN

                             TOTNEBPP(MYATOMI) = TOTNEBPP(MYATOMI) + 1
                             IF (TOTNEBPP(MYATOMI) .GT. MAXDIMPP) THEN
                                CALL ERRORS("neblists","NUMBER OF NEIGHBORS EXCEEDS ARRAY DIMENSION (PP)")
                                RETURN
                             ENDIF

                             NEBPP( 1, TOTNEBPP(MYATOMI), MYATOMI ) = MYATOMJ
                             NEBPP( 2, TOTNEBPP(MYATOMI), MYATOMI ) = PBCX
                             NEBPP( 3, TOTNEBPP(MYATOMI), MYATOMI ) = PBCY
                             NEBPP( 4, TOTNEBPP(MYATOMI), MYATOMI ) = PBCZ

                          ENDIF

                       ENDIF

                       IF (ELECTRO .NE. 0) THEN

                          IF (MAGR2 .GT. MINR .AND. MAGR2 .LT. RCUTCOUL2) THEN

                             TOTNEBCOUL(MYATOMI) = TOTNEBCOUL(MYATOMI) + 1

                             IF (TOTNEBCOUL(MYATOMI) .GT. MAXDIMCOUL) THEN
                                CALL ERRORS("neblists","NUMBER OF NEIGHBORS EXCEEDS ARRAY DIMENSION (COUL)")
                                RETURN
                             ENDIF

                             NEBCOUL( 1, TOTNEBCOUL(MYATOMI), MYATOMI ) = MYATOMJ
                             NEBCOUL( 2, TOTNEBCOUL(MYATOMI), MYATOMI ) = PBCX
                             NEBCOUL( 3, TOTNEBCOUL(MYATOMI), MYATOMI ) = PBCY
                             NEBCOUL( 4, TOTNEBCOUL(MYATOMI), MYATOMI ) = PBCZ

                          ENDIF

                       ENDIF


                    ENDDO
                 ENDDO
              ENDDO

           ENDDO

        ENDDO

     ENDDO

     DEALLOCATE(TOTINCELL, CELLLIST)

  ELSEIF (PBCON .EQ. 0) THEN

     ! Now we're doing building the neighbor lists for gas-phase systems

     DO I = 1, NATS
        DO J = 1, NATS

           RIJ(1) = CR(1,J) - CR(1,I)

           RIJ(2) = CR(2,J) - CR(2,I)

           RIJ(3) = CR(3,J) - CR(3,I)

           MAGR2 = RIJ(1)*RIJ(1) + RIJ(2)*RIJ(2) + RIJ(3)*RIJ(3)

           IF (MAGR2 .GT. MINR .AND. MAGR2 .LT. RCUTTB2) THEN

              TOTNEBTB(I) = TOTNEBTB(I) + 1
              NEBTB( 1, TOTNEBTB(I), I ) = J
              NEBTB( 2, TOTNEBTB(I), I ) = 0
              NEBTB( 3, TOTNEBTB(I), I ) = 0
              NEBTB( 4, TOTNEBTB(I), I ) = 0

           ENDIF

           IF (PPOTON .NE. 0) THEN

              IF (MAGR2 .GT. MINR .AND. MAGR2 .LT. PPMAX2) THEN

                 TOTNEBPP(I) = TOTNEBPP(I) + 1
                 NEBPP( 1, TOTNEBPP(I), I ) = J
                 NEBPP( 2, TOTNEBPP(I), I ) = 0
                 NEBPP( 3, TOTNEBPP(I), I ) = 0
                 NEBPP( 4, TOTNEBPP(I), I ) = 0

              ENDIF

           ENDIF

           IF (ELECTRO .NE. 0) THEN

              IF (MAGR2 .GT. MINR .AND. MAGR2 .LT. RCUTCOUL2) THEN

                 TOTNEBCOUL(I) = TOTNEBCOUL(I) + 1
                 NEBCOUL( 1, TOTNEBCOUL(I), I ) = J
                 NEBCOUL( 2, TOTNEBCOUL(I), I ) = 0
                 NEBCOUL( 3, TOTNEBCOUL(I), I ) = 0
                 NEBCOUL( 4, TOTNEBCOUL(I), I ) = 0

              ENDIF

           ENDIF

        ENDDO
     ENDDO

  ENDIF

  ! Let's get the dimensions of the arrays about right for the next
  ! loop through here

  MAXDIMTB = 2*MAXVAL(TOTNEBTB)
  IF (PPOTON .NE. 0) MAXDIMPP = 2*MAXVAL(TOTNEBPP)
  IF (ELECTRO .NE. 0) MAXDIMCOUL = 2*MAXVAL(TOTNEBCOUL)

  RETURN

END SUBROUTINE NEBLISTS
