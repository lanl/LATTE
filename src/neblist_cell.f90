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
  INTEGER :: AMIALLO
  INTEGER :: XRANGE, YRANGE, ZRANGE
  INTEGER :: MAXNEBTB, MAXNEBPP, MAXNEBCOUL
  INTEGER :: NX, NY, NZ, XCELL, YCELL, ZCELL, NUMCELLS
  INTEGER :: PBCI, PBCJ, PBCK, ATOMJ, CELLID, NEB
  INTEGER, ALLOCATABLE :: TOTPERCELL(:), CELLATOMID(:, :, :)
  REAL(LATTEPREC) :: X, Y, Z, R2
  REAL(LATTEPREC) :: MINCRX, MINCRY, MINCRZ
  REAL(LATTEPREC), PARAMETER :: MINR = 0.01
  REAL(LATTEPREC), SAVE :: MAXCUTTB, MAXCUTPP, MAXCUTCOUL

  IF (PBCON .EQ. 1) CALL PBC

  ! Ensure all coordinates are > 0.0 when building the list.
  ! We will shift them back when we're done

  MINCRX = MINVAL(CR(1,:))
  MINCRY = MINVAL(CR(2,:))
  MINCRZ = MINVAL(CR(3,:))

  DO I = 1, NATS
     CR(1,I) = CR(1,I) - MINCRX
     CR(2,I) = CR(2,I) - MINCRY
     CR(3,I) = CR(3,I) - MINCRZ
  ENDDO

  ! First pass: set up cut-offs

  PRINT*, "1"
  IF (AMIALLO .EQ. 0) THEN

     MAXCUTTB = ZERO
     MAXCUTPP = ZERO

     !$OMP PARALLEL DO DEFAULT(NONE) &
     !$OMP SHARED(NATS, NOINT, ATELE, ELE1, ELE2, BOND, OVERL) &
     !$OMP SHARED(BASISTYPE, PPOTON, NOPPS, POTCOEF, PPELE1, PPELE2) &
     !$OMP PRIVATE(I, J, K) &
     !$OMP REDUCTION(MAX: MAXCUTTB, MAXCUTPP)

     DO I = 1, NATS
        DO J = I, NATS

           DO K = 1, NOINT

              IF ( (ATELE(I) .EQ. ELE1(K) .AND. &
                   ATELE(J) .EQ. ELE2(K)) .OR. &
                   (ATELE(J) .EQ. ELE1(K) .AND. &
                   ATELE(I) .EQ. ELE2(K) )) THEN

                 MAXCUTTB = MAX(BOND(8,K), MAXCUTTB)
                 IF (BASISTYPE .EQ. "NONORTHO") &
                      MAXCUTTB = MAX(OVERL(8,K), MAXCUTTB)

              ENDIF

           ENDDO

           IF (PPOTON .EQ. 1) THEN

              DO K = 1, NOPPS

                 IF ( (ATELE(I) .EQ. PPELE1(K) .AND. &
                      ATELE(J) .EQ. PPELE2(K)) .OR. &
                      (ATELE(J) .EQ. PPELE1(K) .AND. &
                      ATELE(I) .EQ. PPELE2(K)) ) THEN

                    MAXCUTPP = MAX(POTCOEF(10,K), MAXCUTPP)

                 ENDIF

              ENDDO

           ENDIF
        ENDDO
     ENDDO

     !$OMP END PARALLEL DO

     PRINT*, "2"
     MAXCUTTB = MAXCUTTB + SKIN
     MAXCUTPP = MAXCUTPP + SKIN
     MAXCUTCOUL = COULCUT + SKIN

     ! We will figure out good values for these array dimensions later

     MAXDIMTB = 0
     MAXDIMPP = 0
     MAXDIMCOUL = 0

  ELSE

     !
     ! Allow the arrays to grow. It's pointless to let them shrink.
     !

     IF (MAXDIMTB .GT. PREVDIMTB) THEN
        DEALLOCATE(NEBTB)
        ALLOCATE(NEBTB( 4, MAXDIMTB, NATS ))
     ENDIF

     IF (PPOTON .EQ. 1) THEN
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

  ENDIF
  PRINT*, "3"
  ! TB

  IF (PBCON .EQ. 1) THEN

     XRANGE = INT(MAXCUTTB/BOXDIMS(1)) + 1
     YRANGE = INT(MAXCUTTB/BOXDIMS(2)) + 1
     ZRANGE = INT(MAXCUTTB/BOXDIMS(3)) + 1

     NX = INT(BOXDIMS(1)*(ONE + REAL(XRANGE))/MAXCUTTB) + 1
     NY = INT(BOXDIMS(2)*(ONE + REAL(YRANGE))/MAXCUTTB) + 1
     NZ = INT(BOXDIMS(3)*(ONE + REAL(ZRANGE))/MAXCUTTB) + 1

     NUMCELLS = (NX+1) * (NY+1) * (NZ+1)

  ELSE

     XRANGE = 0
     YRANGE = 0
     ZRANGE = 0

     NX = INT(MAXVAL(CR(1,:))/MAXCUTTB) + 1
     NY = INT(MAXVAL(CR(2,:))/MAXCUTTB) + 1
     NZ = INT(MAXVAL(CR(3,:))/MAXCUTTB) + 1

     NUMCELLS = NX * NY * NZ

  ENDIF

  ALLOCATE(TOTPERCELL(NUMCELLS), CELLATOMID(4, 200, NUMCELLS))

  ! Put all the atoms in their cell

  TOTPERCELL = 0

  PRINT*, "4"

  IF (PBCON .EQ. 1) THEN

     DO II = -XRANGE, XRANGE
        DO JJ = -YRANGE, YRANGE
           DO KK = -ZRANGE, ZRANGE

              DO I = 1, NATS

                 IF (CR(1,I) + REAL(II)*BOXDIMS(1) .LT. ZERO .AND. &
                      CR(1,I) + REAL(II)*BOXDIMS(1) .GT. -MAXCUTTB) THEN
                    XCELL = 0
                 ELSEIF (CR(1,I) + REAL(II)*BOXDIMS(1) .GE. ZERO) THEN
                    XCELL = INT((CR(1,I) + REAL(II)*BOXDIMS(1))/MAXCUTTB) + 1
                 ELSE
                    XCELL = -1
                 ENDIF

                 IF (CR(2,I) + REAL(JJ)*BOXDIMS(2) .LT. ZERO .AND. &
                      CR(2,I) + REAL(JJ)*BOXDIMS(2) .GT. -MAXCUTTB) THEN
                    YCELL = 0
                 ELSEIF (CR(2,I) + REAL(JJ)*BOXDIMS(2) .GE. ZERO) THEN
                    YCELL = INT((CR(2,I) + REAL(JJ)*BOXDIMS(2))/MAXCUTTB) + 1
                 ELSE
                    YCELL = -1
                 ENDIF

                 IF (CR(3,I) + REAL(KK)*BOXDIMS(3) .LT. ZERO .AND. &
                      CR(3,I) + REAL(KK)*BOXDIMS(3) .GT. -MAXCUTTB) THEN
                    ZCELL = 0
                 ELSEIF (CR(3,I) + REAL(KK)*BOXDIMS(3) .GE. ZERO) THEN
                    ZCELL = INT((CR(3,I) + REAL(KK)*BOXDIMS(3))/MAXCUTTB) + 1
                 ELSE
                    ZCELL = -1
                 ENDIF

                 IF (XCELL .GE. 0 .AND. YCELL .GE. 0 .AND. ZCELL .GE. 0 &
                      .AND. XCELL .LE. NX .AND. YCELL .LE. NY .AND. &
                      ZCELL .LE. NZ) THEN

                    CELLID = XCELL + (NX+1)*(YCELL-1) + (NX+1)*(NY+1)*(ZCELL-1)
                    CELLID = CELLID + (NX+1) + (NX+1)*(NY+1) + 1

                    TOTPERCELL(CELLID) = TOTPERCELL(CELLID) + 1
                    CELLATOMID(1, TOTPERCELL(CELLID), CELLID) = I
                    CELLATOMID(2, TOTPERCELL(CELLID), CELLID) = II
                    CELLATOMID(3, TOTPERCELL(CELLID), CELLID) = JJ
                    CELLATOMID(4, TOTPERCELL(CELLID), CELLID) = KK


                 ENDIF

              ENDDO
           ENDDO
        ENDDO
     ENDDO

  ELSE

     DO I = 1, NATS

        XCELL = INT(CR(1,I)/MAXCUTTB) + 1
        YCELL = INT(CR(2,I)/MAXCUTTB) + 1
        ZCELL = INT(CR(3,I)/MAXCUTTB) + 1

        CELLID = XCELL + NX*(YCELL-1) + NY*NX*(ZCELL-1)

        TOTPERCELL(CELLID) = TOTPERCELL(CELLID) + 1
        CELLATOMID(1, TOTPERCELL(CELLID), CELLID) = I
        CELLATOMID(2, TOTPERCELL(CELLID), CELLID) = 0
        CELLATOMID(3, TOTPERCELL(CELLID), CELLID) = 0
        CELLATOMID(4, TOTPERCELL(CELLID), CELLID) = 0

     ENDDO

  ENDIF

  PRINT*, "5"

  ! On the first pass through we can use the cell list
  ! to figure out how to dimension the neighbor list arrays

  IF (AMIALLO .EQ. 0) THEN

     IF (PBCON .EQ. 1) THEN

        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP SHARED(NATS, BOXDIMS, CR, TOTPERCELL, CELLATOMID) &
        !$OMP SHARED(NX, NY, NZ, MAXCUTTB) &
        !$OMP PRIVATE(I, XCELL, YCELL, ZCELL, II, JJ, KK, CELLID, J, ATOMJ) &
        !$OMP PRIVATE(PBCI, PBCJ, PBCK, X, Y, Z, R2, NEB)  &
        !$OMP REDUCTION(MAX: MAXDIMTB)

        DO I = 1, NATS

           NEB = 0

           ! Find its cell

           XCELL = INT(CR(1,I)/MAXCUTTB) + 1
           YCELL = INT(CR(2,I)/MAXCUTTB) + 1
           ZCELL = INT(CR(3,I)/MAXCUTTB) + 1

           ! Loop over all neighboring cells including its own cell

           DO KK = ZCELL - 1, ZCELL + 1
              DO JJ = YCELL - 1, YCELL + 1
                 DO II = XCELL - 1, XCELL + 1

                    CELLID = II + (NX+1)*(JJ-1) + (NX+1)*(NY+1)*(KK-1)
                    CELLID = CELLID + (NX+1) + (NX+1)*(NY+1) + 1

                    DO J = 1, TOTPERCELL(CELLID)

                       ATOMJ = CELLATOMID(1, J, CELLID)
                       PBCI = CELLATOMID(2, J, CELLID)
                       PBCJ = CELLATOMID(3, J, CELLID)
                       PBCK = CELLATOMID(4, J, CELLID)

                       X = CR(1,ATOMJ) + REAL(PBCI)*BOXDIMS(1) - CR(1,I)
                       Y = CR(2,ATOMJ) + REAL(PBCJ)*BOXDIMS(2) - CR(2,I)
                       Z = CR(3,ATOMJ) + REAL(PBCK)*BOXDIMS(3) - CR(3,I)

                       R2 = X*X + Y*Y + Z*Z

                       IF (R2 .LT. MAXCUTTB*MAXCUTTB .AND. &
                            R2 .GT. MINR) NEB = NEB + 1

                    ENDDO

                 ENDDO
              ENDDO
           ENDDO

           MAXDIMTB = MAX(NEB, MAXDIMTB)

        ENDDO

        !$OMP END PARALLEL DO

     ELSE ! NO PBC

        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP SHARED(NATS, BOXDIMS, CR, TOTPERCELL, CELLATOMID) &
        !$OMP SHARED(NX, NY, NZ, MAXCUTTB) &
        !$OMP PRIVATE(I, XCELL, YCELL, ZCELL, II, JJ, KK, CELLID, J, ATOMJ) &
        !$OMP PRIVATE(X, Y, Z, R2, NEB)  &
        !$OMP REDUCTION(MAX: MAXDIMTB)

        DO I = 1, NATS

           NEB = 0

           XCELL = INT(CR(1,I)/MAXCUTTB) + 1
           YCELL = INT(CR(2,I)/MAXCUTTB) + 1
           ZCELL = INT(CR(3,I)/MAXCUTTB) + 1

           DO KK = ZCELL - 1, ZCELL + 1
              DO JJ = YCELL - 1, YCELL + 1
                 DO II = XCELL - 1, XCELL + 1

                    IF (II .GE. 1 .AND. JJ .GE. 1 .AND. KK .GE. 1 &
                         .AND. II .LE. NX .AND. JJ .LE. NY .AND. &
                         KK .LE. NZ) THEN


                       CELLID = II + NX*(JJ-1) + NY*NX*(KK-1)

                       DO J = 1, TOTPERCELL(CELLID)

                          ATOMJ = CELLATOMID(1, J, CELLID)

                          X = CR(1,ATOMJ) - CR(1,I)
                          Y = CR(2,ATOMJ) - CR(2,I)
                          Z = CR(3,ATOMJ) - CR(3,I)

                          R2 = X*X + Y*Y + Z*Z

                          IF (R2 .LT. MAXCUTTB*MAXCUTTB .AND. &
                               R2 .GT. MINR) NEB = NEB + 1

                       ENDDO
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO

           MAXDIMTB = MAX(NEB, MAXDIMTB)

        ENDDO

        !$OMP END PARALLEL DO

     ENDIF

     MAXDIMTB = INT(1.5*REAL(MAXDIMTB))

     ALLOCATE( NEBTB(4, MAXDIMTB, NATS) )

  ENDIF
  PRINT*, "6"


  ! Start with atom 1

  TOTNEBTB = 0

  IF (PBCON .EQ. 1) THEN

     !$OMP PARALLEL DO DEFAULT(NONE) &
     !$OMP SHARED(NATS, BOXDIMS, CR, MAXCUTTB, TOTPERCELL, CELLATOMID) &
     !$OMP SHARED(NEBTB, TOTNEBTB, NX, NY, NZ, MAXDIMTB, PBCON) &
     !$OMP PRIVATE(I, XCELL, YCELL, ZCELL, II, JJ, KK, CELLID, J, ATOMJ) &
     !$OMP PRIVATE(PBCI, PBCJ, PBCK, X, Y, Z, R2)

     DO I = 1, NATS

        ! Find its cell

        XCELL = INT(CR(1,I)/MAXCUTTB) + 1
        YCELL = INT(CR(2,I)/MAXCUTTB) + 1
        ZCELL = INT(CR(3,I)/MAXCUTTB) + 1

        !     PRINT*, I, XCELL, YCELL, ZCELL

        ! Loop over all neighboring cells including its own cell

        DO KK = ZCELL - 1, ZCELL + 1
           DO JJ = YCELL - 1, YCELL + 1
              DO II = XCELL - 1, XCELL + 1

                 CELLID = II + (NX+1)*(JJ-1) + (NX+1)*(NY+1)*(KK-1)
                 CELLID = CELLID + (NX+1) + (NX+1)*(NY+1) + 1

                 DO J = 1, TOTPERCELL(CELLID)

                    ATOMJ = CELLATOMID(1, J, CELLID)
                    PBCI = CELLATOMID(2, J, CELLID)
                    PBCJ = CELLATOMID(3, J, CELLID)
                    PBCK = CELLATOMID(4, J, CELLID)

                    X = CR(1,ATOMJ) + REAL(PBCI)*BOXDIMS(1) - CR(1,I)
                    Y = CR(2,ATOMJ) + REAL(PBCJ)*BOXDIMS(2) - CR(2,I)
                    Z = CR(3,ATOMJ) + REAL(PBCK)*BOXDIMS(3) - CR(3,I)

                    R2 = X*X + Y*Y + Z*Z

                    IF (R2 .LT. MAXCUTTB*MAXCUTTB .AND. R2 .GT. MINR) THEN

                       TOTNEBTB(I) = TOTNEBTB(I) + 1

                       IF (TOTNEBTB(I) .GT. MAXDIMTB) THEN
                          CALL ERRORS("neblist_cell","NUMBER OF NEIGHBORS EXCEEDS ARRAY DIMENSION (TB)")
                       ENDIF

                       NEBTB( 1, TOTNEBTB(I), I ) = ATOMJ
                       NEBTB( 2, TOTNEBTB(I), I ) = PBCI
                       NEBTB( 3, TOTNEBTB(I), I ) = PBCJ
                       NEBTB( 4, TOTNEBTB(I), I ) = PBCK

                    ENDIF

                 ENDDO

              ENDDO
           ENDDO
        ENDDO
     ENDDO

     !$OMP END PARALLEL DO

  ELSE ! NO PBC

     !$OMP PARALLEL DO DEFAULT(NONE) &
     !$OMP SHARED(NATS, BOXDIMS, CR, MAXCUTTB, TOTPERCELL, CELLATOMID) &
     !$OMP SHARED(NEBTB, TOTNEBTB, NX, NY, NZ, MAXDIMTB) &
     !$OMP PRIVATE(I, XCELL, YCELL, ZCELL, II, JJ, KK, CELLID, J, ATOMJ) &
     !$OMP PRIVATE(PBCI, PBCJ, PBCK, X, Y, Z, R2)

     DO I = 1, NATS

        XCELL = INT(CR(1,I)/MAXCUTTB) + 1
        YCELL = INT(CR(2,I)/MAXCUTTB) + 1
        ZCELL = INT(CR(3,I)/MAXCUTTB) + 1

        DO KK = ZCELL - 1, ZCELL + 1
           DO JJ = YCELL - 1, YCELL + 1
              DO II = XCELL - 1, XCELL + 1

                 IF (II .GE. 1 .AND. JJ .GE. 1 .AND. KK .GE. 1 &
                      .AND. II .LE. NX .AND. JJ .LE. NY .AND. &
                      KK .LE. NZ) THEN


                    CELLID = II + NX*(JJ-1) + NY*NX*(KK-1)

                    DO J = 1, TOTPERCELL(CELLID)

                       ATOMJ = CELLATOMID(1, J, CELLID)

                       X = CR(1,ATOMJ) - CR(1,I)
                       Y = CR(2,ATOMJ) - CR(2,I)
                       Z = CR(3,ATOMJ) - CR(3,I)

                       R2 = X*X + Y*Y + Z*Z

                       IF (R2 .LT. MAXCUTTB*MAXCUTTB .AND. R2 .GT. MINR) THEN

                          TOTNEBTB(I) = TOTNEBTB(I) + 1

                          IF (TOTNEBTB(I) .GT. MAXDIMTB) THEN
                             CALL ERRORS("neblist_cell","NUMBER OF NEIGHBORS EXCEEDS ARRAY DIMENSION (TB)")
                          ENDIF

                          NEBTB( 1, TOTNEBTB(I), I ) = ATOMJ
                          NEBTB( 2, TOTNEBTB(I), I ) = 0
                          NEBTB( 3, TOTNEBTB(I), I ) = 0
                          NEBTB( 4, TOTNEBTB(I), I ) = 0

                       ENDIF

                    ENDDO
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDDO

     !$OMP END PARALLEL DO

  ENDIF
  PRINT*, "7"
  ! print*, "here"
  DEALLOCATE(TOTPERCELL, CELLATOMID)

  ! Now the list for the pair potential

  IF (PPOTON .EQ. 1) THEN

     IF (PBCON .EQ. 1) THEN

        XRANGE = INT(MAXCUTPP/BOXDIMS(1)) + 1
        YRANGE = INT(MAXCUTPP/BOXDIMS(2)) + 1
        ZRANGE = INT(MAXCUTPP/BOXDIMS(3)) + 1

        NX = INT(BOXDIMS(1)*(ONE + REAL(XRANGE))/MAXCUTPP) + 1
        NY = INT(BOXDIMS(2)*(ONE + REAL(YRANGE))/MAXCUTPP) + 1
        NZ = INT(BOXDIMS(3)*(ONE + REAL(ZRANGE))/MAXCUTPP) + 1

        NUMCELLS = (NX+1) * (NY+1) * (NZ+1)

     ELSE

        XRANGE = 0
        YRANGE = 0
        ZRANGE = 0

        NX = INT(MAXVAL(CR(1,:))/MAXCUTPP) + 1
        NY = INT(MAXVAL(CR(2,:))/MAXCUTPP) + 1
        NZ = INT(MAXVAL(CR(3,:))/MAXCUTPP) + 1

        NUMCELLS = NX * NY * NZ

     ENDIF

     ALLOCATE(TOTPERCELL(NUMCELLS), CELLATOMID(4, NATS, NUMCELLS))

     ! Put all the atoms in their cell

     TOTPERCELL = 0

     IF (PBCON .EQ. 1) THEN

        DO II = -XRANGE, XRANGE
           DO JJ = -YRANGE, YRANGE
              DO KK = -ZRANGE, ZRANGE

                 DO I = 1, NATS

                    IF (CR(1,I) + REAL(II)*BOXDIMS(1) .LT. ZERO .AND. &
                         CR(1,I) + REAL(II)*BOXDIMS(1) .GT. -MAXCUTPP) THEN
                       XCELL = 0
                    ELSEIF (CR(1,I) + REAL(II)*BOXDIMS(1) .GE. ZERO) THEN
                       XCELL = INT((CR(1,I) + REAL(II)*BOXDIMS(1))/MAXCUTPP) + 1
                    ELSE
                       XCELL = -1
                    ENDIF

                    IF (CR(2,I) + REAL(JJ)*BOXDIMS(2) .LT. ZERO .AND. &
                         CR(2,I) + REAL(JJ)*BOXDIMS(2) .GT. -MAXCUTPP) THEN
                       YCELL = 0
                    ELSEIF (CR(2,I) + REAL(JJ)*BOXDIMS(2) .GE. ZERO) THEN
                       YCELL = INT((CR(2,I) + REAL(JJ)*BOXDIMS(2))/MAXCUTPP) + 1
                    ELSE
                       YCELL = -1
                    ENDIF

                    IF (CR(3,I) + REAL(KK)*BOXDIMS(3) .LT. ZERO .AND. &
                         CR(3,I) + REAL(KK)*BOXDIMS(3) .GT. -MAXCUTPP) THEN
                       ZCELL = 0
                    ELSEIF (CR(3,I) + REAL(KK)*BOXDIMS(3) .GE. ZERO) THEN
                       ZCELL = INT((CR(3,I) + REAL(KK)*BOXDIMS(3))/MAXCUTPP) + 1
                    ELSE
                       ZCELL = -1
                    ENDIF

                    IF (XCELL .GE. 0 .AND. YCELL .GE. 0 .AND. ZCELL .GE. 0 &
                         .AND. XCELL .LE. NX .AND. YCELL .LE. NY .AND. &
                         ZCELL .LE. NZ) THEN

                       CELLID = XCELL + (NX+1)*(YCELL-1) + (NX+1)*(NY+1)*(ZCELL-1)
                       CELLID = CELLID + (NX+1) + (NX+1)*(NY+1) + 1

                       TOTPERCELL(CELLID) = TOTPERCELL(CELLID) + 1
                       CELLATOMID(1, TOTPERCELL(CELLID), CELLID) = I
                       CELLATOMID(2, TOTPERCELL(CELLID), CELLID) = II
                       CELLATOMID(3, TOTPERCELL(CELLID), CELLID) = JJ
                       CELLATOMID(4, TOTPERCELL(CELLID), CELLID) = KK


                    ENDIF

                 ENDDO
              ENDDO
           ENDDO
        ENDDO

     ELSE

        DO I = 1, NATS

           XCELL = INT(CR(1,I)/MAXCUTPP) + 1
           YCELL = INT(CR(2,I)/MAXCUTPP) + 1
           ZCELL = INT(CR(3,I)/MAXCUTPP) + 1

           CELLID = XCELL + NX*(YCELL-1) + NY*NX*(ZCELL-1)

           TOTPERCELL(CELLID) = TOTPERCELL(CELLID) + 1
           CELLATOMID(1, TOTPERCELL(CELLID), CELLID) = I
           CELLATOMID(2, TOTPERCELL(CELLID), CELLID) = 0
           CELLATOMID(3, TOTPERCELL(CELLID), CELLID) = 0
           CELLATOMID(4, TOTPERCELL(CELLID), CELLID) = 0

        ENDDO

     ENDIF

     PRINT *, "8"
     !PRINT*, "MAX CELL = ", MAXVAL(TOTPERCELL)
     ! First pass through - figure out dimension of neighborlist array

     IF (AMIALLO .EQ. 0) THEN

        IF (PBCON .EQ. 1) THEN

           PRINT *, "MAX CELL =", MAXVAL(TOTPERCELL)

!!$OMP PARALLEL DO DEFAULT(NONE) &
!!$OMP SHARED(NATS, BOXDIMS, CR, TOTPERCELL, CELLATOMID) &
!!$OMP SHARED(NX, NY, NZ, MAXCUTPP) &
!!$OMP PRIVATE(I, XCELL, YCELL, ZCELL, II, JJ, KK, CELLID, J, ATOMJ) &
!!$OMP PRIVATE(PBCI, PBCJ, PBCK, X, Y, Z, R2, NEB)  &
!!$OMP REDUCTION(MAX: MAXDIMPP)

           DO I = 1, NATS

              NEB = 0

              ! Find its cell

              XCELL = INT(CR(1,I)/MAXCUTPP) + 1
              YCELL = INT(CR(2,I)/MAXCUTPP) + 1
              ZCELL = INT(CR(3,I)/MAXCUTPP) + 1

              ! Loop over all neighboring cells including its own cell

              DO KK = ZCELL - 1, ZCELL + 1
                 DO JJ = YCELL - 1, YCELL + 1
                    DO II = XCELL - 1, XCELL + 1

                       CELLID = II + (NX+1)*(JJ-1) + (NX+1)*(NY+1)*(KK-1)
                       CELLID = CELLID + (NX+1) + (NX+1)*(NY+1) + 1

                       DO J = 1, TOTPERCELL(CELLID)

                          ATOMJ = CELLATOMID(1, J, CELLID)
                          PBCI = CELLATOMID(2, J, CELLID)
                          PBCJ = CELLATOMID(3, J, CELLID)
                          PBCK = CELLATOMID(4, J, CELLID)

                          X = CR(1,ATOMJ) + REAL(PBCI)*BOXDIMS(1) - CR(1,I)
                          Y = CR(2,ATOMJ) + REAL(PBCJ)*BOXDIMS(2) - CR(2,I)
                          Z = CR(3,ATOMJ) + REAL(PBCK)*BOXDIMS(3) - CR(3,I)

                          R2 = X*X + Y*Y + Z*Z

                          IF (R2 .LT. MAXCUTPP*MAXCUTPP .AND. &
                               R2 .GT. MINR) NEB = NEB + 1

                       ENDDO

                    ENDDO
                 ENDDO
              ENDDO

              MAXDIMPP = MAX(NEB, MAXDIMPP)

           ENDDO

!!$OMP END PARALLEL DO

        ELSE ! NO PBC

!!$OMP PARALLEL DO DEFAULT(NONE) &
!!$OMP SHARED(NATS, BOXDIMS, CR, TOTPERCELL, CELLATOMID) &
!!$OMP SHARED(NX, NY, NZ, MAXCUTPP) &
!!$OMP PRIVATE(I, XCELL, YCELL, ZCELL, II, JJ, KK, CELLID, J, ATOMJ) &
!!$OMP PRIVATE(X, Y, Z, R2, NEB)  &
!!$OMP REDUCTION(MAX: MAXDIMPP)

           DO I = 1, NATS

              NEB = 0

              XCELL = INT(CR(1,I)/MAXCUTPP) + 1
              YCELL = INT(CR(2,I)/MAXCUTPP) + 1
              ZCELL = INT(CR(3,I)/MAXCUTPP) + 1

              DO KK = ZCELL - 1, ZCELL + 1
                 DO JJ = YCELL - 1, YCELL + 1
                    DO II = XCELL - 1, XCELL + 1

                       IF (II .GE. 1 .AND. JJ .GE. 1 .AND. KK .GE. 1 &
                            .AND. II .LE. NX .AND. JJ .LE. NY .AND. &
                            KK .LE. NZ) THEN


                          CELLID = II + NX*(JJ-1) + NY*NX*(KK-1)

                          DO J = 1, TOTPERCELL(CELLID)

                             ATOMJ = CELLATOMID(1, J, CELLID)

                             X = CR(1,ATOMJ) - CR(1,I)
                             Y = CR(2,ATOMJ) - CR(2,I)
                             Z = CR(3,ATOMJ) - CR(3,I)

                             R2 = X*X + Y*Y + Z*Z

                             IF (R2 .LT. MAXCUTPP*MAXCUTPP .AND. &
                                  R2 .GT. MINR) NEB = NEB + 1

                          ENDDO
                       ENDIF
                    ENDDO
                 ENDDO
              ENDDO

              MAXDIMPP = MAX(NEB, MAXDIMPP)

           ENDDO

!!$OMP END PARALLEL DO

        ENDIF

        MAXDIMPP = INT(1.5*REAL(MAXDIMPP))

        ALLOCATE( NEBPP(4, MAXDIMPP, NATS) )

     ENDIF

     PRINT*, "maxdimpp = ", MAXDIMPP
     PRINT*, "9"

     ! Start with atom 1

     TOTNEBPP = 0

     IF (PBCON .EQ. 1) THEN

!!$OMP PARALLEL DO DEFAULT(NONE) &
!!$OMP SHARED(NATS, BOXDIMS, CR, MAXCUTPP, TOTPERCELL, CELLATOMID) &
!!$OMP SHARED(NEBPP, TOTNEBPP, NX, NY, NZ, MAXDIMPP) &
!!$OMP PRIVATE(I, XCELL, YCELL, ZCELL, II, JJ, KK, CELLID, J, ATOMJ) &
!!$OMP PRIVATE(PBCI, PBCJ, PBCK, X, Y, Z, R2)

        DO I = 1, NATS

           ! Find its cell

           XCELL = INT(CR(1,I)/MAXCUTPP) + 1
           YCELL = INT(CR(2,I)/MAXCUTPP) + 1
           ZCELL = INT(CR(3,I)/MAXCUTPP) + 1

           ! Loop over all neighboring cells including its own cell

           DO KK = ZCELL - 1, ZCELL + 1
              DO JJ = YCELL - 1, YCELL + 1
                 DO II = XCELL - 1, XCELL + 1

                    CELLID = II + (NX+1)*(JJ-1) + (NX+1)*(NY+1)*(KK-1)
                    CELLID = CELLID + (NX+1) + (NX+1)*(NY+1) + 1

                    DO J = 1, TOTPERCELL(CELLID)

                       ATOMJ = CELLATOMID(1, J, CELLID)
                       PBCI = CELLATOMID(2, J, CELLID)
                       PBCJ = CELLATOMID(3, J, CELLID)
                       PBCK = CELLATOMID(4, J, CELLID)

                       X = CR(1,ATOMJ) + REAL(PBCI)*BOXDIMS(1) - CR(1,I)
                       Y = CR(2,ATOMJ) + REAL(PBCJ)*BOXDIMS(2) - CR(2,I)
                       Z = CR(3,ATOMJ) + REAL(PBCK)*BOXDIMS(3) - CR(3,I)

                       R2 = X*X + Y*Y + Z*Z

                       IF (R2 .LT. MAXCUTPP*MAXCUTPP .AND. R2 .GT. MINR) THEN

                          TOTNEBPP(I) = TOTNEBPP(I) + 1

                          IF (TOTNEBPP(I) .GT. MAXDIMPP) THEN
                             CALL ERRORS("neblist_cell","NUMBER OF NEIGHBORS EXCEEDS ARRAY DIMENSION (PP)")
                          ENDIF

                          NEBPP( 1, TOTNEBPP(I), I ) = ATOMJ
                          NEBPP( 2, TOTNEBPP(I), I ) = PBCI
                          NEBPP( 3, TOTNEBPP(I), I ) = PBCJ
                          NEBPP( 4, TOTNEBPP(I), I ) = PBCK

                       ENDIF

                    ENDDO

                 ENDDO
              ENDDO
           ENDDO
        ENDDO

!!$OMP END PARALLEL DO

     ELSE ! NO PBC

        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP SHARED(NATS, BOXDIMS, CR, MAXCUTPP, TOTPERCELL, CELLATOMID) &
        !$OMP SHARED(NEBPP, TOTNEBPP, NX, NY, NZ, MAXDIMPP) &
        !$OMP PRIVATE(I, XCELL, YCELL, ZCELL, II, JJ, KK, CELLID, J, ATOMJ) &
        !$OMP PRIVATE(PBCI, PBCJ, PBCK, X, Y, Z, R2)

        DO I = 1, NATS

           XCELL = INT(CR(1,I)/MAXCUTPP) + 1
           YCELL = INT(CR(2,I)/MAXCUTPP) + 1
           ZCELL = INT(CR(3,I)/MAXCUTPP) + 1

           DO KK = ZCELL - 1, ZCELL + 1
              DO JJ = YCELL - 1, YCELL + 1
                 DO II = XCELL - 1, XCELL + 1

                    IF (II .GE. 1 .AND. JJ .GE. 1 .AND. KK .GE. 1 &
                         .AND. II .LE. NX .AND. JJ .LE. NY .AND. &
                         KK .LE. NZ) THEN


                       CELLID = II + NX*(JJ-1) + NY*NX*(KK-1)

                       DO J = 1, TOTPERCELL(CELLID)

                          ATOMJ = CELLATOMID(1, J, CELLID)

                          X = CR(1,ATOMJ) - CR(1,I)
                          Y = CR(2,ATOMJ) - CR(2,I)
                          Z = CR(3,ATOMJ) - CR(3,I)

                          R2 = X*X + Y*Y + Z*Z

                          IF (R2 .LT. MAXCUTPP*MAXCUTPP .AND. R2 .GT. MINR) THEN

                             TOTNEBPP(I) = TOTNEBPP(I) + 1

                             IF (TOTNEBPP(I) .GT. MAXDIMPP) THEN
                                CALL ERRORS("neblist_cell","NUMBER OF NEIGHBORS EXCEEDS ARRAY DIMENSION (PP)")
                             ENDIF

                             NEBPP( 1, TOTNEBPP(I), I ) = ATOMJ
                             NEBPP( 2, TOTNEBPP(I), I ) = 0
                             NEBPP( 3, TOTNEBPP(I), I ) = 0
                             NEBPP( 4, TOTNEBPP(I), I ) = 0

                          ENDIF

                       ENDDO
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO
        ENDDO

        !$OMP END PARALLEL DO

     ENDIF

     DEALLOCATE(TOTPERCELL, CELLATOMID)

  ENDIF
  PRINT*, "9"
  ! Coulomb now...

  IF (ELECTRO .EQ. 1) THEN

     IF (PBCON .EQ. 1) THEN

        XRANGE = INT(MAXCUTCOUL/BOXDIMS(1)) + 1
        YRANGE = INT(MAXCUTCOUL/BOXDIMS(2)) + 1
        ZRANGE = INT(MAXCUTCOUL/BOXDIMS(3)) + 1

        NX = INT(BOXDIMS(1)*(ONE + REAL(XRANGE))/MAXCUTCOUL) + 1
        NY = INT(BOXDIMS(2)*(ONE + REAL(YRANGE))/MAXCUTCOUL) + 1
        NZ = INT(BOXDIMS(3)*(ONE + REAL(ZRANGE))/MAXCUTCOUL) + 1

        NUMCELLS = (NX+1) * (NY+1) * (NZ+1)

     ELSE

        XRANGE = 0
        YRANGE = 0
        ZRANGE = 0

        NX = INT(MAXVAL(CR(1,:))/MAXCUTCOUL) + 1
        NY = INT(MAXVAL(CR(2,:))/MAXCUTCOUL) + 1
        NZ = INT(MAXVAL(CR(3,:))/MAXCUTCOUL) + 1

        NUMCELLS = NX * NY * NZ
        PRINT*, NUMCELLS
     ENDIF

     ALLOCATE(TOTPERCELL(NUMCELLS), CELLATOMID(4, NATS, NUMCELLS))

     ! Put all the atoms in their cell

     TOTPERCELL = 0

     IF (PBCON .EQ. 1) THEN

        DO II = -XRANGE, XRANGE
           DO JJ = -YRANGE, YRANGE
              DO KK = -ZRANGE, ZRANGE

                 DO I = 1, NATS

                    IF (CR(1,I) + REAL(II)*BOXDIMS(1) .LT. ZERO .AND. &
                         CR(1,I) + REAL(II)*BOXDIMS(1) .GT. -MAXCUTCOUL) THEN
                       XCELL = 0
                    ELSEIF (CR(1,I) + REAL(II)*BOXDIMS(1) .GE. ZERO) THEN
                       XCELL = INT((CR(1,I) + REAL(II)*BOXDIMS(1))/MAXCUTCOUL) + 1
                    ELSE
                       XCELL = -1
                    ENDIF

                    IF (CR(2,I) + REAL(JJ)*BOXDIMS(2) .LT. ZERO .AND. &
                         CR(2,I) + REAL(JJ)*BOXDIMS(2) .GT. -MAXCUTCOUL) THEN
                       YCELL = 0
                    ELSEIF (CR(2,I) + REAL(JJ)*BOXDIMS(2) .GE. ZERO) THEN
                       YCELL = INT((CR(2,I) + REAL(JJ)*BOXDIMS(2))/MAXCUTCOUL) + 1
                    ELSE
                       YCELL = -1
                    ENDIF

                    IF (CR(3,I) + REAL(KK)*BOXDIMS(3) .LT. ZERO .AND. &
                         CR(3,I) + REAL(KK)*BOXDIMS(3) .GT. -MAXCUTCOUL) THEN
                       ZCELL = 0
                    ELSEIF (CR(3,I) + REAL(KK)*BOXDIMS(3) .GE. ZERO) THEN
                       ZCELL = INT((CR(3,I) + REAL(KK)*BOXDIMS(3))/MAXCUTCOUL) + 1
                    ELSE
                       ZCELL = -1
                    ENDIF

                    IF (XCELL .GE. 0 .AND. YCELL .GE. 0 .AND. ZCELL .GE. 0 &
                         .AND. XCELL .LE. NX .AND. YCELL .LE. NY .AND. &
                         ZCELL .LE. NZ) THEN

                       CELLID = XCELL + (NX+1)*(YCELL-1) + (NX+1)*(NY+1)*(ZCELL-1)
                       CELLID = CELLID + (NX+1) + (NX+1)*(NY+1) + 1

                       TOTPERCELL(CELLID) = TOTPERCELL(CELLID) + 1
                       CELLATOMID(1, TOTPERCELL(CELLID), CELLID) = I
                       CELLATOMID(2, TOTPERCELL(CELLID), CELLID) = II
                       CELLATOMID(3, TOTPERCELL(CELLID), CELLID) = JJ
                       CELLATOMID(4, TOTPERCELL(CELLID), CELLID) = KK


                    ENDIF

                 ENDDO
              ENDDO
           ENDDO
        ENDDO

     ELSE

        DO I = 1, NATS

           XCELL = INT(CR(1,I)/MAXCUTCOUL) + 1
           YCELL = INT(CR(2,I)/MAXCUTCOUL) + 1
           ZCELL = INT(CR(3,I)/MAXCUTCOUL) + 1

           CELLID = XCELL + NX*(YCELL-1) + NY*NX*(ZCELL-1)

           TOTPERCELL(CELLID) = TOTPERCELL(CELLID) + 1
           CELLATOMID(1, TOTPERCELL(CELLID), CELLID) = I
           CELLATOMID(2, TOTPERCELL(CELLID), CELLID) = 0
           CELLATOMID(3, TOTPERCELL(CELLID), CELLID) = 0
           CELLATOMID(4, TOTPERCELL(CELLID), CELLID) = 0

        ENDDO

     ENDIF

     ! First pass through - figure out dimension of neighborlist array

     IF (AMIALLO .EQ. 0) THEN

        IF (PBCON .EQ. 1) THEN

           !$OMP PARALLEL DO DEFAULT(NONE) &
           !$OMP SHARED(NATS, BOXDIMS, CR, TOTPERCELL, CELLATOMID) &
           !$OMP SHARED(NX, NY, NZ, MAXCUTCOUL) &
           !$OMP PRIVATE(I, XCELL, YCELL, ZCELL, II, JJ, KK, CELLID, J, ATOMJ) &
           !$OMP PRIVATE(PBCI, PBCJ, PBCK, X, Y, Z, R2, NEB)  &
           !$OMP REDUCTION(MAX: MAXDIMCOUL)

           DO I = 1, NATS

              NEB = 0

              ! Find its cell

              XCELL = INT(CR(1,I)/MAXCUTCOUL) + 1
              YCELL = INT(CR(2,I)/MAXCUTCOUL) + 1
              ZCELL = INT(CR(3,I)/MAXCUTCOUL) + 1

              ! Loop over all neighboring cells including its own cell

              DO KK = ZCELL - 1, ZCELL + 1
                 DO JJ = YCELL - 1, YCELL + 1
                    DO II = XCELL - 1, XCELL + 1

                       CELLID = II + (NX+1)*(JJ-1) + (NX+1)*(NY+1)*(KK-1)
                       CELLID = CELLID + (NX+1) + (NX+1)*(NY+1) + 1

                       DO J = 1, TOTPERCELL(CELLID)

                          ATOMJ = CELLATOMID(1, J, CELLID)
                          PBCI = CELLATOMID(2, J, CELLID)
                          PBCJ = CELLATOMID(3, J, CELLID)
                          PBCK = CELLATOMID(4, J, CELLID)

                          X = CR(1,ATOMJ) + REAL(PBCI)*BOXDIMS(1) - CR(1,I)
                          Y = CR(2,ATOMJ) + REAL(PBCJ)*BOXDIMS(2) - CR(2,I)
                          Z = CR(3,ATOMJ) + REAL(PBCK)*BOXDIMS(3) - CR(3,I)

                          R2 = X*X + Y*Y + Z*Z

                          IF (R2 .LT. MAXCUTCOUL*MAXCUTCOUL .AND. &
                               R2 .GT. MINR) NEB = NEB + 1

                       ENDDO

                    ENDDO
                 ENDDO
              ENDDO

              MAXDIMCOUL = MAX(NEB, MAXDIMCOUL)

           ENDDO

           !$OMP END PARALLEL DO

        ELSE ! NO PBC

           !$OMP PARALLEL DO DEFAULT(NONE) &
           !$OMP SHARED(NATS, BOXDIMS, CR, TOTPERCELL, CELLATOMID) &
           !$OMP SHARED(NX, NY, NZ, MAXCUTCOUL) &
           !$OMP PRIVATE(I, XCELL, YCELL, ZCELL, II, JJ, KK, CELLID, J, ATOMJ) &
           !$OMP PRIVATE(X, Y, Z, R2, NEB)  &
           !$OMP REDUCTION(MAX: MAXDIMCOUL)

           DO I = 1, NATS

              NEB = 0

              XCELL = INT(CR(1,I)/MAXCUTCOUL) + 1
              YCELL = INT(CR(2,I)/MAXCUTCOUL) + 1
              ZCELL = INT(CR(3,I)/MAXCUTCOUL) + 1

              DO KK = ZCELL - 1, ZCELL + 1
                 DO JJ = YCELL - 1, YCELL + 1
                    DO II = XCELL - 1, XCELL + 1

                       IF (II .GE. 1 .AND. JJ .GE. 1 .AND. KK .GE. 1 &
                            .AND. II .LE. NX .AND. JJ .LE. NY .AND. &
                            KK .LE. NZ) THEN


                          CELLID = II + NX*(JJ-1) + NY*NX*(KK-1)

                          DO J = 1, TOTPERCELL(CELLID)

                             ATOMJ = CELLATOMID(1, J, CELLID)

                             X = CR(1,ATOMJ) - CR(1,I)
                             Y = CR(2,ATOMJ) - CR(2,I)
                             Z = CR(3,ATOMJ) - CR(3,I)

                             R2 = X*X + Y*Y + Z*Z

                             IF (R2 .LT. MAXCUTCOUL*MAXCUTCOUL .AND. &
                                  R2 .GT. MINR) NEB = NEB + 1

                          ENDDO
                       ENDIF
                    ENDDO
                 ENDDO
              ENDDO

              MAXDIMCOUL = MAX(NEB, MAXDIMCOUL)

           ENDDO

           !$OMP END PARALLEL DO

        ENDIF

        MAXDIMCOUL = INT(1.5*REAL(MAXDIMCOUL))

        ALLOCATE( NEBCOUL(4, MAXDIMCOUL, NATS) )

     ENDIF



     ! Start with atom 1

     TOTNEBCOUL = 0

     IF (PBCON .EQ. 1) THEN

        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP SHARED(NATS, BOXDIMS, CR, MAXCUTCOUL, TOTPERCELL, CELLATOMID) &
        !$OMP SHARED(NEBCOUL, TOTNEBCOUL, NX, NY, NZ, MAXDIMCOUL, PBCON) &
        !$OMP PRIVATE(I, XCELL, YCELL, ZCELL, II, JJ, KK, CELLID, J, ATOMJ) &
        !$OMP PRIVATE(PBCI, PBCJ, PBCK, X, Y, Z, R2)

        DO I = 1, NATS

           ! Find its cell

           XCELL = INT(CR(1,I)/MAXCUTCOUL) + 1
           YCELL = INT(CR(2,I)/MAXCUTCOUL) + 1
           ZCELL = INT(CR(3,I)/MAXCUTCOUL) + 1

           ! Loop over all neighboring cells including its own cell

           DO KK = ZCELL - 1, ZCELL + 1
              DO JJ = YCELL - 1, YCELL + 1
                 DO II = XCELL - 1, XCELL + 1

                    CELLID = II + (NX+1)*(JJ-1) + (NX+1)*(NY+1)*(KK-1)
                    CELLID = CELLID + (NX+1) + (NX+1)*(NY+1) + 1

                    DO J = 1, TOTPERCELL(CELLID)

                       ATOMJ = CELLATOMID(1, J, CELLID)
                       PBCI = CELLATOMID(2, J, CELLID)
                       PBCJ = CELLATOMID(3, J, CELLID)
                       PBCK = CELLATOMID(4, J, CELLID)

                       X = CR(1,ATOMJ) + REAL(PBCI)*BOXDIMS(1) - CR(1,I)
                       Y = CR(2,ATOMJ) + REAL(PBCJ)*BOXDIMS(2) - CR(2,I)
                       Z = CR(3,ATOMJ) + REAL(PBCK)*BOXDIMS(3) - CR(3,I)

                       R2 = X*X + Y*Y + Z*Z

                       IF (R2 .LT. MAXCUTCOUL*MAXCUTCOUL .AND. R2 .GT. MINR) THEN

                          TOTNEBCOUL(I) = TOTNEBCOUL(I) + 1

                          IF (TOTNEBCOUL(I) .GT. MAXDIMCOUL) THEN
                             CALL ERRORS("neblist_cell","NUMBER OF NEIGHBORS EXCEEDS ARRAY DIMENSION (COUL)")
                          ENDIF

                          NEBCOUL( 1, TOTNEBCOUL(I), I ) = ATOMJ
                          NEBCOUL( 2, TOTNEBCOUL(I), I ) = PBCI
                          NEBCOUL( 3, TOTNEBCOUL(I), I ) = PBCJ
                          NEBCOUL( 4, TOTNEBCOUL(I), I ) = PBCK

                       ENDIF

                    ENDDO

                 ENDDO
              ENDDO
           ENDDO
        ENDDO

        !$OMP END PARALLEL DO

     ELSE ! NO PBC

        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP SHARED(NATS, BOXDIMS, CR, MAXCUTCOUL, TOTPERCELL, CELLATOMID) &
        !$OMP SHARED(NEBCOUL, TOTNEBCOUL, NX, NY, NZ, MAXDIMCOUL) &
        !$OMP PRIVATE(I, XCELL, YCELL, ZCELL, II, JJ, KK, CELLID, J, ATOMJ) &
        !$OMP PRIVATE(X, Y, Z, R2)

        DO I = 1, NATS

           XCELL = INT(CR(1,I)/MAXCUTCOUL) + 1
           YCELL = INT(CR(2,I)/MAXCUTCOUL) + 1
           ZCELL = INT(CR(3,I)/MAXCUTCOUL) + 1

           DO KK = ZCELL - 1, ZCELL + 1
              DO JJ = YCELL - 1, YCELL + 1
                 DO II = XCELL - 1, XCELL + 1

                    IF (II .GE. 1 .AND. JJ .GE. 1 .AND. KK .GE. 1 &
                         .AND. II .LE. NX .AND. JJ .LE. NY .AND. &
                         KK .LE. NZ) THEN

                       CELLID = II + NX*(JJ-1) + NY*NX*(KK-1)

                       DO J = 1, TOTPERCELL(CELLID)

                          ATOMJ = CELLATOMID(1, J, CELLID)

                          X = CR(1,ATOMJ) - CR(1,I)
                          Y = CR(2,ATOMJ) - CR(2,I)
                          Z = CR(3,ATOMJ) - CR(3,I)

                          R2 = X*X + Y*Y + Z*Z

                          IF (R2 .LT. MAXCUTCOUL*MAXCUTCOUL .AND. R2 .GT. MINR) THEN

                             TOTNEBCOUL(I) = TOTNEBCOUL(I) + 1

                             IF (TOTNEBCOUL(I) .GT. MAXDIMCOUL) THEN
                                CALL ERRORS("neblist_cell","NUMBER OF NEIGHBORS EXCEEDS ARRAY DIMENSION (COUL)")
                             ENDIF

                             NEBCOUL( 1, TOTNEBCOUL(I), I ) = ATOMJ
                             NEBCOUL( 2, TOTNEBCOUL(I), I ) = 0
                             NEBCOUL( 3, TOTNEBCOUL(I), I ) = 0
                             NEBCOUL( 4, TOTNEBCOUL(I), I ) = 0

                          ENDIF

                       ENDDO
                    ENDIF
                 ENDDO
              ENDDO
           ENDDO
        ENDDO

        !$OMP END PARALLEL DO

     ENDIF

     DEALLOCATE(TOTPERCELL, CELLATOMID)

  ENDIF

  DO I = 1, NATS
     CR(1,I) = CR(1,I) + MINCRX
     CR(2,I) = CR(2,I) + MINCRY
     CR(3,I) = CR(3,I) + MINCRZ
  ENDDO

  ! Check whether we need to redimension arrays

  MAXNEBTB = MAXVAL(TOTNEBTB)
  IF (PPOTON .EQ. 1) MAXNEBPP = MAXVAL(TOTNEBPP)
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

  IF (MAXNEBPP .GT. PREVNEBPP) THEN
     PREVNEBPP = MAXNEBPP
     MAXDIMPP = INT(REAL(MAXNEBPP)*1.5)
  ENDIF

  IF (ELECTRO .EQ. 1 .AND. MAXNEBCOUL .GT. PREVNEBCOUL) THEN
     PREVNEBCOUL = MAXNEBCOUL
     MAXDIMCOUL = INT(REAL(MAXNEBCOUL)*1.5)
  ENDIF

  RETURN

END SUBROUTINE NEBLISTS
