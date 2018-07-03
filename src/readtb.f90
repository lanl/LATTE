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

SUBROUTINE READTB

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE UNIVARRAY
  USE MDARRAY
  USE NONOARRAY
  USE SPINARRAY
  USE KSPACEARRAY
  USE LATTEPARSER_LATTE_MOD

  IMPLICIT NONE

  INTEGER :: I, J, K, MAXENTRY, N, NUMENTRY
  CHARACTER(LEN=20) :: HD
  REAL(LATTEPREC) :: TAILPARAMS(6)
  REAL(LATTEPREC) :: JUNK, P, QN, SIG, UN
  REAL(LATTEPREC), ALLOCATABLE :: U(:)

  IF (EXISTERROR) RETURN

  OPEN(UNIT=22,STATUS="OLD", FILE=TRIM(PARAMPATH)//"/electrons.dat")

  READ(22,*) HD, NOELEM

  IF(.NOT.ALLOCATED(WSS))THEN
     ALLOCATE(WSS(NOELEM),WPP(NOELEM),WDD(NOELEM),WFF(NOELEM))
  ENDIF

  ALLOCATE(ELE(NOELEM), BASIS(NOELEM), ATOCC(NOELEM), HES(NOELEM), &
       HEP(NOELEM), HED(NOELEM), HEF(NOELEM), MASS(NOELEM), &
       HUBBARDU(NOELEM))

  READ(22,*) HD, HD, HD, HD, HD, HD, HD, HD, HD, HD, HD, HD, HD

  DO I = 1, NOELEM
     READ(22,*) ELE(I), BASIS(I), ATOCC(I), HES(I), HEP(I), HED(I), HEF(I), &
          MASS(I), HUBBARDU(I), WSS(I), WPP(I), WDD(I), WFF(I)
  ENDDO

  CLOSE(22)

  IF (SCLTYPE .EQ. "EXP") THEN

     IF (BASISTYPE .EQ. "ORTHO") THEN
        OPEN(UNIT=11,STATUS="OLD", FILE=TRIM(PARAMPATH)//"/bondints.ortho")
     ELSE
        OPEN(UNIT=11,STATUS="OLD", FILE=TRIM(PARAMPATH)//"/bondints.nonortho")
     ENDIF

  ELSEIF (SCLTYPE .EQ. "TABLE") THEN

     OPEN(UNIT=11,STATUS="OLD", FILE=TRIM(PARAMPATH)//"/bondints.table")

  ELSE 
     print*, "Choose SCLTYPE either EXP or TABLE"
     STOP
  ENDIF


  READ(11,*) HD, NOINT

  ALLOCATE(ELE1(NOINT), ELE2(NOINT), BTYPE(NOINT))


  IF (SCLTYPE .EQ. "EXP") THEN

     IF (BASISTYPE .EQ. "ORTHO") THEN
        
        ALLOCATE(BOND(14, NOINT), HCUT(NOINT))

        READ(11,*) HD, HD, HD, HD, HD, HD, HD, HD, HD, HD, HD
        
        DO I = 1, NOINT
           
           READ(11,*) ELE1(I), ELE2(I), BTYPE(I), (BOND(J,I), J = 1, 8)
           HCUT(I) = BOND(8,I)

        ENDDO
        
     ELSEIF (BASISTYPE .EQ. "NONORTHO") THEN

        ALLOCATE( BOND(14,NOINT), HCUT(NOINT), OVERL(14,NOINT), SCUT(NOINT) )
        
        READ(11,*) HD, HD, HD, HD, HD, HD, HD, HD, HD, HD, HD, &
             HD, HD, HD, HD, HD, HD, HD, HD
        
        DO I = 1, NOINT

           READ(11,*) ELE1(I), ELE2(I), BTYPE(I), (BOND(J,I), J = 1, 8), &
                (OVERL(K,I), K = 1, 8)
           HCUT(I) = BOND(8,I)
           SCUT(I) = OVERL(8,I)
           
        ENDDO
        
     ENDIF
     
  ELSEIF (SCLTYPE .EQ. "TABLE") THEN

     ! Read in the tables derived from Plato                                     

     MAXENTRY = 0
     DO I = 1, NOINT
        READ(11,*) HD, HD, HD
        READ(11,*) NUMENTRY

        IF (NUMENTRY .GT. MAXENTRY) MAXENTRY = NUMENTRY

        DO J = 1, NUMENTRY
           READ(11,*) JUNK, JUNK, JUNK
        ENDDO
     ENDDO

     REWIND(11)

     ALLOCATE(TABR(MAXENTRY,NOINT), TABH(MAXENTRY, NOINT), TABS(MAXENTRY, NOINT), &
          LENTABINT(NOINT), HSPL(MAXENTRY, NOINT), SSPL(MAXENTRY, NOINT), &
          HCUT(NOINT), SCUT(NOINT))


     TABR = ZERO
     TABH = ZERO
     TABS = ZERO

     HCUT = ZERO
     SCUT = ZERO

     
     READ(11,*) HD, NOINT
     DO I = 1, NOINT
        READ(11,*) ELE1(I), ELE2(I), BTYPE(I)
        READ(11,*) LENTABINT(I)
        DO J = 1, LENTABINT(I)
           READ(11,*) TABR(J,I), TABS(J,I), TABH(J,I)
        ENDDO

        DO J = 1, LENTABINT(I)
           IF (TABR(J,I) .GT. HCUT(I)) HCUT(I) = TABR(J,I)
        ENDDO
        SCUT(I) = HCUT(I)

     ENDDO

     ALLOCATE(U(MAXENTRY))

     ! H first                                                                   

     DO I = 1, NOINT

        N = LENTABINT(I)

        HSPL(1,I) = ZERO
        U(1) = ZERO

        DO J = 2, N-1
           SIG = (TABR(J,I) - TABR(J-1,I))/(TABR(J+1,I) - TABR(J-1,I))
           P = SIG*HSPL(J-1,I) + TWO
           HSPL(J,I) = (SIG - ONE)/P
           U(J) = (SIX*((TABH(J+1,I) - TABH(J,I)) / &
                (TABR(J+1,I) - TABR(J,I)) - (TABH(J,I) - TABH(J-1,I)) &
                /(TABR(J,I) - TABR(J-1,I)))/(TABR(J+1,I)-TABR(J-1,I)) &
                - SIG*U(J-1))/P
        ENDDO

        QN = ZERO
        UN = ZERO

        HSPL(N,I) = (UN - QN*U(N-1))/(QN*HSPL(N-1, I) + ONE)
        
        DO K = N-1, 1, -1
           HSPL(K,I) = HSPL(K,I)*HSPL(K+1,I) + U(K)
        ENDDO

     ENDDO

          ! Now the overlap                                                           

     DO I = 1, NOINT

        N = LENTABINT(I)

        SSPL(1,I) = ZERO
        U(1) = ZERO

        DO J = 2, N-1
           SIG = (TABR(J,I) - TABR(J-1,I))/(TABR(J+1,I) - TABR(J-1,I))
           P = SIG*SSPL(J-1,I) + TWO
           SSPL(J,I) = (SIG - ONE)/P
           U(J) = (SIX*((TABS(J+1,I) - TABS(J,I)) / &
                (TABR(J+1,I) - TABR(J,I)) - (TABS(J,I) - TABS(J-1,I)) &
                /(TABR(J,I) - TABR(J-1,I)))/(TABR(J+1,I)-TABR(J-1,I)) &
                - SIG*U(J-1))/P
        ENDDO

        QN = ZERO
        UN = ZERO

        SSPL(N,I) = (UN - QN*U(N-1))/(QN*SSPL(N-1, I) + ONE)

        DO K = N-1, 1, -1
           SSPL(K,I) = SSPL(K,I)*SSPL(K+1,I) + U(K)
        ENDDO

     ENDDO

     DEALLOCATE(U)

  ENDIF

  CLOSE(11)

  ! If we're doing k-space integration, let's read in the k point mesh
  IF (KON .EQ. 1) THEN

     IF (LATTEINEXISTS) THEN
        CALL PARSE_KMESH("latte.in")
     ELSE
        OPEN(UNIT=11, STATUS="OLD", FILE=TRIM(PARAMPATH)//"/kmesh.in")
        READ(11,*) NKX, NKY, NKZ
        READ(11,*) KSHIFT(1), KSHIFT(2), KSHIFT(3)
        CLOSE (11)
     ENDIF

     NKTOT = NKX*NKY*NKZ

  ENDIF

  IF (SCLTYPE .EQ. "EXP") THEN

     DO I = 1, NOINT
        
        CALL UNIVTAILCOEF(BOND(:,I))
        
        IF (BASISTYPE .EQ. "NONORTHO") CALL UNIVTAILCOEF(OVERL(:,I))
        
     ENDDO
     
  ENDIF

  RETURN

END SUBROUTINE READTB
