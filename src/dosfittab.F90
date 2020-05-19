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

SUBROUTINE DOSFITTAB

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE DIAGARRAY
  USE KSPACEARRAY
  USE SPINARRAY
  USE MDARRAY
  USE UNIVARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, K, JUNK, ACC, II, COUNT, N
  INTEGER :: NSTEP, PICK, CPLOC(1), LOOPTARGET, AB
  INTEGER :: MYELE, MYONSITE
  REAL(LATTEPREC) :: EMIN, EMAX, MERIT, OLDMERIT, RN, ESTEP, EF
  REAL(LATTEPREC) :: P, QN, SIG, UN, HSCLREMEMBER
  REAL(LATTEPREC), ALLOCATABLE :: DFTDOS(:), DFTINTDOS(:)
  REAL(LATTEPREC), ALLOCATABLE :: TBDOS(:), ONSITEBEST(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: HSCL(:), HSCLBEST(:), TABHORIG(:,:), U(:)
  REAL(LATTEPREC), PARAMETER ::  ETA = 0.05, DOSWEIGHT = 1.0D0
  REAL(LATTEPREC) :: ENERGY, INTDOS, NUME, JUNKNUM, MINERR
  REAL(LATTEPREC) :: ONSITE_REMEMBER, BIOROS
  COMPLEX(LATTEPREC) :: CMPARG
  REAL(LATTEPREC), EXTERNAL :: GAUSSRN

  MINERR = 1.0D6

  ! Read in the DOS from a VASP calculation

  OPEN(UNIT=20, STATUS="OLD", FILE="DOSCAR2fit.dat")

  READ(20,*) EMAX, EMIN, NSTEP, EF, JUNKNUM

  EMAX = EMAX - EF
  EMIN = EMIN - EF

  LOOPTARGET  = INT(BNDFIL*REAL(HDIM*NKTOT))
  
!  PRINT*, MAXVAL(LENTABINT)

  ALLOCATE(DFTDOS(NSTEP), DFTINTDOS(NSTEP), TBDOS(NSTEP))
  ALLOCATE(HSCL(NOINT), HSCLBEST(NOINT), TABHORIG(NOINT, MAXVAL(LENTABINT)))
  ALLOCATE(U(MAXVAL(LENTABINT)))
  ALLOCATE(ONSITEBEST(4,NOELEM))

  HSCL = ONE
  HSCLBEST = HSCL
  TABHORIG = TABH
  DO I = 1, NOELEM
     ONSITEBEST(1,I) = HES(I)
     ONSITEBEST(2,I) = HEP(I)
     ONSITEBEST(3,I) = HED(I)
     ONSITEBEST(4,I) = HEF(I)
  ENDDO


  DO I = 1, NSTEP

     READ(20,*) JUNKNUM, DFTDOS(I), DFTINTDOS(I)

  ENDDO

  CLOSE(20)

  ! Initialize stuff

  
  IF (CONTROL .EQ. 1) THEN
     CALL ALLOCATEDIAG
  ELSEIF (CONTROL .EQ. 2 .OR. CONTROL .EQ. 4 .OR. CONTROL .EQ. 5) THEN
     CALL ALLOCATEPURE
  ELSEIF (CONTROL .EQ. 3) THEN
     CALL FERMIALLOCATE
  ENDIF

  IF (ELECTRO .EQ. 1) THEN
     CALL ALLOCATECOULOMB
     CALL INITCOULOMB
  ENDIF

  !
  ! Allocate stuff for building the neighbor lists
  !

  CALL ALLOCATENEBARRAYS

  CALL NEBLISTS(0)  

  IF (BASISTYPE .EQ. "NONORTHO") CALL ALLOCATENONO

  ! Compute DOS from LATTE and compare

  CALL INITRNG

  ACC = 0

  NUME = ZERO
  DO I = 1, NATS
     NUME = NUME + ATOCC(ELEMPOINTER(I))
  ENDDO

  DO II = 1, NFITSTEP


     
     ! Pick one to change
     

     IF (II .GT. 1) THEN
        

        CALL RANDOM_NUMBER(RN)

        BIOROS = RN

        IF (BIOROS .LT. 0.5) THEN

           CALL RANDOM_NUMBER(RN)
           
           PICK = INT(RN*REAL(INT2FIT)) + 1
           
           CALL RANDOM_NUMBER(RN)
           
           CALL RANDOM_NUMBER(RN)
        
           HSCLREMEMBER = HSCL(PICK)
           
           HSCL(PICK) = HSCL(PICK)*(ONE + MCSIGMA*(TWO*RN - ONE))

           DO I = 1, NOINT
              DO J = 1, LENTABINT(I)
                 TABH(J,I) = TABHORIG(J,I)*HSCL(I)
              ENDDO
           ENDDO
           
           ! Re-do the splines
           
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

        ELSE ! Now try changing the onsites

           CALL RANDOM_NUMBER(RN)

           MYELE = INT(REAL(NOELEM)*RN) + 1

           CALL RANDOM_NUMBER(RN)

           MYONSITE = INT(4.0*RN)  + 1

           CALL RANDOM_NUMBER(RN)
           IF (MYONSITE .EQ. 1) THEN
              ONSITE_REMEMBER = HES(MYELE)
              HES(MYELE) = HES(MYELE)*(ONE + MCSIGMA*(TWO*RN - ONE))
           ELSEIF (MYONSITE .EQ. 2) THEN
              ONSITE_REMEMBER = HEP(MYELE)
              HEP(MYELE) = HEP(MYELE)*(ONE + MCSIGMA*(TWO*RN - ONE))
           ELSEIF (MYONSITE .EQ. 3) THEN
              ONSITE_REMEMBER = HED(MYELE)
              HED(MYELE) = HED(MYELE)*(ONE + MCSIGMA*(TWO*RN - ONE))
           ELSEIF (MYONSITE .EQ. 4) THEN
              ONSITE_REMEMBER = HEF(MYELE)
              HEF(MYELE) = HEF(MYELE)*(ONE + MCSIGMA*(TWO*RN - ONE))
           ENDIF
           
        ENDIF
        
     ENDIF
        
     ! Re-build cut-off tails
          
     ! Build H

     CALL KBLDNEWH

     IF (BASISTYPE .EQ. "NONORTHO") CALL KORTHOMYH

     IF (QFIT .EQ. 0) THEN

        CALL KDIAGMYH

     ! Find the chemical potential - we need the looptarget'th lowest 
     ! eigenvalue

        COUNT = 0
        DO I = 1, NKTOT
           DO J = 1, HDIM
              COUNT = COUNT + 1
              CPLIST(COUNT) = KEVALS(J,I)
           ENDDO
        ENDDO
        
        DO I = 1, LOOPTARGET
           CHEMPOT = MINVAL(CPLIST)
           CPLOC = MINLOC(CPLIST)
           CPLIST(CPLOC) = 1.0D12 ! We do this so we don't get this eigenvalue again on the next loop
        ENDDO
        
     ELSE
        
        
        IF (ELECTRO .EQ. 0) THEN
           CALL QNEUTRAL(0, 1)  ! Local charge neutrality
        ELSE
           CALL QCONSISTENCY(0,1) ! Self-consistent charges
        ENDIF
        
        
     ENDIF
     

     KEVALS = KEVALS - CHEMPOT
     
     ! Compute DOS

     TBDOS = ZERO
     
     ESTEP = (EMAX - EMIN)/REAL(NSTEP)
     
     DO I = 1, NSTEP
        
        ENERGY = EMIN + REAL(I-1)*ESTEP
        
        DO K = 1, NKTOT
           DO J = 1, HDIM
              
              CMPARG = ONE/(ENERGY - KEVALS(J,K) + CMPLX(ZERO,ETA))
              
              TBDOS(I) = TBDOS(I) - (ONE/PI)*AIMAG(CMPARG)
              
           ENDDO
        ENDDO
        
     ENDDO
     
     TBDOS = TBDOS/REAL(NKTOT)
     

     ! Normalize the DOS -> integral to Ef = Ne
     
     
     INTDOS = ZERO
     COUNT = 0
     DO I = 1, NSTEP
        
        ENERGY = EMIN + REAL(I-1)*ESTEP
        
        IF (ENERGY .LE. ZERO) THEN
           COUNT = COUNT + 1
           IF (MOD(I,2) .EQ. 0) THEN
              INTDOS = INTDOS + FOUR*TBDOS(I)
           ELSE
              INTDOS = INTDOS + TWO*TBDOS(I)
           ENDIF
        ENDIF
        
     ENDDO
     
     INTDOS = INTDOS - TBDOS(1) - TBDOS(COUNT)
     
     INTDOS = INTDOS*ESTEP/THREE
     
     TBDOS = TBDOS*NUME/INTDOS
     
     DO I = 1, NSTEP
        
        TBDOS(I) = ABS(TBDOS(I) - DFTDOS(I))
        
     ENDDO
     
     INTDOS = ZERO
     COUNT = 0
     DO I = 1, NSTEP
        
        ENERGY = EMIN + REAL(I-1)*ESTEP
        
!        IF (ENERGY .LE. ZERO) THEN
!           COUNT = COUNT + 1
           IF (MOD(I,2) .EQ. 0) THEN
              INTDOS = INTDOS + FOUR*TBDOS(I)
           ELSE
              INTDOS = INTDOS + TWO*TBDOS(I)
           ENDIF
!        ENDIF

     ENDDO
     
     INTDOS = INTDOS - TBDOS(1) - TBDOS(NSTEP)
     
     INTDOS = INTDOS*ESTEP/THREE
     
     MERIT = INTDOS
     
     IF (II .EQ. 1) THEN
        OLDMERIT = MERIT
        MINERR = MERIT
     ENDIF


     IF (MERIT .LT. OLDMERIT) THEN
        
        ACC = ACC + 1
        OLDMERIT = MERIT

        IF (MERIT .LT. MINERR) THEN
           
           MINERR = MERIT
           HSCLBEST = HSCL

           DO I = 1, NOELEM
              ONSITEBEST(1,I) = HES(I)
              ONSITEBEST(2,I) = HEP(I)
              ONSITEBEST(3,I) = HED(I)
              ONSITEBEST(4,I) = HEF(I)
           ENDDO
           
        ENDIF

     ELSE
        
        CALL RANDOM_NUMBER(RN)

        IF ( EXP( -(MERIT - OLDMERIT)*MCBETA) .GT. RN ) THEN
           
           ACC = ACC + 1
           OLDMERIT = MERIT
           
        ELSE
           
           IF (BIOROS .LT. 0.5) THEN
              HSCL(PICK) = HSCLREMEMBER
           ELSE

              IF (MYONSITE .EQ. 1) THEN
                 HES(MYELE) = ONSITE_REMEMBER
              ELSEIF (MYONSITE .EQ. 2) THEN
                 HEP(MYELE) = ONSITE_REMEMBER
              ELSEIF (MYONSITE .EQ. 2) THEN
                 HED(MYELE) = ONSITE_REMEMBER
              ELSEIF (MYONSITE .EQ. 2) THEN
                 HEF(MYELE) = ONSITE_REMEMBER
              ENDIF
              
           ENDIF

        ENDIF
        
     ENDIF
     
     WRITE(*,11) II, ACC, MERIT, OLDMERIT
!     WRITE(*,11) II, ACC, MERIT, OLDMERIT, &
!          (BOND(1,J), J = 1, INT2FIT)
11   FORMAT(2I9, 2F12.6) 
  
  ENDDO
  
  HSCL = HSCLBEST
!  PRINT*, HSCLBEST

  DO I = 1, NOINT
     DO J = 1, LENTABINT(I)
        TABH(J,I) = TABHORIG(J,I)*HSCL(I)
     ENDDO
  ENDDO

  DO I = 1, NOELEM
     HES(I) = ONSITEBEST(1,I)
     HEP(I) = ONSITEBEST(2,I)
     HED(I) = ONSITEBEST(3,I)
     HEF(I) = ONSITEBEST(4,I)
  ENDDO


 
!  TABH = TABHORIG
  ! Re-do the splines                                                                 
  
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
  
  
  CALL KBLDNEWH
  
  IF (BASISTYPE .EQ. "NONORTHO") CALL KORTHOMYH
  
  IF (QFIT .EQ. 0) THEN
     
     CALL KDIAGMYH
     
     ! Find the chemical potential - we need the looptarget'th lowest 
     ! eigenvalue
     
     COUNT = 0
     DO I = 1, NKTOT
        DO J = 1, HDIM
           COUNT = COUNT + 1
           CPLIST(COUNT) = KEVALS(J,I)
        ENDDO
     ENDDO
     
     DO I = 1, LOOPTARGET
        CHEMPOT = MINVAL(CPLIST)
        CPLOC = MINLOC(CPLIST)
        CPLIST(CPLOC) = 1.0D12 ! We do this so we don't get this eigenvalue again on the next loop
     ENDDO
     
  ELSE
     
     IF (ELECTRO .EQ. 0) THEN
        CALL QNEUTRAL(0, 1)  ! Local charge neutrality
     ELSE
        CALL QCONSISTENCY(0,1) ! Self-consistent charges
     ENDIF
     
     
  ENDIF
  
  KEVALS = KEVALS - CHEMPOT
  
  
  ! Build H                                                         
  
  
  ! Compute and normalize the DOS one more time

  TBDOS = ZERO
  
  ESTEP = (EMAX - EMIN)/REAL(NSTEP)
  
  DO I = 1, NSTEP
     
     ENERGY = EMIN + REAL(I-1)*ESTEP
     
     DO K = 1, NKTOT
        DO J = 1, HDIM
              
           CMPARG = ONE/(ENERGY - KEVALS(J,K) + CMPLX(ZERO,ETA))
           
           TBDOS(I) = TBDOS(I) - (ONE/PI)*AIMAG(CMPARG)
           
        ENDDO
     ENDDO
     
  ENDDO
  
  TBDOS = TBDOS/REAL(NKTOT)
  
    
  ! Normalize the DOS -> integral to Ef = Ne
  
  
  INTDOS = ZERO
  COUNT = 0
  DO I = 1, NSTEP
     
     ENERGY = EMIN + REAL(I-1)*ESTEP
     
     IF (ENERGY .LE. ZERO) THEN
        COUNT = COUNT + 1
        IF (MOD(I,2) .EQ. 0) THEN
           INTDOS = INTDOS + FOUR*TBDOS(I)
        ELSE
           INTDOS = INTDOS + TWO*TBDOS(I)
        ENDIF
     ENDIF
     
  ENDDO
  
  INTDOS = INTDOS - TBDOS(1) - TBDOS(COUNT)
  
  INTDOS = INTDOS*ESTEP/THREE
  
  TBDOS = TBDOS*NUME/INTDOS
    
  OPEN(UNIT=20, STATUS="UNKNOWN", FILE="checkdosfit.dat")
  
  DO I = 1, NSTEP
     
     ENERGY = EMIN + REAL(I-1)*ESTEP
     
     WRITE(20,10) ENERGY, TBDOS(I), DFTDOS(I)
     
  ENDDO
!  WRITE(20,*) HSCLBEST
  
12 FORMAT(A2,1X,50F12.6)

  CLOSE(20)
  
10 FORMAT(F12.3, 2G18.9)

  OPEN(UNIT=20, STATUS="UNKNOWN", FILE="bondints.table.opt")

  WRITE(20,'("Noints= ", I3)') NOINT
  DO I = 1, NOINT
     WRITE(20,*) ELE1(I), ELE2(I), BTYPE(I)
     WRITE(20,*) LENTABINT(I)
     DO J = 1, LENTABINT(I)
        WRITE(20,99) TABR(J,I), TABS(J,I), TABH(J,I)
     ENDDO
  ENDDO

99 FORMAT(3F16.8)
     
  CLOSE(20)

  OPEN(UNIT=20, STATUS="UNKNOWN", FILE="electrons.dat.opt")

  WRITE(20,*) "X ", "X ", "X ", "X ", "X ", "X ", "X ", "X ", "X ", "X ", "X ", "X ", "X "

  DO I = 1, NOELEM
     WRITE(20,30) ELE(I), BASIS(I), ATOCC(I), HES(I), HEP(I), HED(I), HEF(I), &
       MASS(I), HUBBARDU(I), WSS(I), WPP(I), WDD(I), WFF(I)
  ENDDO

30 FORMAT(A2,1X,A4,1X,11F12.6)
        
  CLOSE (20)
  
  DEALLOCATE(TBDOS, DFTDOS, DFTINTDOS)
  DEALLOCATE(HSCL, HSCLBEST, TABHORIG, U)
  DEALLOCATE(ONSITEBEST)

  RETURN

END SUBROUTINE DOSFITTAB
