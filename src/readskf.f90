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

SUBROUTINE READSKF
  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE UNIVARRAY
  USE MDARRAY
  USE NONOARRAY
  USE SPINARRAY
  USE SKFPARSER_
  
  INTEGER :: I, J, K, MAXENTRY, N, NUMENTRY
  REAL(LATTEPREC), ALLOCATABLE :: U(:)
  TYPE(SKFPARSER), ALLOCATABLE :: PARSERS(:)
  CHARACTER(100), ALLOCATABLE :: TOKENS(:)
  REAL(8), ALLOCATABLE :: RARRAY(:)
  INTEGER :: JJ
  CHARACTER(3), PARAMETER :: BTYPESF(10) = [ "dds", "ddp", "ddd", "pds", "pdp", "pps", "ppp", "sds", "sps", "sss" ]
  CHARACTER(3), PARAMETER :: BTYPEEF(20) = [ "ffs", "ffp", "ffd", "fff", "dfs", "dfp", "dfd", "dds", "ddp", "ddd", &
                                             "pfs", "pfp", "pds", "pdp", "pps", "ppp", "sfs", "sds", "sps", "sss" ]

	WRITE(*,*) "SKFFILES = ", trim(SKFFILES)
	call FString_split( SKFFILES, TOKENS, "," )
	
	ALLOCATE( PARSERS(SIZE(TOKENS)) )
	
	NOELEM = 0
	NOINT = 0
	MAXENTRY = 0
	DO I=1, SIZE(PARSERS)
		WRITE(*,*) ""
		WRITE(*,*) "Reading SKF file "//TRIM(PARAMPATH)//"/"//TRIM(TOKENS(I))
		WRITE(*,*) "------------------------------------------------"
		call PARSERS(I)%init( TRIM(PARAMPATH)//"/"//TRIM(TOKENS(I)) )
		call PARSERS(I)%show()
		
		IF ( PARSERS(I)%nGridPoints .GT. MAXENTRY ) MAXENTRY = PARSERS(I)%nGridPoints
		NOINT = NOINT + PARSERS(I)%nIntegrals
		
		IF( PARSERS(I)%homoNuclearCase ) NOELEM = NOELEM + 1
	END DO
	
	IF(.NOT.ALLOCATED(WSS))THEN
		ALLOCATE(WSS(NOELEM),WPP(NOELEM),WDD(NOELEM),WFF(NOELEM))
	ENDIF

	ALLOCATE(ELE(NOELEM), BASIS(NOELEM), ATOCC(NOELEM), HES(NOELEM), &
			HEP(NOELEM), HED(NOELEM), HEF(NOELEM), MASS(NOELEM), &
			HUBBARDU(NOELEM))
	
	DO I=1, SIZE(PARSERS)
		IF( PARSERS(I)%homoNuclearCase ) THEN
				ELE(I) = TRIM(PARSERS(I)%symbols(1))

				ATOCC(I) = SUM(PARSERS(I)%f)

				IF( .NOT. PARSERS(I)%extendedFormat ) THEN
					BASIS(I) = "spdf"
					
					HES(I) = PARSERS(I)%E(4)
					HEP(I) = PARSERS(I)%E(3)
					HED(I) = PARSERS(I)%E(2)
					HEF(I) = PARSERS(I)%E(1)
					
!             HUBBARDU(I) = PARSERS(I)%U(4)
!             HUBBARDU(I) = PARSERS(I)%U(3)
!             HUBBARDU(I) = PARSERS(I)%U(2)
					HUBBARDU(I) = PARSERS(I)%U(1) ! Hay que ver que capa seleccionar, probablemente con las occupaciones

					! Toca hacer otro fichero como los .yam de ADF para los Hubbard magneticos
					WSS(I) = PARSERS(I)%W(4)
					WPP(I) = PARSERS(I)%W(3)
					WDD(I) = PARSERS(I)%W(2)
					WFF(I) = PARSERS(I)%W(1)
				ELSE
					BASIS(I) = "spd"
					
					HES(I) = PARSERS(I)%E(3)
					HEP(I) = PARSERS(I)%E(2)
					HED(I) = PARSERS(I)%E(1)
					HEF(I) = 0.0
					
!             HUBBARDU(I) = PARSERS(I)%U(3)
!             HUBBARDU(I) = PARSERS(I)%U(2)
					HUBBARDU(I) = PARSERS(I)%U(1) ! Hay que ver que capa seleccionar, probablemente con las occupaciones
					
					! Toca hacer otro fichero como los .yam de ADF para los Hubbard magneticos
					WSS(I) = PARSERS(I)%W(3)
					WPP(I) = PARSERS(I)%W(2)
					WDD(I) = PARSERS(I)%W(1)
				END IF
				
				MASS(I) = PARSERS(I)%mass
		END IF
	END DO

	ALLOCATE(ELE1(NOINT), ELE2(NOINT), BTYPE(NOINT))
	
	ALLOCATE(TABR(MAXENTRY,NOINT), TABH(MAXENTRY, NOINT), TABS(MAXENTRY, NOINT), &
				LENTABINT(NOINT), HSPL(MAXENTRY, NOINT), SSPL(MAXENTRY, NOINT), &
				HCUT(NOINT), SCUT(NOINT))

	TABR = ZERO
	TABH = ZERO
	TABS = ZERO

	HCUT = ZERO
	SCUT = ZERO
	
	K=1
	DO I=1, SIZE(PARSERS)
		ELE1(K) = TRIM(PARSERS(I)%symbols(1))
		ELE2(K) = TRIM(PARSERS(I)%symbols(1))
		
		DO J=1, PARSERS(I)%nIntegrals
			
			IF( .NOT. PARSERS(I)%extendedFormat ) THEN
				BTYPE(K) = TRIM(BTYPESF(J))
			ELSE
				BTYPE(K) = TRIM(BTYPEEF(J))
			END IF
			
			LENTABINT(K) = PARSERS(I)%nGridPoints
			
			CALL PARSERS(I)%getR( RARRAY )
			TABR(:,K) = RARRAY( 1:SIZE(TABR(:,K)) )/1.88972612456506_8  ! bohrs to angs
			
			CALL PARSERS(I)%getH( J, RARRAY )
			TABH(:,K) = RARRAY( 1:SIZE(TABH(:,K)) )/0.0367493088244753_8  ! Hartree to eV
			
			CALL PARSERS(I)%getS( J, RARRAY )
			TABS(:,K) = RARRAY( 1:SIZE(TABS(:,K)) )
			
			DO JJ = 1, LENTABINT(K)
				IF (TABR(JJ,K) .GT. HCUT(K)) HCUT(K) = TABR(JJ,K)
			ENDDO
			SCUT(K) = HCUT(K)
			
			K=K+1
		END DO
	END DO
	
	IF( ALLOCATED(RARRAY) ) DEALLOCATE(RARRAY)
	IF( ALLOCATED(TOKENS) ) DEALLOCATE(TOKENS)
	IF( ALLOCATED(PARSERS) ) DEALLOCATE(PARSERS)
	
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

END SUBROUTINE READSKF
