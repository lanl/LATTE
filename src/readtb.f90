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

  INTEGER :: I, J, K
  CHARACTER(LEN=20) :: HD
  REAL(LATTEPREC) :: TAILPARAMS(6)
  LOGICAL :: LATTEINEXISTS

  OPEN(UNIT=22,STATUS="OLD", FILE=trim(PARAMPATH)//"/electrons.dat")

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

  IF (BASISTYPE .EQ. "ORTHO") THEN
     OPEN(UNIT=11,STATUS="OLD", FILE=trim(PARAMPATH)//"/bondints.ortho")
  ELSE
     OPEN(UNIT=11,STATUS="OLD", FILE=trim(PARAMPATH)//"/bondints.nonortho")
  ENDIF

  READ(11,*) HD, NOINT

  ALLOCATE(ELE1(NOINT), ELE2(NOINT), BTYPE(NOINT),&
       BOND(14, NOINT))


  IF (BASISTYPE .EQ. "ORTHO") THEN
     
     READ(11,*) HD, HD, HD, HD, HD, HD, HD, HD, HD, HD, HD
      
     DO I = 1, NOINT
        
        READ(11,*) ELE1(I), ELE2(I), BTYPE(I), (BOND(J,I), J = 1, 8)
        
     ENDDO

  ELSE
     
     ALLOCATE( OVERL(14,NOINT) )
 
     READ(11,*) HD, HD, HD, HD, HD, HD, HD, HD, HD, HD, HD, &
          HD, HD, HD, HD, HD, HD, HD, HD
    
     DO I = 1, NOINT
        
        READ(11,*) ELE1(I), ELE2(I), BTYPE(I), (BOND(J,I), J = 1, 8), &
             (OVERL(K,I), K = 1, 8)
        
     ENDDO
     
  ENDIF

  CLOSE(11)

  ! If we're doing k-space integration, let's read in the k point mesh
  IF (KON .EQ. 1) THEN

    INQUIRE( FILE="latte.in", exist=LATTEINEXISTS )
    IF (LATTEINEXISTS) THEN
      CALL PARSE_KMESH("latte.in")
    ELSE
      OPEN(UNIT=11, STATUS="OLD", FILE=trim(PARAMPATH)//"/kmesh.in")
      READ(11,*) NKX, NKY, NKZ
      READ(11,*) KSHIFT(1), KSHIFT(2), KSHIFT(3)
      CLOSE (11)
    ENDIF

    NKTOT = NKX*NKY*NKZ

  ENDIF
     
  DO I = 1, NOINT

     CALL UNIVTAILCOEF(BOND(:,I))

     IF (BASISTYPE .EQ. "NONORTHO") CALL UNIVTAILCOEF(OVERL(:,I))

  ENDDO
  
  RETURN
  
END SUBROUTINE READTB
