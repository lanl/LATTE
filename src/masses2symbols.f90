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

!> Soubroutine to get the Symbols out of the masses of the elements. 
!! \param TYPES atom type index.
!! \param NTYPES Number of types.
!! \param MASSES_IN Masses for every type.
!! \param NATSIN Number of total atoms. 
!! \param SYMBOLS Symbols for every atom. 
!! 
SUBROUTINE MASSES2SYMBOLS(TYPES,NTYPES,MASSES_IN,NATSIN,SYMBOLS)

  USE SETUPARRAY
  USE MYPRECISION 

  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: NTYPES
  INTEGER, INTENT(IN) :: TYPES(NTYPES)
  INTEGER :: NATSIN
    
  REAL(LATTEPREC), INTENT(IN) :: MASSES_IN(NTYPES)
  CHARACTER(LEN=2), INTENT(INOUT) :: SYMBOLS(NATSIN)
  
  INTEGER, PARAMETER :: NZ = 103 
  INTEGER, PARAMETER :: DP = LATTEPREC   

  !> Element symbol
  !!
  CHARACTER(2), PARAMETER :: ELEMENT_SYMBOL(NZ) = [CHARACTER(2) :: &
    "H" ,          "He" ,         "Li" ,         "Be" ,         &
    "B" ,          "C" ,          "N" ,          "O" ,          &
    "F" ,          "Ne" ,         "Na" ,         "Mg" ,         &
    "Al" ,         "Si" ,         "P" ,          "S" ,          &
    "Cl" ,         "Ar" ,         "K" ,          "Ca" ,         &
    "Sc" ,         "Ti" ,         "V" ,          "Cr" ,         &
    "Mn" ,         "Fe" ,         "Co" ,         "Ni" ,         &
    "Cu" ,         "Zn" ,         "Ga" ,         "Ge" ,         &
    "As" ,         "Se" ,         "Br" ,         "Kr" ,         &
    "Rb" ,         "Sr" ,         "Y" ,          "Zr" ,         &
    "Nb" ,         "Mo" ,         "Tc" ,         "Ru" ,         &
    "Rh" ,         "Pd" ,         "Ag" ,         "Cd" ,         &
    "In" ,         "Sn" ,         "Sb" ,         "Te" ,         &
    "I" ,          "Xe" ,         "Cs" ,         "Ba" ,         &
    "La" ,         "Ce" ,         "Pr" ,         "Nd" ,         &
    "Pm" ,         "Sm" ,         "Eu" ,         "Gd" ,         &
    "Tb" ,         "Dy" ,         "Ho" ,         "Er" ,         &
    "Tm" ,         "Yb" ,         "Lu" ,         "Hf" ,         &
    "Ta" ,         "W" ,          "Re" ,         "Os" ,         &
    "Ir" ,         "Pt" ,         "Au" ,         "Hg" ,         &
    "Tl" ,         "Pb" ,         "Bi" ,         "Po" ,         &
    "At" ,         "Rn" ,         "Fr" ,         "Ra" ,         &
    "Ac" ,         "Th" ,         "Pa" ,         "U" ,          &
    "Np" ,         "Pu" ,         "Am" ,         "Cm" ,         &
    "Bk" ,         "Cf" ,         "Es" ,         "Fm" ,         &
    "Md" ,         "No" ,         "Lr"                          &
    ]

  !> Element mass in atomic mass units (1.66 x 10-27 kg)
  !!
  REAL(DP), PARAMETER :: ELEMENT_MASS(NZ) = (/ &
    1.007825032 ,  4.002603254 ,  7.01600455 ,   9.0121822 ,    &
    11.0093054 ,   12.0 ,         14.003074005 , 15.99491462 ,  &
    18.99840322 ,  19.992440175 , 22.989769281 , 23.9850417 ,   &
    26.98153863 ,  27.976926532 , 30.97376163 ,  31.972071 ,    &
    34.96885268 ,  39.962383123 , 38.96370668 ,  39.96259098 ,  &
    44.9559119 ,   47.9479463 ,   50.9439595 ,   51.9405075 ,   &
    54.9380451 ,   55.9349375 ,   58.933195 ,    57.9353429 ,   &
    62.9295975 ,   63.929142 ,    68.925573 ,    73.921177 ,    &
    74.921596 ,    79.916521 ,    78.918337 ,    83.911507 ,    &
    84.911789 ,    87.905612 ,    88.905848 ,    89.904704 ,    &
    92.906378 ,    97.905408 ,    97.907216 ,    101.904349 ,   &
    102.905504 ,   105.903486 ,   106.905097 ,   113.903358 ,   &
    114.903878 ,   119.902194 ,   120.903815 ,   129.906224 ,   &
    126.904473 ,   131.904153 ,   132.905451 ,   137.905247 ,   &
    138.906353 ,   139.905438 ,   140.907652 ,   141.907723 ,   &
    144.912749 ,   151.919732 ,   152.92123 ,    157.924103 ,   &
    158.925346 ,   163.929174 ,   164.930322 ,   165.930293 ,   &
    168.934213 ,   173.938862 ,   174.940771 ,   179.94655 ,    &
    180.947995 ,   183.950931 ,   186.955753 ,   191.96148 ,    &
    192.962926 ,   194.964791 ,   196.966568 ,   201.970643 ,   &
    204.974427 ,   207.976652 ,   208.980398 ,   208.98243 ,    &
    209.987148 ,   222.017577 ,   223.019735 ,   226.025409 ,   &
    227.027752 ,   232.038055 ,   231.035884 ,   238.050788 ,   &
    237.048173 ,   244.064204 ,   243.061381 ,   247.070354 ,   &
    247.070307 ,   251.079587 ,   252.08298 ,    257.095105 ,   &
    258.098431 ,   259.10103 ,    262.10963                     &
    /)

    INTEGER :: I, J
    CHARACTER(LEN=2), ALLOCATABLE :: TYPE_SYMBOLS(:)   
 
    ALLOCATE(TYPE_SYMBOLS(NTYPES))
  
    DO I=1,NTYPES
      DO J=1,NZ
        IF(ABS(MASSES_IN(I) - ELEMENT_MASS(J)) < 0.5) THEN 
          TYPE_SYMBOLS(I) = ELEMENT_SYMBOL(J)
          EXIT
        ENDIF  
      ENDDO       
    ENDDO           

    DO I=1,NATSIN
      SYMBOLS(I) = TYPE_SYMBOLS(TYPES(I))
    ENDDO  

    DEALLOCATE(TYPE_SYMBOLS)    
        
END SUBROUTINE MASSES2SYMBOLS


