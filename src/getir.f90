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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                !
!  Some functions were taken from the qmd-progress project                       !
!  https://github.com/lanl/qmd-progress                                          !
!--------------------------------------------------------------------------------!
!                                                                                !
! Copyright (c) 2016, Los Alamos National Security, LLC                          !
! All rights reserved.                                                           !
!                                                                                !
! Copyright 2016. Los Alamos National Security, LLC. This software was           !
! produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos       !
! National Laboratory (LANL), which is operated by Los Alamos National           !
! Security, LLC for the U.S. Department of Energy. The U.S. Government has       !
! rights to use, reproduce, and distribute this software.  NEITHER THE           !
! GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY,           !
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.     !
! If software is modified to produce derivative works, such modified             !
! software should be clearly marked, so as not to confuse it with the version    !
! available from LANL.                                                           !
!                                                                                !
! Additionally, redistribution and use in source and binary forms, with or       !
! without modification, are permitted provided that the following conditions are !
! met:                                                                           !
! - Redistributions of source code must retain the above copyright               !
! notice, this list of conditions and the following disclaimer.                  !
! - Redistributions in binary form must reproduce the above                      !
! copyright notice, this list of conditions and the following disclaimer in      !
! the documentation and/or other materials provided with the distribution.       !
! - Neither the name of Los Alamos National Security, LLC, Los Alamos            !
! National Laboratory, LANL, the U.S. Government, nor the names of its           !
! contributors may be used to endorse or promote products derived from this      !
! software without specific prior written permission.                            !
!                                                                                !
! THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND             !
! CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT     !
! NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A    !
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL       !
! SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,              !
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,           !
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS          !
! OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON      !
! ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT        !
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF       !
! THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              !
!                                                                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                   !!
!!  Some functions were taken from the Molecule class of the SciFT project           !!
!!-----------------------------------------------------------------------------------!!
!!                                                                                   !!
!!  This file is part of SciFT project                                               !!
!!  https://github.com/nfaguirrec/scift                                              !!
!!  Copyright (c) 2010-2014 Nestor F. Aguirre (nfaguirrec@gmail.com)                 !!
!!                                                                                   !!
!!  Redistribution and use in source and binary forms, with or without               !!
!!  modification, are permitted provided that the following conditions are met:      !!
!!                                                                                   !!
!!  1. Redistributions of source code must retain the above copyright notice, this   !!
!!     list of conditions and the following disclaimer.                              !!
!!  2. Redistributions in binary form must reproduce the above copyright notice,     !!
!!     this list of conditions and the following disclaimer in the documentation     !!
!!     and/or other materials provided with the distribution.                        !!
!!  3. Neither the name of the copyright holders nor the names of its contributors   !!
!!     may be used to endorse or promote products derived from this software         !!
!!     without specific prior written permission.                                    !!
!!                                                                                   !!
!!  The copyright holders provide no reassurances that the source code provided      !!
!!  does not infringe any patent, copyright, or any other intellectual property      !!
!!  rights of third parties.  The copyright holders disclaim any liability to any    !!
!!  recipient for claims brought against recipient by any third party for            !!
!!  infringement of that parties intellectual property rights.                       !!
!!                                                                                   !!
!!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  !!
!!  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    !!
!!  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           !!
!!  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR  !!
!!  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES   !!
!!  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     !!
!!  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      !!
!!  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       !!
!!  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    !!
!!  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     !!
!!                                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!>
!! Returns the center of mass
!!
subroutine getCenterOfMass( centerOfMass )
  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MDARRAY

  implicit none

  real(8) :: centerOfMass(3)
  
  centerOfMass(1) = sum(MASS(ELEMPOINTER(:))*CR(1,:))
  centerOfMass(2) = sum(MASS(ELEMPOINTER(:))*CR(2,:))
  centerOfMass(3) = sum(MASS(ELEMPOINTER(:))*CR(3,:))
  centerOfMass = centerOfMass / sum(MASS(ELEMPOINTER(:)))
  
  centerOfMass = centerOfMass*1.88972612456506_8 ! angs to a.u.
  
end subroutine getCenterOfMass

!!>
!! Gets the inertia tensor and its eigen vectors and eigen values
!!
subroutine getInertiaTensor( inertiaTensor, nRot )
  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MDARRAY

  implicit none
  
  real(8), intent(out) :: inertiaTensor(3,3)
  integer, optional, intent(out) :: nRot
  
  real(8) :: eValsInertiaTensor(3)
  real(8) :: eVecsInertiaTensor(3,3)
  
  integer :: i, j, atom1, atom2
  real(8), allocatable :: geometryInCM(:,:)
  real(8) :: centerOfMass(3)
  real(8) :: massi
  
  integer :: ssign
  
  real(8), allocatable :: workSpace(:)
  integer :: ssize, info
  
  ALLOCATE( geometryInCM(3,NATS) )
  
  call getCenterOfMass( centerOfMass )
  
  do i=1,NATS
    geometryInCM(:,i) = CR(:,i)*1.88972612456506_8 - centerOfMass(:)   ! angs to a.u.
  end do
  
  inertiaTensor = 0.0_8
  do i=1,NATS
    massi = real(MASS(ELEMPOINTER(i)),8)*1822.88853_8 ! amu to a.u.
    inertiaTensor(1,1) = inertiaTensor(1,1) + massi * ( geometryInCM(2,i)**2 + geometryInCM(3,i)**2)
    inertiaTensor(1,2) = inertiaTensor(1,2) - massi * ( geometryInCM(1,i) * geometryInCM(2,i) )
    inertiaTensor(1,3) = inertiaTensor(1,3) - massi * ( geometryInCM(1,i) * geometryInCM(3,i) )
    inertiaTensor(2,2) = inertiaTensor(2,2) + massi * ( geometryInCM(1,i)**2 + geometryInCM(3,i)**2)
    inertiaTensor(2,3) = inertiaTensor(2,3) - massi * ( geometryInCM(2,i) * geometryInCM(3,i) )
    inertiaTensor(3,3) = inertiaTensor(3,3) + massi * ( geometryInCM(1,i)**2 + geometryInCM(2,i)**2)
  end do
  
  inertiaTensor(2,1) =inertiaTensor(1,2)
  inertiaTensor(3,1) =inertiaTensor(1,3)
  inertiaTensor(3,2) =inertiaTensor(2,3)
  
  eVecsInertiaTensor = inertiaTensor
  allocate( workSpace( 3*3*3-1 ) )
  
  ! Compute the eigen values and eigen vectors using the upper elements of the symmetric matrix
  call dsyev( 'V', 'L', 3, eVecsInertiaTensor, 3, eValsInertiaTensor, workSpace, 3*3*3-1, info )
  
  if( info /= 0 ) then
    write(*,*) "### ERROR ### Diagonalizing the inertia tensor"
    stop
  end if
  
  deallocate( workSpace )
  
  !! Checks the determinant's sign
  ssign = eVecsInertiaTensor(1,1)*( &
          eVecsInertiaTensor(2,2)*eVecsInertiaTensor(3,3) &
          -eVecsInertiaTensor(3,2)*eVecsInertiaTensor(2,3)) &
          -eVecsInertiaTensor(1,2)*( &
          eVecsInertiaTensor(2,1)*eVecsInertiaTensor(3,3) &
          -eVecsInertiaTensor(3,1)*eVecsInertiaTensor(2,3)) &
          +eVecsInertiaTensor(1,3)*( &
          eVecsInertiaTensor(2,1)*eVecsInertiaTensor(3,2) &
          -eVecsInertiaTensor(3,1)*eVecsInertiaTensor(2,2))
  
  !! Presers the handedness of the inertia tensor
  if ( ssign < 0.0 ) then
    eVecsInertiaTensor(1,2) = -eVecsInertiaTensor(1,2)
    eVecsInertiaTensor(2,2) = -eVecsInertiaTensor(2,2)
    eVecsInertiaTensor(3,2) = -eVecsInertiaTensor(3,2)
  endif
  
  !! Verifies if the inertia tensor is correct
  if ( 	abs( eVecsInertiaTensor(1,1) ) < 1d-10 .and. &
        abs( eVecsInertiaTensor(2,2) ) < 1d-10 .and. &
        abs( eVecsInertiaTensor(3,3) ) < 1d-10 ) then
    write(*,*) "### ERROR ### Invalid inertia tensor."
    stop
  end if
  
  if( present(nRot) ) then
    nRot = 3
    if( abs( eValsInertiaTensor(1) ) < 10.0 ) then
      nRot = 2 ! @todo Here is important make a second check to verify that molecule is indeed lineal
    end if
    write(*,*) ""
    write(*,*) "nRot = ", nRot
  end if
  
  write(*,"(A,3F20.5)") "Inertia moments: ", eValsInertiaTensor
  write(*,*) ""
  write(*,*) "Inertia tensor: "
  do i=1,3
    write(*,"(3F15.2)") inertiaTensor(:,i)
  end do
  write(*,*) ""
  
  inertiaTensor = eVecsInertiaTensor !<<< @todo OJO esto no es verdad
  
  write(*,*) ""
  write(*,*) "Centered geometry:"
  do atom1=1,NATS
    write(*,"(F8.2,4X,3F8.2)") MASS(ELEMPOINTER(atom1)), geometryInCM(:,atom1)
  end do
  write(*,*) ""
  
  deallocate( geometryInCM )

end subroutine getInertiaTensor

!> 
!! @brief Proyecta un numero dado de elementos sobre el resto de los elementos del espacio vectorial 
!!		y ortogonaliza el espacio resultante mediante un proceso Gram-Schmidt
!!
!! @param this Espacio vectorial
!<
subroutine projectLastElements( vectorialSpace, output )
  implicit none
  real(8) :: vectorialSpace(:,:)
  real(8), allocatable :: output(:,:)
  
  integer :: i
  integer :: last
  real(8) :: squareNorm
  real(8) :: projectionOverOrthogonalizedBasis
  
  last = size(vectorialSpace,dim=2)
  allocate( output( size(vectorialSpace,dim=1), last ) )
  output = vectorialSpace
  
  !!***********************************************************************************
  !! Realiza de ortogonalizacion sobre los last-1 vectores, previamente ortogonalizados.
  !!
  do i=1,last-1
    squareNorm = dot_product( output(:,i), output(:,i) )
    
    projectionOverOrthogonalizedBasis = dot_product( output(:,i), vectorialSpace(:,last) )
    
    if ( squareNorm>1.0D-12 ) then
      output( :, last ) = output( :, last ) - projectionOverOrthogonalizedBasis/sqrt(squareNorm)*output(:,i)
    end if
  end do
  
  squareNorm = dot_product( output(:,last), output(:,last) )	
  
  if ( squareNorm>1.0D-12 ) then
    output( :, last )=output( :, last )/sqrt(squareNorm)
  end if
  
  !!
  !!******************************************************************

end subroutine projectLastElements

!!
!! @brief Orthogonalize the components of a vectorial space
!!
subroutine orthogonalizeLinearVectorialSpace( matrix )
  implicit none
  real(8), allocatable :: matrix(:,:)
  
  interface
    subroutine projectLastElements( vectorialSpace, output )
      implicit none
      real(8) :: vectorialSpace(:,:)
      real(8), allocatable :: output(:,:)
    end subroutine projectLastElements
  end interface

  integer :: i
  integer :: last
  real(8) :: norm
  real(8), allocatable :: output(:,:)
  
  last = size(matrix,dim=2)
  norm=sqrt(dot_product(matrix(:,1),matrix(:,1)))
  
  matrix(:,1)=matrix(:,1)/norm
  
  !!
  !! Realiza de ortogonalizacion consecutiva de cada uno de los vectores
  !! presentes en la matriz
  !!
  do i=2,last
    call projectLastElements( matrix(:,1:i), output )
    matrix(:,1:i) = output
    if( allocated(output) ) deallocate(output)
  end do
  
  !! Reortonormaliza para asegurar la ortonormalizacion
  do i=2,last
    call projectLastElements( matrix(:,1:i), output )
    matrix(:,1:i) = output
    if( allocated(output) ) deallocate(output)
  end do
  
end subroutine orthogonalizeLinearVectorialSpace

!>
!! indexing array so that array(indexes(j)), j=1..n is in
!! ascending numerical order.
!! method is heapsort, see also subroutine hpsort.
!! taken from numerical recipies, p 233.
!!
subroutine rsort( array, indexes )
  real(8), intent(in) :: array(:)
  integer, intent(inout) :: indexes(:)
  
  integer :: i, j, l
  integer :: n
  integer :: id, ir
  real(8) :: value
  
!   if( .not. allocated(array) ) then
!     write(*,*) "Error in Math_sort, array not allocated"
!     stop
!   end if
!   
!   if( .not. allocated(indexes) ) then
!     write(*,*) "Error in Math_sort, indexes not allocated"
!     stop
!   end if
  
  if( size(array) /= size(indexes) ) then
    write(*,*) "Error in Math_sort, array and indexes have different size"
    stop
  end if
  
  n = size(array)
  
  do j=1,n
    indexes(j)=j
  end do
  
  if( n == 1 ) return
  
  l=n/2+1
  ir=n
  
  do while( .true. )
    if( l > 1 ) then
      l = l-1
      id = indexes(l)
      value = array(id)
    else
      id=indexes(ir)
      value=array(id)
      indexes(ir)=indexes(1)
      ir=ir-1
      
      if(ir == 1) then
        indexes(1)=id
        return
      end if
    end if
    
    i = l
    j = 2*l
    
    do while( j <= ir )
      if( j < ir ) then
        if( array(indexes(j)) < array(indexes(j+1)) ) then
          j=j+1
        end if
      end if
      
      if( value < array(indexes(j)) ) then
        indexes(i)=indexes(j)
        i=j
        j=2*j
      else
        j=ir+1
      end if
    end do
    
    indexes(i)=id
  end do
end subroutine rsort

!>
!! @brief  Intercambia dos bloques de columnas especificados por los rangos A y B
!!
!! @warning Coloca el bloque de columnas especificado al inicio de la matriz 
!!		el resto de columnas al final de la misma
!! @warning Actualmente no se soportan rangos intermedios, solamente rangos contains
!!			abiertos que incluyan el elemento terminal
!! @todo Dar soporte a rangos no consecutivos
!<
subroutine swapBlockOfColumns( matrix, rangeSpecification )
  implicit none
  real(8), allocatable :: matrix(:,:)
  integer, intent(in) :: rangeSpecification(2)

  real(8), allocatable :: auxMatrix(:,:)

  allocate( auxMatrix(size(matrix,dim=1), size(matrix,dim=2) ) )
  auxMatrix = matrix

  matrix(:, 1: rangeSpecification(2)-rangeSpecification(1)+1) = auxMatrix(:, rangeSpecification(1):rangeSpecification(2) )
  matrix(:, rangeSpecification(2)-rangeSpecification(1)+2:size(matrix,dim=2) ) = auxMatrix(:,1:rangeSpecification(1)-1)

  deallocate(auxMatrix)
end subroutine swapBlockOfColumns


!>
!! Returns the projector of constants of force to make infinitesimal
!! translations and rotations.
!!
subroutine getForceConstantsProjector( projector, nVib, nRotAndTrans )
  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MDARRAY

  implicit none
  real(8), allocatable :: projector(:,:)
  integer :: nVib, nRot
  integer :: nRotAndTrans

  interface
    subroutine getInertiaTensor( inertiaTensor, nRot )
      real(8), intent(out) :: inertiaTensor(3,3)
      integer, optional, intent(out) :: nRot
    end subroutine getInertiaTensor

    subroutine orthogonalizeLinearVectorialSpace( matrix )
      implicit none
      real(8), allocatable :: matrix(:,:)
    end subroutine orthogonalizeLinearVectorialSpace
    
    subroutine swapBlockOfColumns( matrix, rangeSpecification )
      implicit none
      real(8), allocatable :: matrix(:,:)
      integer, intent(in) :: rangeSpecification(2)
    end subroutine swapBlockOfColumns
    
    subroutine rsort( array, indexes )
      real(8), intent(in) :: array(:)
      integer, intent(inout) :: indexes(:)
    end subroutine rsort
  end interface
  
  integer :: i
  integer :: j
  integer :: index_x
  integer :: index_y
  integer :: index_z
  integer :: aux
  logical :: isNull
  real(8) :: sqrtMass
  real(8) :: coordinatesProyected(3)
  real(8) :: geometryInCM(3)
  real(8) :: centerOfMass(3)
  real(8) :: squareNorm
  real(8) :: projectorNorm(6)
  integer :: indexes(6)
!   type(LinearVectorialSpace) :: spaceOfForceConstants
  real(8) :: inertiaTensor(3,3)

  allocate( projector(3*NATS,3*NATS) )
  
  call getCenterOfMass( centerOfMass )
  call getInertiaTensor( inertiaTensor, nRot )
  
  projector = 0.0_8

  do i=1,NATS

    index_x = 3*i - 2
    index_y = 3*i - 1
    index_z = 3*i

    sqrtMass = sqrt( MASS(ELEMPOINTER(i))*1822.88853_8 ) ! amu to a.u.
    geometryInCM = CR(:,i)*1.88972612456506_8-centerOfMass(:)  ! angs to a.u.

    !!
    !! Projects the cartesian coordinates on the inertia tensor
    !!
    coordinatesProyected(1)=dot_product( geometryInCM, inertiaTensor(1,:) )
    coordinatesProyected(2)=dot_product( geometryInCM, inertiaTensor(2,:) )
    coordinatesProyected(3)=dot_product( geometryInCM, inertiaTensor(3,:) )
    
    projector(index_x,1) = sqrtMass
    projector(index_y,2) = sqrtMass
    projector(index_z,3) = sqrtMass

    projector(index_x,4) = 	(coordinatesProyected(2)*inertiaTensor(1,3) &
                - coordinatesProyected(3)*inertiaTensor(1,2) )/sqrtMass
    projector(index_y,4) = 	(coordinatesProyected(2)*inertiaTensor(2,3) &
                - coordinatesProyected(3)*inertiaTensor(2,2) )/sqrtMass
    projector(index_z,4) = 	(coordinatesProyected(2)*inertiaTensor(3,3) &
                - coordinatesProyected(3)*inertiaTensor(3,2) )/sqrtMass

    projector(index_x,5) = 	(coordinatesProyected(3)*inertiaTensor(1,1) &
                - coordinatesProyected(1)*inertiaTensor(1,3) )/sqrtMass
    projector(index_y,5) = 	(coordinatesProyected(3)*inertiaTensor(2,1) &
                - coordinatesProyected(1)*inertiaTensor(2,3) )/sqrtMass
    projector(index_z,5) = 	(coordinatesProyected(3)*inertiaTensor(3,1) &
                - coordinatesProyected(1)*inertiaTensor(3,3) )/sqrtMass

    projector(index_x,6) =	(coordinatesProyected(1)*inertiaTensor(1,2) &
                - coordinatesProyected(2)*inertiaTensor(1,1) )/sqrtMass
    projector(index_y,6) = 	(coordinatesProyected(1)*inertiaTensor(2,2) &
                - coordinatesProyected(2)*inertiaTensor(2,1) )/sqrtMass
    projector(index_z,6) =	(coordinatesProyected(1)*inertiaTensor(3,2) &
                - coordinatesProyected(2)*inertiaTensor(3,1) )/sqrtMass
                
  end do
  
  write(*,*) ""
  write(*,*) "Very initial projector: "
  do i=1,size(projector,dim=2)
    write(*,"(6F8.2)") projector(i,:)
  end do
  write(*,*) ""

  !! Verfies if the six vectors are actually rotational
  !! and translational normal modes
  nRotAndTrans = 0
  isNull=.false.
  aux = 0
  
  write(*,*) "projectorNorm"
  projectorNorm = 0.0_8
  do i=1,6
    projectorNorm(i) = dot_product( projector(:,i),projector(:,i) )
    write(*,*) i, projectorNorm(i)
  end do
  
  call rsort( projectorNorm, indexes )
  
  write(*,*) "projectorNorm(sorted)"
  do i=1,6
    write(*,*) i, projectorNorm( indexes(i) ), projectorNorm( indexes(i) )/sum(projectorNorm)
  end do
  write(*,*) ""

!   do i=1,6
  do i=1,size(projector(1,:))

    squareNorm = dot_product( projector(:,i),projector(:,i) )

    if ( squareNorm > 1.0d-6 ) then

      projector(:,i) = projector(:,i) / sqrt( squareNorm )
      nRotAndTrans = nRotAndTrans + 1
      
      if ( isNull ) then
        projector(:,i-aux) = projector(:,i)
        projector(:,i) = 0.0_8
      end if
      
      write(*,*) i, "real", squareNorm
    else
      
      isNull = .true.
      aux = aux+1
      
      write(*,*) i, "fake", squareNorm

    end if

  end do
  
  !!
  !!***********************************************************************
  nVib = 3*NATS - nRotAndTrans
  !!Adiciona una serie de vectores asociados a vibraciones con el fin
  !! de completar la matriz de transformacion
  j=1
  do i=nRotAndTrans+1,3*NATS
    projector(j,i)=1.0_8
    j=j+1
  end do
  
  write(*,*) ""
  write(*,*) "Initial projector: "
  do i=1,size(projector,dim=2)
    write(*,"(6F8.2)") projector(i,:)
  end do
  write(*,*) ""

  !! Construye un espacio vectorial lineal con los vectores previemente generados
  !!**********
  !! 	Proyecta los vectores asociados a grados de libertad vibracionales sobre los 
  !! 	asociados a grados de libertad rotacionales y traslacionales. - Los primeros 
  !!	N vectores (nRotAndTrans ),  se asocian a los grados 
  !!	de libertad rotacionales y translacionales el resto se asocian a los grados de libertad
  !!	vibracionales -
  !!**
  call orthogonalizeLinearVectorialSpace( projector )
  
  write(*,*) ""
  write(*,*) "Orthogonalized projector: "
  do i=1,size(projector,dim=2)
    write(*,"(6F8.2)") projector(i,:)
  end do
  write(*,*) ""

  !!
  !! Reordena el proyector colocando los vectores asociados al vibraciones al principio y los
  !! asociados a rotaciones y traslaciones al final
  !!
  call swapBlockOfColumns( projector, [nRotAndTrans+1, 3*NATS] )
  
  write(*,*) ""
  write(*,*) "Orthogonalized & sorted projector: "
  do i=1,size(projector,dim=2)
    write(*,"(6F8.2)") projector(i,:)
  end do
  write(*,*) ""
  write(*,*) "nRotAndTrans = ", nRotAndTrans
  write(*,*) "nVib = ", nVib
  write(*,*) ""

end subroutine getForceConstantsProjector

!>
!!  @brief Calculates the core of the transformation matrix to remove the external
!!         degrees of freedom from the gradient and the hessian
!<
subroutine getTransformationMatrix( output )
  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MYPRECISION
  USE MDARRAY

  implicit none
  
  interface
    subroutine getForceConstantsProjector( projector, nVib, nRotAndTrans )
      implicit none
      real(8), allocatable :: projector(:,:)
      integer :: nVib, nRotAndTrans
    end subroutine getForceConstantsProjector
  end interface 
  
  real(8), allocatable :: output(:,:)
  real(8), allocatable :: forceConstansProjector(:,:)
  integer :: nVib, nRotAndTrans
  
  allocate( output(NATS*3,NATS*3) )
  output = 0.0_8
  call getForceConstantsProjector( forceConstansProjector, nVib, nRotAndTrans )

  output = -1.0_8 * matmul( forceConstansProjector(:,1:nVib+1), &
            transpose(forceConstansProjector(:,1:nVib+1) ) )

end subroutine getTransformationMatrix

  !> Translate to geometric center.
  !! \param coords Coordinates of the system (see system_type).
  !! \param lattice_vectors System lattice vectors.
  !! \param origin (min(x),min(y),min(z)) set as the origin of the system.
  !!
  subroutine prg_translatetogeomcandfoldtobox(coords,lattice_vectors,origin)
    implicit none
    integer                              ::  i
    real(8), allocatable, intent(inout) ::  origin(:),coords(:,:)
    real(8), intent(in)                 ::  lattice_vectors(:,:)
    real(8)                             ::  geomc(3)

    if(.not.allocated(origin)) allocate(origin(3))

    ! Getting the geometric center.
    do i=1,size(coords,dim=2)
       geomc = geomc + coords(:,i)
    enddo

    geomc = geomc/size(coords,dim=2)

    coords(1,:) = coords(1,:) - geomc(1)
    coords(2,:) = coords(2,:) - geomc(2)
    coords(3,:) = coords(3,:) - geomc(3)

    coords(1,:) = mod(coords(1,:)*1000000,1000000*lattice_vectors(1,1))/1000000;
    coords(2,:) = mod(coords(2,:)*1000000,1000000*lattice_vectors(2,2))/1000000;
    coords(3,:) = mod(coords(3,:)*1000000,1000000*lattice_vectors(3,3))/1000000;

    origin = 0.0

    !Shift origin slightly
    origin(1) = -1.0d-1 ; origin(2) = -1.0d-1; origin(3) = -1.0d-1

  end subroutine prg_translatetogeomcandfoldtobox

 !> Wrap around atom i using pbc.
  !! \param coords Coordinates of the system (see system_type).
  !! \param lattice_vectors System lattice vectors.
  !! \param index Index atom to wrap around
  !!
  subroutine prg_wraparound(coords,lattice_vectors,index,verbose)
    implicit none
    integer                              ::  i, nats
    integer, intent(in)                  ::  index
    integer, intent(in), optional        ::  verbose
    real(8), allocatable, intent(inout) ::  coords(:,:)
    real(8), allocatable                ::  origin(:)
    real(8), intent(in)                 ::  lattice_vectors(:,:)

    if(present(verbose) .and. verbose >= 1)write(*,*)"In prg_wraparound ..."

    if(.not.allocated(origin)) allocate(origin(3))

    nats=size(coords,dim=2)

    origin(1) = -coords(1,index) + lattice_vectors(1,1)/2.0d0
    origin(2) = -coords(2,index) + lattice_vectors(2,2)/2.0d0
    origin(3) = -coords(3,index) + lattice_vectors(3,3)/2.0d0

    coords(1,:) = coords(1,:) + origin(1)
    coords(2,:) = coords(2,:) + origin(2)
    coords(3,:) = coords(3,:) + origin(3)

    !$omp parallel do default(none) private(i) &
    !$omp shared(coords,lattice_vectors,nats)
    do i=1,nats
       if(coords(1,i) > lattice_vectors(1,1))coords(1,i)=coords(1,i)-lattice_vectors(1,1)
       if(coords(2,i) > lattice_vectors(2,2))coords(2,i)=coords(2,i)-lattice_vectors(2,2)
       if(coords(3,i) > lattice_vectors(3,3))coords(3,i)=coords(3,i)-lattice_vectors(3,3)
       if(coords(1,i) < 0.0d0)coords(1,i)=coords(1,i)+lattice_vectors(1,1)
       if(coords(2,i) < 0.0d0)coords(2,i)=coords(2,i)+lattice_vectors(2,2)
       if(coords(3,i) < 0.0d0)coords(3,i)=coords(3,i)+lattice_vectors(3,3)
    enddo
    !$end omp parallel do

  end subroutine prg_wraparound

  !> Translate geometric center to the center of the box.
  !! \param coords Coordinates of the system (see system_type).
  !! \param lattice_vectors System lattice vectors.
  !! \param verbose Verbosity level.
  !!
  subroutine prg_centeratbox(coords,lattice_vectors,verbose)
    implicit none
    integer                              ::  i, nats
    integer, intent(in), optional        ::  verbose
    real(8)                             ::  gc(3)
    real(8), allocatable, intent(inout) ::  coords(:,:)
    real(8), intent(in)                 ::  lattice_vectors(:,:)

    if(present(verbose) .and. verbose >= 1)write(*,*)"In prg_centeratbox ..."

    nats=size(coords,dim=2)

    gc= 0.0d0

    !$omp parallel do default(none) private(i) &
    !$omp shared(coords,nats) &
    !$omp reduction(+:gc)
    do i=1,nats
       gc=gc + coords(:,i)
    enddo
    !$omp end parallel do

    gc=gc/real(nats,8)

    !$omp parallel do default(none) private(i) &
    !$omp shared(coords,lattice_vectors,nats, gc)
    do i=1,nats
       coords(1,i) = coords(1,i) + lattice_vectors(1,1)/2.0d0 - gc(1)
       coords(2,i) = coords(2,i) + lattice_vectors(2,2)/2.0d0 - gc(2)
       coords(3,i) = coords(3,i) + lattice_vectors(3,3)/2.0d0 - gc(3)
    enddo
    !$omp end parallel do

  end subroutine prg_centeratbox
  
!!>
!! Prepares the geometry for IR calculation
!!
SUBROUTINE PREPAREIR

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MYPRECISION
  USE MDARRAY

  implicit none
  
  interface
    subroutine prg_wraparound(coords,lattice_vectors,index,verbose)
      implicit none
      integer, intent(in)                  ::  index
      integer, intent(in), optional        ::  verbose
      real(8), allocatable, intent(inout) ::  coords(:,:)
      real(8), intent(in)                 ::  lattice_vectors(:,:)
    end subroutine prg_wraparound
    subroutine prg_centeratbox(coords,lattice_vectors,verbose)
      implicit none
      integer, intent(in), optional        ::  verbose
      real(8), allocatable, intent(inout) ::  coords(:,:)
      real(8), intent(in)                 ::  lattice_vectors(:,:)
    end subroutine prg_centeratbox
  end interface 
  
  call prg_wraparound(CR,BOX,1,1)
  call prg_centeratbox(CR,BOX,1)
  
END SUBROUTINE PREPAREIR  

!!>
!! Calculates the IR spectrum
!!
SUBROUTINE GETIR

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MYPRECISION
  USE MDARRAY

  implicit none
  
  interface
    subroutine getTransformationMatrix( output )
      implicit none
      real(8), allocatable :: output(:,:)
    end subroutine getTransformationMatrix
    subroutine prg_translatetogeomcandfoldtobox(coords,lattice_vectors,origin)
      implicit none
      real(8), allocatable, intent(inout) ::  origin(:),coords(:,:)
      real(8), intent(in)                 ::  lattice_vectors(:,:)
      real(8)                             ::  geomc(3)
    end subroutine prg_translatetogeomcandfoldtobox
    subroutine prg_wraparound(coords,lattice_vectors,index,verbose)
      implicit none
      integer, intent(in)                  ::  index
      integer, intent(in), optional        ::  verbose
      real(8), allocatable, intent(inout) ::  coords(:,:)
      real(8), intent(in)                 ::  lattice_vectors(:,:)
    end subroutine prg_wraparound
    subroutine prg_centeratbox(coords,lattice_vectors,verbose)
      implicit none
      integer, intent(in), optional        ::  verbose
      real(8), allocatable, intent(inout) ::  coords(:,:)
      real(8), intent(in)                 ::  lattice_vectors(:,:)
    end subroutine prg_centeratbox
  end interface 
  
  integer :: i, j, p1, p2, atom1, atom2
  real(8) :: m1, m2
  REAL(8) :: d2Vdp1dp2, HB
  real(8), allocatable :: geom0(:,:)
  
  real(8) :: inertiaTensor(3,3)
  
  real(8), allocatable :: hessian(:,:)
  real(8), allocatable :: eVecsHessian(:,:)
  real(8), allocatable :: eValsHessian(:)
  
  real(8), allocatable :: transformationMatrix(:,:)
  real(8), allocatable :: workSpace(:)
  integer :: ssize, info
  
  IF (EXISTERROR) RETURN
 
  write(*,*) ""
  write(*,*) "Original geometry:"
  do atom1=1,NATS
    write(*,"(F8.2,4X,3F8.2)") MASS(ELEMPOINTER(atom1)), CR(:,atom1)
  end do
  write(*,*) ""
  
  call getTransformationMatrix( transformationMatrix )
  
  write(*,*) " transformationMatrix : "
  do i=1,size(transformationMatrix,dim=2)
    write(*,"(6F8.2)") transformationMatrix(i,:)
  end do
  
!   do i=1,size(transformationMatrix,dim=1)
!     transformationMatrix(i,i) = 1.0_8 + transformationMatrix(i,i)
!   end do

  CALL GETFORCE
  
  write(*,*) " Norm. grad = ", sqrt(sum(FTOT**2))
  
  allocate( geom0(3,NATS) )
  allocate( hessian(3*NATS,3*NATS) )
  
  geom0 = CR
  
  hb = 0.001d0

  p1 = 1
  do atom1=1,NATS; do i=1,3
    m1 = MASS(ELEMPOINTER(atom1))*1822.88853_8 ! amu to a.u.
    
!     write(*,*) MASS(ELEMPOINTER(atom1)), CR(I,atom1)
    
    p2 = 1
    do atom2=1,NATS; do j=1,3
      m2 = MASS(ELEMPOINTER(atom2))*1822.88853_8 ! amu to a.u.
      
      FTOT = 0.0d0
      CR = geom0
      
      !! Second derivatives are calculated as follows
      !! a,b = atom1, atom2
      !! i,j = x, y, or z
      !!
      !! dV(ai,bj)/dai = -F(ai,bj)
      !!
      !! d2V(ai,bj)/daidbj = d( dV(ai,bj)/dai )/dbj = d( -F(ai,bj) )/dbj
      !!
      !! d2V(ai,bj)/daidbj ~ ( -F(ai,bj+hb) + F(ai,bj-hb) )/(2*hb)
      
      CR(J,atom2) = CR(J,atom2) + hb
      CALL GETFORCE
      
      d2Vdp1dp2 = -FTOT(I,atom1)
      
      CR(J,atom2) = CR(J,atom2) - 2*hb
      CALL GETFORCE
      
      d2Vdp1dp2 = ( d2Vdp1dp2 + FTOT(I,atom1) )/2.0d0/hb ! Derivative
      d2Vdp1dp2 = d2Vdp1dp2*0.0102908545816127_8  ! eV/angs^2 to a.u.
      
      d2Vdp1dp2 = d2Vdp1dp2/sqrt(m1*m2) ! Mass weighted derivative
      
      hessian( p1, p2 ) = d2Vdp1dp2
        
      p2 = p2 + 1
    end do; end do
    
    p1 = p1 + 1
  end do; end do
  
  write(*,*) "Hessian:"
  do i=1,3*NATS
    do j=1,3*NATS
      write(*,"(E20.2)",advance="no") hessian(i,j)
    end do
    write(*,*) ""
  end do
  
  hessian = matmul( transpose(transformationMatrix), matmul( hessian, transformationMatrix ) )
  
  CR = geom0 ! Restore original geometry
  
!   write(*,*) hessian
  allocate( eValsHessian(3*NATS) )
  allocate( eVecsHessian(3*NATS,3*NATS) )
  
  write(*,*) "D^T*Hessian*D:"
  do i=1,3*NATS
    do j=1,3*NATS
      write(*,"(E20.2)",advance="no") hessian(i,j)
    end do
    write(*,*) ""
  end do
  
  allocate( workSpace( 3*3*NATS-1 ) )
  eVecsHessian = hessian

  ! Compute the eigen values and eigen vectors using the upper elements of the symmetric matrix
  call dsyev( 'V', 'L', 3*NATS, eVecsHessian, 3*NATS, eValsHessian, workSpace, 3*3*NATS-1, info )
  
  deallocate( workSpace )
  
  if ( info /= 0 ) then
    write(*,*) "### ERROR ### GETIR.dsyev: matrix diagonalization failed"
    stop
  end if
  
  write(*,*) ""
  write(*,*) " Vib. Freqs (cm-1) "
  write(*,*) "-------------------"
  do i=1,3*NATS
    if( eValsHessian(i) < 0.0_8 ) then
      write(*,"(I10,F20.10)") i, -sqrt(abs(eValsHessian(i)))*219474.63068_8  ! a.u. to cm-1  
    else
      write(*,"(I10,F20.10)") i, sqrt(eValsHessian(i))*219474.63068_8  ! a.u. to cm-1  
    end if
  end do
  write(*,*) ""

!   WRITE(*,*)"==========================="
!   WRITE(*,*)"Im stoping in subroutine getir"
!   WRITE(*,*)"Grep for GETIR"
!   WRITE(*,*)"==========================="
  
  deallocate( geom0 )
  deallocate( hessian )
  deallocate( eValsHessian )
  deallocate( eVecsHessian )
  
  STOP

  RETURN

END SUBROUTINE GETIR
