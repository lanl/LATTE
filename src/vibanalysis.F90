!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 2019.  Los Alamos National Security, LLC. This material was    !
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

!!>
!! Returns the center of mass
!!
subroutine getCenterOfMass( centerOfMass )
  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MDARRAY
  USE MYPRECISION

  implicit none

  real(latteprec) :: centerOfMass(3)

  centerOfMass(1) = sum(MASS(ELEMPOINTER(:))*CR(1,:))
  centerOfMass(2) = sum(MASS(ELEMPOINTER(:))*CR(2,:))
  centerOfMass(3) = sum(MASS(ELEMPOINTER(:))*CR(3,:))
  centerOfMass = centerOfMass / sum(MASS(ELEMPOINTER(:)))

  centerOfMass = centerOfMass*1.88972612456506_latteprec ! angs to a.u.

end subroutine getCenterOfMass

!!>
!! Gets the inertia tensor and its eigen vectors and eigen values
!!
subroutine getAxesOfInertia( axesOfInertia, nRot )
  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MDARRAY
  USE MYPRECISION

  implicit none
  
  real(latteprec), parameter :: amu = 1822.88853_latteprec

  real(latteprec), intent(out) :: axesOfInertia(3,3)
  integer, optional, intent(out) :: nRot

  real(latteprec) :: inertiaTensor(3,3)
  real(latteprec) :: inertiaMoments(3)

  integer :: i, j, atom1, atom2
  real(latteprec), allocatable :: geometryInCM(:,:)
  real(latteprec) :: centerOfMass(3)
  real(latteprec) :: massi

  integer :: ssign

  real(latteprec), allocatable :: workSpace(:)
  integer :: ssize, info

  ALLOCATE( geometryInCM(3,NATS) )

  call getCenterOfMass( centerOfMass )

  do i=1,NATS
    geometryInCM(:,i) = CR(:,i)*1.88972612456506_latteprec - centerOfMass(:)   ! angs to a.u.
  end do

  inertiaTensor = 0.0_latteprec
  do i=1,NATS
    massi = real(MASS(ELEMPOINTER(i)),latteprec)*amu ! amu to a.u.
    
    inertiaTensor(1,1) = inertiaTensor(1,1) + massi * ( geometryInCM(2,i)**2 + geometryInCM(3,i)**2)
    inertiaTensor(2,2) = inertiaTensor(2,2) + massi * ( geometryInCM(1,i)**2 + geometryInCM(3,i)**2)
    inertiaTensor(3,3) = inertiaTensor(3,3) + massi * ( geometryInCM(1,i)**2 + geometryInCM(2,i)**2)
    
    inertiaTensor(1,2) = inertiaTensor(1,2) - massi * ( geometryInCM(1,i) * geometryInCM(2,i) )
    inertiaTensor(1,3) = inertiaTensor(1,3) - massi * ( geometryInCM(1,i) * geometryInCM(3,i) )
    inertiaTensor(2,3) = inertiaTensor(2,3) - massi * ( geometryInCM(2,i) * geometryInCM(3,i) )
  end do

  inertiaTensor(2,1) =inertiaTensor(1,2)
  inertiaTensor(3,1) =inertiaTensor(1,3)
  inertiaTensor(3,2) =inertiaTensor(2,3)

  axesOfInertia = inertiaTensor
  allocate( workSpace( 3*3*3-1 ) )

  ! Compute the eigen values and eigen vectors using the upper elements of the symmetric matrix
  call dsyev( 'V', 'L', 3, axesOfInertia, 3, inertiaMoments, workSpace, 3*3*3-1, info )

  if( info /= 0 ) then
    write(*,*) "### ERROR ### Diagonalizing the inertia tensor"
    stop
  end if

  deallocate( workSpace )

!   !! Checks the determinant's sign
!   ssign = axesOfInertia(1,1)*( &
!           axesOfInertia(2,2)*axesOfInertia(3,3) &
!           -axesOfInertia(3,2)*axesOfInertia(2,3)) &
!           -axesOfInertia(1,2)*( &
!           axesOfInertia(2,1)*axesOfInertia(3,3) &
!           -axesOfInertia(3,1)*axesOfInertia(2,3)) &
!           +axesOfInertia(1,3)*( &
!           axesOfInertia(2,1)*axesOfInertia(3,2) &
!           -axesOfInertia(3,1)*axesOfInertia(2,2))
! 
!   !! Preserbs the handedness of the inertia tensor
!   if ( ssign < 0.0 ) then
!     axesOfInertia(1,2) = -axesOfInertia(1,2)
!     axesOfInertia(2,2) = -axesOfInertia(2,2)
!     axesOfInertia(3,2) = -axesOfInertia(3,2)
!   endif

  !! Verifies if the inertia tensor is correct
  if ( 	abs( inertiaTensor(1,1) ) < 1d-10 .and. &
        abs( inertiaTensor(2,2) ) < 1d-10 .and. &
        abs( inertiaTensor(3,3) ) < 1d-10 ) then
    write(*,*) "### ERROR ### Invalid inertia tensor."
	write(*,*) ""
	write(*,*) "inertiaTensor: "
	do i=1,size(inertiaTensor,dim=1)
		do j=1,size(inertiaTensor,dim=2)
			write(*,"(F8.2)",advance="no") inertiaTensor(i,j)
		end do
		write(*,*) ""
	end do
	write(*,*) ""
    stop
  end if

  if( present(nRot) ) then
    nRot = 3
    if( abs( inertiaMoments(1) ) < 10.0 ) then
      nRot = 2 ! @todo Here is important make a second check to verify that molecule is indeed lineal
    end if
    write(*,*) ""
    write(*,*) "nRot = ", nRot
  end if

  write(*,*) ""
  write(*,"(A,3F20.5)") " Principal Moments of Inertia: ", inertiaMoments
  write(*,*) ""
  write(*,*) "Principal Axes of Inertia: "
  do i=1,3
  	write(*,"(3F15.2)") axesOfInertia(:,i)
  end do
  
  write(*,*) ""
  write(*,*) "Coordinates in the Center of Mass:"
  do atom1=1,NATS
  	write(*,"(F8.2,4X,3F8.2)") MASS(ELEMPOINTER(atom1)), geometryInCM(:,atom1)
  end do
  write(*,*) ""

  deallocate( geometryInCM )
  
end subroutine getAxesOfInertia

!>
!! @brief 
!<
subroutine projectLastElements( matrix, output )
  USE MYPRECISION
  implicit none
  real(latteprec) :: matrix(:,:)
  real(latteprec), allocatable :: output(:,:)

  integer :: i
  integer :: last
  real(latteprec) :: squareNorm
  real(latteprec) :: projectionOverOrthogonalizedBasis

  last = size(matrix,dim=2)
  allocate( output( size(matrix,dim=1), last ) )
  output = matrix

  do i=1,last-1
    squareNorm = dot_product( output(:,i), output(:,i) )
    projectionOverOrthogonalizedBasis = dot_product( output(:,i), output(:,last) )

    if ( squareNorm > 1.0d-12 ) then
      output( :, last ) = output( :, last ) - ( projectionOverOrthogonalizedBasis/squareNorm )*output(:,i)
    end if
  end do
end subroutine projectLastElements

!>
!! @brief Orthogonalize the matrix following the Gram–Schmidt method
!<
subroutine orthogonalizeMatrix( matrix )
  USE MYPRECISION
  implicit none
  real(latteprec), allocatable :: matrix(:,:)

  interface
    subroutine projectLastElements( matrix, output )
      USE MYPRECISION
      implicit none
      real(latteprec) :: matrix(:,:)
      real(latteprec), allocatable :: output(:,:)
    end subroutine projectLastElements
  end interface

  integer :: i, j, k
  integer :: last
  real(latteprec) :: norm
  real(latteprec), allocatable :: output(:,:)

  last = size(matrix,dim=2)
  norm=sqrt(dot_product(matrix(:,1),matrix(:,1)))

  matrix(:,1)=matrix(:,1)/norm

  do i=2,last
    call projectLastElements( matrix(:,1:i), output )
    matrix(:,1:i) = output
    if( allocated(output) ) deallocate(output)
  end do
  
!   do i=2,last
!     call projectLastElements( matrix(:,1:i), output )
!     matrix(:,1:i) = output
!     if( allocated(output) ) deallocate(output)
!   end do

end subroutine orthogonalizeMatrix

!>
!! @brief Returns the transformation matrix that transforms from mass-weighted cartesian coordinates
!!        to internal coordinates, where rotation and translation have been separated out.
!!
subroutine getProjector( projector, nVib, nRotAndTrans )
  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MDARRAY
  USE MYPRECISION

  implicit none
  real(latteprec), allocatable :: projector(:,:)
  integer :: nVib, nRot
  integer :: nRotAndTrans

  interface
    subroutine getAxesOfInertia( axesOfInertia, nRot )
      USE MYPRECISION
      real(latteprec), intent(out) :: axesOfInertia(3,3)
      integer, optional, intent(out) :: nRot
    end subroutine getAxesOfInertia

    subroutine orthogonalizeMatrix( matrix )
      USE MYPRECISION
      implicit none
      real(latteprec), allocatable :: matrix(:,:)
    end subroutine orthogonalizeMatrix
  end interface

  integer :: i, j
  integer :: index_x, index_y, index_z
  integer :: aux
  logical :: isNull
  real(latteprec) :: sqrtMass
  real(latteprec) :: coordinatesProyected(3)
  real(latteprec) :: geometryInCM(3)
  real(latteprec) :: centerOfMass(3)
  real(latteprec) :: squareNorm
  real(latteprec) :: axesOfInertia(3,3)

  allocate( projector(3*NATS,3*NATS) )

  call getCenterOfMass( centerOfMass )
  call getAxesOfInertia( axesOfInertia, nRot )

  projector = 0.0_latteprec

  do i=1,NATS

    index_x = 3*i - 2
    index_y = 3*i - 1
    index_z = 3*i

    sqrtMass = sqrt( MASS(ELEMPOINTER(i))*1822.88853_latteprec ) ! amu to a.u.
    geometryInCM = CR(:,i)*1.88972612456506_latteprec-centerOfMass(:)  ! angs to a.u.

    !!
    !! Projects the cartesian coordinates on the inertia axes
    !!
    coordinatesProyected(1) = dot_product( geometryInCM, axesOfInertia(:,1) )
    coordinatesProyected(2) = dot_product( geometryInCM, axesOfInertia(:,2) )
    coordinatesProyected(3) = dot_product( geometryInCM, axesOfInertia(:,3) )

    projector(index_x,1) = sqrtMass
    projector(index_y,2) = sqrtMass
    projector(index_z,3) = sqrtMass

    projector(index_x,4) = ( coordinatesProyected(2)*axesOfInertia(1,3) - coordinatesProyected(3)*axesOfInertia(1,2) )*sqrtMass
    projector(index_y,4) = ( coordinatesProyected(2)*axesOfInertia(2,3) - coordinatesProyected(3)*axesOfInertia(2,2) )*sqrtMass
    projector(index_z,4) = ( coordinatesProyected(2)*axesOfInertia(3,3) - coordinatesProyected(3)*axesOfInertia(3,2) )*sqrtMass

    projector(index_x,5) = ( coordinatesProyected(3)*axesOfInertia(1,1) - coordinatesProyected(1)*axesOfInertia(1,3) )*sqrtMass
    projector(index_y,5) = ( coordinatesProyected(3)*axesOfInertia(2,1) - coordinatesProyected(1)*axesOfInertia(2,3) )*sqrtMass
    projector(index_z,5) = ( coordinatesProyected(3)*axesOfInertia(3,1) - coordinatesProyected(1)*axesOfInertia(3,3) )*sqrtMass

    projector(index_x,6) = ( coordinatesProyected(1)*axesOfInertia(1,2) - coordinatesProyected(2)*axesOfInertia(1,1) )*sqrtMass
    projector(index_y,6) = ( coordinatesProyected(1)*axesOfInertia(2,2) - coordinatesProyected(2)*axesOfInertia(2,1) )*sqrtMass
    projector(index_z,6) = ( coordinatesProyected(1)*axesOfInertia(3,2) - coordinatesProyected(2)*axesOfInertia(3,1) )*sqrtMass

  end do

  if( VERBOSE >= 1 ) then
	write(*,*) ""
	write(*,*) "Initial projector: "
	do i=1,size(projector,dim=1)
		do j=1,size(projector,dim=2)
			write(*,"(F8.2)",advance="no") projector(i,j)
		end do
		write(*,*) ""
	end do
	write(*,*) ""
  end if

  nRotAndTrans = 0
  isNull=.false.
  aux = 0
  
  ! Rotations and translations are detected and placed at the begining of the transformation matrix
  do i=1,size(projector(1,:))

    squareNorm = dot_product( projector(:,i),projector(:,i) )

    if ( squareNorm > 1.0d-4 ) then

      projector(:,i) = projector(:,i) / sqrt( squareNorm )
      nRotAndTrans = nRotAndTrans + 1

      if ( isNull ) then
        projector(:,i-aux) = projector(:,i)
        projector(:,i) = 0.0_latteprec
      end if

      if( VERBOSE >= 1 ) write(*,*) i, "real", squareNorm
    else

      isNull = .true.
      aux = aux+1

      if( VERBOSE >= 1 ) write(*,*) i, "fake", squareNorm

    end if

  end do
  
  if( VERBOSE >= 1 ) then
	write(*,*) ""
	write(*,*) "Normalized Initial projector: "
	do i=1,size(projector,dim=1)
		do j=1,size(projector,dim=2)
			write(*,"(F8.2)",advance="no") projector(i,j)
		end do
		write(*,*) ""
	end do
  end if

  nVib = 3*NATS - nRotAndTrans
  
  !! Some vectors associated with vibrations are added at the end to compleate the transformation matrix
  j=1
  do i=nRotAndTrans+1,3*NATS
    projector(j,i)=1.0_latteprec
    j=j+1
  end do

  if( VERBOSE >= 1 ) then
	write(*,*) ""
	write(*,*) "Projector: "
	do i=1,size(projector,dim=1)
		do j=1,size(projector,dim=2)
			write(*,"(F8.2)",advance="no") projector(i,j)
		end do
		write(*,*) ""
	end do
  end if

  !! A Schmidt orthogonalization is used to generate Nvib=3N-6 (or 3N-5) remaining vectors, which
  !! are orthogonal to the five or six rotational and translational vectors. 
  call orthogonalizeMatrix( projector )

  if( VERBOSE >= 1 ) then
	write(*,*) ""
	write(*,*) "Orthogonalized Projector: "
	do i=1,size(projector,dim=1)
		do j=1,size(projector,dim=2)
			write(*,"(F8.2)",advance="no") projector(i,j)
		end do
		write(*,*) ""
	end do
  end if
  
  write(*,*) ""
  write(*,*) "nRotAndTrans = ", nRotAndTrans
  write(*,*) "nVib = ", nVib
  write(*,*) ""
end subroutine getProjector

!>
!!  @brief Calculates the core of the transformation matrix to remove the external
!!         degrees of freedom from the hessian
!<
subroutine getTransformationMatrix( output, nVib )
  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MYPRECISION
  USE MDARRAY
  USE MYPRECISION

  implicit none
  
  integer :: nVib

  interface
    subroutine getProjector( projector, nVib, nRotAndTrans )
      USE MYPRECISION
      implicit none
      real(latteprec), allocatable :: projector(:,:)
      integer :: nVib, nRotAndTrans
    end subroutine getProjector
  end interface

  real(latteprec), allocatable :: output(:,:)
  real(latteprec), allocatable :: projector(:,:)
  integer :: i, nRotAndTrans

  allocate( output(3*NATS,3*NATS) )
  output = 0.0_latteprec
  call getProjector( projector, nVib, nRotAndTrans )

  output = -1.0_latteprec * matmul( projector, transpose(projector) )
            
  do i=1,size(output,dim=1)
    output(i,i) = 1.0_latteprec + output(i,i)
  end do

end subroutine getTransformationMatrix

!!>
!! Calculates the dipole
!!
SUBROUTINE GETDIPOLEVEC(DIPOLEVEC)

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE COULOMBARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I
  REAL(LATTEPREC), INTENT(OUT) :: DIPOLEVEC(3)
  IF (EXISTERROR) RETURN

  DIPOLEVEC = ZERO

  DO I = 1, NATS

     DIPOLEVEC(1) = DIPOLEVEC(1) + CR(1,I)*DELTAQ(I)
     DIPOLEVEC(2) = DIPOLEVEC(2) + CR(2,I)*DELTAQ(I)
     DIPOLEVEC(3) = DIPOLEVEC(3) + CR(3,I)*DELTAQ(I)

  ENDDO
  
  RETURN

END SUBROUTINE GETDIPOLEVEC

!!>
!! Saves the vibrational analysis in MOLDEN format
!!
subroutine saveMolden( symbols, coords, freqs, intensities, lCart, nVib )
  USE MYPRECISION
  USE CONSTANTS_MOD
  implicit none
  
  real(latteprec), parameter :: cm1 = 4.5563352670953d-06
  real(latteprec), parameter :: angs = 1.88972612456506_latteprec
  
  character(*), allocatable :: symbols(:)
  real(latteprec), allocatable :: coords(:,:)
  real(latteprec), allocatable :: freqs(:)
  real(latteprec), allocatable :: intensities(:)
  real(latteprec), allocatable :: lCart(:,:)
  integer :: nVib
  
  integer :: i, j, k, l, unit
  
  if( len(trim(VIBSAVE)) == 0 ) return
  
  unit = 333
  open( unit=unit, file=trim(VIBSAVE), action="write" )
  
  write(unit,*) "[Molden Format]"
  write(unit,*) "[GEOMETRIES] XYZ"
  write(unit,"(I5)") size(coords,dim=2)
  write(unit,*) ""
  do i=1,size(coords,dim=2)
    write(unit,"(A5,3F15.8)") symbols(i), coords(:,i) ! angs
  end do
  
  write(unit,*) "[FREQ]"
  do i=size(freqs,dim=1)-nVib+1,size(freqs,dim=1)
!   do i=1,size(freqs,dim=1)
    write(unit,"(F15.8)") freqs(i)/cm1 ! a.u. to cm-1
  end do
  
  write(unit,*) "[INT]"
  do i=size(freqs,dim=1)-nVib+1,size(freqs,dim=1)
!   do i=1,size(freqs,dim=1)
    write(unit,"(2F15.8)") intensities(i), 0.0_latteprec
  end do
  
  write(unit,*) "[FR-COORD]"

  do i=1,size(coords,dim=2)
    write(unit,"(A5,3F15.8)") symbols(i), coords(:,i)*angs ! angs to a.u.
  end do
  
  write(unit,*) "[FR-NORM-COORD]"
  l = 1
  do j=size(freqs,dim=1)-nVib+1,size(freqs,dim=1)
!   do j=1,size(lCart,dim=2)
    write(unit,*) "vibration", l
    i = 1
    do while( i<=size(lCart,dim=1) )
      do k=1,3
        write(unit,"(F8.2)",advance="no") lCart(i,j)
        i = i + 1
      end do
      write(unit,*) ""
    end do
    l = l + 1
  end do
  
  close( unit )
end subroutine saveMolden

!!>
!! @brief Carries out the vibrational analysis
!!
SUBROUTINE VIBRATIONAL_ANALYSIS

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MYPRECISION
  USE MDARRAY
  USE MYPRECISION
  USE TIMER_MOD

  implicit none

  interface
    subroutine getTransformationMatrix( output, nVib )
      USE MYPRECISION
      implicit none
      real(latteprec), allocatable :: output(:,:)
      integer :: nVib
    end subroutine getTransformationMatrix
	subroutine saveMolden( symbols, coords, freqs, intensities, lCart, nVib )
	  USE MYPRECISION
	  implicit none
	  character(*), allocatable :: symbols(:)
	  real(latteprec), allocatable :: coords(:,:)
	  real(latteprec), allocatable :: freqs(:)
	  real(latteprec), allocatable :: intensities(:)
	  real(latteprec), allocatable :: lCart(:,:)
	  integer :: nVib
	end subroutine saveMolden
	SUBROUTINE GETDIPOLEVEC(DIPOLEVEC)
	  USE MYPRECISION
	  implicit none
	  REAL(LATTEPREC), intent(out) :: DIPOLEVEC(3)
	end SUBROUTINE GETDIPOLEVEC
  end interface
  
  real(latteprec), parameter :: cm1 = 4.5563352670953d-06
  real(latteprec), parameter :: eV = 0.0367493088244753_latteprec
  real(latteprec), parameter :: amu = 1822.88853_latteprec
  real(latteprec), parameter :: angs = 1.88972612456506_latteprec
  real(latteprec), parameter :: debye = 0.393456_latteprec
  real(latteprec), parameter :: kcalMol = 0.00159360121939816_latteprec

  integer :: i, j, p1, p2, atom1, atom2
  real(latteprec) :: m1, m2
  REAL(latteprec) :: d2Vdp1dp2, HB
  real(latteprec), allocatable :: geom0(:,:)

  real(latteprec) :: inertiaTensor(3,3)

  real(latteprec), allocatable :: hessian(:,:)
  real(latteprec), allocatable :: eVecsHessian(:,:)
  real(latteprec), allocatable :: eValsHessian(:)

  real(latteprec), allocatable :: transformationMatrix(:,:)
  real(latteprec), allocatable :: workSpace(:)
  integer :: ssize, info
  
  integer :: nVib
  real(latteprec), allocatable :: freqs(:)
  real(latteprec), allocatable :: lCart(:,:)
  real(latteprec) :: squareNorm
  
  real(latteprec) :: dipoleVec(3)
  real(latteprec), allocatable :: dMudR(:,:)
  real(latteprec), allocatable :: dMudQ(:,:)
  integer, allocatable :: dMudR_check(:)
  real(latteprec), allocatable :: intensities(:)
  
  real(latteprec) :: mlsi

  IF (EXISTERROR) RETURN
  
  if( NATS < 2 ) return

  if( VERBOSE >= 1 ) then
	write(*,*) ""
	write(*,*) "Original Geometry:"
	do atom1=1,NATS
		write(*,"(F8.2,4X,3F8.2)") MASS(ELEMPOINTER(atom1)), CR(:,atom1)
	end do
	write(*,*) ""
  end if

  mlsi = TIME_MLS()
  call getTransformationMatrix( transformationMatrix, nVib )
  write(*,"(A,F10.2,A)") " Transformation Matrix Calculation", TIME_MLS() - MLSI, "ms"

  if( VERBOSE >= 1 ) then
	write(*,*) ""
	write(*,*) "Transformation Matrix:"
	do i=1,size(transformationMatrix,dim=1)
		do j=1,size(transformationMatrix,dim=2)
		write(*,"(F8.2)", advance="no") transformationMatrix(i,j)
		end do
		write(*,*) ""
	end do
  end if

  CALL GETFORCE

  write(*,*) ""
  write(*,*) "Norm of the Gradient = ", sqrt(sum(FTOT**2)), " eV/A"
  
  if( sqrt(sum(FTOT**2)) > 1d-2 ) then
	write(*,*) "@@@ WARNING @@@ This is not a local minimum."
	write(*,*) "                Consider optimizing the geometry to get a reliable vibrational analysis."
  end if
  write(*,*) ""

  allocate( geom0(3,NATS) )
  allocate( hessian(3*NATS,3*NATS) )
  allocate( dMudR(3,3*NATS) )
  allocate( dMudQ(3,3*NATS) )
  allocate( dMudR_check(3*NATS) )
  
  hessian = 0.0_latteprec
  dMudR = 0.0_latteprec
  dMudQ = 0.0_latteprec
  dMudR_check = 0

  geom0 = CR

  hb = 0.001_latteprec ! angs
  
  mlsi = TIME_MLS()
  
  p1 = 1
  do atom1=1,NATS
    m1 = MASS(ELEMPOINTER(atom1))*amu ! amu to a.u.
    do i=1,3

		p2 = 1
		do atom2=1,NATS
			m2 = MASS(ELEMPOINTER(atom2))*amu ! amu to a.u.
			do j=1,3
			
			if( p2 >= p1 .and. norm2( CR(:,atom2)-CR(:,atom1) ) < 4.0_latteprec ) then
				FTOT = 0.0_latteprec
				CR = geom0
		
				!! Second derivatives of the potential are calculated as follows
				!! a,b = atom1, atom2
				!! i,j = x, y, or z
				!!
				!! dV(ai,bj)/dai = -F(ai,bj)
				!!
				!! d2V(ai,bj)/daidbj = d( dV(ai,bj)/dai )/dbj = d( -F(ai,bj) )/dbj = -d( F(ai,bj) )/dbj
				!!
				!! d2V(ai,bj)/daidbj ~ -( -0.5*F(ai,bj-hb) + 0.0*F(ai,bj) + 0.5*F(ai,bj+hb) )/hb
				!!                   ~  ( 0.5*F(ai,bj-hb) - 0.5*F(ai,bj+hb) )/hb
				!!
				!! First derivatives of the dipole are calculated as follows
				!!
				!! dMu(bj)/dbj = ( -0.5*Mu(bj-hb) + 0.0*Mu(bj) + 0.5*Mu(bj+hb) )/hb
				!!             = ( -0.5*Mu(bj-hb) + 0.5*Mu(bj+hb) )/hb
				CR(J,atom2) = CR(J,atom2) + hb
				CALL GETFORCE
				call GETDIPOLEVEC( dipoleVec )
		
				d2Vdp1dp2 = -0.5_latteprec*FTOT(I,atom1)
				if( dMudR_check(p2)==0 ) dMudR(:,p2) = -0.5_latteprec*dipoleVec
		
				CR(J,atom2) = CR(J,atom2) - 2.0_latteprec*hb
				CALL GETFORCE
				call GETDIPOLEVEC( dipoleVec )
		
				d2Vdp1dp2 = ( d2Vdp1dp2 + 0.5_latteprec*FTOT(I,atom1) )/hb ! Derivative
				if( dMudR_check(p2)==0 ) dMudR(:,p2) = ( dMudR(:,p2) + 0.5_latteprec*dipoleVec )/hb ! Derivative
				
				d2Vdp1dp2 = d2Vdp1dp2*(eV/angs**2)  ! eV/angs^2 to a.u.
				d2Vdp1dp2 = d2Vdp1dp2/sqrt(m1*m2) ! Mass weighted derivative
				
				if( dMudR_check(p2)==0 ) then
					dMudR(:,p2) = dMudR(:,p2)*angs ! angs to a.u.
					dMudR(:,p2) = dMudR(:,p2)/sqrt(m2) ! Mass weighted derivative
				end if
		
				hessian( p1, p2 ) = d2Vdp1dp2
				hessian( p2, p1 ) = d2Vdp1dp2 ! The matrix is symmetric
				
				dMudR_check(p2) = 1 ! Ensures that the derivative is calculated only once
				
				if( VERBOSE >= 2 ) then
					write(*,*) "    atom1=", atom1, "atom2=", atom2, i, j, d2Vdp1dp2
				end if
				call flush(6)
			end if

			p2 = p2 + 1
			end do
		end do

		p1 = p1 + 1
	end do
  end do
  
!   dMudR = 0.0_latteprec
!   p1 = 1
!   do atom1=1,NATS
!     m1 = MASS(ELEMPOINTER(atom1))*amu ! amu to a.u.
!     do i=1,3
! 		CR = geom0
! 		
! 		!! First derivatives of the dipole are calculated as follows
! 		!!
! 		!! dMu(bj)/dbj = ( -0.5*Mu(bj-hb) + 0.0*Mu(bj) + 0.5*Mu(bj+hb) )/hb
! 		!!             = ( -0.5*Mu(bj-hb) + 0.5*Mu(bj+hb) )/hb
! 
! 		CR(i,atom1) = CR(i,atom1) + hb
! 		CALL GETFORCE
! 		call GETDIPOLEVEC( dipoleVec )
! 
! 		dMudR(:,p1) = -0.5_latteprec*dipoleVec
! 
! 		CR(i,atom1) = CR(i,atom1) - 2.0_latteprec*hb
! 		CALL GETFORCE
! 		call GETDIPOLEVEC( dipoleVec )
! 
! 		dMudR(:,p1) = ( dMudR(:,p1) + 0.5_latteprec*dipoleVec )/hb   ! Derivative
! 		
! 		dMudR(:,p1) = dMudR(:,p1)*angs ! angs to a.u.
! 		dMudR(:,p1) = dMudR(:,p1)/sqrt(m1) ! Mass weighted derivative
! 		
! 		p1 = p1 + 1
! 	end do
!   end do

  write(*,*) ""
  write(*,"(A,F10.2,A)") " Hessian Matrix Calculation & Dipole Derivative: ", TIME_MLS() - MLSI, "ms"

  if( VERBOSE >= 1 ) then
    write(*,*) ""
	write(*,*) "Hessian: (a.u.)x10^4"
	do i=1,3*NATS
		do j=1,3*NATS
			write(*,"(F10.4)",advance="no") hessian(i,j)*10000.0
		end do
		write(*,*) ""
	end do
  
    write(*,*) ""
    write(*,"(A,3F10.5,A)") " Dipole= ", dipoleVec, " a.u."
    write(*,*) ""
    write(*,*) "Dipole derivative: a.u."
    do i=1,3*NATS
    	do j=1,3
  		write(*,"(F10.4)",advance="no") dMudR(j,i)
    	end do
    	write(*,*) ""
    end do
  end if
  
!   hessian = matmul( transpose(transformationMatrix), matmul( hessian, transformationMatrix ) )

  CR = geom0 ! Restores original geometry
  CALL GETFORCE
  call GETDIPOLEVEC( dipoleVec )
  
  allocate( eValsHessian(3*NATS) )
  allocate( eVecsHessian(3*NATS,3*NATS) )

  if( VERBOSE >= 1 ) then
    write(*,*) ""
	write(*,*) "D^T*Hessian*D: x10^4"
	do i=1,3*NATS
		do j=1,3*NATS
		write(*,"(F10.4)",advance="no") hessian(i,j)*10000.0
		end do
		write(*,*) ""
	end do
  end if

  allocate( workSpace( 3*3*NATS-1 ) )
  eVecsHessian = hessian

  mlsi = TIME_MLS()
  
  ! Compute the eigen values and eigen vectors using the upper elements of the symmetric matrix
  call dsyev( 'V', 'L', 3*NATS, eVecsHessian, 3*NATS, eValsHessian, workSpace, 3*3*NATS-1, info )
  
  write(*,*) ""
  write(*,"(A,F10.2,A)") " Hessian Diagonalization", TIME_MLS() - MLSI, "ms"

  deallocate( workSpace )

  if ( info /= 0 ) then
    write(*,*) "### ERROR ### VIBANALYSIS.dsyev: matrix diagonalization failed"
    stop
  end if

  allocate( freqs(3*NATS) )
  allocate( intensities(3*NATS) )
  
  write(*,*) ""
  write(*,*) " Low Vib. Freqs (cm-1) "
  write(*,*) "-----------------------"
  do i=1,3*NATS-nVib
    if( eValsHessian(i) < 0.0_latteprec ) then
	  freqs(i) = -sqrt(abs(eValsHessian(i)))
      write(*,"(I10,F20.10)") i, freqs(i)/cm1 ! a.u. to cm-1
    else
      freqs(i) = sqrt(eValsHessian(i))
      write(*,"(I10,F20.10)") i, freqs(i)/cm1 ! a.u. to cm-1
    end if
  end do
  
  write(*,*) ""
  write(*,*) " Vib. Freqs (cm-1) "
  write(*,*) "-------------------"
  do i=3*NATS-nVib+1,3*NATS
    if( eValsHessian(i) < 0.0_latteprec ) then
	  freqs(i) = -sqrt(abs(eValsHessian(i)))
      write(*,"(I10,F20.10)") i, freqs(i)/cm1 ! a.u. to cm-1
    else
      freqs(i) = sqrt(eValsHessian(i))
      write(*,"(I10,F20.10)") i, freqs(i)/cm1 ! a.u. to cm-1
    end if
  end do
  write(*,*) ""
  
  ZPE = 0.0_latteprec
  do i=3*NATS-nVib+1,3*NATS
    if( freqs(i) > 0.0_latteprec ) then
      ZPE = ZPE + freqs(i)
    end if
  end do
  ZPE = ZPE/2.0_latteprec
  
  write(*,*) ""
  write(*,*) " ZPE = ", ZPE/cm1, "cm-1"  ! a.u. to cm-1
  write(*,*) " ZPE = ", ZPE/eV, "eV"  ! a.u. to eV
  write(*,*) " ZPE = ", ZPE/kcalMol, "kcal/mol"  ! a.u. to kcalMol
  write(*,*) ""
  
  ZPE = ZPE/eV ! Stored value in eV
  
  allocate( lCart(3*NATS,3*NATS) )
  lCart = 0.0_latteprec

  p1 = 1
  do atom1=1,NATS
    m1 = MASS(ELEMPOINTER(atom1))*amu ! amu to a.u.
    do i=1,3
    
		p2 = 1
		do atom2=1,NATS; do j=1,3
		
			if( p1 == p2 ) then
				lCart(p1,p2) = 1.0_latteprec/sqrt(m1)
			end if
			
			p2 = p2 + 1
		end do; end do
		
		p1 = p1 + 1
    end do
  end do
  
  if( VERBOSE >= 1 ) then
    write(*,*) ""
	write(*,*) "Mass matrix:"
	do i=1,3*NATS
		do j=1,3*NATS
		write(*,"(6F8.2)",advance="no") lCart(i,j)
		end do
		write(*,*) ""
	end do
  end if

  lCart = matmul( lCart, eVecsHessian )
  
  do i=1,3*NATS
    squareNorm = sqrt( sum(lCart(:,i)**2) )
    if( squareNorm > 1d-9 ) then
		lCart(:,i) = lCart(:,i)/squareNorm
	end if
  end do

  if( VERBOSE >= 1 ) then
	write(*,*) ""
	write(*,*) "Cartesian Displacements:"
	do i=1,3*NATS
		do j=1,3*NATS
		write(*,"(6F8.2)",advance="no") lCart(i,j)
		end do
		write(*,*) ""
	end do
  end if
  
  dMudQ = 0.0_latteprec
  p1 = 1 
  do i=1,3 ! x, y, z
    do j=1,NATS ! Normal mode
		dMudQ(i,p1) = dot_product( dMudR(i,:), eVecsHessian(:,p1) )
		p1 = p1 + 1
	end do
  end do
  
  write(*,*) ""
  write(*,*) "IR intensities:"
  write(*,*) "---------------"
  do i=1,3*NATS
	if( i > 3*NATS-nVib ) then
		intensities(i) = dot_product( dMudQ(:,i), dMudQ(:,i) )/((debye/angs)**2/amu)
	else
		intensities(i) = 0.0_latteprec
	end if
  	write(*,"(I5,F10.5)") i, intensities(i)
  end do
  write(*,*) ""
  
  call saveMolden( ATELE, CR, freqs, intensities, lCart, nVib )

  deallocate( geom0 )
  deallocate( hessian )
  deallocate( eValsHessian )
  deallocate( eVecsHessian )
  deallocate( transformationMatrix )
  deallocate( freqs )
  deallocate( lCart )
  deallocate( dMudR )
  deallocate( dMudQ )
  deallocate( dMudR_check )
  deallocate( intensities )

  RETURN

END SUBROUTINE VIBRATIONAL_ANALYSIS
