!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 2018.  Los Alamos National Security, LLC. This material was    !
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                   !!
!!  Several functions were taken from the String class of the SciFT project          !!
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

!>
!! @brief
!!
module skfParser_
	implicit none
	private
	
	character(3), public, parameter :: elementsSimpleFormat(10) &
		= ["dd0", "dd1", "dd2", "pd0", "pd1", "pp0", "pp1", "sd0", "sp0", "ss0"]
		
	character(3), public, parameter :: elementsExtendedFormat(20) &
		= [ "ff0", "ff1", "ff2", "ff3", "df0", "df1", "df2", "dd0", "dd1", "dd2", &
		    "pf0", "pf1", "pd0", "pd1", "pp0", "pp1", "sf0", "sd0", "sp0", "ss0"]
	
	public :: &
		FString_count, &
		FString_replace, &
		FString_removeTabs, &
		FString_split, &
		FString_toLogical, &
		FString_toInteger, &
		FString_toReal, &
		FString_toComplex, &
		FString_toIntegerArray, &
		FString_toRealArray, &
		FString_fromInteger, &
		FString_fromReal, &
		FString_fromLogical, &
		FString_fromIntegerArray, &
		FString_fromRealArray, &
		skfParser_test
	
	type, public :: skfParser
	
		character(1000) :: iFileName
		logical :: extendedFormat
		logical :: homoNuclearCase
		logical :: spline !< .true. if the repulsive interaction is described via splines
		
		real(8) :: gridDist !< Distance between the grid points of the integral table
		integer :: nGridPoints !< Number of points in the table
		integer :: nIntegrals !< Number of integrals. 20 for "simple format"
		real(8), allocatable :: E(:) !< On-site energies for the angular momenta (d,p,s) "simple format" or (f,d,p,s) "extended format"
		real(8), allocatable :: U(:) !< Hubbard U values for the appropriate angular momenta (d,p,s) "simple format" or (f,d,p,s) "extended format"
		real(8), allocatable :: f(:) !< Occupations (for the neutral atom) in the ground state (d,p,s) "simple format" or (f,d,p,s) "extended format"
		real(8) :: mass !< Mass of the given atom in atomic mass units
		real(8) :: c(2:9), rcut !< Polynomial coefficients and the cutoff radius of the repulsive interaction (zero, if the repulsive is described by splines)
		real(8), allocatable :: H(:,:) !< Integral table containing the DFTB Hamiltonian elements
		real(8), allocatable :: S(:,:) !< Integral table containing the overlap matrix elements
		
		integer :: nInt !< Number of (subsequent) intervals being described by various cubic splines
		real(8) :: cutoff !< The cutoff of the repulsive interaction
		real(8) :: a(1:3) !< Coefficients to describe the exponential part of the repulsive 
		real(8), allocatable :: start(:), end(:) !< start (r0) and end describing the bounds of the distance range
		real(8), allocatable :: c0(:), c1(:), c2(:), c3(:), c4(:), c5(:) !< Coefficients of the spline. c4 and c5 are zero for all, except for last point
		
		character(3) :: symbols(2) ! Loaded from the name of the file
		
		! Non-standard parameters in the sdk format. They are loaded from the skmf metafile
		real(8), allocatable :: W(:) !< Magnetic Hubbard coefficients (d,p,s) "simple format" or (f,d,p,s) "extended format"
		
		contains
			generic :: init => initSkfParser
			generic :: assignment(=) => copySkfParser
			
			procedure :: initSkfParser
			procedure :: copySkfParser
			final :: destroySkfParser
			procedure :: str
			procedure :: show
			
			procedure :: getR
			procedure :: getH
			procedure :: getS
			procedure :: getV
			
			procedure, private :: parseFile
			procedure, private :: loadMetafile
	end type skfParser
	
	contains
	
	!>
	!! @brief Constructor
	!!
	subroutine initSkfParser( this, iFileName )
		class(skfParser) :: this
		character(*), intent(in) :: iFileName
		
		call this%parseFile( iFileName )
	end subroutine initSkfParser
	
	!>
	!! @brief Copy constructor
	!!
	subroutine copySkfParser( this, other )
		class(skfParser), intent(inout) :: this
		class(skfParser), intent(in) :: other
		
		this%iFileName = other%iFileName
	end subroutine copySkfParser
	
	!>
	!! @brief Destructor
	!!
	subroutine destroySkfParser( this )
		type(skfParser) :: this
		
		if( allocated(this%E) ) deallocate(this%E)
		if( allocated(this%U) ) deallocate(this%U)
		if( allocated(this%f) ) deallocate(this%f)
		if( allocated(this%H) ) deallocate(this%H)
		if( allocated(this%S) ) deallocate(this%S)
		
		if( allocated(this%start) ) deallocate(this%start)
		if( allocated(this%end) ) deallocate(this%end)
		if( allocated(this%c0) ) deallocate(this%c0)
		if( allocated(this%c1) ) deallocate(this%c1)
		if( allocated(this%c2) ) deallocate(this%c2)
		if( allocated(this%c3) ) deallocate(this%c3)
		if( allocated(this%c4) ) deallocate(this%c4)
		if( allocated(this%c5) ) deallocate(this%c5)
	end subroutine destroySkfParser
	
	!>
	!! @brief Convert to string
	!!
	function str( this ) result( output )
		class(skfParser) :: this 
		character(:), allocatable :: output
		
		integer :: i
		
		output = ""
		output = trim(output)//"file="//trim(this%iFileName)//new_line('')
		output = trim(output)//"symbols="//trim(this%symbols(1))//","//trim(this%symbols(2))//new_line('')
		output = trim(output)//"extendedFormat="//trim(FString_fromLogical(this%extendedFormat))//new_line('')
		output = trim(output)//"homoNuclearCase="//trim(FString_fromLogical(this%homoNuclearCase))//new_line('')
		output = trim(output)//"gridDist="//trim(adjustl(FString_fromReal(this%gridDist,"(F10.3)")))//new_line('')
		output = trim(output)//"nGridPoints="//trim(FString_fromInteger(this%nGridPoints))//new_line('')
		output = trim(output)//"nIntegrals="//trim(FString_fromInteger(this%nIntegrals))//new_line('')
		
		if( this%homoNuclearCase ) then
			output = trim(output)//"E="//trim(FString_fromRealArray(this%E,"("//FString_fromInteger(size(this%E))//"F10.5)"))//new_line('')
			output = trim(output)//"U="//trim(FString_fromRealArray(this%U,"("//FString_fromInteger(size(this%U))//"F10.5)"))//new_line('')
			output = trim(output)//"f="//trim(FString_fromRealArray(this%f,"("//FString_fromInteger(size(this%f))//"F10.5)"))//new_line('')
			output = trim(output)//"W="//trim(FString_fromRealArray(this%W,"("//FString_fromInteger(size(this%W))//"F10.5)"))//new_line('')
		end if
		
		output = trim(output)//"mass="//trim(FString_fromReal(this%mass,"(F10.5)"))//new_line('')
		output = trim(output)//"c="//trim(FString_fromRealArray(this%c,"("//FString_fromInteger(size(this%c))//"F10.5)"))//new_line('')
		output = trim(output)//"rcut="//trim(FString_fromReal(this%rcut,"(F10.5)"))//new_line('')
		
		do i=1,this%nIntegrals
			if( this%extendedFormat ) then
				output = trim(output)//"H"//trim(elementsExtendedFormat(i))//"="
			else
				output = trim(output)//"H"//trim(elementsSimpleFormat(i))//"="
			end if
			
			output = trim(output)// &
			trim(FString_fromReal(this%H(i,1),"(F12.5)"))// &
			trim(FString_fromReal(this%H(i,2),"(F12.5)"))// &
			trim(FString_fromReal(this%H(i,3),"(F12.5)"))// &
			trim(FString_fromReal(this%H(i,4),"(F12.5)"))// &
			trim(FString_fromReal(this%H(i,5),"(F12.5)"))// &
			"   ..."// &
			trim(FString_fromReal(this%H(i,this%nGridPoints-4),"(F12.5)"))// &
			trim(FString_fromReal(this%H(i,this%nGridPoints-3),"(F12.5)"))// &
			trim(FString_fromReal(this%H(i,this%nGridPoints-2),"(F12.5)"))// &
			trim(FString_fromReal(this%H(i,this%nGridPoints-1),"(F12.5)"))// &
			trim(FString_fromReal(this%H(i,this%nGridPoints),"(F12.5)"))// &
			new_line('')
		end do
		
		do i=1,this%nIntegrals
			if( this%extendedFormat ) then
				output = trim(output)//"S"//trim(elementsExtendedFormat(i))//"="
			else
				output = trim(output)//"S"//trim(elementsSimpleFormat(i))//"="
			end if
			
			output = trim(output)// &
			trim(FString_fromReal(this%S(i,1),"(F12.5)"))// &
			trim(FString_fromReal(this%S(i,2),"(F12.5)"))// &
			trim(FString_fromReal(this%S(i,3),"(F12.5)"))// &
			trim(FString_fromReal(this%S(i,4),"(F12.5)"))// &
			trim(FString_fromReal(this%S(i,5),"(F12.5)"))// &
			"   ..."// &
			trim(FString_fromReal(this%S(i,this%nGridPoints-4),"(F12.5)"))// &
			trim(FString_fromReal(this%S(i,this%nGridPoints-3),"(F12.5)"))// &
			trim(FString_fromReal(this%S(i,this%nGridPoints-2),"(F12.5)"))// &
			trim(FString_fromReal(this%S(i,this%nGridPoints-1),"(F12.5)"))// &
			trim(FString_fromReal(this%S(i,this%nGridPoints),"(F12.5)"))// &
			new_line('')
		end do
		
		output = trim(output)//"spline="//trim(FString_fromLogical(this%spline))//new_line('')
		output = trim(output)//"nInt="//trim(FString_fromInteger(this%nInt))//new_line('')
		output = trim(output)//"cutoff="//trim(FString_fromReal(this%cutoff,"(F10.5)"))//new_line('')
		output = trim(output)//"a="//trim(FString_fromRealArray(this%a,"("//FString_fromInteger(size(this%a))//"F10.5)"))//new_line('')
		
		output = trim(output)//"start="
		output = trim(output)// &
		trim(FString_fromReal(this%start(1),"(F12.5)"))// &
		trim(FString_fromReal(this%start(2),"(F12.5)"))// &
		trim(FString_fromReal(this%start(3),"(F12.5)"))// &
		"   ..."// &
		trim(FString_fromReal(this%start(this%nInt-2),"(F12.5)"))// &
		trim(FString_fromReal(this%start(this%nInt-1),"(F12.5)"))// &
		trim(FString_fromReal(this%start(this%nInt),"(F12.5)"))// &
		new_line('')
		
		output = trim(output)//"end="
		output = trim(output)// &
		trim(FString_fromReal(this%end(1),"(F12.5)"))// &
		trim(FString_fromReal(this%end(2),"(F12.5)"))// &
		trim(FString_fromReal(this%end(3),"(F12.5)"))// &
		"   ..."// &
		trim(FString_fromReal(this%end(this%nInt-2),"(F12.5)"))// &
		trim(FString_fromReal(this%end(this%nInt-1),"(F12.5)"))// &
		trim(FString_fromReal(this%end(this%nInt),"(F12.5)"))// &
		new_line('')
		
		output = trim(output)//"c0="
		output = trim(output)// &
		trim(FString_fromReal(this%c0(1),"(F15.5)"))// &
		trim(FString_fromReal(this%c0(2),"(F15.5)"))// &
		trim(FString_fromReal(this%c0(3),"(F15.5)"))// &
		"   ..."// &
		trim(FString_fromReal(this%c0(this%nInt-2),"(F15.5)"))// &
		trim(FString_fromReal(this%c0(this%nInt-1),"(F15.5)"))// &
		trim(FString_fromReal(this%c0(this%nInt),"(F15.5)"))// &
		new_line('')
		
		output = trim(output)//"c1="
		output = trim(output)// &
		trim(FString_fromReal(this%c1(1),"(F15.5)"))// &
		trim(FString_fromReal(this%c1(2),"(F15.5)"))// &
		trim(FString_fromReal(this%c1(3),"(F15.5)"))// &
		"   ..."// &
		trim(FString_fromReal(this%c1(this%nInt-2),"(F15.5)"))// &
		trim(FString_fromReal(this%c1(this%nInt-1),"(F15.5)"))// &
		trim(FString_fromReal(this%c1(this%nInt),"(F15.5)"))// &
		new_line('')
		
		output = trim(output)//"c2="
		output = trim(output)// &
		trim(FString_fromReal(this%c2(1),"(F15.5)"))// &
		trim(FString_fromReal(this%c2(2),"(F15.5)"))// &
		trim(FString_fromReal(this%c2(3),"(F15.5)"))// &
		"   ..."// &
		trim(FString_fromReal(this%c2(this%nInt-2),"(F15.5)"))// &
		trim(FString_fromReal(this%c2(this%nInt-1),"(F15.5)"))// &
		trim(FString_fromReal(this%c2(this%nInt),"(F15.5)"))// &
		new_line('')
		
		output = trim(output)//"c3="
		output = trim(output)// &
		trim(FString_fromReal(this%c3(1),"(F15.5)"))// &
		trim(FString_fromReal(this%c3(2),"(F15.5)"))// &
		trim(FString_fromReal(this%c3(3),"(F15.5)"))// &
		"   ..."// &
		trim(FString_fromReal(this%c3(this%nInt-2),"(F15.5)"))// &
		trim(FString_fromReal(this%c3(this%nInt-1),"(F15.5)"))// &
		trim(FString_fromReal(this%c3(this%nInt),"(F15.5)"))// &
		new_line('')
		
		output = trim(output)//"c4="
		output = trim(output)// &
		trim(FString_fromReal(this%c4(1),"(F15.5)"))// &
		trim(FString_fromReal(this%c4(2),"(F15.5)"))// &
		trim(FString_fromReal(this%c4(3),"(F15.5)"))// &
		"   ..."// &
		trim(FString_fromReal(this%c4(this%nInt-2),"(F15.5)"))// &
		trim(FString_fromReal(this%c4(this%nInt-1),"(F15.5)"))// &
		trim(FString_fromReal(this%c4(this%nInt),"(F15.5)"))// &
		new_line('')
		
		output = trim(output)//"c5="
		output = trim(output)// &
		trim(FString_fromReal(this%c5(1),"(F15.5)"))// &
		trim(FString_fromReal(this%c5(2),"(F15.5)"))// &
		trim(FString_fromReal(this%c5(3),"(F15.5)"))// &
		"   ..."// &
		trim(FString_fromReal(this%c5(this%nInt-2),"(F15.5)"))// &
		trim(FString_fromReal(this%c5(this%nInt-1),"(F15.5)"))// &
		trim(FString_fromReal(this%c5(this%nInt),"(F15.5)"))// &
		new_line('')
		
	end function str
	
	!>
	!! @brief Show 
	!!
	subroutine show( this, unit )
		class(skfParser) :: this
		integer, optional, intent(in) :: unit
		
		integer :: effunit
		
		effunit = 6
		if( present(unit) ) effunit = unit
		
		write(effunit,"(a)") trim(this%str())
	end subroutine show
	
	!>
	!! @brief Replaces strings like "3*0.15" by "0.15 0.15 0.15"
	!!
	subroutine expandString( str, output )
		character(*), intent(in) :: str
		character(*), intent(out) :: output
		
		integer :: i, j
		character(100), allocatable :: tokens(:), tokens2(:)
		integer :: nTimes
		
		call FString_split( str, tokens, " ," )
		
		output = ""
		do i=1,size(tokens)
			call FString_split( tokens(i), tokens2, "*" )
			
			if( size(tokens2) == 2 ) then
				do j=1,FString_toInteger( tokens2(1) )
					output = trim(output)//" "//trim(tokens2(2))
				end do
			else
				output = trim(output)//" "//trim(tokens(i))
			end if
		end do
		
		deallocate( tokens )
		deallocate( tokens2 )
		
	end subroutine expandString
	
	!>
	!! @brief
	!!
	subroutine parseFile( this, iFileName )
		class(skfParser) :: this
		character(*), intent(in) :: iFileName
		
		integer :: iostat
		integer :: nLine
		character(2000) :: buffer
		character(100), allocatable :: tokens(:)
		real(8), allocatable :: rBuffer(:)
		integer :: i
		integer :: nGridPointsCounter
		logical :: splineSectionLocated
		integer :: nIntCounter
		
		call FString_split( iFileName, tokens, "/-." )
		this%symbols = [ trim(tokens(size(tokens)-2)), trim(tokens(size(tokens)-1)) ]
		
		this%iFileName = iFileName
		
		open( unit=330, file=this%iFileName, status="old", iostat=iostat )
		
		if( iostat /= 0 ) then
		
			write(*, *) "### Error ###: The file ( ", trim(this%iFileName), " ) cannot be open"
			
		else
		
			nLine = 0
			iostat = 1
			splineSectionLocated = .false.
			
			do while( iostat /= -1 )
			
				nLine = nLine + 1
				read( 330, "(a)", iostat=iostat ) buffer
				buffer = FString_removeTabs( buffer )
				
				if( len(trim(adjustl(buffer))) == 0 ) cycle
				
				if( nLine == 1 ) then
				
					call FString_split( buffer, tokens, " ," )
				
					this%extendedFormat = .false.
					if( trim(tokens(1)) == "@" ) this%extendedFormat = .true.
					
					if( this%extendedFormat ) then
						read( 330, "(a)", iostat=iostat ) buffer
						buffer = FString_removeTabs( buffer )
						call FString_split( buffer, tokens, " ," )
					end if
					
					this%nIntegrals = 10
					if( this%extendedFormat ) this%nIntegrals = 20
					
					this%gridDist = FString_toReal( tokens(1) )
					this%nGridPoints = FString_toInteger( tokens(2) )
					nGridPointsCounter = 0
					
					if( allocated(this%H) ) deallocate(this%H)
					if( allocated(this%S) ) deallocate(this%S)
					
					allocate( this%H(this%nIntegrals,this%nGridPoints) )
					allocate( this%S(this%nIntegrals,this%nGridPoints) )
					
					this%H = 0.0_8
					this%S = 0.0_8
					
				else if( nLine == 2 ) then
					
					if( index(buffer, "*") /= 0 ) call expandString( buffer, buffer )
					call FString_split( buffer, tokens, " ," )
					
					if( allocated(this%E) ) deallocate( this%E )
					if( allocated(this%U) ) deallocate( this%U )
					if( allocated(this%f) ) deallocate( this%f )
					if( allocated(this%W) ) deallocate( this%W )
					
					this%homoNuclearCase = .false.
					if( size(tokens) == 10 .or. size(tokens) == 13 ) then
						this%homoNuclearCase = .true.
					
						if( this%extendedFormat ) then
							
							allocate( this%E(4) )
							allocate( this%U(4) )
							allocate( this%f(4) )
							allocate( this%W(4) )
							
							this%E = [ FString_toReal(tokens(1)), FString_toReal(tokens(2)), FString_toReal(tokens(3)), FString_toReal(tokens(4)) ]
							this%U = [ FString_toReal(tokens(6)), FString_toReal(tokens(7)), FString_toReal(tokens(8)), FString_toReal(tokens(9)) ]
							this%f = [ FString_toReal(tokens(10)), FString_toReal(tokens(11)), FString_toReal(tokens(12)), FString_toReal(tokens(13)) ]
							this%W = 0.0_8 ! This is read from the metadata file
							
						else
						
							allocate( this%E(3) )
							allocate( this%U(3) )
							allocate( this%f(3) )
							allocate( this%W(3) )
							
							this%E = [ FString_toReal(tokens(1)), FString_toReal(tokens(2)), FString_toReal(tokens(3)) ]
							this%U = [ FString_toReal(tokens(5)), FString_toReal(tokens(6)), FString_toReal(tokens(7)) ]
							this%f = [ FString_toReal(tokens(8)), FString_toReal(tokens(9)), FString_toReal(tokens(10)) ]
							this%W = 0.0_8 ! This is read from the metadata file
							
						end if
						
						call this%loadMetafile()
						
						nLine = nLine + 1
						read( 330, "(a)", iostat=iostat ) buffer
						buffer = FString_removeTabs( buffer )
						
						if( len(trim(adjustl(buffer))) == 0 ) cycle
						
						if( index(buffer, "*") /= 0 ) call expandString( buffer, buffer )
						call FString_split( buffer, tokens, " ," )
					end if
					
					this%mass = 0.0_8
					if( this%homoNuclearCase ) this%mass = FString_toReal(tokens(1))
					
					this%c(2:9) = 0.0_8
					if( this%homoNuclearCase ) this%c(2:9) = [ FString_toReal(tokens(2)), FString_toReal(tokens(3)), FString_toReal(tokens(4)), &
										   FString_toReal(tokens(5)), FString_toReal(tokens(6)), FString_toReal(tokens(7)), &
										   FString_toReal(tokens(8)), FString_toReal(tokens(9)) ]
															   
					if( abs(sum(this%c)) < 1d-10 ) this%spline = .true.
					
					this%rcut = 0.0_8
					if( this%homoNuclearCase ) this%rcut = FString_toReal(tokens(9))
					
				else if( nLine > 3 .and. nLine <= 3+this%nGridPoints .and. trim(adjustl(buffer)) /= "Spline" ) then
					
					nGridPointsCounter = nGridPointsCounter + 1
					
					if( index(buffer, "*") /= 0 ) call expandString( buffer, buffer )
					call FString_split( buffer, tokens, " ," )
					
					do i=1,this%nIntegrals
						this%H(i,nGridPointsCounter) = FString_toReal(tokens(i))
						this%S(i,nGridPointsCounter) = FString_toReal(tokens(this%nIntegrals+i))
					end do
					
				else if( trim(adjustl(buffer)) == "Spline" ) then
				
					splineSectionLocated = .true.
					
					read( 330, "(a)", iostat=iostat ) buffer
					buffer = FString_removeTabs( buffer )
					call FString_split( buffer, tokens, " ," )
					
					this%nInt = FString_toInteger(tokens(1))
					this%cutoff = FString_toReal(tokens(2))
					nIntCounter = 0
					
					if( allocated(this%start) ) deallocate( this%start )
					if( allocated(this%end) ) deallocate( this%end )
					if( allocated(this%c0) ) deallocate( this%c0 )
					if( allocated(this%c1) ) deallocate( this%c1 )
					if( allocated(this%c2) ) deallocate( this%c2 )
					if( allocated(this%c3) ) deallocate( this%c3 )
					if( allocated(this%c4) ) deallocate( this%c4 )
					if( allocated(this%c5) ) deallocate( this%c5 )
					
					allocate( this%start(this%nInt) )
					allocate( this%end(this%nInt) )
					allocate( this%c0(this%nInt) )
					allocate( this%c1(this%nInt) )
					allocate( this%c2(this%nInt) )
					allocate( this%c3(this%nInt) )
					allocate( this%c4(this%nInt) )
					allocate( this%c5(this%nInt) )
					
					this%start = 0.0_8
					this%end = 0.0_8
					this%c0 = 0.0_8
					this%c1 = 0.0_8
					this%c2 = 0.0_8
					this%c3 = 0.0_8
					this%c4 = 0.0_8
					this%c5 = 0.0_8
					
				else if( splineSectionLocated .and. nIntCounter == 0 ) then
					
					nIntCounter = nIntCounter + 1
					
					if( index(buffer, "*") /= 0 ) call expandString( buffer, buffer )
					call FString_split( buffer, tokens, " ," )
					
					this%a(1:3) = [ FString_toReal(tokens(1)), FString_toReal(tokens(2)), FString_toReal(tokens(3)) ]
					
				else if( splineSectionLocated .and. nIntCounter <= this%nInt ) then
					
					nIntCounter = nIntCounter + 1
					
					if( index(buffer, "*") /= 0 ) call expandString( buffer, buffer )
					call FString_split( buffer, tokens, " ," )
					
! 					write(*,"(I5,4X,A,4X,I5)") nIntCounter, trim(buffer), size(tokens)
					
					this%start(nIntCounter-1) = FString_toReal(tokens(1))
					this%end(nIntCounter-1) = FString_toReal(tokens(2))
					this%c0(nIntCounter-1) = FString_toReal(tokens(3))
					this%c1(nIntCounter-1) = FString_toReal(tokens(4))
					this%c2(nIntCounter-1) = FString_toReal(tokens(5))
					this%c3(nIntCounter-1) = FString_toReal(tokens(6))
					
					if( size(tokens) == 8 ) then
						this%c4(nIntCounter-1) = FString_toReal(tokens(7))
						this%c5(nIntCounter-1) = FString_toReal(tokens(8))
					end if
					
				end if
				
			end do
		end if
		
		if( nGridPointsCounter /= this%nGridPoints ) then
			write(*,"(A)") "@@@ WARNING @@@ Number of points declared in the first line are different"
			write(*,"(A)") "                than the actual number of points in the table (" &
							//trim(FString_fromInteger(this%nGridPoints))//"/="//trim(FString_fromInteger(nGridPointsCounter))//")" 
			write(*,*) ""
		end if
		
		close( unit=330 )
		if( allocated(tokens) ) deallocate( tokens )
		
	end subroutine parseFile
	
	!>
	!! @brief
	!!
	subroutine loadMetafile( this )
		class(skfParser) :: this
		
		character(1000) :: iFileName
		integer :: nLine
		integer :: iostat
		character(2000) :: buffer
		character(100), allocatable :: tokens(:)
		real(8) :: Ws, Wp, Wd, Wf
		
		iFileName = FString_replace( this%iFileName, ".skf", ".skmf" )
		open( unit=331, file=iFileName, status="old", iostat=iostat )
		
		if( iostat /= 0 ) then
			write(*, *) "### Error ###: The file ( ", trim(iFileName), " ) cannot be open"
		else
			this%W = 0.0_8
			
			nLine = 0
			iostat = 1
			do while( iostat /= -1 )
				nLine = nLine + 1
			
				read( 331, "(a)", iostat=iostat ) buffer
				buffer = FString_removeTabs( buffer )
				
				if( len(trim(adjustl(buffer))) == 0 ) cycle
				
				call FString_split( buffer, tokens, " " )
				
				if( .not. this%extendedFormat ) then
					Ws = 0.0_8
					if( trim(tokens(1)) == "MAGNETIC_HUBBARD_S" ) Ws = FString_toReal(tokens(2))
					
					Wp = 0.0_8
					if( trim(tokens(1)) == "MAGNETIC_HUBBARD_P" ) Wp = FString_toReal(tokens(2))
					
					Wd = 0.0_8
					if( trim(tokens(1)) == "MAGNETIC_HUBBARD_D" ) Wd = FString_toReal(tokens(2))
					
					this%W = [ Wd, Wp, Ws ]
				else
				
					Ws = 0.0_8
					if( trim(tokens(1)) == "MAGNETIC_HUBBARD_S" ) Ws = FString_toReal(tokens(2))
					
					Wp = 0.0_8
					if( trim(tokens(1)) == "MAGNETIC_HUBBARD_P" ) Wp = FString_toReal(tokens(2))
					
					Wd = 0.0_8
					if( trim(tokens(1)) == "MAGNETIC_HUBBARD_D" ) Wd = FString_toReal(tokens(2))
					
					Wf = 0.0_8
					if( trim(tokens(1)) == "MAGNETIC_HUBBARD_F" ) Wf = FString_toReal(tokens(2))
					
					this%W = [ Wf, Wd, Wp, Ws ]
				end if
			end do
		
		end if
		
		if( allocated(tokens) ) deallocate( tokens )
		
		close( unit=331 )
	end subroutine loadMetafile
	
	!>
	!! @brief Returns the radial grid
	!!
	subroutine getR( this, array )
		class(skfParser), intent(in) :: this
		real(8), allocatable :: array(:)
		
		integer :: i
		
		if( allocated(array) ) deallocate(array)
		allocate( array(this%nGridPoints) )
		
		do i=1,this%nGridPoints
			array(i) = 0.0 + real(i-1,8)*this%gridDist
		end do
	end subroutine getR
	
	!>
	!! @brief
	!!
	subroutine getH( this, i, array )
		class(skfParser), intent(in) :: this
		integer, intent(in) :: i
		real(8), allocatable :: array(:)
		
		if( allocated(array) ) deallocate(array)
		allocate( array(this%nGridPoints) )
		
		array = this%H( i, : )
	end subroutine getH
	
	!>
	!! @brief
	!!
	subroutine getS( this, i, array )
		class(skfParser), intent(in) :: this
		integer, intent(in) :: i
		real(8), allocatable :: array(:)
		
		if( allocated(array) ) deallocate(array)
		allocate( array(this%nGridPoints) )
		
		array = this%S( i, : )
	end subroutine getS
	
	!>
	!! @brief
	!!
	subroutine getV( this, i, array )
		class(skfParser), intent(in) :: this
		integer, intent(in) :: i
		real(8), allocatable :: array(:)
		
		if( allocated(array) ) deallocate(array)
		allocate( array(this%nGridPoints) )
	end subroutine getV
	
	!>
	!! @brief 
	!!
	function FString_count( str, ref, matchPos, wholeWords, wordSeparators ) result( nMatches )
		character(*), intent(in) :: str
		character(*), intent(in) :: ref
		integer, allocatable, optional, intent(out) :: matchPos(:)
		logical, optional, intent(in) :: wholeWords
		character(*), optional, intent(in) :: wordSeparators
		integer :: nMatches
		
		logical :: effWholeWords
		character(:), allocatable :: effWordSeparators
		
		integer :: pos
		integer, allocatable :: tmpMatchPos(:)
		integer, allocatable :: tmpMatchWholeWord(:) ! 0 or 1
		character(:), allocatable :: strBuffer
		integer :: i, n
		
		effWholeWords = .false.
		if( present(wholeWords) ) effWholeWords = wholeWords
		
		! wordSeparators": "./\\()\"'-:,.;<>~!@#$%^&*|+=[]{}`~?\\.",
		effWordSeparators = ":;@-.,/_~?&=%+#*()[]{} "//achar(9)
		if( present(wordSeparators) ) effWordSeparators = wordSeparators
		
		strBuffer = str
		
		! En el peor de los casos todos los caracteres son ref
		allocate( tmpMatchPos(len(str)) )
		allocate( tmpMatchWholeWord(len(str)) )
		
		tmpMatchPos = 0
		tmpMatchWholeWord = 0
		
		n = 1
		do while( .true. )
			pos = index( strBuffer, ref )
			
			if( pos == 0 ) exit
			
			if( ( &
				  ( index( effWordSeparators, strBuffer(pos-1:pos-1) ) /= 0 .or. pos-1 == 0 ).and. &
				  index( effWordSeparators, strBuffer(pos+len(ref):pos+len(ref)) ) /= 0 &
			) ) then
				tmpMatchWholeWord(n) = 1
			end if
			
			tmpMatchPos(n) = pos + merge( 0, tmpMatchPos(n-1), n<1 )
			n = n + 1
			
			strBuffer = strBuffer(pos+1:)
		end do
		
		nMatches = n-1
		
		if( present(matchPos) ) then
			if( effWholeWords ) then
				if( allocated(matchPos) ) deallocate( matchPos )
				allocate( matchPos(sum(tmpMatchWholeWord)) )
				
				i = 1
				do n=1,nMatches
					if( tmpMatchWholeWord(n) == 1 ) then
						matchPos(i) = tmpMatchPos(n)
						i = i+1
					end if
				end do
				
				nMatches = sum(tmpMatchWholeWord)
			else
				if( allocated(matchPos) ) deallocate( matchPos )
				allocate( matchPos(nMatches) )
				
				matchPos = tmpMatchPos(1:nMatches)
			end if
		end if
		
		deallocate(tmpMatchPos)
	end function FString_count
	
	!>
	!! @brief 
	!!
	function FString_replace( str, before, after, wholeWords, wordSeparators ) result( output )
		character(*), intent(in) :: str
		character(*), intent(in) :: before
		character(*), intent(in) :: after
		logical, optional, intent(in) :: wholeWords
		character(*), optional, intent(in) :: wordSeparators
		character(:), allocatable :: output
		
		integer :: i
		integer :: nMatches
		integer, allocatable :: matchPos(:)
		
		nMatches = FString_count( str, before, matchPos, wholeWords, wordSeparators )
		
		if( nMatches == 0 ) then
			output = str
			return
		end if
		
		output = ""
		do i=1,nMatches+1
			if( i==1 ) then
				output = str(1:matchPos(i)-1)//after
			else if ( i==size(matchPos)+1 ) then
				output = output//str(min(matchPos(i-1)+len(before),len(str)+1):len(str))
			else
				output = output//str(matchPos(i-1)+len(before):matchPos(i)-1)//after
			end if
		end do
		
		deallocate( matchPos )
	end function FString_replace
	
	!>
	!! @brief Replace every occurrence of the character \tab in this string
	!!        for four blank spaces
	!! @todo Hay que ofrecer la opción de seleccionar el tamaño del tab
	!!       por ejemplo subroutine FString_removeTabs( str, tabSize )
	!! @todo Creo que retornar un allocatable es peligroso
	!!
	function FString_removeTabs( str ) result( output )
		character(*), intent(inout) :: str
		character(:), allocatable :: output
		
		output = FString_replace( str, achar(9), "    " )
	end function FString_removeTabs
	
	!>
	!! @brief 
	!!
	subroutine FString_split( str, tokens, delimiters )
		character(*), intent(in) :: str
		character(*), allocatable, intent(out) :: tokens(:)
		character(*), intent(in) :: delimiters
		
		integer, allocatable :: posTokenBegin(:), posTokenEnd(:)
		logical :: isDelimiter, isDelimiterPrev
		integer :: ntokens
		integer :: i, j
		
		if( allocated(tokens) ) then
			deallocate(tokens)
		end if
		
		if( len(trim(str)) == 0 ) then
			allocate( tokens(1) )
			tokens(1) = str
			return
		end if
		
		! En el peor de los casos todos los caracteres son separadores
		allocate( posTokenBegin(len(str)) )
		allocate( posTokenEnd(len(str)) )
		
		posTokenBegin = 0
		posTokenEnd = 0
		
		!! scan
		ntokens = 1
		isDelimiterPrev = .true.
		do i=1,len(str)
! 			write(*,"(3A)",advance="no") str(i:i)
			
			isDelimiter = .false.
			do j=1,len(delimiters)
				if( str(i:i) == delimiters(j:j) ) then
					isDelimiter = .true.
					exit
				end if
			end do
			
			if( isDelimiter .and. .not. isDelimiterPrev ) then
				posTokenEnd( ntokens ) = i-1
				ntokens = ntokens + 1
! 				write(*,*) "    E"
			else if( .not. isDelimiter .and. isDelimiterPrev ) then
				posTokenBegin( ntokens ) = i
! 				write(*,*) "    B"
! 			else
! 				write(*,*) ""
			end if
			
			isDelimiterPrev = isDelimiter
		end do
		
		if( posTokenEnd(ntokens) == 0 .and. .not. isDelimiter ) then
			posTokenEnd( ntokens ) = len( str )
		else
			ntokens = ntokens - 1 
		end if
		
! 		write(*,"(A,<len(str)>I3)") "ntokens = ", ntokens
! 		write(*,"(A,<len(str)>I3)") "  begin = ", posTokenBegin
! 		write(*,"(A,<len(str)>I3)") "    end = ", posTokenEnd
		
		allocate( tokens(ntokens) )
		
		do i=1,ntokens
			tokens(i) = str( posTokenBegin(i):posTokenEnd(i) )
		end do
		
		deallocate( posTokenBegin )
		deallocate( posTokenEnd )
		
	end subroutine FString_split
	
	!>
	!! @brief 
	!!
	function FString_toLogical( str ) result( output )
		character(*), intent(in) :: str
		logical :: output
		
		read( str, * ) output
	end function FString_toLogical
	
	!>
	!! @brief 
	!!
	function FString_toInteger( str ) result( output )
		character(*), intent(in) :: str
		integer :: output
		
		read( str, * ) output
	end function FString_toInteger
	
	!>
	!! @brief 
	!!
	function FString_toReal( str ) result( output )
		character(*), intent(in) :: str
		real(8) :: output
		
		read( str, * ) output
	end function FString_toReal
	
	!>
	!! @brief 
	!!
	function FString_toComplex( str ) result( output )
		character(*), intent(in) :: str
		complex(8) :: output
		
		read( str, * ) output
	end function FString_toComplex
	
	!>
	!! @brief 
	!!
	subroutine FString_toIntegerArray( str, output )
		character(*), intent(in) :: str
		integer, allocatable :: output(:)
		
		character(1000), allocatable :: tokens(:)
		character(1000) :: strBuffer
		integer :: i
		
		call FString_split( trim(adjustl(str)), tokens, "()[]" )
		strBuffer = tokens(1)
		call FString_split( strBuffer, tokens, "," )
		
		if(  allocated(output) ) deallocate( output )
		allocate( output(size(tokens)) )
		
		do i=1,size(tokens)
			read( tokens(i), * ) output(i)
		end do
		
		deallocate(tokens)
	end subroutine FString_toIntegerArray
	
	!>
	!! @brief 
	!!
	subroutine FString_toRealArray( str, output )
		character(*), intent(in) :: str
		real(8), allocatable :: output(:)
		
		character(100000), allocatable :: tokens(:)
		character(100000) :: strBuffer
		integer :: i
		
		call FString_split( trim(adjustl(str)), tokens, "()[]" )
		strBuffer = tokens(1)
		
		if( len_trim(strBuffer) == 0 ) then
			deallocate(tokens)
! 			deallocate(strBuffer)
			return
		end if
		
		call FString_split( strBuffer, tokens, "," )
		
		if(  allocated(output) ) deallocate( output )
		allocate( output(size(tokens)) )
		
		do i=1,size(tokens)
			read( tokens(i), * ) output(i)
		end do
		
		deallocate(tokens)
! 		deallocate(strBuffer)
	end subroutine FString_toRealArray
	
	!>
	!! @brief 
	!!
	function FString_fromInteger( val, format ) result( output )
		integer, intent(in) :: val
		character(*), optional, intent(in) :: format
		character(1000) :: output
		
		character(1000) :: strBuffer
		
		if( present(format) ) then
			write( strBuffer, format ) val
			output = strBuffer
		else
			write( strBuffer, * ) val
			output = trim(adjustl(strBuffer))
		end if
	end function FString_fromInteger
	
	!>
	!! @brief 
	!!
	function FString_fromReal( val, format ) result( output )
		real(8), intent(in) :: val
		character(*), optional, intent(in) :: format
		character(:), allocatable :: output
		
		character(1000) :: strBuffer
		
		if( present(format) ) then
			write( strBuffer, format ) val
			output = strBuffer
		else
			write( strBuffer, * ) val
			output = trim(adjustl(strBuffer))
		end if
	end function FString_fromReal
	
	!>
	!! @brief 
	!!
	function FString_fromLogical( val, format ) result( output )
		logical, intent(in) :: val
		character(*), optional, intent(in) :: format
		character(:), allocatable :: output
		
		character(1000) :: strBuffer
		
		if( present(format) ) then
			write( strBuffer, format ) val
			output = strBuffer
		else
			write( strBuffer, * ) val
			output = trim(adjustl(strBuffer))
		end if
	end function FString_fromLogical
	
	!>
	!! @brief 
	!!
	function FString_fromIntegerArray( val, format ) result( output )
		integer, intent(in) :: val(:)
		character(*), optional, intent(in) :: format
		character(1000) :: output
		
		character(1000) :: strBuffer
		
		if( present(format) ) then
			write( strBuffer, format ) val
		else
			write( strBuffer, * ) val
		end if
		
		output = "( "//trim(adjustl(strBuffer))//" )"
	end function FString_fromIntegerArray
	
	!>
	!! @brief 
	!!
	function FString_fromRealArray( val, format ) result( output )
		real(8), intent(in) :: val(:)
		character(*), optional, intent(in) :: format
		character(:), allocatable :: output
		
		character(1000) :: strBuffer
		
		if( present(format) ) then
			write( strBuffer, format ) val
		else
			write( strBuffer, * ) val
		end if
		
		output = "( "//trim(adjustl(strBuffer))//" )"
	end function FString_fromRealArray
	
	!>
	!! @brief Test method
	!!
	subroutine skfParser_test()
		type(skfParser) :: parser
		
		write(*,*) ""
		write(*,*) "Testing for: extendedFormat=F, homoNuclearCase=T"
		write(*,*) "------------------------------------------------"
		call parser%init( "mio-1-1/C-C.skf" )  ! simple format
		call parser%show()
		write(*,*) ""
		write(*,*) "Testing for: extendedFormat=F, homoNuclearCase=F"
		write(*,*) "------------------------------------------------"
		call parser%init( "mio-1-1/C-H.skf" )  ! simple format
		call parser%show()
! 		write(*,*) ""
! 		write(*,*) "Testing for: extendedFormat=F, homoNuclearCase=T"
! 		write(*,*) "------------------------------------------------"
! 		call parser%init( "tiorg-0-1/Ti-Ti.skf" )  ! simple format
! 		call parser%show()
! 		write(*,*) ""
! 		write(*,*) "Testing for: extendedFormat=F, homoNuclearCase=F"
! 		write(*,*) "------------------------------------------------"
! 		call parser%init( "tiorg-0-1/Ti-O.skf" )  ! simple format
! 		call parser%show()
! 		write(*,*) ""
! 		write(*,*) "Testing for: extendedFormat=T, homoNuclearCase=T"
! 		write(*,*) "------------------------------------------------"
! 		call parser%init( "rare-0-2/EuEu.skf" ) ! extended format
! 		call parser%show()
! 		write(*,*) ""
! 		write(*,*) "Testing for: extendedFormat=T, homoNuclearCase=F"
! 		write(*,*) "------------------------------------------------"
! 		call parser%init( "rare-0-2/EuN.skf" ) ! extended format
! 		call parser%show()
		
	end subroutine skfParser_test
	
end module skfParser_
