
!--------------------------------------------------------------------------
!
!    Copyright (C) 2025  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published
!    by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main
!    directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

! contains compatibility code for old compiler versions

! revision log :
!
! 28.01.2025    ggu     written from scratch
! 08.03.2025    ggu     added compiler information
! 09.03.2025    ggu     added compiler profile
! 09.03.2025    ggu     added compiler profile
! 12.03.2025    clc     fixed bug for INTEL compiler in compiler information

! notes :
!
! find list of defined macros:
!	gfortran -cpp -E -dM - < /dev/null
!	icc -E -dM - < /dev/null
!
! gfortran:
!	__GFORTRAN__		gfortran compiler
!	__VERSION__		version
!	__GNUC__		major version
!
! intel:
!	__INTEL_COMPILER	version

!**************************************************************************

!============================================================
	module mod_compatibility
!============================================================

	implicit none

	logical, parameter :: comp_debug = .false.

        INTERFACE gfindloc
        MODULE PROCEDURE         gfindloc_i1
        END INTERFACE

!============================================================
	contains
!============================================================

	function gfindloc_i1(array,val)

	integer gfindloc_i1
	integer array(:)
	integer val

	integer i,n

	!gfindloc_i1 = findloc(array,val,1)
	!return

!#if defined(__GFORTRAN__) && __GNUC__ == 4
	n = size(array)
	do i=1,n
	  if( array(i) == val ) exit
	end do
	if( i > n ) i = 0
	gfindloc_i1 = i
!#else
!	gfindloc_i1 = findloc(array,val,1)
!#endif

	end function

!************************************************************

	subroutine compiler(string)

	implicit none

	character*(*) string

	integer nvers
	character*40 version

#if defined(__GFORTRAN__)
	if( comp_debug ) write(6,*) 'gfortran compiler'
	nvers = __GNUC__
	version = __VERSION__
	string = 'gfortran (' // trim(version) // ')'
	!write(6,*) 'version = ',trim(version)
	!write(6,*) 'major version = ',nvers
#elif defined(__INTEL_COMPILER)
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ ":" TOSTRING(__LINE__)
	if( comp_debug ) write(6,*) 'intel compiler'
	nvers = __INTEL_COMPILER
	string = 'INTEL (' // TOSTRING(__INTEL_COMPILER) // ')'
	!write(6,*) 'version = ',nvers
	!write(6,*) 'major version = ',nvers/100
#else
	if( comp_debug ) write(6,*) 'compiler not recognized'
	string = 'unknown compiler'
#endif

	end

!************************************************************

	subroutine compiler_profile(string)

	implicit none

	character*(*) string

	integer nvers
	character*40 cprofile

#if defined(SHYFEM_NORMAL)
	cprofile = 'normal'
#elif defined(SHYFEM_CHECK)
	cprofile = 'check'
#elif defined(SHYFEM_SPEED)
	cprofile = 'speed'
#else
	cprofile = 'unknown'
#endif

	if( comp_debug ) write(6,*) 'compiler profile: ',trim(cprofile)
	string = cprofile

	end

!============================================================
	end module mod_compatibility
!============================================================

	subroutine sub_test_compatibility

	use mod_compatibility

	implicit none

	integer nvers
	character*40 version
	character*40 cprofile
	character*80 string

#if defined(__GFORTRAN__)
	write(6,*) 'gfortran compiler'
	nvers = __GNUC__
	version = __VERSION__
	write(6,*) 'version = ',trim(version)
	write(6,*) 'major version = ',nvers
#elif defined(__INTEL_COMPILER)
	write(6,*) 'intel compiler'
	nvers = __INTEL_COMPILER
	write(6,*) 'version = ',nvers
	write(6,*) 'major version = ',nvers/100
#else
	write(6,*) 'compiler not recognized'
#endif

#if defined(SHYFEM_NORMAL)
	cprofile = 'normal'
#elif defined(SHYFEM_CHECK)
	cprofile = 'check'
#elif defined(SHYFEM_SPEED)
	cprofile = 'speed'
#else
	cprofile = 'unknown'
#endif
	write(6,*) 'compiler profile: ',trim(cprofile)
	call exit(0)

	call compiler(string)
	write(6,*) 'compiler: ',trim(string)

	end

!************************************************************
!	program main_test_compatibility
!	call sub_test_compatibility
!	end program
!************************************************************

