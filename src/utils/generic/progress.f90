
!--------------------------------------------------------------------------
!
!    Copyright (C) 2024  Ian Bush
!    Copyright (C) 2026  Georg Umgiesser
!
!    This file is part of SHYFEM.
!
!    SHYFEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SHYFEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SHYFEM. Please see the file COPYING in the main directory.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Contributions to this file can be found below in the revision log.
!
!--------------------------------------------------------------------------

! implement a progress bar in fortran
!
! some pieces have been copied from the internt site:
! stackoverflow.com/questions/78002119/how-to-print-a-progress-bar-in-fortran
!
! revision log :
!
! 26.01.2026    ggu     implemented from the internet

!********************************************************************

!====================================================================
	module mod_progress_bar
!====================================================================

	implicit none

	integer, save :: lmax  = 0
	character*80, save :: pfull = ' '
	character*80, save :: pact

!====================================================================
	contains
!====================================================================

	subroutine progress_bar_init(string)

! initialize progress bar

	implicit none

	character*(*) string

	lmax = 80 - 11 - len(string)

	pfull = ' '
	pfull = repeat('|',lmax)

	end

!********************************************************************

	subroutine progress_bar_finalize

! finalize progress bar

	implicit none

	write(6,*)

	end

!********************************************************************

	subroutine progress_bar_print(string,progress) 

! print progress bar

	implicit none

	character*(*) string		!string before progress bar
	real, intent(in) :: progress	!progress [0-1]

	integer :: lpad
    
	lpad = int(progress * (lmax + 1))
	lpad = min(lpad,lmax)
	pact(1:lmax) = ' '
	pact(1:lpad) = pfull(1:lpad)
    
	write(6,'(a, a, f5.1, a2)', advance="no") & 
     &	        char(13), string//' ', progress*100, '% '
	write(6, "(A)", advance="no") "["
	write(6, "(A)", advance="no")  pact(1:lmax)
	write(6, "(A)", advance="no") "]"
    
	end subroutine progress_bar_print

!====================================================================
	end module mod_progress_bar
!====================================================================

!********************************************************************
!********************************************************************
!********************************************************************
! testing routines
!********************************************************************
!********************************************************************
!********************************************************************

	subroutine test_progress_bar

! testing the progress bar

	use mod_progress_bar

	implicit none

	integer, parameter :: ndim = 500
	integer i
	real progress
	double precision load
	character*40 string

	string = 'my progress '
	string(40:40) = '+'
	call progress_bar_init(string)

	do i=1,ndim
	  call run_load(300000,load)
	  progress = i / float(ndim)
	  call progress_bar_print(string,progress)
	end do

	call progress_bar_finalize

	end

!********************************************************************

	subroutine run_load(n,load)

! artificial loading of the cpu

	implicit none

	integer n
	double precision load

	integer i
	double precision fact
	
	load = 0.
	fact = 0.9

	do i=1,n
	  load = load + exp(fact)
	end do

	end

!********************************************************************
!	program main_print_progress_bar
!	call test_progress_bar
!	end
!********************************************************************

