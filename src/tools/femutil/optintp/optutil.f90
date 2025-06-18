
!--------------------------------------------------------------------------
!
!    Copyright (C) 2003-2004,2009,2014-2015,2018-2019  Georg Umgiesser
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

! utilities for optimal interpolation
!
! revision log :
!
! 12.01.2025	ggu	written from scratch
!
!****************************************************************

	subroutine get_options(nmax,string,ns,values)

! reads maximum n values from string
! returns number of values in ns and values in values
! if error ns = -1

	implicit none

	integer, intent(in)	 :: nmax
	character*80, intent(in) :: string
	integer, intent(out)	 :: ns
	real, intent(out)	 :: values(nmax)

	integer iscanf

	ns = iscanf(string,values,0)
	if( ns < 0 ) then		!read error
	  return
	else if( ns > nmax ) then	!error - too many values
	  ns = -1
	  return
	end if

	ns = iscanf(string,values,nmax)

	end subroutine

!****************************************************************

	subroutine mima_flag(n,array,flag,amin,amax)

	implicit none

	integer n
	real array(n)
	real flag
	real amin,amax

	amin = minval(array,array/=flag)
	amax = maxval(array,array/=flag)

	end

!****************************************************************

