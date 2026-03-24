
!--------------------------------------------------------------------------
!
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

! revision log :
!
! 20.03.2025	ggu	copied from shyfem_sub.f90

!**************************************************************************

!==================================================================
	module mod_shyfem
!==================================================================

	implicit none

        logical, save :: bdebug_shyfem	= .false.
        logical, save :: bdebout	= .false.
        logical, save :: bmpirun	= .false.

        logical, save :: bverbose	= .false.
        logical, save :: bquiet		= .false.
        logical, save :: bsilent	= .false.

        logical, save :: bfirst		= .true.

!==================================================================
        end module mod_shyfem
!==================================================================

