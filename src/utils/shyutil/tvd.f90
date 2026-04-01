
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015,2019  Georg Umgiesser
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
! 10.07.2015	ggu	changed VERS_7_1_50
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 03.09.2024	ggu	new data structures
! 23.09.2024	ggu	new array itvdup()

!**************************************************************************

!==================================================================
        module mod_tvd
!==================================================================

        implicit none

        integer, private, save :: nel_tvd = 0

	integer, save :: itvd_type = 0

        real, allocatable, save :: xtvdup(:,:,:)	!x-coordinate
        real, allocatable, save :: ytvdup(:,:,:)	!y-coordinate
        integer, allocatable, save :: ietvdup(:,:,:)	!internal element number
        integer, allocatable, save :: ieetvdup(:,:,:)	!external elem number
        integer, allocatable, save :: iatvdup(:,:,:)	!domain of element
        integer, allocatable, save :: ltvdup(:,:,:)	!lmax of element
        integer, allocatable, save :: itvdup(:,:,:)	!index into array

!==================================================================
	contains
!==================================================================

        subroutine mod_tvd_init(nel)

        integer nel

        if( nel == nel_tvd ) return

        if( nel_tvd > 0 ) then
          deallocate(xtvdup)
          deallocate(ytvdup)
          deallocate(ietvdup)
          deallocate(ieetvdup)
          deallocate(iatvdup)
          deallocate(ltvdup)
          deallocate(itvdup)
        end if

        nel_tvd = nel

        if( nel == 0 ) return

        allocate(xtvdup(3,3,nel))
        allocate(ytvdup(3,3,nel))
        allocate(ietvdup(3,3,nel))
        allocate(ieetvdup(3,3,nel))
        allocate(iatvdup(3,3,nel))
        allocate(ltvdup(3,3,nel))
        allocate(itvdup(3,3,nel))

	!write(6,*) 'mod_tvd_init successfully called ',nel

        end subroutine mod_tvd_init

!==================================================================
        end module mod_tvd
!==================================================================

