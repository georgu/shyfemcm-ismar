
!--------------------------------------------------------------------------
!
!    Copyright (C) 2024  Georg Umgiesser
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

! routines to easily reallocate arrays
!
! revision log :
!
! 11.12.2024    ggu     written from scratch

!===============================================================
	module mod_reallocate
!===============================================================

	implicit none

        INTERFACE mem_reallocate
        MODULE PROCEDURE   mem_reallocate_i
        END INTERFACE

        INTERFACE mem_extend
        MODULE PROCEDURE   mem_extend_i
        END INTERFACE

!===============================================================
	contains
!===============================================================

	subroutine mem_reallocate_i(nnew,array)

	integer nnew
	integer, allocatable :: array(:)

	integer nold
	integer, allocatable :: aux(:)

	nold = size(array)

	allocate(aux(nnew))
	aux(1:nold) = array(1:nold)
        call move_alloc(aux,array)

	end

!***************************************************************

	subroutine mem_extend_i(array,extend)

	integer, allocatable :: array(:)
	integer, allocatable ::  extend(:)

	integer nold,nnew,ntot
	integer, allocatable :: aux(:)

	nold = size(array)
	nnew = size(extend)
	ntot = nold + nnew

	allocate(aux(ntot))
	aux(1:nold) = array(1:nold)
        call move_alloc(aux,array)
	array(nold+1:ntot) = extend(1:nnew)

	end

!===============================================================
	end module mod_reallocate
!===============================================================

