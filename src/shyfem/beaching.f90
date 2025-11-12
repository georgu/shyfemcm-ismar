
!--------------------------------------------------------------------------
!
!    Copyright (C) 2003,2008-2009,2011,2015-2019  Georg Umgiesser
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

! routines for simulating beaching of tracers
!
! revision log :
!
! 12.11.2025    ggu     created from scratch
!
!**************************************************************

!==============================================================
	module mod_beaching
!==============================================================

	logical, save :: bbeach = .false.
	logical, save :: bbeach_debug = .true.
	integer, save, allocatable :: beach_node(:)

!==============================================================
	end module mod_beaching
!==============================================================

!**************************************************************

	subroutine beaching_init

	use basin
	use mod_beaching
	use shympi
	use mod_geom

	implicit none

	integer ie,ii,k,ia
	integer nvar,ivar,id
	logical b2d
	double precision dtime
	character*10 type

	logical is_open_boundary_node

	integer, allocatable :: node_area_code(:)
	real, allocatable :: value(:)

	if( .not. bbeach ) return

	allocate(beach_node(nkn))
	allocate(node_area_code(nkn))
	allocate(value(nkn))
	node_area_code = -1

	do ie=1,nel
	  ia = iarv(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    node_area_code(k) = max(node_area_code(k),ia)
	  end do
	end do

	call shympi_exchange_2d_node(node_area_code)

	beach_node = 0

	do k=1,nkn
	  if( node_area_code(k) < 2 ) then
	    if( is_material_boundary_node(k) ) then
	      if( .not. is_open_boundary_node(k) ) then
	        beach_node(k) = 1
	      end if
	    end if
	  end if
	end do
	
	call shympi_exchange_2d_node(beach_node)

	if( .not. bbeach_debug ) return

	nvar = 1
	ivar = 75		!general index
	b2d = .true.
	type = 'beachi'
	dtime = 0.
	value = beach_node

	call shyfem_init_scalar_file(type,nvar,b2d,id)
	call shy_write_scalar2d(id,type,dtime,nvar,ivar,value)
	call shy_close_output_file(id)

	end

!**************************************************************

