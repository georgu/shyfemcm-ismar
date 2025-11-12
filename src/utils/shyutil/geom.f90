
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
! 24.07.2015	ggu	changed VERS_7_1_82
! 16.12.2015	ggu	changed VERS_7_3_16
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 13.04.2022	ggu	new array iboundv (indicator of boundary node)
! 26.04.2022	ggu	dxv, dyv eliminated
! 07.11.2024	ggu	introduced buselink (must be false)

!**************************************************************************

!==================================================================
	module mod_geom
!==================================================================

	implicit none

	integer, private, save  :: nkn_geom = 0
	integer, private, save  :: nel_geom = 0
	integer, private, save  :: ngr_geom = 0
	integer, private, save  :: nlk_geom = 0

	integer, save :: maxlnk = 0
	logical, parameter :: buselink = .false.

!	ilinkv(1) == 0
!	ilinkv(k+1) == max entries in lenkv and linkv
!	ibase = ilinkv(k)
!	n = ilinkv(k+1)
!	n - ibase -> nodes around k
!	if( lenkv(n) == 0 ) boundary node

	integer, allocatable, save :: ilinkv(:)  !pointer into arrays below
	integer, allocatable, save :: lenkv(:)   !element numbers
	integer, allocatable, save :: lenkiiv(:) !vertex number of node in elem
	integer, allocatable, save :: linkv(:)   !node numbers

	integer, allocatable, save :: iboundv(:) !1: material bound
	integer, allocatable, save :: ieltv(:,:)
	integer, allocatable, save :: kantv(:,:)

	integer, allocatable, save :: nlist(:,:)
	integer, allocatable, save :: elist(:,:)

!==================================================================
	contains
!==================================================================

	subroutine mod_geom_init(nkn,nel,ngr)

	integer nkn,nel,ngr

	integer nlk

        if( ngr == ngr_geom .and. nel == nel_geom .and. nkn == nkn_geom ) return

        if( ngr > 0 .or. nel > 0 .or. nkn > 0 ) then
          if( ngr == 0 .or. nel == 0 .or. nkn == 0 ) then
            write(6,*) 'ngr,nel,nkn: ',ngr,nel,nkn
            stop 'error stop mod_geom_init: incompatible parameters'
          end if
        end if

        if( nkn_geom > 0 ) then
	  if( buselink ) then
            deallocate(ilinkv)
            deallocate(lenkv)
            deallocate(lenkiiv)
            deallocate(linkv)
	  end if
          deallocate(iboundv)
          deallocate(ieltv)
          deallocate(kantv)
          deallocate(nlist)
          deallocate(elist)
        end if

	nlk = 3*nel + 2*nkn
	maxlnk = ngr

        ngr_geom = ngr
        nel_geom = nel
        nkn_geom = nkn
        nlk_geom = nlk

	if( nkn == 0 ) return

	if( buselink ) then
          allocate(ilinkv(nkn+1))
          allocate(lenkv(nlk))
          allocate(lenkiiv(nlk))
          allocate(linkv(nlk))
	end if
        allocate(iboundv(nkn))
        allocate(ieltv(3,nel))
        allocate(kantv(2,nkn))
        allocate(nlist(0:ngr,nkn))
        allocate(elist(0:ngr,nkn))

	end subroutine mod_geom_init

!------------------------------------------------------------------

	subroutine mod_geom_release

	call mod_geom_init(0,0,0)

	end subroutine mod_geom_release

!------------------------------------------------------------------

	pure function is_material_boundary_node(k)

	logical is_material_boundary_node
	integer, intent(in) :: k

	is_material_boundary_node = ( iboundv(k) /= 0 )

	end function is_material_boundary_node

!------------------------------------------------------------------

	pure function n_nodes_around(k)

	integer n_nodes_around
	integer, intent(in) :: k

	n_nodes_around = nlist(0,k)

	end function n_nodes_around

!------------------------------------------------------------------

	pure function n_elems_around(k)

	integer n_elems_around
	integer, intent(in) :: k

	n_elems_around = elist(0,k)

	end function n_elems_around

!==================================================================
	end module mod_geom
!==================================================================

