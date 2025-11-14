
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
! 13.11.2025    ggu     added burrying
!
!**************************************************************

!==============================================================
	module mod_beaching
!==============================================================

	logical, save :: bbeach = .true.
	logical, save :: bbeach_debug = .true.
	integer, save :: iudbg = 889

	integer, save, allocatable :: beach_node(:)
	real, save, allocatable :: beach_value(:,:)
	real, save, allocatable :: burry_value(:,:)
	real, save :: beach_rate = 0.3
	real, save :: burry_rate = 0.1
        double precision, save :: da_beach(5)

!==============================================================
	end module mod_beaching
!==============================================================

!**************************************************************

	subroutine beaching_init

	use basin
	use mod_beaching
	use shympi
	use mod_geom
	use mod_conz

	implicit none

	integer ie,ii,k,ia
	integer nvar,ivar,id
	logical b2d
	double precision dtime
	character*10 type

	logical is_open_boundary_node
        logical has_output_d

	integer, allocatable :: node_area_code(:)
	real, allocatable :: value(:)

	if( .not. bbeach ) return

	nvar = iconz
	b2d = .true.

	if( .not. bbeach_debug ) iudbg = 0

	if( beach_rate + burry_rate > 1. ) then
	  write(6,*) 'sum of beach_rate + burry_rate > 1'
	  write(6,*) beach_rate,burry_rate,beach_rate+burry_rate
	  stop 'error stop beaching_init: rate too high'
	end if

	allocate(beach_node(nkn))
	allocate(beach_value(nkn,nvar))
	allocate(burry_value(nkn,nvar))
	allocate(node_area_code(nkn))
	allocate(value(nkn))
	node_area_code = -1
        beach_value = 0
        burry_value = 0

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

        call init_output_d('itmcon','idtcon',da_beach)
        if( has_output_d(da_beach) ) then
          call shyfem_init_scalar_file('beach',nvar,b2d,id)
          da_beach(4) = id
          call shyfem_init_scalar_file('burry',nvar,b2d,id)
          da_beach(5) = id
        end if

	if( .not. bbeach_debug ) return

	type = 'beachi'
	dtime = 0.
	value = beach_node

	nvar = 1
	ivar = 75		!general index

	call shyfem_init_scalar_file(type,nvar,b2d,id)
	call shy_write_scalar2d(id,type,dtime,nvar,ivar,value)
	call shy_close_output_file(id)

	end

!**************************************************************

        subroutine beaching_run

        use basin
        use levels
	use mod_beaching
	use mod_conz
	use shympi

        implicit none

	logical bwrite
        integer k,ivar,id,iv,nvar,i,lmax
        real cb,db,c1
        real cu,du,cl
        double precision dtime
        double precision cbsum,cusum
	real, allocatable :: cmax1(:),cmaxl(:),bmax(:),umax(:)
	real cmax10,cmaxl0,bmax0,umax0
        character*20 aline

        logical next_output_d

	if( .not. bbeach ) return

	bwrite = iudbg > 0 .and. my_id == 0

        nvar = iconz
        cbsum = 0.
        cusum = 0.

        do k=1,nkn
	  lmax = ilhkv(k)
          do iv=1,nvar
            c1 = conzv(1,k,iv)
            cl = conzv(lmax,k,iv)
            cb = beach_value(k,iv)
            cu = burry_value(k,iv)
            db = beach_rate * c1
            du = burry_rate * cl
            if( beach_node(k) > 0 ) then
              c1 = c1 - db
              cb = cb + db
              conzv(1,k,iv) = c1
              beach_value(k,iv) = cb
	    end if
            cl = cl - du
            cu = cu + du
            conzv(lmax,k,iv) = cl
            burry_value(k,iv) = cu
            cbsum = cbsum + cb
            cusum = cusum + cu
          end do
        end do

	if( iudbg > 0 ) then
	  allocate(cmax1(nvar),cmaxl(nvar),bmax(nvar),umax(nvar))
	  do iv=1,nvar
	    cmax1(iv) = maxval(conzv(1,:,iv))
	    cmaxl(iv) = maxval(conzv(lmax,:,iv))
	    bmax(iv) = maxval(beach_value(:,iv))
	    umax(iv) = maxval(burry_value(:,iv))
	    cmax1(iv) = shympi_max(cmax1(iv))	!FIXME - too slow
	    cmaxl(iv) = shympi_max(cmaxl(iv))
	    bmax(iv) = shympi_max(bmax(iv))
	    umax(iv) = shympi_max(umax(iv))
	  end do
	end if

	if( bwrite ) then
          call get_act_dtime(dtime)
          call get_act_timeline(aline)

          write(iudbg,*) 'bsum: ',aline,cbsum,cusum
	  do iv=1,nvar
	   write(iudbg,*) iv,bmax(iv),umax(iv)
	  end do
	  do iv=1,nvar
	   write(iudbg,*) iv,cmax1(iv),cmaxl(iv)
	  end do
	  do k=1,0,nkn/10
	    lmax = ilhkv(k)
	    write(iudbg,*) conzv(1,k,4),conzv(lmax,k,4),conzv(1,k,1)
	  end do
	  flush(iudbg)
	end if

        if( next_output_d(da_beach) ) then
          id = nint(da_beach(4))
          do iv=1,nvar
            ivar = 300 + iv
            call shy_write_scalar_record(id,dtime,ivar,1,beach_value(:,iv))
          end do
          id = nint(da_beach(5))
          do iv=1,nvar
            ivar = 300 + iv
            call shy_write_scalar_record(id,dtime,ivar,1,burry_value(:,iv))
          end do
        end if

        end

!**************************************************************

