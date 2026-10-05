
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

! implements a fast search method for point in element
!
! revision log :
!
! 01.10.2024	ggu	fast_find started
! 03.10.2024	ggu	fast_find debugged
! 03.10.2024	ggu	small bug fixes

!**************************************************************************

!==================================================================
        module mod_fast_find
!==================================================================

        implicit none

        logical, private, save :: bdebug = .false.
        logical, private, save :: bplot = .false.

        integer, private, save :: nsize = 0
        integer, private, save :: idmax = 0
        integer, private, save :: nbx,nby
        real, private, save :: box_size = 0
	real, private, save :: xbmin,ybmin,xbmax,ybmax

	integer, save, allocatable :: n_elem_list(:)
	integer, save, allocatable :: elem_list(:,:)

!==================================================================
	contains
!==================================================================

        subroutine mod_fast_find_init

	write(6,*) 'mod_fast_find_init successfully called '

        end subroutine mod_fast_find_init

!******************************************************************

	subroutine fast_find_initialize(nin)

! sets up data structure

	use basin
	use stack

	integer, intent(in) :: nin

	integer nx,ny
	integer ix,iy
	integer ixmin,iymin,ixmax,iymax
	integer id,idback
	integer ie,ii,k
	integer ns,ne,nemax,netot,neaver,nezero
	real x(3),y(3)
	real xmin,ymin,xmax,ymax
	real xemin,yemin,xemax,yemax
	real dx,dy
	real dbx,dby
	real, parameter :: eps = 0.01

!	------------------------------------------------
!	find box sizes and limits
!	------------------------------------------------

	if( idmax /= 0 ) call fast_find_finalize

	if( bdebug ) write(6,*) 'starting seting up fast_find routine'

	nsize = nin
	if( nin <= 0 ) nsize = 10		!default

	call bas_get_minmax(xmin,ymin,xmax,ymax)

	dx = xmax - xmin
	dy = ymax - ymin
	xbmin = xmin - eps*dx
	xbmax = xmax + eps*dx
	ybmin = ymin - eps*dy
	ybmax = ymax + eps*dy
	dbx = dx
	dby = dy

	if( dx > dy ) then
	  box_size = (xbmax-xbmin)/nsize
	  nbx = nsize
	  nby = 1 + (ybmax-ybmin) / box_size
	  ybmax = ybmin + nby*box_size
	else
	  box_size = (ybmax-ybmin)/nsize
	  nby = nsize
	  nbx = 1 + (xbmax-xbmin) / box_size
	  xbmax = xbmin + nbx*box_size
	end if

	if( bdebug ) then
	  write(6,*) 'box: ',nsize,nbx,nby
	  write(6,*) 'box size: ',dx,dy,box_size
	  write(6,*) 'minmax: ',xbmin,xbmax,ybmin,ybmax
	end if

	if( bplot ) call fast_find_plot

!	------------------------------------------------
!	initialize stack for all boxes
!	------------------------------------------------

	id = 0
	do iy=1,nby
	  do ix=1,nbx
	    id = id + 1
	    call stack_init(idback)
	    if( id /= idback ) then
	      !this happens because we assumed that no other routine
	      !has already called stack ... make it more general
	      stop 'error stop fast_find_initialize: id/=idback'
	    end if
	  end do
	end do
	idmax = id

	if( bdebug ) write(6,*) 'idmax: ',idmax

!	------------------------------------------------
!	insert elements
!	------------------------------------------------

	do ie=1,nel

          do ii=1,3
            k = nen3v(ii,ie)
            x(ii) = xgv(k)
            y(ii) = ygv(k)
          end do

          xemin = min(x(1),x(2),x(3))
          yemin = min(y(1),y(2),y(3))
          xemax = max(x(1),x(2),x(3))
          yemax = max(y(1),y(2),y(3))
	  ixmin = 1 + (xemin-xbmin)/box_size
	  ixmax = 1 + (xemax-xbmin)/box_size
	  iymin = 1 + (yemin-ybmin)/box_size
	  iymax = 1 + (yemax-ybmin)/box_size

	  do iy=iymin,iymax
	    do ix=ixmin,ixmax
	      id = (iy-1)*nbx + ix
	      if( id < 1 .or. id > idmax ) then
		write(6,*) 'id,idmax: ',id,idmax
		stop 'error stop fast_find_initialize: internal error (1)'
	      end if
	      call stack_push(id,ie)
	    end do
	  end do

        end do

!	------------------------------------------------
!	convert from stack to element list
!	------------------------------------------------

	nemax = 0
	netot = 0
	neaver = 0
	nezero = 0
	do id=1,idmax
	  ne = stack_fill(id)
	  nemax = max(nemax,ne)
	  netot = netot + ne
	  if( ne == 0 ) nezero = nezero + 1
	  !write(6,*) id,ne
	end do
	neaver = netot / idmax

	if( bdebug ) then
	  write(6,'(a,6i8)') 'stack ne: ',idmax,nemax,netot,neaver,nezero,nel
	end if

	allocate(n_elem_list(idmax))
	allocate(elem_list(nemax,idmax))

	do id=1,idmax
	  call stack_get_entries(id,ns,elem_list(:,id))
	  n_elem_list(id) = ns
	  ne = stack_fill(id)
	  if( ne /= ns ) then
	    write(6,*) 'n,ne: ',ns,ne
	    stop 'error stop fast_find_initialize: internal error (w)'
	  end if
	end do

	do id=1,idmax
	  call stack_delete(id)
	end do

	if( bdebug ) write(6,*) 'finished seting up fast_find routine'

!	------------------------------------------------
!	end of routine
!	------------------------------------------------

	end

!******************************************************************

	subroutine fast_find_finalize

! delete data structure

	use basin

	if( idmax == 0 ) return

	nsize = 0
	idmax = 0
	deallocate(n_elem_list)
	deallocate(elem_list)

	end

!******************************************************************

	subroutine fast_find_search(x,y,ie)

! finds unique element for x,y - if not found returns ie==0

	use basin

	real x,y
	integer ie

	integer ix,iy,id
	integer ifound,n,iee,i,ie_ext
	integer, parameter :: iedim = 100
	integer ie_found(iedim)

	logical in_element

	if( idmax == 0 ) then
	  stop 'error stop fast_find_search: fast_find not initialized'
	end if

	ix = 1 + (x-xbmin)/box_size
	iy = 1 + (y-ybmin)/box_size
	id = (iy-1)*nbx + ix

	if( id < 1 .or. id > idmax ) then
	  write(6,*) 'id,idmax: ',id,idmax
	  stop 'error stop fast_find_search: internal error (1)'
	end if
	
	n = n_elem_list(id)
	ifound = 0

	do i=1,n
	  iee = elem_list(i,id)
	  if( in_element(iee,x,y) ) then
	    ifound = ifound + 1
	    if( ifound > iedim) stop 'error stop fast_find_search: iedim'
	    ie_found(ifound) = iee
	  end if
	end do

	if( ifound == 0 ) then
	  ie = 0
	else if( ifound == 1 ) then
	  ie = ie_found(1)
	else
	  ie = 0
	  ie_ext = 0
	  do i=1,ifound
	    iee = ie_found(i)
	    if( ipev(iee) > ie_ext ) then
	      ie = iee
	      ie_ext = ipev(iee)
	    end if
	  end do
	end if

	end

!******************************************************************

	subroutine fast_find_plot

	implicit none

	integer n,l,ix,iy
	real x,y,x1,x2,y1,y2

	n = 0
	l = 0

	!write(6,*) 'grid ',nbx,nby

	do ix=0,nbx
	  x = xbmin + ix*box_size
	  !write(6,*) ix,x
	  n = n + 1
	  y1 = ybmin
	  write(78,1000) 1,n,0,x,y1
	  n = n + 1
	  y2 = ybmax
	  write(78,1000) 1,n,0,x,y2
	  l = l + 1
	  write(78,2000) 3,l,0,2,n-1,n
	end do

	write(78,*)

	do iy=0,nby
	  y = ybmin + iy*box_size
	  !write(6,*) iy,y
	  n = n + 1
	  x1 = xbmin
	  write(78,1000) 1,n,0,x1,y
	  n = n + 1
	  x2 = xbmax
	  write(78,1000) 1,n,0,x2,y
	  l = l + 1
	  write(78,2000) 3,l,0,2,n-1,n
	end do

	return
 1000	format(i1,2i8,2e16.8)
 2000	format(i1,5i8)
	end

!==================================================================
        end module mod_fast_find
!==================================================================

