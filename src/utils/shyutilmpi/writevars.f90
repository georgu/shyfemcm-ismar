
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

! routines for writing selected nodes to file
!
! revision log :
!
! 24.09.2026    ggu     written from scratch
!
! notes :
!
! only the first layer (surface) is written to file
! works also for MPI
!
! typical call from shyfem:
!
!	call writevars_init(nnodes,nodes)	!just call it once at start
!
!	call writevars_do(iunit,n,nlvdi,var)	!write on every time step
!	
!	nnodes		total number of nodes to be written
!	nodes(nnodes)	external node numbers
!	iunit		unit number of the file (different for different vars)
!	n		horizonatl extension of grid (typically nkn)
!	nlvdi		vertical extension of grid (typically nlv)
!	var(nlvdi,nkn)	variables to be written
!
! example:
!
!   at simulation start:
!
!	nnodes = 3
!	nodes(1) = 100				!give nodes with external nums
!	nodes(2) = 1000
!	nodes(3) = 1500
!	call writevars_init(nnodes,nodes)	!just call it once at start
!
!   in time loop:
!
!	call writevars_do(200,nkn,nlv,tempv)	!writes temperature values
!
! please note that for different variables different file units have to be used

!=================================================================
	module mod_writevars
!=================================================================

	logical, parameter :: bflush = .false.	!flush after every write

	character*80 format
	integer, save :: ndim = 0
	integer, save, allocatable :: inodes(:)
	integer, save, allocatable :: enodes(:)
	integer, save, allocatable :: node_id(:)

!=================================================================
	end module mod_writevars
!=================================================================

	subroutine writevars_init(nnodes,nodes)

	use shympi
	use mod_writevars

	implicit none

	integer nnodes
	integer nodes(nnodes)

	logical berror,bmaster
	integer i,node
	integer, allocatable :: found(:)

	integer ipint

	if( ndim > 0 ) return			!already initialized

	ndim = nnodes
	allocate(inodes(ndim),enodes(ndim),node_id(ndim))
	allocate(found(ndim))
	inodes = 0
	node_id = 0
	found = 0
	
	write(format,'(a,i3,a)') '(a20,',ndim,'f12.4)'

	enodes = nodes			!save external nodes

	do i=1,ndim
	  node = ipint(nodes(i))
	  if( node == 0 ) cycle
	  if( id_node(node) /= my_id ) node = 0	!take node from my_id
	  inodes(i) = node
	  node_id(i) = id_node(node)
	  found(i) = 1
	end do

	call shympi_gather_and_sum(node_id)
	call shympi_gather_and_sum(found)

	bmaster = shympi_is_master()
	berror = .false.

	do i=1,ndim
	  if( bmaster ) then
	    write(6,*) 'writevars_init: ',i,enodes(i),node_id(i),found(i)
	  end if
	end do

	do i=1,ndim
	  if( found(i) /= 1 ) then
	    if( bmaster ) write(6,*) '*** no such node: ',enodes(i),found(i)
	    berror = .true.
	  end if
	end do

	if( berror .and. bmaster ) then
	  stop 'error stop writevars_init: unknown node(s)'
	end if

	call shympi_barrier

	end

!*****************************************************************

	subroutine writevars_do(iunit,n,nlvdi,var)

	use shympi
	use mod_writevars

	implicit none

	integer iunit			!unit to write variable
	integer n			!horizontal dimension of var (nkn)
	integer nlvdi			!vertical dimension of var (nlv)
	real var(nlvdi,n)		!variable to be written

	integer i,node
	real values(ndim)
	character*20 aline

!-----------------------------------------------------------------
! check if initialized
!-----------------------------------------------------------------

	if( ndim == 0 ) then
	  write(6,*) 'routines writevars are not initialized'
	  stop 'error stop writevars_do: not initialized'
	end if

!-----------------------------------------------------------------
! normal call - get values and write to file
!-----------------------------------------------------------------

	values = 0.
	do i=1,ndim
	  node = inodes(i)
	  if( node > 0 ) values(i) = var(1,node)
	end do

	call shympi_gather_and_sum(values)

	if( shympi_is_master() ) then
	  call get_act_timeline(aline)
	  write(iunit,format) aline,values
	  if( bflush ) flush(iunit)
	end if

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

	end

!*****************************************************************

