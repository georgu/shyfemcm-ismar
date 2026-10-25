
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017-2019  Georg Umgiesser
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

!------------------------------------------------------------------
!
! connectivity routines
!
! revision log :
!
! 19.05.2020    ccf     started from scratch
! 12.04.2022    ggu     adapted
! 25.10.2026    ggu     new calls to floodfill routines (nodes not yet finished)
!
! notes :
!
! module mod_connect must be initialized
!
! usage :
!
! floodfill_elem_all(nel,nareas,iareas)		!fills whole domain
! floodfill_elem(nel,acode,ac,nareas,iareas)	!fills domain with code ac
!
!------------------------------------------------------------------

!==================================================================
        module mod_flood
!==================================================================

        implicit none

!==================================================================
        contains
!==================================================================

!==================================================================
        end module mod_flood
!==================================================================

!******************************************************************
!******************************************************************
!******************************************************************
! element routines
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine floodfill_elem_all(nel,nareas,iareas)

! fills whole domain and returns nareas and iareas(:)

	implicit none

	integer, intent(in) :: nel		!total number of elements
	integer, intent(out) :: nareas		!total number of domains found
	integer, intent(out) :: iareas(nel)	!domain number for each element

	integer ie,ic,icolor
	integer acode(nel)	

	nareas = 0
	iareas = 0
	acode = 0		!make sure to fill whole domain
	icolor = 0

	do

	  do ie=1,nel
	    if( iareas(ie) == 0 ) exit
	  end do

	  !write(6,*) 'checking element ',ie,icolor

	  if( ie > nel ) exit
	  icolor = icolor + 1

	  call fill_elem_area(nel,acode,ie,icolor,iareas)

	end do

	nareas = icolor

	ic = count( iareas == 0 )
	if( ic > 0 ) stop 'error stop floodfill_elem_all: internal error'

	end

!******************************************************************

	subroutine floodfill_elem(nel,acode,ac,nareas,iareas)

! fills domain with element area code ac and returns nareas and iareas(:)

	implicit none

	integer, intent(in) :: nel		!total number of elements
	integer, intent(in) :: acode(nel)	!area code of elements
	integer, intent(in) :: ac		!area code to fill
	integer, intent(out) :: nareas		!number of connected domains
	integer, intent(out) :: iareas(nel)	!domain number for each element

	integer ie,icolor

	nareas = 0
	iareas = 0

	do

	  do ie=1,nel
	    if( acode(ie) == ac .and. iareas(ie) == 0 ) exit
	  end do

	  if( ie > nel ) exit
	  nareas = nareas + 1

	  icolor = nareas
	  call fill_elem_area(nel,acode,ie,icolor,iareas)

	end do

	end

!******************************************************************

	subroutine fill_elem_area(nel,acode,ie,icolor,aaux)

! internal routine - not to be called directly
!
! flags aaux with icolor if can be reached from element ie and has acode of ie
! aaux must be 0 at entry
!
! at return, all elements filled have aaux(:) == icolor

	use queue
	use mod_connect

	implicit none

	integer, intent(in) :: nel		!total number of elements
	integer, intent(in) :: acode(nel)	!area code of elements
	integer, intent(in) :: ie		!element to start with
	integer, intent(in) :: icolor		!color to insert into aaux
	integer, intent(inout) :: aaux(nel)	!array to show filling

	!aaux(ie)==0 if not yet filled, >0 otherwise

	integer ac,ii,ieo,ien,id
	integer, allocatable :: ecv(:,:)

	ac = acode(ie)

	call queue_init(id)
	call queue_enqueue(id,ie)

	allocate(ecv(3,nel))
	call connect_get_ecv(nel,ecv)

	do

	  if( .not. queue_dequeue(id,ieo) ) exit

	  if( aaux(ieo) > 0 ) cycle
	  aaux(ieo) = icolor

	  do ii=1,3
	    ien = ecv(ii,ieo)
	    if( ien > 0 ) then
	      if( acode(ien) == ac .and. aaux(ien) == 0 ) then
		call queue_enqueue(id,ien)
	      end if
	    end if
	  end do

	end do

	call queue_delete(id)

	end

!******************************************************************
!******************************************************************
!******************************************************************
! node routines
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine floodfill_node(nkn,acode,ac,nareas)

! fills domain with element area code ac and returns nareas and iareas(:)

	implicit none

	integer nkn		!total number of nodes
	integer acode(nkn)	!area code of nodes
	integer ac		!area code to look for
	integer nareas		!number of connected domains with area code ac

	integer ainfo(2,1)

	nareas = 0
	call floodfill_node_info(nkn,acode,ac,nareas,ainfo)

	end

!******************************************************************

	subroutine floodfill_node_info(nkn,acode,ac,nareas,ainfo)

	implicit none

	integer nkn		!total number of nodes
	integer acode(nkn)	!area code of nodes
	integer ac		!area code to look for
	integer nareas		!number of connected domains with area code ac
	integer ainfo(2,nareas)	!info on domains of color ac (return)

	integer k,kstart
	integer ndim
	integer nnode
	integer aaux(nkn)	!aux array to remember what has been filled

	ndim = nareas

	nareas = 0
	aaux = 0
	kstart = 1

	do

	  do k=kstart,nkn
	    if( acode(k) == ac .and. aaux(k) == 0 ) exit
	  end do

	  if( k > nkn ) return

	  call fill_node_area(nkn,acode,k,aaux,nnode)
	  kstart = k

	  nareas = nareas + 1
	  if( ndim > 0 ) then
	    if( nareas > ndim ) goto 99
	    ainfo(1,nareas) = k
	    ainfo(2,nareas) = nnode
	  end if

	end do

	return
   99	continue
	stop 'error stop floodfill_node_info: iareas>ndim'
	end

!******************************************************************

	subroutine floodfill_node_color(nkn,acode,acnew,k)

	implicit none

	integer nkn		!total number of nodes
	integer acode(nkn)	!area code of nodes
	integer acnew		!new area code to color area
	integer k		!node to start from

	integer nnode
	integer aaux(nkn)	!aux array to remember what has been filled

	aaux = 0
	call fill_node_area(nkn,acode,k,aaux,nnode)
	where( aaux > 0 ) acode = acnew

	end

!******************************************************************

	subroutine fill_node_area(nkn,acode,k,aaux,nnode)

! flags aaux with 1 if it can be reached from node k and has acode of k

	use queue
	use mod_connect

	implicit none

	integer nkn		!total number of nodes
	integer acode(nkn)	!area code of nodes
	integer k		!node to start with
	integer aaux(nkn)	!aux array to show what has been filled
	integer nnode		!number of nodes in area (return)

	integer ac,i,n,ko,kn,id
	integer ngr
	integer, allocatable :: nlist(:,:)
	integer, allocatable :: elist(:,:)

	nnode = 0
	ac = acode(k)

	call queue_init(id)
	call queue_enqueue(id,k)

	call connect_get_grade(ngr)
	allocate(nlist(0:ngr,nkn))
	allocate(elist(0:ngr,nkn))
	call connect_get_lists(nkn,ngr,nlist,elist)

	do

	  if( .not. queue_dequeue(id,ko) ) exit

	  aaux(ko) = aaux(ko) + 1
	  if( aaux(ko) > 1 ) cycle
	  nnode = nnode + 1

	  n = nlist(0,ko)
	  do i=1,n
	    kn = nlist(i,ko)
	    if( acode(kn) == ac .and. aaux(kn) == 0 ) then
	      call queue_enqueue(id,kn)
	    end if
	  end do

	end do

	call queue_delete(id)

	end

!******************************************************************

