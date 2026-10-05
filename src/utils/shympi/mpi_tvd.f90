
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

! tvd routines
!
! revision log :
!
! 02.08.2024	ggu	started mpi tvd
! 28.08.2024	ggu	first part finding elements finished
! 03.09.2024	ggu	more on mpi-tvd
! 10.09.2024	ggu	data structures localized
! 27.09.2024	ggu	first running version... still to be cleaned
! 05.10.2024	ggu	introduced quad_tree_search() - feature complete
! 07.10.2024	ggu	tvd mpi finished
!
! notes :
!
! tvd_mpi_init		initializes the mpi tvd variables
! tvd_mpi_prepare	prepares mpi variables for every call to scalar advect.
!
! still to do :
!
! prepare upwind c for more than one variable

!==================================================================
        module shympi_tvd
!==================================================================

	implicit none

	logical, parameter :: bmpitvd = .true. !handles tvd

	integer, parameter :: ind_x = 1
	integer, parameter :: ind_y = 2
	integer, parameter :: ind_ia_this = 3
	integer, parameter :: ind_j_this = 4
	integer, parameter :: ind_ii_this = 5
	integer, parameter :: ind_k_this = 6
	integer, parameter :: ind_ie_this = 7
	integer, parameter :: ind_iee_this = 8
	integer, parameter :: ind_lmax_this = 9
	integer, parameter :: ind_ia_found = 10
	integer, parameter :: ind_ie_found = 11
	integer, parameter :: ind_iee_found = 12
	integer, parameter :: ind_lmax_found = 13

	integer, parameter :: nlist = 13
	integer, parameter :: nlist_type = 2

	integer, save :: n_tvd_receive_total,n_tvd_send_total
	integer, save :: n_my_array_list,n_my_receive_list,n_my_send_list
	integer, save, allocatable :: my_receive_list(:)
	integer, save, allocatable :: my_send_list(:)
	real, save, allocatable :: my_array_list(:,:)
	double precision, save, allocatable :: my_xi_list(:,:)
	real, save, allocatable :: values_tvd_aux(:,:)
	real, save, allocatable :: values_tvd_local(:,:)
	real, save, allocatable :: values_tvd_remote(:,:)
	real, save, allocatable :: tvd_buffer_in(:)
	real, save, allocatable :: tvd_buffer_out(:)
	real, save, allocatable :: tvd_mbuffer_in(:)
	real, save, allocatable :: tvd_mbuffer_out(:)

	integer, save :: nblocks_multi
	integer, save, allocatable :: exchange_multi_blocks(:,:)

! btvddebug regulates writing of all debug files:
!	fort.34*, fort.50[01], tvd_debug.txt
! btvdassert regulates assesing vital relationships between variables
! iu_debug writes debug output to fort.500
! iuunf_debug writes debug output to fort.501
! btvddebug could also be set in files concentration.f90 and tvd_admin.f90
!
! for production runs set btvddebug and btvdassert to .false.

	logical, save :: btvddebug = .false.	!writes tvd_debug.txt
	logical, save :: btvdassert = .false.	!asserts important values
	logical, save :: bcheck_multi = .false.	!checks multi exchange
	logical, save :: bdebout = .false.	!allows for writing to fort.500
	integer, save :: iu_debug = 500		!writes fort.500 if > 0
	integer, save :: iuunf_debug = 501	!writes fort.501 if > 0
	integer, save :: ifreq_debug = 1	!frequency of debug in fort.500
	real, save, allocatable :: tvd_debug_aux(:,:,:)
	logical, save, allocatable :: blevel(:)

!==================================================================
        contains
!==================================================================

	subroutine jii2k(j,ii,k)

	integer j,ii,k

	integer, parameter :: mat(3,3) = reshape((/0,1,2,3,0,4,5,6,0/),(/3,3/))
	
	k = mat(j,ii)

	end subroutine

!******************************************************************

	subroutine k2jii(k,j,ii)

	integer j,ii,k

	integer, parameter :: mat(2,6) = &
     &		reshape((/2,1,3,1,1,2,3,2,1,3,2,3/),(/2,6/))
	
	j = mat(1,k)
	ii = mat(2,k)

	end subroutine

!==================================================================
        end module shympi_tvd
!==================================================================

	subroutine tvd_mpi_init

! this sets up all data structures needed for tvd mpi

	use basin
	use levels
	use shympi
	use shympi_tvd
	use mod_tvd
	use mod_quad_tree

	implicit none

	logical bdebug,bfound
	integer ie,ii,j,i,ie_new,ix,n,ind,iaf,ian,k
	integer my_ia,ia,id,inlist,nfound,nchanged
	integer ie_ext,iee_ext,iee_old,ie_int,iee,ieq
	integer ilist,maxlist,maxnlist
	integer ia_found,ia_needed,ie_local,iudb
	integer n_tvd_r,n_tvd_s,n_tvd_a,n_tvd
	integer lmax
	real x,y
	real xmin,ymin,xmax,ymax

	real, allocatable :: raux(:,:)
	real, allocatable :: rlist(:,:)
	real, allocatable :: newlist(:,:)
	real, allocatable :: rlists(:,:,:)
	real, allocatable :: newlists(:,:,:)
	integer, allocatable :: count(:)
	integer, allocatable :: index(:)
	integer, allocatable :: n_tvd_receive(:)
	integer, allocatable :: n_tvd_send(:)
	integer, allocatable :: n_tvd_all(:)
	integer, allocatable :: n_tvd_index(:)
	integer, allocatable :: tvd_receive(:,:)
	integer, allocatable :: tvd_send(:,:)
	integer, allocatable :: tvd_index(:,:)
	integer, allocatable :: matrix(:,:)
	integer, allocatable :: n_tvd_receive_from(:)
	integer, allocatable :: n_tvd_send_to(:)
	integer :: n_tvd_tot
	integer :: na,nr,ns
	double precision xi(3),xd,yd

!----------------------------------------------------------
! starting collection of tvd information
!----------------------------------------------------------

	if( .not. bmpi ) return
	if( .not. bmpitvd ) return		!do not handle mpi tvd

	write(6,*) 'configuring tvd with mpi ',my_id

	bdebug = btvddebug
	my_ia = my_id + 1
	iudb = 340 + my_ia

	ilist = 0
	allocate(raux(nlist,nel*6))
	raux = 0.

!----------------------------------------------------------
! setting up initial list with all information of this domain
!----------------------------------------------------------

	do ie=1,nel
	  ie_ext = ipev(ie)
	  do ii=1,3
	    do j=1,3
	      if( ii == j ) cycle
	      ilist = ilist + 1
	      call jii2k(j,ii,k)
	      raux(ind_x,ilist) = xtvdup(j,ii,ie)
	      raux(ind_y,ilist) = ytvdup(j,ii,ie)
	      raux(ind_ia_this,ilist) = my_ia		!needed in this domain
	      raux(ind_j_this,ilist) = j
	      raux(ind_ii_this,ilist) = ii
	      raux(ind_k_this,ilist) = k
	      raux(ind_ie_this,ilist) = ie		!internal local element
	      raux(ind_iee_this,ilist) = ipev(ie)	!external local element
	      raux(ind_lmax_this,ilist) = ilhv(ie)	!lmax
	    end do
	  end do
	end do

	call shympi_barrier

!----------------------------------------------------------
! find biggest list information and exchange
!----------------------------------------------------------

	write(6,*) 'exchanging list information ',my_id

	maxlist = shympi_max(ilist)

	write(6,*) 'array sizes: ',ilist,maxlist,my_id
	if( bdebug ) then
	  write(iudb,*) 'debug output from tvd_mpi_handle for area ',my_ia
	  write(iudb,*) 'array sizes: ',ilist,maxlist,my_ia
	end if

	allocate(rlist(nlist,maxlist))
	allocate(rlists(nlist,maxlist,n_threads))
	rlist = 0.

	rlist(:,1:ilist) = raux(:,1:ilist)
	call shympi_gather(nlist,rlist,rlists)

	!------------------------------------------------
	! in rlist/rlists is information on all upwind points
	!------------------------------------------------

	call shympi_barrier
	write(6,*) 'all information collected ',my_id

!----------------------------------------------------------
! look for missing information needed in other domains
!----------------------------------------------------------

	call bas_get_minmax(xmin,ymin,xmax,ymax)

	nfound = 0

	deallocate(raux)
	allocate(raux(nlist,maxlist*n_threads))
	raux = 0.
	inlist = 0

	do ia=1,n_threads
	  if( ia == my_ia ) cycle
	  do i=1,maxlist
	    ie_local = nint(rlists(ind_ie_this,i,ia))
	    ia_needed = nint(rlists(ind_ia_this,i,ia))
	    ie_ext = nint(rlists(ind_iee_this,i,ia))
	    ie = 0
	    iee_ext = 0
	    if( ie_local == 0 ) exit			!end of list
	    x = rlists(ind_x,i,ia)
	    y = rlists(ind_y,i,ia)
	    if( x < xmin .or. x > xmax ) cycle
	    if( y < ymin .or. y > ymax ) cycle
	    call quad_tree_search(x,y,ieq)
	    ie = ieq
	    !call find_unique_element(x,y,ie)
	    if( ie /= ieq ) stop 'error stop tvd_mpi_init: ie/=ieq'
	    if( ie > 0 ) then				!element found
	      nfound = nfound + 1
	      iee_ext = ipev(ie)
	      rlists(ind_ia_found,i,ia) = my_ia		!found for this domain
	      rlists(ind_ie_found,i,ia) = ie		!internal element number
	      rlists(ind_iee_found,i,ia) = iee_ext	!external element number
	      rlists(ind_lmax_found,i,ia) = ilhv(ie)	!lmax
	      inlist = inlist + 1
	      raux(:,inlist) = rlists(:,i,ia)
	    end if
	  end do
	end do

	!------------------------------------------------
	! in raux are points found in this domain needed in other domains
	!------------------------------------------------

	call shympi_barrier

!----------------------------------------------------------
! exchange information between domains
!----------------------------------------------------------

	maxnlist = shympi_max(inlist)

	if( bdebug ) write(iudb,*) 'new list info ',inlist,maxnlist,my_ia
	
	call shympi_barrier

	allocate(newlist(nlist,maxnlist))
	allocate(newlists(nlist,maxnlist,n_threads))
	newlist = 0.
	newlists = 0.
	newlist(:,1:inlist) = raux(:,1:inlist)

	call shympi_gather(nlist,newlist,newlists)

	!------------------------------------------------
	! in newlist/newlists are all points that need exchange
	!------------------------------------------------

	call shympi_barrier
	write(6,*) 'all new information collected ',my_id

!----------------------------------------------------------
! finding unique element across domains
!----------------------------------------------------------

	nchanged = 0

	do ia=1,n_threads
	  if( ia == my_ia ) cycle
	  do i=1,maxnlist
	    ie = nint(newlists(ind_ie_found,i,ia))
	    if( ie == 0 ) exit				 !end of list
	    ia_found = nint(newlists(ind_ia_found,i,ia)) !found in this domain
	    ia_needed = nint(newlists(ind_ia_this,i,ia)) !needed in this domain
	    ie_int = nint(newlists(ind_ie_found,i,ia))	 !in this element found
	    iee_ext = nint(newlists(ind_iee_found,i,ia)) !in this element found
	    lmax = nint(newlists(ind_lmax_found,i,ia))   !lmax of remote
	    if( ia_found /= ia_needed ) then
	     if( my_ia == ia_needed ) then
	      nchanged = nchanged + 1
	      j = nint(newlists(ind_j_this,i,ia))
	      ii = nint(newlists(ind_ii_this,i,ia))
	      ie = nint(newlists(ind_ie_this,i,ia))
	      iee_old = ieetvdup(j,ii,ie)
	      if( iee_ext > iee_old ) then		 !always take highest
	        ieetvdup(j,ii,ie) = iee_ext
	        ietvdup(j,ii,ie) = ie_int
	        iatvdup(j,ii,ie) = ia_found
	        ltvdup(j,ii,ie) = lmax
	      end if
	     end if
	    end if
	  end do
	end do

	call shympi_barrier
	write(6,*) 'all new information elaborated ',nchanged,my_id

!----------------------------------------------------------
! prepare data exchange
!----------------------------------------------------------

	deallocate(raux)
	allocate(raux(nlist,9*nel))
	raux = 0.

	ilist = 0
	do ie=1,nel
	  do ii=1,3
	    do j=1,3
	        if( ii == j ) cycle
	        ia = iatvdup(j,ii,ie)		!area where point found
		if( ia == my_ia ) cycle
		ilist = ilist + 1
		call jii2k(j,ii,k)
		raux(ind_x,ilist) = xtvdup(j,ii,ie)
		raux(ind_y,ilist) = ytvdup(j,ii,ie)
		raux(ind_ia_this,ilist) = my_ia	!needed in this domain
		raux(ind_j_this,ilist) = j
		raux(ind_ii_this,ilist) = ii
		raux(ind_k_this,ilist) = k
		raux(ind_ie_this,ilist) = ie		!internal local element
		raux(ind_iee_this,ilist) = ipev(ie)	!external local element
	        raux(ind_ia_found,ilist) = iatvdup(j,ii,ie)!area where found
	        raux(ind_ie_found,ilist) = ietvdup(j,ii,ie)!internal remote
	        raux(ind_iee_found,ilist) = ieetvdup(j,ii,ie)!external remote
	        raux(ind_lmax_this,ilist) = ilhv(ie) 	!lmax local
	        raux(ind_lmax_found,ilist) = ltvdup(j,ii,ie)!lmax remote
	    end do
	  end do
	end do
	  
	!------------------------------------------------
	! in raux is all available information
	!------------------------------------------------

	call assert_list_coherence(ilist,raux,.true.)
	call shympi_barrier

!------------------------------------------------
! prepare count and index for sorting
!------------------------------------------------

	allocate(count(n_threads))
	allocate(index(n_threads))
	count = 0

	do i=1,ilist
	  ia = nint(raux(ind_ia_found,i))
	  count(ia) = count(ia) + 1
	end do
	call tvd_assert('count(my_ia) /= 0',count(my_ia) == 0)

	ix = 0
	do ia=1,n_threads
	  index(ia) = ix
	  ix = ix + count(ia)
	end do
	call tvd_assert('ix /= ilist',ix == ilist)

	maxlist = shympi_max(ilist)

	if( bdebug ) then
	  write(iudb,*) 'new array sizes: ',ilist,maxlist,my_ia
	  write(iudb,*) 'count: ',count
	  write(iudb,*) 'index: ',index
	end if

	deallocate(rlist)
	deallocate(rlists)
	allocate(rlist(nlist,maxlist))
	allocate(rlists(nlist,maxlist,n_threads))

!------------------------------------------------
! sort raux and copy to rlist
! after this in rlist are the same entries as in raux, but sorted by ia
!------------------------------------------------

	if( bdebug .and. ilist > 0 ) then
	  write(iudb,*) 'not sorted list: ',my_ia,ilist,maxlist
	  call write_array_header(iudb)
	  do i=1,ilist
	    call write_array_debug(iudb,i,raux(:,i))
	  end do
	end if

	rlist = 0.

	do i=1,ilist
	  ia = nint(raux(ind_ia_found,i))
	  ix = index(ia) + 1
	  rlist(:,ix) = raux(:,i)
	  index(ia) = ix
	end do

	if( bdebug .and. ilist > 0 ) then
	  write(iudb,*) 'sorted list: ',my_ia,ilist,maxlist
	  call write_array_header(iudb)
	  do i=1,ilist
	    call write_array_debug(iudb,i,rlist(:,i))
	  end do
	end if

	call assert_list_coherence(ilist,rlist,.true.)

	call shympi_gather(nlist,rlist,rlists)

!------------------------------------------------
! compute total number of points to receive or send
!------------------------------------------------

	allocate(n_tvd_receive(n_threads))
	allocate(n_tvd_send(n_threads))
	allocate(n_tvd_all(n_threads))

	n_tvd_receive = 0
	n_tvd_send = 0
	n_tvd_all = 0

	do ia=1,n_threads
	  do i=1,maxlist
	    ie = nint(rlists(ind_iee_found,i,ia))
	    if( ie == 0 ) exit				!end of list
	    ia_found = nint(rlists(ind_ia_found,i,ia))	!found in this domain
	    if( ia_found == ia ) cycle
	    if( ia_found == my_ia ) then
	      n_tvd_send(ia) = n_tvd_send(ia) + 1
	    else
	      n_tvd_receive(ia) = n_tvd_receive(ia) + 1
	    end if
	    n_tvd_all(ia) = n_tvd_all(ia) + 1
	  end do
	end do

	n_tvd_r = maxval(n_tvd_receive)
	n_tvd_s = maxval(n_tvd_send)
	n_tvd_a = maxval(n_tvd_all)
	n_tvd = max(n_tvd_r,n_tvd_s)

	write(6,*) 'dims: ',n_tvd_r,n_tvd_s,n_tvd_a,n_tvd,my_ia

!------------------------------------------------
! set up receive and send list index
!------------------------------------------------

	allocate(tvd_receive(n_tvd_r,n_threads))
	allocate(tvd_send(n_tvd_s,n_threads))
	allocate(tvd_index(n_tvd,n_threads))
	allocate(n_tvd_index(n_threads))
	allocate(matrix(n_threads,n_threads))
	allocate(n_tvd_send_to(n_threads))
	allocate(n_tvd_receive_from(n_threads))
	n_tvd_receive = 0
	n_tvd_send = 0
	n_tvd_index = 0
	tvd_receive = 0
	tvd_send = 0
	tvd_index = 0
	matrix = 0
	n_tvd_send_to = 0
	n_tvd_receive_from = 0
	n_tvd_tot = 0
	n_my_array_list = 0
	n_my_receive_list = 0
	n_my_send_list = 0

	do ia=1,n_threads
	  do i=1,maxlist
	    ie = nint(rlists(ind_ie_this,i,ia))
	    if( ie == 0 ) exit				!end of list
	    ia_needed = nint(rlists(ind_ia_this,i,ia))  !needed in this domain
	    ia_found = nint(rlists(ind_ia_found,i,ia))	!found in this domain
	    if( ia_found == ia_needed ) cycle
	    if( ia_found == my_ia ) then
	      n = n_tvd_send(ia) + 1
	      if( n > n_tvd_s ) stop 'error stop tvd_mpi_init: n > n_tvd_s'
	      tvd_send(n,ia) = i
	      n_tvd_send(ia) = n
	      n = n_tvd_index(ia) + 1
	      if( n > n_tvd ) stop 'error stop tvd_mpi_init: n > n_tvd'
	      tvd_index(n,ia) = i
	      n_tvd_index(ia) = n
	      n_tvd_send_to(ia_needed) =  n_tvd_send_to(ia_needed) + 1
	    else if( ia_needed == my_ia ) then
	      n = n_tvd_receive(ia) + 1
	      if( n > n_tvd_r ) stop 'error stop tvd_mpi_init: n > n_tvd_r'
	      tvd_receive(n,ia) = i
	      n_tvd_receive(ia) = n
	      n = n_tvd_index(ia) + 1
	      if( n > n_tvd ) stop 'error stop tvd_mpi_init: n > n_tvd'
	      tvd_index(n,ia) = i
	      n_tvd_index(ia) = n
	      n_tvd_receive_from(ia_found) =  n_tvd_receive_from(ia_found) + 1
	    end if
	    if( ia_found == my_ia .or. ia_needed == my_ia ) then
	      n_my_array_list = n_my_array_list + 1
	      if( ia_found == my_ia ) then
	        n_my_send_list = n_my_send_list + 1
	      else
	        n_my_receive_list = n_my_receive_list + 1
	      end if
	    end if
	    matrix(ia_found,ia_needed) = matrix(ia_found,ia_needed) + 1
	    if( ia_needed /= ia_found ) then
	      n_tvd_tot = n_tvd_tot + 1
	    end if
	  end do
	end do

	n_tvd_receive_total = sum( n_tvd_receive )
	n_tvd_send_total = sum( n_tvd_send )

	call tvd_assert('list sizes' &
     &		,n_my_receive_list + n_my_send_list == n_my_array_list)

!------------------------------------------------
! finally set up receive, send, and array list
! my_receive_list and my_send_list are indices into my_array_list
! in my_array_list all points that have to be exchanged are contained
!------------------------------------------------

	allocate(my_receive_list(n_my_receive_list))
	allocate(my_send_list(n_my_send_list))
	allocate(my_array_list(nlist,n_my_array_list))
	my_receive_list = 0
	my_send_list = 0
	my_array_list = 0.

	nr = 0
	ns = 0
	na = 0
	do ia=1,n_threads
	  do i=1,maxlist
	    ie = nint(rlists(ind_ie_found,i,ia))
	    if( ie == 0 ) exit				!end of list
	    ia_needed = nint(rlists(ind_ia_this,i,ia))	!needed in this domain
	    ia_found = nint(rlists(ind_ia_found,i,ia))	!found in this domain
	    if( ia_found == ia_needed ) cycle
	    if( ia_needed == my_ia ) then
	      na = na + 1
	      my_array_list(:,na) = rlists(:,i,ia)
	      nr = nr + 1
	      my_receive_list(nr) = na
	    else if( ia_found == my_ia ) then
	      na = na + 1
	      my_array_list(:,na) = rlists(:,i,ia)
	      ns = ns + 1
	      my_send_list(ns) = na
	    end if
	  end do
	end do

	call tvd_assert('na /= n_my_array_list',na == n_my_array_list)
	call tvd_assert('nr /= n_my_receive_list',nr == n_my_receive_list)
	call tvd_assert('ns /= n_my_send_list',ns == n_my_send_list)

!------------------------------------------------
! setting up xi values and itvdup array
!------------------------------------------------

	allocate(my_xi_list(3,n_my_send_list))
	my_xi_list = 0.

	if( bdebug ) then
	  write(iudb,*) '====================================='
	  write(iudb,*) 'sending xi list'
	  write(iudb,*) '====================================='
	end if

	do i=1,n_my_send_list
	  na = my_send_list(i)
	  xd = my_array_list(ind_x,na)
	  yd = my_array_list(ind_y,na)
	  ie = my_array_list(ind_ie_found,na)
	  iee = my_array_list(ind_iee_found,na)
	  call xy2xi(ie,xd,yd,xi)
	  my_xi_list(:,i) = xi(:)
	  if( bdebug ) write(iudb,'(i6,2e14.6,3e14.6)') iee,xd,yd,xi
	end do

	itvdup(:,:,:) = 0
	do i=1,n_my_receive_list
	  na = my_receive_list(i)
	  j = my_array_list(ind_j_this,na)
	  ii = my_array_list(ind_ii_this,na)
	  ie = my_array_list(ind_ie_this,na)
	  itvdup(j,ii,ie) = i
	end do

	call shympi_barrier
	write(6,*) 'finished configuring tvd with mpi ',my_id
	flush(6)

	call assert_list_coherence(n_my_array_list,my_array_list,.true.)

!------------------------------------------------
! beyond here only debug output
!------------------------------------------------

	if( .not. bdebug ) return

	write(iudb,*) '====================================='
	write(iudb,*) 'my list:'
	write(iudb,*) '====================================='
	write(iudb,*) 'my array list of area: ',my_ia,n_my_array_list
	call write_array_header(iudb)
	do i=1,n_my_array_list
	  call write_array_debug(iudb,i,my_array_list(:,i))
	end do
	write(iudb,*) 'my receive list of area: ',my_ia,n_my_receive_list
	call write_array_header(iudb)
	do i=1,n_my_receive_list
	  na = my_receive_list(i)
	  call write_array_debug(iudb,i,my_array_list(:,na))
	end do
	write(iudb,*) 'my send list of area: ',my_ia,n_my_send_list
	call write_array_header(iudb)
	do i=1,n_my_send_list
	  na = my_send_list(i)
	  call write_array_debug(iudb,i,my_array_list(:,na))
	end do

	write(iudb,*) '====================================='
	write(iudb,*) 'send/receive index start ',my_ia
	write(iudb,*) n_tvd_r,n_tvd_s,n_tvd

	write(iudb,*) 'receive max: ',n_tvd_r
	write(iudb,*) 'receive dim: ',n_tvd_receive
	write(iudb,*) 'receive tot: ',n_tvd_receive_total
	do ia=1,n_threads
	  n = n_tvd_receive(ia)
	  if( n == 0 ) cycle
	  write(iudb,*) 'receiving from... n = ',n
	  call write_array_header(iudb)
	  do i=1,n
	    ind = tvd_receive(i,ia)
	    call write_array_debug(iudb,i,rlists(:,ind,ia))
	  end do
	end do
	write(iudb,*) 'send max: ',n_tvd_s
	write(iudb,*) 'send dim: ',n_tvd_send
	write(iudb,*) 'send tot: ',n_tvd_send_total
	do ia=1,n_threads
	  n = n_tvd_send(ia)
	  if( n == 0 ) cycle
	  write(iudb,*) 'sending to... n = ',n
	  call write_array_header(iudb)
	  do i=1,n
	    ind = tvd_send(i,ia)
	    call write_array_debug(iudb,i,rlists(:,ind,ia))
	  end do
	end do
	do ia=1,n_threads
	  n = n_tvd_index(ia)
	  if( ia == my_ia ) then
	    if(any(tvd_receive(1:n,ia)/=tvd_index(1:n,ia))) then
		write(6,*) 'receive,index: ',n_tvd_index(ia),n_tvd_receive(ia)
		stop 'error stop tvd_mpi_init: receive'
	    end if
	  else
	    if(any(tvd_send(1:n,ia)/=tvd_index(1:n,ia))) then
		write(6,*)'send,index: ', n_tvd_index(ia),n_tvd_send(ia)
		stop 'error stop tvd_mpi_init: send'
	    end if
	  end if
	end do
	write(iudb,*) '====================================='
	write(iudb,*) 'matrix exchange in domain ',my_ia
	do iaf=1,n_threads
	  do ian=1,n_threads
	    n = matrix(iaf,ian)
	    if( n == 0 ) cycle
	    write(iudb,'(i5,a,i5,a,i5)') iaf,' -> ',ian,' : ',n
	  end do
	end do
	write(iudb,*) '====================================='
	write(iudb,*) 'matrix all'
	write(iudb,'(a,10i5)') ' from/to:',(i,i=1,n_threads)
	do iaf=1,n_threads
	  write(iudb,'(i5,4x,10i5)') iaf,(matrix(iaf,ian),ian=1,n_threads)
	end do
	write(iudb,*) '====================================='
	do ia=1,n_threads
	  write(iudb,9) my_ia,' receives from ',ia,' : ',n_tvd_receive_from(ia)
	end do
	do ia=1,n_threads
	  write(iudb,9) my_ia,' sends to      ',ia,' : ',n_tvd_send_to(ia)
	end do
	write(iudb,*) '====================================='

	call shympi_barrier

!----------------------------------------------------------
! all info exchanged
!----------------------------------------------------------

	call shympi_barrier
	call flush(iudb)

!----------------------------------------------------------
! end of routine
!----------------------------------------------------------

	return
    9	format(i5,a,i5,a,i10)
	end

!******************************************************************
!******************************************************************
!******************************************************************

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine allocate_tvd_arrays(nvals)

! allocates global arrays

	use shympi
	use shympi_tvd

	implicit none

	integer nvals
	integer lmax

	if( allocated(values_tvd_local) ) return

	lmax = nlv_global

	allocate( values_tvd_aux(lmax,n_tvd_receive_total) )
	allocate( values_tvd_local(lmax,n_tvd_receive_total) )
	allocate( values_tvd_remote(lmax,n_tvd_send_total) )
	allocate( tvd_buffer_in(lmax) )
	allocate( tvd_buffer_out(lmax) )
	allocate( tvd_mbuffer_in(lmax*n_my_array_list) )
	allocate( tvd_mbuffer_out(lmax*n_my_array_list) )

	end

!******************************************************************

	subroutine tvd_mpi_run(values)

! this routine must be called before transport and diffusion of scalars
! handles preparation of remote values and exchange

	use basin
	use levels
	use shympi
	use shympi_tvd

	implicit none

	real values(nlvdi,nkn)

	integer, save :: nvals = 1		!needed later for multi-var

	if( .not. bmpi ) return
	if( .not. bmpitvd ) return		!do not handle mpi tvd

	call allocate_tvd_arrays(nvals)		!executed only first time

	call tvd_prepare_remote(values)		!prepare remote values

	if( bcheck_multi ) call tvd_exchange_single

	call tvd_exchange_prepare_multi		!executed only first time
	call tvd_exchange_multi			!exchange values between domains

	!call tvd_exchange_debug

	end

!******************************************************************

	subroutine tvd_prepare_remote(values)

! prepares remote values to be exchanged later

	use basin
	use levels
	use shympi
	use shympi_tvd

	implicit none

	real values(nlvdi,nkn)

	integer i,na,ie,ii,k,l,ia,iee
	integer lmax,lmax_remote,lmax_this
	integer my_ia
	double precision c,cacum
	double precision xi(3)

	my_ia = my_id + 1

	do i=1,n_my_send_list
	  na = my_send_list(i)
	  ia = nint(my_array_list(ind_ia_found,na))
	  ie = nint(my_array_list(ind_ie_found,na))
	  iee = nint(my_array_list(ind_iee_found,na))
	  lmax = ilhv(ie)
	  lmax_remote = nint(my_array_list(ind_lmax_found,na))
	  lmax_this = nint(my_array_list(ind_lmax_this,na))
	  if( btvdassert ) then
	    call tvd_assert('tvd_prepare_remote: (1)',ia==my_ia)
	    call tvd_assert('tvd_prepare_remote: (2)',ie<=nel)
	    call tvd_assert('tvd_prepare_remote: (3)',lmax==lmax_remote)
	  end if
	  xi(:) = my_xi_list(:,i)
	  values_tvd_remote(:,i) = 0.
	  do l=1,lmax
	    cacum = 0.
	    do ii=1,3
	      k = nen3v(ii,ie)
	      c = values(l,k)
	      cacum = cacum + c * xi(ii)
	    end do
	    values_tvd_remote(l,i) = cacum
	  end do
	end do

	end

!******************************************************************

	subroutine tvd_exchange_single

! exchanges upwind values - one by one - should be used only for debug

	use shympi
	use shympi_tvd

	implicit none

	logical bdebug
	integer i,na,n,nin,nout,iudb
	integer my_ia
	integer id_from,id_to
	integer lmax,lmax_remote

	bdebug = .false.	!leave this at .false.

	my_ia = my_id + 1
	iudb = 350 + my_ia
	nin = 0
	nout = 0

	lmax = nlv_global
	values_tvd_aux = 0.

	do i=1,n_my_array_list
	  na = i
	  id_from = nint(my_array_list(ind_ia_found,na)) - 1
	  id_to = nint(my_array_list(ind_ia_this,na)) - 1

	  n = lmax

	  if( id_from == my_id ) then
	    nout = nout + 1
            tvd_buffer_in(1:n) = values_tvd_remote(1:n,nout)
	  end if

	  if( bdebug ) then
	    write(iudb,*) 'exchanging: ',my_ia,id_from+1,id_to+1,n &
     &			,values_tvd_remote(1,nout)
	  end if

          call shympi_receive(id_from,id_to,n,tvd_buffer_in,tvd_buffer_out)

	  if( id_to /= my_id ) cycle
	  nin = nin + 1
          values_tvd_aux(1:n,nin) = tvd_buffer_out(1:n)

	  if( bdebug ) then
	    write(iudb,*) 'exchanged: ',my_ia,nin,values_tvd_aux(1,nin)
	  end if
	end do

	if( nin /= n_my_receive_list .or. nout /= n_my_send_list ) then
	  write(6,*) my_id,my_ia
	  write(6,*) nin,n_my_receive_list
	  write(6,*) nout,n_my_send_list
	  stop 'error stop tvd_exchange: internal (1)'
	end if

	end

!******************************************************************

	subroutine tvd_exchange_multi

! exchanges upwind values - multi block exchange - more efficient than single

	use shympi
	use shympi_tvd

	implicit none

	logical bdebug,bexchange
	integer i,na,n,nin,nout,iudb,nbuf,ib,id_to_old,nfill,ip
	integer my_ia
	integer id_from,id_to
	integer ia,ia_min,ia_max,ia_tot,iblock
	integer lmax,lmax_remote,l
	integer ipmin,ipmax

	bdebug = .false.	!leave this at .false.

	my_ia = my_id + 1
	iudb = 350 + my_ia
	nin = 0
	nout = 0
	nbuf = 0
	nfill = 0

	lmax = nlv_global
	values_tvd_local = 0.

	do iblock=1,nblocks_multi
	  id_from = exchange_multi_blocks(1,iblock)
	  id_to   = exchange_multi_blocks(2,iblock)
	  ia_min  = exchange_multi_blocks(3,iblock)
	  ia_max  = exchange_multi_blocks(4,iblock)
	  ia_tot  = exchange_multi_blocks(5,iblock)
	  lmax    = exchange_multi_blocks(6,iblock)
	  nbuf = 0
	  do ia=ia_min,ia_max
	    nbuf = nbuf + 1
	    if( id_from == my_id ) then
	      nout = nout + 1
	      ipmax = nbuf * lmax
	      ipmin = ipmax - lmax + 1
              tvd_mbuffer_in(ipmin:ipmax) = values_tvd_remote(1:lmax,nout)
	    end if
	  end do

	  nfill = nbuf * lmax
          call shympi_receive(id_from,id_to &
     &			,nfill,tvd_mbuffer_in,tvd_mbuffer_out)

	  if( id_to /= my_id ) cycle

	  nbuf = nfill / lmax
	  do ib=1,nbuf
	    nin = nin + 1
	    ipmax = ib * lmax
	    ipmin = ipmax - lmax + 1
            values_tvd_local(1:lmax,nin) = tvd_mbuffer_out(ipmin:ipmax)
	  end do
	end do

	if( .not. bcheck_multi ) return

	if( any(values_tvd_local /= values_tvd_aux ) ) then
	  write(6,*) 'nin: ',nin,n_tvd_receive_total,lmax
	  do i=1,nin
	    do l=1,n
	      write(6,*) i,l,values_tvd_aux(l,i),values_tvd_local(l,i),my_id
	    end do
	  end do
	  stop 'error stop tvd_exchange_multi: exchanged values differ in multi'
	end if

	end

!******************************************************************

	subroutine tvd_exchange_prepare_multi

! prepares data structures to exchange multi values

	use shympi
	use shympi_tvd

	implicit none

	logical bnew
	integer iblock
	integer id_from,id_to
	integer id_from_old,id_to_old
	integer i,na
	integer my_ia,iudb
	integer ntot
	integer lmax,lmax_found,lmax_this
	integer, save :: icall = 0

	if( icall > 0 ) return
	icall = icall + 1

	my_ia = my_id + 1
	iudb = 340 + my_ia

	iblock = 0
	id_from_old = -1
	id_to_old = -1

	do i=1,n_my_array_list
	  na = i
	  id_from = nint(my_array_list(ind_ia_found,na)) - 1
	  id_to = nint(my_array_list(ind_ia_this,na)) - 1
	  bnew = .false.
	  if( id_from /= id_from_old ) bnew = .true.
	  if( id_to /= id_to_old ) bnew = .true.
	  if( bnew ) iblock = iblock + 1
	  id_from_old = id_from
	  id_to_old = id_to
	end do

	nblocks_multi = iblock
	allocate(exchange_multi_blocks(6,nblocks_multi))
	exchange_multi_blocks = 0

	iblock = 0
	id_from_old = -1
	id_to_old = -1

	do i=1,n_my_array_list
	  na = i
	  id_from = nint(my_array_list(ind_ia_found,na)) - 1
	  id_to = nint(my_array_list(ind_ia_this,na)) - 1
	  lmax_found = nint(my_array_list(ind_lmax_found,na))
	  lmax_this = nint(my_array_list(ind_lmax_this,na))
	  bnew = .false.
	  if( id_from /= id_from_old ) bnew = .true.
	  if( id_to /= id_to_old ) bnew = .true.
	  if( bnew ) then
	    iblock = iblock + 1
	    exchange_multi_blocks(1,iblock) = id_from
	    exchange_multi_blocks(2,iblock) = id_to
	    exchange_multi_blocks(3,iblock) = i
	  end if
	  exchange_multi_blocks(4,iblock) = i
	  exchange_multi_blocks(5,iblock) = exchange_multi_blocks(5,iblock) + 1
	  lmax = exchange_multi_blocks(6,iblock)
	  lmax = max(lmax,lmax_found,lmax_this)
	  exchange_multi_blocks(6,iblock) = lmax
	  id_from_old = id_from
	  id_to_old = id_to
	end do

	ntot = 0
	do iblock=1,nblocks_multi
	  ntot = ntot + exchange_multi_blocks(5,iblock)
	end do
	if( ntot /= n_my_array_list ) then
	  stop 'error stop tvd_exchange_prepare_multi: ntot (internal error)'
	end if

	if( .not. btvddebug ) return

	write(iudb,*) '-----------------------------'
	write(iudb,*) 'exchange blocks ',my_ia
	write(iudb,*) '-----------------------------'

	do iblock=1,nblocks_multi
	  write(iudb,'(7i8)') iblock,exchange_multi_blocks(:,iblock)
	end do

	!write(6,*) 'total values to be exchanged: ',ntot,n_my_array_list,my_id
	write(iudb,*) 'total values to be exchanged: ',ntot,n_my_array_list

	end

!******************************************************************

	subroutine tvd_get_upwind(l,ipoint,cu)

! gets upwind value stored in array

	use basin
	use levels
	use shympi
	use shympi_tvd

	implicit none

	integer l,ipoint
	real cu

        cu = values_tvd_local(l,ipoint)

	end
	
!******************************************************************
!******************************************************************
!******************************************************************
! from here check and debug routines
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine assert_list_coherence(nmax,list,bconsolidated)

	use basin
	use shympi
	use shympi_tvd

	implicit none

	integer nmax
	real list(nlist,nmax)
	logical bconsolidated

	integer i,j,ii,ie,my_ia
	integer ia_found,ia_needed

	if( .not. btvdassert ) return

	my_ia = my_id + 1

	do i=1,nmax
	  ie = nint(list(ind_ie_this,i))
	  ia_found = nint(list(ind_ia_found,i))
	  ia_needed = nint(list(ind_ia_this,i))
	  if( .not. bconsolidated .and. ie == 0 ) cycle		!end of data
	  j = nint(list(ind_j_this,i))
	  ii = nint(list(ind_ii_this,i))
	  call tvd_assert('j out of bounds',j>=0.and.j<=3)
	  call tvd_assert('ii out of bounds',ii>=0.and.ii<=3)
	  call tvd_assert('ie<=0 ',ie>0.)
	  if( ia_needed == my_ia ) then
	    call tvd_assert('ie out of bounds',ie<=nel)
	  end if
	  call tvd_assert('ia_found out of bounds' &
     &			,ia_found >= 1 .or. ia_found <= n_threads)
	  call tvd_assert('ia_needed out of bounds' &
     &			,ia_needed >= 1 .or. ia_needed <= n_threads)
	end do

	end

!******************************************************************

	subroutine tvd_exchange_debug

! not generally usefull - only use in real debugging situations

	use shympi
	use shympi_tvd

	implicit none

	integer i,na,n,iee
	integer nin,nout
	integer my_ia,iaf
	integer iudb
	integer id_from,id_to
	integer lmax,lmax_remote
	real val,val0

	my_ia = my_id + 1
	nin = 0
	nout = 0
	iudb = 340 + my_ia

	write(iudb,*) '----------------------- starting exchange ',my_id

	do i=1,n_my_array_list
	  na = i
	  id_from = nint(my_array_list(ind_ia_found,na)) - 1
	  id_to = nint(my_array_list(ind_ia_this,na)) - 1
	  lmax = nint(my_array_list(ind_lmax_found,na))
	  iee = nint(my_array_list(ind_iee_found,na))
	  iaf = nint(my_array_list(ind_ia_found,na))
	  val = iaf*iee

	  if( id_from == my_id ) nout = nout + 1
	  n = lmax
          tvd_buffer_in(1:n) = val
          tvd_buffer_out = 0.
	  write(iudb,*) 'exchanging ',id_from,id_to,my_id,iaf,val
          call shympi_receive(id_from,id_to,n,tvd_buffer_in,tvd_buffer_out)
	  if( id_to /= my_id ) cycle
	  nin = nin + 1
          values_tvd_local(1:n,nin) = tvd_buffer_out(1:n)
	end do

	if( nin /= n_my_receive_list .or. nout /= n_my_send_list ) then
	  write(6,*) my_id,my_ia
	  write(6,*) nin,n_my_receive_list
	  write(6,*) nout,n_my_send_list
	  stop 'error stop tvd_exchange_debug: internal (1)'
	end if

	write(iudb,*) 'finished exchange: ',nin,nout
	write(iudb,*) 'debug exchange: ',nin,n_my_receive_list

	do i=1,n_my_receive_list
	  na = my_receive_list(i)
	  id_from = nint(my_array_list(ind_ia_found,na)) - 1
	  id_to = nint(my_array_list(ind_ia_this,na)) - 1
	  lmax = nint(my_array_list(ind_lmax_found,na))
	  iee = nint(my_array_list(ind_iee_found,na))
	  iaf = nint(my_array_list(ind_ia_found,na))
	  val = values_tvd_local(1,i)
	  val0 = iaf*iee
	  write(iudb,*) id_from,id_to,val,val0,my_id
	end do

	end

!******************************************************************
!******************************************************************
!******************************************************************
! debug routines
! set btvddebug = .true. in concentration.f90 and tvd_admin.f90
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine tvd_debug_initialize(dtime,what,isact)

	use basin
	use levels
	use shympi
	use shympi_tvd

	implicit none

	double precision dtime
	character*(*) what
	integer isact
	integer l,lmax,lmax_form,lmax_unf
	integer itime
	integer, save :: icall = 0
	character*80 what_aux

	if( .not. btvddebug ) return
	if( iu_debug <= 0 ) return

	if( icall == 0 ) then
	  allocate( tvd_debug_aux(nlvdi,3,nel) )
	  allocate( blevel(nlv_global) )
	end if
	tvd_debug_aux = 0.

	icall = icall + 1

	if( ifreq_debug < 0 ) then
	  bdebout = .false.
	else if( ifreq_debug > 0 ) then
	  bdebout = .false.
	  if( mod(icall,ifreq_debug) == 0 ) bdebout = .true.
	else
	  bdebout = .true.
	end if

	!bdebout = .false.
	!if( dtime == 1200 .and. what == 'temp' ) bdebout = .true.

	lmax_unf = 0
	lmax_form = 0
	if( bdebout ) then	!we write the output
	  lmax_unf = nlv_global
	  lmax = 0
	  blevel = .false.
	  do l=1,lmax_unf,lmax_unf/2
	    blevel(l) = .true.
	    lmax = lmax + 1
	  end do
	  lmax_form = lmax
	end if

	if( shympi_is_master() ) then
	  if( icall == 0 ) write(iu_debug,*) 'written by tvd_debug_initialize'
	  write(iu_debug,*) 'time = ',dtime,isact,lmax_form,'  ',trim(what)
	  flush(iu_debug)
	  what_aux(1:80) = what
	  itime = nint(dtime)
	  write(iuunf_debug) dtime,isact,lmax_unf,what_aux
	  flush(iuunf_debug)
	end if

	end

!******************************************************************

	subroutine tvd_debug_accum(ie,l,conu)

	use shympi
	use shympi_tvd

	implicit none

	integer ie,l
	real conu(3)

	if( .not. btvddebug ) return
	if( iu_debug <= 0 ) return

	tvd_debug_aux(l,:,ie) = conu(:)

	end

!******************************************************************

	subroutine tvd_debug_finalize

	use basin
	use levels
	use shympi
	use shympi_tvd

	implicit none

	logical bout
	integer l,lmax,ie,ie_ext
	real, allocatable :: aux(:,:),auxg(:,:)
	integer, save :: icall = 0
	integer, save :: nrec = 0
	integer, save :: nrec_form = 0
	integer, save :: nrec_unf = 0

	if( .not. btvddebug ) return
	if( iu_debug <= 0 ) return

	bout = .false.		!if false checks how may records are written
	bout = .true.		!if false checks how may records are written
	icall = icall + 1

	if( .not. bdebout ) return

	allocate(aux(3,nel))
	allocate(auxg(3,nel_global))

	lmax = nlv_global

	do l=1,lmax
	  auxg = 0.
	  aux = 0.
	  if( l <= nlv ) aux(:,:) = tvd_debug_aux(l,:,:)
	  call shympi_l2g_array(3,aux,auxg)
	  if( shympi_is_master() ) then
	    nrec_unf = nrec_unf + 1
	    if( blevel(l) ) then
	      nrec_form = nrec_form + 1
	      write(iu_debug,*) 'level = ',l,nel_global,nrec_form,icall
	    end if
	    write(iuunf_debug) l,nel_global,nrec_unf,icall
	    if( .not. bout ) cycle
	    do ie=1,nel_global
	      ie_ext = ip_ext_elem(ie)
	      if( blevel(l) ) then
	        write(iu_debug,*) ie,ie_ext,auxg(:,ie)
	      end if
	      write(iuunf_debug) ie,ie_ext,auxg(:,ie)
	    end do
	  end if
	end do

	if( shympi_is_master() ) then
	  flush(iu_debug)
	  flush(iuunf_debug)
	end if

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine write_tvd_debug(nel)

	use shympi
	use shympi_tvd
	use mod_tvd

	implicit none

	integer nel

	integer ie,n,j,ii,iu,ie_ext
	integer, allocatable :: iedebug(:,:)
	integer, allocatable :: ieout(:,:)
	integer, allocatable :: iadebug(:,:)
	integer, allocatable :: iaout(:,:)
	character*80 file

	if( .not. btvddebug ) return

	iu = 583
	file='tvd_debug.txt'

	call shympi_barrier
	call flush(6)

	write(6,*) 'collecting tvd debug information...',my_id

	allocate(iedebug(9,nel))
	allocate(ieout(9,nel_global))
	allocate(iadebug(9,nel))
	allocate(iaout(9,nel_global))
	iedebug = 0
	iadebug = 0
	ieout = 0
	iaout = 0

	do ie=1,nel
	  n = 0
	  do ii=1,3
	    do j=1,3
	      n = n + 1
	      iedebug(n,ie) = ieetvdup(j,ii,ie)
	    end do
	  end do
	  iadebug(:,ie) = reshape(ieetvdup(:,:,ie),(/9/))
	  if( any( iedebug(:,ie) /= iadebug(:,ie) ) ) then
	    stop 'error stop write_tvd_debug: reshape not working'
	  end if
	  iadebug(:,ie) = reshape(iatvdup(:,:,ie),(/9/))
	end do

	call shympi_l2g_array(9,iedebug,ieout)
	call shympi_l2g_array(9,iadebug,iaout)

	write(6,*) 'finished collecting tvd debug information...',my_id

	if( shympi_is_master() ) then

	open(iu,file=file,status='unknown',form='formatted')
	write(iu,*) 'debug information written by write_tvd_debug'
	do ie=1,nel_global
	  ie_ext = ip_ext_elem(ie)
	  write(iu,'(11i7)') ie,ie_ext,ieout(:,ie)
	  !write(iu,'(9i7)') iaout(:,ie)
	end do
	close(iu)
	write(6,*) 'tvd debug information written to file ',trim(file)

	end if

	write(6,*) 'finished tvd debug...',my_id
	call shympi_barrier
	flush(6)

	end

!******************************************************************

	subroutine tvd_assert(text,bassert)

	implicit none

	character*(*) text
	logical bassert

	if( bassert ) return

	write(6,*) 'assertion violated: ',trim(text)
	stop 'error stop tvd_assert: assertion violated'

	end

!******************************************************************

	subroutine write_array_debug(iudb,i,array)

	use shympi_tvd

	implicit none

	integer iudb,i
	real array(nlist)

	if( nlist_type == 1 ) then
	  write(iudb,'(10i6)') i,nint(array(3:10))
	else if( nlist_type == 2 ) then
	  write(iudb,'(i4,11i6)') i,nint(array(3:nlist))
	end if

	end

!******************************************************************

	subroutine write_array_header(iudb)

	use shympi_tvd

	implicit none

	integer iudb

	if( nlist_type == 1 ) then
	  write(iudb,*) '    i   ian     j    ii    ie   iee   iaf   ief  ieef'
	else if( nlist_type == 2 ) then
	  write(iudb,*) '  i   ian     j    ii     k    ie   iee  lmax' &
     &		//	'   iaf   ief  ieef lmaxf'
	end if

	end

!******************************************************************
!******************************************************************
!******************************************************************
! test routines
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine test_tvd_convert

	use shympi_tvd

	implicit none

	integer k,j,ii,m

	write(6,*) 'testing tvd convert routines...'

	do k=1,6
	  call k2jii(k,j,ii)
	  call jii2k(j,ii,m)
	  if( k /= m ) stop 'error stop text_tvd_convert: k/=m'
	  write(6,*) j,ii,k,m
	end do

	end

!******************************************************************
!	program main_test_tvd
!	call test_tvd_convert
!	end
!******************************************************************

