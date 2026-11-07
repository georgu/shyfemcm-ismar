
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

!
! connectivity routines
!
! revision log :
!
! 19.05.2020	ccf	started from scratch
! 12.04.2022	ggu	adapted
! 07.11.2024	ggu	is ready now (completely auto-sufficient)
!
!****************************************************************

!==================================================================
	module mod_connect
!==================================================================

	implicit none

	private

	integer, save :: nkn = 0
	integer, save :: nel = 0
	integer, save :: ngr = 0

	integer, save, allocatable :: nen3v(:,:)

	integer, save, allocatable :: nlist(:,:)
	integer, save, allocatable :: elist(:,:)
	integer, save, allocatable :: ecv(:,:)
	integer, save, allocatable :: bound(:)
	integer, save, allocatable :: kant(:,:)

	integer, save :: ierrors = 0
	integer, save, allocatable :: kerrors(:)
	integer, save, allocatable :: iperrors(:)

	public :: connect_init
	public :: connect_get_grade
	public :: connect_get_grades
	public :: connect_get_lists
	public :: connect_get_ecv
	public :: connect_get_kant
	public :: connect_get_bound
	public :: connect_errors
	public :: connect_release

! nk = nlist(0,k) is total number of neighboring nodes of node k
! nlist(1:nk,k) is list of neighboring nodes of node k
! ek = elist(0,k) is total number of connected elems of node k
! elist(1:ek,k) is list of connected elems of node k
! ecv(3,e) is list of elements connected to elem e
! bound(k) is indicator if k is boundary node
! kant(:,k) are neigboring boundary nodes of k if k is boundary node

! calling sequence:
!
!	call connect_init(nkn,nel,nen3v,ierr)
!
!	call connect_get_grade(ngr)
!	call connect_get_grades(nkn,ngrade,egrade,ngr)
!	call connect_get_lists(nkn,ngr,nlist,elist)
!       call connect_get_ecv(nel,ecv)
!       call connect_get_kant(nkn,kant)
!       call connect_get_bound(nkn,bound)
!
!	call connect_errors(ierr,kerr)
!	call connect_release

!==================================================================
	contains
!==================================================================

	subroutine connect_internal_allocate(nkn_l,nel_l,ngr_l,nen3v_l)

	integer nkn_l,nel_l,ngr_l
	integer nen3v_l(3,nel_l)

	logical brealloc

	brealloc = .true.
	if( ngr > 0 ) then
	  if( nkn == nkn_l .and. nel == nel_l ) brealloc = .false.
	end if
	brealloc = .true.

	if( .not. brealloc ) return

	call connect_internal_deallocate

	nkn = nkn_l
	nel = nel_l
	ngr = ngr_l

	allocate(nen3v(3,nel))
	allocate(nlist(0:ngr,nkn))
	allocate(elist(0:ngr,nkn))
	allocate(bound(nkn))
	allocate(ecv(3,nel))
	allocate(kant(2,nkn))

	ierrors = 0
	allocate(kerrors(nkn))
	allocate(iperrors(nkn))
	kerrors = 0
	iperrors = 0

	nen3v = nen3v_l

	end

!******************************************************************

	subroutine connect_internal_deallocate

	if( ngr == 0 ) return

	ngr = 0
	nkn = 0
	nel = 0

	deallocate(nen3v)
	deallocate(nlist)
	deallocate(elist)
	deallocate(bound)
	deallocate(ecv)
	deallocate(kant)

	ierrors = 0
	deallocate(kerrors)
	deallocate(iperrors)

	end

!******************************************************************

	subroutine connect_internal_init(nkn_l,nel_l,nen3v_l,ierr)

	integer nkn_l,nel_l
	integer nen3v_l(3,nel_l)
	integer ierr

	integer ngr_l

	if( ngr > 0 ) call connect_internal_deallocate

	call make_grade(nkn_l,nel_l,nen3v_l,ngr_l,ierr)
	if( ierr /= 0 ) return
	call connect_internal_allocate(nkn_l,nel_l,ngr_l,nen3v_l)

	call make_ne_list
	call make_errors
	call sort_ne_lists
	call make_bound
	call make_kant
	call make_ecv

	end

!******************************************************************

	subroutine connect_internal_parameters(nkn_l,nel_l,ngr_l)

	integer nkn_l,nel_l,ngr_l

	nkn_l = nkn
	nel_l = nel
	ngr_l = ngr

	end

!******************************************************************

	subroutine connect_internal_check(text,nkn_l,nel_l,ngr_l)

	character*(*) text
	integer nkn_l,nel_l,ngr_l

	logical berror

	if( ngr == 0 ) then
	  write(6,*) trim(text)
	  stop 'error stop connect_internal_check: no initialization'
	end if

	berror = .false.

	if( nkn_l > 0 .and. nkn_l /= nkn ) berror = .true.
	if( nel_l > 0 .and. nel_l /= nel ) berror = .true.
	if( ngr_l > 0 .and. ngr_l /= ngr ) berror = .true.

	if( berror ) then
	  write(6,*) trim(text)
	  write(6,*) 'nkn: ',nkn_l,nkn
	  write(6,*) 'nel: ',nel_l,nel
	  write(6,*) 'ngr: ',ngr_l,ngr
	  stop 'error stop connect_internal_check: incompatible parameters'
	end if

	end

!******************************************************************

	subroutine connect_get_element_index(nen3v_l)

	integer nen3v_l(3,nel)

	nen3v_l = nen3v

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine make_errors

	implicit none

	integer k

	ierrors = 0

	do k=1,nkn
	  if( nlist(0,k) - elist(0,k) > 1 ) then
	    ierrors = ierrors + 1
	    kerrors(ierrors) = k
	    iperrors(k) = 1
	  end if
	end do

	end

!******************************************************************

	subroutine make_grade(nkn,nel,nen3v,ngr,ierr)

! makes grade ngr - also computes ngrade and egrade

	implicit none

	integer nkn,nel
	integer nen3v(3,nel)
	integer ngr		!maximum grade (return)
	integer ierr

	integer ie,ii,k,ii1,k1,ngre
	integer idiff,nbnd
	integer, allocatable :: ngrade(:)
	integer, allocatable :: egrade(:)
	integer, allocatable :: nalist(:,:)

	ierr = 0

	allocate(ngrade(nkn))
	allocate(egrade(nkn))

	egrade = 0

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    egrade(k) = egrade(k) + 1
	  end do
	end do

	ngr = 1 + maxval(egrade)
	allocate(nalist(2*ngr,nkn))
	!write(6,*) 'ngr estimate: ',ngr
	
	ngrade = 0
	nalist = 0

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    ii1 = 1 + mod(ii,3)
	    k1 = nen3v(ii1,ie)
	    call insert_in_list(k,k1)
	    call insert_in_list(k1,k)
	  end do
	end do

	ngr = maxval(ngrade)
	ngre = maxval(egrade)

	nbnd = 0
	do k=1,nkn
	  idiff = abs( ngrade(k) - egrade(k) )
	  if( idiff == 1 ) then
	    nbnd = nbnd + 1
	  else if( idiff /= 0 ) then
	    write(6,*) 'error in grade: ',k,ngrade(k),egrade(k)
	    ierr = ierr + 1
	  end if
	end do

	contains

	subroutine insert_in_list(k1,k2)
	integer k1,k2
	integer n,i
	n = ngrade(k1)
	if( any( nalist(1:n,k1) == k2 ) ) return
	n = n + 1
	if( n > 2*ngr ) stop 'error stop make_grade: ngr'
	nalist(n,k1) = k2
	ngrade(k1) = n
	end subroutine insert_in_list

	end

!******************************************************************

	subroutine make_ne_list

	implicit none

	integer ie,ii,ii1,k,k1

	nlist = 0
	elist = 0

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    ii1 = 1 + mod(ii,3)
	    k1 = nen3v(ii1,ie)
	    call insert_item_in_list(k,k1,nlist)
	    call insert_item_in_list(k1,k,nlist)
	    call insert_item_in_list(k,ie,elist)
	  end do
	end do

	end

!*******************************************************************

	subroutine sort_ne_lists

! sorts element and node list
!
! start from boundary node if exists
! handles multiple border nodes (inserting 0 elements)

	implicit none

	logical bbound,bextra
	integer ie,k,kn,kb,iee,ii
	integer i,j,nn,ne
	integer ien,ipe,ipk,naux,kk
	integer in,ib
	integer el(ngr),elaux(ngr)
	integer nl(ngr),nlaux(ngr)
	character*80 text

!------------------------------------------------------------
! first pass - sort elements
!
! deals with multi border nodes
! inserts zero elements in element list if elements are not connected
!------------------------------------------------------------

	do k=1,nkn
	  ne = elist(0,k)
	  nn = nlist(0,k)
	  el(1:ne) = elist(1:ne,k)
	  ipe = 0
	  elaux = 0
	  bbound = ( nn /= ne )
	  bextra = ( nn > ne+1 )
	  ie = el(1)
	  if( bbound ) then
	    call look_for_first_element(k,ie)	!finds first element
	    call exchange_item(k,ie,1,el)
	  end if
	  ipe = ipe + 1
	  elaux(ipe) = ie
	  do i=2,ne
	    ie = el(i-1)
	    ien = enext(ie,k)
	    if( ien == 0 ) then
	      ipe = ipe + 1
	      elaux(ipe) = 0			!insert zero element
	      ien = el(i)
	      call look_for_first_element(k,ien)!finds first element
	    end if
	    call exchange_item(k,ien,i,el)
	    ipe = ipe + 1
	    elaux(ipe) = ien
	  end do
	  if( bbound ) then			!ipe is filling of element list
	    if( ipe+1 /= nn ) goto 91
	  else
	    if( ipe /= nn ) goto 91
	  end if
	  ne = ipe
	  elist(0,k) = ne
	  elist(1:ne,k) = elaux(1:ne)
	  if( bextra ) then
	    write(6,*) 'node inconsistency (elements): ',k
	    write(6,*) ne,nn
	    write(6,*) elist(1:ne,k)
	  end if
	end do

!------------------------------------------------------------
! second pass - sort nodes
!
! nodes cannot have zero values
!------------------------------------------------------------

	do k=1,nkn
	  ne = elist(0,k)
	  nn = nlist(0,k)
	  el(1:ne) = elist(1:ne,k)
	  nl(1:nn) = nlist(1:nn,k)
	  ipk = 0
	  bbound = ( nn /= ne )
	  bextra = .false.
	  do i=1,ne
	    ie = el(i)
	    if( ie == 0 ) then
	      bextra = .true.
	      kk = kbhnd(el(i-1),k)
	      ipk = ipk + 1
	      nlaux(ipk) = kk
	      call exchange_item(k,kk,ipk,nl)
	      cycle
	    end if
	    kk = knext(ie,k)
	    ipk = ipk + 1
	    nlaux(ipk) = kk
	    call exchange_item(k,kk,ipk,nl)
	  end do
	  if( bbound ) then		!add last node
	    ie = el(ne)
	    kk = kbhnd(ie,k)
	    ipk = ipk + 1
	    nlaux(ipk) = kk
	    call exchange_item(k,kk,ipk,nl)
	  end if
	  if( ipk /= nn ) goto 92
	  nlist(1:nn,k) = nlaux(1:nn)
	  if( bextra ) then
	    write(6,*) 'node inconsistency (nodes): ',k
	    write(6,*) ne,nn
	    write(6,*) nlist(1:nn,k)
	    write(6,*) elist(1:ne,k)
	    write(6,*) 'continue processing...'
	  end if
	end do

!------------------------------------------------------------
! third pass - check
!------------------------------------------------------------

	text = 'checking lists'

	do k=1,nkn
	  ne = elist(0,k)
	  nn = nlist(0,k)
	  el(1:ne) = elist(1:ne,k)
	  nl(1:nn) = nlist(1:nn,k)
	  bbound = ( nn /= ne )
	  !naux = count(el(1:ne)==0)
	  naux = ne
	  if( bbound ) naux = naux + 1
	  !if( bbound .and. ne+naux+1 /= nn ) goto 97
	  !if( bbound .and. ne+1 /= nn ) goto 97
	  if( naux /= nn ) goto 97
	  if( any(nl(1:nn)==0) ) goto 97
	  do i=1,ne-1
	    ie = el(i)
	    iee = el(i+1)
	    if( ie == 0 .or. iee == 0 ) cycle
	    kb = kbhnd(ie,k)
	    kn = knext(iee,k)
	    if( kn /= kb ) goto 98
	  end do
	  if( bbound ) cycle
	  kn = knext(el(1),k)
	  kb = kbhnd(el(ne),k)
	  if( kn /= kb ) goto 98
	end do

!------------------------------------------------------------
! end of routine
!------------------------------------------------------------

	return
   91	continue
	write(6,*) 'inconsistency in element list: ',k,nn,ne
	write(6,*) bbound,ipe
	write(6,*) elaux(1:ne)
	write(6,*) elist(0,k),nlist(0,k)
	write(6,*) elist(1:ne,k)
	write(6,*) nlist(1:nn,k)
	stop 'error stop sort_ne_lists: inconsistency'
   92	continue
	write(6,*) 'inconsistency in node list: ',k,nn,ne
	write(6,*) nlaux(1:ne)
	stop 'error stop sort_ne_lists: inconsistency'
   97	continue
	write(6,*) 'phase: ',trim(text)
	write(6,*) 'zeros in list: ',k
	write(6,*) nn,ne,naux
	call info('checking nodes',k)
	stop 'error stop sort_ne_lists: zeros'
   98	continue
	write(6,*) 'phase: ',trim(text)
	write(6,*) 'k,i: ',k,i
	write(6,*) 'kn,kb not consistent: ',ie,iee,kn,kb
	call info('checking nodes',k)
	stop 'error stop sort_ne_lists: inconsistency'
   99	continue
	write(6,*) 'phase: ',trim(text)
	write(6,*) 'error in list'
	write(6,*)  'ie,iee: ',ie,iee
	call info('on error',k)
	stop 'error stop sort_ne_lists: generic'

	contains

	subroutine look_for_first_element(k,ie)
	integer k,ie
	integer iee
	do
	  iee = ebhnd(ie,k)
	  if( iee == 0 ) exit
	  ie = iee
	end do
	end

	subroutine insert_item(k,ie,j)
	integer k,ie,j
	integer kn
	kn = knext(ie,k)
	call exchange_item(k,ie,j,el)
	call exchange_item(k,kn,j,nl)
	end

	end

!*******************************************************************

	subroutine make_bound

! makes bound - indicator of boundary nodes

	implicit none

	integer k,n

	bound = 0

	do k=1,nkn
	  n = nlist(0,k) - elist(0,k)
	  if( n == 1 ) then
	  !if( n > 0 ) then
	    bound(k) = 1
	  else if( n /= 0 ) then
	    write(6,*) k,n
	    stop 'error stop make_bound: error in connectivity'
	  end if
	end do

	end

!*******************************************************************

	subroutine make_kant

! makes kant - neighboring boundary nodes

	implicit none

	integer k,nk

	kant = 0

	do k=1,nkn
	  if( bound(k) == 0 ) cycle	!bound must already be set up
	  if( iperrors(k) /= 0 ) cycle	!cannot handle irregular boundary node
	  nk = nlist(0,k)
	  kant(1,k) = nlist(1,k)
	  kant(2,k) = nlist(nk,k)
	end do

	end

!*******************************************************************

	subroutine make_ecv

! makes neighbor information of elements

	implicit none

	logical bbound
	integer k,ne,nn,i
	integer ie,ii,iee
	integer el(ngr)

	ecv = 0

	do k=1,nkn
	  ne = elist(0,k)
	  nn = nlist(0,k)
	  el(1:ne) = elist(1:ne,k)
	  do i=1,ne-1
	    call ecv_insert(el(i),el(i+1))
	  end do
	  if( nn == ne ) call ecv_insert(el(ne),el(1))
	end do

! check neighbor list

	do ie=1,nel
	  do ii=1,3
	    iee = ecv(ii,ie)
	    bbound = is_edge_on_boundary(ii,ie)
	    if( iee == 0 ) then
	      if( bbound ) cycle	!edge on boundary... ok
	    else if( any(ecv(:,iee)==ie) ) then
	      cycle				!ie in neighbor list... ok
	    end if
	    write(6,*) ie,ii,iee,bbound
	    write(6,*) ie,ecv(:,ie)
	    if( iee > 0 ) write(6,*) iee,ecv(:,iee)
	    stop 'error stop make_ecv: checking ecv'
	  end do
	end do

	contains

	subroutine get_edge_nodes(ii,ie,k1,k2)
	integer ii,ie,k1,k2
	integer ii1,ii2
	ii1 = 1+mod(ii,3)
	ii2 = 1+mod(ii1,3)
	k1 = nen3v(ii1,ie)
	k2 = nen3v(ii2,ie)
	end

	function is_edge_on_boundary(ii,ie)
	integer ii,ie
	logical is_edge_on_boundary
	integer ii1,ii2,k1,k2
	call get_edge_nodes(ii,ie,k1,k2)
	is_edge_on_boundary = .false.
	if( elist(0,k1) == nlist(0,k1) ) return
	if( elist(0,k2) == nlist(0,k2) ) return
	is_edge_on_boundary = .true.
	end

	subroutine ecv_insert(ie1,ie2)
	integer ie1,ie2
	integer ib,in
	if( ie1 == 0 .or. ie2 == 0 ) return
	ib = ibhnd(ie2,k)
	in = inext(ie1,k)
	ecv(in,ie1) = ie2
	ecv(ib,ie2) = ie1
	end

	end

!*******************************************************************
!******************************************************************
!******************************************************************
!******************************************************************
! local private subroutines
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine info(text,k)
	character*(*) text
	integer k
	integer i,ie,nn,ne
	nn = nlist(0,k)
	ne = elist(0,k)
	write(6,*) '----------------------'
	write(6,*) trim(text)
	write(6,*)  'k',k
	write(6,*)  'nn,ne',nn,ne
	write(6,*)  'norig',nlist(1:nn,k)
	write(6,*)  'eorig',elist(1:ne,k)
	!write(6,*)  'nl',nl(1:nn)
	!write(6,*)  'el',el(1:ne)
	write(6,*)  'element list:'
	do i=1,ne
	  ie = elist(i,k)
	  write(6,*) ie,nen3v(:,ie)
	end do
	write(6,*) '----------------------'
	end

	function ebhnd(ie,k)
	integer ebhnd
	integer ie,k
	integer i,iee,kn,ne
	kn = knext(ie,k)
	ne = elist(0,k)
	do i=1,ne
	  iee = elist(i,k)
	  if( kbhnd(iee,k) == kn ) exit
	end do
	if( i > ne ) iee = 0
	ebhnd = iee
	end

	function enext(ie,k)
	integer enext
	integer ie,k
	integer i,iee,kn,ne
	kn = kbhnd(ie,k)
	ne = elist(0,k)
	do i=1,ne
	  iee = elist(i,k)
	  if( knext(iee,k) == kn ) exit
	end do
	if( i > ne ) iee = 0
	enext = iee
	end

	function kbhnd(ie,k)
	integer kbhnd
	integer ie,k
	integer ii
	if( ie == 0 ) stop 'error stop kbhnd: ie==0'
	do ii=1,3
	  if( nen3v(ii,ie) == k ) exit
	end do
	if( ii > 3 ) stop 'error stop kbhnd: node not found'
	ii = 1 + mod(ii+1,3)
	kbhnd = nen3v(ii,ie)
	end

	function knext(ie,k)
	integer knext
	integer ie,k
	integer ii
	if( ie == 0 ) stop 'error stop knext: ie==0'
	do ii=1,3
	  if( nen3v(ii,ie) == k ) exit
	end do
	if( ii > 3 ) stop 'error stop knext: node not found'
	ii = 1 + mod(ii,3)
	knext = nen3v(ii,ie)
	end
	
	function ibhnd(ie,k)
	integer ibhnd
	integer ie,k
	integer ii
	if( ie == 0 ) stop 'error stop ibhnd: ie==0'
	do ii=1,3
	  if( nen3v(ii,ie) == k ) exit
	end do
	if( ii > 3 ) stop 'error stop ibhnd: node not found'
	ii = 1 + mod(ii+1,3)
	ibhnd = ii
	end

	function inext(ie,k)
	integer inext
	integer ie,k
	integer ii
	if( ie == 0 ) stop 'error stop inext: ie==0'
	do ii=1,3
	  if( nen3v(ii,ie) == k ) exit
	end do
	if( ii > 3 ) stop 'error stop inext: node not found'
	ii = 1 + mod(ii,3)
	inext = ii
	end

	subroutine exchange_item(k,item,jp,list)
	integer k,item,jp
	integer list(ngr)
	integer j,iaux
	j = locate_in_array(item,list)
	if( j == 0 ) then
	  write(6,*) 'looking for: ',item
	  call info('looking for item in exchange_item',k)
	  stop 'error stop exchange_item: cannot find item'
	end if
	if( j == jp ) return
	iaux = list(j)
	list(j) = list(jp)
	list(jp) = iaux
	end

	function locate_in_array(value,array)
	integer value,array(:)
	integer locate_in_array
	integer n,i
	n = size(array)
	do i=1,n
	  if( array(i) == value ) exit
	end do
	locate_in_array = i
	if( i > n ) locate_in_array = 0
	end

	subroutine insert_item_in_list(k,item,list)
	integer k,item
	integer list(0:ngr,nkn)
	integer n,i
	n = list(0,k)
	if(k==0) write(6,1000) 1,item,n,list(1:n,k)
	do i=1,n
	  if( list(i,k) == item ) return
	end do
	n = n + 1
	if( n > ngr ) then
	  write(6,*) n,ngr,k,item
	  write(6,*) list(0:ngr,k)
	  stop 'error stop insert_item_in_list: n>ngr'
	end if
	list(n,k) = item
	list(0,k) = n
	if(k==0) write(6,1000) 3,item,n,list(1:n,k)
 1000	format(15i5)
	end subroutine insert_item_in_list

!*******************************************************************
!*******************************************************************
!*******************************************************************
! public routines - to be called from outside
!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine connect_init(nkn_l,nel_l,nen3v_l,ierr)

	implicit none

	integer, intent(in) :: nkn_l,nel_l
	integer, intent(in) :: nen3v_l(3,nel_l)
	integer, intent(out) :: ierr

	call connect_internal_init(nkn_l,nel_l,nen3v_l,ierr)

	end

!*******************************************************************

	subroutine connect_release

	implicit none

	call connect_internal_deallocate

	end

!*******************************************************************

	subroutine connect_errors(ierr,kerr)

! fills kerr with error nodes (up to ierr)

	implicit none

	integer ierr
	integer kerr(ierr)

	integer i,k,ndim

	ndim = ierr

	if( ierrors == 0 ) return

	write(6,*) 'connection errors found: ',ierrors
	do i=1,ierrors
	  k = kerrors(i)
	  write(6,*) 'connection error in node ',k
	  if( i <= ndim ) kerr(i) = k
	end do

	if( ndim > 0 .and. ndim < ierrors ) then
	  write(6,*) 'a maximum of ',ndim,'errors have been returned'
	end if

	if( ndim > 0 ) ierr = min(ierrors,ndim)

	end

!*******************************************************************

	subroutine connect_get_grade(ngr_l)

	implicit none

	integer, intent(out) :: ngr_l

	ngr_l = ngr

	end

!*******************************************************************

	subroutine connect_get_grades(nkn_l,ngrade,egrade,ngr_l)

	implicit none

	integer, intent(in) :: nkn_l
	integer, intent(out) :: ngrade(nkn_l)
	integer, intent(out) :: egrade(nkn_l)
	integer, intent(out) :: ngr_l

	call connect_internal_check('get_grades',nkn_l,0,0)

	ngrade(:) = nlist(0,:)
	egrade(:) = elist(0,:)
	ngr_l = ngr

	end

!*******************************************************************

	subroutine connect_get_lists(nkn_l,ngr_l,nlist_l,elist_l)

	implicit none

	integer, intent(in) :: nkn_l
	integer, intent(in) :: ngr_l
	integer, intent(out) :: nlist_l(0:ngr_l,nkn_l)
	integer, intent(out) :: elist_l(0:ngr_l,nkn_l)

	call connect_internal_check('get_lists',nkn_l,0,ngr_l)

	nlist_l = nlist
	elist_l = elist

	end

!*******************************************************************

	subroutine connect_get_bound(nkn_l,bound_l)

	implicit none

	integer, intent(in) :: nkn_l
	integer, intent(out) :: bound_l(nkn_l)

	call connect_internal_check('get_bound',nkn_l,0,0)

	bound_l = bound

	end

!*******************************************************************

	subroutine connect_get_ecv(nel_l,ecv_l)

	implicit none

	integer, intent(in) :: nel_l
	integer, intent(out) :: ecv_l(3,nel_l)

	call connect_internal_check('get_ecv',0,nel_l,0)

	ecv_l = ecv

	end

!*******************************************************************

	subroutine connect_get_kant(nkn_l,kant_l)

	implicit none

	integer, intent(in) :: nkn_l
	integer, intent(out) :: kant_l(2,nkn_l)

	call connect_internal_check('get_kant',nkn_l,0,0)

	kant_l = kant

	end

!*******************************************************************

!==================================================================
	end module mod_connect
!==================================================================

