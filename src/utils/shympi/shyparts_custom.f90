
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

! revision log :
!
! 19.05.2020	ccf	started from scratch
! 12.04.2022	ggu	adapted
! 29.10.2023	ggu	ready for custom partitioning
! 10.04.2024	ggu	restructured for sdda algorithm
!
!****************************************************************

        subroutine do_partition(nkn,nel,nen3v,nparts,npart,epart)

! shyparts dummy routine

        implicit none

        integer nkn,nel
        integer nen3v(3,nel)
        integer nparts
        integer npart(nkn)
        integer epart(nel)

	real pquality

	npart = 0
	epart = 0
	write(6,*) 'using SDDA algorithm for partitioning'
	write(6,*) 'running do_custom with np = ',nparts

        call do_sdda(nkn,nel,nen3v,nparts,npart,epart)

	!call info_partition(nparts,npart,pquality)

	end

!*******************************************************************

	subroutine check_partition(npart,epart,bdebug,ierr1,ierr2)

	use basin

	implicit none

        integer npart(nkn)
        integer epart(nel)
	logical bdebug
	integer ierr1,ierr2

	ierr1 = 0
	ierr2 = 0

	end

!*******************************************************************

        subroutine do_sdda(nkn,nel,nen3v,nparts,npart,epart)

! shyparts custom routine

	use mod_color

        implicit none

        integer nkn,nel
        integer nen3v(3,nel)
        integer nparts
        integer npart(nkn)
        integer epart(nel)

	logical bdebug
	integer nroots
	integer ic,icnew,nic
	integer np,ngr,ngrmin,k,ncol,nmax,knew,i,itype
	integer kroots(nparts+1)
	integer, allocatable :: ngrade(:),egrade(:)
	integer, allocatable :: ncolor(:)
	integer, allocatable :: matrix(:,:)
	integer, allocatable :: dist(:,:)

	bdebug = .false.
	np = nparts

	allocate(ngrade(nkn),egrade(nkn))
	allocate(ncolor(nkn))
	allocate(dist(nkn,0:nparts))

!------------------------------------------------------------
! initialize connectivity routines
!------------------------------------------------------------

	call connect_init(nkn,nel,nen3v)
	call connect_get_grades(nkn,ngrade,egrade,ngr)

	ic = 1
	icnew = 2
	npart = ic
	epart = 0
	call divide_domain(ic,icnew,nkn,npart)
	call make_epart(nkn,nel,nen3v,npart,epart)
	call exchange_nodes_between_two(ic,icnew,nkn,npart,nel,nen3v)
!	call write_partition_to_grd('domain_divide',bdebug &
!     &		,np,npart,epart)

	return

!------------------------------------------------------------
! find first node with lowest grade
!------------------------------------------------------------

	ngrmin = minval(ngrade)
	do k=1,nkn
	  if( ngrade(k) == ngrmin ) exit
	end do
	if( k > nkn ) stop 'error stop do_sdda: k>nkn'

	nroots = 1
	kroots(nroots) = k
	npart = 0
	epart = 0

!------------------------------------------------------------
! find other source nodes
!------------------------------------------------------------

	write(6,*) 'start with root node ',kroots(nroots)
	do nroots=1,nparts-1
	  dist(:,nroots) = -1
	  call compute_distance_from_root(kroots(nroots),nkn,dist(:,nroots))
	  call compute_total_distance(nroots,nkn,dist)
	  call choose_new_source_node(nroots,nkn,dist,knew)
	  kroots(nroots+1) = knew
	  np = nroots
	  npart = dist(:,0)
	  call write_partition_to_grd('domain_dist',bdebug &
     &		,np,npart,epart)
	end do

	  np = nparts
	  npart = dist(:,0)
	  call write_partition_to_grd('domain_dist',bdebug &
     &		,np,npart,epart)

	write(6,*) 'source nodes found: ',nroots

	open(100,file='rootnodes.grd',status='unknown',form='formatted')
	do i=1,nroots
	  k = kroots(i)
	  write(6,'(i3,i6,20i4)') i,k,dist(k,:)
	end do

	itype = 7
	call write_single_nodes('rootnodes.grd',nroots,kroots,itype)

	npart = -1
	do i=1,nparts
	  k = kroots(i)
	  npart(k) = i
	end do

	call fill_domain(nkn,npart)
	call make_epart(nkn,nel,nen3v,npart,epart)
	np = nparts
	call write_partition_to_grd('domain_custom',bdebug &
     &		,np,npart,epart)
!	call write_partition_to_grd('domain_color',bdebug &
!     &		,np,ncolor,epart)

        call exchange_nodes(np,nkn,npart,nel,nen3v)

	return

!------------------------------------------------------------
! loop over partitioning
!------------------------------------------------------------

	do np=1,nparts
	  write(6,*) 'calling ffill ',np
	  call ffill(nkn,nel,nen3v,ngrade,np,kroots,npart)
	  call make_epart(nkn,nel,nen3v,npart,epart)
	  call check_connected_domains(nkn,nel,npart,epart)
	  call make_node_matrix(nel,nen3v,nkn,npart,nmax,matrix)
	  call color_graph(nkn,npart,nmax,matrix,ncol,ncolor)
	  !call color_nodes(nel,nen3v,nkn,npart,ncol,ncolor)
	  call release_color_matrix(matrix)
	  write(6,*) 'ncol = ',ncol
	  write(6,*) 'npart min/max: ',minval(npart),maxval(npart)
	  write(6,*) 'epart min/max: ',minval(epart),maxval(epart)
	  call write_partition_to_grd('domain_custom',bdebug &
     &		,np,npart,epart)
	  call write_partition_to_grd('domain_color',bdebug &
     &		,np,ncolor,epart)
	end do

!------------------------------------------------------------
! end of routine
!------------------------------------------------------------

	end

!*******************************************************************

	subroutine fill_domain(nkn,npart)

	use mod_connect

	implicit none

	integer nkn
	integer npart(nkn)

	integer i,k,iloop,ncol,ic,kk,nk
	integer, allocatable :: newcolor(:)

	allocate(newcolor(nkn))

	iloop = 0
	do
	  iloop = iloop + 1
	  ncol = 0
	  newcolor = npart
	  do k=1,nkn
	    ic = npart(k)
	    if( ic < 0 ) cycle
	    nk = nlist(0,k)
	    do i=1,nk
	      kk = nlist(i,k)
	      if( npart(kk) == -1 ) then
	        newcolor(kk) = ic
	        ncol = ncol + 1
	      end if
	    end do
	  end do
	  npart = newcolor
	  !write(6,*) 'in fill loop: ',iloop,ncol
	  if( ncol == 0 ) exit
	end do

	write(6,*) 'fill loop: ',iloop

	end

!*******************************************************************

	subroutine compute_distance_from_root(kroot,nkn,dist)

! compute distance from kroot on all nodes that have dist(k) == -1

	use mod_connect

	implicit none

	integer kroot           ! source node to start from
	integer nkn		! total nodes in domain
	integer dist(nkn)	! distance from node kroot

	integer idist,k,nk,i,kk,ndist,maxdist
	integer, allocatable :: newdist(:)

	idist = 0
	dist(kroot) = 0
	allocate(newdist(nkn))

	do
	  ndist = 0
	  idist = idist + 1
	  newdist = dist
	  do k=1,nkn
	    if( dist(k) < 0 ) cycle
	    nk = nlist(0,k)
	    do i=1,nk
	      kk = nlist(i,k)
	      if( dist(kk) == -1 ) then
	        newdist(kk) = idist
	        ndist = ndist + 1
	      end if
	    end do
	  end do
	  dist = newdist
	  !write(6,*) 'in loop: ',idist,ndist
	  if( ndist == 0 ) exit
	end do

	maxdist = maxval(dist(:))
	write(6,*) 'maximum dist: ',kroot,idist,maxdist

	end

!*******************************************************************

	subroutine compute_total_distance(nroots,nkn,dist)

	implicit none

	integer nroots
	integer nkn
	integer dist(nkn,0:nroots)

	integer k,maxdist

	do k=1,nkn
	  dist(k,0) = sum( dist(k,1:nroots) )
	end do

	maxdist = maxval(dist(:,0))
	write(6,*) 'maximum total distance: ',nroots,maxdist

	end

!*******************************************************************

	subroutine choose_new_source_node(nroots,nkn,dist,knew)

	use mod_sort

	implicit none

	integer nroots
	integer nkn
	integer dist(nkn,0:nroots)
	integer knew

	logical bsortdist,bsortrange
	integer mindist,maxdist,ic,i,ind
	integer range,minrange,maxrange,dmin,dmax
	integer k,krange,rdist
	integer, allocatable :: totdist(:), totrange(:)
	integer, allocatable :: dindex(:)
	integer, allocatable :: rindex(:)

	bsortdist = .true.
	bsortrange = .false.

	allocate(totdist(nkn))
	allocate(totrange(nkn))
	allocate(dindex(nkn))
	allocate(rindex(nkn))

	totdist = dist(:,0)
	do k=1,nkn
	  dmin = minval(dist(k,1:nroots))
	  dmax = maxval(dist(k,1:nroots))
	  totrange(k) = dmax - dmin
	end do

	call sort(nkn,totdist,dindex)
	call sort(nkn,totrange,rindex)
	call sort_invert(nkn,rindex)

	mindist = minval(totdist)
	maxdist = maxval(totdist)
	ic = count( dist(:,0) == maxdist )
	if( ic <= 0 ) stop 'error stop choose_new_source_node: ic'
	minrange = minval(totrange)
	maxrange = maxval(totrange)

	write(6,*) 'min/max dist: ',mindist,maxdist
	write(6,*) 'min/max range: ',minrange,maxrange

	write(6,*) 'sorted by distance:'
	do i=1,nkn
	  ind = dindex(i)
	  !if( totdist(ind) > maxdist - 5 ) then
	    !write(6,*) i,totdist(ind),totrange(ind)
	  !end if
	end do

	write(6,*) 'sorted by range:'
	do i=1,nkn
	  ind = rindex(i)
	  !if( totdist(ind) > maxdist - 5 ) then
	    !write(6,*) i,totdist(ind),totrange(ind)
	  !end if
	end do

	if( bsortdist ) then
	  minrange = nkn
	  do k=1,nkn
	    if( totdist(k) /= maxdist ) cycle
	    range = totrange(k)
	    if( range < minrange ) then
	      minrange = range
	      krange = k
	    end if
	  end do
	else if( bsortrange ) then
	  maxdist = 0
	  do k=1,nkn
	    if( totrange(k) /= minrange ) cycle
	    rdist = totdist(k)
	    if( rdist > maxdist ) then
	      maxdist = rdist
	      krange = k
	    end if
	  end do
	else
	  stop 'error stop choose_new_source_node: no sort algo'
	end if

	knew = krange

	write(6,*) 'choose: ',nroots,maxdist,minrange,knew

	end

!*******************************************************************

	subroutine exchange_nodes(np,nkn,npart,nel,nen3v)

	implicit none

	integer np
	integer nkn
	integer npart(nkn)
	integer nel
	integer nen3v(3,nel)

	integer ie,ii,iii,k,kk
	integer ic,icc
	integer ia,minc,maxc
	integer, allocatable :: icount(:)
	integer, allocatable :: counts(:)
	integer, allocatable :: ind(:)
	integer, allocatable :: matrix(:,:)

	allocate(icount(np))
	allocate(counts(np))
	allocate(ind(np))
	allocate(matrix(np,np))

	icount = 0
	matrix = 0

	do k=1,nkn
	  ic = npart(k)
	  icount(ic) = icount(ic) + 1
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    iii = 1+mod(ii,3)
	    kk = nen3v(iii,ie)
	    ic = npart(k)
	    icc = npart(kk)
	    if( ic == icc ) cycle
	    matrix(ic,icc) = matrix(ic,icc) + 1
	    matrix(icc,ic) = matrix(icc,ic) + 1
	  end do
	end do

	do ia=1,np
	  ic = icount(ia)
	  write(6,*) ia,ic
	end do

	minc = minval(icount)
	ic = findloc(icount,minc,1)
	write(6,*) 'exchanging: ',ic,minc
	counts(:) = matrix(ic,:)
	maxc = maxval(counts)
	icc = findloc(counts,maxc,1)
	write(6,*) 'with: ',icc,maxc

	do ia=1,np
	  if( counts(ia) /= 0 ) ind(ia) = 1
	end do

	!stop

	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine exchange_nodes_between_two(ic,icnew,nkn,npart,nel,nen3v)

	use queue
	use mod_connect
	use mod_flood

	implicit none

	integer ic
	integer icnew
	integer nkn
	integer npart(nkn)
	integer nel
	integer nen3v(3,nel)

	integer id
	integer ic1,ic2,nc1,nc2,nic1,nic2
	integer k,kk,i,ie,ii,ne,nk,ieo
	integer nfill
	integer iloop
	integer ics(3)
	integer, allocatable :: cinfo(:,:)
	integer icc,kcc,ncc
	integer :: ainfo(2,nkn)

        call queue_init(id)

	nc1 = count( npart == ic )
	nc2 = count( npart == icnew )

	if( nc1 > nc2 ) then		! convert ic to icnew
	  ic1 = ic
	  ic2 = icnew
	else				! convert icnew to ic
	  ic1 = icnew
	  ic2 = ic
	end if

	nc1 = count( npart == ic1 )
	nc2 = count( npart == ic2 )
	write(6,*) 'ic1 -> ic2',ic1,ic2,nc1,nc2

	do k=1,nkn
	  if( npart(k) /= ic1 ) cycle
	  call enqueue_if_color(id,ic2,k,nkn,npart,nel,nen3v)
	end do
	
	iloop = 0
	nfill = queue_fill(id)
	write(6,*) 'queue filled: ',nfill

	write(6,*) '           iloop   nfill     nc1     nc2' // &
     &			'    nic1    nic2'

        do
	  if( nc2 > nc1 ) exit
	  nfill = queue_fill(id)
	  if( nfill <= 0 ) exit
	  if( mod(iloop,100) == 0 ) then
	    nic1 = nkn
	    call floodfill_node_info(nkn,npart,ic1,nic1,ainfo)
	    nic2 = 1
	    !call floodfill_node(nkn,npart,ic2,nic2)
	    write(6,'(a,6i8)') 'in loop: ',iloop,nfill,nc1,nc2,nic1,nic2
	    if( nic1 > 1 ) then
	      ncc = minval(ainfo(2,:))
	      icc = findloc(ainfo(2,:),ncc,1)
	      kcc = ainfo(1,icc)
	      call floodfill_node_color(nkn,npart,ic2,kcc)
	      nc1 = nc1 - ncc
	      nc2 = nc2 + ncc
	      nic1 = nic1 - 1
	      write(6,*) 'detached area found: ',ncc
	      write(6,'(a,6i8)') 'in loop: ',iloop,nfill,nc1,nc2,nic1,nic2
	      iloop = iloop - 1
	    end if
	  end if
          if( .not. queue_dequeue(id,ieo) ) exit
	  do ii=1,3
	    kk = nen3v(ii,ieo)
	    ics(ii) = npart(kk)
	  end do
	  ic = count( ics == ic1 )
	  if( ic == 0 ) then			!already changed
	    ! nothing
	  else if( ic == 1 ) then		!change ic1 to ic2
	    !iloop = 0
	    nc1 = nc1 - 1
	    nc2 = nc2 + 1
	    ii = findloc(ics,ic1,1)
	    kk = nen3v(ii,ieo)
	    npart(kk) = ic2
	    call enqueue_if_color(id,ic1,kk,nkn,npart,nel,nen3v)
	  else					!put in queue again
	    if( ii < 4 ) call queue_enqueue(id,ieo)
	  end if
	  iloop = iloop + 1
	end do

        call queue_delete(id)

	write(6,*) 'ics: ',nc1,nc2
	nc1 = count( npart == ic1 )
	nc2 = count( npart == ic2 )
	write(6,*) 'ics: ',nc1,nc2

	call floodfill_node(nkn,npart,ic1,nic1)
	call floodfill_node(nkn,npart,ic2,nic2)

	write(6,*) 'connected: ',nic1,nic2
	write(6,*) 'iloop = ',iloop

	end

!*******************************************************************

	subroutine enqueue_if_color(id,ic,k,nkn,npart,nel,nen3v)

	use mod_connect
	use queue

	integer id,k,ic
	integer nkn
	integer npart(nkn)
	integer nel
	integer nen3v(3,nel)

	integer ne,i,ie,ii,kk

	  ne = elist(0,k)
	  do i=1,ne
	    ie = elist(i,k)
	    do ii=1,3
	      kk = nen3v(ii,ie)
	      if( npart(kk) == ic ) exit
	    end do
	    if( ii < 4 ) call queue_enqueue(id,ie)
	  end do

	end

!*******************************************************************

	subroutine divide_domain(ic,icnew,nkn,npart)

	use mod_connect

	implicit none

	integer ic
	integer icnew
	integer nkn
	integer npart(nkn)

	integer k,kk,i,nk,ng
	integer kroot,kdist
	integer maxdist,mingrade,nc
	integer, allocatable :: grade(:)
	integer, allocatable :: dist(:)
	integer, allocatable :: newpart(:)

!------------------------------------------------------
! find minimum grade of domain with color ic
!------------------------------------------------------

	allocate( grade(nkn) )
	allocate( dist(nkn) )
	allocate( newpart(nkn) )

	grade = nkn

	do k=1,nkn
	  if( npart(k) /= ic ) cycle
	  nk = nlist(0,k)
	  ng = 0
	  do i=1,nk
	    kk = nlist(i,k)
	    if( npart(kk) == ic ) ng = ng + 1
	  end do
	  grade(k) = ng
	end do

	mingrade = minval( grade )
	kroot = findloc(grade,mingrade,1)

!------------------------------------------------------
! find kroot and kdist
!------------------------------------------------------

	dist = -2
	where( npart == ic ) dist = -1

	call compute_distance_from_root(kroot,nkn,dist)

	maxdist = maxval( dist )

	mingrade = nkn
	do k=1,nkn
	  if( dist(k) == maxdist ) mingrade = min(mingrade,grade(k))
	end do
	do k=1,nkn
	  if( dist(k) == maxdist .and. grade(k) == mingrade ) exit
	end do
	if( k > nkn ) stop 'error stop divide_domain: no node'

	kdist = k
	write(6,*) kroot,kdist,maxdist,mingrade

!------------------------------------------------------
! find fill domain ic from kroot and kdist
!------------------------------------------------------

	newpart = -2
	where( npart == ic ) newpart = -1
	newpart(kroot) = ic
	newpart(kdist) = icnew

	call fill_domain(nkn,newpart)

	where( newpart /= -2 ) npart = newpart

!------------------------------------------------------
! end of routine - in npart are new domains ic and icnew
!------------------------------------------------------

	end

!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine ffill(nkn,nel,nen3v,ngrade,np,kroots,npart)

! shyparts custom routine

        implicit none

        integer nkn,nel
        integer nen3v(3,nel)
        integer ngrade(nkn)
        integer np
	integer kroots(np+1)
        integer npart(nkn)

	integer ie,ii,k,is,ncol,ic
	integer i,nchange,ia,dmax,id,iloop
	integer kn,ng
	integer kk(3)

	integer color1(nkn)
	integer color2(nkn)
	integer idist(nkn)

	color1 = 0
	idist = 0

	do i=1,np
	  k = kroots(i)
	  color1(k) = i
	  idist(k) = 1
	end do
	color2 = color1

	iloop = 0

	do
	  iloop = iloop + 1
	  id = iloop + 1
	  nchange = 0
	  color1 = color2
	  do ie=1,nel
	    ncol = 0
	    is = 0
	    do ii=1,3
	      k = nen3v(ii,ie)
	      kk(ii) = k
	      if( color1(k) /= 0 ) then
		ncol = ncol + 1
		is = is + ii
	      end if
	    end do
	    if( ncol > 0 .and. ncol < 3 ) then
	      nchange = nchange + 1
	      if( ncol == 1 ) then
	        ic = color1(kk(is))
	      else
		is = 6 - is
		ia = 1+mod(is,3)	!use any color of the two available
	        ic = color1(kk(ia))
	      end if
	      if( ic <= 0 ) stop 'error stop ffill: internal error (1)'
	      do ii=1,3
		k = kk(ii)
		if( color1(k) == 0 ) then
		  color2(k) = ic
		  idist(k) = id
		end if
	      end do
	    end if
	  end do
	  if( nchange == 0 ) exit
	end do

	npart = color1
	dmax = maxval(idist)

	kn = 0
	ng = nkn
	do k=1,nkn
	  if( idist(k) == dmax ) then
	    if( ngrade(k) < ng ) then
	      ng = ngrade(k)
	      kn = k
	    end if
	  end if
	end do
	if( kn == 0 ) stop 'error stop ffill: kn==0'

	kroots(np+1) = kn

	end

!*******************************************************************

	subroutine make_epart(nkn,nel,nen3v,npart,epart)

! make epart from npart - only elements fully in domain are colored, other -1

	implicit none

	integer nkn,nel
	integer nen3v(3,nel)
	integer npart(nkn)		!nodal partition [1-np]
	integer epart(nel)		!elemental partition [1-np] (return)

	integer ie,ii,k,ic,ics(3)

	epart = -1

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    ics(ii) = npart(k)
	  end do
	  ic = ics(1)
	  if( all( ics == ic ) ) epart(ie) = ic
	end do

	end

!*******************************************************************

	subroutine check_connected_domains(nkn,nel,npart,epart)

	use mod_flood

	implicit none

	integer nkn,nel
	integer npart(nkn)
	integer epart(nel)
	
	integer nmax,i,iae

	nmax = maxval(npart)
	write(6,*) 'stats: ',nmax
	call domain_stats(nmax,nel,epart)

	do i=0,nmax
	  call floodfill_elem(nel,epart,i,iae)
	  write(6,*) 'domain by elems: ',i,iae
	  if( iae > 1 ) goto 99
	  call floodfill_elem(nkn,npart,i,iae)
	  write(6,*) 'domain by nodes: ',i,iae
	  if( iae > 1 ) goto 99
	end do

	return
   99	continue
	write(6,*) 'domain is not connected: ',i,iae
	stop 'error stop check_connected_domains: not connected'
	end

!*******************************************************************

	subroutine domain_stats(nmax,n,part)

	implicit none

	integer nmax,n
	integer part(n)

	integer i,ic,nmin
	integer, allocatable :: count(:)

	nmin = -1
	allocate(count(nmin:nmax))

	count = 0

	do i=1,n
	  ic = part(i)
	  if( ic < nmin .or. ic > nmax ) goto 99
	  count(ic) = count(ic) + 1
	end do

	write(6,*) '----------------------'
	do i=nmin,nmax
	  write(6,*) i,count(i)
	end do
	write(6,*) '----------------------'

	return
   99	continue
	write(6,*) i,ic
	stop 'error stop domain_stats: value not possible'
	end

!*******************************************************************

