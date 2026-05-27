
!--------------------------------------------------------------------------
!
!    Copyright (C) 2013,2015,2017-2019  Georg Umgiesser
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

c revision log :
c
c 07.05.2013	ggu	copied from bastreat
c 14.05.2013	ggu	subroutines cleaned
c 10.10.2015	ggu	changed VERS_7_3_2
c 31.03.2017	ggu	changed VERS_7_5_24
c 22.02.2018	ggu	changed VERS_7_5_42
c 14.02.2019	ggu	changed VERS_7_5_56
c 16.02.2019	ggu	changed VERS_7_5_60
c 10.04.2021	ggu	better error handling and info output
c 10.11.2021	ggu	avoid warning for stack size
c 10.02.2022	ggu	better error message for not connected domain
c 16.02.2022	ggu	new routine basboxgrd()
c 09.03.2022	ggu	write also file index_sections.grd for the sections
c 10.05.2022	ggu	write element types on external element numbers
c 10.05.2022	ggu	new routine sort_multiple_sections() for predictability
c 18.05.2022	ggu	some more checks in sort_multiple_sections()
c 08.06.2022	ggu	save external nodes in basboxgrd for grd_write_item()
c 16.06.2022	ggu	write coloring errors to grd files
c 13.07.2022	ggu	forgot setting listold(0,:)
c 20.12.2024	ggu	general read for index file
c 22.05.2026	ggu	new routines for regularizing edges and connectivity
c 27.05.2026	ggu	finalized swapping elements to resolve non-connectivity
c
c****************************************************************

        subroutine basbox

c reads grid with box information and writes index file boxes.txt

	use mod_geom
	use mod_depth
	use evgeom
	use basin
	use basutil

	implicit none

	!integer ndim

        character*60 line
	integer node,nit
	integer mode,np,n,niter,i
        integer ner,nco,nknh,nelh,nli
	integer nlidim,nlndim
	integer idepth
	integer nlkdi
	logical bstop

	logical, parameter :: bnew = .true.	!new version of boxes.txt
	integer, parameter :: nbxdim = 100	!max number of boxes
	integer, parameter :: nlbdim = 250	!max number of boundary nodes

	integer nbox
	integer nblink(nbxdim)
	integer, allocatable :: boxinf(:,:,:)
	integer neib(nbxdim)
	integer iaux1(2,nlbdim)
	integer iaux2(nlbdim)

c-----------------------------------------------------------------
c check if we have read a bas file and not a grd file
c-----------------------------------------------------------------

        if( .not. breadbas .and. .not. bnew ) then
          write(6,*) 'for -box we need a bas file'
          stop 'error stop basbox: need a bas file'
        end if

c-----------------------------------------------------------------
c handle boxes
c-----------------------------------------------------------------

	allocate(boxinf(3,nlbdim,nbxdim))

	call check_box_connection

	call check_boxes(nbxdim,nbox,nblink)
	call handle_boxes(nbxdim,nlbdim,nbox,nblink,boxinf)

c-----------------------------------------------------------------
c write index file boxes.txt
c-----------------------------------------------------------------

	open(69,file='boxes.txt',form='formatted',status='unknown')

	call write_boxes(nbxdim,nlbdim,nbox,nblink,boxinf,bnew)
	call sort_boxes(nbxdim,nlbdim,nbox,nblink,boxinf,neib
     +			,iaux1,iaux2)

	close(69)
 
	call make_line_boxes

c-----------------------------------------------------------------
c end of routine
c-----------------------------------------------------------------

	write(6,*) 
	write(6,*) 'finished writing files'
	write(6,*) 'the index is in file boxes.txt'

	return
	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine sort_boxes(nbxdim,nlbdim,nbox,nblink,boxinf,neib
     +				,iaux1,iaux2)

c creates boxinf that contains information on sections
c boxinf(ii,n,ib)
c	ib is box number
c	n is an index of found intersections
c	ii=1	first node
c	ii=2	second node
c	ii=3	box number of neibor
c there are nblink(ib) node pairs in boxinf, so n=1,nblink(ib)

	use basin

	implicit none

	integer nbxdim,nlbdim,nbox
	integer nblink(nbxdim)
	integer boxinf(3,nlbdim,nbxdim)
	integer neib(nbxdim)
	integer iaux1(2,nlbdim)
	integer iaux2(nlbdim)
	integer ib,n,i,ii,ibn,nn,k
	integer jfill,j,ibb,nf,ns,nt,ntt
	integer id

	integer kantv(2,nkn)

	nt = 0
	ntt = 0
	id = 0
	do ib=1,nbox
	  neib(ib) = 0
	end do

	do k=1,nkn
	  kantv(1,k) = 0
	  kantv(2,k) = 0
	end do

	do ib=1,nbox
	  n = nblink(ib)
	  if( n .gt. 0 ) then
	    !write(6,*) '******** ',ib,n
	    do i=1,n
	      ibn = boxinf(3,i,ib)	!neibor box
	      neib(ibn) = neib(ibn) + 1
	    end do
	    !do ibn=1,nbox
	    do ibn=ib+1,nbox
	      nn = neib(ibn)
	      if( nn .gt. 0 ) then
		!write(6,*) '--- ',ib,ibn,nn
		neib(ibn) = 0
		jfill = 0
		do j=1,n
	          ibb = boxinf(3,j,ib)	!neibor box
		  if( ibn .eq. ibb ) then
		    jfill = jfill + 1
		    iaux1(1,jfill) = boxinf(1,j,ib)
		    iaux1(2,jfill) = boxinf(2,j,ib)
		  end if
		end do
		if( nn .ne. jfill ) then
			write(6,*) nn,jfill
			stop 'internal error'
		end if
		do j=1,jfill
		  !write(6,*) j,iaux1(1,j),iaux1(2,j)
		end do
		call sort_section(nlbdim,jfill,iaux1,nf,iaux2,nkn,kantv)
		!write(6,*) nf
		!write(6,*) (ipv(iaux2(i)),i=1,nf)

		call count_sections(nf,iaux2,ns)
		call invert_list(nf,iaux2)
		call sort_multiple_sections(ns,nf,iaux2,ipv)
		!write(66,*) ns,nf,ib,ibn
		!write(66,*) (ipv(iaux2(i)),i=1,nf)

		nt = nt + nf - ns + 1	!real nodes, no zeros
		ntt = ntt + nf + 1
		call write_section(nf,iaux2,id,ib,ibn,ipv)

	      end if
	    end do
	  end if
	end do

	write(69,*) 0,0,0,0
	write(6,*) 'total number of sections: ',id
	write(6,*) 'total number of nodes in sections: ',nt
	write(6,*) 'total number of needed nodes in sections: ',ntt

	end

c*******************************************************************

	subroutine sort_multiple_sections(ns,nf,list,ipv)

! sorts multiple section from one box to another to be predictable

	implicit none

	integer ns,nf
	integer list(nf)
	integer ipv(*)

	integer, allocatable :: listold(:,:)
	integer, allocatable :: listnew(:,:)
	integer, allocatable :: listaux(:,:)

	logical bdebug
	integer kint,kext,kmax
	integer is,ie,isect,n,i,id,nn
	integer nns,nnn
	integer ismax,nsect

	if( ns <= 1 ) return	!just one section

	bdebug = .true.
	bdebug = .false.

	call list_get_number_of_blocks(nf,list,nns,nnn)
	if( nns /= ns ) stop 'error stop sort_multiple_sections: (1)'

	allocate(listold(0:nf,ns))
	allocate(listnew(0:nf,ns))
	allocate(listaux(0:nf,ns))

	listold = 0
	listnew = 0
	listaux = 0

	if( bdebug ) then
	write(6,*) '==========================='
	write(6,*) '---------- listorig -----------'
	write(6,*) ns,nf
	write(6,*) list(1:nf)
	end if

	id = 0
	i = 1
	do while( i .le. nf )
	  is = i
	  do while( i .le. nf .and. list(i) .gt. 0 )	!FIXME - read past end of array
	    i = i + 1
	  end do
	  ie = i - 1
	  id = id + 1
	  nn = ie - is + 1
	  listold(1:nn,id) = list(is:ie)
	  listold(0,id) = nn
	  i = i + 1
	end do

	call list_split_blocks(ns,nf,nf,list,listaux)

	if( bdebug ) then
	write(6,*) '---------- listold -----------'
	do isect=1,ns
	  write(6,*) isect,listold(0,isect)
	  write(6,*) listold(1:nf,isect)
	end do
	end if

	if( any( listaux /= listold ) ) then
	  write(6,*) '---------- listaux -----------'
	  do isect=1,ns
	    write(6,*) isect,listaux(0,isect)
	    write(6,*) listaux(1:nf,isect)
	  end do
	  stop 'error stop sort_multiple_sections: listaux/=listold (3)'
	end if
	if( id /= ns ) then
	  stop 'error stop sort_multiple_sections: (2)'
	end if

	nsect = ns
	do while( nsect > 0 )
	  ismax = 0
	  kmax = 0
	  do is=1,nsect
	    kint = listold(1,is)
	    kext = ipv(kint)
	    if( kext > kmax ) then
	      kmax = kext
	      ismax = is
	    end if
	  end do
	  listnew(:,nsect) = listold(:,ismax)
	  listold(:,ismax) = listold(:,nsect)
	  nsect = nsect - 1
	end do

	if( bdebug ) then
	write(6,*) '---------- listnew -----------'
	do isect=1,ns
	  write(6,*) listnew(1:nf,isect)
	end do
	end if

	is = 1
	do isect=1,ns
	  do i=1,nf
	    if( listnew(i,isect) == 0 ) exit
	  end do
	  n = i - 1
	  ie = is + n - 1
	  list(is:ie) = listnew(1:n,isect)
	  list(ie+1) = 0
	  is = ie + 2
	end do
	  
	if( bdebug ) then
	write(6,*) '---------- listfinal -----------'
	write(6,*) list(1:nf)
	end if

	end

c*******************************************************************

	subroutine sort_section(nlbdim,jfill,iaux1,nf,iaux2,nkn,kantv)

c given a list of points sorts the nodes to be continuous
c deals with sections that are not connected
c inserts between these connections a 0

	implicit none

	integer nlbdim,jfill
	integer iaux1(2,nlbdim)
	integer nf
	integer iaux2(nlbdim)
	integer nkn
	integer kantv(2,nkn)

	integer k,k1,k2,j,kk,kn,kstop

	nf = 0
	do k=1,nkn	!only test... can be deleted
	  if( kantv(1,k) .ne. 0 ) goto 99
	  if( kantv(2,k) .ne. 0 ) goto 99
	end do

	do j=1,jfill
	  k1 = iaux1(1,j)
	  k2 = iaux1(2,j)
	  kantv(1,k1) = k2
	  kantv(2,k2) = k1
	end do

	do k=1,nkn
	  if( kantv(2,k) .ne. 0 ) then
	    kk = k
	    kstop = k
	    do while( kantv(2,kk) .gt. 0 .and. kantv(2,kk) .ne. kstop )
	      kk = kantv(2,kk)
	    end do
	    do while( kk .gt. 0 )
	      nf = nf + 1
	      iaux2(nf) = kk
	      kn = kantv(1,kk)
	      kantv(1,kk) = 0
	      kantv(2,kk) = 0
	      kk = kn
	    end do
	    nf = nf + 1
	    iaux2(nf) = 0	!insert 0 between non-continuous connections
	  end if
	end do

	nf = nf - 1	!skip trailing 0

	do k=1,nkn	!only test... can be deleted
	  if( kantv(1,k) .ne. 0 ) goto 99
	  if( kantv(2,k) .ne. 0 ) goto 99
	end do

	return
   99	continue
	write(6,*) k,kantv(1,k),kantv(2,k)
	stop 'error stop sort_section: internal error'
	end

c*******************************************************************

	subroutine count_sections(n,list,ns)

	implicit none

	integer n
	integer list(n)
	integer ns

	integer i

	ns = 1
	do i=1,n
	  if( list(i) .eq. 0 ) ns = ns + 1
	end do

	end

c*******************************************************************

	subroutine write_section(n,list,id,ib,ibn,ipv)

! nodes of section are written as external numbers

	implicit none

	integer n		!total number of nodes in section
	integer list(n)		!node numbers of section (internal)
	integer id		!id of section (consecutive)
	integer ib		!first box number (from-box)
	integer ibn		!second box number (to-box)
	integer ipv(*)		!external node numbers

	integer i,j,is,ie,nn

	!write(6,*) 
	i = 1
	do while( i .le. n )
	  is = i
	  do while( list(i) .gt. 0 )	!FIXME - read past end of array
	    i = i + 1
	  end do
	  ie = i - 1
	  id = id + 1
	  nn = ie - is + 1
	  !write(68,*) id,nn,ib,ibn
	  !write(68,*) (ipv(list(j)),j=is,ie)
	  write(69,*) id,nn,ib,ibn
	  write(69,'((8i9))') (ipv(list(j)),j=is,ie)
	  i = i + 1
	end do

	end

c*******************************************************************

	subroutine invert_list(n,list)

	implicit none

	integer n
	integer list(n)

	integer i,j,iaux

	do i=1,n/2
	  j=n+1-i
	  iaux = list(i)
	  list(i) = list(j)
	  list(j) = iaux
	end do

	end

c*******************************************************************

	subroutine write_boxes(nbxdim,nlbdim,nbox,nblink,boxinf,bnew)

! in new version external element numbers are written

	use basin

	implicit none

	integer nbxdim,nlbdim,nbox
	integer nblink(nbxdim)
	integer boxinf(3,nlbdim,nbxdim)
	logical bnew				!new version of boxes.txt

	integer ib,n,i,ii,ie
	integer nu

	integer, parameter :: idbox = 473226
	integer, parameter :: ftype = 9
	integer, parameter :: nversbox = 3

	nu = 0

	do ib=1,nbox
	  n = nblink(ib)
	  if( n .gt. 0 ) then
	    nu = nu + 1
	    !write(6,*) ib,n
	    do i=1,n
	      !write(6,*) (boxinf(ii,i,ib),ii=1,3)
	    end do
	  end if
	end do

	if( bnew ) then
	  write(69,*) idbox,ftype,nversbox
	  write(69,*) nel,nbox,nu
	  do ie=1,nel
	    write(69,*) ipev(ie),iarv(ie)
	  end do
	else
	  write(69,*) nel,nbox,nu
	  write(69,'((8i9))') (iarv(ie),ie=1,nel)
	end if

	end

c*******************************************************************

	subroutine check_boxes(nbxdim,nbox,nblink)

	use mod_geom
	use basin

	implicit none

	integer nbxdim,nbox
	integer nblink(nbxdim)

	integer ib,ie,ien,ii,i1,i2
	integer ia,n,nb


	nbox = 0
	do ib=1,nbxdim
	  nblink(ib) = 0
	end do
	
	do ie=1,nel
	  ia = iarv(ie)
	  if( ia .le. 0 ) goto 99
	  if( ia .gt. nbxdim ) goto 99
	  nblink(ia) = nblink(ia) + 1
	  nbox = max(nbox,ia)
	end do

	nb = 0
	do ib=1,nbxdim
	  if( nblink(ib) .gt. 0 ) nb = nb + 1
	end do
	
	write(6,*) 'total number of boxes: ',nbox,nb

	return
   99	continue
	write(6,*) 'ia,nbxdim: ',ia,nbxdim
	stop 'error stop check_boxes: nbxdim'
	end

c*******************************************************************

	subroutine handle_boxes(nbxdim,nlbdim,nbox,nblink,boxinf)

c sets up box index

	use mod_geom
	use basin

	implicit none

	integer nbxdim			!box dimension
	integer nlbdim			!boundary node dimension 
	integer nbox			!total number of boxes (return)
	integer nblink(nbxdim)		!total number of boundary nodes for box
	integer boxinf(3,nlbdim,nbxdim)	!info for boundary nodes of box

	integer ib,ie,ien,ii,i1,i2
	integer ia,ian,n
	integer nboxes(nbxdim)

	nbox = 0
	nblink = 0
	nboxes = 0
	
	ian = 0

	do ie=1,nel
	  ia = iarv(ie)				!element type - is box number
	  if( ia .le. 0 ) goto 99
	  if( ia .gt. nbxdim ) goto 99
	  nboxes(ia) = nboxes(ia) + 1
	  nbox = max(nbox,ia)
	  do ii=1,3
	    ien = ieltv(ii,ie)
	    ian = 0
	    if( ien .gt. 0 ) ian = iarv(ien)
	    if( ian .gt. nbxdim ) goto 99
	    if( ien .gt. 0 .and. ian .ne. ia ) then  !neighbor is diff box
	      i1 = mod(ii,3) + 1
	      i2 = mod(ii+1,3) + 1
	      n = nblink(ia) + 1
	      if( n .gt. nlbdim ) goto 98
	      boxinf(1,n,ia) = nen3v(i1,ie)
	      boxinf(2,n,ia) = nen3v(i2,ie)
	      boxinf(3,n,ia) = ian
	      nblink(ia) = n
	    end if
	  end do
	end do

	write(6,*) 'box info:'
	!write(6,*) ' box number    elements   bnd-nodes'
	!do ia=1,nbox
	!  write(6,*) ia,nboxes(ia),nblink(ia)
	!end do
	call write_box_info(nbox,nboxes)

	call handle_regularize_edges
	stop

	return
   98	continue
	write(6,*) 'too many boundary nodes in box ',ia
	write(6,*) 'n,nlbdim: ',n,nlbdim
	write(6,*) 'increase nlbdim and recompile'
	stop 'error stop handle_boxes: nlbdim'
   99	continue
	write(6,*) 'too many boxes...'
	write(6,*) 'ia,ian,nbxdim: ',ia,ian,nbxdim
	write(6,*) 'increase nbxdim and recompile'
	stop 'error stop handle_boxes: nbxdim'
	end

c*******************************************************************
c*******************************************************************
c*******************************************************************

	subroutine check_box_connection

c checks if all boxes are connected

	use basin
	use mod_geom

	implicit none

	integer ie,k,ii
	integer i,j,nc,ic,nt,nnocol,icerror,icn,ico,ntot
	integer ijmin,ijmax
	integer icol,icolmax
	integer ierr,ierr_swap
	integer nmin,nmax
	integer icolor(nel)
	integer icon(nel)
	integer list(3,nel)	!insert not connected areas in this list

	integer ieext,ipext

	icerror = 0
	ic = 0
	ierr = 0
	icolmax = 0
	do ie=1,nel
	  icon(ie) = 0
	  icolor(ie) = 0
	end do

	do ie=1,nel
	  if( icon(ie) .eq. 0 ) then
	    ic = ic + 1
	    call color_box_area(ie,icon,icol,nc)
	    list(1,ic) = icol
	    list(2,ic) = nc
	    list(3,ic) = ie
	    !write(6,*) icol,nc
	    if( icol .gt. nel ) goto 99
	    if( icolor(icol) .ne. 0 ) then
	      !write(6,*) '*** area not connected: ',icol
	      ierr = ierr + 1
	    end if
	    icolor(icol) = icolor(icol) + nc
	    icolmax = max(icolmax,icol)
	  end if
	end do

! list(1,i)	area code
! list(2,i)	total number of elements
! list(3,i)	start element index

	if( ierr > 0 ) then	!here error treatment for not connected areas
	 do i=1,ic
	  icol = list(1,i)
	  do j=i+1,ic
	    if( list(1,j) == icol ) then
	        if( list(2,i) > list(2,j) ) then
		  ijmax = i
		  ijmin = j
		else
		  ijmax = j
		  ijmin = i
		end if
		ie = list(3,ijmin)
		nc = list(2,ijmin)
		ico = list(1,ijmin)
		ierr_swap = 1
		if( nc == 1 ) then
		  write(6,*) 'single not connected element found: ',ico
		  call exchange_one_element(ie,ierr_swap)
		  if( ierr_swap == 0 ) then
		    icn = iarv(ie)
		    write(6,*) 'swapping colors works ',ico,' -> ',icn
		  else
		    write(6,*) 'could not automatically resolve connections'
		  end if
		end if
		if( ierr_swap == 0 ) cycle
		icerror = icerror + 1
		write(6,*) '==================================='
		write(6,*) 'not connected area found ',icol
	        write(6,*) 'area            elements'
     +		// '   elem number (int)'
     +		// '   elem number (ext)'
		write(6,1123) list(:,ijmax),ieext(list(3,ijmax))
		write(6,1123) list(:,ijmin),ieext(list(3,ijmin))
	        write(6,*) 'not connected area contains element (int/ext)'
     +			,ie,ieext(ie)
	        write(6,*) 'coordinates of element are: '
		do ii=1,3
		  k = nen3v(ii,ie)
		  write(6,*) ii,ipext(k),xgv(k),ygv(k)
		end do
	        call write_grd_error(icerror,list(3,i),list(3,j))
		write(6,*) '==================================='
	    end if
	  end do
	 end do
	end if

	nnocol = count( icon == 0 )
	if( nnocol > 0 ) goto 96

	nt = 0
	ic = 0
	icol = 0
	nmin = nel
	nmax = 0
	do ie=1,nel
	  if( icolor(ie) .gt. 0 ) then
	    icol = ie
	    ic = ic + 1
	    nc = icolor(ie)
	    nmax = max(nmax,nc)
	    nmin = min(nmin,nc)
	    nt = nt + nc
	  end if
	end do
	
	if( ic /= icolmax ) goto 97

	write(6,*) 
	write(6,*) 'number of boxes: ',ic
	write(6,*) 'maximum box index: ',icol
	write(6,*) 'largest area contains element: ',nmax
	write(6,*) 'smallest area contains element: ',nmin
	write(6,*) 

	if( ierr /= 0 ) then
	  call write_box_info(ic,icolor)
	end if

	if( nel .ne. nt ) goto 98

	call write_grd_from bas('swapboxes.grd')

	if( ierr .gt. 0 ) then
	  write(6,*) 'There are non connected areas: ',icerror
	  write(6,*) 'files are in error_*.grd'
	  write(6,*) 'new boxes are in swapboxes.grd'
	  stop 'error stop check_box_connection: errors'
	else
	  write(6,*) 'all areas are connected'
	end if

 1123	format(i5,3i20)
	return
   96	continue
	write(6,*) 'some elements could not be colored: ',nnocol
	write(6,*) 'the reason for this is unknown'
	stop 'error stop check_box_connection: could not color'
   97	continue
	write(6,*) 'ic,icolmax: ',ic,icolmax
	stop 'error stop check_box_connection: area numbers have holes'
   98	continue
	write(6,*) 'nel,nt: ',nel,nt
	stop 'error stop check_box_connection: internal error'
   99	continue
	write(6,*) 'icol = ',icol
	stop 'error stop check_box_connection: color code too high'
	end

c*******************************************************************

	subroutine color_box_area(iestart,icon,icol,nc)

c colors all accessible elements starting from iestart
c colors only elements with same area code
c uses this area code to color the elements
c area code 0 is not allowed !!!!

	use basin
	use mod_geom

	implicit none

	integer iestart
	integer icon(nel)
	integer icol		!color used (return)
	integer nc		!total number of elements colored (return)

	integer ip,ien,ii,ie
	integer list(nel)	!stack to hold still to color elements

	nc = 0
	ip = 1
	list(ip) = iestart
	icol = iarv(iestart)
	if( icol .le. 0 ) goto 98

	do while( ip .gt. 0 ) 
	  ie = list(ip)
	  if( icon(ie) .le. 0 ) nc = nc + 1
	  icon(ie) = icol
	  !write(6,*) icol,ip,ie
	  ip = ip -1
	  do ii=1,3
	    ien = ieltv(ii,ie)
	    if( ien .gt. 0 ) then
	      if( icon(ien) .eq. 0 .and. iarv(ien) .eq. icol ) then
	        ip = ip + 1
	        list(ip) = ien
	      end if
	    end if
	  end do
	end do

	return
   98	continue
	write(6,*) 'cannot color with icol < 1 : ',icol
	stop 'error stop color_box_area: iarv'
   99	continue
	write(6,*) ip,ie,ien,icol,icon(ien)
	stop 'error stop color_box_area: internal error (1)'
	end

!*******************************************************************

	subroutine check_element_connection(iecheck,ierr)

! checks if color of element is unique (connected)

	use basin

	implicit none

	integer iecheck		!check this element
	integer ierr		!error code: 0 if ok, else not connected

	integer icol,nc,ie
	integer icon(nel)

	icon = 0
	call color_box_area(iecheck,icon,icol,nc)
	
	ierr = 0

	do ie=1,nel
	  if( iarv(ie) == icol .and. icon(ie) == 0 ) then
	    !write(6,*) 'area not connected: ',ie,icol,iarv(ie),icon(ie)
	    ierr = ierr + 1
	  end if
	end do

	!if( ierr > 0 ) then
	!  write(6,*) 'initial color was: ',iarv(iecheck),nc,ierr
	!end if

	end

!*******************************************************************

	subroutine exchange_one_element(ie,ierr)

	use basin
	use mod_geom

	implicit none

	integer ie
	integer ierr

	logical bw
	integer ii,i
	integer icol,ii_start,icol_new,iii,ien,icoln,ncol
	integer ie_start,ie_try,ii_aux
	integer iens(3)
	integer icols(3)

	ierr = 1	!fake error

	icol = iarv(ie)

!	--------------------------------------
!	find colors of neigboring elements
!	exclude areas that are already not connected
!	--------------------------------------

	icols = 0
	do ii=1,3
	  iens(ii) = ieltv(ii,ie)
	  if( iens(ii) > 0 ) then
	    icols(ii) = iarv(iens(ii))
	    call check_element_connection(iens(ii),ierr)
	    if( ierr /= 0 ) icols(ii) = 0
	  end if
	end do

	if( any(icols==icol) ) goto 99	!this would mean it is connected

!	--------------------------------------
!	find element with predominant color and start from this one
!	--------------------------------------

	ncol = 2
	if( icols(1) == icols(2) ) then
	  ii_start = 1
	  if( icols(1) == icols(3) ) ncol = 1
	else if( icols(1) == icols(3) ) then
	  ii_start = 1
	else if( icols(2) == icols(3) ) then
	  ii_start = 2
	else					!all colors are different
	  ii_start = 1
	  ncol = 3
	end if

!	--------------------------------------
!	loop over colors of element and neighbors
!	--------------------------------------

	ii_start = 1 + mod(ii_start+1,3)	!ii_start-1

	ierr = 0

	bw = ( icol == 33 )
	bw = .false.

	if( bw ) then
	  write(6,*) '======================================='
	  write(6,*) 'swapping color icol: ',icol
	  write(6,*) '======================================='
	  write(6,*) ncol,' different colors in neibors found'
	  do ii=1,3
	    write(6,*) ii,iens(ii),icols(ii)
	  end do
	end if

	ierr = -1

	do ii=1,3
	  ii_start = 1 + mod(ii_start,3)
	  icol_new = icols(ii_start)
	  if( icol_new == 0 ) cycle
	  iarv(ie) = icol_new
	  if( bw ) write(6,*) 'trying color: ',icol_new
	  call check_element_connection(ie,ierr)
	  if( bw ) write(6,*) 'checking: ',0,icol,ierr
	  if( ierr /= 0 ) cycle	!old elem still not connected
	  do iii=1,3			!now check the surrounding areas
	    ien = iens(iii)
	    icoln = icols(iii)
	    if( ien == 0 ) cycle
	    if( icoln == 0 ) cycle
	    call check_element_connection(ien,ierr)
	    if( bw ) write(6,*) 'checking: ',iii,icoln,ierr
	    if( ierr /= 0 ) exit	!area not connected
	  end do
	  if( ierr == 0 ) exit		!this color works
	end do

	if( ierr == 0 ) then
	  if( bw ) write(6,*) 'swapping colors works: ',icol,icol_new
	else
	  iarv(ie) = icol		!set back old color
	  if( bw ) write(6,*) 'cannot automatically reconnect area...'
	end if

	return
   99	continue
	write(6,*) 'area is connected: ',ie,ipev(ie)
	stop 'error stop exchange_one_element: internal error (1)'
	end

!*******************************************************************
!*******************************************************************
!*******************************************************************
! here create grd file from index file
!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine basboxgrd

	use basin
	use basutil

	implicit none

	integer ios
	integer ne,nbox,nmax
	integer ie,ia,itot
	integer iline,itype
	integer is,nnodes,ib1,ib2
	integer k,kext,i
	integer nout,nvers
	integer, allocatable :: iboxes(:)
	integer, allocatable :: extnodes(:),nodes(:),used(:)
	real x,y
	real, parameter :: flag = -999.
	character*80 file

	integer ipint

	write(6,*) 're-creating grd file from bas and index.txt'

!-----------------------------------------------------------
! open file
!-----------------------------------------------------------

        if( .not. breadbas ) then
          write(6,*) 'for -boxgrd we need a bas file'
          stop 'error stop basboxgrd: need a bas file'
        end if

	if( index_file == ' ' ) then
	  write(6,*) 'no index file given...'
	  stop 'error stop basboxgrd: no index file given'
	end if

	open(1,file=index_file,status='old',form='formatted',iostat=ios)
	if( ios /= 0 ) then
	  write(6,*) 'cannot open file: ',trim(index_file)
	  stop 'error stop basboxgrd: no such index file'
	end if
	write(6,*) 'reading index file...'

!-----------------------------------------------------------
! read box header info
!-----------------------------------------------------------

	call boxes_read_header(ne,nbox,nmax,nvers)

!-----------------------------------------------------------
! read box information and write index.grd
!-----------------------------------------------------------

	write(6,*) 'reading areas of index file...'

	allocate(iboxes(ne))
	call boxes_read_boxes(ne,iboxes,nbox,nmax,nvers)
	
	iarv = iboxes

        call basin_to_grd
        call grd_write('index.grd')
        write(6,*) 'The basin has been written to index.grd'

!-----------------------------------------------------------
! read section information and write index_sections.grd
!-----------------------------------------------------------

	nout = 2
	file = 'index_sections.grd'
	open(nout,file=file,status='unknown',form='formatted')
	write(6,*) 'reading sections of index file...'

	allocate(extnodes(nkn),nodes(nkn),used(nkn))
	used = 0
	iline = 0

	do
	  call boxes_read_section(nnodes,nodes,extnodes,ib1,ib2)
	  if( nnodes == 0 ) exit
	  iline = iline + 1
	  itype = iline
	  do i=1,nnodes
	    k = nodes(i)
	    kext = extnodes(i)
	    if( used(k) /= 0 ) cycle	!do not write nodes more than once
	    used(k) = 1
	    x = xgv(k)
	    y = ygv(k)
	    call grd_write_node(nout,kext,0,x,y,flag)
	  end do
          call grd_write_item(nout,3,iline,itype,nnodes,extnodes,flag)
	end do

	close(nout)
	close(1)
        write(6,*) 'The sections have been written to ',trim(file)

!-----------------------------------------------------------
! end of routine
!-----------------------------------------------------------

	return
	end

!*******************************************************************

	subroutine boxes_read_header(ne,nbox,nmax,nvers)

	use basin

	implicit none

	integer ne,nbox,nmax,nvers

	integer id,ftype
        integer, parameter :: idbox = 473226            !id for box files

	read(1,*) id,ftype,nvers

        if( id /= idbox .or. ftype /= 9 ) then  !old read: nel,nbox,nb
          ne = id
          nbox = ftype
          nmax = nvers
	  nvers = 2
        else                                    !new read: id,ftype,nvers
          read(1,*) ne,nbox,nmax
	end if

        if( ne .ne. nel ) goto 99
	if( nvers < 2 .or. nvers > 3 ) goto 98

	return
   98   continue
        write(6,*) 'nvers = ',nvers
        stop 'error stop boxes_read_header: cannot read version'
   99   continue
        write(6,*) nel,ne
        stop 'error stop boxes_read_header: nel/=nelbox'
	end

!*******************************************************************

	subroutine boxes_read_boxes(ne,iboxes,nbox,nmax,nvers)

	use basin

	implicit none

	integer ne
	integer iboxes(ne)
	integer nbox,nmax
	integer nvers

	integer ie,ie_int,ie_ext,ia,ic,itot
	integer, allocatable :: icount(:)

	integer ieint

	allocate(icount(0:nmax))

        iboxes = -1

	if( nvers <= 2 ) then
	  read(1,*) (iboxes(ie),ie=1,ne)
	else
          do ie=1,ne
            read(1,*) ie_ext,ia
            ie_int = ieint(ie_ext)
            if( ie_int > 0 ) iboxes(ie_int) = ia
          end do
	end if

        ic = count( iboxes == -1 )
        if( ic > 0 ) goto 99

	icount = 0
	do ie=1,ne
	  ia = iboxes(ie)
	  if( ia < 0 .or. ia > nmax ) goto 98
	  icount(ia) = icount(ia) + 1
	end do

	write(6,'(a)') '        ia     count'
	itot = 0
	do ia=0,nmax
	  write(6,'(2i10)') ia,icount(ia)
	  itot = itot + icount(ia)
	end do
	write(6,*) 'total: ',itot
	if( itot /= ne ) goto 97

	return
   97	continue
        stop 'error stop boxes_read_boxes: itot /= ne'
   98	continue
        write(6,*) 'area code out of bounds: ',ia
        write(6,*) 'possible between: ',0,nmax
        stop 'error stop boxes_read_boxes: erroneous area code'
   99	continue
        write(6,*) 'some elements are not part of a box: ',ic
        stop 'error stop boxes_read_boxes: incomplete box information'
	end

!*******************************************************************

	subroutine boxes_read_section(nnodes,nodes,extnodes,ib1,ib2)

	use basin

	implicit none

	integer nnodes
	integer ib1,ib2
	integer nodes(nkn)
	integer extnodes(nkn)

	integer is,i
	integer, save :: iline = 0
	integer ipint

	read(1,*) is,nnodes,ib1,ib2
	if( is == 0 ) return
	iline = iline + 1
	if( is /= iline ) goto 99
	read(1,*) (extnodes(i),i=1,nnodes)
	do i=1,nnodes
	  nodes(i) = ipint(extnodes(i))	!nodes are external numbers!!!!
	end do

	if( iline == 1 ) then
	  write(6,'(a)') '   section     nodes   ib_from     ib_to'
	end if
	write(6,'(4i10)') is,nnodes,ib1,ib2

	return
   99	continue
	write(6,*) is,iline
	stop 'error stop boxes_read_section: is /= iline'
	end

!*******************************************************************

	subroutine write_grd_error(icerror,ie1,ie2)

	use basin

	implicit none

	integer icerror,ie1,ie2

	integer i,k,ie
	integer nc1,nc2,icol1,icol2,ic1,ic2,ic0
	integer inext(nkn)
	integer ieext(nel)
	integer intype(nkn)
	integer ietype(nel)
	integer icon(nel)
	character*80 text,file,string

	text  = 'coloring'
	write(string,'(i3)') icerror
	do i=1,3
	  if( string(i:i) == ' ' ) string(i:i) = '0'
	end do
	file = 'error_' // trim(string) // '.grd'
	write(6,*) 'writing file: ',trim(file)

        do k=1,nkn
          inext(k) = ipv(k)
        end do
        do ie=1,nel
          ieext(ie) = ipev(ie)
        end do

	icon = 0
	!iarv = 1
	call color_box_area(ie1,icon,icol1,nc1)
	ic1 = count( icon == icol1 )
	where( icon == icol1 ) icon = icol1+1
	!iarv = 2
	call color_box_area(ie2,icon,icol2,nc2)
	ic2 = count( icon == icol2 )
	write(6,*) 'coloring 1: ',ie1,nc1,icol1,ic1,nel
	write(6,*) 'coloring 2: ',ie2,nc2,icol2,ic2,nel
	ic0 = count( icon == 0 )
	write(6,*) 'coloring 0: ',ic0,nel
	where( icon == 0 ) icon = -1

	intype = 0
	ietype = icon

	call write_grd_file(file,text,nkn,nel,xgv,ygv,nen3v
     +                  ,inext,ieext,intype,ietype)

	end

!*****************************************************************

	subroutine make_line_boxes

! creates lines of boxes from box type of elements

	use basin
	use mod_geom

	implicit none

	logical binsert,bdebug,bfound,bterminal
	integer ie,ii,ia,ien,ian,nbox,nmax,ntot,it
	integer i,ibase,nc,ic,ncc,iline,itot,icum
	integer ii1,ii2,k1,k2
	integer istart,iend
	integer iu,iu1
	integer, allocatable :: ncount(:)
	integer, allocatable :: icount(:)
	integer, allocatable :: ibound(:,:)
	integer, allocatable :: istartend(:,:)
	integer, allocatable :: icouples(:,:)
	integer, allocatable :: iorder(:)
	integer, allocatable :: line_nodes(:)
	integer, allocatable :: iuse(:)
	character*10 string
	character*80 header
	character*80 file1,file2

	bdebug = .true.
	bdebug = .false.
	bterminal = .false.

	write(6,*) 'creating edge list of boxes'

	file1 = 'line_boxes.grd'
	file2 = 'line_boxes.txt'

! count number of boxes

	nbox = 0
	do ie=1,nel
	  ia = iarv(ie)
	  if( ia > nbox ) nbox = ia
	end do
	if( nbox /= maxval(iarv) ) stop 'internal error 7'
	if( 1 /= minval(iarv) ) stop 'internal error 8'

	allocate(ncount(nbox))
	allocate(icount(0:nbox))
	ncount = 0

! count number of edges for each box

	do ie=1,nel
	  ia = iarv(ie)
	  do ii=1,3
	    ien = ieltv(ii,ie)
	    if( ien > 0 ) then
	      ian = iarv(ien)
	      if( ia /= ian ) then
		ncount(ia) = ncount(ia) + 1
	        !write(6,*) ' c   ',ie,ien,ia,ian
	      end if
	    else
	      ncount(ia) = ncount(ia) + 1	!is boundary side
	    end if
	  end do
	end do

! get cumulative index for insertion of edges

	nmax = 0
	ntot = 0
	write(6,*) 'nbox = ',nbox
	icount(0) = 0
	do ia=1,nbox
	  nmax = max(nmax,ncount(ia))
	  ntot = ntot + ncount(ia)
	  icount(ia) = icount(ia-1) + ncount(ia)
	  write(6,*) ia,ncount(ia)
	end do

	write(6,*) 'nmax = ',nmax
	write(6,*) 'ntot = ',ntot

! retrieve edges (node numbers) and put into icouples

	ncount = 0
	allocate(icouples(2,ntot))

	do ie=1,nel
	  ia = iarv(ie)
	  do ii=1,3
	    binsert = .false.
	    ien = ieltv(ii,ie)
	    if( ien > 0 ) then
	      ian = iarv(ien)
	      if( ia /= ian ) then
	        binsert = .true.
	      end if
	    else
	      binsert = .true.
	    end if
	    if( binsert ) then		! edge found
	      ncount(ia) = ncount(ia) + 1
	      it = icount(ia-1) + ncount(ia)
	      ii1 = 1 + mod(ii,3)
	      ii2 = 1 + mod(ii1,3)
	      k1 = nen3v(ii1,ie)
	      k2 = nen3v(ii2,ie)
	      icouples(1,it) = k1
	      icouples(2,it) = k2
	    end if
	  end do
	end do

! run over boxes and compute ordered list of nodes
! ibound is aux array that contains unordered edges of box

	allocate(ibound(2,nmax))
	allocate(istartend(2,nmax))
	allocate(iorder(nmax))
	allocate(line_nodes(nmax))
	allocate(iuse(nkn))
	iuse = 0

	iu = 88
	open(iu,file=file1,status='unknown',form='formatted')
	open(iu+1,file=file2,status='unknown',form='formatted')
	write(iu+1,*) nbox

	header = '   ibox  iline istart   iend   itot   icum     nc'
	if( bterminal ) write(6,'(a)') header

	do ia=1,nbox
	  ibase = icount(ia-1)
	  nc = ncount(ia)
	  ibound(:,1:nc) = icouples(:,ibase+1:ibase+nc)

	  call order_edges(nc,ibound,iorder,iline,istartend)

	  icum = 0
	  do i=1,iline
	    itot = istartend(2,i) - istartend(1,i) + 1
	    icum = icum + itot
	    if( bterminal ) then
	      write(6,'(7i7)') ia,i,istartend(:,i),itot,icum,nc
	    end if
	  end do
	  if( icum /= nc ) stop 'error stop intern 11'

	  istart = istartend(1,1)
	  iend = istartend(2,1)
	  if( iline > 1 ) then
	    call find_outer_line(iline,istartend,istart,iend)
	  end if
	  itot = iend - istart + 1
	  if( iline > 1 .and. bterminal ) then
	    write(6,*) 'outerline: ',istart,iend,itot
	  end if
	  line_nodes(1:itot) = iorder(istart:iend)

	  call write_line(iu,ia,itot,line_nodes,nkn,iuse)
	end do

	write(6,*) 'lines have been written to files:'
	write(6,*) trim(file1)
	write(6,*) trim(file2)

	end

!*****************************************************************

	subroutine order_edges(ncc,ibound,iorder,iline,istartend)

! takes a list of edges (ibound) and creates continuous ordered list
! more than one line can be found (islands)
! in iline is total number of lines found
! in istartend are indices of start and end of lines in iorder

	implicit none

	integer ncc		  ! total number of edges
	integer ibound(2,ncc)	  ! edges
	integer iorder(ncc)	  ! ordered nodes
	integer iline		  ! total number of lines found (return)
	integer istartend(2,ncc)  ! start and end of line in iorder (return)

	logical bdebug,bfound	
	integer icum,nc,i
	integer ic,icstart,icend,ictot

	bdebug = .false.

	if( bdebug ) then
	  write(6,*) 'original'
	  do i=1,nc
	    write(6,*) i,ibound(:,i)
	  end do
	  write(6,*) 'ordering'
	end if

	iorder = 0
	icum = 0

	nc = ncc
	ic = 0
	iline = 0

	do while( nc > 0 ) 

	  ic = ic + 1
	  icstart = ic
	  iorder(ic) = ibound(2,1)
	  if( bdebug ) write(6,*) ic,nc,ibound(:,1)
	  ibound(:,1) = ibound(:,nc)
	  nc = nc - 1

	  do while( nc > 0 )
	    bfound = .false.
	    do i=1,nc
	      if( ibound(1,i) == iorder(ic) ) then
		bfound = .true.
	        ic = ic + 1
		if( bdebug ) write(6,*) ic,nc,ibound(:,i)
	        if( ic > ncc ) stop 'error stop intern 4'
	        iorder(ic) = ibound(2,i)
	  	ibound(:,i) = ibound(:,nc)
		exit
	      end if
	    end do
	    if( bfound ) then
	      nc = nc - 1
	    else
	      exit
	    end if
	  end do

	  iline = iline + 1
	  icend = ic
	  ictot = icend - icstart + 1
	  icum = icum + ictot
	  istartend(1,iline) = icstart
	  istartend(2,iline) = icend
	  !write(6,'(7i6)') ncc,nc,ic,ictot,icum,iline
	end do

	end

!*****************************************************************

	subroutine find_outer_line(iline,istartend,istart,iend)

! returns line with maximum nodes
!
! could not be outer line
! real program should test inclusion of lines into others

	implicit none

	integer iline
	integer istartend(2,iline)
	integer istart,iend

	integer il,imax,i,itot

	il = 0
	imax = 0
	do i=1,iline
	  itot = istartend(2,i) - istartend(1,i) + 1
	  if( itot > imax ) then
	    imax = itot
	    il = i
	  end if
	end do

	istart = istartend(1,il)
	iend = istartend(2,il)

	!itot = istartend(2,il) - istartend(1,il) + 1
	!write(6,*) 'outerline: ',il,itot

	end

!*****************************************************************

	subroutine write_line(iu,ia,itot,line_nodes,n,iuse)

	use basin

	implicit none

	integer iu
	integer ia
	integer itot
	integer line_nodes(itot)
	integer n
	integer iuse(n)

	integer iu1,i,k
	integer is,ie
	real x,y

	iu = 88
	iu1 = iu + 1

	write(iu1,*) ia,itot

	do i=1,itot
	  k = line_nodes(i)
	  x = xgv(k)
	  y = ygv(k)
	  if( iuse(k) == 0 ) then
	    write(iu,'(i1,2i8,2e16.8)') 1,k,ia,x,y
	  end if
	  iuse(k) = iuse(k) + 1
	  write(iu1,*) i,x,y
	end do

	write(iu,'(i1,3i8)') 3,ia,ia,itot+1	!first node twice
	
	ie = 0
	do
	  is = ie + 1
	  if( is > itot ) exit
	  ie = ie + 10
	  if( ie > itot ) ie = itot
	  write(iu,'(10i8)') line_nodes(is:ie)
	  !write(6,*) 'line: ',ia,is,ie
	end do
	write(iu,'(10i8)') line_nodes(1:1)	!close line

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine handle_regularize_edges

	implicit none

	integer itot_start,itot_end
	character*80 file

	call regularize_edges(itot_start)
	call regularize_edges(itot_end)
	write(6,*) 'regularizing elements before/after: '
     &		,itot_start,itot_end
	call check_box_connection

	file = 'newboxes.grd'
	write(6,*) 'writing files with new boxes: ',trim(file)
	call write_grd_from bas(file)

	end

!*****************************************************************

	subroutine regularize_edges(itot)

	use basin
	use mod_geom

	implicit none

	integer ie,ia,ic,ii,ien
	integer ian,ians(3)
	integer ii1,ii2,ii3
	integer itot

	itot = 0

	do ie=1,nel
	  ia = iarv(ie)
	  ians = -1
	  do ii=1,3
	    ien = ieltv(ii,ie)
	    if( ien <= 0 ) cycle		! boundary
	    ians(ii) = iarv(ien)
	  end do
	  ic = count( ians == ia )
	  if( ic == 0 ) goto 99			! error - cannot be 0
	  if( ic > 1 ) cycle			! more than one of same color
	  ii1 = findloc(ians,ia,1)
	  if( ii1 == 0 ) goto 99		! error - cannot be 0
	  ii2 = 1 + mod(ii1,3)
	  ii3 = 1 + mod(ii2,3)
	  if( ians(ii2) /= ians(ii3) ) cycle	! other colors are different
	  ian = ians(ii2)
	  if( ian == -1 ) cycle			! do not handle boundaries
	  !write(6,*) 'regularize: ',ia,ians
	  iarv(ie) = ian
	  itot = itot + 1
	end do

	return
   99	continue
	write(6,*) 'error in element ie ',ie,ipev(ie)
	write(6,*) 'color: ',ia
	write(6,*) 'color of neigbors: ',ians
	stop 'error stop straigthen_edges: color error'
	end

!*****************************************************************

	subroutine write_grd_from bas(file)

	use basin

	implicit none

	character*(*) file

	integer inext(nkn)
	integer intype(nkn)
	integer ieext(nel)
	integer ietype(nel)
	character*80 text

	text = 'regularized'

        inext = ipv
        ieext = ipev

	intype = 0
	ietype = iarv

	call write_grd_file(file,text,nkn,nel,xgv,ygv,nen3v
     +                  ,inext,ieext,intype,ietype)

	end

!*****************************************************************

	subroutine write_box_info(ic,icolor)

	use basin

	implicit none

	integer ic
	integer icolor(ic)

	integer ntot,icol
	integer ncolumns,ncc
	integer ic1,ic2,ic3

	write(6,*) 'box number and number of elements: '
	ntot = 0
	do icol=1,ic
	  !write(6,*) icol,icolor(icol)
	  ntot = ntot + icolor(icol)
	end do
	write(6,*) 'total: ',nel,ntot

	ncolumns = 3
	ncc = ic / ncolumns
	if( ncc * ncolumns /= ic ) ncc = ncc + 1
	!write(6,*) 'ncc = ',ncc

	write(6,'(a)') '     box   color     box   color     box   color'
	do ic1=1,ncc
	  ic2 = ncc + ic1
	  ic3 = 2*ncc + ic1
	  if( ic3 <= ic ) then
	    write(6,1000) ic1,icolor(ic1),ic2,icolor(ic2),ic3,icolor(ic3)
	  else
	    write(6,1000) ic1,icolor(ic1),ic2,icolor(ic2)
	  end if
	end do

	return
 1000	format(6i8)
	end

!*****************************************************************

