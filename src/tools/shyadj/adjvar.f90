
!--------------------------------------------------------------------------
!
!    Copyright (C) 2003,2010,2013,2015,2018-2019  Georg Umgiesser
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

!  description :
! 
!  various utility routines for adj
! 
!  contents :
! 
!  subroutine smooth_grid(npass,omega,nkn,nel,nen3v,nbound,xgv,ygv,dx,dy,ic)
!                smoothing of internal nodes
!  subroutine chkgrd
!                check of grid
! 
!  revision log :
! 
!  01.01.2003	ggu	written
!  19.05.2003	ggu	some more utility routines
!  10.03.2010	ggu	area computation changed in checkarea (bug in 64 bit)
!  23.03.2010	ggu	changed v6.1.1
!  25.10.2013	ggu	changed VERS_6_1_68
!  19.01.2015	ggu	changed VERS_7_1_2
!  10.07.2015	ggu	changed VERS_7_1_50
!  24.07.2015	ggu	changed VERS_7_1_82
!  30.07.2015	ggu	changed VERS_7_1_83
!  12.10.2015	ggu	changed VERS_7_3_3
!  19.10.2015	ggu	changed VERS_7_3_6
!  18.12.2018	ggu	changed VERS_7_5_52
!  21.05.2019	ggu	changed VERS_7_5_62
!  23.02.2026   ggu     checks to avoid negative areas
!  27.02.2026   ggu     checks grade index for not unique nodes
! 
! **************************************************************

	subroutine smooth_grid(npass,omega)

!  smoothing of internal nodes

	use mod_adj_grade
	use basin

	implicit none

	integer npass
	real omega

	integer n,k,ie,ii,nks,ipos
	integer nos,nb
	integer ic(nkn)
	real dx(nkn)
	real dy(nkn)
	real xm,ym
	real amax
	integer nosmooth(nkn)

	write(6,*) 'smoothing grid '

	nosmooth = 0
	do k=1,nkn
	  nks = ngrade(k)
          call check_angles(k,nks,ngri(:,k),amax,ipos)
	  if( amax > 180. ) nosmooth(k) = 1
	end do
	nos = count(nosmooth==1)
	nb = count(nbound/=0)
	write(6,*) 'cannot smooth nodes: ',nos,nb,nkn

	call checkarea(0,'before smoothing')

	do n=1,npass

	  !if( mod(n,10) .eq. 0 ) write(6,*) 'smoothing... pass ',n
	  if( mod(n,1) .eq. 0 ) write(6,*) 'smoothing... pass ',n

!  initialize

	  do k=1,nkn
	    dx(k) = 0.
	    dy(k) = 0.
	    ic(k) = 0
	  end do

!  loop over elements -> accumulate

	  do ie=1,nel
	    xm = 0.
	    ym = 0.
	    do ii=1,3
	      k = nen3v(ii,ie)
	      xm = xm + xgv(k)
	      ym = ym + ygv(k)
	    end do
	    do ii=1,3
	      k = nen3v(ii,ie)
	      if( nbound(k) .eq. 0 .and. nosmooth(k) == 0 ) then
		ic(k) = ic(k) + 2
		dx(k) = dx(k) + 3. * xgv(k) - xm
		dy(k) = dy(k) + 3. * ygv(k) - ym
	      end if
	    end do
	  end do

!  change coordinate

	  do k=1,nkn
	    if( ic(k) .gt. 0 ) then
		xgv(k) = xgv(k) - omega * dx(k) / ic(k)
		ygv(k) = ygv(k) - omega * dy(k) / ic(k)
	    end if
	  end do

!  check area

	  call checkarea(0,'during smoothing')

	end do

	write(6,*) 'smoothing finished - passes ',npass
	call checkarea(0,'after smoothing')

	end

! **************************************************************

	subroutine mkstatic(nkn,ianv,nbound)

!  marks nodes as static (not moveable)

	use mod_adj_static

	implicit none

	integer nkn
	integer ianv(1)
        integer nbound(1)

        integer k

        do k=1,nkn
          if( ianv(k) .eq. iastatic .and. nbound(k) .eq. 0 ) then
	    write(6,*) '*** marking node ',k,' (internal) as static'
            nbound(k) = nbstatic
          end if
        end do

        end

! **************************************************************

	subroutine chkgrd(text)

!  check of grid
! 
!  checks consistency of node and grade index

	use mod_adj_grade
	use mod_adj_static
	use basin

	implicit none

	character*(*) text

	integer iaux(nkn)

	logical bstop,bverb
	integer n,k,ie,ii,kk,i,ke
	integer i1,i2,k1,k2
        integer nb,ic

	integer nextgr

	bverb = .true.
	bverb = .false.

	bstop = .false.

! --------------------------------------------------
!  check element index for unknown nodes
! --------------------------------------------------

	if( text .ne. ' ' ) write(6,*) 'chkgrd: '//trim(text)

        if( bverb ) write(6,*) 'checking for unknown nodes ...'

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( k .le. 0 .or. k .gt. nkn ) then
		bstop =.true.
		write(6,*) 'chkgrd (1): ',ie,ii,k
	    end if
	  end do
	end do

	if( bstop ) then
		call wr0grd
		write(6,*) 'nkn,nel: ',nkn,nel
		stop 'error stop chkgrd (1)'
	end if

! --------------------------------------------------
!  check grades for strange grades
! --------------------------------------------------

        if( bverb ) write(6,*) 'checking for strange grades ...'

	do k=1,nkn
	  n = ngrade(k)
	  if( n .le. 0 ) then
		bstop =.true.
		write(6,*) 'chkgrd (2): ',k,n
	  end if
	end do

	if( bstop ) then
		call wr0grd
		write(6,*) 'nkn,nel: ',nkn,nel
		stop 'error stop chkgrd (2)'
	end if

! --------------------------------------------------
!  check grade index for unknown nodes
! --------------------------------------------------

        if( bverb ) write(6,*) 'checking grade index for nodes ...'

	do k=1,nkn
	  n = ngrade(k)
	  do i=1,n
	    kk = ngri(i,k)
	    if( kk .le. 0 .or. kk .gt. nkn ) then
		bstop =.true.
		write(6,*) 'chkgrd (3): ',k,n,i,kk
	    end if
	  end do
	end do

	if( bstop ) then
		call wr0grd
		write(6,*) 'nkn,nel: ',nkn,nel
		stop 'error stop chkgrd (3)'
	end if

! --------------------------------------------------
!  check grade index for not unique nodes
! --------------------------------------------------

        if( bverb ) write(6,*) 'checking grade index for nodes ...'

	do k=1,nkn
	  n = ngrade(k)
	  do i=1,n
	    kk = ngri(i,k)
	    ic = count(kk==ngri(i+1:n,k))
	    if( ic > 0 ) then
		bstop =.true.
		write(6,*) 'chkgrd (7): ',k,n,i,kk
		write(6,*) 'chkgrd nodes: ',ngri(1:n,k)
	    end if
	  end do
	if( k == 2460 ) then
	  write(6,*) 'testing not unique: ',k,n
	  write(6,*) ngri(1:n,k)
	end if
	end do

	if( bstop ) then
		call wr0grd
		write(6,*) 'nkn,nel: ',nkn,nel
		stop 'error stop chkgrd (7): not unique nodes'
	end if

! --------------------------------------------------
!  check consistency of grade index
! --------------------------------------------------

        if( bverb ) write(6,*) 'consistency check  ...'

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    i1 = mod(ii,3) + 1
	    k1 = nen3v(i1,ie)
	    i2 = mod(i1,3) + 1
	    k2 = nen3v(i2,ie)

	    kk = nextgr(k,k1,ngrdi,ngrade,ngri)
	    if( kk .ne. k2 ) then
		bstop =.true.
		write(6,*) 'chkgrd (4): ',ie,k,k1,k2,kk
	    end if
	  end do
	end do

	if( bstop ) then
		call wr0grd
		write(6,*) 'nkn,nel: ',nkn,nel
		stop 'error stop chkgrd (4)'
	end if

! --------------------------------------------------
!  is grade still ok?
! --------------------------------------------------

        if( bverb ) write(6,*) 'checking for final grades ...'

	do k=1,nkn
	  iaux(k) = 0
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    iaux(k) = iaux(k) + 1
	  end do
	end do

	do k=1,nkn
          nb = nbound(k)
	  if( nb .ne. 0 .and. nb .ne. nbstatic ) then	!boundary node
	    if( iaux(k) + 1 .ne. ngrade(k) ) then
		write(6,*) 'chkgrd (5a): ',k,nbound(k),iaux(k),ngrade(k)
		ke = ipv(k)
		write(6,*) '   (problem in boundary node ',ke,' )'
		bstop =.true.
	    end if
	  else				!internal node
	    if( iaux(k) .ne. ngrade(k) ) then
		write(6,*) 'chkgrd (5b): ',k,nbound(k),iaux(k),ngrade(k)
		ke = ipv(k)
		write(6,*) '   (problem in internal node ',ke,' )'
		bstop =.true.
	    end if
	  end if
	end do

	if( bstop ) then
		call wr0grd
		write(6,*) 'nkn,nel: ',nkn,nel
		stop 'error stop chkgrd (5)'
	end if

! --------------------------------------------------
!  checking area
! --------------------------------------------------

	call checkarea(0,text)

! --------------------------------------------------
!  end of routine
! --------------------------------------------------

        if( bverb ) write(6,*) 'check ok ...'

	end

! *******************************************************

	subroutine nshift(ip,n,nga,iau)

!  shifts array - index ip will be first index
!  iau is auxiliary array

	implicit none

	integer ip,n
	integer nga(n),iau(n)

	integer i,ipa

	ipa = ip - 1 
	do i=1,n
	  ipa = mod(ipa,n) + 1
	  iau(i) = nga(ipa)
	end do

	do i=1,n
	  nga(i) = iau(i)
	end do

	end

! *******************************************************

	subroutine checkarea(k,text)

!  check if area is positive

	use basin

        implicit none

	integer k
	character*(*) text

	logical bstop
	integer ie,ii,i1,i2,k1,k2
        integer ieext
	real x1,y1,x2,y2,x3,y3
	real aj				!is twice the area
	character*80 string

	real areat

	bstop = .false.
	string = text

	do ie=1,nel
	  x1 = xgv(nen3v(1,ie))
	  x2 = xgv(nen3v(2,ie))
	  x3 = xgv(nen3v(3,ie))
	  y1 = ygv(nen3v(1,ie))
	  y2 = ygv(nen3v(2,ie))
	  y3 = ygv(nen3v(3,ie))
	  aj = areat(x1,y1,x2,y2,x3,y3)
	  if( aj .le. 0. ) then
            call eint2ext(ie,ieext)
	    write(6,*) 'element ',ie,' (extern ',ieext &
     &                          ,') has negative area...'
	    call elem_info(ie)
	    bstop = .true.
	  end if
        end do

	if( bstop ) then
	    write(6,*) 'error while checking: ',trim(string)
	    if( k /= 0 ) write(6,*) 'error while handling node ',k
	    call wr0grd
	    write(6,*) 'nkn,nel: ',nkn,nel
	    stop 'error stop checkarea'
	end if

	end

! *******************************************************

	subroutine elem_info(ie)

	use mod_adj_grade
	use basin

        implicit none

	integer ie,ieext,ii,k

        call eint2ext(ie,ieext)

	write(6,*) 'info on element ie = ',ie,' extern = ',ieext

	do ii=1,3
	  k = nen3v(ii,ie)
	  call node_info(k)
	end do

	if( bplot_error ) call plot_element(ie)

	end

! *******************************************************

	subroutine node_info(k)

!  writes info on node

	use mod_adj_grade
	use basin

        implicit none

	integer k,kext

	integer i

	if( k .le. 0 ) return

        call nint2ext(k,kext)

	write(6,*) 'info on node k = ',k,' extern = ',kext
	write(6,*) ngrade(k),nbound(k),xgv(k),ygv(k)
	write(6,*) nkn,nel,ngrdi
	write(6,*) (ngri(i,k),i=1,ngrade(k))

	end

! ****************************************************************
! ****************************************************************
! ****************************************************************

	subroutine make_grd_name(text,file)

	implicit none

	character*(*) text,file

	integer i
	integer, save :: igrd = 0
	character*5 num

	igrd = igrd + 1

	write(num,'(i5)') igrd
	num = adjustr(num)
	do i=1,5
	  if( num(i:i) == ' ' ) num(i:i) = '0'
	end do

	file = 'plot_' // text // '_' // num // '.grd'

	write(6,*) 'plotting to file ',trim(file)

	end 

! ****************************************************************

	subroutine make_unique(n,list)

	implicit none

	integer n
	integer list(n)

	integer i,nn

	call sort(n,list)

	nn = 1
	do i=2,n
	  if( list(i) == list(i-1) ) cycle
	  nn = nn + 1
	  if( i == nn ) cycle
	  list(nn) = list(i)
	end do

	n = nn

	end

! ****************************************************************

	subroutine plot_nodes(n,list)

	use basin

	implicit none

	integer n
	integer list(n)

	integer i,kk
	integer nn,ne
	integer nodes(3*nkn)
	integer elems(3*nkn)
	character*80 file

	nn = 0
	ne = 0
	do i=1,n
	  kk = list(i)
	  call add_nodes_and_elements(kk,nn,nodes,ne,elems)
	end do
	write(6,*) 'found total nodes and elems: ',nn,ne

	call make_unique(nn,nodes)
	call make_unique(ne,elems)
	write(6,*) 'found unique nodes and elems: ',nn,ne

	call make_grd_name('node',file)
	call write_partial_grid(file,nn,nodes,ne,elems)

	end

! ****************************************************************

	subroutine plot_node(kk)

	implicit none

	integer kk

	integer n
	integer list(1)

	n = 1
	list(1) = kk
	call plot_nodes(n,list)

	end

! ****************************************************************

	subroutine plot_elements(n,list)

	use basin

	implicit none

	integer n
	integer list(n)

	integer k,ii,ie,i
	integer nn,ne
	integer nodes(3*nkn)
	integer elems(3*nkn)
	character*80 file

	nn = 0
	ne = 0
	do i=1,n
	  ie = list(i)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    call add_nodes_and_elements(k,nn,nodes,ne,elems)
	  end do
	end do
	write(6,*) 'found total nodes and elems: ',nn,ne
	  
	call make_unique(nn,nodes)
	call make_unique(ne,elems)
	write(6,*) 'found unique nodes and elems: ',nn,ne

	call make_grd_name('elem',file)
	call write_partial_grid(file,nn,nodes,ne,elems)

	end

! ****************************************************************

	subroutine plot_element(ie)

	implicit none

	integer ie

	integer n
	integer list(1)

	n = 1
	list(1) = ie
	call plot_elements(n,list)

	end

! ****************************************************************

	subroutine write_partial_grid(file,nn,nodes,ne,elems)

	use basin

	implicit none

	character*(*) file
	integer nn,ne
	integer nodes(nn)
	integer elems(ne)

	integer i,k,ie,ii
	integer naux(3)

	open(1,file=file,status='unknown',form='formatted')

	do i=1,nn
	  k = nodes(i)
	  write(1,1000) 1,ipv(k),0,xgv(k),ygv(k)
	end do

	do i=1,ne
	  ie = elems(i)
	  do ii=1,3
	    naux(ii) = ipv(nen3v(ii,ie))
	  end do
	  write(1,2000) 2,ipev(ie),0,3,naux(:)
	end do

	close(1)

	return
 1000	format(i1,i8,i5,2f14.5)
 2000	format(i1,i8,2i5,3i8)
	end

! ****************************************************************

	subroutine add_nodes_and_elements(k,nn,nodes,ne,elems)

! adds nodes and elements to list

	use basin

	implicit none

	integer k,nn,ne
	integer nodes(nkn)
	integer elems(nel)

	integer ie,ii,kk

	do ie=1,nel
	  do ii=1,3
	    kk = nen3v(ii,ie)
	    if( kk == k ) then
	      nodes(nn+1:) = nen3v(:,ie)
	      nn = nn + 3
	      ne = ne + 1
	      elems(ne) = ie
	    end if
	  end do
	end do

	end

! ****************************************************************

	subroutine find_nodes_and_elements(k,nn,nodes,ne,elems)

! creates nodes and elements list (list is initialized to 0)

	use basin

	implicit none

	integer k,nn,ne
	integer nodes(nkn)
	integer elems(nel)

	integer ie,ii,kk

	nn = 0
	ne = 0
	do ie=1,nel
	  do ii=1,3
	    kk = nen3v(ii,ie)
	    if( kk == k ) then
	      nodes(nn+1:) = nen3v(:,ie)
	      nn = nn + 3
	      ne = ne + 1
	      elems(ne) = ie
	    end if
	  end do
	end do

	end

! ****************************************************************

	function rangle(k1,k2,k3)

!  gives angle between nodes
! 
!  rangle	angle [degrees] ,     rangle < 180 => right turn
!  k1,k2,k3	node numbers

	use basin

	implicit none

	real rangle
	integer k1,k2,k3

	real x1,y1,x2,y2,x3,y3
	real angle

	x1 = xgv(k1)
	y1 = ygv(k1)
	x2 = xgv(k2)
	y2 = ygv(k2)
	x3 = xgv(k3)
	y3 = ygv(k3)

	rangle = angle(x1,y1,x2,y2,x3,y3)

	end

! ***********************************************************

        subroutine check_angles(k,n,list,amax,ipos)

! computes maximum angle of node list around k (must be sorted)

        implicit none

        integer k		!central node
        integer n
        integer list(n)
        real amax		!maximum angle (return)
        integer ipos		!position of maximum angle in list (return)

        integer i,ibefore,iafter
        real a1,a2
        real angle(n)

        real rangle

        ipos = 0
        amax = 0.
        do i=1,n
          ibefore = i - 1
          if( ibefore == 0 ) ibefore = n
          iafter = i + 1
          if( iafter > n ) iafter = 1
          a1 = rangle(k,list(i),list(iafter))
          if( a1 > 180.) a1 = 360. - a1
          a2 = rangle(k,list(i),list(ibefore))
          if( a2 > 180.) a2 = 360. - a2
          angle(i) = a1 + a2
          if( angle(i) > amax ) then
            amax = angle(i)
            ipos = i
          end if
        end do

        !write(6,*) k,list
        !write(6,*) k,angle

        end

! ************************************************************

        subroutine nint2ext(kint,kext)

	use mod_adj_grade
	use basin

        implicit none

        integer kint,kext

        kext = ipv(kint)

        end

! ************************************************************

        subroutine eint2ext(ieint,ieext)

	use mod_adj_grade
	use basin

        implicit none

        integer ieint,ieext

        ieext = ipev(ieint)

        end

! ************************************************************

