
!--------------------------------------------------------------------------
!
!    Copyright (C) 2003,2010,2015,2018-2019  Georg Umgiesser
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
!  revision log :
! 
!  01.01.2003	ggu	written
!  23.03.2010	ggu	changed v6.1.1
!  19.01.2015	ggu	changed VERS_7_1_2
!  24.07.2015	ggu	changed VERS_7_1_82
!  30.07.2015	ggu	changed VERS_7_1_83
!  12.10.2015	ggu	changed VERS_7_3_3
!  18.12.2018	ggu	changed VERS_7_5_52
!  21.05.2019	ggu	changed VERS_7_5_62
!  23.02.2026	ggu	checks to avoid negative areas
!  06.03.2026	ggu	completely restructured
! 
!  description :
! 
!  4 grade routines
! 
!  contents :
! 
!  subroutine elim4(nkn,nel,ngrddi,ngrade,ngri,nen3v)	eliminates grade=4 nodes
!  subroutine elim4n(k,nel,ngrddi,ngrade,ngri,nen3v)	eliminates node
!  subroutine el4to2(k,k1,k2,ieind,ienew)		four elements to two
!  subroutine unifel(k,ik1,ik2,inew)			unifies two elements
! 
! ***********************************************************

	subroutine elimlow

!  eliminates grade=4 nodes (and less)

	use mod_adj_grade
	use mod_progress_bar
	use basin

	implicit none

	logical bok,bprog
	integer k,n,nc
	integer itot,ielim
	real perc

	bprog = bquiet .and. .not. bsilent
	bprog = .not. bverbose .and. .not. bsilent

	if( .not. bquiet ) write(6,*) 'eliminating 3 grades...'
	if( bprog ) call progress_bar_init('elim3')

	itot = 0
	ielim = 0

	do k=1,nkn
          perc = k/float(nkn)
          if( bprog ) call progress_bar_print(k,nkn)
	  if( nbound(k) .eq. 0 ) then
	    n = ngrade(k)
	    if( n .eq. 3 ) then
	      call elim3(k,bok)
	      itot = itot + 1
	      if( bok ) ielim = ielim + 1
	    end if
	  end if
	end do

	if( bprog ) call progress_bar_finalize
	if( .not. bsilent ) write(6,*) 'nodes eliminated: ',ielim,' of ',itot

	if( .not. bquiet ) write(6,*) 'eliminating 4 grades...'
	if( bprog ) call progress_bar_init('elim4')

	itot = 0
	ielim = 0

	do k=1,nkn
          perc = k/float(nkn)
          if( bprog ) call progress_bar_print(k,nkn)
	  if( nbound(k) .eq. 0 ) then
	    n = ngrade(k)
	    if( n .eq. 4 ) then
	      call elim4(k,bok)
	      itot = itot + 1
	      if( bok ) ielim = ielim + 1
	    end if
	  end if
	end do

	if( bprog ) call progress_bar_finalize
	if( .not. bsilent ) write(6,*) 'nodes eliminated: ',ielim,' of ',itot

	end

! ***********************************************************

	subroutine elim3(k,bok)

!  eliminates node and all attached elements

	use mod_adj_grade
	use basin

	implicit none

	integer k
	logical bok

	integer i,n,kk,ie
	integer neibs(ngr),ngneib(ngr)

	integer ifindel

	if( k .gt. nkn ) return

	n = ngrade(k)
	if( n .ne. 3 ) stop 'error stop elim3: not grade 3'

	bok = .true.
	if( bverbose ) write(6,*) 'elim3 ',k,n

	do i=1,n
	  neibs(i) = ngri(i,k)
	  ngneib(i) = ngrade(neibs(i))
	end do

!  delete grade from neibors

	do i=1,n
	  kk = neibs(i)
	  call delgr(kk,k,ngrdi,ngrade,ngri)
	end do

!  delete elements

	ie = ifindel(neibs(1),neibs(2),k)
	if( ie .eq. 0 ) stop 'error stop elim3: internal error (1)'
	call setele(ie,neibs(1),neibs(2),neibs(3),nen3v)

	ie = ifindel(neibs(2),neibs(3),k)
	if( ie .eq. 0 ) stop 'error stop elim3: internal error (1)'
	call delele(ie)

	ie = ifindel(neibs(3),neibs(1),k)
	if( ie .eq. 0 ) stop 'error stop elim3: internal error (1)'
	call delele(ie)

!  delete node

	call delnod(k)

	end

! ***********************************************************

	subroutine elim4(k,bok)

!  eliminates node

	use mod_adj_grade
	use basin

	implicit none

	integer k
	logical bok

	integer i,n,k1,k2
	integer ipos,ipos1,ipos2,iposa
	integer neibs(ngr),ngneib(ngr)

	integer ie,ii,ip
	integer kk
	integer ielem(4), ieind(3,4), ienew(3,2)
	real amax

	if( k .gt. nkn ) return

	bok = .false.

	n = ngrade(k)

	do i=1,n
	  neibs(i) = ngri(i,k)
	  ngneib(i) = ngrade(neibs(i))
	end do

	call check_angles(k,n,neibs,amax,iposa)
	if( amax < 180 ) iposa = 0
	if( iposa > 2 ) iposa = iposa - 2

!  we have two solutions -> eliminate 1/3 or 2/4 connection
! 	can eliminate only if grade on both nodes is at least 6

	ipos1 = 0
	if( ngneib(1) .ge. 6 .and. ngneib(3) .ge. 6 ) ipos1 = 1
	ipos2 = 0
	if( ngneib(2) .ge. 6 .and. ngneib(4) .ge. 6 ) ipos2 = 1

	if( ipos1*ipos2 .eq. 0 ) then	!at least one not possible
	  if( ipos1+ipos2 .eq. 0 ) then	!none possible
		ipos = 0
	  else				!just one possible
		if( ipos1 .gt. 0 ) then
			ipos = 1
		else
			ipos = 2
		end if
	  end if
	  if( iposa > 0 .and. ipos /= iposa ) ipos = 0	!cannot eliminate
	  if( iposa > 0 ) ipos = 0	!cannot eliminate
	else				!both possible
	  ipos1 =  ngneib(1) + ngneib(3)
	  ipos2 =  ngneib(2) + ngneib(4)
	if( iposa > 0 ) write(6,*) 'B ipos1,ipos2: ',ipos1,ipos2
	  if( ipos1 .eq. ipos2 ) then	!most equilibrated: 7/7 better than 6/8
		ipos1 = ngneib(1) * ngneib(3)
		ipos2 = ngneib(2) * ngneib(4)
	  end if
	if( iposa > 0 ) write(6,*) 'C ipos1,ipos2: ',ipos1,ipos2
	  if( iposa .gt. 0 ) then
		ipos = iposa
		ipos = 0
	  else if( ipos1 .gt. ipos2 ) then
		ipos = 1
	  else
		ipos = 2
	  end if
	end if

	if( ipos .eq. 0 ) then
		if( bverbose ) write(6,*) 'Cannot eliminate node: ',k
	else if( ipos1 .eq. ipos2 ) then
! 		write(6,*) 'Both solutions equivalent; first chosen'
	else
! 		write(6,*) 'Best solutions: ',ipos
	end if

	if( bdebug ) then
	  write(6,*) 'poss: ',ipos1,ipos2
	  write(6,*) 'neibs: ',(neibs(i),i=1,n)
	  write(6,*) 'grade: ',(ngneib(i),i=1,n)
	end if

	if( ipos .eq. 0 ) return

!  now we have the information -> eliminate node

	bok = .true.

	if( bverbose ) write(6,*) 'elim4: ',k,n,ipos

!  get element numbers and index -> ielem, ieind

	i = 0
	do ie=1,nel
	  do ii=1,3
	    kk = nen3v(ii,ie)
	    if( kk .eq. k ) then
		i = i + 1
		ielem(i) = ie
	    end if
	  end do
	end do
	if( i .ne. 4 ) then
	  write(6,*) k,i
	  write(6,*) (ielem(k),k=1,i)
	  stop 'error stop elim4: internal error (2)'
	end if
	do i=1,4
	  ie = ielem(i)
	  do ii=1,3
	    ieind(ii,i) = nen3v(ii,ie)
	  end do
	end do

!  make new elements (locally)

	k1 = neibs(ipos)
	k2 = neibs(ipos+2)
	call el4to2(k,k1,k2,ieind,ienew)

	if( bdebug ) then
	write(6,*) 'old elements:'
	do i=1,4
	  write(6,*) ielem(i),(ieind(ii,i),ii=1,3)
	end do
	write(6,*) 'new elements:'
	do i=1,2
	  write(6,*) (ienew(ii,i),ii=1,3)
	end do
	end if

!  copy new elements -> still need to copy hev,...

	do i=1,2
	  ie = ielem(i)
	  call setele(ie,ienew(1,i),ienew(2,i),ienew(3,i),nen3v)
	end do

!  nodes that change grade -> adjust

	do ip=ipos,4,2
	  kk = neibs(ip)
	  call delgr(kk,k,ngrdi,ngrade,ngri)
	end do

!  nodes that do not change grade -> exchange information

	ipos1 = mod(ipos,2)+1
	ipos2 = ipos1 + 2
	k1 = neibs(ipos1)
	k2 = neibs(ipos2)

	call exchgr(k1,k,k2,ngrdi,ngrade,ngri)
	call exchgr(k2,k,k1,ngrdi,ngrade,ngri)

!  eliminate other two elements
!  here we must count downward, otherwise, if the last element
!  to be eliminated is number nel, we will have a bug

	do i=4,3,-1
	  ie = ielem(i)
	  call delele(ie)
	end do

!  eliminate node

	call delnod(k)

	end
	
! ***********************************************************

	subroutine el4to2(k,k1,k2,ieind,ienew)

!  four elements to two

	implicit none

	integer k,k1,k2
	integer ieind(3,4), ienew(3,2)

	integer iip,i,ii
	integer ip(4)

	iip = 0
	do i=1,4
	  do ii=1,3
	    if( ieind(ii,i) .eq. k1 ) then
		iip = iip + 1
		ip(iip) = i
	    end if
	  end do
	end do
	if( iip .ne. 2 ) then
		stop 'error stop el4to2: internal error (5)'
	end if

	call unifel(k,ieind(1,ip(1)),ieind(1,ip(2)),ienew(1,1))

	iip = 0
	do i=1,4
	  do ii=1,3
	    if( ieind(ii,i) .eq. k2 ) then
		iip = iip + 1
		ip(iip) = i
	    end if
	  end do
	end do
	if( iip .ne. 2 ) then
		stop 'error stop el4to2: internal error (5)'
	end if

	call unifel(k,ieind(1,ip(1)),ieind(1,ip(2)),ienew(1,2))

	end

! ***********************************************************

	subroutine unifel(k,ik1,ik2,inew)

!  unifies two elements with common node k

	implicit none

	integer k
	integer ik1(3),ik2(3),inew(3)

	integer ikp1,ikp2,in1,in2,ip1,ip2
	integer ii

!  first find node k in elements

	ikp1 = 0
	ikp2 = 0
	do ii=1,3
	  if( ik1(ii) .eq. k ) ikp1 = ii
	  if( ik2(ii) .eq. k ) ikp2 = ii
	end do

	if( ikp1 .eq. 0 .or. ikp2 .eq. 0 ) then
	  write(6,*) 'unifel : ',k
	  write(6,*) (ik1(ii),ii=1,3)
	  write(6,*) (ik2(ii),ii=1,3)
	  stop 'error stop unifel: no node k'
	end if

	in1 = mod(ikp1,3) + 1
	ip1 = mod(ikp1+1,3) + 1
	in2 = mod(ikp2,3) + 1
	ip2 = mod(ikp2+1,3) + 1
	if( ik1(in1) .eq. ik2(ip2) ) then
	  inew(1) = ik1(ip1)
	  inew(2) = ik2(in2)
	  inew(3) = ik1(in1)	!common node
	else if( ik1(ip1) .eq. ik2(in2) ) then
	  inew(1) = ik2(ip2)
	  inew(2) = ik1(in1)
	  inew(3) = ik1(ip1)	!common node
	else
	  write(6,*) (ik1(ii),ii=1,3)
	  write(6,*) (ik2(ii),ii=1,3)
	  stop 'error stop unifel: no second common node'
	end if

	end

! ***********************************************************

