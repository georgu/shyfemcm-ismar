
!--------------------------------------------------------------------------
!
!    Copyright (C) 2003,2010,2012,2015,2018-2019  Georg Umgiesser
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
!  30.03.2012	ggu	changed VERS_6_1_51
!  19.01.2015	ggu	changed VERS_7_1_2
!  21.05.2015	ggu	changed VERS_7_1_11
!  24.07.2015	ggu	changed VERS_7_1_82
!  30.07.2015	ggu	changed VERS_7_1_83
!  12.10.2015	ggu	changed VERS_7_3_3
!  18.12.2018	ggu	changed VERS_7_5_52
!  21.05.2019	ggu	changed VERS_7_5_62
!  23.02.2026   ggu     checks to avoid negative areas
!  06.03.2026   ggu     completely restructured
! 
!  description :
! 
!  5-7-5 grade routines
! 
!  contents :
! 
!  subroutine elim57(nkn,nel,ngrddi,ngrade,nbound,ngri,nen3v)
! 			eliminates special 575 connections (shell)
!  subroutine elim575(k,nkn,nel,ngrddi,ngrade,nbound,ngri,nen3v)
! 			eliminates 5-7-5 connections
! 
! ***********************************************************

	subroutine elim_5_7_5

!  eliminates 5-7-5 grades

	use mod_adj_grade
        use mod_progress_bar
	use basin

	implicit none

	logical bok,bprog
        integer k,n
	integer itot,ielim
	real perc
	character*10 string

        bprog = bquiet .and. .not. bsilent
        bprog = .not. bverbose .and. .not. bsilent

        if( .not. bquiet ) write(6,*) 'eliminating 5-7-5 grades... '
	if( bprog ) call progress_bar_init('elim575')

	itot = 0
	ielim = 0

        do k=1,nkn
          perc = k/float(nkn)
          if( bprog ) call progress_bar_print(k,nkn)
          n = ngrade(k)
          if( n .eq. 7 .and. nbound(k) .eq. 0 ) then
            call elim575(k,bok)
            itot = itot + 1
            if( bok ) ielim = ielim + 1
	    !call chkgrd('checking in elim575')
          end if
        end do

        if( bprog ) call progress_bar_finalize

	end

! ***********************************************************

	subroutine elim575(k,bok)

!  eliminates 5-7-5 connections
! 
!  create one new node and two new elements

	use mod_adj_grade
	use basin

	implicit none

	integer k
	logical bok

        integer n,i,nc,ii
	integer ie,kk,iii,ks,nks
	integer ip1,ip2,ipos
	integer ip
	integer ng,idp
	integer ngav(ngrdi)	!we do not need 0 index
	integer ngrv(ngrdi)
	integer nbav(ngrdi)
	integer iau(ngrdi)
	real x,y,xm,ym
	real amax

	bok = .false.

	if( k .gt. nkn ) return

	if( bdebug ) write(6,*) 'elim575 new node: ',k

!  make list

        n = ngrade(k)
	if( n > ngrdi ) stop 'error stop elim575: ngrdi'
	do i=1,n
	  ngav(i) = ngri(i,k)
	end do

	do i=1,n
	  ngrv(i) = ngrade(ngav(i))
	  nbav(i) = 0
	  if( nbound(ngav(i)) .ne. 0 ) then
	    ngrv(i) = 6
	    nbav(i) = 1
	  end if
	end do

!  check if exchange is possible

	nc = 0
	do i=1,n
	  ng = ngrv(i)
	  if( ng .eq. 5 ) nc = nc + 1
	end do

	if( nc .ne. 2 ) return	!if not exactly 2 cannot proceed

        ks = k
        nks = ngrade(ks)
        call check_angles(ks,nks,ngri(:,ks),amax,ipos)
        if( amax > 180 ) then
          if( bverbose ) write(6,*) 'cannot eliminate... ks angle > 180: ',k
          return
        end if

!  find out distance of 5 grades

	nc = 0
	ip1 = 0
	ip2 = 0
	do i=1,n
	  ng = ngrv(i)
	  if( ng .eq. 5 ) then
	    nc = nc + 1
	    if( nc .eq. 1 ) then
		ip1 = i
	    else
		ip2 = i
	    end if
	  end if
	end do

	if( bdebug ) then
	  do i=1,n
	    write(6,*) ngav(i),ngrv(i),nbav(i)
	  end do
	end if

	idp = ip2 - ip1
	if( idp .le. 2 .or. idp .ge. 5 ) return

	if( bverbose ) write(6,*) 'elim575: ',k,ip1,ip2,idp

	bok = .true.

	if( bdebug ) then
	  write(6,*) ngav(ip1),ngav(ip2)
	  write(6,'(7i10)') (ngav(i),i=1,7)
	  write(6,'(7i10)') (ngrv(i),i=1,7)
	  write(6,'(7i10)') (nbav(i),i=1,7)
	end if

!  reorder node list
!  node 1 is a 5-grade, and node 5 is a 5-grade

	ip = ip1
	if( idp .eq. 3 ) ip = ip2
	call nshift(ip,n,ngav,iau)
	call nshift(ip,n,ngrv,iau)
	call nshift(ip,n,nbav,iau)

	if( bdebug ) then
	  write(6,'(7i10)') (ngav(i),i=1,7)
	  write(6,'(7i10)') (ngrv(i),i=1,7)
	  write(6,'(7i10)') (nbav(i),i=1,7)
	end if

!  new node 

	call newnod(nkn)

!  substitute new node for old one in node index

	do ie=1,nel
	  do ii=1,3
	    kk = nen3v(ii,ie)
	    if( kk .eq. k ) then
	      do iii=1,3
		if( nen3v(iii,ie) .eq. ngav(2) ) nen3v(ii,ie) = nkn
		if( nen3v(iii,ie) .eq. ngav(4) ) nen3v(ii,ie) = nkn
	      end do
	    end if
	  end do
	end do

!  new elements

	call newele(nel)
	call setele(nel,ngav(1),nkn,k,nen3v)
	call newele(nel)
	call setele(nel,ngav(5),k,nkn,nen3v)

!  adjust grade index of old node (5 grade)

	call delgr(k,ngav(2),ngrdi,ngrade,ngri)
	n = ngrade(k)
	call delgr(k,ngav(3),ngrdi,ngrade,ngri)
	n = ngrade(k)
	call delgr(k,ngav(4),ngrdi,ngrade,ngri)
	n = ngrade(k)
	call insgr(k,ngav(1),nkn,ngrdi,ngrade,ngri)
	n = ngrade(k)

!  adjust grade index of new node (6 grade)

	ngri(:,nkn) = 0
	do i=1,5
	  ngri(i,nkn) = ngav(i)
	end do
	n = 6
	ngri(n,nkn) = k
	ngrade(nkn) = n
	call resort_index(n,ngri(:n,nkn))

!  adjust grade index of 5-5 nodes

	n = ngrade(k)
	call insgrb(ngav(1),k,nkn,ngrdi,ngrade,ngri)
	n = ngrade(k)
	call insgr(ngav(5),k,nkn,ngrdi,ngrade,ngri)
	n = ngrade(k)

!  substitute new node in grade index of nodes close to new node

	call exchgr(ngav(2),k,nkn,ngrdi,ngrade,ngri)
	call exchgr(ngav(3),k,nkn,ngrdi,ngrade,ngri)
	call exchgr(ngav(4),k,nkn,ngrdi,ngrade,ngri)

	if( bdebug ) then
	  call prgr(k,ngrdi,ngrade,ngri)
	  call prgr(nkn,ngrdi,ngrade,ngri)
	  call prgr(ngav(1),ngrdi,ngrade,ngri)
	  call prgr(ngav(5),ngrdi,ngrade,ngri)
	end if

!  adjust coordinates

	xm = 0.5 * ( xgv(ngav(6)) + xgv(ngav(7)) )
	ym = 0.5 * ( ygv(ngav(6)) + ygv(ngav(7)) )
	x = xgv(ngav(3))
	y = ygv(ngav(3))

	xgv(k) = xm + (1./3.) * ( x - xm )
	ygv(k) = ym + (1./3.) * ( y - ym )
	xgv(nkn) = xm + (2./3.) * ( x - xm )
	ygv(nkn) = ym + (2./3.) * ( y - ym )

	end

! *******************************************************

	subroutine node_debug(k,nkn,nel,nen3v,xgv,ygv)

	integer nen3v(3,nel)
	real xgv(nkn),ygv(nkn)

	iu = 79

	write(iu,*) k,nkn,nel
	do i=1,nel
	  do ii=1,3
	    write(iu,*) nen3v(ii,i)
	  end do
	end do
	do i=1,nkn
	  write(iu,*) xgv(i),ygv(i)
	end do

	end

! *******************************************************

