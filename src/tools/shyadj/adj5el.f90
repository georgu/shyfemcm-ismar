
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

!  description :
! 
!  5-5 grade routines
! 
!  contents :
! 
!  subroutine elim5(nkn,nel,ngrddi,ngrade,nbound,ngri,nen3v)
! 			eliminates low grades
!  subroutine elim55(k,nkn,nel,ngrddi,ngrade,nbound,ngri,nen3v)
! 			eliminates 5-5 connections
! 
!  revision log :
! 
!  01.01.2003	ggu	written
!  23.03.2010	ggu	changed v6.1.1
!  19.01.2015	ggu	changed VERS_7_1_2
!  21.05.2015	ggu	changed VERS_7_1_11
!  10.07.2015	ggu	changed VERS_7_1_50
!  24.07.2015	ggu	changed VERS_7_1_82
!  30.07.2015	ggu	changed VERS_7_1_83
!  11.10.2015	ggu	bug fix: fused node was not moved
!  18.12.2018	ggu	changed VERS_7_5_52
!  21.05.2019	ggu	changed VERS_7_5_62
!  23.02.2026   ggu     checks to avoid negative areas
!  27.02.2026   ggu     some more checks
! 
! ***********************************************************

	subroutine elim_5_5

!  eliminates 5-5 grades
! 
!  the two nodes with 5-5 are fused together
!  the two elements attached to these nodes are deleted

	use mod_adj_grade
	use basin

	implicit none

        integer k,n
	integer ielim
	logical bok

        write(6,*) 'eliminating grades for grade 5... '

	ielim = 0

        do k=1,nkn
          n = ngrade(k)
          if( n .eq. 5 .and. nbound(k) .eq. 0 ) then
            call elim55(k,bok)
	    if( bok ) ielim = ielim + 1
	    if( bok ) call chkgrd('checking inside 5 grade')
	  end if
        end do

        write(6,*) ielim, ' nodes have been changed'

	end

! ***********************************************************

	subroutine elim55(k,bok)

!  eliminates 5-5 connections

	use mod_adj_grade
	use basin

	implicit none

	integer k
	logical bok

	logical berr
        integer n,i,nc,nmax,ii,ks,nks
	integer ie1,ie2
	integer ip1,ip2
	integer np,nt,nn
	integer nval,ip
	integer ipos,ng,nb
	real amax
	integer ngav(0:ngrdi+1)
	integer ngrv(0:ngrdi+1)
	integer nbav(0:ngrdi+1)
	integer naux(0:ngrdi+1)
	integer iplist(ngrdi)

	integer ifindel

	berr = .false.
	bok = .false.

	if( k .gt. nkn ) return

	if( bdebug ) then
	  write(6,*) '==============================================='
	  write(6,*) 'debug of new node: ',k
	end if

!  make circular list
! 
!  ngav 	node numbers around k
!  ngrv		grades of node numbers around k
!  nbav		boundary flag for nodes around k

        n = ngrade(k)
	ngav(0) = ngri(n,k)
	do i=1,n
	  ngav(i) = ngri(i,k)
	end do
	ngav(n+1) = ngri(1,k)

	do i=0,n+1
	  ngrv(i) = ngrade(ngav(i))
	  nbav(i) = 0
	  if( nbound(ngav(i)) .ne. 0 ) then	!boundary node
	    ngrv(i) = 6				!fake perfect grade
	    nbav(i) = 1
	  end if
	end do

!  check if exchange is possible

	nc = 0		!how many of this nmax value
	nmax = 0	!maximum sum of grades -> must be at least 3
	ip = 0		!pointer to node in list that has been chosen
	do i=1,n
	  np = ngrv(i-1)
	  nt = ngrv(i)
	  nn = ngrv(i+1)

	  nval = np + nn - n - nt

	  if( nval .gt. nmax ) then
	    nc = 1
	    iplist(nc) = i
	    nmax = nval
	  else if( nval .eq. nmax ) then
	    nc = nc + 1
	    iplist(nc) = i
	  end if
	end do

	if( nmax .lt. 3 ) return

	if( nc > 1 ) then
	  do i=1,nc
	    ip = iplist(i)
	    ks = ngav(ip)
	    ng = ngrv(ip)
	    nb = nbav(ip)
	    naux(1:ng) = ngri(1:ng,ks)
	    write(6,*) '5-5 before: ',ks,ng,nb
	    write(6,*) '5-5 naux: ',naux(1:ng)
	    call check_angles(ks,ng,naux(1:n),amax,ipos)
	    write(6,*) '5-5 amax: ',ks,amax
	  end do
	end if
!        do i=1,n
!          neibs(i) = ngri(i,k)
!          ngneib(i) = ngrade(neibs(i))
!        end do

	ip = 0
	do i=1,nc
	  ip = iplist(i)
	  if( nbav(ip) == 0 ) exit	!take the first non boundary node
	end do

	if( i > nc ) return		!no possible node

	call check_angles(k,n,ngav(1:n),amax,ipos)
	if( amax > 180 ) then
	  write(6,*) 'cannot eliminate... k angle > 180'
	  return
	end if

	write(6,*) k,n,nmax,nc,ip,amax

	!if( k == 2821 ) call plot_node(k)
	!if( k == 2834 ) berr = .true.

!  nc gives number of occurences of this value of nmax ...
!  ip is the pointer to the node to be exchanged
!  we know that it is not a boundary node, so we can shift it
! 
!  k is eliminated, ngav(ip) is retained

	ks = ngav(ip)		! node to be shifted
	nks = ngrade(ks)
	!if( nks /= 5 ) stop 'error stop: nks/=5'
	call check_angles(ks,nks,ngri(:,ks),amax,ipos)
	if( amax > 180 ) then
	  write(6,*) 'cannot eliminate... ks angle > 180'
	  return
	end if

	if( berr ) then
	  write(6,*) 'writing grid for k=2834'
	  call plot_node(k)
	  call plot_node(ks)
	  call plot_nodes(2,(/k,ks/))
	end if

	if( bdebug ) then
	    write(6,*) 'exchanging with node ... ',ks
	    write(6,'(7i10)') (ngav(i),i=0,n+1)
	    write(6,'(7i10)') (ngrv(i),i=0,n+1)
	    write(6,'(7i10)') (nbav(i),i=0,n+1)
	    call plosno(k)
	    call plosno(ngav(ip))
	end if

!  find elements that have to be deleted

	ie1 = ifindel(k,ngav(ip),ngav(ip+1))
	ie2 = ifindel(k,ngav(ip-1),ngav(ip))

	if( ie1 .eq. 0 .or. ie2 .eq. 0 ) then
	  stop 'error stop elim55: internal error (2)'
	end if

	if( bdebug ) then
	  write(6,*) 'elements to be deleted... ',ie1,ie2
	  write(6,*) ie1,k,ngav(ip),ngav(ip+1)
	  write(6,*) (nen3v(ii,ie1),ii=1,3)
	  write(6,*) ie2,k,ngav(ip-1),ngav(ip)
	  write(6,*) (nen3v(ii,ie2),ii=1,3)
	  call plosel2(ie1,ie2)
	end if

!  delete elements

	if( ie1 .gt. ie2 ) then		!to avoid bug
	  call delele(ie1)
	  call delele(ie2)
	else
	  call delele(ie2)
	  call delele(ie1)
	end if

	if( bdebug ) then
	  write(6,*) 'grade index befor manipulation:'
	  call prgr(k,ngrdi,ngrade,ngri)
	  call prgr(ngav(ip),ngrdi,ngrade,ngri)
	  call prgr(ngav(ip-1),ngrdi,ngrade,ngri)
	  call prgr(ngav(ip+1),ngrdi,ngrade,ngri)
	end if

!  new coordinates for node

	xgv(ks) = 0.5 * ( xgv(k) + xgv(ks) )
	ygv(ks) = 0.5 * ( ygv(k) + ygv(ks) )

!  substitute all occurrences of k with ks

	call subnod(k,ks)

	if( bdebug ) then
	  write(6,*) 'after substitution...'
	  call prgr(ngav(ip),ngrdi,ngrade,ngri)
	  call prgr(ngav(ip-1),ngrdi,ngrade,ngri)
	  call prgr(ngav(ip+1),ngrdi,ngrade,ngri)
	end if

!  adjourn grade (delete) for nodes ip-1, ip+1

	call delgr(ngav(ip-1),ngav(ip),ngrdi,ngrade,ngri)
	call delgr(ngav(ip+1),ngav(ip),ngrdi,ngrade,ngri)

	if( bdebug ) then
	  write(6,*) 'after deleting ip-1,ip+1...'
	  call prgr(ngav(ip),ngrdi,ngrade,ngri)
	  call prgr(ngav(ip-1),ngrdi,ngrade,ngri)
	  call prgr(ngav(ip+1),ngrdi,ngrade,ngri)
	end if

!  adjourn grade index for ip and delete node k finally

	call delgr(ngav(ip),ngav(ip),ngrdi,ngrade,ngri)
	call delnod(k)
	call subval(n+2,ngav(0),nkn+1,k)	!if nkn is in ngav

	if( bdebug ) then
	  write(6,*) 'after deleting ip...'
	  call prgr(ngav(ip),ngrdi,ngrade,ngri)
	end if

	ip1 = mod(ip+2,n)
	ip2 = mod(ip+3,n)
	call insgr(ngav(ip),ngav(ip+1),ngav(ip1),ngrdi,ngrade,ngri)
	call insgr(ngav(ip),ngav(ip1),ngav(ip2),ngrdi,ngrade,ngri)

	if( bdebug ) then
	  write(6,*) 'grade index after manipulation:'
	  call prgr(ngav(ip),ngrdi,ngrade,ngri)
	  call prgr(ngav(ip-1),ngrdi,ngrade,ngri)
	  call prgr(ngav(ip+1),ngrdi,ngrade,ngri)
	end if

	if( bdebug ) then
	  call plosno(ngav(ip))
	  call plosno(ngav(ip-1))
	  call plosno(ngav(ip+1))
	end if

	if( bdebug ) then
	  write(6,*) 'end of debug of node ',k
	  write(6,*) '==============================================='
	end if

	if( berr ) write(6,*) 'calling checkarea'
	if( bcheck ) call checkarea(k,' ')
	if( berr ) write(6,*) 'finished calling checkarea'
	! should only check elements around ks -> we need element index

	bok = .true.

	end

! ***********************************************************

