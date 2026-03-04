
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
!  grading routines
! 
!  contents :
! 
!  subroutine maxgrd(nkn,nel,nen3v,ngrade,nmax)
! 			computes maximum grade (except boundary)
!  subroutine statgrd(nkn,nmax,nsg,ngrade)
! 			computes statistics on grade
!  subroutine mkbound(nkn,nel,ngrddi,nen3v,ngrade,ngri)
! 			marks boundary nodes and determines exact grade
!  subroutine sortgr(n,iv)
! 			compacts array iv (throws out double numbers)
!  subroutine compact(n,iv)
! 			compacts array iv (throws out double numbers)
! 
!  function nextgr(k,kgr,ngrddi,ngrade,ngri)
! 	               finds next node in grade index
!  subroutine exchgr(k,kold,knew,ngrddi,ngrade,ngri)
! 			exchanges node in grade index
!  subroutine delgr(k,kold,ngrddi,ngrade,ngri)
! 	      	       deletes grade in index (kold from index of node k)
!  subroutine insgr(k,kafter,knew,ngrddi,ngrade,ngri)
! 	               inserts grade in index 
!  subroutine prgr(k,ngrddi,ngrade,ngri)
! 	               prints grade index
! 
!  revision log :
! 
!  01.01.2003	ggu	written
!  19.05.2003	ggu	some changes in stats
!  23.03.2010	ggu	changed v6.1.1
!  25.10.2013	ggu	changed VERS_6_1_68
!  24.07.2015	ggu	changed VERS_7_1_82
!  30.07.2015	ggu	changed VERS_7_1_83
!  18.12.2018	ggu	changed VERS_7_5_52
!  21.05.2019	ggu	changed VERS_7_5_62
! 
! ***********************************************************

	subroutine maxgrd(nkn,nel,nen3v,nmax)

!  computes maximum grade (except boundary -> is one too low)

	implicit none

	integer nkn,nel,nmax
	integer nen3v(3,nel)

	integer ngrade(nkn)

	integer k,ie,ii

	do k=1,nkn
	  ngrade(k) = 0
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    ngrade(k) = ngrade(k) + 1
	  end do
	end do

	nmax = 0
	do k=1,nkn
	  if( ngrade(k) .gt. nmax ) nmax = ngrade(k)
	end do

	end

! ***********************************************************

	subroutine setgrd(nkn,nel,nen3v,ngrade)

!  computes maximum grade (except boundary -> is one too low)

	implicit none

	integer nkn,nel
	integer nen3v(3,1)
	integer ngrade(1)

	integer k,ie,ii

	do k=1,nkn
	  ngrade(k) = 0
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    ngrade(k) = ngrade(k) + 1
	  end do
	end do

	end

! ***********************************************************

	subroutine statgrd(text,nkn,nmax,ngrade,nbound)

!  computes statistics on grade

	use mod_adj_static

	implicit none

	character*(*) text
	integer nkn,nmax
	integer ngrade(nkn)
	integer nbound(nkn)

	integer nsa(0:nmax)
	integer nsi(0:nmax)
	integer nsb(0:nmax)

	integer i,k,n,nb
	integer naux,ndiff
	integer nin,nbn,nto,nis

	do i=0,nmax
	  nsa(i) = 0
	  nsi(i) = 0
	  nsb(i) = 0
	end do

	nin = 0			! total number of internal nodes
	nbn = 0			! total number of boundary nodes
	nto = 0			! total grade

	do k=1,nkn
	  n = ngrade(k)
	  nb = nbound(k)	! = 1 if on boundary
	  if( nb .eq. nbstatic ) nb = 0	!NBSTATIC
	  nto = nto + n
	  if( n .lt. 0 ) stop 'error stop statgrd: internal error'
	  nsa(n) = nsa(n) + 1
	  if( nb .ne. 0 ) then
	    nsb(n) = nsb(n) + 1
	  else
	    nsi(n) = nsi(n) + 1
	  end if

	  if( nb .ne. 0 ) then
		nbn = nbn + 1
	  else
		nin = nin + 1
	  end if
	end do

! 	ntheo = 6*nin + 4*nbn + 6*nis - 6	!theoretical grade

	naux = 6*nin + 4*nbn
	ndiff = nto-naux	! must be divisible by 6, -6 for no island
	nis = ndiff/6 + 1	! total number of islands

	write(6,'(a,a)') 'statistics on grades: ',text
	write(6,'(5x,a)') '  internal  boundary   islands    total grade'
	write(6,'(5x,3i10,i15)') nin,nbn,nis,nto
	write(6,'(3x,a)') 'grade       all  internal  boundary'
	do i=0,nmax
	  write(6,'(3x,i5,3i10)') i,nsa(i),nsi(i),nsb(i)
	end do

	if( mod(ndiff,6) .ne. 0 .and. nbn .gt. 0 ) then	!not on first call
	  write(6,*) 'error in element layout'
	  write(6,*) 'ndiff = ',ndiff
	  write(6,*) 'nis = ',nis
	  write(6,*) 'nto = ',nto
	  write(6,*) 'naux = ',naux
	  write(6,*) 'nin = ',nin
	  write(6,*) 'nbn = ',nbn
	  stop 'error stop statgrd: element layout'
	end if

	end

! ***********************************************************

	subroutine mkbound(nkn,nel,ngrddi,nen3v,ngrade,nbound,ngri)

!  marks boundary nodes and determines exact grade

	implicit none

	integer nkn,nel,ngrddi
	integer nen3v(3,nel)
	integer ngrade(nkn)
	integer nbound(nkn)
	integer ngri(2*ngrddi,nkn)

	integer k,ie,ii,n,nb
	integer i1,i2
	integer i

	!write(6,*) 'grading nodes...'

	ngrade = 0
	nbound = 0
	ngri = 0

	do ie=1,nel
	  do ii=1,3
	    i1 = mod(ii,3) + 1
	    i2 = mod(i1,3) + 1
	    k  = nen3v(ii,ie)
	    n = ngrade(k)
! 	    nodes are inserted in counter-clockwise sense
	    if( n + 2 .gt. 2*ngrddi ) goto 99
	    ngri(n+1,k) = nen3v(i1,ie)
	    ngri(n+2,k) = nen3v(i2,ie)
	    ngrade(k) = ngrade(k) + 2
	  end do
	end do

	nb = 0
	do k=1,nkn
	  n = ngrade(k)
	  call sortgr(ngrddi,n,ngri(1,k))	!compact ngri and compute n
	  if( ngrade(k)/2 .ne. n ) then		!boundary node
	    nb = nb + 1
	    nbound(k) = 1
	  end if
	  ngrade(k) = n
	  call resort_index(n,ngri(:,k))
	end do

	!write(6,*) 'grading nodes done ... boundary nodes: ',nb,nkn

	return
   99	continue
	write(6,*) 'internal error: grade too high: ',n+2,2*ngrddi
	stop 'error stop mkbound: ngrddi'
	end

! ***********************************************************

	subroutine sortgr(ngrddi,n,iv)

!  compacts array iv (throws out double numbers)

	implicit none

	integer ngrddi
	integer n
	integer iv(n)

        integer ndim
	logical bnew
	integer ib,ie,new
	integer i
	integer iaux(2*ngrddi)      !aux array

!  initialize

	ndim = 2*ngrddi
	ib = ndim/2
	ie = ib + 1

	iaux(ib) = iv(1)
	iaux(ie) = iv(2)
	iv(1) = 0
	iv(2) = 0

!  look forward

	bnew = .true.
	do while( bnew )
	  new = 0
	  do i=1,n,2		!only odd nodes
	    if( iv(i) .eq. iaux(ie) ) new = i
	  end do
	  if( new .eq. 0 ) then
	    bnew = .false.
	  else
	    ie = ie + 1
            if( ie .gt. ndim ) stop 'error stop sortgr: ndim'
	    iaux(ie) = iv(new+1)
	    iv(new) = 0
	    iv(new+1) = 0
	  end if
	end do

!  look backward

	bnew = .true.
	do while( bnew )
	  new = 0
	  do i=2,n,2		!only even nodes
	    if( iv(i) .eq. iaux(ib) ) new = i
	  end do
	  if( new .eq. 0 ) then
	    bnew = .false.
	  else
	    ib = ib - 1
            if( ib .lt. 1 ) stop 'error stop sortgr: ndim'
	    iaux(ib) = iv(new-1)
	    iv(new) = 0
	    iv(new-1) = 0
	  end if
	end do

!  adjust for internal node

	if( iaux(ie) .eq. iaux(ib) ) ie = ie - 1

!  copy to original array

	n = 0
	do i=ib,ie
	  n = n + 1
	  iv(n) = iaux(i)
	end do

	end

! ***********************************************************

	subroutine compact(n,iv)

!  compacts array iv (throws out double numbers)
! 
!  not used

	implicit none

	integer n
	integer iv(n)

	integer ivold,i,idiff,nmax

!  set to 0 all double nodes

	ivold = iv(1)
	do i=2,n
	  if( ivold .eq. iv(i) ) then
	    iv(i) = 0
	  else
	    ivold = iv(i)
	  end if
	end do

!  copy remaining nodes to front of array (leave no holes)

	idiff = 0
	nmax = 0
	do i=1,n
	  if( iv(i) .eq. 0 ) then
	    idiff = idiff + 1
	  else
	    iv(i-idiff) = iv(i)
	    nmax = nmax + 1
	  end if
	end do

!  delete useless entries

	do i=nmax+1,n
	  iv(i) = 0
	end do

	n = nmax

	end

! ***********************************************************

	subroutine resort_index(n,index)

	implicit none

	integer n
	integer index(n)

	integer i,imin
	integer iaux(n)

	imin = minloc(index,1)
	if( imin == 1 ) return

	do i=1,n
	  iaux(i) = index(imin)
	  imin = imin + 1
	  if( imin > n ) imin = 1
	end do

	index = iaux

	end

! ***********************************************************
! ***********************************************************
! ***********************************************************

	function nextgr(k,kgr,ngrddi,ngrade,ngri)

!  finds next node in grade index
! 
!  in grade index of node k find node after kgr
!  does not take into account boundary nodes

	implicit none

	integer nextgr
	integer k,kgr,ngrddi
	integer ngrade(1)
	integer ngri(2*ngrddi,1)

	integer i,n,kn,kext

	nextgr = 0

	n = ngrade(k)

	do i=1,n
	  kn = ngri(i,k)
	  if( kn .eq. kgr) goto 1
	end do

	write(6,*) k,kgr,n,ngrddi
        call nint2ext(k,kext)
        write(6,*) 'int/ext k: ',k,kext
        call nint2ext(kgr,kext)
        write(6,*) 'int/ext kgr: ',kgr,kext
        write(6,*) (ngri(i,k),i=1,n)
	stop 'error stop nextgr: no such node'

    1	continue
	i = i + 1
	if( i .gt. n ) i = 1
	nextgr = ngri(i,k)

	end

! **************************************************************

	subroutine exchgr(k,kold,knew,ngrddi,ngrade,ngri)

!  exchanges node in grade index

	implicit none

	integer k,kold,knew,ngrddi
	integer ngrade(1)
	integer ngri(2*ngrddi,1)

	integer n,i

	n = ngrade(k)

	do i=1,n
	  if( ngri(i,k) .eq. kold ) exit
	end do

	if( i > n ) then
	  write(6,*) k,kold,knew
	  stop 'error stop exchgr: no such node'
	end if

	ngri(i,k) = knew
	call resort_index(n,ngri(:n,k))

	end

! **************************************************************

	subroutine delgr(k,kold,ngrddi,ngrade,ngri)

!  deletes grade in index (kold from index of node k)

	implicit none

	integer k,kold,ngrddi
	integer ngrade(1)
	integer ngri(2*ngrddi,1)

	integer n,i,idiff

	n = ngrade(k)
	ngrade(k) = ngrade(k) - 1

	idiff = 0
	do i=1,n
	  if( idiff .eq. 0 .and. ngri(i,k) .eq. kold ) then
	    idiff = 1
	  else
	    ngri(i-idiff,k) = ngri(i,k)
	  end if
	end do
	ngri(n,k) = 0
	n = n - 1

	if( idiff .eq. 0 ) then
	  write(6,*) k,kold
	  stop 'error stop delgr: no such node'
	end if

	call resort_index(n,ngri(:n,k))

	end

! **************************************************************

	subroutine insgrb(k,kbefor,knew,ngrddi,ngrade,ngri)

!  inserts grade in index (knew in index of node k befor node kbefor)

	implicit none

	integer k,kbefor,knew,ngrddi
	integer ngrade(1)
	integer ngri(2*ngrddi,1)

	integer n,i

	n = ngrade(k)
	ngrade(k) = ngrade(k) + 1

	do i=n,1,-1			!backwards
	  if( ngri(i,k) .eq. kbefor ) exit
	  ngri(i+1,k) = ngri(i,k)
	end do

	if( i < 1 ) then
	  write(6,*) k,kbefor,knew
	  write(6,*) 'No such node : ',kbefor
	  stop 'error stop insgrb: no such node'
	end if

	ngri(i+1,k) = ngri(i,k)
	ngri(i,k) = knew
	n = n + 1
	call resort_index(n,ngri(:n,k))

	end 

! **************************************************************

	subroutine insgr(k,kafter,knew,ngrddi,ngrade,ngri)

!  inserts grade in index (knew in index of node k after node kafter)

	implicit none

	integer k,kafter,knew,ngrddi
	integer ngrade(1)
	integer ngri(2*ngrddi,1)

	integer n,i

	n = ngrade(k)
	ngrade(k) = ngrade(k) + 1

	do i=n,1,-1			!backwards
	  if( ngri(i,k) .eq. kafter ) exit
	  ngri(i+1,k) = ngri(i,k)
	end do

	if( i < 1 ) then
	  write(6,*) 'insgr: ',k,kafter,knew
	  write(6,*) 'No such node : ',kafter
	  stop 'error stop insgr'
	end if

	ngri(i+1,k) = knew
	n = n + 1
	call resort_index(n,ngri(:n,k))

	end 

! **************************************************************

	subroutine prgr(k,ngrddi,ngrade,ngri)

!  prints grade index

	implicit none

	integer k,ngrddi
	integer ngrade(1)
	integer ngri(2*ngrddi,1)

	integer n,i

	n = ngrade(k)

	write(6,*) 'grade index : ',k,n
	write(6,*) (ngri(i,k),i=1,n)

	end 

! **************************************************************

