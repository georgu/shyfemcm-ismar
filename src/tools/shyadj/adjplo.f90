
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

!  description :
! 
!  plotting routines
! 
!  contents :
! 
!  subroutine plobas(nkn,nel,nen3v,ngrade,xgv,ygv)
! 			plots basin
!  subroutine pltgrd(n,k,k1,k2,xgv,ygv)
! 			plots grade in color
!  subroutine pltsgrd(igr,nkn,nel,nen3v,ngrade,xgv,ygv)
! 			plots only special grade
!  subroutine write_grid(file,nkn,nel,xgv,ygv,nen3v)
! 			writes quick and dirty results to file
! 
!  revision log :
! 
!  01.01.2003	ggu	written
!  19.05.2003	ggu	plot some more info on debug plot
!  23.03.2010	ggu	changed v6.1.1
!  30.03.2012	ggu	changed VERS_6_1_51
!  19.01.2015	ggu	changed VERS_7_1_2
!  10.07.2015	ggu	changed VERS_7_1_50
!  24.07.2015	ggu	changed VERS_7_1_82
!  30.07.2015	ggu	changed VERS_7_1_83
!  12.10.2015	ggu	changed VERS_7_3_3
!  18.12.2018	ggu	changed VERS_7_5_52
!  21.05.2019	ggu	changed VERS_7_5_62
!  06.03.2026   ggu     completely restructured
! 
! ***********************************************************

	subroutine plobas

!  plots basin

	use mod_adj_grade
	use basin

	implicit none

	integer ie,ii,k,n,i1,i2
	real xmin,xmax,ymin,ymax

	call mima(xgv,nkn,xmin,xmax)
	call mima(ygv,nkn,ymin,ymax)

	if( bverbose ) write(6,*) 'start plotting basin...'

	call qstart

	call qworld(xmin,ymin,xmax,ymax)
	call qrcfy

	do ie=1,nel
	  k = nen3v(3,ie)
	  call qmove(xgv(k),ygv(k))
	  do ii=1,3
	    k = nen3v(ii,ie)
! 	    if( k .gt. 0 ) call qplot(xgv(k),ygv(k))	!HACK
	    call qplot(xgv(k),ygv(k))
	  end do
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    n = ngrade(k)
	    if( n .ne. 6 .and. n .gt. 0 ) then
		i1 = mod(ii,3) + 1
		i2 = mod(i1,3) + 1
		call pltgrd(n,k,nen3v(i1,ie),nen3v(i2,ie),xgv,ygv)
	    end if
	  end do
	end do

	call qend

	if( bverbose ) write(6,*) 'finished plotting basin...'

	end

! ***********************************************************

	subroutine pltgrd(n,k,k1,k2,xgv,ygv)

!  plots grade in color

	implicit none

	integer n,k,k1,k2
	real xgv(1),ygv(1)
	real x(3),y(3)
	real col

	if( n .lt. 6 ) col = 0.35
	if( n .lt. 5 ) col = 0.15
	if( n .gt. 6 ) col = 0.65
	if( n .gt. 7 ) col = 0.85
	if( n .gt. 8 ) col = 0.95

	x(1) = xgv(k)
	y(1) = ygv(k)
	x(2) = xgv(k) + 0.3 * (xgv(k1)-xgv(k))
	y(2) = ygv(k) + 0.3 * (ygv(k1)-ygv(k))
	x(3) = xgv(k) + 0.3 * (xgv(k2)-xgv(k))
	y(3) = ygv(k) + 0.3 * (ygv(k2)-ygv(k))

	call qhue(col)
	call qafill(3,x,y)

	end

! ***********************************************************

	subroutine pltsgrd(igr,nkn,nel,nen3v,ngrade,xgv,ygv)

!  plots only special grade

	implicit none

	integer igr
	integer nkn,nel
	integer nen3v(3,1)
	integer ngrade(1)
	real xgv(1),ygv(1)

	integer k,n

	do k=1,nkn
	  n = ngrade(k)
	  if( n .ne. igr ) then
	    ngrade(k) = ngrade(k) - 1000
	  end if
	end do
	call plobas
	do k=1,nkn
	  n = ngrade(k)
	  if( n .ne. igr ) then
	    ngrade(k) = ngrade(k) + 1000
	  end if
	end do

	end

! ***********************************************************

	subroutine wr0grd

!  writes quick and dirty results from basin to file

	implicit none

	call write_grid('new0.grd')

	end

! ***********************************************************

	subroutine write_grid(file)

!  writes quick and dirty results to file

	implicit none

	character*(*) file

	call basin_to_grd
	call grd_write(file)

	end

! ***********************************************************
! ***********************************************************
! ***********************************************************

	subroutine plosel2(ie1,ie2)

!  plots element

	use mod_adj_grade
	use basin

	implicit none

	integer ie1,ie2

	character*11 line
	integer ie,ii,k,n
	integer iev(2),iee,i,kk,id
	real xmin,xmax,ymin,ymax
	real dx,dy
	real x,y,xc,yc

	integer ialfa
	iev(1) = ie1
	iev(2) = ie2

	xmin = 1.e+35
	xmax = -xmin
	ymin = xmin
	ymax = -ymin
	x = 0.
	y = 0.

	do iee=1,2
	  ie = iev(iee)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    n = ngrade(k)
	    do i=1,n
	      kk = ngri(i,k)
	      x = xgv(kk)
	      y = ygv(kk)
	      if( x .gt. xmax ) xmax = x
	      if( x .lt. xmin ) xmin = x
	      if( y .gt. ymax ) ymax = y
	      if( y .lt. ymin ) ymin = y
	    end do
	  end do
	end do

	dx = xmax - xmin
	dy = ymax - ymin
	xmin = xmin - 0.2*dx
	xmax = xmax + 0.2*dx
	ymin = ymin - 0.2*dy
	ymax = ymax + 0.2*dy

	call qstart

	call qworld(xmin,ymin,xmax,ymax)
	call qrcfy

	do iee=1,2
	  ie = iev(iee)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    xc = xgv(k)
	    yc = ygv(k)
	    id = ialfa(float(k),line,-1,-1)
	    call qtext(xc,yc,line)
	    n = ngrade(k)
	    do i=1,n
	      kk = ngri(i,k)
	      x = xgv(kk)
	      y = ygv(kk)
	      call qline(xc,yc,x,y)
	      id = ialfa(float(kk),line,-1,-1)
	      call qtext(x,y,line)
	    end do
	  end do
	end do

	call qend

	end

! ***********************************************************

	subroutine plosno(k)

!  plots node and neighborhood

	use mod_adj_grade
	use basin

	implicit none

	integer k

	integer n,ip,i,kk
	real xmin,xmax,ymin,ymax

	call mnmx0(k,xmin,xmax,ymin,ymax)

	call qstart

	call qworld(xmin,ymin,xmax,ymax)
	call qrcfy

        n = ngrade(k)

	call plonn(k)

	do i=1,n
	   kk = ngri(i,k)
	   call plosno0(kk)
	   call plonn(kk)
	end do

	call qend

	end

! ***********************************************************

	subroutine plosno0(k)

!  plots one node

	use mod_adj_grade
	use basin

	implicit none

	integer k

	integer n,ip,i,kk
	real x,y,xc,yc

	xc = xgv(k)
	yc = ygv(k)
        n = ngrade(k)

	do i=1,n
	   kk = ngri(i,k)
	   x = xgv(kk)
	   y = ygv(kk)
	   call qline(xc,yc,x,y)
	end do

	end

! ***********************************************************

	subroutine plonn(k)

!  plots node number

	use mod_adj_grade
	use basin

	implicit none

	integer k

	character*11 line
	logical isbound
	integer id,idnew,idtot
	integer ng
	real x,y
	integer ialfa

	line = ' '

	isbound = nbound(k) .eq. 1
        ng = ngrade(k)
	id = ialfa(float(k),line,-1,-1)
	idtot = id

	if( isbound ) then
	  line(id+1:id+4) = ' (B)'
	  idtot = id + 4
	else if( ng .ne. 6 ) then
	  line(id+1:id+2) = ' ('
	  idnew = ialfa(float(ng),line(id+3:),-1,-1)
	  idtot = id + 3 + idnew
	  line(idtot:idtot) = ')'
	end if

	x = xgv(k)
	y = ygv(k)
	call qtext(x,y,line(1:idtot))

	end

! ******************************************************

	subroutine mnmx0(k,xmin,xmax,ymin,ymax)

!  computes min/max of (x,y) of node k and neighbors

	use mod_adj_grade

	implicit none

	integer k
	real xmin,xmax,ymin,ymax

	integer n,ip,i,kk
	real dx,dy

	xmin = 1.e+35
	xmax = -xmin
	ymin = xmin
	ymax = -ymin

        n = ngrade(k)

	call mnmx(k,xmin,xmax,ymin,ymax)

	do i=1,n
	   kk = ngri(i,k)
	   call mnmx(kk,xmin,xmax,ymin,ymax)
	end do

	dx = xmax - xmin
	dy = ymax - ymin
	xmin = xmin - 0.2*dx
	xmax = xmax + 0.2*dx
	ymin = ymin - 0.2*dy
	ymax = ymax + 0.2*dy

	end

! ******************************************************

	subroutine mnmx(k,xmin,xmax,ymin,ymax)

!  computes min/max of (x,y) of node k
!  xmin... must be already initialized

	use mod_adj_grade
	use basin

	implicit none

	integer k
	real xmin,xmax,ymin,ymax

	integer n,ip,i,kk
	real x,y

        n = ngrade(k)

	do i=1,n
	   kk = ngri(i,k)
	   x = xgv(kk)
	   y = ygv(kk)
	   if( x .gt. xmax ) xmax = x
	   if( x .lt. xmin ) xmin = x
	   if( y .gt. ymax ) ymax = y
	   if( y .lt. ymin ) ymin = y
	end do

	end

! ******************************************************

