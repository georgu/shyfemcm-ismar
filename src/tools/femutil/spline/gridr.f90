
!--------------------------------------------------------------------------
!
!    Copyright (C) 2002,2004-2005,2011-2012  Georg Umgiesser
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

! smoothing program
!
! contents :
!
! revision log :
! 
! 01.02.2002    ggu     handle periodic line
! 06.08.2002    ggu     read also depth value
! 23.09.2004    ggu     adapted for malta, bug fix
! 19.10.2005    ggu     documentation and description
! 02.12.2011    ggu     bug fix in intpdep() and reduce()
! 02.12.2011    ggu     use depth also for smoothing (change in distxy())
! 27.05.2012    ggu     renamed ndim to nsdim in smooth()
! 12.01.2026    ggu     large arrays transformed to allocatable
! 12.02.2026    ggu     new routine despike, bug fix for retriv
!
! description :
!
! This program smothes and reduces the number of points on a line
!
! The program ask for the GRD file that contains the lines
! It also asks for the smoothing parameter (sigma) and reduction
!	parameter (reduce)
! The units of these values must be in the units of the line coordinates.
! 
! Smoothing:
!
! The program uses a gaussian curve as the smoothing operator, with
! its standard deviation given by sigma. Therefore, if (x,y) coordinates
! of the lines are in meters and sigma=100, the lengthscale for the 
! smoothing is about 100 meters, which means that points closer than 100
! meters to the point will have much more effect than points further away
! than 100 m.
!
! Reduction:
!
! The program uses the parameter reduce to decide if a point should
! be kept or not. If reduce=100, than all points closer than 100 m
! to the reference point are discarded. At the end a line is left with
! a distance between the points of minimum 100 m.
!
! Adaption of resolution:
!
! If the points of the line have no depth values, then the 
! smoothing and reduction is working as described above. If, however,
! depth values do exists, then they are used to change the resolution
! of the line points. For a value of 2, twice the points in the reduction
! algorithm are retained. In other words, the value for reduce is changed
! locally to a value of 100/2=50, and only points closer as 50 are
! eliminated. Line points without depth receive an implicit value of 1.
! The value -1 indicates that this point should never be eliminated 
! nor should it be smoothed.
!
! To do:
!
! the fixed points should be handled better (stop smoothing, reducing)
! adaption of reduce on short lines, look at line globally
!
!***********************************************************************

	program gridr

! smoothing program
!
! same as gridf but smooth on length of line, not points

	use basin
	use grd
	use mod_depth

	implicit none

	integer ner
	integer nco,nk,ne,nl,nne,nnl
	logical bstop
	character*80 file

	logical bperiod
	integer nt,nll
	integer i

	real, allocatable :: xt(:)
	real, allocatable :: yt(:)
	real, allocatable :: ht(:)

	integer l
	integer nline
	integer nnode
	real sigma
	real reduce
	real rspike

	call shyfem_copyright('gridr - smoothing of lines')

!-----------------------------------
! sigma		smoothing parameter (size of gaussian kernel)
! reduce	reduction of points in line (in meters)
!-----------------------------------
	sigma = 2.0
	sigma = 1.5
	sigma = 500.
	sigma = 200.
	reduce = 4.0
	reduce = 200.
!----------------------------------- hakata
	sigma = 200.
	reduce = 200.
!----------------------------------- circle
	sigma = 0.0
	reduce = 0.1
!----------------------------------- lido
	sigma = 0.
	reduce = 100.
!----------------------------------- curonian lagoon
	sigma = 200.
	reduce = 200.
!----------------------------------- brasile
	sigma = 0.005
	reduce = 0.005
!----------------------------------- malta
	sigma = 0.1
	reduce = 0.2
	sigma = 0.1 * 245.
	reduce = 0.2 * 245.
!----------------------------------- hue
	sigma = 30.
	reduce = 300.
!-----------------------------------

	call handle_command_line(file,sigma,reduce,rspike)

!------------------------------------------------------

	nline = 0
	nnode = 0

        call grd_read(file)
        call grd_get_params(nk,ne,nl,nne,nnl)

        write(6,'(a,5i7)') 'grid parameters: ',nk,ne,nl
        if( nk .le. 0 ) stop

        call grd_to_basin
	call mod_depth_init(nkn,nel)
	call grd_get_nodal_depth(hkv)

	write(6,*) 'nodes             : ',nk
	write(6,*) 'elements          : ',ne
	write(6,*) 'lines             : ',nl
	write(6,*) 'nodes in elements : ',nne
	write(6,*) 'nodes in lines    : ',nnl

	allocate(xt(nnl),yt(nnl),ht(nnl))

	call insnod(nkn,ipv)

	open(99,file='smooth.grd',status='unknown',form='formatted')
	open(98,file='reduce.grd',status='unknown',form='formatted')
	open(97,file='despike.grd',status='unknown',form='formatted')

	do l=1,nl
	  nll = nnl
	  nline = ipplv(l)
	  call extrli(l,nl,ipplv,ialv,ipntlv,inodlv,xgv,ygv,hkv &
     &				,xt,yt,ht,nll,nt)
	  call mkperiod(xt,yt,nll,bperiod)
	  call mkstats(nl,nll,xt,yt,ht,bperiod)
	  call intpdep(nline,ht,nll,bperiod)
	  call smooth(sigma,xt,yt,ht,nll,bperiod)
	  call wrline(99,nline,nnode,nll,xt,yt,ht,nt,bperiod)
	  call reduce_points(reduce,xt,yt,ht,nll)
	  call wrline(98,nline,nnode,nll,xt,yt,ht,nt,bperiod)
	  call despike(rspike,xt,yt,ht,nll)
	  call wrline(97,nline,nnode,nll,xt,yt,ht,nt,bperiod)
	end do

	write(6,*) 'routine finished...'

	close(99)
	close(98)
	close(97)

	write(6,*) 'files smooth.grd, reduce.grd, and despike.grd written'

	end

!**********************************************************

	subroutine extrli(l,nli,iplv,ialrv,ipntlv,inodlv,xgv,ygv,hkv &
     &		,xt,yt,ht,nl,nt)

! extracts line from list

	implicit none

	integer l		!actual line
	integer nli		!total number of lines
	integer iplv(1)
	integer ialrv(1)
	integer ipntlv(0:1)
	integer inodlv(1)
	real xgv(1)
	real ygv(1)
	real hkv(1)
	integer nl		!on entry dim, on return number of nodes in line
	real xt(1)
	real yt(1)
	real ht(1)
	integer nt		!type of line

	integer nvert,i,ibase
	integer node,ier,k

	integer retriv

	if( l .gt. nli ) stop 'error stop extrli : l > nli'

	nvert = ipntlv(l) - ipntlv(l-1)
	ibase = ipntlv(l-1)

	if( nvert .gt. nl ) goto 98
	nt = ialrv(l)

	!write(6,*) 'extracting line ',l,iplv(l),nvert,nt

	do i=1,nvert
	    node = inodlv(ibase+i)
	    !ier = retriv(node,k)
	    k = node
	    if( ier .lt. 0 ) goto 99
	    if( ier .eq. 0 ) then
	      write(6,*) 'node not found: ',node
		stop
	    end if
	    xt(i) = xgv(k)
	    yt(i) = ygv(k)
	    ht(i) = hkv(k)
	end do

	nl = nvert

	return
   98	continue
	write(6,*) '*** Dimension error for nl'
	write(6,*) 'Number of vertices in line: nvert = ',nvert
	write(6,*) 'Dimension for vertices:        nl = ',nl
	write(6,*) 'Please adjust dimension of ndim'
	stop 'error stop extrli : dimension nl'
   99	continue
	stop 'error stop extrli: hash routines'
	end

!**********************************************************

	subroutine mkperiod(xt,yt,nl,bperiod)

! decides if line is periodic or not

	implicit none

	real xt(1)
	real yt(1)
	integer nl
	logical bperiod

	bperiod = .false.

	if( xt(1) .eq. xt(nl) .and. yt(1) .eq. yt(nl) ) then
	  bperiod = .true.
	  nl = nl - 1
	end if

	end

!**********************************************************

	subroutine prlixy(nli,iplv,ialrv,ipntlv,inodlv,xgv,ygv)

! print info on lines

	implicit none

	integer nli
	integer iplv(1)
	integer ialrv(1)
	integer ipntlv(0:1)
	integer inodlv(1)
	real xgv(1)
	real ygv(1)

	integer l,nvert,i,ibase
	integer node,ier,k

	integer retriv

	do l=1,nli
	  nvert = ipntlv(l) - ipntlv(l-1)
	  ibase = ipntlv(l-1)
	  write(6,*) 'line : ',l,iplv(l),ialrv(l),nvert
	  do i=1,nvert
	    node = inodlv(ibase+i)
	    !ier = retriv(node,k)
	    k = node
	    if( ier .lt. 0 ) goto 99
	    write(6,*) i,node,xgv(k),ygv(k)
	  end do
	end do

	return
   99	continue
	stop 'error stop prlixy: hash routines'
	end

!**********************************************************

	subroutine prline(nli,iplv,ialrv,ipntlv,inodlv)

! print info on lines

	implicit none

	integer nli
	integer iplv(1)
	integer ialrv(1)
	integer ipntlv(0:1)
	integer inodlv(1)

	integer l,nvert,i,ibase

	do l=1,nli
	  nvert = ipntlv(l) - ipntlv(l-1)
	  ibase = ipntlv(l-1)
	  write(6,*) 'line : ',l,iplv(l),ialrv(l),nvert
	  write(6,*) (inodlv(ibase+i),i=1,nvert)
	end do

	end

!**********************************************************

	subroutine insnod(nkn,ipv)

! inserts nodes into hash table

	implicit none

	integer nkn
	integer ipv(nkn)

	integer k,ier
	integer insert

	call init

	do k=1,nkn
	  ier = insert(ipv(k),k)
	  if( ier .lt. 0 ) goto 99
	end do

	return
   99	continue
	stop 'error stop insnod: hash routines'
	end

!********************************************************

	subroutine distxy(ndim,nl,xt,yt,ht,dxy)

! sets up distance between nodes

! dxy(i) -> distance to next node (i,i+1)
! dxy(i-1) -> distance to previous node (i,i-1)

	implicit none

	integer ndim
	integer nl
	real xt(1)
	real yt(1)
	real ht(1)
	real dxy(-ndim:2*ndim)

	integer i
	double precision xo,yo,xn,yn
	double precision dx,dy,dist
	double precision ho,hn,h

	xn = xt(nl)
	yn = yt(nl)
	hn = ht(nl)

	do i=1,nl
	  xo = xn
	  yo = yn
	  ho = hn
	  xn = xt(i)
	  yn = yt(i)
	  hn = ht(i)

	  dx = xn - xo
	  dy = yn - yo
	  dist = sqrt( dx*dx + dy*dy )
	  if( hn .gt. 0. ) dist = dist * hn

	  dxy(i-1) = dist
	  dxy(i-1+nl) = dist
	  dxy(i-1-nl) = dist
	end do

	dxy(2*nl) = dxy(0)

	end
	  
!********************************************************
!********************************************************
!********************************************************

	subroutine smooth(sigma,xt,yt,ht,nl,bperiod)

! smoothing

	implicit none

	real sigma
	real xt(1)
	real yt(1)
	real ht(1)
	integer nl
	logical bperiod

	integer nsdim
	parameter(nsdim=40000)

	integer i
	integer ngk
	real dist,dx,dy,distot
	real gk(-nsdim:nsdim)
	real, allocatable :: raux(:)
	real, allocatable :: dxy(:)
	real, allocatable :: htt(:)

	if( nl .gt. nsdim ) goto 99

        if( sigma .le. 0. ) return

	allocate(raux(-nl:2*nl))
	allocate(dxy(-nl:2*nl))
	allocate(htt(nl))

! set up dxy

	htt = 0.
	call distxy(nsdim,nl,xt,yt,htt,dxy)

	distot = 0.
	do i=1,nl-1
	  distot = distot + dxy(i)
	end do

	dx = xt(1) - xt(nl)
	dy = yt(1) - yt(nl)
	dist = sqrt ( dx**2 + dy**2 )
	write(6,*) 'distance : ',nl,dist,distot,dist/distot,bperiod

! create smoothing kernel - only use with regular point spacing

!	ngk = ndim
!	call gkernel(sigma,nsdim,ngk,gk)
!	write(6,*) 'Gaussian kernel : ',sigma,ngk

! smooth x and y

	call grsmooth(nsdim,nl,sigma,xt,raux,dxy,ht,bperiod)
	call grsmooth(nsdim,nl,sigma,yt,raux,dxy,ht,bperiod)

	return
   99	continue
	write(6,*) '*** Dimension error for nsdim'
	write(6,*) 'Number of vertices in line: nl = ',nl
	write(6,*) 'Dimension for vertices:  nsdim = ',nsdim
	write(6,*) 'Please adjust dimension of nsdim in routine smooth'
	stop 'error stop smooth: nsdim'
	end

!********************************************************

	subroutine gkernel(sigma,ndim,ngk,gk)

! gaussian smoothing - creates kernel to be used with gsmooth

	real sigma
	integer ndim
	integer ngk
	real gk(-ndim:ndim)

	integer i
	real pi,eps,a,b,x,val

        pi = 4.*atan(1.)
        eps = 1.e-7

        a = 1. / ( sqrt(2.*pi) * sigma )
        b = 1. / ( 2 * sigma * sigma )

	do i=0,ndim
          x = i
          val = a * exp( -x**2 * b )
	  gk(i) = val
	  gk(-i) = val
	  if( val .lt. eps ) goto 1
        end do
    1	continue
	ngk = i-1

!	do i=-ngk,ngk
!	  write(6,*) i,gk(i)
!	end do

	end

!********************************************************

	subroutine gsmooth(ndim,nl,rt,raux,ngk,gk)

! gaussian smoothing - only good for regular point spacing

	integer ndim
	integer nl
	real rt(1)
	real raux(-ndim:2*ndim)
	integer ngk
	real gk(-ndim:ndim)

	integer i,j
	real val

! set up circular array

	do i=1,nl
	  val = rt(i)
	  raux(i) = val
	  raux(i+nl) = val
	  raux(i-nl) = val
	end do

! convolution

	do i=1,nl
	  val = 0.
	  do j=-ngk,ngk
	    val = val + raux(i+j) * gk(j)
	  end do
	  rt(i) = val
	end do

	end

!********************************************************

	subroutine grsmooth(ndim,nl,sigma,rt,raux,dxy,ht,bperiod)

! gaussian smoothing - good for irregular point spacing

	integer ndim
	integer nl
	real sigma
	real rt(1)
	real raux(-ndim:2*ndim)
	real dxy(-ndim:2*ndim)
	real ht(1)
	logical bperiod

	integer i,j
	real val
	real rkern,tkern
	real rintv

! set up circular array

	do i=1,nl
	  val = rt(i)
	  raux(i) = val
	  raux(i+nl) = val
	  raux(i-nl) = val
	end do

! set up gaussian kernel parameters

        pi = 4.*atan(1.)
        eps = 1.e-7

        a = 1. / ( sqrt(2.*pi) * sigma )
        b = - 1. / ( 2 * sigma * sigma )

! convolution

	do i=1,nl

	  rintv = 0.5 * ( dxy(i-1) + dxy(i) )
	  rkern = a
	  weight = rkern * rintv
	  tkern = weight
	  val = raux(i) * weight
	  jmin = i - nl
	  jmax = i + nl
	  if( .not. bperiod ) then
	    jmin = max(2,jmin)
	    jmax = min(nl-1,jmax)
	  end if
!	  write(6,*) i,i,weight,rkern,tkern,val

	  rkern = 1.
	  dist = 0.
	  j = i
	  do while( rkern .gt. eps .and. j .lt. jmax )
	    j = j + 1
	    dist = dist + dxy(j-1)
	    rintv = 0.5 * ( dxy(j-1) + dxy(j) )
	    rkern = a * exp( b * dist * dist )
	    weight = rkern * rintv
	    tkern = tkern + weight
	    val = val + raux(j) * weight
!	    write(6,*) i,j,weight,rkern,tkern,val,dist
	  end do

!	  write(6,*) i,j-i,nl,dist,rkern

	  rkern = 1.
	  dist = 0.
	  j = i
	  do while( rkern .gt. eps .and. j .gt. jmin )
	    j = j - 1
	    dist = dist + dxy(j)
	    rintv = 0.5 * ( dxy(j-1) + dxy(j) )
	    rkern = a * exp( b * dist * dist )
	    weight = rkern * rintv
	    tkern = tkern + weight
	    val = val + raux(j) * weight
!	    write(6,*) i,j,weight,rkern,tkern,val,dist
	  end do

!	  write(6,*) i,i-j,nl,dist,rkern

!	  write(6,*) '*** ',i,val,rt(i),tkern

	  if( tkern .gt. 0. .and. ht(i) .ge. 0. ) then
	    rt(i) = val / tkern
	  else
	    rt(i) = raux(i)
	  end if

	end do

	end

!********************************************************

	subroutine wrline(iunit,nline,nnode,nl,xt,yt,ht,nt,bperiod)

! write line

	integer iunit
	integer nline
	integer nnode
	integer nl
	real xt(1)
	real yt(1)
	real ht(1)
	integer nt
	logical bperiod

	integer i
	integer ntot

	if( nl .le. 0 ) return

	do i=1,nl
	  write(iunit,'(i1,2i10,3e14.6)') 1,nnode+i,nt,xt(i),yt(i),ht(i)
	end do
	  
!	nline = nline + 1

	ntot = nl
	if( bperiod ) ntot = ntot + 1

	write(iunit,'(i1,3i10)') 3,nline,nt,ntot
	write(iunit,'(7i10)') (nnode+i,i=1,nl)
	if( bperiod ) then
	  write(iunit,*) nnode+1
	end if

	nnode = nnode + nl

	end

!********************************************************

	subroutine despike(ratio,xt,yt,ht,nl)

! reduces points in line

	implicit none

	real ratio
	real xt(nl)
	real yt(nl)
	real ht(nl)
	integer nl

	integer i,nnew,ispike
	real dist1,dist2,dist3
	real dx1,dy1,dx2,dy2,dx3,dy3
	real h,hmin,hmax
	real, parameter :: flag = -999.

! very simplicistic approach

	if( ratio <= 0 ) return

	hmin = minval(ht)
	hmax = maxval(ht)
	write(6,*) 'despike: hmin,hmax: ',nl,hmin,hmax

	do i=2,nl-1
	  dx1 = xt(i) - xt(i-1)
	  dy1 = yt(i) - yt(i-1)
	  dx2 = xt(i+1) - xt(i)
	  dy2 = yt(i+1) - yt(i)
	  dx3 = xt(i+1) - xt(i-1)
	  dy3 = yt(i+1) - yt(i-1)
	  dist1 = sqrt( dx1*dx1 + dy1*dy1 )
	  dist2 = sqrt( dx2*dx2 + dy2*dy2 )
	  dist3 = sqrt( dx3*dx3 + dy3*dy3 )
	  if( (dist1+dist2)/dist3 > ratio ) then
	    ht(i) = flag
	    write(6,*) i,dist1+dist2,dist3,(dist1+dist2)/dist3
	  end if
	end do

	ispike = count( ht == flag )
	write(6,*) ispike,' spikes found'

	nnew = 0
	do i=1,nl
	  if( ht(i) /= flag ) then
	    nnew = nnew + 1
	    xt(nnew) = xt(i)
	    yt(nnew) = yt(i)
	    ht(nnew) = ht(i)
	  end if
	end do

!	if( nnew .lt. 3 ) stop 'error stop reduce: line too short...'
	if( nnew .lt. 3 ) then
	  write(6,*) 'line too short -> eliminated'
	  nnew = 0
	end if

	write(6,*) 'despike: from ',nl,'  to ',nnew,'  nodes'
	nl = nnew

	end

!********************************************************

	subroutine reduce_points(reduce,xt,yt,ht,nl)

! reduces points in line

	implicit none

	real reduce
	real xt(nl)
	real yt(nl)
	real ht(nl)
	integer nl

	integer i,nnew
	real rr,rtot
	real dist,dx,dy
	real h,hmin,hmax

! very simplicistic approach

	if( reduce <= 0 ) return

	hmin = minval(ht)
	hmax = maxval(ht)
	write(6,*) 'hmin,hmax: ',hmin,hmax

	rr = reduce
	rtot = 1.

	if( nl .lt. 3 * rr ) then
	  rr = nl / 3.
	end if

	nnew = 0
	rtot = 0.
	if( ht(1) .lt. 0. .or. ht(1) .ge. 0. ) then	!bug fix 23.09.2004
	    nnew = nnew + 1
	    xt(nnew) = xt(1)
	    yt(nnew) = yt(1)
	    ht(nnew) = ht(1)
	end if
	do i=2,nl
	  dx = xt(i) - xt(i-1)
	  dy = yt(i) - yt(i-1)
	  h = ht(i)
	  dist = sqrt( dx*dx + dy*dy )
	  !write(6,*) dist,h,dist/h,rtot,reduce
	  !if( h .gt. 0. ) dist = dist * h
	  rtot = rtot + dist
	  if( rtot .gt. reduce .or. h .lt. 0. ) then
	    nnew = nnew + 1
	    xt(nnew) = xt(i)
	    yt(nnew) = yt(i)
	    ht(nnew) = ht(i)
	    rtot = 0.
	  end if
	end do

!	if( nnew .lt. 3 ) stop 'error stop reduce: line too short...'
	if( nnew .lt. 3 ) then
	  write(6,*) 'line too short -> eliminated'
	  nnew = 0
	end if

	write(6,*) 'reduce: from ',nl,'  to ',nnew,'  nodes'
	nl = nnew

	end

!********************************************************

	  subroutine intpdep(nline,ht,nl,bperiod)

! interpolates depth values for points in line
!
! negative values are not interpolated and are left alone

	implicit none

      integer nline
	real ht(1)
	integer nl
	logical bperiod

	integer i
	integer ifirst,ilast,inext,nonzero
	integer nval
	real value,vfirst,vlast
	real aux,hflag
	logical bdebug

	bdebug = .true.
	bdebug = nline .eq. 9
	bdebug = nline .eq. 6
	bdebug = .false.

        hflag = -990.

!------------------------------------------------------------
! look for values greater than 0
!------------------------------------------------------------

	nonzero = 0
	ifirst = 0
	ilast = 0

	do i=1,nl
	  if( ht(i) .gt. 0. ) then
	    nonzero = nonzero + 1
	    if( nonzero .eq. 1 ) ifirst = i
	    ilast = i
	  end if
	end do

      if( bdebug ) write(6,*) 'nz: ',nl,nonzero,ifirst,ilast

!------------------------------------------------------------
! if no value or 1 value found -> set everything constant
!------------------------------------------------------------

	if( nonzero .le. 1 ) then
	  if( nonzero .eq. 0 ) then
	    value = 1.
	  else
	    value = ht(ifirst)
	  end if
	  do i=1,nl
	    if( ht(i) .ge. 0. .or. ht(i) .lt. hflag ) ht(i) = value
	  end do
	  return
	end if

!------------------------------------------------------------
! interpolate extreme values
!------------------------------------------------------------

	vfirst = ht(ifirst)
	vlast = ht(ilast)
	nval = nl - ilast + ifirst
	aux = (vfirst-vlast)/nval
	inext = 0

	if( bdebug ) write(6,*) bperiod
	if( bdebug ) write(6,*) ifirst,ilast,nval,vfirst,vlast

	i = 1
	do while( ht(i) .lt. 0 .and. ht(i) .gt. hflag )
	  i = i + 1
	end do
	if( bperiod ) then
	  inext = nl - ilast + i
	  value = vlast + inext * aux
	else
	  value = ht(ifirst)
	end if
	if( bdebug ) write(6,*) 'first: ',i,inext,value,ht(i)
	ht(i) = value
	ifirst = i

	i = nl
	do while( ht(i) .lt. 0 .and. ht(i) .gt. hflag )
	  i = i - 1
	end do
	if( bperiod ) then
	  inext = i - ilast
	  value = vlast + inext * aux
	else
	  value = ht(ilast)
	end if
	if( bdebug ) write(6,*) 'last: ',i,inext,value,ht(i)
	ht(i) = value
	ilast = i

!------------------------------------------------------------
! interpolate central values
!------------------------------------------------------------

	do while( ifirst .lt. ilast )
	  do i=ifirst+1,nl
	    if( ht(i) .gt. 0. ) goto 1
	  end do
	  stop 'error stop intpdep: impossible branch (1)'
    1	  continue
	  inext = i
	  nval = inext - ifirst
	  if( bdebug ) write(6,*) ifirst,inext,ilast,ht(ifirst),ht(inext)
	  do i=ifirst+1,inext-1
	    value = ht(ifirst) + (i-ifirst) * (ht(inext)-ht(ifirst))/nval
	    if( ht(i) .ge. 0. .or. ht(i) .lt. hflag ) ht(i) = value
	    if( bdebug ) write(6,*) i,value,ht(i)
	  end do
	  ifirst = inext
	end do

!------------------------------------------------------------
! do some sanity check
!------------------------------------------------------------

	do i=1,nl
	  if( ht(i) .eq. 0 ) then
	    write(6,*) i,nl,ht(i)
	    stop 'error stop intpdep: values not interpolated'
	  end if
	  if( bdebug ) write(6,*) i,ht(i)
	end do

!------------------------------------------------------------
! end of routine
!------------------------------------------------------------

	end

!********************************************************

	subroutine mkstats(nl,nll,xt,yt,ht,bperiod)

	implicit none

	integer nl
	integer nll
	real xt(nll)
	real yt(nll)
	real ht(nll)
	logical bperiod

	integer i,ndim
	real dmin,dmax,dmed,d
	real, allocatable :: dxy(:)
	real, allocatable :: htt(:)

	ndim = nll
	allocate(dxy(-ndim:2*ndim))
	allocate(htt(ndim))

	htt = 0
	call distxy(ndim,nll,xt,yt,htt,dxy)

	dmin = 1.e+10
	dmax = 0.
	dmed = 0.

	do i=1,nll
	  d = dxy(i)
	!write(6,*) xt(i),yt(i),dxy(i)
	  dmin = min(dmin,d)
	  dmax = max(dmax,d)
	  dmed = dmed + d
	end do
	dmed = dmed / nll

	write(6,*) nl,nll,dmin,dmax,dmed

	end

!********************************************************

	subroutine handle_command_line(file,sigma,reduce,rspike)

	use clo

	implicit none

	character*(*) file
	real sigma
	real reduce
	real rspike

	integer nfile

        call clo_init('smooth','grd-file','1.0')

        call clo_add_info('smoothes line and reduces points')
        call clo_add_option('sigma sigma',0. &
     &                    ,'standard deviation for smoothing is sigma')
        call clo_add_option('reduce dmin',0. &
     &                    ,'elimination of points with distance < dmin')
        call clo_add_option('despike rspike',0. &
     &                    ,'elimination of points with spike ratio > rspike')

        call clo_parse_options(1)       !expecting 1 file

        call clo_get_option('sigma',sigma)
        call clo_get_option('reduce',reduce)
        call clo_get_option('despike',rspike)

        nfile = clo_number_of_files()
        if( nfile > 0 ) call clo_get_file(1,file)

	write(6,*) 'file name: ',trim(file)
	write(6,*) 'sigma: ',sigma
	write(6,*) 'reduce: ',reduce
	write(6,*) 'rspike: ',rspike

	end

!********************************************************

