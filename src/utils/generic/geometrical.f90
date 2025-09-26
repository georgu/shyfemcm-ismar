
!--------------------------------------------------------------------------
!
!    Copyright (C) 1999,2001,2010-2011,2014-2015,2014-2015  Georg Umgiesser
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

! geometrical routines
!
! contents :
!
! function areat(x1,y1,x2,y2,x3,y3)
! function areapoly(n,x,y)
!
! subroutine centert(x1,y1,x2,y2,x3,y3,xc,yc)
! subroutine centerpoly(n,x,y,xc,yc)
!
! function left(x1,y1,x2,y2,x3,y3)
! function lefton(x1,y1,x2,y2,x3,y3)
! function collinear(x1,y1,x2,y2,x3,y3)
!
! function between(x1,y1,x2,y2,x3,y3)
! function betweens(x1,y1,x2,y2,x3,y3)
! function angle(x1,y1,x2,y2,x3,y3)
!
! subroutine dist_point_to_line(px,py,ax,ay,bx,by,qx,qy,d)
!
! function segsegint(x1,y1,x2,y2,x3,y3,x4,y4,xi,yi)
! function parallelint(x1,y1,x2,y2,x3,y3,x4,y4,xi,yi)
!
! function isconvex(n,x,y)
! function inconvex(n,x,y,x0,y0)
! function polyinpoly(n1,x1,y1,n2,x2,y2)
!
! function inpoly(n,x,y,x0,y0)		!shell for inpoly0 or inpoly1
! function inpoly1(n,x,y,x0,y0)		!fast algorithm
! function inpoly0(n,x,y,x0,y0)		!using winding number
! function winding(n,x,y,x0,y0)
!
! subroutine c2p_ocean(u,v,s,d)
!
! subroutine c2p(u,v,s,d)		converts cartesian to polar coordinates
! subroutine p2c(s,d,u,v)		converts polar to cartesian coordinates
! subroutine polar2xy(n,speed,dir,uv,vv)	as p2c but for arrays
! subroutine xy2polar(n,uv,vv,speed,dir)	as c2p but for arrays
!
! revision log :
!
! 18.03.1999	ggu	written from scratch (based on Joseph O-Rourke, 1998)
! 28.03.1999	ggu	bug fix : area2 -> areat (area was 2x)
! 28.03.1999	ggu	new routine polyinpoly
! 29.03.1999	ggu	segsegint is working with double precision internally
! 19.11.1999	ggu	new routines inpoly and winding
! 06.02.2001	ggu	new routine angle
! 22.03.2010	ggu	new routine inpoly() has been implemented
! 05.07.2010	ggu	new routine dist_point_to_line(), bug fix in inpoly1()
! 22.07.2010	ggu	changed VERS_6_1_9
! 14.07.2011	ggu	changed VERS_6_1_27
! 21.02.2014	ggu	new routine centert and centerpoly
! 07.03.2014	ggu	changed VERS_6_1_72
! 14.09.2015	ggu	changed VERS_7_2_2
! 07.10.2017	ggu	new routine c2p_ocean()
! 31.08.2018	ggu	changed VERS_7_5_49
! 25.10.2018	ggu	changed VERS_7_5_51
! 18.12.2018	ggu	changed VERS_7_5_52
! 14.02.2019	ggu	changed VERS_7_5_56
! 16.02.2019	ggu	changed VERS_7_5_60
! 23.09.2025	ggu	more routines to convert cartesian to polar and vicev.
!
!**************************************************************
!**************************************************************
!**************************************************************

	function areat(x1,y1,x2,y2,x3,y3)

! computes area of triangle

	implicit none

	real areat
	real x1,y1,x2,y2,x3,y3

	areat = 0.5 * ( (x2-x1) * (y3-y1) - (x3-x1) * (y2-y1) )

	end

!**************************************************************

	function areapoly(n,x,y)

! computes area of polygon - first point must be different from last point

	implicit none

	real areapoly
	integer n
	real x(n),y(n)

	integer i
	real sum
	real x1,y1,xo,yo,xn,yn

	real areat

	sum = 0.

	x1 = x(1)
	y1 = y(1)
	xn = x(2)
	yn = y(2)

	do i=3,n

	  xo = xn
	  yo = yn
	  xn = x(i)
	  yn = y(i)

	  sum = sum + areat(x1,y1,xo,yo,xn,yn)

	end do

	areapoly = sum

	end

!**************************************************************
!**************************************************************
!**************************************************************

	subroutine centert(x1,y1,x2,y2,x3,y3,xc,yc)

! computes center of triangle

	implicit none

	real x1,y1,x2,y2,x3,y3
	real xc,yc

	xc = (x1+x2+x3) / 3
	yc = (y1+y2+y3) / 3

	end

!**************************************************************

	subroutine centerpoly(n,x,y,xc,yc)

! computes center of polygon (center of gravity)

	implicit none

	integer n
	real x(n),y(n)
	real xc,yc

	integer i
	real x1,y1,xo,yo,xn,yn
	real area
	double precision areasum,xd,yd

        real areat

	xc = 0.
	yc = 0.

	if( n .le. 0 ) return

        areasum = 0.
	xd = 0.
	yd = 0.

        x1 = x(1)
        y1 = y(1)
        xn = x(2)
        yn = y(2)

        do i=3,n

          xo = xn
          yo = yn
          xn = x(i)
          yn = y(i)

          area = areat(x1,y1,xo,yo,xn,yn)

	  areasum = areasum + area
	  xd = xd + area*(x1+xo+xn)
	  yd = yd + area*(y1+yo+yn)

        end do

	xc = xd / (3.*areasum)
	yc = yd / (3.*areasum)

	end

!**************************************************************
!**************************************************************
!**************************************************************

	function left(x1,y1,x2,y2,x3,y3)

! left turn ?

	implicit none

	logical left
	real x1,y1,x2,y2,x3,y3

	real areat

	left = areat(x1,y1,x2,y2,x3,y3) .gt. 0.

	end

!**************************************************************

	function lefton(x1,y1,x2,y2,x3,y3)

! left turn or straight ?

	implicit none

	logical lefton
	real x1,y1,x2,y2,x3,y3

	real areat

	lefton = areat(x1,y1,x2,y2,x3,y3) .ge. 0.

	end

!**************************************************************

	function collinear(x1,y1,x2,y2,x3,y3)

! left turn ?

	implicit none

	logical collinear
	real x1,y1,x2,y2,x3,y3

	real areat

	collinear = areat(x1,y1,x2,y2,x3,y3) .eq. 0.

	end

!**************************************************************
!**************************************************************
!**************************************************************

	function between(x1,y1,x2,y2,x3,y3)

! is 3 between 1 and 2 ?

	implicit none

	logical between
	real x1,y1,x2,y2,x3,y3

	logical b1,b2
	logical collinear

	if( .not. collinear(x1,y1,x2,y2,x3,y3) ) then	!not on one line
	  between = .false.
	else if( x1 .ne. x2 ) then			!check on x coordinate
	  b1 = x1 .le. x3 .and. x3 .le. x2
	  b2 = x1 .ge. x3 .and. x3 .ge. x2
	  between = b1 .or. b2
	else						!check on y coordinate
	  b1 = y1 .le. y3 .and. y3 .le. y2
	  b2 = y1 .ge. y3 .and. y3 .ge. y2
	  between = b1 .or. b2
	end if

	end

!**************************************************************

	function betweens(x1,y1,x2,y2,x3,y3)

! is 3 between 1 and 2 ? (we already know that points are collinear)

	implicit none

	logical betweens
	real x1,y1,x2,y2,x3,y3

	logical b1,b2

	if( x1 .ne. x2 ) then				!check on x coordinate
	  b1 = x1 .le. x3 .and. x3 .le. x2
	  b2 = x1 .ge. x3 .and. x3 .ge. x2
	  betweens = b1 .or. b2
	else						!check on y coordinate
	  b1 = y1 .le. y3 .and. y3 .le. y2
	  b2 = y1 .ge. y3 .and. y3 .ge. y2
	  betweens = b1 .or. b2
	end if

	end

!**************************************************************

	function angle(x1,y1,x2,y2,x3,y3)

! angle between 1-2-3 in degrees 
!
! a < 180 => right turn
! a = 180 => straight
! a > 180 => left turn

	implicit none

	real angle
	real x1,y1,x2,y2,x3,y3

	double precision null,full,pi,rad
	parameter( null = 0. , full = 360. )
	parameter( pi = 3.14159265 , rad = 180. / pi )

	double precision norm,aux,ang
	double precision dx1,dy1,dx2,dy2

        dx1=x1-x2
        dy1=y1-y2
        dx2=x3-x2
        dy2=y3-y2
 
        norm = sqrt( (dx1*dx1+dy1*dy1) * (dx2*dx2+dy2*dy2) )
        aux = ( dx1*dx2 + dy1*dy2 ) / norm
 
        ang = acos( aux ) * rad
        if( dx1*dy2 - dy1*dx2 .lt. null ) ang = full - ang

	angle = ang

	end

!**************************************************************
!**************************************************************
!**************************************************************

        subroutine dist_point_to_line(px,py,ax,ay,bx,by,qx,qy,d)

! computes distance from point P to line given by A - B
! returns point Q on line closest to P and its distance d

        implicit none

        real px,py              ! single point P
        real ax,ay              ! point A defining line
        real bx,by              ! point B defining line
        real qx,qy              ! closest point Q on line to P (return)
        real d                  ! distance between Q and P (return)

        real rx,ry,sx,sy
        real lambda

        rx = bx - ax
        ry = by - ay
        sx = px - ax
        sy = py - ay

        lambda = (rx*sx + ry*sy) / (rx*rx + ry*ry)

        qx = ax + lambda * rx
        qy = ay + lambda * ry

        d = sqrt( (qx-px)**2 + (qy-py)**2 )

        end

!**************************************************************
!**************************************************************
!**************************************************************

	function segsegint(x1r,y1r,x2r,y2r,x3r,y3r,x4r,y4r,xi,yi)

! intersects two lines ( 1-2 and 3-4 ) and returns intersection in (xi,yi)
!
! return code :
!			0	no intersection
!			1	proper intersection
!			2	segments collinearly overlap
!			3	endpoint of one segment is on other
!
! test for .ne. 0 to find out on any intersection

	implicit none

	integer segsegint
	real x1r,y1r,x2r,y2r,x3r,y3r,x4r,y4r,xi,yi

	real*8 x1,y1,x2,y2,x3,y3,x4,y4
	real*8 denom,nums,numt
	real*8 s,t
	logical bsi,bti,bso,bto
	integer parallelint

	logical bdebug
	bdebug = .false.

! convert to double precision

	x1 = x1r
	x2 = x2r
	x3 = x3r
	x4 = x4r
	y1 = y1r
	y2 = y2r
	y3 = y3r
	y4 = y4r

! compute denominator

	denom = (x1-x2)*(y4-y3) + (x4-x3)*(y2-y1)

! if denominator segments are parallel -> handle seperately

	if( denom .eq. 0. ) then
	  segsegint = parallelint(x1r,y1r,x2r,y2r,x3r,y3r,x4r,y4r,xi,yi)
	  return
	end if

! compute intersection

	nums = x1*(y4-y3) + x3*(y1-y4) + x4*(y3-y1)
	s = nums / denom

	if( nums .eq. 0.0 .or. nums .eq. denom ) segsegint = 3
	if( s .eq. 0.0 .or. s .eq. 1.0 ) segsegint = 3

	numt = - ( x1*(y3-y2) + x2*(y1-y3) + x3*(y2-y1) )
	t = numt / denom

	if( numt .eq. 0.0 .or. numt .eq. denom ) segsegint = 3
	if( t .eq. 0.0 .or. t .eq. 1.0 ) segsegint = 3

! see if segments intersect

	bsi = 0.0 .lt. s .and. s .lt. 1.0
	bti = 0.0 .lt. t .and. t .lt. 1.0

	bso = 0.0 .gt. s .or. s .gt. 1.0
	bto = 0.0 .gt. t .or. t .gt. 1.0

	if( bsi .and. bti ) then
	  segsegint = 1
	else if( bso .or. bto ) then
	  segsegint = 0
	else
	  if( segsegint .ne. 3 ) then
		stop 'error stop segsegint: internal error'
	  end if
	end if

! compute point of intersection

	xi = x1 + s * ( x2 - x1 )
	yi = y1 + s * ( y2 - y1 )

	if( bdebug ) then
	  write(6,*) '---------- segsegint -------------'
	  write(6,*) x1,y1,x2,y2
	  write(6,*) x3,y3,x4,y4
	  write(6,*) denom,nums,numt,s,t
	  write(6,*) bsi,bti,bso,bto
	  write(6,*) segsegint,xi,yi
	  write(6,*) '----------------------------------'
	end if

	end

!**************************************************************

	function parallelint(x1,y1,x2,y2,x3,y3,x4,y4,xi,yi)

! intersects two parallel lines ( 1-2 and 3-4 ) and returns intersection
!
! return code :
!			0	no intersection
!			2	segments collinearly overlap

	implicit none

	integer parallelint
	real x1,y1,x2,y2,x3,y3,x4,y4,xi,yi

	logical collinear,betweens

	if( .not. collinear(x1,y1,x2,y2,x3,y3) ) then
	  parallelint = 0
	else if( betweens(x1,y1,x2,y2,x3,y3) ) then
	  xi = x3
	  yi = y3
	  parallelint = 2
	else if( betweens(x1,y1,x2,y2,x4,y4) ) then
	  xi = x4
	  yi = y4
	  parallelint = 2
	else if( betweens(x3,y3,x4,y4,x1,y1) ) then
	  xi = x1
	  yi = y1
	  parallelint = 2
	else if( betweens(x3,y3,x4,y4,x2,y2) ) then
	  xi = x2
	  yi = y2
	  parallelint = 2
	else
	  parallelint = 0
	end if

	end

!**************************************************************
!**************************************************************
!**************************************************************

	function isconvex(n,x,y)

! is polygon convex ?

	implicit none

	logical isconvex
	integer n
	real x(n),y(n)

	integer i
	real xl,yl,xm,ym,xn,yn
	logical lefton

	isconvex = .false.

	xm = x(n-1)
	ym = y(n-1)
	xn = x(n)
	yn = y(n)

	do i=1,n
	  xl = xm
	  yl = ym
	  xm = xn
	  ym = yn
	  xn = x(i)
	  yn = y(i)

	  if( .not. lefton(xl,yl,xm,ym,xn,yn) ) then
!	    write(6,*) 'isconvex: ',i-1,xm,ym
	    return
	  end if
	end do

	isconvex = .true.

	end

!**************************************************************

	function inconvex(n,x,y,x0,y0)

! checks if point (x0,y0) is in convex polygon (border is inside)

	implicit none

        logical inconvex !true if (x0,y0) is inside closed line (return)
	integer n        !total number of points in line
	real x(n),y(n)   !coordinates of line
	real x0,y0       !coordinates of point to check

	integer i
	real xm,ym,xn,yn
	logical lefton

	inconvex = .false.

	xn = x(n)
	yn = y(n)

	do i=1,n
	  xm = xn
	  ym = yn
	  xn = x(i)
	  yn = y(i)

	  if( .not. lefton(xm,ym,xn,yn,x0,y0) ) return
	end do

	inconvex = .true.

	end

!**************************************************************

	function polyinpoly(n1,x1,y1,n2,x2,y2)

! polygon 1 in convex polygon 2 ? (border is inside)

	implicit none

	logical polyinpoly
	integer n1,n2
	real x1(n1),y1(n1)
	real x2(n2),y2(n2)

	integer i
	logical inconvex

	polyinpoly = .false.

	do i=1,n1
	  if( .not. inconvex(n2,x2,y2,x1(i),y1(i)) ) return
	end do

	polyinpoly = .true.

	end

!**************************************************************
!**************************************************************
!**************************************************************

	function inpoly(n,x,y,x0,y0)

! shell for in-polygon check

	implicit none

        logical inpoly   !true if (x0,y0) is inside closed line (return)
	integer n        !total number of points in line
	real x(n),y(n)   !coordinates of line
	real x0,y0       !coordinates of point to check

	logical inpoly0	!old routine -> use classical winding number
	logical inpoly1	!new routine -> uses new algorithm from Dan Sunday

	if( n .gt. 0 ) then
	  inpoly = inpoly1(n,x,y,x0,y0)		!use new by default
	else if( n .lt. 0 ) then
	  inpoly = inpoly0(-n,x,y,x0,y0)
	else
	  inpoly = .false.
	end if

	end

!**************************************************************

	function inpoly1(n,x,y,x0,y0)

! checks if point (x0,y0) is in general polygon (border ?)
!
! fast algorithm from Dan Sunday
!
! http://softsurfer.com/Archive/algorithm_0103/algorithm_0103.htm

	implicit none

        logical inpoly1  !true if (x0,y0) is inside closed line (return)
	integer n        !total number of points in line
	real x(n),y(n)   !coordinates of line
	real x0,y0       !coordinates of point to check

	integer i,iw
	real xo,yo,xn,yn
	real area

	real areat

	iw = 0
	xo=x(n)
	yo=y(n)

	do i=1,n
	  xn=x(i)
	  yn=y(i)
	  area = areat(xo,yo,xn,yn,x0,y0)

	  if( yo .le. y0 ) then
	    if( yn .gt. y0 ) then
	      if( area .gt. 0. ) iw = iw + 1
	    end if
	  else
	    if( yn .le. y0 ) then
	      if( area .lt. 0. ) iw = iw - 1
	    end if
	  end if

	  xo = xn
	  yo = yn
	end do

	inpoly1 = iw .ne. 0

	end

!**************************************************************

	function inpoly0(n,x,y,x0,y0)

! checks if point (x0,y0) is in general polygon (border is inside)
!
! use classical winding number

	implicit none

        logical inpoly0  !true if (x0,y0) is inside closed line (return)
	integer n        !total number of points in line
	real x(n),y(n)   !coordinates of line
	real x0,y0       !coordinates of point to check

	integer i
	real xm,ym,xn,yn
	real xmin,xmax,ymin,ymax
	logical bdebug

	logical isconvex,inconvex
	integer winding

	bdebug = .false.

	if(bdebug) write(6,*) 'debug is set in inpoly...'

!-------------------------------------------------------
! initialize inpoly
!-------------------------------------------------------

	inpoly0 = .false.

	if( n .le. 0 ) return

!-------------------------------------------------------
! maybe polygon is convex ?
!-------------------------------------------------------

	if( isconvex(n,x,y) ) then
	  inpoly0 = inconvex(n,x,y,x0,y0)
	  return
	end if

!-------------------------------------------------------
! no it is not -> get min/max and make first decisions
!-------------------------------------------------------

	xmin = x(1)
	xmax = x(1)
	ymin = y(1)
	ymax = y(1)

	do i=2,n
	  xmin = min(xmin,x(i))
	  xmax = max(xmax,x(i))
	  ymin = min(ymin,y(i))
	  ymax = max(ymax,y(i))
	end do

	if( bdebug ) then
	  write(6,*) n,x0,y0
	  write(6,*) xmin,xmax,ymin,ymax
	end if

	if( x0 .lt. xmin .or. x0 .gt. xmax ) return
	if( y0 .lt. ymin .or. y0 .gt. ymax ) return

!-------------------------------------------------------
! not yet decided -> proceed with heavy algorithm
!-------------------------------------------------------

	if( winding(n,x,y,x0,y0) .eq. 0 ) return

!-------------------------------------------------------
! ok, it is inside
!-------------------------------------------------------

	inpoly0 = .true.

!-------------------------------------------------------
! end of routine
!-------------------------------------------------------

	end

!**************************************************************

	function winding(n,x,y,x0,y0)

! compute winding number for x0/y0 (integer)

	implicit none

	integer winding
	integer n
	real x(n),y(n)
	real x0,y0

	integer i
	logical bdebug
	real xm,ym,xn,yn
	double precision dxn,dyn,dxm,dym
	double precision scal,cros,modu,arg
	double precision totangle,angle
	real pi2

	bdebug = .false.

	pi2 = 2. * 4. * atan(1.)
	totangle = 0.

	xn=x(n)
	yn=y(n)

	do i=1,n
	  xm=x(i)
	  ym=y(i)

	  dxn = xn-x0
	  dyn = yn-y0
	  dxm = xm-x0
	  dym = ym-y0

	  scal = dxn * dxm + dyn * dym
	  cros = dxn * dym - dxm * dyn
	  modu = sqrt( (dxn**2 + dyn**2) * (dxm**2 + dym**2) )
	  arg = scal / modu

	  if( arg .gt.  1. ) arg =  1.
	  if( arg .lt. -1. ) arg = -1.

	  angle = acos( arg )	! arg in [-1,+1], acos in [0,pi]

	  if( cros .lt. 0. ) angle = -angle

	  totangle = totangle + angle

	  xn = xm
	  yn = ym
	end do

	winding = nint( totangle / pi2 )

	if( bdebug ) then
	  write(6,*) 'winding... ',totangle,winding
	end if

	end

!**************************************************************
!**************************************************************
!**************************************************************

	subroutine c2p_ocean(u,v,s,d)

! converts cartesian to polar coordinates
! oceanographic convention: 0 -> current to north

	real u,v,s,d

	call c2p(u,v,s,d)

        d = d + 180.
        if( d > 360. ) d = d - 360.

	end

!**************************************************************

	subroutine c2p(u,v,s,d)

! converts cartesian to polar coordinates
!
! meteorological convention: 0 -> wind from north
! if oceanographic convention is needed: d=d+180; if(d>360) d=d-360

	implicit none

	real u,v
	real s,d

	real a
        real, parameter :: rad = 45 / atan(1.)

!----------------------------------------------------------------
! compute speed
!----------------------------------------------------------------

	s = sqrt( u**2 + v**2 )

!----------------------------------------------------------------
! compute mathematical angle
!----------------------------------------------------------------

	if( u .eq. 0. ) then
	  if( v .gt. 0. ) then
	    a = 90.
	  else
	    a = -90.
	  end if
	else
	  a = rad * atan(v/u)
	end if

	if( u .lt. 0. ) a = a + 180.

!----------------------------------------------------------------
! convert to wind directions
!----------------------------------------------------------------

	d = 270 - a

	if( d .gt. 360. ) d = d - 360.
	if( d .lt. 0.   ) d = d + 360.

!----------------------------------------------------------------
! end of routine
!----------------------------------------------------------------

	end

!**************************************************************

	subroutine p2c(s,d,u,v)

! converts polar to cartesian coordinates

	implicit none

	real s,d,u,v

	real a
        real, parameter :: rad = atan(1.) / 45.

        a = d
        a = 90. - a + 180.
        do while( a .lt. 0. )
          a = a + 360.
        end do
        a = mod(a,360.)

        u = s * cos( rad * a )
        v = s * sin( rad * a )

	end

!**************************************************************

        subroutine polar2xy(n,speed,dir,uv,vv)

        implicit none

        integer n
        real speed(n), dir(n)
        real uv(n), vv(n)

        integer i

        do i=1,n
	  call p2c(speed(i),dir(i),uv(i),vv(i))
        end do

        end

!**************************************************************

        subroutine xy2polar(n,uv,vv,speed,dir)

        implicit none

        integer n
        real uv(n), vv(n)
        real speed(n), dir(n)

        integer i

        do i=1,n
	  call c2p(uv(i),vv(i),speed(i),dir(i))
	end do

	end

!**************************************************************
!**************************************************************
!**************************************************************
! test routines for ploar conversion
!**************************************************************
!**************************************************************
!**************************************************************

	subroutine test_polar

	implicit none

	integer, parameter :: ndim = 100000
	real, parameter :: eps = 1.e-6

	integer i,ierr
	real u,v,uu,vv,s,d,dmax
	real du,dv,dd

	ierr = 0
	dmax = 0.

	call test_polar_local(1.,0.)
	call test_polar_local(1.,1.)
	call test_polar_local(-1.,1.)

	do i=1,ndim
	  call random_number(u)
	  call random_number(v)
	  u = 2*u - 1.
	  v = 2*v - 1.
	  call c2p(u,v,s,d)
	  call p2c(s,d,uu,vv)
	  du = abs(u-uu)
	  dv = abs(v-vv)
	  dd = max(du,dv)
	  dmax = max(dmax,dd)
	  if( dd > eps ) then
	    ierr = ierr + 1
	    write(6,'(6f10.2,e14.6)') u,v,s,d,uu,vv,dd
	  end if
	  !write(6,'(6f10.2,e14.6)') u,v,s,d,uu,vv,dd
	end do

	write(6,*) eps,dmax,ierr,ndim

	end

!**************************************************************

	subroutine test_polar_local(u,v)

	implicit none

	real u,v

	real uu,vv,s,d

	call c2p(u,v,s,d)
	call p2c(s,d,uu,vv)

	write(6,'(6f10.2)') u,v,s,d,uu,vv
	
	end

!**************************************************************
!	program main_test_polar
!	call test_polar
!	end
!**************************************************************
	
