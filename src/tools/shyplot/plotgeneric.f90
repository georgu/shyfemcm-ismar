
!--------------------------------------------------------------------------
!
!    Copyright (C) 2001,2003,2010,2014,2018-2019  Georg Umgiesser
!    Copyright (C) 2012  Debora Bellafiore
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

!  general utilities for plot routines
! 
!  revision log :
! 
!  16.10.2001	ggu	new routine fpfeil()
!  03.09.2003	ggu	bug fix XBXT in fpfeil
!  05.10.2003	ggu	bug fix S_IS_0 in pfeil
!  23.03.2010	ggu	changed v6.1.1
!  26.03.2010	ggu	new routine fcmpfeil for distorted scales
!  29.09.2010	ggu	changed VERS_6_1_12
!  16.03.2012	dbf	bug fix in fcmpfeil (d was used but not set)
!  30.10.2014	ggu	changed VERS_7_0_4
!  18.12.2018	ggu	changed VERS_7_5_52
!  13.03.2019	ggu	changed VERS_7_5_61
!  21.05.2019	ggu	changed VERS_7_5_62
!  13.11.2025	ggu	new routine poctogon() (not yet finished)
! 
! **************************************************

	subroutine pbox(xmin,ymin,xmax,ymax)

!  plots box around rectangle with actual color

	implicit none

	real xmin,ymin,xmax,ymax

        call qline(xmin,ymin,xmin,ymax)
        call qline(xmin,ymax,xmax,ymax)
        call qline(xmax,ymax,xmax,ymin)
        call qline(xmax,ymin,xmin,ymin)

	end

! *************************************************************

        subroutine pwbox(mode,xmin,ymin,xmax,ymax)

!  plots box with white background (erasing)
! 
!  mode != 0 -> plot outline

        implicit none

        integer mode
        real xmin,ymin,xmax,ymax

        call qwhite(1)
        call qgray(1.)
        call qrfill(xmin,ymin,xmax,ymax)
        call qwhite(0)
        call qgray(0.)

        if( mode .ne. 0 ) call pbox(xmin,ymin,xmax,ymax)

        end

! **************************************************

	subroutine pcircle(x,y,r)

!  plots circle around (x,y) with radius r

	implicit none

	real x,y,r

	call qarc(x,y,r,0.,360.)

	end

! **************************************************

	subroutine pcirclefill(x,y,r)

!  fill circle around (x,y) with radius r with actual color

	implicit none

	real x,y,r

	call qarcf(x,y,r,0.,360.)

	end

! **************************************************

	subroutine pcross(x,y,r)

!  plots cross at (x,y) with size r

	implicit none

	real x,y,r

	real dxy
        real, parameter :: sq2 = 1./sqrt(2.)

	dxy = r * sq2

	call qline(x-dxy,y-dxy,x+dxy,y+dxy)
	call qline(x+dxy,y-dxy,x-dxy,y+dxy)

	end 

! **************************************************

	subroutine make_octogon(x,y,r,xx,yy)

! makes coordinates of octogon

	implicit none

	real x,y
	real r
	real xx(8),yy(8)

        real, parameter :: sq2 = 1./sqrt(2.)
        real rsq2

        rsq2 = r * sq2

        xx(1) = x + r
        yy(1) = y
        xx(2) = x + rsq2
        yy(2) = y + rsq2
        xx(3) = x 
        yy(3) = y + r
        xx(4) = x - rsq2
        yy(4) = y + rsq2
        xx(5) = x - r
        yy(5) = y
        xx(6) = x - rsq2
        yy(6) = y - rsq2
        xx(7) = x 
        yy(7) = y - r
        xx(8) = x + rsq2
        yy(8) = y - rsq2

	end 

! **************************************************

	subroutine poctogon(x,y,r)

! plots octogon with center at (x,y) and with size r

	implicit none

	real x,y,r

        real xx(8),yy(8)

	call make_octogon(x,y,r,xx,yy)
	call qpoly(8,xx,yy)

	end 

! **************************************************

	subroutine poctogonfill(x,y,r)

! fills octogon with center at (x,y) and with size r with actual color

	implicit none

	real x,y,r

        real xx(8),yy(8)

	call make_octogon(x,y,r,xx,yy)
	call qafill(8,xx,yy)

	end 

! **************************************************

        subroutine pfeil(x,y,u,v,s)

!  draws arrow
! 
!  x,y           coordinates of arrow base
!  u,v           size of arrow
!  s             scaling of vector

        implicit none

        real dl,cos,sin
        parameter (dl=0.1)
        parameter (cos=.866025,sin=.5)

        real x,y,u,v,s
        real xu,yv,d
        real x1,y1,x2,y2

        if(u.eq.0..and.v.eq.0.) return
        if(s.eq.0.) return                      !S_IS_0

        xu=x+u*s
        yv=y+v*s

        d=0.2*s

        x1=xu-d*(cos*u+sin*v)
        y1=yv-d*(-sin*u+cos*v)
        x2=xu-d*(cos*u-sin*v)
        y2=yv-d*(sin*u+cos*v)

        call qmove(x,y)
        call qplot(xu,yv)
        call qplot(x1,y1)
        call qplot(x2,y2)
        call qplot(xu,yv)

        end

! **************************************************

        subroutine fpfeil(xb,yb,xt,yt,s)

!  draws full arrow
! 
!  xb,yb         coordinates of arrow base
!  xt,yt         coordinates of arrow tip
!  s             size of tip (in fraction, 1 -> whole length)

        implicit none

        real xb,yb,xt,yt
        real s

        real cos,sin
        parameter (cos=.866025,sin=.5)  !cos(30), sin(30)

        real u,v
        real d
        real xs,ys
        real x1,y1,x2,y2
        real x3(3),y3(3)

        if( xb .eq. xt .and. yb .eq. yt ) return	!bug fix XBXT

        u = xt - xb
        v = yt - yb

        d=s

        x1=xt-d*(cos*u+sin*v)
        y1=yt-d*(-sin*u+cos*v)
        x2=xt-d*(cos*u-sin*v)
        y2=yt-d*(sin*u+cos*v)

        xs = 0.5 * ( x1 + x2 )   !this is the point on the bottom of the
        ys = 0.5 * ( y1 + y2 )   !...arrow, lying on the line

        call qmove(xb,yb)
        !call qplot(xt,yt)      !if line is too thick, this is not nice
        call qplot(xs,ys)

        x3(1) = xt
        y3(1) = yt
        x3(2) = x1
        y3(2) = y1
        x3(3) = x2
        y3(3) = y2

        call qafill(3,x3,y3)

        !call qplot(x1,y1)
        !call qplot(x2,y2)
        !call qplot(xt,yt)

        end

! **************************************************

        subroutine fcmpfeil(xb,yb,xt,yt,s)

!  draws full arrow - adjusts for distortion
! 
!  xb,yb         coordinates of arrow base
!  xt,yt         coordinates of arrow tip
!  s             size of tip (in fraction, 1 -> whole length)
! 		... 0: only line, no arrow tip, <0: do not plot anything
! 		...good value is 0.3

        implicit none

        real xb,yb,xt,yt
        real s

        real cos,sin
        parameter (cos=.866025,sin=.5)  !cos(30), sin(30)

        real u,v
        real d,dd,rl,ddd
        real xs,ys,xn,yn,xa,ya
        real x1,y1,x2,y2
        real x3(3),y3(3)
	real xcm,ycm

! 	---------------------------------------------
! 	check if we have to plot anything
! 	---------------------------------------------

        if( xb .eq. xt .and. yb .eq. yt ) return	!bug fix XBXT
	if( s .lt. 0 ) return				!negative -> no arrow

! 	---------------------------------------------
! 	first only line is plotted
! 	---------------------------------------------

        u = xt - xb
        v = yt - yb

        d = s
	xs = xt - d*u
	ys = yt - d*v

        call qmove(xb,yb)
        call qplot(xs,ys)

	if( d .le. 0 ) return

! 	---------------------------------------------
! 	we still have to plot the tip 
! 	(x1,y1) and (x2,y2) are lateral coords of tip
! 	---------------------------------------------

	call qcm(xcm,ycm)
	dd = xcm/ycm

! 	---------------------------------------------
! 	do all computations in plotting coordinates
! 	---------------------------------------------

	v = v * dd	!account for distortion

	xn = -v		!normal vector
	yn = u

	rl = 0.5 * d * sqrt( (u*u+v*v) / (xn*xn+yn*yn) )

	x1 = rl * xn
	y1 = rl * yn
	x2 = rl * xn
	y2 = rl * yn

! 	---------------------------------------------
! 	now get final points in world coordinates
! 	---------------------------------------------

	x1 = xs + x1
	y1 = ys + y1/dd
	x2 = xs - x2
	y2 = ys - y2/dd

! 	---------------------------------------------
! 	prepare and plot tip
! 	---------------------------------------------

        x3(1) = xt
        y3(1) = yt
        x3(2) = x1
        y3(2) = y1
        x3(3) = x2
        y3(3) = y2

        call qafill(3,x3,y3)

! 	---------------------------------------------
! 	end of routine
! 	---------------------------------------------

        end

! **************************************************

        subroutine convert_wind(s,d,u,v)

        implicit none

        real s,d,u,v

        real dir
        real pi,rad
        parameter(pi=3.14159,rad=pi/180.)

        dir = d
        dir = 90. - dir + 180.
        do while( dir .lt. 0. )
          dir = dir + 360.
        end do
        dir = mod(dir,360.)

        u = s*cos(rad*dir)
        v = s*sin(rad*dir)

        end subroutine convert_wind

! **************************************************

