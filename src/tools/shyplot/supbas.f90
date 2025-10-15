
!--------------------------------------------------------------------------
!
!    Copyright (C) 1999-2001,2003-2005,2009-2020  Georg Umgiesser
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

! routines for plotting basin
! 
! contents :
! 
! subroutine basinit			internal initialization
! subroutine determine_frame_dimensions finds dimension of frame
! subroutine determine_typical_length_scale   computes length scale for basin
! 
! subroutine bash(mode)			plots basin (shell)
! 
! subroutine plot_basin(mode)		plots basin
! subroutine plot_basin_number(mode)	plots basin with node/element numbers
! subroutine bpixel			writes the pixel size as a comment 
! subroutine plot_frame(mode)		plots frame around plot
!
! subroutine basmima(mode)		computes min/max of basin
! subroutine setbas(x0,y0,x1,y1)	sets min/max of basin to plot
! subroutine getbas(x0,y0,x1,y1)	gets min/max of basin to plot
! subroutine plot_boundline		plots boundary line for basin
! 
! subroutine handle_spherical		handles spherical coordinates
! subroutine spherical_fact(fact,afact) computes factors for spherical coords
!
! subroutine handle_reg_grid		handles plotting of regular grid
! subroutine compute_reg_grid(bverb,dreg,imicro) find best regular grid
! subroutine plot_reg_grid(ngrid,dist,gray) plots regular grid on top of basin
!
! subroutine handle_bw_frame		handles labeling of regular grid
! subroutine plot_bw_frame(imicro,x0,y0,dx,dy,xmin,xmax,ymin,ymax)
!					plot black/white frame around grid
! subroutine plot_box(ic,x1,y1,x2,y2)	plots rectangular box
! 
! subroutine setreg_grid(regpar)	sets regular grid values
! subroutine plot_regular_points	plots regular points
! subroutine plot_obs_points(mode,file,text) plots special points (obs/meteo)
! subroutine plot_plus(x,y,dx,dy)	plots a plus char
! subroutine plot_cross(x,y,dx,dy)	plots a cross char
! subroutine plot_islands(cgray)	plots islands in gray
! subroutine plot_outer_area(n,xa,ya)	plots outer area
!
!  revision log :
! 
!  16.02.1999	ggu	bpixel: write bounding box to ps file (as comment)
!  09.02.2000	ggu	use inboxdim to compute box to plot
!  12.06.2000	ggu	get gray value for basin mode 3 from str file
!  11.02.2001	ggu	routine to compute typical length scale
!  21.08.2003	ggu	occupy is called elsewhere
!  16.12.2004	ggu	changed reggrid to plot regular grid
!  02.03.2005	ggu	in bash: bug fix -> get size of grid if not given
!  12.06.2009	ggu	new routines to handle spherical coords. & regular grid
!  15.06.2009	ggu	call to reggrid() changeed -> pass in gray value
!  14.09.2009	ggu	new routine divdist()
!  22.02.2010	ggu	new routine bw_frame() to plot bw scale around plot
!  23.03.2010	ggu	changed v6.1.1
!  09.04.2010	ggu	bug fix in frac_pos() -> maybe compiler error
!  20.12.2010	ggu	changed VERS_6_1_16
!  17.05.2011	ggu	new routine basin_number()
!  31.05.2011	ggu	changed VERS_6_1_23
!  30.03.2012	ggu	changed VERS_6_1_51
!  30.08.2012	ggu	new routines to automatically label spherical grid
!  12.09.2012	ggu	changed VERS_6_1_57
!  24.10.2012	ggu	bug in labelling non spherical grid (returned -1)
!  02.05.2013	ggu	handle fact in spherical coords
!  02.05.2013	ggu	meteo point plotting (plot_meteo_points())
!  13.06.2013	ggu	bug fix in spherical_fact() -> set fact to 1
!  19.06.2013	ggu	changed VERS_6_1_66
!  13.12.2013	ggu	new mode=4 for plotting gray grid over scalar variable
!  28.01.2014	ggu	changed VERS_6_1_71
!  30.05.2014	ggu	new metpnt for meteo points, imicro computed
!  18.07.2014	ggu	changed VERS_7_0_1
!  13.10.2014	ggu	changed VERS_7_0_2
!  26.11.2014	ggu	changed VERS_7_0_7
!  05.12.2014	ggu	changed VERS_7_0_8
!  23.12.2014	ggu	changed VERS_7_0_11
!  15.01.2015	ggu	changed VERS_7_1_1
!  19.01.2015	ggu	changed VERS_7_1_3
!  10.02.2015	ggu	also plot other points, also regular points
!  26.02.2015	ggu	changed VERS_7_1_5
!  05.05.2015	ggu	changed VERS_7_1_10
!  10.07.2015	ggu	changed VERS_7_1_50
!  17.07.2015	ggu	changed VERS_7_1_80
!  20.07.2015	ggu	changed VERS_7_1_81
!  12.10.2015	ggu	fix in rround() to handle rmaster==0 case
!  19.02.2016	ggu	changed VERS_7_5_2
!  10.06.2016	ggu	changed VERS_7_5_13
!  17.06.2016	ggu	changed VERS_7_5_15
!  27.06.2016	ggu	changed VERS_7_5_16
!  30.09.2016	ggu	changed VERS_7_5_18
!  25.05.2017	ggu	changed VERS_7_5_28
!  11.07.2017	ggu	changed VERS_7_5_30
!  26.09.2017	ggu	changed VERS_7_5_32
!  14.11.2017	ggu	changed VERS_7_5_36
!  22.02.2018	ggu	changed VERS_7_5_42
!  19.04.2018	ggu	changed VERS_7_5_45
!  06.07.2018	ggu	changed VERS_7_5_48
!  18.12.2018	ggu	changed VERS_7_5_52
!  13.03.2019	ggu	changed VERS_7_5_61
!  21.05.2019	ggu	changed VERS_7_5_62
!  13.02.2020	ggu	rounding routines into new file subround.f
!  10.11.2021   ggu     avoid warning for stack size
!  14.10.2025   ggu     big overhaul of routines
! 
!  notes :
! 
!  reg_grid is regular grid plotted over basin (reggrd and reggry)
!
!  rules for regular bw_frame around plot:
! 
!  if reggrd > 0
!	-> use this value
!  else if reggrd == 0		
!	-> do not plot frame
!  else if reggrd == -1
!	-> automatically determine
!  end if
!
!  if no legend is requested
!	-> do not plot frame
!  else if legend is requested
! 	if spherical		
!		-> plot bw_frame and no scale/north
! 	else			
!		-> plot scale/north and no bw_frame
!	end if
!  end if
! 
! *************************************************************

	module mod_bash

	implicit none

	logical, save :: bbverb = .true.
	logical, save :: breg = .false.		!regular grid is available
	logical, save :: bbasin = .false.	!basin has been read
	logical, save :: ball = .false.		!plot whole regular grid

	real, save :: regp(7)
	real, save :: xmin,ymin,xmax,ymax	!plotting area

	end module mod_bash

! *************************************************************

	subroutine bash_set_regpar(regpar)

	use mod_bash

	implicit none

	real regpar(7)

	breg = .true.
	regp = regpar

	end

! *************************************************************

	subroutine bash_set_regular(bregular)

	use mod_bash

	implicit none

	logical bregular

	breg = bregular

	end

! *************************************************************

	subroutine bash_get_regular_minmax(xpmin,ypmin,xpmax,ypmax)

	use mod_bash

	implicit none

	real xpmin,ypmin,xpmax,ypmax

	integer nx,ny
	real dx,dy

	if( .not. breg ) then
	  stop 'error stop bash_get_regular_minmax: not regular'
	end if

        nx = nint(regp(1))
        ny = nint(regp(2))
        xpmin = regp(3)
        ypmin = regp(4)
        dx = regp(5)
        dy = regp(6)
	xpmax = xmin + (nx-1)*dx
	ypmax = ymin + (ny-1)*dy

	end

! *************************************************************

	subroutine bash_set_regall(bregall)

	use mod_bash

	implicit none

	logical bregall

	ball = bregall

	end

! *************************************************************

	subroutine bash_verbose(bverbose)

	use mod_bash

	implicit none

	logical bverbose

	bbverb = bverbose

	end

! *************************************************************
! *************************************************************
! *************************************************************

	subroutine basinit

!  internal initialization

	use mod_bash

	implicit none

	logical, save :: binit = .false.

	if( binit ) return

	binit = .true.

	call determine_frame_dimensions
	call determine_typical_length_scale

	end

! **************************************************************

	subroutine determine_frame_dimensions

! finds dimension of frame

	use mod_bash
	use basin

	implicit none

	integer nx,ny
	real dx,dy
	real x0,y0,x1,y1
	logical inboxdim

	if( basin_has_read_basin() ) then
	  bbasin = .true.
	  call basmima(0)		!compute exact dimensions as default
	  if( breg .and. ball ) then
	    call bash_get_regular_minmax(xmin,ymin,xmax,ymax)
	  end if
	else if ( breg ) then
	  call bash_get_regular_minmax(xmin,ymin,xmax,ymax)
	else
	  stop 'error stop determine_frame_dimensions: cannot find min/max'
	end if

	if( inboxdim(' ',x0,y0,x1,y1) ) then	!size of plotting is given
	  call setbas(x0,y0,x1,y1)
	end if
	  
	call getbas(xmin,ymin,xmax,ymax)

	end

! *************************************************************

	subroutine determine_typical_length_scale

!  computes typical length scale for basin
! 
!  dxygrd is used first
!  then typls is used if given
!  else it is computed from grid

	use basin, only : nkn,nel,ngr,mbw
	use mod_bash

	implicit none

	logical bdebug
	integer ie
	real area,ao
	real dist,typls			!typical length scale
	real fact,afact
	real dxygrd
	real dx,dy
	double precision acu

	real getpar,aomega_elem

	bdebug = .false.

	if( bbasin ) then
	  acu = 0.
	  do ie=1,nel
	    ao = aomega_elem(ie)
	    acu = acu + ao
	  end do
	  area = 24. * acu / nel	!a little bit bigger than one element
	  call spherical_fact(fact,afact)  !correct for spherical coordinates
	else if( breg ) then
	  area = (xmax-xmin)*(ymax-ymin)
	  dx = regp(5)
	  dy = regp(6)
	  area = dx*dy
	  afact = 100			!empirical
	else
	  stop 'error stop determine_typical_length_scale: impossible'
	end if

	dist = sqrt(area/afact)

	typls = getpar('typls')
	dxygrd = getpar('dxygrd')
	if( dxygrd .gt. 0. ) then
	  typls = dxygrd
	else if( typls .gt. 0. ) then
	  typls = typls
	else
	  typls = dist
	end if
	!typls = 0.25
	!typls = 0.012
	call putpar('typls',typls)

	if( bdebug ) then
	write(6,*) '=============================='
	write(6,*) 'ggguuu breg: ',breg,bbasin
	write(6,*) 'ggguuu area: ',area,afact,dxygrd
	write(6,*) 'ggguuu typls: ',typls,dist
	write(6,*) 'ggguuu xyminmax: ',xmin,ymin,xmax,ymax
	write(6,*) '=============================='
	end if

	if( bbverb ) then
	  write(6,*) 'typical length scale for basin: ',dist
	  write(6,*) 'typical length scale used     : ',typls
	end if

	end

! **************************************************************
! **************************************************************
! **************************************************************

	subroutine bash(mode)

! plots basin (shell)
!
! mode	
!	0: only scaling, annotation, and bw_frame
!	1: net  
!	2: boundary line, frame, reg_grid, and legend
!	3: net in gray (for bathymetry - use bgray)
!	4: net in gray (for scalar and velocities - use bsgray)
!
! call first with mode == 0 (scaling)
! then call with desired mode
! must be called with mode == 2 before closing the plot

	implicit none

	integer mode

	integer ifreg
	real x0,y0,x1,y1
	real x0leg,y0leg,x1leg,y1leg
	real dxygrd,x,y
	real cgray
	character*80 file
	real getpar
	logical bverb
	logical inboxdim
        logical is_spherical
	logical is_box_given

	bverb = .false.

!--------------------------------------
! initializing
!--------------------------------------

	cgray = getpar('cislnd')	!gray scale for islands

	call basinit

!--------------------------------------
! if mode == 0 -> do scaling
!--------------------------------------

	if( mode .eq. 0 ) then

	  !call basmima(0)
	  !if( inboxdim(' ',x0,y0,x1,y1) ) then	!size of plotting is given
	  !  call setbas(x0,y0,x1,y1)
	  !end if

	  !call getbas(x0,y0,x1,y1)		!bug fix 2.3.2005

	  !--------------------------------------
	  ! prepare regular grid
	  !--------------------------------------

	  dxygrd = getpar('dxygrd')
	  x = x0 + dxygrd/2.
	  y = y0 + dxygrd/2.
	  call setgeo(x,y,dxygrd,dxygrd,-999.)

	  call annote		!annotation
	  call plot_basin(0)	!scaling
	  call handle_bw_frame
	  call plot_islands(cgray)

	  return

	end if

!--------------------------------------
! plot
!--------------------------------------

        call qcomm('Plotting basin')
        call plot_basin(mode)

!--------------------------------------
! end of plot - now only labeling
!--------------------------------------

	if( mode .ne. 2 ) return

	!--------------------------------------
	! boundary line and islands
	!--------------------------------------

	call plot_boundline

	!--------------------------------------
	! frame
	!--------------------------------------

	call plot_frame(0)
	call handle_reg_grid

	!--------------------------------------
	! user defined legend
	!--------------------------------------

	call legplo

	!--------------------------------------
	! legend (north and scale)
	!--------------------------------------

	!if( inboxdim('leg',x0,y0,x1,y1) ) then	!write legend
	if( is_box_given('leg') ) then	!write legend
          if( is_spherical() ) then
	    if( bverb ) then
              write(6,*) 'coordinates are spherical' // &
     &			' ...no north and scale written'
	    end if
            return
	  else
	    if( inboxdim('leg',x0leg,y0leg,x1leg,y1leg) ) then
	      call scale_legend(x0leg,y0leg,x1leg,y1leg)
	    end if
          end if
	end if

	!--------------------------------------
	! special output
	!--------------------------------------

	call getfnm('metpnt',file)
	if( file .ne. ' ' ) then
	  call plot_obs_points(1,file,'meteo points')
	end if

	call getfnm('obspnt',file)
	if( file .ne. ' ' ) then
	  call plot_obs_points(2,file,'observation points')
	end if

	ifreg = nint(getpar('ifreg'))
	if( ifreg .gt. 0 ) then
	  call plot_regular_points
	end if

!--------------------------------------
! end of routine
!--------------------------------------

	end

! *************************************************************

	subroutine plot_basin(mode)

!  plots basin
! 
!  mode	
! 	0: only scaling  
! 	1: net  
! 	2: boundary  
! 	3: net in gray (for bathymetry - use bgray)
! 	4: net in gray (for scalar and velocities - use bsgray)

	use mod_geom
	use basin
	use mod_bash
	use mod_hydro_plot

	implicit none

	integer mode

	integer ie,kn,ii,k
	real gray

	real getpar

	call basinit

	if( mode .eq. 0 ) then				!scaling
	  call qworld(xmin,ymin,xmax,ymax)
	  call qrcfy
	  call handle_spherical
	  call bpixel
	else if( mode .eq. 2 ) then			!plot boundary
	  call qgray(0.)
	  gray = getpar('bbgray')
	  if( gray < 1. ) then
	   call qgray(gray)
	   do k=1,nkn
	    if( kantv(1,k) .ne. 0 ) then
	      kn=kantv(1,k)
	      if( kn .gt. k ) call qline(xgv(k),ygv(k),xgv(kn),ygv(kn))
	      kn=kantv(2,k)
	      if( kn .gt. k ) call qline(xgv(k),ygv(k),xgv(kn),ygv(kn))
	    end if
	   end do
	  end if
	else if( mode .gt. 0 .and. mode .le. 4 ) then	!plot grid
	  call qgray(0.)
	  if( mode .eq. 3 ) then
	    gray = getpar('bgray')
	    if( gray .lt. 0. ) return			!do not plot net
	    call qgray(gray)
	  end if
	  if( mode .eq. 4 ) then
	    gray = getpar('bsgray')
	    if( .not. basin_has_read_basin() ) return	!no basin read
	    if( gray .lt. 0. ) return			!do not plot net
	    call qgray(gray)
	  end if

	  do ie=1,nel
	   if( bplot(ie) ) then
	    kn = nen3v(3,ie)
	    call qmove(xgv(kn),ygv(kn))
	    do ii=1,3
	      kn = nen3v(ii,ie)
	      call qplot(xgv(kn),ygv(kn))
	    end do
	   end if
	  end do
	else
	  write(6,*) 'mode = ',mode
	  stop 'error stop basin: internal error'
	end if

	call qgray(0.)

	end

! *************************************************************

	subroutine plot_basin_number(mode)

!  plots basin with node and element numbers
! 
!  mode	1: node number   2: element number
! 	positive: external    negative: internal

	use mod_geom
	use basin
	use mod_bash

	implicit none

	integer mode

	integer ie,kn,ii,k,ke,iee,i
	real gray,x,y
	character*10 s

	integer ialfa,ipext,ieext

	write(6,*) '======================================='
	write(6,*) '======================================='
	write(6,*) '======================================='
	write(6,*) 'debug print of node and element numbers'
	write(6,*) 'please remove in production code'
	write(6,*) '======================================='
	write(6,*) '======================================='
	write(6,*) '======================================='

	call qgray(0.)
	call qtxts(12)

	if( abs(mode) .eq. 1 ) then

	  do k=1,nkn
	    ke = k
	    if( mode .gt. 0 ) ke = ipext(k)
	    x = xgv(k)
	    y = ygv(k)
	    i = ialfa(float(ke),s,-1,-1)
	    call qtext(x,y,s)
	  end do

	else if( abs(mode) .eq. 2 ) then

	  call qtxtcc(0,0)
	  do ie=1,nel
	    x = 0.
	    y = 0.
	    do ii=1,3
	      k = nen3v(ii,ie)
	      x = x + xgv(k)
	      y = y + ygv(k)
	    end do
	    x = x / 3.
	    y = y / 3.
	    iee = ie
	    if( mode .gt. 0 ) iee = ieext(ie)
	    i = ialfa(float(iee),s,-1,-1)
	    call qtext(x,y,s)
	  end do
	  call qtxtcc(-1,-1)
	else

	  write(6,*) mode
	  stop 'error stop plot_basin_number: mode'

	end if

	end

! *****************************************************************

	subroutine bpixel

!  writes the pixel size as a comment to PS file

	implicit none

	character*80 line
	real fact
	real xvmin,yvmin,xvmax,yvmax
	real dxy
	integer idxy
	integer ixmin,iymin,ixmax,iymax

	call qgetvp(xvmin,yvmin,xvmax,yvmax)

	fact = 28.3	!factor to translate from PS dimensions to pixel
	dxy = 1.	!translation of plot
	idxy = 5	!leave some space around plot

	ixmin = fact * ( xvmin + dxy ) - idxy
	iymin = fact * ( yvmin + dxy ) - idxy
	ixmax = fact * ( xvmax + dxy ) + idxy
	iymax = fact * ( yvmax + dxy ) + idxy

! 	write(6,*) 'message from basin() ...'
! 	write(6,*) '****** ',xvmin,yvmin,xvmax,yvmax
! 	write(6,*) '****** ',ixmin,iymin,ixmax,iymax

	write(line,'(a,4i8)') 'InternalBoundingBox' &
     &				,ixmin,iymin,ixmax,iymax
	call qcomm(line)

	end

! *****************************************************************

	subroutine plot_frame(mode)

!  plots frame around plot
! 
!  0: just box around plot (only mode allowed)

	use mod_bash

	implicit none

	integer mode
	integer iaux
	real pxareg,pyareg,pxereg,pyereg
	real plon0,plat0,dlon,dlat
	real alonmin,alatmin,alonmax,alatmax
	real alon,alat
	real x,y

	call basinit

	pxareg = xmin
	pyareg = ymin
	pxereg = xmax
	pyereg = ymax

	call qcomm('Plotting frame around basin')
	call qgray(0.)
	call qmove(pxareg,pyareg)
	call qplot(pxereg,pyareg)
	call qplot(pxereg,pyereg)
	call qplot(pxareg,pyereg)
	call qplot(pxareg,pyareg)

	if( mode .eq. 0 ) return

	stop 'error stop plot_frame: only mode == 0 is allowed'
	end

! *************************************************************

	subroutine basmima(mode)

!  computes min/max of basin
! 
!  mode	0: exact dimensions  1: larger dimensions

	use basin
	use mod_bash

	implicit none

	integer mode

	integer k
	real x,y
	real dist

	real rdist,rround

	!call basinit
	if( .not. basin_has_read_basin() ) then
	  xmin = 0.
	  ymin = 0.
	  xmax = 0.
	  ymax = 0.
	  return
	end if

	xmin = xgv(1)
	xmax = xgv(1)
	ymin = ygv(1)
	ymax = ygv(1)
	do k=2,nkn
	  x = xgv(k)
	  y = ygv(k)
	  if( x .gt. xmax ) xmax = x
	  if( x .lt. xmin ) xmin = x
	  if( y .gt. ymax ) ymax = y
	  if( y .lt. ymin ) ymin = y
	end do

	if( mode .eq. 0 ) return

	dist = rdist(xmin,ymin,xmax,ymax)

	xmin = rround(xmin,dist,-1)
	ymin = rround(ymin,dist,-1)
	xmax = rround(xmax,dist,+1)
	ymax = rround(ymax,dist,+1)

	end

! *************************************************************

	subroutine setbas(x0,y0,x1,y1)

!  sets min/max of basin to plot

	use mod_bash

	implicit none

	real x0,y0,x1,y1

	call basinit

	xmin = x0
	ymin = y0
	xmax = x1
	ymax = y1

	end

! *************************************************************

	subroutine getbas(x0,y0,x1,y1)

!  gets min/max of basin to plot

	use mod_bash

	implicit none

	real x0,y0,x1,y1

	x0 = xmin
	y0 = ymin
	x1 = xmax
	y1 = ymax

	end

! *************************************************************

	subroutine plot_boundline

!  plots boundary line for basin

	implicit none

	integer iflag,i
	real x,y
	character*80 bndlin
	integer, save :: npoints = 0
	real, save, allocatable :: xx(:), yy(:)
	integer, save, allocatable :: ifl(:)

	logical is_grd_file

	call basinit

	call getfnm('bndlin',bndlin)
	if( bndlin .eq. " " ) return

	if( npoints == 0 ) then
	  call read_all_lines(bndlin,npoints,xx,yy,ifl)
	  if( npoints <= 0 ) goto 99
	  allocate(xx(npoints),yy(npoints),ifl(npoints))
	  call read_all_lines(bndlin,npoints,xx,yy,ifl)
	end if

	call qcomm('plotting boundary line')
	call qgray(0.)

	do i=1,npoints
	  x = xx(i)
	  y = yy(i)
	  iflag = ifl(i)
	  if( iflag .eq. 1 ) then
	    call qmove(x,y)
	  else
	    call qplot(x,y)
	  end if
	end do

	return
   99	continue
	write(6,*) 'error reading boundary line file: ',trim(bndlin)
	stop 'error stop plot_boundline'
	end

! **************************************************************

	subroutine handle_spherical

!  handles spherical coordinates

	use mod_bash

	implicit none

	real fact,afact
	logical is_spherical

	call spherical_fact(fact,afact)

	call qfact(fact,1.0)

	if( bbverb ) then
          if( is_spherical() ) then
	    write(6,*) 'Using factor for spherical coordinates: ',fact
	  else
	    write(6,*) 'Using factor for coordinates: ',fact
	  end if
	end if

	end

! **************************************************************

	subroutine spherical_fact(fact,afact)

!  computes factors for for spherical coordinates

	implicit none

	real fact	!factor in x direction
	real afact	!area factor between geo and cartesian coordinates

	real geomile,onedeg
	parameter( geomile = 1855.325 , onedeg = 60.*geomile )

	real y,pi,rad
	real x0,y0,x1,y1

	real getpar
	logical is_spherical

	fact = 1.
	afact = 1.

        if( .not. is_spherical() ) return	!only for spherical

	call getbas(x0,y0,x1,y1)
	y = 0.5*(y0+y1)
	pi = 4.*atan(1.)
	rad = pi/180.

	fact = cos(y*rad)
	afact = fact * onedeg * onedeg

	end

! **************************************************************
! **************************************************************
! **************************************************************

	subroutine handle_reg_grid

!  handles plotting of regular grid

	use mod_bash

	implicit none

	integer ngrid
	real reggrd
	real reggry
	real dreg
	integer imicro

	real getpar

	if( bbverb ) write(6,*) 'starting handle_reg_grid...'

	reggrd = getpar('reggrd')
	reggry = getpar('reggry')

	call compute_reg_grid(bbverb,reggrd,imicro)

	if( reggrd > 0 .and. reggry < 1. ) then
	  ngrid = 0
	  call plot_reg_grid(ngrid,reggrd,reggry)
	end if

	if( bbverb ) write(6,*) 'ending handle_reg_grid...'

	end

! **************************************************************

	subroutine compute_reg_grid(bverb,dreg,imicro)

!  tries to find best regular grid spacing value
!
! is called both in handle_reg_grid() and handle_bw_frame()

	implicit none

	logical bverb
	real dreg		!is reggrd in calling programs
	integer imicro

	real x0,y0,x1,y1
	real dx,dy,dxy
	real dsreg

	logical is_spherical,is_box_given
	real rnext,rnextsub

	dsreg = 0
	imicro = 0

        if( .not. is_spherical() ) then		!only for non-spherical
	  if( dreg .lt. 0. ) dreg = 0.
	end if

	if( dreg == 0 ) return

	call getbas(x0,y0,x1,y1)

	dx = x1 - x0
	dy = y1 - y0
	dxy = min(dx,dy)	!use smaller side

	dxy = dxy/4.		!around 4 grid lines

	if( dreg < 0 ) dreg = rnext(dxy,1)	!compute only if not given
	if( dreg == 0 ) return
	dsreg = rnextsub(dreg)
	imicro = nint(dreg/dsreg)

	if( bverb ) then
	  write(6,*) 'reg,micro: ',dreg,dsreg,imicro
	  write(6,*) 'new reggrd: ',dreg,x0,x1,y0,y1
	end if

	end

! *************************************************************

	subroutine plot_reg_grid(ngrid,dist,gray)

!  plots regular grid on top of basin
! 
!  if ngrid > 0                          use ngrid
!  if ngrid <= 0 and dist > 0.           use dist
!  otherwise                             do nothing

	use mod_bash

	implicit none

        integer ngrid           !divide into ngrid parts
	real dist               !regular distance of grid
	real gray		!gray value

	real xdmin,ydmin,xdmax,ydmax
	real x,y,x0,y0,dx,dy
	integer n,nx,ny

	real rround

	call basinit

	xdmin = rround(xmin,dist,-1)
	xdmax = rround(xmax,dist,+1)
	ydmin = rround(ymin,dist,-1)
	ydmax = rround(ymax,dist,+1)

        if( ngrid .gt. 0 ) then
          nx = ngrid - 1
          ny = ngrid - 1
          x0 = xmin
          y0 = ymin
          dx = (xmax-xmin)/ngrid
          dy = (ymax-ymin)/ngrid
        else if( dist .gt. 0. ) then
	  nx = nint((xdmax-xdmin)/dist)
	  ny = nint((ydmax-ydmin)/dist)
          x0 = xdmin
          y0 = ydmin
          dx = dist
          dy = dist
        else
          return
        end if

	call qcomm('regular grid')
	call qgray(gray)

	do n=0,nx
	  x = x0 + n * dx
	  if( x .gt. xmin .and. x .lt. xmax ) then
	    call qline(x,ymin,x,ymax)
	  end if
	end do

	do n=0,ny
	  y = y0 + n * dy
	  if( y .gt. ymin .and. y .lt. ymax ) then
	    call qline(xmin,y,xmax,y)
	  end if
	end do

	call qgray(0.)

	end

! **************************************************************
! **************************************************************
! **************************************************************

	subroutine handle_bw_frame

!  handles labeling of frame
	
	use plot_fonts
	use mod_bash

	implicit none

	integer nx,ny,n,i,nc
	integer imicro
	integer nxymax
	real reggrd
	real xvmin,yvmin,xvmax,yvmax
	real xdmin,ydmin,xdmax,ydmax
	real x0,y0,dx,dy
	real dist,x,y
	real size,ftext,eps, fextra_space, mm_per_point
	character*10 string
	logical bdebug

	real getpar
	real rround
	integer ialfa

	bdebug = .true.
	bdebug = .false.

	eps = 1.e-5

	call basinit

	if( bbverb ) write(6,*) 'starting handle_bw_frame...'

	!fs_bw_frame (font size) is defined in module plot_fonts

	fextra_space = 0.4   !factor for space between text and outer line
	mm_per_point = 0.352777778              !mm per point of font size
	size = (fs_bw_frame * mm_per_point)/10.
	size = size + size * fextra_space   !space around plot for labeling [cm]													
	ftext  = 2.5	     !factor to shift text vertically (from inner line)	
	nxymax = 50	     !not more than these number of reg grids

	reggrd = getpar('reggrd')
	imicro = nint(getpar('regdst'))

	if( bdebug ) write(6,*) 'reggrd,imicro (1): ',reggrd,imicro
	call compute_reg_grid(bbverb,reggrd,imicro)	!check if automatic
	if( bdebug ) write(6,*) 'reggrd,imicro (2): ',reggrd,imicro

	if( reggrd .eq. 0. ) return

	call plot_frame(0)

	call qgetvp(xvmin,yvmin,xvmax,yvmax)

	xvmin = xvmin + size
	yvmin = yvmin + size
	xvmax = xvmax - size
	yvmax = yvmax - size

	call qsetvpnc(xvmin,yvmin,xvmax,yvmax)

!  here labeling

	dist = reggrd
	call frac_pos(dist,nc)
	if( nc .eq. 0 ) nc = -1
	if( bbverb ) write(6,*) 'handle_bw_frame: ',dist,nc,imicro

	xdmin = rround(xmin,dist,-1)
	xdmax = rround(xmax,dist,+1)
	ydmin = rround(ymin,dist,-1)
	ydmax = rround(ymax,dist,+1)

	nx = nint((xdmax-xdmin)/dist)
	ny = nint((ydmax-ydmin)/dist)
        x0 = xdmin
        y0 = ydmin
        dx = dist
        dy = dist

	!if( nx .gt. nxymax .or. ny .gt. nxymax ) goto 99
        !DWNH bug

	call qcomm('labeling bw_frame')
	call qfont('Times-Roman')
	call qgray(0.0)

	!call qtxts(9) !Font size defined in module plot_fonts now
	call qtxts(fs_bw_frame)

	do n=0,nx
	 x = x0 + n * dx
	 !if( x .ge. xmin .and. x .le. xmax ) then
	 if( x-xmin .ge. -eps .and. x-xmax .le. eps ) then
	  i = ialfa(x,string,nc,-1)
	  call qtxtcr(0.,-ftext)
	  call qtext(x,ymax,string(1:i))
	  call qtxtcr(0.,+ftext)
	  call qtext(x,ymin,string(1:i))
	 end if
	end do

	do n=0,ny
	 y = y0 + n * dy
	 !if( y .ge. ymin .and. y .le. ymax ) then
	 if( y-ymin .ge. -eps .and. y-ymax .le. eps ) then
	  i = ialfa(y,string,nc,-1)
	  call qtxtr(90.)
	  call qtxtcr(0.,-ftext)
	  call qtext(xmin,y,string(1:i))
	  call qtxtr(-90.)
	  call qtxtcr(0.,-ftext)
	  call qtext(xmax,y,string(1:i))
	 end if
	end do

	call qtxtcr(-1.,-1.)
	call qtxtr(0.)
	call qgray(0.)

	call plot_bw_frame(imicro,x0,y0,dx,dy,xmin,xmax,ymin,ymax)

	call qsetvp(xvmin,yvmin,xvmax,yvmax)

	if( bbverb ) write(6,*) 'ending handle_bw_frame...'

	return
   99	continue
	write(6,*) 'nx,ny: ',nx,ny
	stop 'error stop handle_bw_frame: nx,ny too high'
	end

! **************************************************************

	subroutine plot_bw_frame(imicro,x0,y0,dx,dy,xmin,xmax,ymin,ymax)

!  plot black/white frame around geographical grid

	implicit none

	integer imicro
	real x0,y0,dx,dy,xmin,xmax,ymin,ymax

	integer i
	real xcm,ycm
	real hx,hy
	real ddx,ddy
	real x,y,xx0,yy0
	real xn,yn
	real fact

! ------------------------------------------------------------
!  set parameters
! ------------------------------------------------------------

	if( dx .le. 0. .or. dy .le. 0. ) goto 97

	call qcm(xcm,ycm)
	hx = 0.1*xcm
	hy = 0.1*ycm

	fact = 1.5

! 	here we should define how many boxes inbetween labels -> empirical

	if( imicro .le. 0 ) then
	  ddx = dx/2.
	  if( ddx .gt. fact*xcm ) ddx = dx/4.
	  ddy = dy/2.
	  if( ddy .gt. fact*ycm ) ddy = dy/4.
	else
	  ddx = dx/imicro
	  ddy = dy/imicro
	end if

! ------------------------------------------------------------
!  plot along x-axis
! ------------------------------------------------------------

	x = x0
	do while( x .lt. xmin )
	  x = x + dx
	end do
	xx0 = x

	x = x0
	if( imicro .ge. 0 ) x = 2.*xmax
	do while( x .le. xmax )		!plot ticks on x axis
	  call qline(x,ymax,x,ymax+hx)
	  call qline(x,ymin,x,ymin-hx)
	  x = x + dx
	end do

	x = xx0
	if( imicro .lt. 0 ) x = xmax
	i = 0
	do while( x .lt. xmax )
	  i = mod(i+1,2)
	  xn = min(x+ddx,xmax)
	  call plot_box(i,x,ymax,xn,ymax+hy)
	  call plot_box(i,x,ymin,xn,ymin-hy)
	  x = x + ddx
	end do

	x = xx0
	if( imicro .lt. 0 ) x = xmin
	i = 1
	do while( x .gt. xmin )
	  i = mod(i+1,2)
	  xn = max(x-ddx,xmin)
	  call plot_box(i,xn,ymax,x,ymax+hy)
	  call plot_box(i,xn,ymin,x,ymin-hy)
	  x = x - ddx
	end do

! ------------------------------------------------------------
!  plot small corner boxes
! ------------------------------------------------------------

	if( imicro .ge. 0 ) then
	  call plot_box(0,xmin-hx,ymax,xmin,ymax+hy)
	  call plot_box(0,xmax,ymax,xmax+hx,ymax+hy)
	  call plot_box(0,xmin-hx,ymin-hy,xmin,ymin)
	  call plot_box(0,xmax,ymin-hy,xmax+hx,ymin)
	end if

! ------------------------------------------------------------
!  plot along y-axis
! ------------------------------------------------------------

	y = y0
	do while( y .lt. ymin )
	  y = y + dy
	end do
	yy0 = y

	y = y0
	if( imicro .ge. 0 ) y = 2.*ymax
	do while( y .le. ymax )		!plot ticks on y axis
	  call qline(xmin,y,xmin-hx,y)
	  call qline(xmax,y,xmax+hx,y)
	  y = y + dy
	end do

	y = yy0
	if( imicro .lt. 0 ) y = ymax
	i = 0
	do while( y .lt. ymax )
	  i = mod(i+1,2)
	  yn = min(y+ddy,ymax)
	  call plot_box(i,xmin-hx,y,xmin,yn)
	  call plot_box(i,xmax,y,xmax+hx,yn)
	  y = y + ddy
	end do

	y = yy0
	if( imicro .lt. 0 ) y = ymin
	i = 1
	do while( y .gt. ymin )
	  i = mod(i+1,2)
	  yn = max(y-ddy,ymin)
	  call plot_box(i,xmin-hx,yn,xmin,y)
	  call plot_box(i,xmax,yn,xmax+hx,y)
	  y = y - ddy
	end do

! ------------------------------------------------------------
!  end of routine
! ------------------------------------------------------------

	return
   97	continue
	write(6,*) 'dx,dy: ',dx,dy
	stop 'error stop plot_bw_frame: dx,dy'
	end

! **************************************************************

	subroutine plot_box(ic,x1,y1,x2,y2)

! plots rectangular box
!
! if ic == 1 fills rectangle

	implicit none

	integer ic
	real x1,y1,x2,y2

	if( ic .eq. 1 ) then
	  call qrfill(x1,y1,x2,y2)
	end if

	call qline(x1,y1,x2,y1)
	call qline(x2,y1,x2,y2)
	call qline(x2,y2,x1,y2)
	call qline(x1,y2,x1,y1)

	end

! **************************************************************
! **************************************************************
! **************************************************************

	subroutine plot_regular_points

! plots regular points

	use mod_bash

	implicit none

	integer nx,ny
	real ddx,ddy,dxy
	real x0,y0,x1,y1,dx,dy,flag
	real fact,afact
	real x,y
	real fplus,eps

	fplus = 0.3		!relative size of plus sign
	fplus = 0.1		!relative size of plus sign
	fplus = 0.05		!relative size of plus sign
	eps = 1.e-5

        nx = nint(regp(1))
        ny = nint(regp(2))
        x0 = regp(3)
        y0 = regp(4)
        dx = regp(5)
        dy = regp(6)
        flag = regp(7)

	if( dx <= 0. .or. dy <= 0. ) return
	if( nx <= 0. .or. ny <= 0. ) return

	write(6,*) 'plotting regular points: ',flag
	write(6,*) x0,y0,dx,dy
	write(6,*) xmin,ymin,xmax,ymax

	x1 = x0 + (nx-1)*dx
	y1 = y0 + (ny-1)*dy

	do while( x0 + dx < xmin )
	  x0 = x0 + dx
	end do
	do while( y0 + dy < ymin )
	  y0 = y0 + dy
	end do
	do while( x1 - dx > xmax )
	  x1 = x1 - dx
	end do
	do while( y1 - dy > ymax )
	  y1 = y1 - dy
	end do

	nx = (xmax-xmin)/dx
	ny = (ymax-ymin)/dy
	!nx = max(nx,20)
	!ny = max(ny,20)
	ddx = fplus*(xmax-xmin)/(nx-1)
	ddy = fplus*(ymax-ymin)/(ny-1)
	dxy = min(ddx,ddy)
	call spherical_fact(fact,afact)
	ddx = dxy/fact
	ddy = dxy

	write(6,*) nx,ny,ddx,ddy

	call qgray(0.)

	y = y0
	do while( y <= y1 + eps )
	  x = x0
	  do while( x <= x1 + eps )
	    call plot_plus(x,y,ddx,ddy)
	    x = x + dx
	  end do
	  y = y + dy
	end do

	end

! **************************************************************

	subroutine plot_obs_points(mode,file,text)

! plots special points from observation/meteo file
! 
! name of coords.dat has to passed into the subroutine
! if sea_land.dat exists it is used, otherwise we can do without

	implicit none

	integer mode
	character*(*) file
	character*(*) text

	integer ndim
	parameter (ndim=30000)

	integer nx,ny,nz,n
	integer idum,i,ios
	real xx,yy
	real fact,afact
	real xmin,xmax,ymin,ymax
	real dx,dy,dxy
	real fplus
	real, allocatable, save :: x(:),y(:),rf(:)
	save n,dx,dy,dxy

	integer icall
	save icall
	data icall / 0 /

	if( icall .eq. -1 ) return

	if( icall .eq. 0 ) then
	  allocate(x(ndim))
	  allocate(y(ndim))
	  allocate(rf(ndim))
	  open(1,file=file,status='old',form='formatted',err=88)

	  if( mode .eq. 1 ) then
	    read(1,*) nx,ny,nz
	    n = nx*ny
	    if( n .gt. ndim ) goto 99
	    do i=1,n
	      read(1,*) idum,idum,x(i),y(i)
	    end do
	  else if( mode .eq. 2 ) then
	    i = 0
	    do
	      !read(1,*,iostat=ios) xx,yy
	      read(1,*,iostat=ios) idum,idum,xx,yy
	      if( ios < 0 ) exit
	      if( ios > 0 ) goto 98
	      i = i + 1
	      if( i .gt. ndim ) goto 99
	      x(i) = xx
	      y(i) = yy
	    end do
	    n = i
	    nx = 10
	    ny = 10
	  else
	    write(6,*) 'mode = ',mode
	    stop 'error stop plot_obs_points: mode not recognized'
	  end if

	  close(1)

	  call mima(x,n,xmin,xmax)
	  call mima(y,n,ymin,ymax)

	  fplus = 0.3
	  fplus = 0.1
	  dx = fplus*(xmax-xmin)/(nx-1)
	  dy = fplus*(ymax-ymin)/(ny-1)
	  dxy = min(dx,dy)
	  call spherical_fact(fact,afact)
	  dx = dxy/fact
	  dy = dxy

	  do i=1,nx*ny
	    rf(i) = 1.
	  end do
	  open(1,file='sea_land.dat',status='old',form='formatted',err=77)
	  read(1,*)
	  read(1,*) nx,ny,nz
	  if( n .ne. nx*ny ) goto 97
	  read(1,*) (rf(i),i=1,n)
	  close(1)
   77	  continue
	  icall = 1
	end if

	write(6,*) 'plotting obs points'
	call qcomm('plotting obs points')
	call qgray(0.)
	!dx = 0.02
	!dy = dx
	do i=1,n
	  if( rf(i) .gt. 0 ) then
	    call qgray(0.)
	    call plot_plus(x(i),y(i),dx,dy)
	  else
	    call qgray(0.5)
	    call plot_cross(x(i),y(i),dx,dy)
	  end if
	end do

	return
   88	continue
	icall = -1
	write(6,*) 'message from plot_obs_points:'
	write(6,*) 'no coords.dat file ... cannot plot coordinates'
	write(6,*) 'file: ',file
	write(6,*) 'trying to plot '//trim(text)
	stop 'error stop plot_obs_points: file'
   98	continue
	write(6,*) 'file: ',trim(file)
	write(6,*) 'mode,irec: ',mode,i
	write(6,*) 'trying to plot '//trim(text)
	stop 'error stop plot_obs_points: read error'
   99	continue
	write(6,*) n,ndim
	write(6,*) 'trying to plot '//trim(text)
	stop 'error stop plot_obs_points: ndim'
   97	continue
	write(6,*) nx,ny,n,ndim
	write(6,*) 'trying to plot '//trim(text)
	stop 'error stop plot_obs_points: different length'
	end

! **************************************************************

	subroutine plot_plus(x,y,dx,dy)

! plots a plus char

	implicit none

	real x,y,dx,dy

	call qline(x,y-dy,x,y+dy)
	call qline(x-dx,y,x+dx,y)

	end

! **************************************************************

	subroutine plot_cross(x,y,dx,dy)

! plots a cross char

	implicit none

	real x,y,dx,dy

	call qline(x-dx,y-dy,x+dx,y+dy)
	call qline(x-dx,y+dx,x+dx,y-dx)
	
	end

! **************************************************************

	subroutine plot_islands(cgray)

! plots islands in gray

	use mod_geom
	use basin
	use mod_bash


	implicit none

	real cgray	!color
	
	logical bouter
	integer k,kstart,kk,kn,ko
	integer n,nis
	real x0,y0,area,cg
	real xa(nkn),ya(nkn)

	real areapoly

	if( cgray < 0 ) return

	cg = cgray
	bouter = .false.		! plot external island
	if( cg >= 2 ) then
	  cg = cg - 2
	  bouter = .true.
	end if

	call qgray(cg)

	x0 = 0.5*(xmax-xmin)
	y0 = 0.5*(ymax-ymin)

	nis = 0

	do k=1,nkn
	  if( kantv(1,k) > 0 ) then
	    kstart = k
	    kk = k
	    n = 0
	    nis = nis + 1
	    do
	      kn = kantv(1,kk)
	      kantv(1,kk) = -kantv(1,kk)
	      ko = kk
	      kk = kn
	      if( kn == 0 ) stop 'error stop plot_islands: internal error'
	      if( kn < 0 ) exit
	      n = n + 1
	      xa(n) = xgv(kk)
	      ya(n) = ygv(kk)
	    end do
	    area = areapoly(n,xa,ya)
	    if( kstart /= ko ) then
		stop 'error stop plot_islands: internal error'
	    end if
	    !write(69,*) 'debug: ',kk,kn,ko,kstart
	    !write(69,*) 'island found: ',nis,n,area
	    write(6,*) 'island found: ',nis,n,area
	    if( area < 0 ) then	!real island
	      call qafill(n,xa,ya)
	    else if( bouter ) then
	      call plot_outer_area(n,xa,ya)
	    end if
	  end if
	end do

	do k=1,nkn
	  if( kantv(1,k) < 0 ) kantv(1,k) = -kantv(1,k)
	  if( kantv(2,k) < 0 ) kantv(2,k) = -kantv(2,k)
	end do

	call qgray(0.)

	end

! **************************************************************

	subroutine plot_outer_area(n,xa,ya)

! plots outer area

	use mod_bash

	implicit none

	integer n
	real xa(n),ya(n)

	integer i,ilow
	real dx,dy,xlow,ylow
	real xn(n+7),yn(n+7)

	dx = xmax - xmin
	dy = ymax - ymin

!	-------------------------------------------
!	look for lowest y value in ya
!	-------------------------------------------

	ilow = 0
	ylow = ya(1) + 1.
	do i=1,n
	  if( ya(i) < ylow ) then
	    ylow = ya(i)
	    ilow = i
	  end if
	end do
	xlow = xa(ilow)

!	-------------------------------------------
!	in (xlow,ylow) are coordinates for lowest y value
!	-------------------------------------------

!	-------------------------------------------
!	copy points until (xlow,ylow)
!	-------------------------------------------

	do i=1,ilow
	  xn(i) = xa(i)
	  yn(i) = ya(i)
	end do
	
!	-------------------------------------------
!	add points to comprise whole domain
!	-------------------------------------------

	i = ilow

	i = i + 1
	xn(i) = xlow
	yn(i) = ymin - dy

	i = i + 1
	xn(i) = xmin - dx
	yn(i) = ymin - dy

	i = i + 1
	xn(i) = xmin - dx
	yn(i) = ymax + dy

	i = i + 1
	xn(i) = xmax + dx
	yn(i) = ymax + dy

	i = i + 1
	xn(i) = xmax + dx
	yn(i) = ymin - dy

	i = i + 1
	xn(i) = xlow
	yn(i) = ymin - dy

!	-------------------------------------------
!	add remaining points
!	-------------------------------------------

	do i=ilow,n
	  xn(i+7) = xa(i)
	  yn(i+7) = ya(i)
	end do

!	-------------------------------------------
!	fill the area (outer area of basin)
!	-------------------------------------------

	call qafill(n+6,xn,yn)

!	-------------------------------------------
!	end of routine remaining points
!	-------------------------------------------

	end

! **************************************************************

