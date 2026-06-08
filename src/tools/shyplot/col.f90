
!--------------------------------------------------------------------------
!
!    Copyright (C) 2003-2004,2009-2011,2013-2019  Georg Umgiesser
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

!  include file for colors
! 
!  revision log :
! 
!  05.08.2003	ggu	changed adjcoltab, new cfast, init_cfast
!  20.08.2003	ggu	some more comments
!  21.08.2003	ggu	colini is called in colsetup
!  04.10.2004	ggu	adjust amin and decimal point in colauto
!  21.04.2009	ggu	allow for constant values (make 2 isolines)
!  23.03.2010	ggu	changed v6.1.1
!  09.12.2011	ggu	changed VERS_6_1_38
!  13.06.2013	ggu	small bug fix for computation of dval with rnext
!  01.10.2013	ggu	bug fix for valmin/valmax too close
!  25.10.2013	ggu	changed VERS_6_1_68
!  18.07.2014	ggu	changed VERS_7_0_1
!  06.05.2015	ggu	logarithmic scale introduced (ilog, blog)
!  21.05.2015	ggu	changed VERS_7_1_11
!  16.11.2015	ggu	changed VERS_7_3_14
!  18.12.2015	ggu	changed VERS_7_3_17
!  30.05.2016	ggu	changed VERS_7_5_11
!  15.06.2016	ggu	adjcoltab deleted
!  30.09.2016	ggu	changed VERS_7_5_18
!  13.04.2017	ggu	changed VERS_7_5_25
!  11.07.2017	ggu	changed VERS_7_5_30
!  07.12.2017	ggu	changed VERS_7_5_40
!  24.01.2018	ggu	changed VERS_7_5_41
!  19.04.2018	ggu	changed VERS_7_5_45
!  31.08.2018	ggu	changed VERS_7_5_49
!  18.12.2018	ggu	changed VERS_7_5_52
!  21.05.2019	ggu	changed VERS_7_5_62
!  20.09.2022	ggu	in colauto() use absolute value in case val is constant
!  18.09.2024	ggu	new parameter rfaccol
! 
! **********************************************************
! 
! ---------------------------------------------------------------------
! 
!  typical usage:
! 
!  $color
!          x0col = 0.1  y0col = 0.65
!          x1col = 0.45 y1col = 0.7
!   
!          legcol = 'Bathymetry [m]'
!          ndccol = -1  faccol = 1
!          icolor = 1
!   
!          colmin = 0   colmax = 0.666
!          valmin = 0   valmax = 300
!          niso = 50   nisomx = 70
!          nctick = 5
!  $end
! 
!  other variables:
! 
! 	 isoval =    1.   3.   5.
! 	 color  = 0.3  0.4  0.5  0.6
! 
! **********************************************************

!==========================================================
	module color
!==========================================================

	implicit none

        integer, parameter :: isodim = 256   !max number of isolines
        integer, parameter :: coldim = 1024  !max number of color table entries

        integer, save :: isopar          !dimension of arrays ( = isodim )
        integer, save :: isoanz          !number of isolines used
        integer, save :: nisord          !number of isolines read
        integer, save :: ncolrd          !number of colors read
        integer, save :: icauto          !automatic computation of iso-values
        integer, save :: nriso           !number of single isolines to plot
        integer, save :: iusear          !use array looking up values for colors

        real, save :: fnull              !flag for value not to plot

        real, save :: fiso(isodim)       !array with iso-values
        real, save :: ciso(isodim+1)     !array with colors
        real, save :: riso(isodim)       !array with iso-values of isolines

        integer, save :: icmax = 0       !filling of colortable
        real, save :: coltab(3,coldim)   !custom colortable

        character*80, save :: colfil = ' '

	integer, save :: icoltab_def = 1 !default color table (gray values)
	real, save :: color_def = 0. !default color 

!==========================================================
	end module color
!==========================================================

	subroutine colsetup

	use color

!  sets up color table after read from $color section

	logical bdebug
	integer niso,nisomx,nisodf
	integer icolor,ipllog
	real rfiso
	real colmin,colmax,valmin,valmax,dval
	logical bisord,bcolrd,bfunc,blog

	real getpar

	bdebug = .true.
	bdebug = .false.

! -------------------------------------------------
!  initialize value and check for error
! -------------------------------------------------

	call colini
	call colchk			!finishes set up of color section

	icauto = 0

	colmin = getpar('colmin')	!min color [0..1]
	colmax = getpar('colmax')	!max color [0..1]
	valmin = getpar('valmin')	!min isovalue
	valmax = getpar('valmax')	!max isovalue
	dval   = getpar('dval')		!distance of isovalues
	rfiso  = getpar('rfiso')	!interpolating function
	ipllog = nint(getpar('ipllog'))	!logarithmic scale
	niso   = nint(getpar('niso'))	!total number of isovalues to use
	nisomx = nint(getpar('nisomx'))	!maximum number of isovalues allowed
	nisodf = nint(getpar('nisodf'))	!default number of isovalues to use
	icolor = nint(getpar('icolor'))	!color table

	if( bdebug ) then
	  write(6,*) 'colsetup'
	  write(6,*) 'colmin/max: ',colmin,colmax
	  write(6,*) 'dval,valmin/max: ',dval,valmin,valmax
	  write(6,*) 'niso,nisomx,nisodf: ',niso,nisomx,nisodf
	  write(6,*) 'nisord,ncolrd: ',nisord,ncolrd
	end if

	if( colmin .eq. colmax ) then
	  write(6,*) 'colmin/max : ',colmin,colmax
	  write(6,*) 'values must be different...'
	  stop 'error stop colsetup: colmin/max'
	end if

	bisord = nisord .gt. 0		!isovalues have been read
	bcolrd = ncolrd .gt. 0		!colors have been read
	bfunc = rfiso .ne. 0		!power scale
	blog = ipllog .ne. 0		!logarithmic scale
	iusear = 0			!use array for lookup
	if( bfunc .or. blog ) iusear = 1

! -----------------------------------------------------
!  see what has been read -> compute niso
! -----------------------------------------------------

	if( bisord .or. bcolrd ) then

! 	  --------------------------------------------
! 	  at least one type of values has been read
! 	  --------------------------------------------

	  if( bisord .and. bcolrd ) then
	    if( nisord + 1 .ne. ncolrd ) then
		write(6,*) 'Number of isolines and colors read'
		write(6,*) '  are not compatible. There must be'
		write(6,*) '  one color more than isovalues.'
		write(6,*) 'nisord, ncolrd : ',nisord,ncolrd
		stop 'error stop colsetup: nisord/ncolrd'
	    end if
	  end if

! 	  --------------------------------------------
! 	  set niso
! 	  --------------------------------------------

	  if( bisord ) niso = nisord
	  if( bcolrd ) niso = ncolrd - 1

! 	  --------------------------------------------
! 	  check computed number of isovalues
! 	  --------------------------------------------

	  if( niso .le. 0 ) then
	    write(6,*) 'Cannot use only one color.'
	    stop 'error stop colsetup: niso'
	  else if( niso .eq. 1 ) then
	    if( .not. bisord .and. valmax .ne. valmin ) then
	      write(6,*) 'One isovalue only possible for valmin=valmax'
	      stop 'error stop colsetup: niso'
	    end if
	  end if

! 	  --------------------------------------------
! 	  if we cannot determine isovalues -> automatic
! 	  --------------------------------------------

	  if( .not. bisord .and. valmax .eq. valmin ) icauto = 1

	else if( valmin .ne. valmax ) then

! 	  --------------------------------------------
! 	  isovalues can be determined
! 	  --------------------------------------------

	  if( dval .ne. 0. ) then
	    niso = nint((valmax-valmin)/dval) + 1
	    if( niso .le. 0 ) then
	      write(6,*) 'Error computing niso.'
	      write(6,*) 'If dval < 0 then valmax < valmin.'
	      stop 'error stop colsetup: niso'
	    end if
	  else
	    if( niso .le. 0 ) niso = nisodf
	  end if

	else ! if( valmin .eq. valmax ) then

! 	  --------------------------------------------
! 	  here we can do only automatic evaluation
! 	  --------------------------------------------

	  if( dval .ne. 0. ) then
	    niso = 0
! 	  else
! 	    use given value or 0
	  end if

	  icauto = 1

	end if

! -----------------------------------------------------
!  check computed niso value and store
! -----------------------------------------------------

	if( niso .gt. nisomx ) then
	  write(6,*) 'Value of niso > nisomx.'
	  write(6,*) 'niso,nisomx : ',niso,nisomx
	  write(6,*) 'Either reduce iso-values or increase nisomx.'
	  stop 'error stop colsetup: nisomx'
	end if

	call putpar('niso',float(niso))
	isoanz = niso

! -----------------------------------------------------
!  compute all other values that can be set up
! -----------------------------------------------------

	if( icauto .eq. 0 ) then

! 	  --------------------------------------------
! 	  no automatic evaluation -> set up everything
! 	  --------------------------------------------

	  if( niso .le. 0 ) stop 'internal error colsetup (1)'

	  if( .not. bcolrd ) then
	    call mkval(niso+1,ciso,colmin,colmax,0.,0)
	  end if
	  if( .not. bisord ) then
	    if( bfunc .or. blog ) then
	      call mkval(-niso,fiso,valmin,valmax,rfiso,ipllog)
	    else
	      call mkval(niso,fiso,valmin,valmax,dval,ipllog)
	    end if
	  end if

	else ! if( icauto .ne. 0 ) then

! 	  --------------------------------------------
! 	  automatic evaluation -> set up as much as possible
! 	  --------------------------------------------

	  if( niso .gt. 0 ) then

	    if( .not. bcolrd ) then
	      call mkval(niso+1,ciso,colmin,colmax,0.,0)
	    end if
	    if( bisord ) stop 'internal error colsetup (2)'
	    if( valmin .ne. valmax ) stop 'internal error colsetup (3)'

	  end if

	end if

! -----------------------------------------------------
!  debug output
! -----------------------------------------------------

	if( bdebug ) then
	  write(6,*) 'colsetup: ',icauto
	  write(6,*) 'colmin/max: ',colmin,colmax
	  write(6,*) 'dval,valmin/max: ',dval,valmin,valmax
	  write(6,*) 'niso,nisomx,nisodf: ',niso,nisomx,nisodf
	end if

! -----------------------------------------------------
!  end of routine
! -----------------------------------------------------

	end

! **********************************************************

	subroutine colauto(valmin,valmax)

	use color

!  sets up color table in automatic mode

	implicit none

	real valmin,valmax		!min/max for array of values

	logical bdebug,bfunc,blog
	integer i
	integer niso,nisomx,nisodf
	integer icolor
	real colmin,colmax,dval
	real amin,amax
        real fact,rfact,val,fdec,eps
        integer ndec,idec,ipllog
	real rfiso,ddval

	real getpar
	real rnext,rnexta

	bdebug = .true.
	bdebug = .false.

	amin = 0.
	amax = 0.

! -------------------------------------------------
!  check if we are in automatic mode
! -------------------------------------------------

	if( bdebug ) then
	  write(6,*) 'colauto: ',icauto
	  if( icauto .eq. 0 ) then
	    write(6,*) 'colors: ',isoanz+1
	    write(6,*) (ciso(i),i=1,isoanz+1)
	    write(6,*) 'isovals: ',isoanz
	    write(6,*) (fiso(i),i=1,isoanz)
	  end if
	end if

	if( icauto .lt. 0 ) stop 'error stop colauto: no initialization'
	if( icauto .eq. 0 ) return

! -------------------------------------------------
!  initialize values
! -------------------------------------------------

	colmin = getpar('colmin')	!min color [0..1]
	colmax = getpar('colmax')	!max color [0..1]
	dval   = getpar('dval')		!distance of isovalues
	niso   = nint(getpar('niso'))	!total number of isovalues
	rfiso  = getpar('rfiso')	!interpolating function
	ipllog = nint(getpar('ipllog'))	!logarithmic scale
	nisomx = nint(getpar('nisomx'))	!maximum number of isovalues allowed
	nisodf = nint(getpar('nisodf'))	!default number of isovalues to use
	icolor = nint(getpar('icolor'))	!color table

	bfunc = rfiso .ne. 0
	blog = ipllog .ne. 0
	iusear = 0			!use array for lookup
	if( bfunc .or. blog ) iusear = 1

	if( bdebug ) then
	  write(6,*) 'colmin/max: ',colmin,colmax
	  write(6,*) 'dval,valmin/max: ',dval,valmin,valmax
	  write(6,*) 'niso,nisomx,nisodf: ',niso,nisomx,nisodf
	end if

! -------------------------------------------------
!  set up data structures
! -------------------------------------------------

	if( valmin > valmax ) then	!no good values
	  valmin = 0.
	  valmax = 0.
	end if

	if( niso .gt. 0 ) then

! 	  ----------------------------------------
! 	  we already know how many isolines to draw (and color table)
! 	  ----------------------------------------

	  dval = 0.
	  if( bfunc .or. blog ) then
	    call mkval(-niso,fiso,valmin,valmax,rfiso,ipllog)
	  else
	    call mkval(niso,fiso,valmin,valmax,dval,ipllog)
	  end if

	else

! 	  ----------------------------------------
! 	  we must compute everything -> niso and dval first
! 	  ----------------------------------------

	  if( dval .ne. 0. ) then
	    niso = nint((valmax-valmin)/dval) + 1
	  else			!find a good dval
	    niso = nisodf
	    dval = (valmax-valmin)/niso
	    if( dval .ne. 0. ) then
	      dval = rnext(dval,+3)	!round dval upwards
	    end if
	    if( dval .ne. 0. ) then	!see if really rounding
	      niso = nint((valmax-valmin)/dval) + 1
	    else
	      dval = 1.
	      niso = 1
	    end if
	  end if

! 	  ----------------------------------------
! 	  be sure values computed make sense
! 	  ----------------------------------------

	  ddval = max(abs(valmax),abs(valmin))
	  eps = 0.
	  if( ddval > 0. ) eps = dval/ddval
	  if( eps .lt. 1.e-4 ) then
	    dval = 1.
	    niso = 1
	  end if

! 	  ----------------------------------------
! 	  check niso if in range
! 	  ----------------------------------------

	  if( niso .le. 0 ) then
	    write(6,*) 'Error computing niso.'
	    write(6,*) dval,valmin,valmax,niso
	    write(6,*) 'If dval < 0 then valmax < valmin.'
	    stop 'error stop colauto: niso'
	  else if( niso .gt. nisomx ) then
	    write(6,*) 'Value of niso > nisomx.'
	    write(6,*) 'niso,nisomx : ',niso,nisomx
	    write(6,*) 'Either increase dval or increase nisomx.'
	    stop 'error stop colauto: nisomx'
	  end if

! 	  ----------------------------------------
! 	  now we have a good niso and dval -> get amin
! 	  ----------------------------------------

	  amin = valmin / dval
	  amin = rnexta(amin,+4)
	  amin = amin * dval
          do while( amin .gt. valmin )          !find starting point
            amin = amin - dval
          end do
          if( amin + niso*dval .lt. amax ) then
            amin = amin + dval
          end if
          amax = amin + (niso-1)*dval   !this is maximum isoline
	  amax = valmax

! 	  ----------------------------------------
! 	  special case - constant array
! 	  ----------------------------------------

	  if( valmin .eq. valmax ) then
	    niso = 2
	    dval = abs(valmax) / 5.
	    if( dval .eq. 0 ) dval = 1.
	    amin = valmax - 0.5*dval
	    amax = amin + dval
	  end if

! 	  ----------------------------------------
! 	  adjust decimal point
! 	  ----------------------------------------

          ndec = getpar('ndccol')
          fact = getpar('faccol')
          rfact = getpar('rfaccol')
	  if( fact == 1. .and. rfact /= 1. ) fact = 1. / rfact
	  
          val = fact * dval
          fdec = log10(val)
          if( fdec .lt. 0. ) then
            idec = 1 + ifix(-fdec)      !1 for 0.1-0.9
            call putpar('ndccol',float(idec))
            !write(6,*) 'decimal point adjusted: ',fact,dval,ndec,idec
          end if

! 	  ----------------------------------------
! 	  set up arrays
! 	  ----------------------------------------

	  if( bfunc .or. blog ) then
	    call mkval(-niso,fiso,valmin,valmax,rfiso,ipllog)
	  else
	    call mkval(niso,fiso,valmin,valmax,dval,ipllog)
	  end if
	  call mkval(niso+1,ciso,colmin,colmax,0.,0)

	end if

	isoanz = niso

	if( bdebug ) then
	  write(6,*) 'colauto final...'
	  write(6,*) 'colmin/max: ',colmin,colmax
	  write(6,*) 'dval,valmin/max: ',dval,valmin,valmax
	  write(6,*) 'niso,nisomx,nisodf: ',niso,nisomx,nisodf
	  write(6,*) 'amin,amax: ',amin,amax
	  write(6,*) 'colors: ',niso+1
	  write(6,*) (ciso(i),i=1,niso+1)
	  write(6,*) 'isovals: ',niso
	  write(6,*) (fiso(i),i=1,niso)
	  write(6,*) 'isoanz: ',isoanz
	end if

	end

! **********************************************************

	subroutine mkval(nn,array,amin,amax,da,ilog)

!  makes array of values
! 
!  if nn is negative in da is power of function to use
!  if nn is negative and ilog is not zero - use logarithmic scale

	implicit none

	integer nn
	real array(abs(nn))
	real amin,amax
	real da
	integer ilog		!logarithmic scale

	logical bfunc,bdebug
	integer i,n
	real dval,r,rn,val

	bdebug = .true.
	bdebug = .false.

	bfunc = nn .lt. 0
	n = abs(nn)
	rn = n
	  
	if( n .le. 0 ) then
	  return
	else if( n .eq. 1 ) then
	  dval = 0.
	else
	  dval = da
	  if( dval .eq. 0. ) dval = (amax-amin) / (n-1)
	end if

	if( bfunc ) then			!non-linear division
	  if( ilog .ne. 0 ) then		!logarithmic division
	    write(6,*) 'using logarithmic scale... '
	    if( amin .le. 0 ) then
	      write(6,*) 'amin,amax: ',amin,amax
	      write(6,*) 'need amin > 0 for logarithmic scale'
	      stop 'error stop mkval: amin <= 0'
	    end if
	    do i=1,n
              r = float(i-1)/(n-1)
              val = amin * (amax/amin)**r
	      array(i) = val
	    end do
	  else if( da .ne. 0. ) then			!polynomial division
	    write(6,*) 'using polynomial scale... ',da
	    dval = amax-amin
	    do i=1,n
	      if( da .gt. 0 ) then
	        array(i) = amin+dval*((i-1)/(n-1.))**da
	      else
	        array(i) = amin+dval*((i/rn)**da)/((1./rn)**da)
	      end if
	    end do
	  else
	    write(6,*) 'nn,da,ilog: ',nn,da,ilog
	    stop 'error stop mkval: error in parameters'
	  end if
	else					!linear division
	  do i=1,n
	    array(i) = amin + (i-1) * dval
	  end do
	end if

	if( bdebug ) then
	  write(6,*) 'mkval: ',n,da
	  write(6,*) (array(i),i=1,n)
	end if

	end

! **********************************************************

	subroutine colset_reg(nmin,nmax)

	use color

	implicit none

	integer nmin,nmax

	logical bldebug
	integer narea,ncol,niso,i
	real colmin,colmax,valmin,valmax,dval

	bldebug = .true.
	bldebug = .false.

	narea = nmax - nmin + 1
	ncol = narea
	niso = narea - 1

	colmin = 0.05
	colmax = 0.95
	dval = 0.
	valmin = nmin + 0.5
	valmax = nmax - 0.5

	call mkval(narea,ciso,colmin,colmax,0.,0)
	call mkval(niso,fiso,valmin,valmax,dval,0)

	if( .not. bldebug ) return

	write(6,*) 'colset_reg: ',narea
	write(6,*) 'colors: ',ncol
	write(6,*) (ciso(i),i=1,ncol)
	write(6,*) 'isolines: ',niso
	write(6,*) (fiso(i),i=1,niso)

	end

! **********************************************************
! **********************************************************
! **********************************************************

	subroutine read_color_table(cname,imap,berr)

	use color

	implicit none

	character*(*) cname
	integer imap
	logical berr

	logical bdebug
	integer icolor,i
	character*80 cfile

	bdebug = .false.
	cfile = colfil

	call read_colormap(cfile,cname,imap,coldim,coltab,berr)
	call set_max_color(coldim,coltab,icmax)

	if( bdebug ) then
	  write(6,*) 'coldim: ',coldim,icmax
	  do i=1,coldim
	    write(6,*) i,coltab(:,i)
	  end do
	end if

	if( berr ) then
	  !stop 'error stop admin_color_table: reading colormap'
	end if

	end

! **********************************************************

	subroutine admin_color_table

	use color

	implicit none

	logical berr,bdebug
	integer icolor,i,imap
	character*80 cfile,cname

	bdebug = .false.
	icmax = 0
	imap = 0

	call getfnm('coltab',cname)
	call getfnm('colfil',cfile)

	!write(6,*) 'colfil = ',trim(cfile)
	!write(6,*) 'coltab = ',trim(cname)

	if( cname == ' ' .or. cfile == ' ' ) return

	write(6,*) 'initializing color table ',trim(cname) &
     &		,' from file ',trim(cfile)

	call read_colormap(cfile,cname,imap,coldim,coltab,berr)
	if( berr ) then
	  write(6,*) 'colfil: ',trim(cfile)
	  write(6,*) 'coltab: ',trim(cname)
	  stop 'error stop admin_color_table: reading colormap'
	end if
	call set_max_color(coldim,coltab,icmax)

	if( bdebug ) then
	  write(6,*) 'coldim: ',coldim,icmax
	  do i=1,coldim
	    write(6,*) i,coltab(:,i)
	  end do
	end if

	write(6,*) 'using colormap: ',icmax,'  ',trim(cname)

	icolor = 8
	call set_color_table(icolor)
	call set_default_color_table(icolor)
	call write_color_table

	end

! **********************************************************

	subroutine set_max_color(coldim,coltab,icmax)

	implicit none

	integer coldim,icmax
	real coltab(3,coldim)

	integer ic

	do ic=1,coldim
	  if( coltab(1,ic) < 0. ) exit
	end do

	icmax = ic - 1

	end

! **********************************************************

	function has_color_table()

	use color

	implicit none

	logical has_color_table

	has_color_table = icmax > 0

	end

! **********************************************************

!******************************************************************

	subroutine read_colormap(file,cname,imap,ncdim,cdata,berr)

	implicit none

	character*(*) file
	character*(*) cname	!check if empty, read if set
	integer imap
	integer ncdim
	real cdata(3,ncdim)
	logical berr

	logical bend,bread,bfound,bwrite
	logical bincolormap
	integer iend,nc,nmap,ios
	real col(3)
	character*80 line,name,aux

	nmap = 0
	nc = 0
	bincolormap = .false.		!inside colormap
	bwrite = .false.		!write info to terminal?
	bread = .false.			!have to read data?
	bfound = .false.		!have found right colormap?
	berr = .true.			!have encountered error?
	iend = 0
	bend = .false.

	if( cname == ' ' .and. imap == 0 ) bwrite = .true.
	if( cname /= ' ' .and. imap /= 0 ) then
	  write(6,*) 'cannot give cname and imap: ',trim(cname),imap
	  berr = .true.
	end if

	cdata = -1.

	open(1,file=file,status='old',form='formatted',iostat=ios)
	if( ios /= 0 ) then
	  write(6,*) 'cannot read colormap ',trim(file)
	  berr = .true.
	  return
	end if

	do 
	  read(1,'(a)',iostat=ios) line
	  if( ios /= 0 ) exit
	  if( line == ' ' ) cycle
	  line = adjustl(line)
	  if( line(1:1) == '#' ) cycle
	  if( line(1:1) == '_' ) then
	    bincolormap = .true.
	    nmap = nmap + 1
	    nc = 0
	    iend = scan(line(2:),'_')
	    name = line(2:iend)
	    if( name == cname .or. imap == nmap ) then
	      bread = .true.
	      bfound = .true.
	      if( cname == ' ' ) cname = name
	      if( imap == 0 ) imap = nmap
	    end if
	    iend = scan(line,'[') + 1
	    line = adjustl(line(iend:))
	    !write(6,*) 'start of colormap: ',trim(name),'  ',trim(line)
	    !if( bread ) write(6,*) '...inserting'
	  end if
	  if( bincolormap ) then
	    bend = .false.
	    iend = len_trim(line)
	    if( line(iend-1:iend) == ']]' ) then
	      bend = .true.
	    else if( line(iend-1:iend) == '],' ) then
	      !still in colormap
	    else
	      write(6,*) 'error reading colormap ',trim(line)
	      exit
	    end if
	  end if
	  aux = line(2:iend-2)
          line = aux
	  !write(6,*) 'scanning line: ',trim(line)
	  read(line,*) col
	  !write(6,*) 'scanned line: ',col
	  nc = nc + 1
	  if( bread ) then
	    if( nc > ncdim ) then
	      write(6,*) 'colortable too big: ',ncdim
	      berr = .true.
	      return
	    end if
	    cdata(:,nc) = col
	  end if

	  if( bend ) then
	    bincolormap = .false.
	    bread = .false.
	    if( bwrite ) then
	      write(6,*) 'colormap: ',nmap,nc,'  ',trim(name)
	    end if
	  end if
	end do

	if( bincolormap ) then
	  write(6,*) 'no end of colormap ',trim(file)
	  berr = .true.
	  return
	end if

	if( ios > 0 ) then
	  write(6,*) 'error reading colormap ',trim(file)
	  berr = .true.
	  return
	end if

	close(1)

	berr = .not. bfound .and. .not. bwrite

	end

!******************************************************************

	subroutine color_table_file_init(file)

	use color

	implicit none

	character*(*) file

	character*80 cfile

	cfile = file
	if( cfile == ' ' ) then
	  call get_standard_color_table(cfile)
	end if
	colfil = cfile

	write(6,*) 'using color table file: ',trim(cfile)

	end

!******************************************************************

	subroutine get_standard_color_table(file)

	implicit none

	character*(*) file

	character*80 cfile,dir
	logical filex

	cfile = 'colormap.dat'

	file = cfile
	if( filex(file) ) return

	call get_environment_variable('SHYFEMDIR',dir)
	if( dir /= ' ' ) then
	  file = trim(dir) // '/femplot/color/' // trim(cfile)
	  if( filex(file) ) return
	end if

	call get_environment_variable('HOME',dir)
	if( dir /= ' ' ) then
	  file = trim(dir) // '/shyfem/femplot/color/' // trim(cfile)
	  if( filex(file) ) return
	end if

	file = ' '

	end

!******************************************************************

