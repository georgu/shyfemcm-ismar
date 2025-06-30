
!--------------------------------------------------------------------------
!
!    Copyright (C) 2003-2004,2009,2014-2015,2017-2019  Georg Umgiesser
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

! revision log :
!
! 20.08.2003	ggu	new laplacian interpolation
! 02.09.2003	ggu	some comments, write to .dat file
! 30.10.2003	ggu	subroutine prepare_bc_l included in this file
! 04.03.2004	ggu	writes also number of variables (1)
! 11.03.2009	ggu	bug fix -> declare hev() here
! 28.01.2014	ggu	changed VERS_6_1_71
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 07.02.2015	ggu	OI finished
! 26.02.2015	ggu	changed VERS_7_1_5
! 17.07.2015	ggu	changed VERS_7_1_53
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 14.09.2015	ggu	OI adapted to new modular structure
! 29.09.2015	ggu	move some vars to custimize to the beginning of file
! 26.09.2017	ggu	changed VERS_7_5_32
! 26.03.2018	ggu	completely restructured for model driven obs
! 03.04.2018	ggu	changed VERS_7_5_43
! 16.02.2019	ggu	changed VERS_7_5_60
! 20.02.2025	ggu	completely restructured
! 18.03.2025	ggu	write new complementary files
! 12.06.2025	ggu	minor changes
!
! notes :
!
! format of observation file content:
!
! header_record_1
! data_record_1
! header_record_2
! data_record_2
! etc...
!
! where header_record is:
!
! time nobs
!
! and data record is:
!
! 1     x_1     y_1     val_1     [error_val_1     [std_val_1]]
! 2     x_2     y_2     val_2     [error_val_2     [std_val_2]]
! ...
! n     x_n     y_n     val_n     [error_val_n     [std_val_n]]
! ...
! nobs  x_nobs  y_nobs  val_nobs  [error_val_nobs  [std_val_nobs]]
!
! legend:
!
! time		time in seconds (relative time) or date (yyyy-mm-dd[::hh:MM:ss])
! nobs		number of observations in data record
! x,y		coordinates of observation
! val		value of observation
! error_val	standard deviation of observation errors
! std_val	standard deviation of observation influence
!
!****************************************************************

!================================================================
	module mod_optintp
!================================================================

	logical, save :: bsilent
	logical, save :: bquiet
	logical, save :: bgeo
	logical, save :: blimit
	logical, save :: bback = .false.
	logical, save :: bfile = .false.

	real, parameter :: flag = -999

	real, save :: rl = flag
	real, save :: drl = flag
	real, save :: rlmax = flag
	real, save :: rr = flag
	real, save :: ss = flag
	real, save :: dx = 0.
	real, save :: dy = 0.
	real, save :: xmin = 0.
	real, save :: xmax = 0.
	real, save :: ymin = 0.
	real, save :: ymax = 0.
	real, save :: tmin = 0.
	real, save :: tmax = 0.
	real, save :: backvalue = flag
	double precision, save :: atime0 = 0.

	character*80, save :: sdxy = ' '
	character*80, save :: sminmax = ' '
	character*80, save :: stau = ' '
	character*80, save :: varstring = 'unknown'
	character*80, save :: date0 = ' '
	character*80, save :: backfile = ' '

!================================================================
	end module mod_optintp
!================================================================

        program optintp

! optimal interpolation interpolation

	use mod_optintp
	use mod_depth
	use evgeom
	use basin
	use clo

	implicit none

	integer nobdim
	parameter (nobdim = 4)

	real xobs(nobdim)
	real yobs(nobdim)
	real zobs(nobdim)
	real bobs(nobdim)
	real errors(nobdim)
	real rra(nobdim)
	real ssa(nobdim)
	real rla(nobdim)
	real rlmaxa(nobdim)
	real rlact(nobdim)
	integer iuse(nobdim)
	integer iuse_orig(nobdim)

	real, allocatable :: xback(:)
	real, allocatable :: yback(:)
	real, allocatable :: zback(:)
	real, allocatable :: zanal(:)
	real, allocatable :: zweight(:)
	real, allocatable :: ztau(:)
	real, allocatable :: zweight_old(:)
	real, allocatable :: ztau_old(:)
	real, allocatable :: zaux(:)
	real, allocatable :: hd(:)
	integer, allocatable :: ilhkv(:)

	integer np
	integer iufem,iuobs,iuweight,iutau,iu
	integer jmax,j
	double precision dtime,atime,atime_old
	real hlv(1)
	real regpar(7)
	integer date,time
	character*80 string,format
	character*20 aline

	character*80 file,obsfile
	logical bcart,breg
	logical bfem,bloop,bweight
	logical bdebug
	logical bchanged
	integer k,ie,n,ndim
	integer nobs,nback,nobs_valid
        integer ilev,ivar,irec,nrec
	integer index
	integer iexcl
	integer nx,ny
	real dxy
	real x0,y0,x1,y1
	real zmin,zmax
	real zomin,zomax
	real rmse,rmse1,rmse2

	logical is_spherical
	logical dts_is_atime

	bdebug = .false.
	bfem = .false.
	bfem = .true.
	bloop = .true.
	bloop = .false.

!--------------------------------------------------------------
! parameters and command line options
!--------------------------------------------------------------

	call optintp_init(obsfile)
	call success_init

!-----------------------------------------------------------------
! set some variables
!-----------------------------------------------------------------

	!ivar = 85		!what is in the observations - ice
	!string = 'ice cover [0-1]'
	!date = 19970101

	!ivar = 12		!what is in the observations - temperature
	!string = 'temperature [C]'
	!date = 20100101

	!ivar = 11		!what is in the observations - temperature
	!string = 'salinity [psu]'
	!date = 20100101

	!ivar = 1		!what is in the observations - water level
	!string = 'water level [m]'
	!date = 20130101

	nrec = 10
	nrec = 5
	nrec = 0		!do not read more than this number of records

!-----------------------------------------------------------------
! read in basin
!-----------------------------------------------------------------

	dxy = max((xmax-xmin),(ymax-ymin))

	if( rl == flag ) rl = 0.1*dxy
	if( rlmax == flag ) rlmax = 3.*rl
	if( rr == flag ) rr = 0.01
	if( ss == flag ) ss = 100.*rr

	rra = rr
	ssa = ss
	rla = rl
	rlmaxa = rlmax
	bweight = ( 0 < tmin .and. tmin <= tmax )

	if( .not. bquiet ) then
	  write(6,*) 'parameters used:'
	  write(6,*) 'variable:  ',trim(varstring)
	  write(6,*) 'xmin/xmax: ',xmin,xmax
	  write(6,*) 'ymin/ymax: ',ymin,ymax
	  write(6,*) 'dxy      : ',dxy
	  write(6,*) 'rl:        ',rl
	  write(6,*) 'rlmax:     ',rlmax
	  write(6,*) 'drl:       ',drl
	  write(6,*) 'rr:        ',rr
	  write(6,*) 'ss:        ',ss
	  write(6,*) 'dx,dy:     ',dx,dy
	  write(6,*) 'tmin,tmax: ',tmin,tmax
	  write(6,*) 'limit:     ',blimit
	  write(6,*) 'weight:    ',bweight
	end if

!-----------------------------------------------------------------
! set up ev and background grid
!-----------------------------------------------------------------

	regpar = 0.
	nback = 0
	call setup_background_size(dx,dy,xmin,ymin,xmax,ymax &
     &			,nback,regpar)

	allocate(xback(nback),yback(nback),zback(nback),zanal(nback))
	allocate(zweight(nback),ztau(nback))
	allocate(zweight_old(nback),ztau_old(nback))
	allocate(zaux(nback))

	zback = flag
	zanal = flag

	call setup_background_coords(nback,regpar,xback,yback)

	call read_background_values(backfile,nback,zback)

	if( bback ) then
	  zback = backvalue
	  bobs = backvalue
	end if

	nx = nint(regpar(1))
	ny = nint(regpar(2))
	x0 = regpar(3)
	y0 = regpar(4)
	x1 = x0 + nx*dx
	y1 = y0 + ny*dy

	breg = dx > 0.
	if( .not. breg ) bloop = .false.

	if( .not. bquiet ) then
	  write(6,*) 'background grid:'
	  write(6,*) 'nback: ',nback
	  if( breg ) then
	    write(6,*) 'nx,ny: ',nx,ny
	    write(6,*) 'dx,dy: ',dx,dy
	    write(6,*) 'x0,y0: ',x0,y0
	    write(6,*) 'x1,y1: ',x1,y1
	  end if
	end if

	format = 'formatted'
	np = nback

!-----------------------------------------------------------------
! write single fem files
!-----------------------------------------------------------------

	!if( bfem ) then
	!if( .false. ) then
	if( .true. ) then
	  iu = 3
	  atime = 0.

	  string = varstring
	  file = 'back.fem'
	  zaux = zback
	  open(iu,file=file,status='unknown',form='formatted')
	  call write_fem_record(iu,atime,regpar,string,np,zaux)
	  close(iu)

	  string = 'weight'
	  file = 'weight.fem'
	  zaux = 0.
	  open(iu,file=file,status='unknown',form='formatted')
	  call write_fem_record(iu,atime,regpar,string,np,zaux)
	  close(iu)

	  string = 'time scale tau [s]'
	  file = 'tau.fem'
	  zaux = tmax
	  open(iu,file=file,status='unknown',form='formatted')
	  call write_fem_record(iu,atime,regpar,string,np,zaux)
	  close(iu)
	end if

!-----------------------------------------------------------------
! open files
!-----------------------------------------------------------------

	irec = 0
	!drl = 0.05		!test different rl, distance drl
	jmax = 0
	if( drl > 0. ) jmax = 5
	if( jmax > 0 .and. .not. bquiet ) then
	  write(6,*) 'trying multiple values for rl: ',rl,drl,jmax
	end if

	iuobs = 20
	open(iuobs,file=obsfile,status='old',form='formatted')
	if( .not. bquiet ) write(6,*) 'file with observations: ',trim(obsfile)

	iufem = 21
	open(iufem,file='optintp.fem',status='unknown',form=format)

	if( bweight ) then
	  iuweight = 22
	  open(iuweight,file='optweight.fem',status='unknown',form=format)
	  iutau = 23
	  open(iutau,file='opttau.fem',status='unknown',form=format)
	end if

	iuse = 1
	nobs = nobdim
	call read_use_file(nobs,iuse_orig)	!what observations to use

!-----------------------------------------------------------------
! read observations and interpolate
!-----------------------------------------------------------------

	if( bloop .and. .not. bquiet ) then
	  write(6,*) '                time     ' // &
     &		'min/max interpol  min/max observed   rmse'
	end if

	iexcl = 0

!-----------------------------------------------------------------
! loop over observations
!-----------------------------------------------------------------

	do

	if ( nrec > 0 .and. irec == nrec ) exit		!maximum records reached

	nobs = nobdim
	rra = rr
	rla = rl
	call read_observations(iuobs &
     &				,atime,nobs,xobs,yobs,zobs,rra,rla)
	if( nobs <= 0 ) exit
	irec = irec + 1
	where( iuse_orig == 0 ) zobs = flag
	where( zobs == flag ) iuse = 0
	nobs_valid = count( iuse(1:nobs) == 1 )

	call dts_format_abs_time(atime,aline)
	if( .not. bsilent ) then
	  write(6,*) 'observation read: ',trim(aline),irec,nobs,nobs_valid
	end if

	!call mima(zobs,nobs,zomin,zomax)
	call mima_flag(nobs,zobs,flag,zomin,zomax)
	!write(6,*) 'observations min/max: ',zomin,zomax

	if( irec == 1 ) then
	  call write_obs_grid(nobs,xobs,yobs,zobs)
	end if

!-----------------------------------------------------------------
! interpolate
!-----------------------------------------------------------------

	do j=-jmax,jmax			!this loops over rl if jmax>0
	  rlact = rla * ( 1. + j*drl )
          call opt_intp(nobs,xobs,yobs,zobs,bobs &
     &                  ,nback,bback,xback,yback,zback &
     &                  ,rlact,rlmaxa,ssa,rra,zanal)

	  call mima(zanal,nback,zmin,zmax)
	  if( blimit ) call limit_values(nback,zanal,zomin,zomax)

	  atime = atime + j		!just to distinguish between time steps

	  if( bfem ) then
	    string = varstring
	    call write_fem_record(iufem,atime,regpar,string,np,zanal)
	  end if
	  call obs_have_changed(nobs,xobs,yobs,iuse,bchanged)
	  !if( bchanged ) write(99,*) aline,irec,bchanged
	  if( bweight .and. irec > 1 .and. bchanged ) then
	    if( .not. bquiet ) write(6,*) 'observations have changed: ',aline
	    string = 'weight'
	    call write_fem_record(iuweight,atime_old,regpar,string,np &
     &						,zweight_old)
	    string = 'time scale tau [s]'
	    call write_fem_record(iutau,atime_old,regpar,string,np,ztau_old)
	  end if
	  if( bweight .and. bchanged ) then
	    call make_weight(nobs,xobs,yobs,zobs,rlact,nback,xback,yback &
     &						,zweight)
	    string = 'weight'
	    call write_fem_record(iuweight,atime,regpar,string,np,zweight)
	    call make_tau(tmin,tmax,nback,zweight,ztau)
	    string = 'time scale tau [s]'
	    call write_fem_record(iutau,atime,regpar,string,np,ztau)
	  end if
	  if( bloop ) then
	    call compute_rmse(nobs,xobs,yobs,zobs,nx,ny,regpar &
     &				,zanal,0,rmse1)
	    write(67,*) atime,nobs
            call opt_intp_loop(nobs,xobs,yobs,zobs,bobs &
     &                  ,nback,bback,xback,yback,zback,regpar &
     &                  ,rlact,rlmaxa,ssa,rra,zanal,rmse2,errors)
	    !write(6,*) rmse1,rmse2
	    write(66,*) atime,nobs
	    write(66,*) errors(1:nobs)
	  end if

	  atime_old = atime
	  zweight_old = zweight
	  ztau_old = ztau

	  call dts_format_abs_time(atime,aline)
	  if( .not. bquiet ) then
	    if( bfem .and. .not. bloop ) then	!multiple records in time
	      write(6,3000) 'min/max: ',aline,zmin,zmax,zomin,zomax
 1000	      format(a,a,4f12.2)
	    else if( bloop ) then
	      write(6,3000) 'min/max: ',aline &
     &				,zmin,zmax,zomin,zomax,rmse1,rmse2
 3000	      format(a,a,4f9.2,2x,2f9.3)
	  end if
	end if

	end do

!-----------------------------------------------------------------
! end of loop over observations
!-----------------------------------------------------------------

	end do

!-----------------------------------------------------------------
! final weights and tau
!-----------------------------------------------------------------

	if( .not. bquiet ) then
	  write(6,*) 'final observations have changed: ',aline
	end if

	string = 'weight'
	call write_fem_record(iuweight,atime_old,regpar,string,np &
     &						,zweight_old)
	string = 'time scale tau [s]'
	call write_fem_record(iutau,atime_old,regpar,string,np,ztau_old)

!-----------------------------------------------------------------
! final message
!-----------------------------------------------------------------

	if( .not. bsilent ) then
	  write(6,*) 'total number of records treated: ',irec
	  if( irec == nrec ) then
	    write(6,*) 'limiting records treated to: ',nrec
	  end if
	end if

	call success_final

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

	end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine setup_background_size(dx,dy,xmin,ymin,xmax,ymax &
     &			,nback,regpar)

! computes the limits and the size of the background grid

	use mod_optintp, only : bquiet

	implicit none

	real, intent(in) ::	 dx,dy
	real, intent(in) ::	 xmin,ymin,xmax,ymax
	integer, intent(out) ::	 nback			!size of background grid
	real, intent(out) ::	 regpar(7)		!values of grid

	integer k,i,j
	integer nx,ny
	real x0,y0,x1,y1
	real xdiff,ydiff,flag

	call create_reg(dx,dy,xmin,ymin,xmax,ymax,regpar)
	nx = nint(regpar(1))
	ny = nint(regpar(2))
	nback = nx*ny

	if( nx > 0. .and. .not. bquiet ) then
	  call printreg(regpar)
	end if

	return
   95	continue
	write(6,*) x0,xmin,xmax,x1
	write(6,*) y0,ymin,ymax,y1
	stop 'error stop setup_background: internal error (3)'
   96	continue
	write(6,*) 'dx,dy: ',dx,dy
	write(6,*) 'xmin,xmax: ',xmin,xmax
	write(6,*) 'ymin,ymax: ',ymin,ymax
	stop 'error stop setup_background: error in parameters'
	end

!****************************************************************

	subroutine setup_background_coords(nback,regpar &
     &			,xback,yback)

	implicit none

	integer nback
	real regpar(7)
	real xback(nback)
	real yback(nback)

	integer n,i,j
	integer nx,ny
	real x0,y0,x1,y1,dx,dy

	nx = regpar(1)
	ny = regpar(2)
	x0 = regpar(3)
	y0 = regpar(4)
	dx = regpar(5)
	dy = regpar(6)

	x1 = x0 + nx*dx
	y1 = y0 + ny*dy

	n = 0
	do j=0,ny-1
	  do i=0,nx-1
	    n = n + 1
	    if( n > nback ) goto 99
	    xback(n) = x0 + i*dx
	    yback(n) = y0 + j*dy
	  end do
	end do

	if( nback .ne. n ) goto 99

	return
   99	continue
	write(6,*) nx,ny,nx*ny,nback,n
	stop 'error stop setup_background_coords: internal error (1)'
	end

!****************************************************************

	subroutine read_background_size(file,nback,regpar)

! reads grid and returns characteristics

	implicit none

	character*(*) file
	integer nback
	real regpar(7)

	integer iformat,iunit
	integer nvers,np,lmax,nvar,ntype,ierr
	integer datetime(2)
	real hlv(1)
	double precision dtime

	call fem_file_read_open(file,0,iformat,iunit)

        call fem_file_read_params(iformat,iunit,dtime,nvers,np,lmax    &
                           &              ,nvar,ntype,datetime,ierr)

	if( ierr /= 0 ) goto 99
	if( lmax > 1 ) goto 99
	if( nvar > 1 ) goto 99

	call fem_file_read_2header(iformat,iunit,ntype,lmax,hlv,regpar,ierr)

	if( ierr /= 0 ) goto 99

	close(iunit)
	nback = np

	return
   99	continue
	write(6,*) 'errors in background file'
	write(6,*) 'ierr,lmax,nvar: ',ierr,lmax,nvar
	stop 'error stop read_background_size: parameter errors'
	end

!****************************************************************

	subroutine read_background_values(file,nback,zback)

! reads grid and returns characteristics

	implicit none

	character*(*) file
	integer nback
	real zback(nback)

	integer iformat,iunit
	integer nlvddi
	integer nvers,np,lmax,nvar,ntype,ierr
	integer datetime(2)
	real regpar(7)
	real hlv(1)
	real hd(1)
	integer ilhkv(1)
	double precision dtime
	character*80 string

	if( file == ' ' ) return

	nlvddi = 1

	call fem_file_read_open(file,0,iformat,iunit)

        call fem_file_read_params(iformat,iunit,dtime,nvers,np,lmax    &
                           &              ,nvar,ntype,datetime,ierr)

	if( ierr /= 0 ) goto 99
	if( lmax > 1 ) goto 99
	if( nvar > 1 ) goto 99
	if( np /= nback ) goto 99

	call fem_file_read_2header(iformat,iunit,ntype,lmax,hlv,regpar,ierr)

	if( ierr /= 0 ) goto 99

        call fem_file_read_data(iformat,iunit,nvers,np,lmax   &
                  &        ,string,ilhkv,hd,nlvddi,zback,ierr)

	close(iunit)

	return
   99	continue
	write(6,*) 'errors in background file'
	write(6,*) 'ierr,lmax,nvar: ',ierr,lmax,nvar
	write(6,*) 'np,nback: ',np,nback
	stop 'error stop read_background_values: parameter errors'
	end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine limit_values(nback,zback,zomin,zomax)

! limits computed values to min/max of observations

	implicit none

	integer nback
	real zback(nback)
	real zomin,zomax

	integer i

	do i=1,nback
	  zback(i) = max(zback(i),zomin)
	  zback(i) = min(zback(i),zomax)
	end do

	end

!****************************************************************

	subroutine make_weight(nobs,xobs,yobs,zobs,rl,nback,xback,yback &
     &					,zweight)

	use mod_optintp, only : flag

	implicit none

	integer nobs
	real xobs(nobs)
	real yobs(nobs)
	real zobs(nobs)
	real rl(nobs)
	integer nback
	real xback(nback)
	real yback(nback)
	real zweight(nback)

	integer j,i
	real xi,yi,zi,xj,yj
	real dist2,z
	real rl2(nobs)

	rl2 = rl**2

        do j=1,nback
          xj = xback(j)
          yj = yback(j)
	  z = 0
          do i=1,nobs
            xi = xobs(i)
            yi = yobs(i)
            zi = zobs(i)
	    if( zi == flag ) cycle
            dist2 = (xi-xj)**2 + (yi-yj)**2
            z = z + exp( -dist2/rl2(i) )
          end do
	  z = min(1.,z)
          zweight(j) = z
        end do

	end

!****************************************************************

	subroutine make_tau(tmin,tmax,nback,zweight,ztau)

	implicit none

	real tmin,tmax
	integer nback
	real zweight(nback)
	real ztau(nback)

	integer j
	real z,t

        do j=1,nback
	  z = zweight(j)
	  t = tmin + (1.-z)*(tmax-tmin)
	  ztau(j) = t
	end do

	end

!******************************************************************

	subroutine write_reg_grid(nx,ny,x0,y0,dx,dy)

	implicit none

	integer nx,ny
	real x0,y0
	real dx,dy

	integer iu
	integer ix,iy,i
	real x,y

	iu = 9
	open(iu,file='regfile.grd',form='formatted',status='unknown')

	i = 0
	do iy=0,ny
	  do ix=0,nx
	    i = i + 1
	    x = x0 + ix*dx
	    y = y0 + iy*dy
	    write(iu,*) 1,i,0,x,y
	  end do
	end do

	close(iu)

	end

!******************************************************************

	subroutine write_box_grid(nx,ny,x0,y0,dx,dy)

	implicit none

	integer nx,ny
	real x0,y0
	real dx,dy

	integer iu
	integer ix,iy,in,il
	real x,y
	real x1,x2,y1,y2

	iu = 9
	open(iu,file='boxfile.grd',form='formatted',status='unknown')

	in = 0
	il = 0

	do iy=0,ny
	    il = il + 1
	    x1 = x0
	    x2 = x0 + nx*dx
	    y = y0 + iy*dy
	    write(iu,*) 1,in+1,0,x1,y
	    write(iu,*) 1,in+2,0,x2,y
	    write(iu,*) 3,il,0,2,in+1,in+2
	    in = in + 2
	end do
	do ix=0,nx
	    il = il + 1
	    y1 = y0
	    y2 = y0 + ny*dy
	    x = x0 + ix*dx
	    write(iu,*) 1,in+1,0,x,y1
	    write(iu,*) 1,in+2,0,x,y2
	    write(iu,*) 3,il,0,2,in+1,in+2
	    in = in + 2
	end do

	close(iu)

	end

!******************************************************************

	subroutine write_obs_grid(nobs,xobs,yobs,zobs)

	implicit none

	integer nobs
	real xobs(nobs)
	real yobs(nobs)
	real zobs(nobs)

	integer iu
	integer i
	real x,y,z

	iu = 9
	open(iu,file='obsfile.grd',form='formatted',status='unknown')

	do i=1,nobs
	    x = xobs(i)
	    y = yobs(i)
	    z = zobs(i)
	    write(iu,*) 1,i,3,x,y,z
	end do

	close(iu)

	end

!******************************************************************

	subroutine compute_rmse(nobs,xobs,yobs,zobs,nx,ny,regpar &
     &					,za,iexcl,rmse)

	implicit none

	integer nobs
	real xobs(nobs)
	real yobs(nobs)
	real zobs(nobs)
	integer nx,ny
	real regpar(7)
	real za(nx,ny)
	integer iexcl		!excluded obs - compute error for this
	real rmse

	integer i,ix,iy,ntot
	integer istart,iend
	real dx,dy,x0,y0
	real x,y,z
	real rx,ry,xx,yy,zm
	real zz(4)
	real flag
	double precision rrr

	real rbilin

	x0 = regpar(3)
	y0 = regpar(4)
	dx = regpar(5)
	dy = regpar(6)
	flag = regpar(7)

	rrr = 0.
	istart = 1
	iend = nobs
	ntot = nobs
	if( iexcl > 0 ) then
	  istart = iexcl
	  iend = iexcl
	  ntot = 1
	end if

	do i=istart,iend
	  x = xobs(i)
	  y = yobs(i)
	  z = zobs(i)
	  ix = 1 + (x-x0)/dx
	  if( ix >= nx ) ix = nx - 1
	  xx = x0 + (ix-1)*dx
	  rx = (x - xx)/dx
	  iy = 1 + (y-y0)/dy
	  if( iy >= ny ) iy = ny - 1
	  yy = y0 + (iy-1)*dy
	  ry = (y - yy)/dy
	  zz = (/za(ix,iy),za(ix+1,iy),za(ix,iy+1),za(ix+1,iy+1)/)
	  zm = rbilin(zz,rx,ry,flag)
	!write(6,*) i,ix,iy,rx,ry
	!write(6,*) zz,zm,z
	  if( iexcl > 0 ) write(67,*) i,z,zm
	  rrr = rrr + (zm-z)**2
	end do

	rmse = sqrt( rrr / ntot )

	end

!******************************************************************

        subroutine opt_intp_loop(nobs,xobs,yobs,zobs,bobs &
     &                  ,nback,bback,xback,yback,zback,regpar &
     &                  ,rlact,rlmax,ssa,rra,zanal,rmse,errors)

	implicit none

	integer nobs
	real xobs(nobs)
	real yobs(nobs)
	real zobs(nobs)
	real bobs(nobs)
        integer nback           !size of background field
        logical bback           !background values are given
        real xback(nback)       !x-coordinates of background
        real yback(nback)       !y-coordinates of background
        real zback(nback)       !values of background
	real regpar(7)
        real rlact(nobs)        !length scale for covariance
        real rlmax(nobs)        !max radius to be considered
        real ssa(nobs)          !std of background field
        real rra(nobs)          !std of observation error matrix
        real zanal(nback)       !analysis on return
	real rmse
	real errors(nobs)

	integer i
	integer nx,ny
	real error,rrsave
	double precision rrr

	rrr = 0.
	nx =  nint(regpar(1))
	ny =  nint(regpar(2))

	do i=1,nobs
	  rrsave = rra(i)
	  rra(i) = 1000. * rra(i)
          call opt_intp(nobs,xobs,yobs,zobs,bobs &
     &                  ,nback,bback,xback,yback,zback &
     &                  ,rlact,rlmax,ssa,rra,zanal)
	  call compute_rmse(nobs,xobs,yobs,zobs,nx,ny,regpar &
     &					,zanal,i,error)
	  errors(i) = error
	  rrr = rrr + error**2
	  rra(i) = rrsave
	end do

	rmse = sqrt( rrr / nobs )

	end

!******************************************************************

	subroutine read_use_file(nobs,iuse)

! read use.txt that indicates what observations to use

	implicit none

	integer nobs
	integer iuse(nobs)

	integer ios,i,j,iu,n

	iuse = 1

	iu = 1
	open(iu,file='use.txt',status='old',form='formatted',iostat=ios)
	if( ios /= 0 ) return

	read(iu,*) n
	if( n /= nobs ) then
	  write(6,*) n,nobs
	  stop 'error stop read_use_file: nobs'
	end if

	do i=1,nobs
	  read(iu,*) j,iuse(i)
	  if( i /= j ) stop 'error stop read_use_file: i/=j'
	end do

	close(iu)

	end

!****************************************************************

	subroutine write_fem_record(iu,atime,regpar,string,np,z)

! writes one record of fem file

	use mod_optintp

	implicit none

	integer iu
	double precision atime
	real regpar(7)
	character*(*) string
	integer np
	real z(np)
	
	integer iformat,nvers,lmax,nvar,ntype,nlvdi
	real hlv(1)
	integer ilhkv(1)
	real hd(1)
	double precision dtime
	integer datetime(2)
	integer date,time

	logical dts_is_atime

	iformat = 1
	nvers = 0
	lmax = 1
	nvar = 1
	ntype = 11
	nlvdi = 1
	hlv = 10000.
	ilhkv = 1
	hd = flag

	if( dts_is_atime(atime) ) then
	  call dts_from_abs_time(date,time,atime)
	  datetime = (/date,time/)
	  dtime = 0
	else
	  datetime = 0
	  dtime = atime
	end if

	if( .not. bquiet ) write(6,*) 'writing fem file with ',trim(string)

	call fem_file_write_header(iformat,iu,dtime &
     &                          ,nvers,np,lmax &
     &                          ,nvar,ntype &
     &                          ,nlvdi,hlv,datetime,regpar)
        call fem_file_write_data(iformat,iu &
     &                          ,nvers,np,lmax &
     &                          ,string &
     &                          ,ilhkv,hd &
     &                          ,nlvdi,z)

	end

!****************************************************************

	subroutine obs_have_changed(nobs,x,y,u,bchanged)

	implicit none

	integer nobs
	real x(nobs)
	real y(nobs)
	integer u(nobs)
	logical bchanged

	real, save, allocatable :: xold(:)
	real, save, allocatable :: yold(:)
	integer, save, allocatable :: uold(:)

	if( .not. allocated(uold) ) then
	  allocate(xold(nobs),yold(nobs),uold(nobs))
	  xold = x
	  yold = y
	  uold = u
	  bchanged = .true.
	  return
	end if

	bchanged = .false.

	if( any(x/=xold) ) goto 99
	if( any(y/=yold) ) goto 99
	if( any(u/=uold) ) bchanged = .true.

	xold = x
	yold = y
	uold = u

	return
   99	continue
	write(6,*) 'coordinates have changed:'
	write(6,*) x
	write(6,*) xold
	write(6,*) y
	write(6,*) yold
	stop 'error stop obs_have_changed: x/y have changed'
	end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine parse_time_line(line,atime,nobs)

! parses time header of observations
!
! on return atime is either absolute time (because string given or date0 set)
! or relative time because time is number and date0 is not set

	use mod_optintp

	implicit none

	character*(*) line
	double precision atime
	integer nobs

	integer ns,ierr
	double precision d(10)
	character*80 strings(2)

	logical dts_is_atime
	integer iscans,iscand

	!write(6,*) trim(line)

	ns = iscans(line,strings,2)
	if( ns /= 2 ) goto 99		!we must have two values
	
	call dts_string2time(strings(1),atime,ierr)
	if( ierr /= 0 ) goto 99		!cannot parse time variable
	if( .not. dts_is_atime(atime) ) atime = atime + atime0

	ns = iscand(strings(2),d,1)
	if( ns /= 1 ) goto 99		!no nobs on line
	nobs = nint(d(1))

	return
   99	continue
	write(6,*) 'cannot parse line: ',trim(line)
	stop 'error stop parse_time_line: cannot parse'
	end

!****************************************************************

	subroutine read_observations(iuobs &
     &				,atime,nobs,xobs,yobs,zobs,rra,rla)

! reads one set of observations (one time record)
!
! file format is given in header of this file

	use mod_optintp, only: flag

	implicit none

	integer iuobs			!unit to be read from
	double precision atime		!relative time
	integer nobs			!number of observations
	real xobs(nobs)			!x-coordinate
	real yobs(nobs)			!y-coordinate
	real zobs(nobs)			!value of observation
	real rra(nobs)			!error of observation
	real rla(nobs)			!std of observation

	character*80 line
	logical bdebug
	integer k,kn,ndim,n,ianz,ios,i,j
	integer, save :: nobs_orig = 0
	!double precision dtime,atime
	real f(10)

	integer ipint,iscanf

	bdebug = .false.

	ndim = nobs
	nobs = 0

	if( bdebug ) then
	  write(6,*) '...reading observations from unit :',iuobs
	end if

	!-----------------------------------------------
	! loop over empty lines to find next time header
	!-----------------------------------------------

	do
	  read(iuobs,'(a)',iostat=ios) line
	  if( ios < 0 ) return
	  if( ios > 0 ) goto 95
	  if( line /= ' ' ) exit		!header found
	end do
	call parse_time_line(line,atime,nobs)

	if( nobs > ndim ) goto 96

	if( bdebug ) write(6,*) 'atime = ',atime

	!-----------------------------------------------
	! loop over observations
	!-----------------------------------------------

	xobs = 0.
	yobs = 0.
	zobs = flag

	do i=1,nobs
	  read(iuobs,'(a)',iostat=ios) line
	  if( ios /= 0 ) goto 94
	  ianz = iscanf(line,f,6)
	  if( ianz == 0 ) cycle			!empty line
	  if( ianz < 0 ) goto 92		!error parsing
	  if( ianz < 4 ) goto 91		!too few values
	    j = nint(f(1))
	    if( i /= j ) goto 93		!input is garbled
	    xobs(i) = f(2)
	    yobs(i) = f(3)
	    zobs(i) = f(4)
	    if( ianz >= 5 .and. f(5) > 0. ) rra(i) = f(5)	!optional rr
	    if( ianz >= 6 .and. f(6) > 0. ) rla(i) = f(6)	!optional rl
	end do

	if( nobs_orig == 0 ) nobs_orig = nobs
	if( nobs /= nobs_orig ) goto 90		!number of observations changed

	!-----------------------------------------------
	! end of routine
	!-----------------------------------------------

	return
   90	continue
	write(6,*) 'inconsistent nobs: ',nobs_orig,nobs
	write(6,*) 'number of observations must always be the same'
	stop 'error stop read_observations: nobs_orig /= nobs'
   91	continue
	write(6,*) 'error reading data record: ',ianz
	write(6,*) 'need at least 4 values on line'
	stop 'error stop read_observations: ianz < 4'
   92	continue
	write(6,*) 'error reading data record: ',ianz
	write(6,*) 'line: ',trim(line)
	stop 'error stop read_observations: cannot parse line'
   93	continue
	write(6,*) 'error reading data record: ',i,j
	stop 'error stop read_observations: read error in data record'
   94	continue
	write(6,*) 'error reading data line'
	stop 'error stop read_observations: read error in data line'
   95	continue
	write(6,*) 'error reading header record'
	stop 'error stop read_observations: read error in header line'
   96	continue
	write(6,*) nobs,ndim
	stop 'error stop read_observations: nobs>ndim'
	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine optintp_init(file)

	use clo
	use mod_optintp

	implicit none

	character*(*) file

	integer ns,ierr
	real f(4)

!-------------------------------------------------------------
! initialize and set options
!-------------------------------------------------------------

        call clo_init('optintp','obs-file','1.0')

        call clo_add_info('Optimal interpolation for sparse data')
        call clo_add_extra('  data are interpolated onto a regular field')
        call clo_add_option('quiet',.false.,'do not be verbose')
        call clo_add_option('silent',.false.,'be silent')
        call clo_add_option('geo',.false. &
     &		,'coordinates are spherical')
        call clo_add_option('limit',.false. &
     &		,'limit values to min/max of observations')
        call clo_add_option('rl #',flag &
     &		,'set length scale for covariance matrix')
        call clo_add_option('drl #',flag &
     &		,'try multiple values of rl with percentage step drl')
        call clo_add_option('rlmax #',flag &
     &		,'maximum distance of nodes to be considered')
        call clo_add_option('rr #',flag &
     &		,'std of observation errors')
        call clo_add_option('ss #',flag &
     &		,'std of background field')
        call clo_add_option('dxy dx[,dy]',' ' &
     &		,'dx,dy for regular field')
        call clo_add_option('minmax xmin,ymin,xmax,ymax',' ' &
     &		,'min/max for regular field')
        call clo_add_option('backvalue val',flag &
     &		,'use val as average background value')
        call clo_add_option('backfile file',' ' &
     &		,'use file as background grid')
        call clo_add_option('tau tmin,tmax',' ' &
     &		,'tmin,tmax for time scale')
        call clo_add_option('date0 date',' ' &
     &		,'date for time 0 if relative time')
        call clo_add_option('variable varstring',' ' &
     &		,'name of variable to be written to file')
	call clo_add_extra('defaults for undefined values:')
	call clo_add_extra('  rl=1/10 of basin size, rlmax=10*rl')
	call clo_add_extra('  rr=0.01  ss=100*rr')
	call clo_add_extra('  Default for dx is 0 (use FEM grid)')
	call clo_add_extra('Format for date is yyyy-mm-dd[::hh:MM:ss]')

!-------------------------------------------------------------
! get options
!-------------------------------------------------------------

        call clo_parse_options(1)  !expecting (at least) 1 file after options

        call clo_get_option('silent',bsilent)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('geo',bgeo)
        call clo_get_option('limit',blimit)
        call clo_get_option('rl',rl)
        call clo_get_option('drl',drl)
        call clo_get_option('rlmax',rlmax)
        call clo_get_option('rr',rr)
        call clo_get_option('ss',ss)
        call clo_get_option('dxy',sdxy)
        call clo_get_option('minmax',sminmax)
        call clo_get_option('backvalue',backvalue)
        call clo_get_option('backfile',backfile)
        call clo_get_option('tau',stau)
        call clo_get_option('date0',date0)
        call clo_get_option('variable',varstring)

	!call clo_info()

	call clo_get_file(1,file)

!-------------------------------------------------------------
! elaborate options
!-------------------------------------------------------------

	if( bsilent ) bquiet = .true.

	if( backfile /= ' ' ) then		!background file given
	else
	  call get_options(2,sdxy,ns,f)
	  if( ns < 1 .or. ns > 2 ) then
	    stop 'error stop: no dx,dy given or error'
	  end if
	  dx = f(1)
	  dy = dx
	  if( ns == 2 ) dy = f(2)
	  if( dx*dy == 0 ) stop 'error stop: dx or dy is 0'

	  call get_options(4,sminmax,ns,f)
	  if( ns < 1 .or. ns > 4 ) then
	    stop 'error stop: no minmax given or error'
	  end if
	  xmin = f(1)
	  ymin = f(2)
	  xmax = f(3)
	  ymax = f(4)
	  if( xmin >= xmax .or. ymin >= ymax ) then
	    stop 'error stop: min >= max in minmax'
	  end if
	end if

	if( stau /= ' ' ) then
	  call get_options(2,stau,ns,f)
	  if( ns < 1 .or. ns > 2 ) then
	    stop 'error stop: option tau needs 2 values tmin,tmax'
	  end if
	  tmin = f(1)
	  tmax = f(2)
	  if( .not. bquiet ) write(6,*) 'tmin,tmax: ',tmin,tmax
	end if

	if( date0 /= ' ' ) then
	  call dts_string2time(date0,atime0,ierr)
	  if( ierr /= 0 ) then
	    write(6,*) 'cannot parse date0: ',trim(date0)
	    stop 'error stop: error in date0'
	  end if
	end if

	if( backfile /= ' ' .and. backvalue /= flag ) then
	  write(6,*) 'only one of backfile or backvalue can be given'
	  write(6,*) 'backfile : ',trim(backfile)
	  write(6,*) 'backvalue : ',backvalue
	  stop 'error stop: error in backfile and backvalue'
	end if

	bback = ( backvalue /= flag )
	bfile = ( backfile /= ' ' )

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end subroutine

!******************************************************************

	subroutine success_init

	implicit none

	integer ios
	character*80, save :: file = "success.txt"

	open(unit=1, iostat=ios, file=file, status='old')
	if (ios == 0) close(1, status='delete')

	end subroutine

!******************************************************************

	subroutine success_final

	implicit none

	integer ios
	character*80, save :: file = "success.txt"

	open(1,file=file,status='new',form='formatted')
	write(1,*) 'optintp finished with success'
	close(1)

	end subroutine

!******************************************************************

