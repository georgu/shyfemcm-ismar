
!--------------------------------------------------------------------------
!
!    Copyright (C) 1999,2005,2007,2009-2019  Georg Umgiesser
!    Copyright (C) 2017  Christian Ferrarin
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

! assimilation routines
!
! revision log :
!
! 17.10.2025    ggu     new routine scalar_nudging()
! 21.10.2025    ggu     new routine scalar_assimilation()
! 24.10.2025    ggu     first draft of assimilation ready
! 27.10.2025    ggu     pass 3d array, change of names
! 30.10.2025    ggu     bug fixes
! 07.11.2025    ggu     introduced iuse
!
!****************************************************************

!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine scalar_assimilation(zback)

	use basin
	use levels
	use shympi
	use intp_fem_file

	implicit none

	real zback(nlvdi,nkn)

	real, save, allocatable :: xobs(:), yobs(:)
	real, save, allocatable :: zobss(:,:)
	double precision, save, allocatable :: times(:)
	real, save, allocatable :: bobs(:), zobs(:), zaux(:), zaux2(:)
	real, save, allocatable :: zanal(:)
	real, save, allocatable :: zback2d(:)
	real, save, allocatable :: rl(:),rlmax(:),seo(:),seb(:)
	real, save, allocatable :: iuse(:)
	logical bback
	integer iu
	integer nintp,nfirst,i
	integer nback,nlast,nval
	integer, save :: nvar,idobs
	integer, save :: icall = 0
	integer, save :: iact = 0
	integer, save :: nobs = 0
	integer, save :: nmax = 0
	real rl0,rlmax0,seo0,seb0
	real, parameter :: flag = -999.
	double precision dtime,atime0
	double precision, save :: tstart,tend
	character*80 file_coords
	character*80 file_obs
	character*20 aline
	logical, save :: bdebug = .true.
	logical bwrite

	logical is_time_absolute
	double precision rd_intp_neville

	if( icall == -1 ) return

!---------------------------------------------------------------
! check basics and set parameters
!---------------------------------------------------------------

	if( nlvdi > 1 ) then
	  write(6,*) 'nlvdi = ',nlvdi
	  stop 'error stop scalar_assimilation: only 2d ready'
	end if

	call get_act_dtime(dtime)

	nintp = 2
	nback = nkn
	bback = .true.

	rl0 = 1000.
	rlmax0 = 3.*rl0

	seo0 = 0.1		!standard error observations
	seb0 = 100.*seo0	!standard error background

	file_obs = 'obsall.txt'
	file_coords = 'obs_coords.txt'

	bwrite = bdebug .and. my_id == 0

!---------------------------------------------------------------
! initialize at first call
!---------------------------------------------------------------

	if( icall == 0 ) then

	  nobs = 0
	  iu = 1
	  open(iu,file=file_coords,status='old',form='formatted')
	  call read_coords(iu,nobs,xobs,yobs)	!only count
	  if( nobs == 0 ) icall = -1
	  if( icall == -1 ) return
	  allocate(xobs(nobs),yobs(nobs),zobs(nobs))
	  allocate(bobs(nobs),zaux(nobs),zaux2(nobs))
	  rewind(iu)
	  call read_coords(iu,nobs,xobs,yobs)	!read coords
	  close(iu)

	  nmax = 0
	  iu = 1
	  open(iu,file=file_obs,status='old',form='formatted')
	  call read_timeseries(iu,nmax,nval,times,zobss)	!only count
	  if( nmax == 0 ) icall = -1
	  if( icall == -1 ) return
	  allocate(times(nmax),zobss(nmax,nval))
	  zobss = 0
	  rewind(iu)
	  call read_timeseries(iu,nmax,nval,times,zobss)	!read obs
	  close(iu)

	  if( is_time_absolute(times(1)) ) then
	    call get_absolute_ref_time(atime0)
	    times = times - atime0		!convert to simulation time
	  end if

	  tstart = times(1)
	  tend = times(nmax)

	  if( nobs /= nval ) then
	    write(6,*) 'nobs,nval: ',nobs,nval
	    write(6,*) 'number of coordinates and observations are different'
	    stop 'error stop scalar_assimilation: nobs/=nval'
	  end if

	  write(6,*) 'observations have been read: ',nobs

	  nvar = nobs
	  call iff_ts_init(dtime,file_obs,nintp,nvar,idobs)
	!write(6,*) idobs

	  if( bdebug ) then
	    if( my_id == 0 ) write(666,*) 'starting optintp: ',nobs
	  end if

	  allocate(zanal(nkn))
	  allocate(zback2d(nkn))

	  allocate(rl(nobs),rlmax(nobs),seo(nobs),seb(nobs))
	  allocate(iuse(nobs))
	  iuse = 1
	  rl = rl0
	  rlmax = rlmax0
	  seo = seo0
	  seb = seb0
	end if

	icall = icall + 1
	iuse(5) = 0

	!call iff_ts_intp(idobs,dtime,zaux2)
	if( bdebug ) call get_timeline(dtime,aline)
	if( iff_file_has_data(idobs,dtime) ) then
	  if( bwrite ) write(666,*) 'file has data  ',aline,dtime
	else
	  if( bwrite ) write(666,*) 'file has no data  ',aline,dtime
	  return
	end if

!---------------------------------------------------------------
! see if in observation window
!---------------------------------------------------------------

	if( dtime < tstart .or. dtime > tend ) return

	if( bdebug ) then
	  if( my_id == 0 ) then
	    write(666,*) 'doing optintp',nmax,dtime
	  end if
	end if

!---------------------------------------------------------------
! interpolate observations in time
!---------------------------------------------------------------

	call iff_ts_intp(idobs,dtime,zobs)
	where( iuse == 0 ) zobs = flag

!---------------------------------------------------------------
! compute background values at observation points
!---------------------------------------------------------------

	zback2d = zback(1,:)
	call get_bobs(nobs,xobs,yobs,zback2d,bobs)

!---------------------------------------------------------------
! do optimal interpolation
!---------------------------------------------------------------

	call opt_intp(nobs,xobs,yobs,zobs,bobs                          &
     &                  ,nback,bback,xgv,ygv,zback2d                  &
     &                  ,rl,rlmax,seb,seo,zanal)

!---------------------------------------------------------------
! copy analysis back to background
!---------------------------------------------------------------

!	zobs is old way of observation interpolation
!	zaux2 is new way of observation interpolation
!	bobs is background at observation points
!	zaux is analysis at obervation points

	if( bdebug ) then
	  call get_bobs(nobs,xobs,yobs,zanal,zaux)
	  if( my_id == 0 ) then
	    write(666,*) 'final: '
	    do i=1,nobs
	      write(666,*) i,zobs(i),bobs(i),zaux(i)
	    end do
	  end if
	  !if( icall > 10 ) stop
	end if

	zback(1,:) = zanal

!---------------------------------------------------------------
! end of routine
!---------------------------------------------------------------

	end

!*******************************************************************

	subroutine scalar_nudging(dt,scal,sobs,rtau)

! does nudging of scalars

	use basin
	use levels

	implicit none

	real, intent(in) :: dt			!time step
	real, intent(inout) :: scal(nlvdi,nkn)	!scalar
	real, intent(in) :: sobs(nlvdi,nkn)	!observation of scalar
	real, intent(in) :: rtau(nlvdi,nkn)	!inverse of time scale tau

	integer k,l,lmax
	real r,rr

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    r = dt * rtau(l,k)
	    rr = 1./(1.+r)
	    scal(l,k) = rr*scal(l,k) + (1.-rr)*sobs(l,k)
	  end do
	end do

	end

!*******************************************************************

	subroutine read_coords(iu,nobs,xobs,yobs)

	implicit none

	integer iu
	integer nobs
	real xobs(nobs), yobs(nobs)

	integer i,j
	integer ios
	real x,y
	character*80 name

	i = 0
	do
	  read(iu,*,iostat=ios) j,x,y
	  if( ios < 0 ) exit
	  if( ios > 0 ) goto 99
	  i = i + 1
	  if( i /= j ) goto 98
	  if( nobs == 0 ) cycle
	  if( i > nobs) goto 97
	  xobs(i) = x
	  yobs(i) = y
	end do

	nobs = i

	return
   97	continue
	write(6,*) 'i,nobs: ',i,nobs
	stop 'error stop read_coords: i>nobs'
   98	continue
	write(6,*) 'i,j: ',i,j
	stop 'error stop read_coords: i/=j'
   99	continue
	inquire(iu,name=name)
	write(6,*) 'reading file ',trim(name)
	write(6,*) 'read error close to line ',i
	stop 'error stop read_coords: read error'
	end

!*******************************************************************

	subroutine get_bobs(nobs,xobs,yobs,zback,bobs)

! get value of the background grid on observation points

	use basin
	use levels
	use shympi

	implicit none

	integer nobs
	real xobs(nobs)
	real yobs(nobs)
	real zback(nlvdi,nkn)
	real bobs(nobs)

	integer i,ie,n,ii,k
	integer, save :: icall = 0
	integer, save, allocatable :: ies(:)	!element where obs point is in
	integer iaux(nobs)
	real zp
	real z(nobs)

	if( icall == 0 ) then
	  allocate(ies(nobs))
	  ies = 0
	  iaux = 0

	  do i=1,nobs
	    call find_element(xobs(i),yobs(i),ie)
	    if( ie <= 0 ) cycle
	    if( ie > nel_unique ) cycle
	    ies(i) = ie
	    iaux(i) = 1
	  end do

	  call shympi_gather_and_sum(iaux)
	  n = sum(iaux)
	  if( any(iaux==0) .or. n /= nobs ) then
	    write(6,*) n,nobs
	    write(6,*) iaux
	    write(6,*) ies
	    write(6,*) 'some elements could not be found'
	    stop 'error stop get_bobs: elements not found'
	  end if
	end if

	icall = icall + 1

	bobs = 0.
	do i=1,nobs
	  ie = ies(i)
	  if( ie == 0 ) cycle
	  do ii=1,3
	    k = nen3v(ii,ie)
	    z(ii) = zback(1,k)
	  end do
	  call femintp(ie,z,xobs(i),yobs(i),zp)
	  bobs(i) = zp
	end do

	call shympi_gather_and_sum(bobs)

	end

!*******************************************************************

