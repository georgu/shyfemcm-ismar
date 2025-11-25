
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
! 10.11.2025    ggu     general assimilation of scalar variables
! 14.11.2025    ggu     set up debug file
! 23.11.2025    ggu     bug fix for mode=1
!
!****************************************************************

!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine scalar_assimilation(idass,zback)

	use basin
	use levels
	use shympi
	use intp_fem_file
	use mod_assimil_admin

	implicit none

	integer idass
	real zback(nlvdi,nkn)

	real, allocatable :: xobs(:), yobs(:)
	real, allocatable :: bobs(:), zobs(:), zaux(:)
	real, save, allocatable :: zanal(:)
	real, save, allocatable :: zback2d(:)
	logical bback
	logical basstime
	integer iu
	integer, save :: iu666
	integer nintp,nfirst,i
	integer mode
	integer nback,nlast,nval,ivar,iusec,inofc
	integer, save :: nvar,idobs
	integer, save :: icall = 0
	integer, save :: iact = 0
	integer, save :: nobs = 0
	real tau,rtau,dt,tacum
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
	integer ifemop
	double precision rd_intp_neville

	if( icall == -1 ) return
	if( idass == 0 ) return

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

	bwrite = bdebug .and. my_id == 0

!---------------------------------------------------------------
! initialize at first call
!---------------------------------------------------------------

	if( icall == 0 ) then		!global initialization

	  nobs = passim(idass)%nobs
	  if( nobs <= 0 ) icall = -1
	  if( icall == -1 ) return

	  allocate(zanal(nkn))
	  allocate(zback2d(nkn))

	  iu666 = ifemop('ggu','formatted','unknown')
	  if( iu666 <= 0 ) stop 'iu666'

	end if

	if( .not. passim(idass)%binit ) then

	  file_obs = passim(idass)%fileobs
	  file_coords = passim(idass)%filecoords

	  nobs = 0
	  iu = 1
	  open(iu,file=file_coords,status='old',form='formatted')
	  call read_coords(iu,nobs,xobs,yobs)	!only count
	  if( nobs == 0 ) icall = -1
	  if( icall == -1 ) return
	  allocate(xobs(nobs),yobs(nobs))
	  rewind(iu)
	  call read_coords(iu,nobs,xobs,yobs)	!read coords
	  close(iu)

	  if( nobs /= passim(idass)%nobs ) then
	    write(6,*) 'nobs: ',nobs,passim(idass)%nobs
	    stop 'error stop scalar_assimilation: nobs incompatible'
	  end if

	  allocate(passim(idass)%xobs(nobs))
	  allocate(passim(idass)%yobs(nobs))
	  allocate(passim(idass)%ies(nobs))
	  passim(idass)%xobs = xobs
	  passim(idass)%yobs = yobs
	  passim(idass)%ies = -1		!still to be initialized
	  deallocate(xobs,yobs)

	  write(6,*) 'observations have been read: ',nobs

	  nvar = nobs
	  call iff_ts_init(dtime,file_obs,nintp,nvar,idobs)
	  passim(idass)%idobs = idobs

	  if( nobs /= nvar ) then
	    write(6,*) 'nobs,nval: ',nobs,nval
	    write(6,*) 'number of coordinates and observations are different'
	    stop 'error stop scalar_assimilation: nobs/=nval'
	  end if

	  if( bwrite ) write(iu666,*) 'starting optintp: ',nobs

	  passim(idass)%binit = .true.
	end if

	icall = icall + 1

	nobs = passim(idass)%nobs
	if( nobs <= 0 ) return
	ivar = passim(idass)%ivar
	idobs = passim(idass)%idobs

!---------------------------------------------------------------
! see if in observation window
!---------------------------------------------------------------

	if( bdebug ) then
	  call get_timeline(dtime,aline)
	  if( iff_file_has_data(idobs,dtime) ) then
	    if( bwrite ) write(iu666,*) 'file has data  ',aline,dtime,ivar
	  else
	    if( bwrite ) write(iu666,*) 'file has no data  ',aline,dtime,ivar
	    return
	  end if
	end if

	if( bwrite ) write(iu666,*) 'doing optintp',dtime

!---------------------------------------------------------------
! see if assimilation time step has arrived
!---------------------------------------------------------------

	call get_timestep(dt)
	tau = passim(idass)%tau
	tacum = passim(idass)%tacum
	mode = passim(idass)%mode

	if( mode == 1 ) then
	  tacum = tacum + dt
	  basstime = ( tacum >= tau )	!do assimilation at this time
	  if( basstime ) tacum = tacum - tau
	  passim(idass)%tacum = tacum
	  if( bwrite ) then
	    write(iu666,*) 'checking assimilation in time',basstime,tacum
	  end if
	  if( .not. basstime ) return
	end if

	if( bwrite ) write(iu666,*) 'doing assimilation in time',dtime

!---------------------------------------------------------------
! allocate arrays
!---------------------------------------------------------------

	allocate(xobs(nobs),yobs(nobs))
	allocate(bobs(nobs),zobs(nobs),zaux(nobs))

	xobs = passim(idass)%xobs
	yobs = passim(idass)%yobs

!---------------------------------------------------------------
! interpolate observations in time
!---------------------------------------------------------------

	idobs = passim(idass)%idobs
	call iff_ts_intp(idobs,dtime,zobs)
	where( passim(idass)%iuse == 0 ) zobs = flag

!---------------------------------------------------------------
! compute background values at observation points
!---------------------------------------------------------------

	zback2d = zback(1,:)
	call get_bobs(nobs,xobs,yobs,zback2d,bobs,passim(idass)%ies)

!---------------------------------------------------------------
! do optimal interpolation
!---------------------------------------------------------------

	!write(6,*) passim(idass)%sigback
	!write(6,*) passim(idass)%sigobs
	call opt_intp(nobs,xobs,yobs,zobs,bobs                          &
     &                  ,nback,bback,xgv,ygv,zback2d			&
     &                  ,passim(idass)%corlength			&
     &                  ,passim(idass)%maxlength			&
     &                  ,passim(idass)%sigback				&
     &                  ,passim(idass)%sigobs				&
     &                  ,zanal)

!---------------------------------------------------------------
! write debug message to file
!---------------------------------------------------------------

!	zobs is observations interpolated at actual time
!	bobs is background at observation points
!	zaux is analysis at obervation points

	if( bdebug ) then
	  call get_bobs(nobs,xobs,yobs,zanal,zaux,passim(idass)%ies)
	  if( my_id == 0 ) then
	    iusec = sum(passim(idass)%iuse)
	    inofc = count(zobs/=flag)
	    write(iu666,*) 'final optimal interpolation: ',nobs,iusec,inofc,ivar
	    do i=1,nobs
	      if( zobs(i) == flag ) then
		if( bobs(i) /= zaux(i) ) then
		  write(6,*) i,bobs(i),zaux(i)
		  stop 'error stop scalar_assimilation: bobs/=zanal'
		end if
		cycle
	      end if
	      write(iu666,*) i,zobs(i),bobs(i),zaux(i)
	    end do
	  end if
	end if

!---------------------------------------------------------------
! if needed apply nudging
!---------------------------------------------------------------

	if( mode == 0 .and. tau > 0. ) then
	  if( bwrite ) write(iu666,*) 'doing nudging'
	  rtau = 1./tau
	  call scalar_nudging_0(dt,zback2d,zanal,rtau)
	else if( mode == 1 ) then
	  zback2d = zanal
	else
	  write(6,*) 'error in mode: ',mode,tau
	  stop 'error stop scalar_assimilation: error in mode'
	end if

!---------------------------------------------------------------
! copy analysis back to background
!---------------------------------------------------------------

	zback(1,:) = zback2d

!---------------------------------------------------------------
! end of routine
!---------------------------------------------------------------

	end

!*******************************************************************

	subroutine scalar_nudging(dt,scal,sobs,rtau)

! does nudging of scalars - implicit - rtau matrix version

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

	subroutine scalar_nudging_0(dt,scal,sobs,rtau)

! does nudging of scalars - implicit - rtau scalar version

	use basin
	use levels

	implicit none

	real, intent(in) :: dt			!time step
	real, intent(inout) :: scal(nlvdi,nkn)	!scalar
	real, intent(in) :: sobs(nlvdi,nkn)	!observation of scalar
	real, intent(in) :: rtau		!inverse of time scale tau

	integer k,l,lmax
	real r,rr

	r = dt * rtau
	rr = 1./(1.+r)

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
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

	subroutine get_bobs(nobs,xobs,yobs,zback,bobs,ies)

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
	integer ies(nobs)	!element where obs point is in

	integer i,ie,n,ii,k
	integer, save :: icall = 0
	integer iaux(nobs)
	real zp
	real z(nobs)

	if( ies(1) == -1 ) then
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

        subroutine assimil_ts(tarray,sarray)

	use basin
	use levels
	use shympi
	use mod_assimil_admin

        implicit none

	real tarray(nlvdi,nkn)
	real sarray(nlvdi,nkn)

	integer, save :: icall = 0
	integer, save :: idtemp = 0
	integer, save :: idsalt = 0

	if( icall == -1 ) return

	if( icall == 0 ) then
	  idtemp = get_id_from_var(12)
	  idsalt = get_id_from_var(11)
	end if

	icall = icall + 1
	
	if( idtemp /= 0 ) call scalar_assimilation(idtemp,tarray)
	if( idsalt /= 0 ) call scalar_assimilation(idsalt,sarray)

        end

!*******************************************************************

        subroutine assimil_init_ids(nvar,ivars,ids)

! initializes array ids

	use mod_assimil_admin

	implicit none

	integer, intent(in) 	:: nvar
	integer, intent(in) 	:: ivars(nvar)
	integer, intent(inout) 	:: ids(3,nvar)

	integer id,iv,ivar
	integer idmax

	ids = -1
	idmax = 0
	do iv=1,nvar
	  id = get_id_from_var(ivars(iv))
	  if( id > 0 ) then
	    idmax = idmax + 1
	    ids(1,idmax) = id
	    ids(2,idmax) = iv
	    ids(3,idmax) = ivars(iv)
	  end if
	end do

	end

!*******************************************************************

        subroutine assimil_3dvars(nvar,ivars,array,ids)

! goes through all variables and checks if assimilation has to be carried out
!
! ids on first call must be set to 0

	use basin
	use levels
	use shympi
	use mod_assimil_admin

        implicit none

	integer, intent(in) 	:: nvar
	integer, intent(in) 	:: ivars(nvar)
	real, intent(in) 	:: array(nlvdi,nkn,nvar)
	integer, intent(inout) 	:: ids(3,nvar)

	integer i,id,iv,ivar
	integer idmax

	if( ids(1,1) == -1 ) return

	if( ids(1,1) == 0 ) call assimil_init_ids(nvar,ivars,ids)

	do i=1,nvar
	  id = ids(1,i)
	  if( id < 0 ) exit		!no more variables to assimilate
	  iv = ids(2,i)
	  ivar = ids(3,i)
	  write(6,*) 'doing assimilation for ivar = ',ivar
	  call scalar_assimilation(id,array(:,:,iv))
	end do

        end

!*******************************************************************

        subroutine assimil_2dvars(nvar,ivars,array,ids)

! goes through all variables and checks if assimilation has to be carried out
!
! ids on first call must be set to 0

	use basin
	use levels
	use shympi
	use mod_assimil_admin

        implicit none

	integer, intent(in) 	:: nvar
	integer, intent(in) 	:: ivars(nvar)
	real, intent(in) 	:: array(nkn,nvar)
	integer, intent(inout) 	:: ids(2,nvar)

	integer i,id,iv,ivar
	integer idmax

	if( ids(1,1) == -1 ) return

	if( ids(1,1) == 0 ) call assimil_init_ids(nvar,ivars,ids)

	do i=1,nvar
	  id = ids(1,i)
	  if( id < 0 ) exit		!no more variables to assimilate
	  iv = ids(2,i)
	  ivar = ivars(iv)
	  write(6,*) 'doing assimilation for ivar = ',ivar
	  call scalar_assimilation(id,array(:,iv))
	end do

        end

!*******************************************************************


