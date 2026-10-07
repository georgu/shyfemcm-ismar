
!--------------------------------------------------------------------------
!
!    Copyright (C) 2013-2015,2017-2020  Georg Umgiesser
!    Copyright (C) 2015  Christian Ferrarin
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

! routines for offline data handling
!
! revision log :
!
! 13.06.2013	ggu	new routines written from scratch
! 17.06.2013	ggu	eliminated compiler warnings
! 25.03.2014	ggu	new offline (for T/S)
! 19.12.2014	ggu	changed VERS_7_0_10
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 01.04.2015	ggu	changed VERS_7_1_7
! 06.05.2015	ccf	write offline to .off file
! 06.05.2015	ccf	read offline from offlin file in section name
! 21.05.2015	ggu	changed VERS_7_1_11
! 10.07.2015	ggu	changed VERS_7_1_50
! 13.07.2015	ggu	changed VERS_7_1_51
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 05.11.2015	ggu	revisited and checked
! 09.11.2015	ggu	changed VERS_7_3_13
! 16.11.2015	ggu	changed VERS_7_3_14
! 29.03.2017	ggu	bug fix - input file opened on unit 1
! 05.12.2017	ggu	changed VERS_7_5_39
! 12.11.2018	ggu	linear arrays introduced
! 18.12.2018	ggu	changed VERS_7_5_52
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 04.07.2019	ggu	solved problem for vertical velocity (1 index off)
! 17.02.2020	ggu	femtime eliminated
! 28.04.2020	ggu	restructured and taken out of suboff.f
! 09.10.2020	ggu	added comments and error checking
! 13.10.2020	ggu	default values for T/S introduced
! 06.04.2021	ggu	better error message in off_check_vertical()
! 07.03.2024	ggu	off write routine prepared for use with MPI, new vers 4
! 13.03.2024	ggu	turbulence implemented
! 15.03.2024	ggu	time to iitime, more on turbulence
! 10.04.2024	ggu	finsihed offline, bug fix when reading ze
!
! contents :
!
! mod_offline_init()	initializes offline module
! off_write_record()	writes one record to offline file
! off_read_record()	reads one record from offline file
! off_read_header()	reads header and sets nlv
! off_peek_header()	peeks into header and sets nlv
! off_init_vertical()	initializes vertical indices
! off_check_vertical()	checks if vertical indices have been initialized
! off_error()		exits with error
! off_next_record()	checks info on next record
!
! calling sequence for writing:
!
!	mod_offline_init()
!	do
!	  off_write_record()
!	end do
!
! calling sequence for reading:
!
!	off_peek_header()
!	mod_offline_init()
!	do
!	  off_read_record()
!	end do
!
! ioffline:
!		1	hydro
!		2	T/S
!		4	turbolence
!		7	all of the above
!
! notes :
!
!	reading not yet ready for mpi
!	should use dtime and not it
!
!****************************************************************

!==================================================================
	module mod_offline
!==================================================================

	implicit none

	integer, parameter :: nintp = 4		!2 (linear) or 4 (cubic)

	integer, parameter :: nvers_max = 5	!last version
	integer, parameter :: ioff_max = 7	!maximum we can do

	logical, save :: bwhydro = .true.	!write hydro results
	logical, save :: bwts = .true.		!write TS results
	logical, save :: bwturb = .true.	!write turbulence results

	integer, save :: nvers = 0		!actual version
	integer, save :: ioffline = 0
	integer, save :: idtoff,itmoff,itoff
	!double precision, save :: idtoff,itmoff,itoff
	integer, save :: iread = 0		!what is in file
	integer, save :: iwrite = 0		!what has been written to file
	integer, save :: iuoff			!unit to read/write
	integer, save :: icall = 0
	integer, save :: ioff_mode = 0		!read/write mode
	logical, save :: bfirst = .true.
	logical, save :: bdebug = .false.

	integer, save :: nkn_off = 0		!values of local domain
	integer, save :: nel_off = 0
	integer, save :: nlv_off = 0

	double precision, save :: dtr = 0.
	integer, save :: iitime(nintp)

	integer, save :: idef = 0		!use default values for T/S
	real, save :: tdef = 0
	real, save :: sdef = 0

	double precision, save, allocatable :: ut(:,:,:)
	double precision, save, allocatable :: vt(:,:,:)
	double precision, save, allocatable :: ze(:,:,:)
	double precision, save, allocatable :: wn(:,:,:)
	double precision, save, allocatable :: zn(:,:)
	double precision, save, allocatable :: sn(:,:,:)
	double precision, save, allocatable :: tn(:,:,:)
	double precision, save, allocatable :: vd(:,:,:)	!visv
	double precision, save, allocatable :: dd(:,:,:)	!difv

	integer, save, allocatable :: ile(:)
	integer, save, allocatable :: ilk(:)
	logical, save :: bvinit = .false.	!il files initialized

!==================================================================
	contains
!==================================================================

	subroutine mod_offline_init(nk,ne,nl)

! initializes offline module

	integer nk,ne,nl

        if( nk == nkn_off .and. ne == nel_off &
     &          .and. nl == nlv_off ) return

        if( nk > 0 .or. ne > 0 .or. nl > 0 ) then
          if( nk == 0 .or. ne == 0 .or. nl == 0 ) then
            write(6,*) 'nk,ne,nl: ',nk,ne,nl
            stop 'error stop mod_offline_init: incompatible parameters'
          end if
        end if

	if( nkn_off > 0 ) then
	  deallocate(ut)
	  deallocate(vt)
	  deallocate(ze)
	  deallocate(wn)
	  deallocate(zn)
	  deallocate(sn)
	  deallocate(tn)
	  deallocate(vd)
	  deallocate(dd)
	  deallocate(ile)
	  deallocate(ilk)
	end if

	nkn_off = nk
	nel_off = ne
	nlv_off = nl

	allocate(ut(nl,ne,nintp))
	allocate(vt(nl,ne,nintp))
	allocate(ze(3,ne,nintp))
	allocate(wn(0:nl,nk,nintp))
	allocate(zn(nk,nintp))
	allocate(sn(nl,nk,nintp))
	allocate(tn(nl,nk,nintp))
	allocate(vd(0:nl,nk,nintp))
	allocate(dd(0:nl,nk,nintp))
	allocate(ile(ne))
	allocate(ilk(nk))

	ut = 0.
	vt = 0.
	ze = 0.
	wn = 0.
	zn = 0.
	sn = 0.
	tn = 0.
	vd = 0.
	dd = 0.

	end subroutine mod_offline_init

!==================================================================
	end module mod_offline
!==================================================================

!****************************************************************
	
	subroutine off_write_record(iu,it)

! writes one record to offline file

	use mod_offline
	use shympi
	use mod_trace_point

	implicit none

	integer iu
	integer it

	logical bwrite,bvers4,bvers5,b3d
	integer ie,ii,k,i
	integer nlin,nlink,nline,nlinv
	integer iunit
	integer nkn,nel,nlv
	integer nkng,nelg,nlvg
	double precision dtime
	
	integer, allocatable :: ilkg(:), ileg(:)
	double precision, save, allocatable :: wnaux(:,:)
	double precision, save, allocatable :: dlin(:)
	double precision, save, allocatable :: delemg(:,:)
	double precision, save, allocatable :: dnodeg(:,:)
	double precision, save, allocatable :: dezg(:,:)
	double precision, save, allocatable :: dkzg(:)

	if( .not. bvinit ) then
	  stop 'error stop off_write_record: bvinit is false'
	else if( nkn_off <= 0  ) then
	  stop 'error stop off_write_record: offline not initialized'
	end if

!----------------------------------------------------------
! initialize
!----------------------------------------------------------

	bvers4 = ( nvers_max > 3 )
	bvers5 = ( nvers_max == 5 )
	bwrite = shympi_is_master()

	iunit = iu
	nkn = nkn_off
	nel = nel_off
	nlv = nlv_off

	nkng = nkn_global
	nelg = nel_global
	nlvg = nlv_global

	b3d = ( nlvg > 1 )

!----------------------------------------------------------
! set up auxiliary arrays
!----------------------------------------------------------

	allocate(ilkg(nkng),ileg(nelg))
	call shympi_l2g_array(ilk,ilkg)
	call shympi_l2g_array(ile,ileg)
        call count_linear(nlvg,nkng,1,ilkg,nlink)
        call count_linear(nlvg,nelg,1,ileg,nline)
	nlin = max(nlink,nline)

        if( .not. allocated(dlin) ) allocate(dlin(nlin))
        if( .not. allocated(dnodeg) ) allocate(dnodeg(nlvg,nkng))
        if( .not. allocated(delemg) ) allocate(delemg(nlvg,nelg))
        if( .not. allocated(dezg) ) allocate(dezg(3,nelg))
        if( .not. allocated(dkzg) ) allocate(dkzg(nkng))
        if( .not. allocated(wnaux) ) allocate(wnaux(nlv,nkn))

!----------------------------------------------------------
! write header
!----------------------------------------------------------

	!write(6,*) 'offline nlv: ',my_id,nlv,nlvg

	dtime = it

	bwturb = ( bwturb .and. b3d )	! do turb only for 3d
	if( nvers_max == 3 ) bwturb = .false.

	iwrite = 0
	if( bwhydro ) iwrite = iwrite + 1
	if( bwts ) iwrite = iwrite + 2
	if( bwturb ) iwrite = iwrite + 4

	if( bwrite ) then
	  write(iunit) it,nkng,nelg,nvers_max
	  if( bvers4 ) write(iunit) nlvg,dtime		! new version 4
	  if( bvers4 ) write(iunit) nlink,nline		! new version 4
	  if( bvers5 ) write(iunit) iwrite		! new version 5
	  write(iunit) (ileg(ie),ie=1,nelg)
	  write(iunit) (ilkg(k),k=1,nkng)
	end if

!----------------------------------------------------------
! write currents
!----------------------------------------------------------

	if( bwhydro ) then
	  call trace_point('offline writing hydro')
          nlin = nline
	  call shympi_l2g_array(ut(:,:,1),delemg)
          call dvals2linear(nlvg,nelg,1,ileg,delemg,dlin,nlin)
          if( bwrite ) write(iunit) (dlin(i),i=1,nlin)
	  call shympi_l2g_array(vt(:,:,1),delemg)
          call dvals2linear(nlvg,nelg,1,ileg,delemg,dlin,nlin)
          if( bwrite ) write(iunit) (dlin(i),i=1,nlin)
	end if

!----------------------------------------------------------
! write water levels and vertical velocities
!----------------------------------------------------------

	if( bwhydro ) then
	  call trace_point('offline writing levels')
          nlin = nlink
	  call shympi_l2g_array(3,ze(:,:,1),dezg)
	  if( bwrite ) write(iunit) ((dezg(ii,ie),ii=1,3),ie=1,nelg)
          wnaux(1:nlv,:) = wn(1:nlv,:,1)
	  call shympi_l2g_array(wnaux,dnodeg)
          call dvals2linear(nlvg,nkng,1,ilkg,dnodeg,dlin,nlin)
          if( bwrite ) write(iunit) (dlin(i),i=1,nlin)
	  call shympi_l2g_array(zn(:,1),dkzg)
	  if( bwrite ) write(iunit) (dkzg(k),k=1,nkng)
	end if

!----------------------------------------------------------
! write T/S
!----------------------------------------------------------

	if( bwts ) then
	  call trace_point('offline writing T/S')
          nlin = nlink
	  call shympi_l2g_array(sn(:,:,1),dnodeg)
          call dvals2linear(nlvg,nkng,1,ilkg,dnodeg,dlin,nlin)
          if( bwrite ) write(iunit) (dlin(i),i=1,nlin)
	  call shympi_l2g_array(tn(:,:,1),dnodeg)
          call dvals2linear(nlvg,nkng,1,ilkg,dnodeg,dlin,nlin)
          if( bwrite ) write(iunit) (dlin(i),i=1,nlin)
	end if

!----------------------------------------------------------
! write turbulence
!----------------------------------------------------------

	if( bwturb ) then
	  call trace_point('offline writing turbulence')
          nlin = nlink
          wnaux(1:nlv,:) = vd(1:nlv,:,1)
	  call shympi_l2g_array(wnaux,dnodeg)
          call dvals2linear(nlvg,nkng,1,ilkg,dnodeg,dlin,nlin)
          if( bwrite ) write(iunit) (dlin(i),i=1,nlin)
          wnaux(1:nlv,:) = dd(1:nlv,:,1)
	  call shympi_l2g_array(wnaux,dnodeg)
          call dvals2linear(nlvg,nkng,1,ilkg,dnodeg,dlin,nlin)
          if( bwrite ) write(iunit) (dlin(i),i=1,nlin)
	end if

!----------------------------------------------------------
! end of routine
!----------------------------------------------------------

	end

!****************************************************************

	subroutine off_read_record(iu,ig,it,ierr)

! reads one record from offline file
!
! this is still not working with mpi

	use mod_offline
	use shympi
	use mod_trace_point

	implicit none

	integer iu,ig
	integer it
	integer ierr

	logical bread
	integer ie,ii,k,i
	integer nlin,nlink,nline
	integer iunit
	integer nkng,nelg,nlvg
	integer nlinkaux,nlineaux
	integer nkn,nel,nlv
	
	integer, save, allocatable :: ileaux(:)
	integer, save, allocatable :: ilkaux(:)

	integer, allocatable :: ilkg(:), ileg(:)
	double precision, save, allocatable :: wnaux(:,:)
	double precision, save, allocatable :: dlin(:)
	double precision, save, allocatable :: delemg(:,:)
	double precision, save, allocatable :: dnodeg(:,:)
	double precision, save, allocatable :: dezg(:,:)
	double precision, save, allocatable :: dkzg(:)

	logical off_has_record

	if( nkn_off <= 0  ) then
	  stop 'error stop off_read_record: offline not initialized'
	end if

!----------------------------------------------------------
! initialize
!----------------------------------------------------------

	nkn = nkn_off		! these are values of local domain
	nel = nel_off
	nlv = nlv_off

	iunit = iu

!----------------------------------------------------------
! read header
!----------------------------------------------------------

	call off_read_header(iunit,it,nkng,nelg,nlvg &
     &			,nlinkaux,nlineaux,iread,ierr)
	if( ierr > 0 ) goto 95
	if( ierr < 0 ) goto 98
	if( nkn_global .ne. nkng .or. nel_global .ne. nelg ) goto 97
	if( nlv_global .ne. nlvg ) goto 97

	iitime(ig) = it

!----------------------------------------------------------
! read vertical indices
!----------------------------------------------------------

	if( .not. allocated(ileg) ) allocate(ileg(nelg))
	if( .not. allocated(ilkg) ) allocate(ilkg(nkng))
	if( .not. allocated(ileaux) ) allocate(ileaux(nel))
	if( .not. allocated(ilkaux) ) allocate(ilkaux(nkn))
	read(iunit,iostat=ierr) (ileg(ie),ie=1,nelg)
	call off_error(ierr,it,'reading ile')
	read(iunit,iostat=ierr) (ilkg(k),k=1,nkng)
	call off_error(ierr,it,'reading ilk')

	call shympi_g2l_array(ilkg,ilkaux)
	call shympi_g2l_array(ileg,ileaux)
	if( .not. bvinit ) then
	  call off_init_vertical(nkn,nel,ileaux,ilkaux)
	end if

	call off_check_vertical('element',nel,ileaux,ile)
	call off_check_vertical('node',nkn,ilkaux,ilk)

!----------------------------------------------------------
! set up auxiliary arrays
!----------------------------------------------------------

        call count_linear(nlvg,nkng,1,ilkg,nlink)
        call count_linear(nlvg,nelg,1,ileg,nline)
	if( nlinkaux > 0 .and. nlink /= nlinkaux ) goto 91
	if( nlineaux > 0 .and. nline /= nlineaux ) goto 91
	nlin = max(nlink,nline)

        if( .not. allocated(dlin) ) allocate(dlin(nlin))
        if( .not. allocated(dnodeg) ) allocate(dnodeg(nlvg,nkng))
        if( .not. allocated(delemg) ) allocate(delemg(nlvg,nelg))
        if( .not. allocated(dezg) ) allocate(dezg(3,nelg))
        if( .not. allocated(dkzg) ) allocate(dkzg(nkng))
        if( .not. allocated(wnaux) ) allocate(wnaux(nlv,nkn))

!----------------------------------------------------------
! read currents
!----------------------------------------------------------

	bread = off_has_record(1,iread)

	if( bread ) then
	  call trace_point('offline reading hydro')
	  nlin = nline
          read(iunit,iostat=ierr) (dlin(i),i=1,nlin)
	  call off_error(ierr,it,'reading ut')
          call dlinear2vals(nlvg,nelg,1,ileg,delemg,dlin,nlin)
	  call shympi_g2l_array(delemg,ut(:,:,ig))
          read(iunit,iostat=ierr) (dlin(i),i=1,nlin)
	  call off_error(ierr,it,'reading vt')
          call dlinear2vals(nlvg,nelg,1,ileg,delemg,dlin,nlin)
	  call shympi_g2l_array(delemg,vt(:,:,ig))
	end if

!----------------------------------------------------------
! read water levels and vertical velocities
!----------------------------------------------------------

	bread = off_has_record(1,iread)

	if( bread ) then
	  call trace_point('offline reading levels')
	  nlin = nlink
	  read(iunit,iostat=ierr) ((dezg(ii,ie),ii=1,3),ie=1,nelg)
	  call off_error(ierr,it,'reading ze')
	  call shympi_g2l_array(3,dezg,ze(:,:,ig))

          read(iunit,iostat=ierr) (dlin(i),i=1,nlin)
	  call off_error(ierr,it,'reading wn')
          call dlinear2vals(nlvg,nkng,1,ilkg,dnodeg,dlin,nlin)
	  call shympi_g2l_array(dnodeg,wnaux)
	  wn(1:nlv,:,ig) = wnaux(1:nlv,:)
	  !wn(0,:,ig) = 0.

	  read(iunit,iostat=ierr) (dkzg(k),k=1,nkng)
	  call off_error(ierr,it,'reading zn')
	  call shympi_g2l_array(dkzg,zn(:,ig))
	end if

!----------------------------------------------------------
! read T/S
!----------------------------------------------------------

	bread = off_has_record(2,iread)

	if( bread ) then
	  call trace_point('offline reading T/S')
	  nlin = nlink
          read(iunit,iostat=ierr) (dlin(i),i=1,nlin)
	  call off_error(ierr,it,'reading sn')
          call dlinear2vals(nlvg,nkng,1,ilkg,dnodeg,dlin,nlin)
	  call shympi_g2l_array(dnodeg,sn(:,:,ig))
          read(iunit,iostat=ierr) (dlin(i),i=1,nlin)
	  call off_error(ierr,it,'reading tn')
          call dlinear2vals(nlvg,nkng,1,ilkg,dnodeg,dlin,nlin)
	  call shympi_g2l_array(dnodeg,tn(:,:,ig))
	end if

!----------------------------------------------------------
! write turbulence
!----------------------------------------------------------

	bread = off_has_record(4,iread)

	if( bread ) then
	  call trace_point('offline reading turbulence')
	  nlin = nlink
          read(iunit,iostat=ierr) (dlin(i),i=1,nlin)
	  call off_error(ierr,it,'reading vd')
          call dlinear2vals(nlvg,nkng,1,ilkg,dnodeg,dlin,nlin)
	  call shympi_g2l_array(dnodeg,wnaux)
	  vd(1:nlv,:,ig) = wnaux(1:nlv,:)
          read(iunit,iostat=ierr) (dlin(i),i=1,nlin)
	  call off_error(ierr,it,'reading dd')
          call dlinear2vals(nlvg,nkng,1,ilkg,dnodeg,dlin,nlin)
	  call shympi_g2l_array(dnodeg,wnaux)
	  dd(1:nlv,:,ig) = wnaux(1:nlv,:)
	end if

!----------------------------------------------------------
! if needed set default values for T/S
!----------------------------------------------------------

	if( idef /= 0 ) then
	  sn = sdef
	  tn = tdef
	end if

!----------------------------------------------------------
! end of routine
!----------------------------------------------------------

	!call off_debug_var

	ierr = 0

	return
   91	continue
	write(6,*) 'nlink: ',nlink,nlinkaux
	write(6,*) 'nline: ',nline,nlineaux
	stop 'error stop off_read: inconsistency in nlink/nline'
   95	continue
	stop 'error stop off_read: read error'
   97	continue
	write(6,*) 'nkng,nkn_global: ',nkng,nkn_global
	write(6,*) 'nelg,nel_global: ',nelg,nel_global
	write(6,*) 'nlvg,nlv_global: ',nlvg,nlv_global
	stop 'error stop off_read: parameter mismatch'
   98	continue
	!write(6,*) 'EOF encountered: ',iu,ig
	ierr = -1
	return
	!stop 'error stop off_read: EOF encountered'
   99	continue
	write(6,*) iu,ig
	stop 'error stop off_read: error reading record'
	end

!****************************************************************

	subroutine off_read_header(iu,it,nkn,nel,nlv,nlink,nline &
     &			,itype,ierr)

! reads header and sets nlv

	use mod_offline

	implicit none

	integer iu
	integer it,nkn,nel,nlv
	integer nlink,nline
	integer itype			! what has been read
	integer ierr

	logical bvers4,bvers5
	integer ie,k,n
	integer nlve,nlvk
	double precision dtime

	integer, allocatable :: il(:)

	ierr = 0

	read(iu,err=99,end=98) it,nkn,nel,nvers
	if( nvers < 3 ) goto 96
	if( nvers > nvers_max ) goto 96
	bvers4 = ( nvers >= 4 )
	bvers5 = ( nvers >= 5 )

	nlink = 0
	nline = 0
	itype = 3
	if( bvers4 ) then
	  read(iu,err=97,end=97) nlv,dtime
	  read(iu,err=97,end=97) nlink,nline
	  if( bvers5 ) then
	    read(iu,err=97,end=97) itype
	  end if
	end if

	if( bvers4 ) return

! determine nlv - this should be deleted once we write nlv to file

	n = max(nkn,nel)
	allocate(il(n))
	read(iu,iostat=ierr) (il(ie),ie=1,nel)
	call off_error(ierr,it,'off_read_header: reading file')
	nlve = maxval(il(1:nel))
	read(iu,iostat=ierr) (il(k),k=1,nkn)
	call off_error(ierr,it,'off_read_header: reading file')
	nlvk = maxval(il(1:nkn))

	nlv = max(nlve,nlvk)

	return
   96	continue
	write(6,*) 'nvers: ',nvers
	write(6,*) 'allowed nvers: >= 3 and <= ',nvers_max
	stop 'error stop off_read_header: nvers'
   97	continue
	write(6,*) iu
	stop 'error stop off_read_header: error reading second header'
   98	continue
	write(6,*) 'EOF encountered: ',iu
	ierr = -1
	return
   99	continue
	write(6,*) iu
	stop 'error stop off_read_header: error reading header'
	end

!****************************************************************

	subroutine off_peek_header(iu,it,nkn,nel,nlv,nlink,nline &
     &					,itype,ierr)

! peeks into header and sets nlv

	use mod_offline

	implicit none

	integer iu
	integer it,nkn,nel,nlv
	integer nlink,nline
	integer itype
	integer ierr

	call off_read_header(iu,it,nkn,nel,nlv,nlink,nline &
     &				,itype,ierr)

	rewind(iu)

	end

!****************************************************************

	subroutine off_peek_next_record(iu,it,ierr)

! checks info on next record

	implicit none

	integer iu,it,ierr

	integer nknaux,nelaux,nv
	integer nlvaux
	double precision dtime

	read(iu,err=99,end=98) it,nknaux,nelaux,nv

	if( nv >= 4 ) then
	  read(iu,err=97,end=97) nlvaux,dtime
	  backspace(iu)
	else
	  dtime = it
	end if

	backspace(iu)
	ierr = 0

	return
   97	continue
	write(6,*) iu,nv
	stop 'error stop off_peek_next_record: error reading 2 record'
   98	continue
	it = 0
	ierr = -1
	return
   99	continue
	write(6,*) iu
	stop 'error stop off_peek_next_record: error reading record'
	end

!****************************************************************

	subroutine off_init_vertical(nkn,nel,ilhv,ilhkv)

! initializes vertical indices

	use mod_offline

	implicit none

	integer nkn,nel
	integer ilhv(nel)
	integer ilhkv(nkn)

	if( nkn /= nkn_off ) goto 99
	if( nel /= nel_off ) goto 99

	ile = ilhv
	ilk = ilhkv

	bvinit = .true.

	return
   99	continue
	write(6,*) 'nkn,nel: ',nkn,nel
	write(6,*) 'nkn_off,nel_off: ',nkn_off,nel_off
	stop 'error stop off_init_vertical: non compatible'
	end

!****************************************************************

	subroutine off_check_vertical(text,n,ilaux,il)

! checks if vertical indices have been initialized

	implicit none

	character*(*) text
	integer n
	integer ilaux(n)
	integer il(n)

	integer i

	do i=1,n
	  if( il(i) .le. 0 ) il(i) = ilaux(i)
	  if( il(i) .ne. ilaux(i) ) goto 99
	end do

	return
   99	continue
	write(6,*) 'checking layers for ',trim(text)
	write(6,*) 'different number of layers found in'
	write(6,*) trim(text),i,' (internal number)'
	write(6,*) 'layers: ',il(i),ilaux(i)
	stop 'error stop off_check_vertical: not compatible'
	end 

!****************************************************************
!****************************************************************
!****************************************************************

        function is_time_for_offline()

! this checks if it is time for offline read/write

        use mod_offline

        implicit none

        logical is_time_for_offline    !it is offline time

	integer it
	double precision dtime

	is_time_for_offline = .false.

	if( ioff_mode == 0 ) return

        call get_act_dtime(dtime)
        it = nint(dtime)

	if( ioff_mode == 1 ) then			!write
	  if( mod(it-itmoff,idtoff) /= 0 ) return
	else if( ioff_mode == 2 ) then			!read
	  write(6,*) 'iitime: ',it,iitime
	  if( count( iitime==it ) == 0 ) return
	else
	  stop 'error stop is_time_for_offline: ioff_mode'
	end if

	is_time_for_offline = .true.

	end

!****************************************************************

	function off_has_record(iwant,iread)

	implicit none

	logical off_has_record
	integer iwant		! what we want
	integer iread		! what we have read

	off_has_record = ( mod(iread/iwant,2) /= 0 )

	end

!****************************************************************

        subroutine is_offline(itype,boff)

! this checks if the requested offline data (itype) is wanted and available
!
! itype: 1 hydro, 2 T/S, 4 turb, combinations are possible: 3,7
! itype == 0 -> any offline
!
! iwant is what we want for offline
! iread is what has been read

        use mod_offline

        implicit none

        integer itype    !should we use this offline data?
        logical boff    !data is available and should be used (return)

        integer iwant

        iwant = ioffline                !this is what we want (from idtoff)

        if( itype < 1 .and. itype > 4 ) then         	!check itype
          write(6,*) 'value for itype not allowed: ',itype
          stop 'error stop is_offline: itype'
        else if( iwant .le. 0 ) then        		!no offline
          boff = .false.
        else
          boff = mod(iwant/itype,2) .ne. 0
        end if

        end

!****************************************************************

        subroutine can_do_offline

! checks if all offline data needed is in file

        use mod_offline

        implicit none

        logical bwant,bread
        integer iwant,i

        iwant = ioffline        !this is what we want

        i = 1
        do while( i .le. 4 )
          bwant = mod(iwant/i,2) .ne. 0
          bread = mod(iread/i,2) .ne. 0
          if( bwant .and. .not. bread ) goto 99
          i = i * 2
        end do

        return
   99   continue
        write(6,*) 'iread = ',iread,'  iwant = ',iwant
        write(6,*) 'type = ',i
        write(6,*) 'offline data requested has not been read'
        stop 'error stop can_do_offline: no such data'
        end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine off_debug_var

	use mod_offline
	use basin
	use shympi

	implicit none

	integer nelg,nkng
	integer ie_int,ie_ext
	integer ik_int,ik_ext
	integer i,ie,ik,iu,lmax,it
	double precision dtime
	character*80 name

	integer, save, allocatable :: ies(:)
	integer, save, allocatable :: iks(:)
	integer, save, allocatable :: iue(:)
	integer, save, allocatable :: iuk(:)
	integer, save :: nie = 0
	integer, save :: nik = 0
	integer, save :: icall_local = 0

	integer ieint, ipint
	integer ieext, ipext

	nelg = nel_global
	nkng = nkn_global

        call get_act_dtime(dtime)
        it = nint(dtime)

	if( icall_local == 0 ) then

	  allocate(ies(nel))
	  allocate(iks(nkn))
	  allocate(iue(nel))
	  allocate(iuk(nkn))

	  iu = 300

	  do ie_ext=1,nelg,nelg/20
	    ie_int = ieint(ie_ext)
	    if( ie_int > 0 .and. shympi_is_unique_elem(ie_int) ) then
	      nie = nie + 1
	      ies(nie) = ie_int
	      iu = iu + 1
	      iue(nie) = iu
	      call make_name_with_number('elem',ie_ext,'aux',name)
	      open(iu,file=name,status='unknown',form='formatted')
	    end if
	  end do

	  do ik_ext=1,nkng,nkng/20
	    ik_int = ipint(ik_ext)
	    if( ik_int > 0 .and. shympi_is_unique_node(ik_int) ) then
	      nik = nik + 1
	      iks(nik) = ik_int
	      iu = iu + 1
	      iuk(nik) = iu
	      call make_name_with_number('node',ik_ext,'aux',name)
	      open(iu,file=name,status='unknown',form='formatted')
	    end if
	  end do

	  icall_local = 1
	end if

	do i=1,nie
	  ie = ies(i)
	  ie_ext = ieext(ie)
	  iu = iue(i)
	  lmax = ile(ie)
	  write(iu,*) it,dtime
	  write(iu,*) ie_ext,lmax
	  write(iu,*) ze(:,ie,:)
	  write(iu,*) ut(1:lmax,ie,:)
	  write(iu,*) vt(1:lmax,ie,:)
	end do

	do i=1,nik
	  ik = iks(i)
	  ik_ext = ipext(ik)
	  iu = iuk(i)
	  lmax = ilk(ik)
	  write(iu,*) it,dtime
	  write(iu,*) ik_ext,lmax
	  write(iu,*) zn(ik,:)
	  write(iu,*) wn(1:lmax,ik,:)
	  write(iu,*) sn(1:lmax,ik,:)
	  write(iu,*) tn(1:lmax,ik,:)
	  write(iu,*) vd(1:lmax,ie,:)
	  write(iu,*) dd(1:lmax,ie,:)
	end do

	end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine off_error(ierr,it,text)

! exits with error

	implicit none

	integer ierr,it
	character*(*) text

	if( ierr == 0 ) return 

	write(6,*) '*** error in module offline'
	write(6,*) 'ierr = ',ierr
	write(6,*) 'it = ',it
	write(6,*) trim(text)
	stop 'error stop off_error: module offline'

	end

!****************************************************************

	subroutine off_set_default_ts(t,s)

	use mod_offline

	real t,s

	idef = 1
	tdef = t
	sdef = s

	end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine off_visv_d_debug(it,nlvddi,nkn,a3d)

	use shympi

	implicit none

	integer it
	integer nlvddi,nkn
	double precision a3d(nlvddi,nkn)

	logical bdebug
	integer ik_ext, ik_int, iudbg
	integer ipint

	ik_ext = 1152
	ik_int = ipint(ik_ext)
	bdebug = ik_int > 0
	iudbg = 750 + my_id

	if( .not. bdebug ) return

	write(iudbg,*) 'it: ',it,ik_ext,ik_int,nlvddi
	write(iudbg,*) 'vvv: ',a3d(:,ik_int)

	flush(iudbg)
	
	end

!****************************************************************

	subroutine off_zenv_d_debug(it,nel,zenv)

	use shympi

	implicit none

	integer it
	integer nel
	double precision zenv(3,nel)

	logical bdebug
	integer ie_ext, ie_int, iudbg
	integer ieint

	ie_ext = 4309
	ie_int = ieint(ie_ext)
	bdebug = ie_int > 0
	iudbg = 730 + my_id

	if( .not. bdebug ) return

	write(iudbg,*) 'it: ',it,ie_ext,ie_int
	write(iudbg,*) 'zzz: ',zenv(:,ie_int)

	flush(iudbg)
	
	end

!****************************************************************

	subroutine off_zenv_r_debug(it,nel,zenv)

	use shympi

	implicit none

	integer it
	integer nel
	real zenv(3,nel)

	logical bdebug
	integer ie_ext, ie_int, iudbg
	integer ieint

	ie_ext = 4309
	ie_int = ieint(ie_ext)
	bdebug = ie_int > 0
	iudbg = 740 + my_id

	if( .not. bdebug ) return

	write(iudbg,*) 'it: ',it,ie_ext,ie_int
	write(iudbg,*) 'zzz: ',zenv(:,ie_int)

	flush(iudbg)
	
	end

!****************************************************************
!****************************************************************
!****************************************************************
!
!	program off_main
!	call off_test
!	end
!
!****************************************************************

