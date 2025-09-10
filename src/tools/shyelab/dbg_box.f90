
!--------------------------------------------------------------------------
!
!    Copyright (C) 2024  Georg Umgiesser
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

! debug of box files
!
! revision log :
!
! 14.12.2024    ggu     written from scratch

!==========================================================================
        module mod_dbg_box
!==========================================================================

        implicit none

        logical, save :: bsilent = .false.      !be silent
        logical, save :: bquiet = .false.       !be quiet
        logical, save :: bverbose = .false.     !be verbose
        logical, save :: bstop = .true.         !stop on error
        logical, save :: bnostop = .false.      !do not stop on differences
        logical, save :: bnodiff = .true.       !do not show differences
        !logical, save :: bsummary = .true.      !only show summary
        !logical, save :: bbalance = .true.      !balance time step
        double precision, save :: maxdiff = 0.  !max difference allowed

	integer, save :: ntime = 0
	real, save :: mdiff = 0.

	integer, parameter :: id_box = 473226	!id of box file
	integer, save :: nvers			!version of box file
	integer, save :: ftype			!type of box file

	real, allocatable, save :: box_vals(:,:)
	real, allocatable, save :: box_vals2(:,:)
	integer, allocatable, save :: isect(:,:)
	integer, allocatable, save :: isect2(:,:)
	real, allocatable, save :: sect_vals(:,:)
	real, allocatable, save :: sect_vals2(:,:)
	integer, allocatable, save :: nblayers(:)
	integer, allocatable, save :: nblayers2(:)
	integer, allocatable, save :: nslayers(:)
	integer, allocatable, save :: nslayers2(:)
	real, allocatable, save :: box_2d_vals(:,:)
	real, allocatable, save :: box_2d_vals2(:,:)
	real, allocatable, save :: box_3d_vals(:,:,:)
	real, allocatable, save :: box_3d_vals2(:,:,:)
	real, allocatable, save :: box_3d_sect(:,:,:)
	real, allocatable, save :: box_3d_sect2(:,:,:)

	INTERFACE
          SUBROUTINE read_box_boxes(iu,nbox,nvars,box_vals)
	    implicit none
	    integer, intent(in) :: iu
	    integer, intent(out) :: nbox
	    integer, intent(out) :: nvars
	    real, allocatable, intent(out) :: box_vals(:,:)
          END SUBROUTINE read_box_boxes      
	END INTERFACE

	INTERFACE
	  SUBROUTINE read_box_sections(iu,nsect,isect,sect_vals)
	    implicit none
	    integer, intent(in) :: iu
	    integer, intent(out) :: nsect
	    integer, allocatable, intent(out) :: isect(:,:)
	    real, allocatable, intent(out) :: sect_vals(:,:)
          END SUBROUTINE read_box_sections      
	END INTERFACE

	INTERFACE
	  subroutine read_box_3d_boxes(iu,nbox,nvars,nvars2d,nblmax,nblayers &
     &			,box_2d_vals,box_3d_vals)
	  implicit none
	  integer, intent(in) :: iu
	  integer, intent(out) :: nbox
	  integer, intent(out) :: nvars
	  integer, intent(out) :: nvars2d
	  integer, intent(out) :: nblmax
	  integer, allocatable, intent(out) :: nblayers(:)
	  real, allocatable, intent(out) :: box_2d_vals(:,:)
	  real, allocatable, intent(out) :: box_3d_vals(:,:,:)
	  end subroutine read_box_3d_boxes
	END INTERFACE

	INTERFACE
	  subroutine read_box_3d_sections(iu,nsect,nslmax,nslayers &
     &			,isect,box_3d_sect)
	  implicit none
	  integer, intent(in) :: iu
	  integer, intent(out) :: nsect
	  integer, intent(out) :: nslmax
	  integer, allocatable, intent(out) :: nslayers(:)
	  integer, allocatable, intent(out) :: isect(:,:)
	  real, allocatable, intent(out) :: box_3d_sect(:,:,:)
	  end subroutine read_box_3d_sections
	END INTERFACE

!==========================================================================
        end module mod_dbg_box
!==========================================================================

        program dbg_box

        use clo
        use mod_dbg_box
        use mod_error_stop

        implicit none

        integer nc,ierr
        character*80 routine

        call dbg_box_init

	ierr = 0

        nc = clo_number_of_files()

        if( nc == 0 ) then
          call clo_usage
        else if( nc == 1 ) then
          call read_dbg_box_file
        else if( nc == 2 ) then
          call compare_dbg_box_files(ierr)
        else
          write(6,*) 'nc = ',nc
	  call error_stop(routine,'wrong number of files')
        end if

        if( ierr > 0 ) call error_stop(77)
	call success_stop

        end

!**************************************************************************

        subroutine read_dbg_box_file

! reads one file and outputs info

        use clo
        use mod_dbg_box
        use mod_error_stop

        implicit none

        integer nc
        integer ftype1
	integer nbox,nvars
        integer nsect
	integer nvars2d
        integer nblmax,nslmax
	integer ierr
	integer iu,iu1
        integer ios
        double precision dtime
        character*20 aline
        character*60 name_one,text
        character*80 title,femver
        character*80 routine

	routine = 'read_dbg_box_file'

!--------------------------------------------------
! open file(s)
!--------------------------------------------------

        call clo_get_file(1,name_one)

	iu1 = 1
	iu = iu1
        open(iu1,file=name_one,status='old',form='formatted',iostat=ios)

        if( ios /= 0 ) then
          write(6,*) 'no such file: ',trim(name_one)
	  call error_stop(routine,'error opening file')
        end if

        if( .not. bquiet ) write(6,*) 'opened file: ',trim(name_one)

	call read_box_general_header(iu,ftype1)
        if( .not. bquiet ) write(6,*) 'file type: ',ftype1

!--------------------------------------------------
! start of time loop
!--------------------------------------------------

	ntime = 0
	nbox = 0
	nsect = 0
	nblmax = 0
	nslmax = 0

	do

	  call read_box_time(iu,dtime,aline,ierr)
	  if( ierr /= 0 ) exit
          if( .not. bquiet ) write(6,*) 'time: ',dtime,aline
	  ntime = ntime + 1

	  if( ftype == 1 ) then					!meteo
	    call read_box_boxes(iu,nbox,nvars,box_vals)
	  else if( ftype == 2 ) then				!2d
	    call read_box_boxes(iu,nbox,nvars,box_vals)
	    call read_box_sections(iu,nsect,isect,sect_vals)
	  else if( ftype == 3 ) then				!3d
	    call read_box_3d_boxes(iu,nbox,nvars,nvars2d,nblmax,nblayers &
     &			,box_2d_vals,box_3d_vals)
	    call read_box_3d_sections(iu,nsect,nslmax,nslayers &
     &			,isect,box_3d_sect)
	  else if( ftype == 4 ) then				!vertical
	    call read_box_3d_boxes(iu,nbox,nvars,nvars2d,nblmax,nblayers &
     &			,box_2d_vals,box_3d_vals)
	  else
	    stop 'cannot yet handle'
	  end if

	end do

!--------------------------------------------------
! end of time loop
!--------------------------------------------------

        if( bverbose ) then
	  write(6,*) 'time records read: ',ntime
	  write(6,*) 'file type = ',ftype
	  write(6,*) 'nbox      =  ',nbox
	  write(6,*) 'nvars     = ',nvars
	  write(6,*) 'nsect     = ',nsect
	  write(6,*) 'nblmax    = ',nblmax
	  write(6,*) 'nslmax    = ',nslmax
	end if

!--------------------------------------------------
! end of routine
!--------------------------------------------------

	return
   97	continue
	write(6,*) 'on unit = ',iu
	write(6,*) 'error reading record: ',ierr
	call error_stop(routine,'error reading record')
   98	continue
	write(6,*) 'on unit = ',iu
	write(6,*) 'error reading header2: ',ierr
	call error_stop(routine,'error reading header2')
   99	continue
	write(6,*) 'on unit = ',iu
	write(6,*) 'error reading header: ',ierr
	call error_stop(routine,'error reading header')
        end

!**************************************************************************

        subroutine compare_dbg_box_files(ierr_recs)

! reads one file and outputs info

        use clo
        use mod_dbg_box
        use mod_error_stop

        implicit none

	integer ierr_recs

        integer nc
        integer nbox,nbox2
        integer nsect,nsect2
        integer ftype1,ftype2
        integer nvars,nvars2
        integer nvars2d,nvars2d2
        integer nblmax,nblmax2
        integer nslmax,nslmax2
	integer iu,iu1,iu2
        integer ios,ierr,is,ib
        double precision dtime,dtime2
        character*20 aline,aline2
        character*60 name,name_one,name_two,text
        character*80 title,title2,femver,header,routine

	routine = 'compare_dbg_flx_files'
	header = ' time                iv    ivar      is       l    lmax'

	ierr_recs = 0

!--------------------------------------------------
! open file(s)
!--------------------------------------------------

        call clo_get_file(1,name_one)

	iu1 = 1
	iu = iu1
	name = name_one
        open(iu,file=name,status='old',form='formatted',iostat=ios)
        if( ios /= 0 ) goto 91
        if( .not. bquiet ) write(6,*) 'opened file: ',trim(name)

        call clo_get_file(2,name_two)

	iu2 = 2
	iu = iu2
	name = name_two
        open(iu,file=name,status='old',form='formatted',iostat=ios)
        if( ios /= 0 ) goto 91
        if( .not. bquiet ) write(6,*) 'opened file: ',trim(name)

!--------------------------------------------------
! read header
!--------------------------------------------------

	call read_box_general_header(iu1,ftype1)
	call read_box_general_header(iu2,ftype2)
	if( ftype1 /= ftype2 ) goto 99
        if( .not. bquiet ) write(6,*) 'file type: ',ftype1

!--------------------------------------------------
! start of time loop
!--------------------------------------------------

	ntime = 0
	nbox = 0
	nsect = 0
	nblmax = 0
	nslmax = 0

	do

	  call read_box_time(iu1,dtime,aline,ierr)
	  if( ierr /= 0 ) exit
	  call read_box_time(iu2,dtime2,aline2,ierr)
	  if( ierr /= 0 ) exit
	  if( dtime /= dtime2 ) goto 98
	  if( aline /= aline2 ) goto 98
          if( .not. bquiet ) write(6,*) 'time: ',dtime,aline
	  ntime = ntime + 1

	  if( ftype == 1 ) then					!meteo

	    call read_box_boxes(iu1,nbox,nvars,box_vals)
	    call read_box_boxes(iu2,nbox2,nvars2,box_vals2)
	    if( nbox /= nbox2 ) goto 97
	    if( nvars /= nvars2 ) goto 97
	    call compare_2d_values('box',nbox,nvars &
     &				,box_vals,box_vals2,ierr_recs)

	  else if( ftype == 2 ) then				!2d

	    call read_box_boxes(iu1,nbox,nvars,box_vals)
	    call read_box_boxes(iu2,nbox2,nvars2,box_vals2)
	    if( nbox /= nbox2 ) goto 97
	    if( nvars /= nvars2 ) goto 97
	    call compare_2d_values('box',nbox,nvars &
     &				,box_vals,box_vals2,ierr_recs)
	    call read_box_sections(iu1,nsect,isect,sect_vals)
	    call read_box_sections(iu2,nsect2,isect2,sect_vals2)
	    if( nsect /= nsect2 ) goto 96
	    if( any(isect /= isect2 ) ) goto 95
	    call compare_2d_values('section',nbox,nvars &
     &				,sect_vals,sect_vals2,ierr_recs)

	  else if( ftype == 3 ) then				!3d

	    call read_box_3d_boxes(iu1,nbox,nvars,nvars2d,nblmax,nblayers &
     &			,box_2d_vals,box_3d_vals)
	    call read_box_3d_boxes(iu2,nbox2,nvars2,nvars2d2,nblmax2,nblayers2 &
     &			,box_2d_vals2,box_3d_vals2)
	    if( nbox /= nbox2 ) goto 97
	    if( nvars /= nvars2 ) goto 97
	    if( nblmax /= nblmax2 ) goto 94
	    if( any(nblayers /= nblayers2 ) ) goto 93
	    call compare_2d_values('box',nbox,3 &
     &				,box_2d_vals,box_2d_vals2,ierr_recs)
	    call compare_3d_values('box',nbox,nblmax,nvars &
     &				,box_3d_vals,box_3d_vals2,ierr_recs)
	    call read_box_3d_sections(iu1,nsect,nslmax,nslayers &
     &			,isect,box_3d_sect)
	    call read_box_3d_sections(iu2,nsect2,nslmax2,nslayers2 &
     &			,isect2,box_3d_sect2)
	    if( nsect /= nsect2 ) goto 96
	    if( nslmax /= nslmax2 ) goto 92
	    if( any(nslayers /= nslayers2 ) ) goto 91
	    if( any(isect /= isect2 ) ) goto 95
	    call compare_3d_values('section',nsect,nslmax,3 &
     &				,box_3d_sect,box_3d_sect2,ierr_recs)

	  else if( ftype == 4 ) then				!3d

	    call read_box_3d_boxes(iu1,nbox,nvars,nvars2d,nblmax,nblayers &
     &			,box_2d_vals,box_3d_vals)
	    call read_box_3d_boxes(iu2,nbox2,nvars2,nvars2d2,nblmax2,nblayers2 &
     &			,box_2d_vals2,box_3d_vals2)
	    if( nbox /= nbox2 ) goto 97
	    if( nvars /= nvars2 ) goto 97
	    if( nblmax /= nblmax2 ) goto 94
	    if( any(nblayers /= nblayers2 ) ) goto 93
	    call compare_3d_values('box',nbox,nblmax,nvars &
     &				,box_3d_vals,box_3d_vals2,ierr_recs)

	  else

	    stop 'cannot yet handle'

	  end if

	end do

!--------------------------------------------------
! end of time loop
!--------------------------------------------------

        if( bverbose ) then
	  write(6,*) 'file type = ',ftype
	  write(6,*) 'nbox      =  ',nbox
	  write(6,*) 'nvars     = ',nvars
	  write(6,*) 'nsect     = ',nsect
	  write(6,*) 'nblmax    = ',nblmax
	  write(6,*) 'nslmax    = ',nslmax
	end if

	call box_finalize(ierr_recs)

!--------------------------------------------------
! end of routine
!--------------------------------------------------

	return
   91	continue
	do is=1,nsect
	  write(6,*) is,nslayers(ib),nslayers2(ib)
	end do
	call error_stop(routine,'incompatible nslayers')
   92	continue
	write(6,*) 'nslmax,nslmax2: ',nslmax,nslmax2
	call error_stop(routine,'incompatible 3d section record')
   93	continue
	do ib=1,nbox
	  write(6,*) ib,nblayers(ib),nblayers2(ib)
	end do
	call error_stop(routine,'incompatible nblayers')
   94	continue
	write(6,*) 'nblmax,nblmax2: ',nblmax,nblmax2
	call error_stop(routine,'incompatible 3d box record')
   95	continue
	do is=1,nsect
	  write(6,*) is,isect(:,is),isect2(:,is)
	end do
	call error_stop(routine,'incompatible isect')
   96	continue
	write(6,*) 'nsect,nsect2: ',nsect,nsect2
	call error_stop(routine,'incompatible section record')
   97	continue
	write(6,*) 'nbox,nbox2: ',nbox,nbox2
	write(6,*) 'nvars,nvars2: ',nvars,nvars2
	call error_stop(routine,'incompatible box record')
   98	continue
	write(6,*) 'dtime,dtime2: ',dtime,dtime2
	write(6,*) 'aline,aline2: ',aline,aline2
	call error_stop(routine,'incompatible time record')
   99	continue
	write(6,*) 'ftype1,ftype2: ',ftype1,ftype2
	call error_stop(routine,'non compatible files')
        end

!**************************************************************************
!**************************************************************************
!**************************************************************************

	subroutine read_box_general_header(iu,ftype_local)

	use mod_dbg_box
	use mod_error_stop

	implicit none

	integer iu
	integer ftype_local

	integer, parameter :: nvers_min = 6
	integer, parameter :: nvers_max = 6

	integer id,ios
	character*80, save :: routine = 'read_box_general_header'

	read(iu,*,iostat=ios)
	if( ios /= 0 ) goto 98
	read(iu,'(2x,3i10)') id,ftype,nvers

	ftype_local = ftype

	if( id /= id_box ) goto 99
	if( nvers < nvers_min ) goto 99
	if( nvers > nvers_max ) goto 99

	if( ftype == 1 ) return
	if( ftype == 2 ) return
	if( ftype == 3 ) return
	if( ftype == 4 ) return
	
	write(6,*) 'ftype = ',ftype
	call error_stop(routine,'ftype not supported')
	
	return
   98	continue
	write(6,*) 'error reading general header: file is empty'
	call error_stop(routine,'error reading general header')
   99	continue
	write(6,*) 'error reading general header:'
	write(6,*) 'id,id_box: ',id,id_box
	write(6,*) 'nvers,nvers_min,nvers_max: ',nvers,nvers_min,nvers_max
	call error_stop(routine,'error reading general header')
	end

!**************************************************************************

	subroutine read_box_time(iu,dtime,aline,ierr)

	!use mod_dbg_box
	use mod_error_stop

	implicit none

	integer iu
	double precision dtime
	character*20 aline
	integer ierr

	integer ios
	character*80, save :: routine = 'read_box_time'

	ierr = -1

	read(iu,*,iostat=ios)

	if( ios < 0 ) return
	if( ios > 0 ) goto 99
	
	read(iu,*) dtime,aline
	ierr = 0

	return
   99	continue
	write(6,*) 'error reading time header:'
	write(6,*) 'ierr: ',ierr
	call error_stop(routine,'error reading time header')
	end

!**************************************************************************

	subroutine read_box_boxes_header(iu,nbox,nvars)

	implicit none

	integer iu
	integer nbox,nvars

	read(iu,*)
	read(iu,*) nbox,nvars

	end

!**************************************************************************

	subroutine read_box_boxes(iu,nbox,nvars,box_vals)

	use mod_error_stop

	implicit none

	integer, intent(in) :: iu
	integer, intent(out) :: nbox
	integer, intent(out) :: nvars
	real, allocatable, intent(out) :: box_vals(:,:)

	integer ib,i
	character*80, save :: routine = 'read_box_boxes'

	read(iu,*)
	read(iu,*) nbox,nvars

	if( .not. allocated(box_vals) ) allocate(box_vals(nvars,nbox))
	
	read(iu,*)
	do ib=1,nbox
	  read(iu,*) i,box_vals(:,ib)
	  if( ib /= i ) goto 99
	end do

	return
   99	continue
	write(6,*) 'ib,i: ',ib,i
	call error_stop(routine,'error reading box data')
	end

!**************************************************************************

	subroutine read_box_sections(iu,nsect,isect,sect_vals)

	use mod_error_stop

	implicit none

	integer, intent(in) :: iu
	integer, intent(out) :: nsect
	integer, allocatable, intent(out) :: isect(:,:)
	real, allocatable, intent(out) :: sect_vals(:,:)

	integer is,i
	character*80, save :: routine = 'read_box_sections'

	read(iu,*)
	read(iu,*) nsect

	if( .not. allocated(isect) ) allocate(isect(2,nsect))
	if( .not. allocated(sect_vals) ) allocate(sect_vals(3,nsect))
	
	read(iu,*)
	do is=1,nsect
	  read(iu,*) i,isect(:,is),sect_vals(:,is)
	  if( is /= i ) goto 99
	end do

	return
   99	continue
	write(6,*) 'is,i: ',is,i
	call error_stop(routine,'error reading section data')
	end

!**************************************************************************

	subroutine read_box_3d_boxes(iu,nbox,nvars,nvars2d,nblmax,nblayers &
     &			,box_2d_vals,box_3d_vals)

	use mod_error_stop

	implicit none

	integer, intent(in) :: iu
	integer, intent(out) :: nbox
	integer, intent(out) :: nvars
	integer, intent(out) :: nvars2d		!how many 2d variables
	integer, intent(out) :: nblmax
	integer, allocatable, intent(out) :: nblayers(:)
	real, allocatable, intent(out) :: box_2d_vals(:,:)
	real, allocatable, intent(out) :: box_3d_vals(:,:,:)

	integer ib,l,i
	integer nvers,ftype
	character*80, save :: routine = 'read_box_3d_boxes'

	read(iu,*)
	read(iu,*) nbox,nblmax,nvars

	call get_box_params(nvers,ftype)
	nvars2d = 3
	if( ftype == 4 ) nvars2d = 0

	if( .not. allocated(nblayers) ) allocate(nblayers(nbox))
	if( .not. allocated(box_2d_vals) ) allocate(box_2d_vals(nvars2d,nbox))
	if( .not. allocated(box_3d_vals) ) then
	  allocate(box_3d_vals(nvars,nblmax,nbox))
	end if

	do ib=1,nbox
	  read(iu,*)
	  read(iu,*) i,nblayers(ib),box_2d_vals(:,ib)
	  if( nblayers(ib) > nblmax ) goto 97
	  if( ib /= i ) goto 99
	  read(iu,*)
	  box_3d_vals(:,:,ib) = 0.
	  do l=1,nblayers(ib)
	    read(iu,*) i,box_3d_vals(:,l,ib)
	    if( l /= i ) goto 98
	  end do
	end do

	return
   97	continue
	write(6,*) 'ib,nblayers(ib),nblmax: ',ib,nblayers(ib),nblmax
	call error_stop(routine,'error reading box nblayers')
   98	continue
	write(6,*) 'l,i: ',l,i
	call error_stop(routine,'error reading box level data')
   99	continue
	write(6,*) 'ib,i: ',ib,i
	call error_stop(routine,'error reading box data')
	end

!**************************************************************************

	subroutine read_box_3d_sections(iu,nsect,nslmax,nslayers &
     &			,isect,box_3d_sect)

	use mod_error_stop

	implicit none

	integer, intent(in) :: iu
	integer, intent(out) :: nsect
	integer, intent(out) :: nslmax
	integer, allocatable, intent(out) :: nslayers(:)
	integer, allocatable, intent(out) :: isect(:,:)
	real, allocatable, intent(out) :: box_3d_sect(:,:,:)

	integer is,l,i
	integer ifrom,ito,nl
	character*80, save :: routine = 'read_box_3d_sections'

	read(iu,*)
	read(iu,*) nsect,nslmax

	if( .not. allocated(nslayers) ) allocate(nslayers(nsect))
	if( .not. allocated(isect) ) allocate(isect(2,nsect))
	if( .not. allocated(box_3d_sect) ) then
	  allocate(box_3d_sect(3,nslmax,nsect))
	end if

	do is=1,nsect
	  read(iu,*)
	  read(iu,*) i,isect(:,is),nslayers(is)
	  if( nslayers(is) > nslmax ) goto 97
	  if( is /= i ) goto 99
	  read(iu,*)
	  box_3d_sect(:,:,is) = 0.
	  do l=1,nslayers(is)
	    read(iu,*) i,box_3d_sect(:,l,is)
	    if( l /= i ) goto 98
	  end do
	end do

	return
   97	continue
	write(6,*) 'is,nslayers(ib),nslmax: ',is,nslayers(is),nslmax
	call error_stop(routine,'error reading section nslayers')
   98	continue
	write(6,*) 'l,i: ',l,i
	call error_stop(routine,'error reading section level data')
   99	continue
	write(6,*) 'is,i: ',is,i
	call error_stop(routine,'error reading section data')
	end

!**************************************************************************

	subroutine compare_2d_values(what,nbox,nvars,vals2d,vals2d2,ierr_recs)

	use mod_dbg_box
	use mod_error_stop

	implicit none

	character*(*) what
	integer nbox,nvars
	real vals2d(nvars,nbox)
	real vals2d2(nvars,nbox)
	integer ierr_recs

	logical bwrite
	integer ib,iv
	real eps
	character*80, save :: routine = 'compare_2d_values'

	eps = maxdiff
	bwrite = .not. bnodiff .and. .not. bquiet

	if( all(abs(vals2d-vals2d2) <= eps ) ) return

	if( bwrite ) write(6,*) 'differences found: ',trim(what)

	do ib=1,nbox
	  do iv=1,nvars
	    if( abs(vals2d(iv,ib)-vals2d2(iv,ib)) > eps ) then
	      mdiff = max(mdiff,abs(vals2d(iv,ib)-vals2d2(iv,ib)))
	      ierr_recs = ierr_recs + 1
	      if( bwrite ) write(6,*) iv,ib,vals2d(iv,ib),vals2d2(iv,ib)
	    end if
	  end do
	end do

	if( ierr_recs > 0 .and. bstop ) call box_finalize(ierr_recs)

	end

!**************************************************************************

	subroutine compare_3d_values(what,nbox,lmax,nvars &
     &				,vals3d,vals3d2,ierr_recs)

	use mod_dbg_box
	use mod_error_stop

	implicit none

	character*(*) what
	integer nbox,lmax,nvars
	real vals3d(nvars,lmax,nbox)
	real vals3d2(nvars,lmax,nbox)
	integer ierr_recs

	logical bwrite
	integer ib,iv,l
	real eps
	character*80, save :: routine = 'compare_3d_values'

	eps = maxdiff
	bwrite = .not. bnodiff .and. .not. bquiet

	if( all(abs(vals3d-vals3d2) <= eps ) ) return

	if( bwrite ) write(6,*) 'differences found: ',trim(what)

	do ib=1,nbox
	 do l=1,lmax
	  do iv=1,nvars
	    if( abs(vals3d(iv,l,ib)-vals3d2(iv,l,ib)) > eps ) then
	      mdiff = max(mdiff,abs(vals3d(iv,l,ib)-vals3d2(iv,l,ib)))
	      ierr_recs = ierr_recs + 1
	      if( bwrite ) write(6,*) iv,l,ib,vals3d(iv,l,ib),vals3d2(iv,l,ib)
	    end if
	  end do
	 end do
	end do

	if( ierr_recs > 0 .and. bstop ) call box_finalize(ierr_recs)

	end

!**************************************************************************

	subroutine box_finalize(ierr_recs)

	use mod_dbg_box
	use mod_error_stop

	implicit none

	integer ierr_recs

	character*80, save :: routine = 'box_finalize'

	if( .not. bsilent ) then
          write(6,*) 'time records read: ',ntime
          if( ierr_recs > 0 ) then
            write(6,*) 'total number of errors: ',ierr_recs
            write(6,*) 'maximum error: ',mdiff
	    call error_stop(routine,'differences found')
          else
            write(6,*) 'no differences found'
	    call success_stop
          end if
	end if

	end

!**************************************************************************

	subroutine get_box_params(nvers_local,ftype_local)

	use mod_dbg_box

	implicit none

	integer nvers_local
	integer ftype_local

	nvers_local = nvers
	ftype_local = ftype

	end

!**************************************************************************
!**************************************************************************
!**************************************************************************

        subroutine dbg_box_init

        use clo
        use mod_dbg_box

        implicit none

        logical baux
        character*80 version

        version = '1.0'

        call clo_init('dbg_box_init','box-file(s)',trim(version))

        call clo_add_info('reads box files')

        call clo_add_sep('general options:')
        call clo_add_option('silent',.false.,'be silent')
        call clo_add_option('quiet',.false.,'be quiet')
        call clo_add_option('verbose',.false.,'be verbose')
        call clo_add_option('nodiff',.false.,'do not show differences')
        call clo_add_option('nostop',.false.,'do not stop at error')
        !call clo_add_option('summary',.false.,'do only summary')
        !call clo_add_option('balance',.false.,'balance time records')
        call clo_add_option('maxdiff',0.,'maximum tolerated difference')

        call clo_parse_options

        call clo_get_option('silent',bsilent)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('verbose',bverbose)
        call clo_get_option('nodiff',bnodiff)
        call clo_get_option('nostop',baux)
        !call clo_get_option('summary',bsummary)
        !call clo_get_option('balance',bbalance)
        call clo_get_option('maxdiff',maxdiff)

        if( baux ) bstop = .false.
        if( bsilent ) bquiet = .true.
        if( bquiet ) bverbose = .false.

        end

!**************************************************************************


