
!--------------------------------------------------------------------------
!
!    Copyright (C) 2025  Georg Umgiesser
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

! check dbg files
!
! revision log :
!
! 05.10.2025    ggu     last changes

!============================================================
	module mod_check_shympi_check
!============================================================

	logical, save :: bsilent = .false.
	logical, save :: bquiet = .false.
	logical, save :: bverbose = .false.
	logical, save :: bnodiff = .false.
	logical, save :: bstop = .false.

	logical, save :: bminmax = .true.
	logical, save :: bwrite = .false.
	logical, save :: bover = .false.
	real, save :: rthresh = 0.

!============================================================
	contains
!============================================================

!*************************************************************

	subroutine check_shympi_check_sub

! simplified version of check_shympi_debug
!
! file format is different from usual dbg files

	use clo

	implicit none

	logical btwo
	integer ierr,iu,iu2,ios
	integer l,ll,lmax
	integer ll2,lmax2
	integer nrec,icall,isact,nsize
	integer nrec2,icall2,isact2,nsize2
	logical belem,belem2
	integer itime
	integer nc
	double precision dtime,dtime2
	character*80 what,what2
	character*80 name_one,name_two
	integer, allocatable :: ies(:,:)
	real, allocatable :: vals(:)
	integer, allocatable :: ies2(:,:)
	real, allocatable :: vals2(:)

	iu = 1
	iu2 = 2
	btwo = .false.
	call check_shympi_check_init

	nc = clo_number_of_files()
	if( nc == 0 ) stop 'no file given'
	if( nc > 2 ) stop 'error stop: cannot handle more than 2 files'
	call clo_get_file(1,name_one)
	open(iu,file=name_one,status='old',form='unformatted',iostat=ios)
	if( ios /= 0 ) stop 'error stop check500: cannot open file'
	write(6,*) 'file opened: ',trim(name_one)
	if( nc == 2 ) then
	  btwo = .true.
	  call clo_get_file(2,name_two)
	  open(iu2,file=name_two,status='old',form='unformatted',iostat=ios)
	  if( ios /= 0 ) stop 'error stop check500: cannot open file 2'
	  write(6,*) 'file opened: ',trim(name_two)
	end if

	do
	  call read_time_header(iu,dtime,isact,nsize,lmax,belem,what,ierr)
	  if( ierr /= 0 ) exit
	  if( btwo ) then
	    call read_time_header(iu2,dtime2,isact2,nsize2 &
     &					,lmax2,belem2,what2,ierr)
	    if( ierr /= 0 ) exit
	    if( dtime /= dtime2 ) goto 99
	    if( isact /= isact2 ) goto 99
	    if( nsize /= nsize2 ) goto 99
	    if( lmax /= lmax2 ) goto 99
	    if( lmax /= lmax2 ) goto 99
	    if( belem .neqv. belem2 ) goto 99
	    if( what /= what2 ) goto 99
	  end if
	  write(6,*) dtime,isact,nsize,lmax,trim(what)
	  if( .not. btwo ) call write_data_header
	  do l=1,lmax
	    call read_level_header(iu,ll,nsize,nrec,icall)
	    if( l /= ll ) stop 'error stop check500: l/=ll'
	    if( btwo ) then
	      call read_level_header(iu2,ll2,nsize2,nrec2,icall2)
	      if( ll /= ll2 ) goto 98
	      if( nsize /= nsize2 ) goto 98
	      if( nrec /= nrec2 ) goto 98
	      if( icall /= icall2 ) goto 98
	    end if
	    allocate(ies(2,nsize),vals(nsize))
	    call read_data(iu,nsize,ies,vals)
	    if( btwo ) then
	      allocate(ies2(2,nsize),vals2(nsize))
	      call read_data(iu2,nsize,ies2,vals2)
	      call compare_files(l,nsize,ies,vals,ies2,vals2)
	      deallocate(ies2,vals2)
	    else
	      call write_data(l,nsize,nrec,icall,ies,vals)
	    end if
	    deallocate(ies,vals)
	  end do
	end do

	if( ierr > 0 ) stop 'error stop check500: time record read error'

	stop
   98	continue
	write(6,*) 'level records are not comparable:'
	write(6,*) ll,nsize,nrec,icall
	write(6,*) ll2,nsize2,nrec2,icall2
	stop 'error stop 98'
   99	continue
	write(6,*) 'time records are not comparable:'
	write(6,*) dtime,isact,nsize,lmax,belem,trim(what)
	write(6,*) dtime2,isact2,nsize2,lmax2,belem2,trim(what2)
	stop 'error stop 99'
	end

!*************************************************************

	subroutine compare_files(l,nsize,ies,vals,ies2,vals2)

	implicit none

	integer l
	integer nsize
	integer ies(2,nsize)
	real vals(nsize)
	integer ies2(2,nsize)
	real vals2(nsize)

	integer i,ierror
	integer, save :: ierror_max = 10

	ierror = 0

	do i=1,nsize
	  if( any(ies(:,i) /= ies2(:,i) ) ) then
	    ierror = ierror + 1
	    if( ierror > ierror_max ) cycle
	    write(6,*) l,ies(:,i),ies2(:,i)
	  end if
	end do
	if( ierror > 0 ) then
	  write(6,*) ierror,' errors found, max shown ',ierror_max
	  stop 'error stop: error in ies'
	end if

	if( all( vals == vals2 ) ) return

	write(6,*) ' ierror       l       i   i_int   i_ext' // &
     &			'              val1              val2'

	do i=1,nsize
	  if( vals(i) /= vals2(i) ) then
	    ierror = ierror + 1
	    if( ierror > ierror_max ) cycle
	    write(6,1000) ierror,l,i,ies(:,i),vals(i),vals2(i)
	  end if
	end do
	if( ierror > 0 ) then
	  if( ierror > ierror_max ) then
	    write(6,*) '*** errors found: ',ierror,' max shown ',ierror_max
	  else
	    write(6,*) '*** errors found: ',ierror
	  end if
	  stop 'error stop: error in vals'
	end if

 1000	format(5i8,2e18.8)
	end

!*************************************************************

	subroutine read_time_header(iu,dtime,isact,nsize,lmax,belem,what,ierr)

	implicit none

	integer ierr,iu
	integer lmax
	integer isact
	integer nsize
	logical belem
	double precision dtime
	character*80 what

	integer ios

	what = 'unknown'
	read(iu,iostat=ios) dtime,isact,nsize,lmax,belem,what
	!write(6,*) ios,dtime,isact,nsize,lmax,belem,what
	ierr = ios

	end

!*************************************************************

	subroutine read_level_header(iu,ll,nsize,nrec,icall)

	implicit none

	integer iu
	integer ll,nsize,nrec,icall

	read(iu) ll,nsize,nrec,icall

	end

!*************************************************************

	subroutine read_data(iu,nsize,ies,vals)

        implicit none

	integer iu
	integer nsize
	integer ies(2,nsize)
	real vals(nsize)

	integer i

	read(iu) ies(2,:)
	read(iu) vals(:)

	do i=1,nsize
	  ies(1,i) = i
	end do

	end

!*************************************************************

	subroutine write_data_header

	if( bover ) then
	  write(6,*) '      l       i   i_int   i_ext' // &
     &			'            rmin            rmax' // &
     &			'   nover'
	else
	  write(6,*) '      l       i   i_int   i_ext' // &
     &			'            rmin            rmax'
	end if

	end

!*************************************************************

	subroutine write_data(l,nsize,nrec,icall,ies,vals)

	implicit none

	integer l,nsize
	integer nrec,icall
	integer ies(2,nsize)
	real vals(nsize)

	integer i,nover,iover
	real rmin,rmax,r

	rmin = vals(1)
	rmax = vals(1)
	nover = 0
	do i=1,nsize
	  rmin = min(rmin,vals(i))
	  rmax = max(rmax,vals(i))
	  r = vals(i)
	  if( abs(r) > rthresh ) then
	    nover = nover + 1
	  end if
	end do

	if( bminmax ) then
	  if( bover ) then
	    write(6,1000) l,nsize,nrec,icall,rmin,rmax,nover
	  else
	    write(6,1000) l,nsize,nrec,icall,rmin,rmax
	  end if
	else
	  write(6,*) l,nsize,nrec,icall
	end if

	if( .not. bwrite .or. nover == 0 ) return

	write(6,*) '  iover       l       i   i_int   i_ext' // &
     &			'             val'

	iover = 0
	do i=1,nsize
	  r = vals(i)
	  if( abs(r) > rthresh ) then
	    iover = iover + 1
	    write(6,1100) iover,l,i,ies(:,i),r
	  end if
	end do

	return
 1000	format(4i8,2e16.8,i8)
 1100	format(5i8,2e16.8,i8,e18.8)
	end

!*************************************************************

        subroutine check_shympi_check_init

        use clo

        implicit none

        logical baux
        character*80 version

        version = '2.0'

        call clo_init('check_shympi_check','dbg-file(s)',trim(version))

        call clo_add_info('checks files fort.501')

        call clo_add_sep('general options:')
        call clo_add_option('silent',.false.,'be silent')
        call clo_add_option('quiet',.false.,'be quiet')
        call clo_add_option('nodiff',.false.,'do not show differences')
        call clo_add_option('verbose',.false.,'be verbose')
        call clo_add_option('nostop',.false.,'do not stop at error')
        call clo_add_option('write',.false.,'writes values over threshold')
        call clo_add_option('rthresh',0.,'threshold for writing values')

        call clo_parse_options

        call clo_get_option('silent',bsilent)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('nodiff',bnodiff)
        call clo_get_option('verbose',bverbose)
        call clo_get_option('nostop',baux)
        call clo_get_option('write',bwrite)
        call clo_get_option('rthresh',rthresh)

        if( baux ) bstop = .false.
        if( bsilent ) bquiet = .true.
        if( bquiet ) bverbose = .false.

	if( rthresh > 0. ) bover = .true.

	end

!*************************************************************

!============================================================
	end module mod_check_shympi_check
!============================================================

	program check_shympi_check
	use mod_check_shympi_check
	call check_shympi_check_sub
	end

!*************************************************************

