
!--------------------------------------------------------------------------
!
!    Copyright (C) 2005,2010,2012,2014-2015,2018-2020  Georg Umgiesser
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

! info on restart file

! revision log :
!
! 01.05.2005	ggu	written from scratch
! 23.03.2010	ggu	changed v6.1.1
! 28.09.2010	ggu	changed VERS_6_1_11
! 14.02.2012	ggu	changed VERS_6_1_44
! 29.08.2012	ggu	changed VERS_6_1_56
! 26.11.2014	ggu	changed VERS_7_0_7
! 05.06.2015	ggu	changed VERS_7_1_12
! 05.11.2015	ggu	changed VERS_7_3_12
! 06.07.2018	ggu	changed VERS_7_5_48
! 16.02.2019	ggu	changed VERS_7_5_60
! 03.05.2019	ggu	adapted
! 21.05.2019	ggu	changed VERS_7_5_62
! 09.03.2020	ggu	prepared for mercury restart
! 20.03.2020    ggu     adjusted for new routine calls
! 27.03.2021    ggu     new option -checkval
! 14.04.2021    ggu     bug fix - atime was integer
! 05.10.2025    ggu     new routines to compute difference of rst files

!******************************************************************

!==================================================================
	module mod_rstinf
!==================================================================

	logical, save :: bread = .false.
	logical, save :: bcheckval = .false.
	logical, save :: bquiet = .false.
	logical, save :: bverbose = .false.
	integer, save :: nerror = 0

!==================================================================
	end module mod_rstinf
!==================================================================

	program rstinf

	use clo
	use mod_rstinf

	implicit none

	integer nc,ierr

	call rst_init

        nc = clo_number_of_files()

        if( bverbose ) write(6,*) 'running rstinf...'

        if( nc == 0 ) then
          call clo_usage
        else if( nc == 1 ) then
          call rst_info
        else if( nc == 2 ) then
          call rst_diff(ierr)
        else
          write(6,*) 'nc = ',nc
          stop 'error stop rstinf: wrong number of files'
        end if

        if( ierr > 0 ) then
          if( ierr == 99 ) ierr = 100   !terrible hack - FIXME
          call exit(ierr)
        else
          call exit(99)
        end if

        end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine rst_info

	use clo
	use mod_rstinf

	implicit none

	integer iunit,it,nvers,nrec,nknr,nelr,nlvr,iflag,ierr,ic
	integer nread
	double precision atime
	double precision atime_anf
	double precision atime_end
	character*20 aline
	character*80 file
	character*80 title1
	character*80 title2
	character*80 rstfile

!-------------------------------------------------------------------
! create title strings
!-------------------------------------------------------------------

	title1 = 'version nrec       nkn       nel       nlv' // &
     &                  '     iflag'
!                 12345678901234567890123456789012345678901234567890123
	title2 = '   irec                         atime     date'

!-------------------------------------------------------------------
! initialize and open file
!-------------------------------------------------------------------

	nread = 0
	iunit = 1
	file = ' '

        call clo_get_file(1,rstfile)

	bread = bcheckval

	if( bread ) call open_for_read(rstfile)

	open(iunit,file=rstfile,status='old',form='unformatted')

!-------------------------------------------------------------------
! loop on records
!-------------------------------------------------------------------

	do

	  call rst_skip_record(iunit,atime,nvers,nrec &
     &					,nknr,nelr,nlvr,iflag,ierr)
	  if( ierr .ne. 0 ) exit

	  if( nread == 0 ) then
            write(6,1000) trim(title1)
            write(6,1010) nvers,nrec,nknr,nelr,nlvr,iflag
            write(6,*)
            write(6,1001) trim(title2)
	    atime_anf = atime
	  end if

	  nread = nread + 1
	  call dts_format_abs_time(atime,aline)
	  write(6,1011) nread,atime,aline
	  atime_end = atime

	end do

	if( ierr > 0 ) stop 'error stop rstinf: error reading record'
	if( nread == 0 ) stop 'error stop rstinf: no data in file'

!-------------------------------------------------------------------
! final message
!-------------------------------------------------------------------

        write(6,*)
	write(6,*) 'Number of records read: ',nread
	call dts_format_abs_time(atime_anf,aline)
	write(6,*) 'Initial time in file:   ',atime_anf,aline
	call dts_format_abs_time(atime_end,aline)
	write(6,*) 'Final time in file:     ',atime_end,aline
	write(6,*)
        write(6,1000) trim(title1)
        write(6,1010) nvers,nrec,nknr,nelr,nlvr,iflag
        write(6,*)

	call write_flags(iflag)

!-------------------------------------------------------------------
! end of routine
!-------------------------------------------------------------------

	stop
 1000	format(a)
 1001	format(a)
 1010	format(i7,i5,5i10)
 1011	format(i7,f30.2,5x,a20)
	end

!******************************************************************

	subroutine open_for_read(file)

	use basin
	use levels
	use mod_hydro
	use mod_geom_dynamic
	use mod_ts
	use mod_hydro_vel

	implicit none

	character*(*) file

	integer iunit
	double precision atime
	integer nvers,nrec
	integer nk,ne,nl,iflag,ierr

	iunit = 1
	open(iunit,file=file,status='old',form='unformatted')
	call rst_skip_record(iunit,atime,nvers,nrec &
     &					,nk,ne,nl,iflag,ierr)

	call basin_init(nk,ne)
	call levels_init(nk,ne,nl)
	call mod_hydro_init(nk,ne,nl)
	call mod_geom_dynamic_init(nk,ne)
	call mod_ts_init(nk,nl)
	call mod_hydro_vel_init(nk,ne,nl)

	close(iunit)

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine rst_diff(ierr)

        use clo
        use mod_rstinf

	implicit none

	integer ierr

	integer iunit,it,nvers,nrec,nknr,nelr,nlvr,iflag,ic
	integer nread
	integer iu1,iu2
	double precision atime
	double precision atime_anf
	double precision atime_end
	character*20 aline
	character*80 file
	character*80 rstfile1,rstfile2

!-------------------------------------------------------------------
! initialize and open file
!-------------------------------------------------------------------

	nread = 0
	iunit = 1
	file = ' '
	iu1 = 1
	iu2 = 2

        call clo_check_files(2)
        call clo_get_file(1,rstfile1)
        call clo_get_file(2,rstfile2)

	open(iu1,file=rstfile1,status='old',form='unformatted')
	open(iu2,file=rstfile2,status='old',form='unformatted')

!-------------------------------------------------------------------
! loop on records
!-------------------------------------------------------------------

	do

	  call check_record(iu1,iu2,atime,ierr)
	  if( ierr /= 0 ) exit
	  nread = nread + 1
	  if( nread == 1 ) atime_anf = atime
	  atime_end = atime

	end do

!-------------------------------------------------------------------
! final message
!-------------------------------------------------------------------

	if( ierr > 0 ) then
	  write(6,*) 'read error while reading restart file'
	  stop 'error stop: read error'
	end if

        write(6,*)
	write(6,*) 'Number of records read: ',nread
	write(6,*) 'Differences found: ',nerror
	call dts_format_abs_time(atime_anf,aline)
	write(6,*) 'Initial time in file:   ',atime_anf,aline
	call dts_format_abs_time(atime_end,aline)
	write(6,*) 'Final time in file:     ',atime_end,aline
	write(6,*)

!-------------------------------------------------------------------
! end of routine
!-------------------------------------------------------------------

	stop
 1000	format(a)
 1001	format(a)
 1010	format(i7,i5,5i10)
 1011	format(i7,f30.2,5x,a20)
	end

!******************************************************************

	subroutine check_record(iu1,iu2,atime,ierr)

	implicit none

	integer iu1,iu2,ierr
	double precision atime

	integer nvers
	integer nkn,nel,nlv
	integer ival,i
	logical b3d

	call compare_header(iu1,iu2,nvers,atime,nkn,nel,nlv,ierr)
	if( ierr /= 0 ) return
	b3d = nlv > 1
	write(6,*) 'time = ',atime,nkn,nel,nlv

	if( b3d ) then
	  call check_real('hlv',iu1,iu2,1,nlv)
	  call check_integer('ilhv',iu1,iu2,1,nel)
	  call check_integer('ilhkv',iu1,iu2,1,nkn)
	end if

	call check_integer('iwegv',iu1,iu2,1,nel)
	call check_real('znv',iu1,iu2,1,nkn)
	call check_real('zenv',iu1,iu2,3,nel)
	call check_real('utlnv',iu1,iu2,nlv,nel)
	call check_real('vtlnv',iu1,iu2,nlv,nel)

	call check_real('hm3v',iu1,iu2,3,nel)

	call get_integer('ibarcl',iu1,iu2,ival)
	if( ival > 0 ) then
	  call check_real('saltv',iu1,iu2,nlv,nkn)
	  call check_real('tempv',iu1,iu2,nlv,nkn)
	  call check_real('rhov',iu1,iu2,nlv,nkn)
	end if

	call get_integer('iturb',iu1,iu2,ival)
	if( ival > 0 ) then
	  call check_double('numv_gotm',iu1,iu2,(nlv+1),nkn)
	  call check_double('nuhv_gotm',iu1,iu2,(nlv+1),nkn)
	  call check_double('tken_gotm',iu1,iu2,(nlv+1),nkn)
	  call check_double('eps_gotm',iu1,iu2,(nlv+1),nkn)
	  call check_double('rls_gotm',iu1,iu2,(nlv+1),nkn)
	end if

	call get_integer('iconz',iu1,iu2,ival)
	if( ival > 0 ) then
	  do i=1,ival
	    call check_real('conz',iu1,iu2,nlv,nkn)
	  end do
	end if

	call get_integer('nlv-1',iu1,iu2,ival)
	if( ival > 0 ) then
	  call check_real('wlnv',iu1,iu2,(nlv+1),nkn)
	end if

	call get_integer('ieco',iu1,iu2,ival)
	if( ival > 0 ) then
	  stop 'error stop: ieco not ready'
	end if

	call get_integer('imerc',iu1,iu2,ival)
	if( ival > 0 ) then
	  stop 'error stop: imerc not ready'
	end if

	call get_integer('ibfm',iu1,iu2,ival)
	if( ival > 0 ) then
	  stop 'error stop: ibfm not ready'
	end if

	end

!******************************************************************

	subroutine get_integer(text,iu1,iu2,ival)

	implicit none

	character*(*) text
	integer iu1,iu2,ival

	integer ival1,ival2

	read(iu1) ival1
	read(iu2) ival2

	if( ival1 /= ival2 ) then
	  write(6,*) 'values are different: ',trim(text)
	  stop 'error stop'
	end if

	ival = ival1

	end

!******************************************************************

	subroutine check_real(text,iu1,iu2,nv,nh)

	use mod_rstinf

	implicit none

	character*(*) text
	integer iu1,iu2,nv,nh

	integer n,ic
	real, allocatable :: a1(:)
	real, allocatable :: a2(:)

	n = nv*nh
	allocate(a1(n),a2(n))

	read(iu1) a1
	read(iu2) a2

	if( any(a1/=a2) ) then
	  ic = count(a1/=a2)
	  write(6,*) 'records are different: ',trim(text),ic,n
	  nerror = nerror + 1
	end if

	end

!******************************************************************

	subroutine check_double(text,iu1,iu2,nv,nh)

	use mod_rstinf

	implicit none

	character*(*) text
	integer iu1,iu2,nv,nh

	integer n,ic
	double precision, allocatable :: a1(:)
	double precision, allocatable :: a2(:)

	n = nv*nh
	allocate(a1(n),a2(n))

	read(iu1) a1
	read(iu2) a2

	if( any(a1/=a2) ) then
	  ic = count(a1/=a2)
	  write(6,*) 'records are different: ',trim(text),ic,n
	  if( bverbose ) call show_double(n,nv,nh,a1,a2)
	  nerror = nerror + 1
	end if

	end

!******************************************************************

	subroutine check_integer(text,iu1,iu2,nv,nh)

	use mod_rstinf

	implicit none

	character*(*) text
	integer iu1,iu2,nv,nh

	integer n,ic
	integer, allocatable :: a1(:)
	integer, allocatable :: a2(:)

	n = nv*nh
	allocate(a1(n),a2(n))

	read(iu1) a1
	read(iu2) a2

	if( any(a1/=a2) ) then
	  ic = count(a1/=a2)
	  write(6,*) 'records are different: ',trim(text),ic,n
	  nerror = nerror + 1
	end if

	end

!******************************************************************

	subroutine show_double(n,nv,nh,a1,a2)

	implicit none

	integer n,nv,nh
	double precision a1(n),a2(n)

	integer i,iv,ih
	integer ic,imax

	ic = 0
	imax = 20

	do i=1,n
	  if( a1(i) /= a2(i) ) then
	    ic = ic + 1
	    if( ic > imax ) exit
	    iv = 1+mod(i-1,nv)
	    ih = 1+(i-1)/nv
	    write(6,*) i,iv,ih,a1(i),a2(i),abs(a1(i)-a2(i))
	  end if
	end do

	end
	
!******************************************************************

	subroutine read_header(iunit,nvers,date,time,atime,nkn,nel,nlv,ierr)

	use mod_restart

	implicit none

	integer iunit,nvers,date,time,nkn,nel,nlv,ierr
	integer id,ignore
	double precision atime

        read(iunit,iostat=ierr) id,nvers,ignore
	if( ierr /= 0 ) return
        read(iunit) date,time
        read(iunit) atime
        read(iunit) nkn,nel,nlv

	if( id /= idfrst ) stop 'error stop: not a restart file'
	if( ignore /= 1 ) stop 'error stop: ignore /= 1'
	if( nvers < 17 ) stop 'error stop: nvers < 17'

	end

!******************************************************************

	subroutine compare_header(iu1,iu2,nvers,atime,nkn,nel,nlv,ierr)

	use mod_restart

	implicit none

	integer iu1,iu2,nvers,nkn,nel,nlv,ierr
	double precision atime

	integer nvers1,nvers2
	integer date1,time1,nkn1,nel1,nlv1
	integer date2,time2,nkn2,nel2,nlv2
	double precision atime1,atime2

	call read_header(iu1,nvers1,date1,time1,atime1,nkn1,nel1,nlv1,ierr)
	if( ierr /= 0 ) return
	call read_header(iu2,nvers2,date2,time2,atime2,nkn2,nel2,nlv2,ierr)
	if( ierr /= 0 ) return

	if( nvers1 /= nvers2 ) goto 99
	if( date1 /= date2 ) goto 99
	if( time1 /= time2 ) goto 99
	if( atime1 /= atime2 ) goto 99
	if( nkn1 /= nkn2 ) goto 99
	if( nel1 /= nel2 ) goto 99
	if( nlv1 /= nlv2 ) goto 99

	nvers = nvers1
	atime = atime1
	nkn = nkn1
	nel = nel1
	nlv = nlv1

	return
   99	continue
	write(6,*) nvers1,date1,time1,atime1,nkn1,nel1,nlv1
	write(6,*) nvers2,date2,time2,atime2,nkn2,nel2,nlv2
	stop 'error stop: error comparing header'
	end

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine rst_init

        use clo
        use mod_rstinf

        implicit none

        call shyfem_copyright('rstinf - info on restart file')

        call clo_init('rstinf','rstfile [rstfile2]','1.3')

        call clo_add_info('returns info on records of restart file')

        call clo_add_option('checkval',.false.,'check NaNs in file')
        call clo_add_option('quiet',.false.,'be quiet')
        call clo_add_option('verbose',.false.,'be verbose')

        call clo_add_com('if two files given checks for differences')

        call clo_parse_options

        call clo_get_option('checkval',bcheckval)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('verbose',bverbose)
 
        end

!******************************************************************

