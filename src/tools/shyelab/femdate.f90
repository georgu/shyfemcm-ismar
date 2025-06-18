
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015-2019  Georg Umgiesser
!    Copyright (C) 2017  Marco Bajo
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

! elaborates fem files
!
! revision log :
!
! 14.01.2015	ggu	adapted from feminf
! 20.05.2015	ggu	use bhuman to convert to human readable time
! 05.06.2015	ggu	iextract to extract nodal value
! 05.11.2015	ggu	new option chform to change format
! 04.10.2016	ggu	output flags now similar to shyelab
! 05.10.2016	ggu	allow for expansion of regular grid
! 11.10.2016	ggu	introduced flag for min/max/med computation
! 31.10.2016	ggu	new flag condense (bcondense)
! 16.05.2017	ggu&mbj	better handling of points to extract
! 30.01.2018	ggu	written with new fem_util module
! 22.02.2018	ggu	changed VERS_7_5_42
! 16.02.2019	ggu	changed VERS_7_5_60
! 27.01.2022	ggu	minor changes
! 16.03.2022	ggu	femadd newly written
!
!******************************************************************

	program femdate

! changes date and time in fem file

	use clo
	use fem_util

	implicit none

	character*80 name,string,infile
	integer nfile,i,ierr,iformat
	integer nvar,nvar0,nrecs
	integer ianz
	double precision atime,atime0
	logical bdebug
	logical bextend
	logical bverb,bquiet,bsilent
	logical bunform
	character*20 aline
	character*80 sextend,s(2)
	character*80 stime
	double precision astart,aend
	type(femfile_type), allocatable :: ffinfo(:)
	type(femfile_type) :: ffiout
	type(femrec_type), allocatable :: finfo(:)
	type(femrec_type) :: fout
	type(femrec_type) :: fextra

	integer iscans

	bdebug = .true.
	bdebug = .false.

!--------------------------------------------------------------
! set command line options
!--------------------------------------------------------------

	call clo_init('femadd','fem-files','1.0')

	call clo_add_info('adds vars of multiple fem-files into one')

        call clo_add_sep('output options')

        call clo_add_option('verb',.false.,'be more verbose')
        call clo_add_option('silent',.false.,'be silent')
        call clo_add_option('quiet',.false.,'do not write time records')

        call clo_add_sep('action')

	call clo_add_option('date0 time',' ','change time in records')

        call clo_add_extra('time is YYYY-MM-DD[::hh:mm:ss]')

!--------------------------------------------------------------
! parse command line options
!--------------------------------------------------------------

	call clo_parse_options(1)  !expecting (at least) 1 file after options

!--------------------------------------------------------------
! get command line options
!--------------------------------------------------------------

	call clo_get_option('verb',bverb)
	call clo_get_option('quiet',bquiet)
	call clo_get_option('silent',bsilent)

	if( bsilent ) bquiet = .true.

	call clo_get_option('date0',stime)

	if( stime == ' ' ) then
	  stop 'error stop femdate: no date/time given'
	end if

!--------------------------------------------------------------
! set parameters
!--------------------------------------------------------------

	nfile = clo_number_of_files()

!--------------------------------------------------------------
! open all files
!--------------------------------------------------------------

	if( nfile < 1 ) then
	  write(6,*) 'No file given... exiting'
	  stop 'error stop femdate: no files'
	else if( nfile > 1 ) then
	  write(6,*) 'Can only handle one file at at time... exiting'
	  stop 'error stop femdate: too many files'
	end if

	allocate(ffinfo(nfile))
	allocate(finfo(nfile))

	call femutil_init_record(finfo(1))
        call clo_get_file(1,infile)
	call femutil_open_for_read(infile,0,ffinfo(1),ierr)
	if( ierr /= 0 ) goto 99

	iformat = 1
	call femutil_open_for_write('out.fem',iformat,ffiout)

!--------------------------------------------------------------
! loop on files and read data
!--------------------------------------------------------------

	nrecs = 0
	nvar0 = 0
	nvar = 0

	do

	  i = 1
	  call femutil_read_record(ffinfo(i),finfo(i),ierr)
	  if( ierr < 0 ) exit
	  if( ierr /= 0 ) goto 98
	  call femutil_get_time(finfo(i),atime)

	  nvar = finfo(i)%nvar
	  if( nvar0 == 0 ) nvar0 = nvar
	  if( nvar /= nvar0 ) goto 95

	  nrecs = nrecs + 1

	  fout = finfo(i)
	  call time_change(stime,fout)
	  call femutil_get_time(fout,atime)

          call dts_format_abs_time(atime,aline)
	  if( .not. bquiet ) write(6,*) atime,'  ',aline
	  call femutil_write_record(ffiout,fout)

	  if( nrecs == 5 ) exit
	end do

	if( .not. bsilent ) then
	  write(6,*) 'total number of records treated: ',nrecs
	  write(6,*) 'output written to file out.fem'
	end if

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	stop
   95	continue
	write(6,*) 'nvar is changing: ',nvar,nvar0
	stop 'error stop femdate: nvar not constant'
   96	continue
	write(6,*) 'files are not compatible'
	stop 'error stop femdate: not compatible'
   97	continue
	write(6,*) 'times are not compatible: ',atime0,atime
	stop 'error stop femdate: time error'
   98	continue
	write(6,*) 'error reading record ',i,ierr
	stop 'error stop femdate: read error'
   99	continue
	write(6,*) 'error opening file ',infile
	stop 'error stop femdate: opening error'
        end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine time_change(stime,fin)

	use fem_util

	implicit none

	character*(*) stime
	type(femrec_type) :: fin

	integer ierr
	integer date,time
	double precision atime,dtime,atime0

	!write(6,*) 'datetime: ',fin%datetime
	!write(6,*) 'dtime: ',fin%dtime
	!write(6,*) 'atime: ',fin%atime

	call dts_string2time(stime,atime0,ierr)
	dtime = fin%dtime
	atime = atime0 + dtime
	
	if( ierr /= 0 ) stop 'error stop time_change: error converting time'
	call dts_from_abs_time(date,time,atime)

	fin%datetime = (/date,time/)
	fin%dtime = 0.
	fin%atime = atime

	end

!*****************************************************************

