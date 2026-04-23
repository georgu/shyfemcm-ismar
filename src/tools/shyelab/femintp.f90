
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

! elaborates and interpolates fem files
!
! revision log :
!
! 23.04.2026	ggu	copied from femadd
!
!******************************************************************

	program femintp

! interpolates to different points values in a fem file (staggered values)

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
	logical bx,by
	character*20 aline
	character*80 sextend,s(2)
	character*80 tstart,tend
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

	call clo_init('femintp','file1.fem file2.fem','1.0')

	call clo_add_info('interpolates points in fem1 to points in fem2')

        call clo_add_sep('output options')

        call clo_add_option('verb',.false.,'be more verbose')
        call clo_add_option('silent',.false.,'be silent')
        call clo_add_option('quiet',.false.,'do not write time records')

	call clo_add_sep('actions')
        call clo_add_option('bintpx',.false.,'interpolate in x')
        call clo_add_option('bintpy',.false.,'interpolate in y')

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

	call clo_get_option('bintpx',bx)
	call clo_get_option('bintpy',by)

!--------------------------------------------------------------
! set parameters
!--------------------------------------------------------------

	nfile = clo_number_of_files()

	if( nfile < 1 ) then
	  write(6,*) 'No file given... exiting'
	  stop 'error stop femadd: no files'
	else if( nfile /= 2 ) then
	  write(6,*) 'Need exactly two files'
	  stop 'error stop femadd: nothing to add'
	end if

	if( bdebug ) then
	  write(6,*) nfile
	  write(6,*) bunform
	end if

!--------------------------------------------------------------
! open all files
!--------------------------------------------------------------

	allocate(ffinfo(nfile))
	allocate(finfo(nfile))

	do i=1,nfile
	  call femutil_init_record(finfo(i))
          call clo_get_file(i,infile)
	  call femutil_open_for_read(infile,0,ffinfo(i),ierr)
	  if( ierr /= 0 ) goto 99
	end do

	iformat = 1
	call femutil_open_for_write('out.fem',iformat,ffiout)

!--------------------------------------------------------------
! loop on files and read data
!--------------------------------------------------------------

	nrecs = 0
	nvar0 = 0

	do

	!------------------------------------------------------
	! read new records from files
	!------------------------------------------------------

	nvar = 0
	do i=1,nfile
	  call femutil_read_record(ffinfo(i),finfo(i),ierr)
	  !if( i == 1 .and. ierr < 0 ) exit
	  if( ierr < 0 ) exit
	  if( ierr /= 0 ) goto 98
	  call femutil_get_time(finfo(i),atime)
	  if( i == 1 ) atime0 = atime
	  if( atime /= atime0 ) goto 97
	  !if( .not. femutil_is_compatible(finfo(1),finfo(i)) ) goto 96
	  nvar = finfo(i)%nvar
	  if( nvar0 == 0 ) nvar0 = nvar
	  if( nvar /= nvar0 ) goto 95
	end do

	if( ierr < 0 ) exit
	nrecs = nrecs + 1

	!------------------------------------------------------
	! interpolate data
	!------------------------------------------------------

	call femutil_interpolate_recs(nfile,finfo,fout,bx,by)

	!------------------------------------------------------
	! write to output file
	!------------------------------------------------------

        call dts_format_abs_time(atime,aline)
	if( .not. bquiet ) write(6,*) atime,'  ',aline
	call femutil_write_record(ffiout,fout)

	end do

!--------------------------------------------------------------
! end of loop on files
!--------------------------------------------------------------

	if( .not. bsilent ) then
	  write(6,*) 'total number of records treated: ',nrecs
	  write(6,*) 'output written to file out.fem'
	end if

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	stop
   95	continue
	write(6,*) '*** nvar is changing: ',nvar,nvar0
	stop 'error stop femadd: nvar not constant'
   96	continue
	write(6,*) '*** files are not compatible'
	stop 'error stop femadd: not compatible'
   97	continue
	write(6,*) '*** times are not compatible: ',atime0,atime
	stop 'error stop femadd: time error'
   98	continue
	write(6,*) '*** error reading record ',i,ierr
	stop 'error stop femadd: read error'
   99	continue
	write(6,*) '*** error opening file ',infile
	stop 'error stop femadd: opening error'
        end

!*****************************************************************
!*****************************************************************
!*****************************************************************


!*****************************************************************

