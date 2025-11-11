
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

! administration of assimilation routines
!
! revision log :
!
! 10.11.2025    ggu     general assimilation of scalar variables
!
!****************************************************************

!=============================================================
	module mod_assimil_admin
!=============================================================

	implicit none

	integer, parameter, private :: ndim = 100
	integer, save, private :: idlast = 0

	type, private :: assim
	  logical :: binit = .false.
	  integer :: num
	  integer :: nobs
	  integer :: ivar
	  integer :: idobs
	  integer, allocatable :: iuse(:)
	  real, allocatable :: corlength(:)
	  real, allocatable :: maxlength(:)
	  real, allocatable :: sigobs(:)
	  real, allocatable :: sigback(:)
	  character*80      :: filecoords
	  character*80      :: fileobs
	  real, allocatable :: xobs(:)
	  real, allocatable :: yobs(:)
	  integer, allocatable :: ies(:)
	end type assim

	type(assim), save, allocatable :: passim(:)

!=============================================================
	contains
!=============================================================

	subroutine assimil_init

	integer, save :: icall = 0

	if( icall /= 0 ) return
	allocate(passim(ndim))
	icall = 1

	end

!*************************************************************

	function is_num_in_list(num)

	logical is_num_in_list
	integer num

	integer i

	is_num_in_list = .true.

	do i=1,idlast
	  if( passim(i)%num == num ) return
	end do

	is_num_in_list = .false.

	end

!*************************************************************

	function is_var_in_list(ivar)

	logical is_var_in_list
	integer ivar

	integer i

	is_var_in_list = .true.

	do i=1,idlast
	  if( passim(i)%ivar == ivar ) return
	end do

	is_var_in_list = .false.

	end

!*************************************************************

	function get_id_from_var(ivar)

	integer get_id_from_var
	integer ivar

	integer i

	get_id_from_var = 0

	do i=1,idlast
	  if( passim(i)%ivar == ivar ) then
	    get_id_from_var = i
	    return
	  end if
	end do

	end

!*************************************************************

	function get_new_id()

	integer get_new_id

	idlast = idlast + 1

	if( idlast > ndim ) then
	  write(6,*) 'too many assimil sections'
	  write(6,*) 'please increase ndim in assimil_admin.f90'
	  stop 'error stop get_new_id: idlast > ndim'
	end if

	passim(idlast)%binit = .false.

	get_new_id = idlast

	end

!=============================================================
	end module mod_assimil_admin
!=============================================================

	subroutine read_assimil(num)

	use mod_assimil_admin

	implicit none

	integer num

	call assimil_init

	if( num <= 0 ) then
	  write(6,*) 'section assimil ',num,' not allowed'
	  write(6,*) 'must specify the number of the nudging section'
	  stop 'error stop read_assimil: num <= 0 or not given'
	end if

	if( is_num_in_list(num) ) then
	  write(6,*) 'section assimil ',num,' is already defined'
	  stop 'error stop read_assimil: section not unique'
	end if

	call assimil_init_section(num)

	end

!*************************************************************

	subroutine check_assimil

	use mod_assimil_admin

	implicit none


	end

!*************************************************************

	subroutine assimil_init_section(num)

	use mod_assimil_admin
	use para
	use nls

	implicit none

	integer num

	integer ivar,nobs
	integer n,id

	real getpar

        call sctpar('assimil')
        call sctfnm('assimil')

	call addpar('ivar',-1.)
	call addpar('nobs',0.)

	call para_add_array_value('corlength',0.)
	call para_add_array_value('maxlength',0.)
	call para_add_array_value('sigobs',0.)
	call para_add_array_value('sigback',0.)
	call para_add_array_value('iuse',1.)

	call addfnm('filecoords',' ')
	call addfnm('fileobs',' ')

	call nls_read_namelist('assimil')		!reads all parameters

	nobs = nint(getpar('nobs'))
	ivar = nint(getpar('ivar'))
	if( ivar < 0 .or. nobs <= 0 ) goto 99
	if( is_var_in_list(ivar) ) goto 98

	id = get_new_id()	

	allocate(passim(id)%corlength(nobs))
	allocate(passim(id)%maxlength(nobs))
	allocate(passim(id)%sigobs(nobs))
	allocate(passim(id)%sigback(nobs))
	allocate(passim(id)%iuse(nobs))

	passim(id)%nobs = nobs
	passim(id)%ivar = ivar

	call para_set_array_value('corlength',nobs,passim(id)%corlength)
	call para_set_array_value('maxlength',nobs,passim(id)%maxlength)
	call para_set_array_value('sigobs',nobs,passim(id)%sigobs)
	call para_set_array_value('sigback',nobs,passim(id)%sigback)
	call para_set_array_value('iuse',nobs,passim(id)%iuse)

	call getfnm('filecoords',passim(id)%filecoords)
	call getfnm('fileobs',passim(id)%fileobs)

	if( any( passim(id)%corlength<=0 ) ) goto 97
	if( any( passim(id)%sigobs<=0 ) ) goto 97
	if( any( passim(id)%sigback<=0 ) ) goto 97
	where( passim(id)%maxlength <= 0. )
	  passim(id)%maxlength = 3.*passim(id)%corlength
	end where
	if( passim(id)%filecoords == ' ' ) goto 96
	if( passim(id)%fileobs == ' ' ) goto 96

	call delete_section('assimil')

	return
   96	continue
	write(6,*) 'both filecoords and fileobs must be given'
	write(6,*) 'filecoords = ',trim(passim(id)%filecoords)
	write(6,*) 'fileobs    = ',trim(passim(id)%fileobs)
	stop 'error stop assimil_init_section: file names missing'
   97	continue
	write(6,*) 'some parameters have value 0 which is not allowed'
	write(6,*) 'corlength: ',passim(id)%corlength
	write(6,*) 'sigobs: ',passim(id)%sigobs
	write(6,*) 'sigback: ',passim(id)%sigback
	stop 'error stop assimil_init_section: error in parameters'
   98	continue
	write(6,*) 'ivar has already be specified: ',ivar
	write(6,*) 'can only specify one section with this ivar'
	stop 'error stop assimil_init_section: ivar not unique'
   99	continue
	write(6,*) 'nobs,ivar: ',nobs,ivar
	stop 'error stop assimil_init_section: nobs or ivar not set'
	end

!*************************************************************

	subroutine assimil_ts0

	implicit none

	end

!*************************************************************

