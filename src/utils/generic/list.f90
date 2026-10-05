
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017,2019  Georg Umgiesser
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

! list utility routines
!
! revision log :
!
! 05.10.2024	ggu	file copied from stack

!===============================================================
	module list
!===============================================================

	implicit none

	private

        type :: entry
          integer :: top
          integer :: max
          integer :: type
          double precision, allocatable :: array(:)
          character*80, allocatable :: string(:)
        end type entry

        integer, parameter ::     no_type = 0
        integer, parameter ::  value_type = 1
        integer, parameter :: string_type = 2

        integer, parameter :: empty_error = 1
        integer, parameter ::  type_error = 2
        integer, parameter ::  size_error = 3

        integer, save :: idlast = 0
        integer, save :: ndim = 0
        integer, parameter :: ndim_first = 10
        type(entry), save, allocatable :: pentry(:)

	public :: list_init		!call list_init(id)
	public :: list_delete		!call list_delete(id)
	public :: list_add		!call list_add(id,value)
	public :: list_get		!logical list_get(id,value)
	public :: list_remove		!logical list_remove(id,value)
	public :: list_get_entries	!call list_get_entries(id,n,values)
	public :: list_has_entry	!logical list_has_entry(id,value)
	public :: list_fill		!integer list_fill(id)
	public :: list_clear		!call list_clear(id)
	public :: list_is_empty		!logical list_is_empty(id)
	public :: list_info		!call list_info(id)

        INTERFACE list_add
        MODULE PROCEDURE         list_add_d &
     &                          ,list_add_r &
     &                          ,list_add_i
        END INTERFACE

        INTERFACE list_get
        MODULE PROCEDURE         list_get_d &
     &                          ,list_get_r &
     &                          ,list_get_i
        END INTERFACE

        INTERFACE list_remove
        MODULE PROCEDURE         list_remove_d &
     &                          ,list_remove_r &
     &                          ,list_remove_i
        END INTERFACE

        INTERFACE list_get_entries
        MODULE PROCEDURE         list_get_entries_d &
     &                          ,list_get_entries_r &
     &                          ,list_get_entries_i
        END INTERFACE

        INTERFACE list_has_entry
        MODULE PROCEDURE         list_has_entry_d &
     &                          ,list_has_entry_r &
     &                          ,list_has_entry_i
        END INTERFACE

!===============================================================
	contains
!===============================================================

        subroutine list_init_alloc

        type(entry), allocatable :: paux(:)

        if( ndim == 0 ) then
          ndim = ndim_first
          allocate(pentry(ndim))
        else
          ndim = ndim*2
          allocate(paux(ndim))
          paux(1:ndim/2) = pentry(1:ndim/2)
          call move_alloc(paux,pentry)
        end if

        end subroutine list_init_alloc

!******************************************************************

        subroutine list_init_new_id(id)

        integer id

        idlast = idlast + 1
        if( idlast > ndim ) then
          call list_init_alloc
        end if
        id = idlast

        call list_init_id(id)

        end subroutine list_init_new_id

!******************************************************************

        subroutine list_init_id(id)

        integer id

        if( id > ndim ) then
          stop 'error stop list_init_id: ndim'
        end if

        pentry(id)%top = 0
        pentry(id)%max = 0
        pentry(id)%type = 0

	if( allocated(pentry(id)%array) ) deallocate(pentry(id)%array)
	if( allocated(pentry(id)%string) ) deallocate(pentry(id)%string)

        end subroutine list_init_id

!******************************************************************

	subroutine list_error(id,error)

	integer id,error

	if( error == empty_error ) then
	  write(6,*) 'list: ',id
	  stop 'error stop list: list is empty'
	else if( error == type_error ) then
	  write(6,*) 'list: ',id
	  write(6,*) 'type: ',pentry(id)%type
	  stop 'error stop list: variable is of wrong type'
	else if( error == size_error ) then
	  write(6,*) 'list: ',id
	  write(6,*) 'nfill: ',pentry(id)%top
	  stop 'error stop list: list holds more entries than array size'
	else
	  stop 'error stop list: internal error (1)'
	end if

	end subroutine list_error

!******************************************************************

	subroutine realloc_double(n,value)

	integer n
	double precision, allocatable :: value(:)

	integer nsize
	double precision, allocatable :: daux(:)

	if( n == 0 ) then
	  n = 10
	  allocate(value(n))
	else
	  nsize = min(n,size(value))
          allocate(daux(n))
          daux(1:nsize) = value(1:nsize)
          call move_alloc(daux,value)
	end if

	end subroutine realloc_double

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine list_init(id)
	integer id
        call list_init_new_id(id)
	end subroutine list_init

	subroutine list_delete(id)
	integer id
        call list_init_id(id)
	if( id == idlast ) idlast = idlast - 1
	end subroutine list_delete

!--------------------

	subroutine list_add_i(id,value)
	integer id
	integer value
	call list_add_d(id,dble(value))
	end subroutine list_add_i

	subroutine list_add_r(id,value)
	integer id
	real value
	call list_add_d(id,dble(value))
	end subroutine list_add_r

	subroutine list_add_d(id,value)
	integer id
	double precision value
	integer n
	if( pentry(id)%top >= pentry(id)%max ) then
	  n = 2 * pentry(id)%max
	  call realloc_double(n,pentry(id)%array)
	  pentry(id)%max = n
	end if
	if( pentry(id)%type == no_type ) then
	  pentry(id)%type = value_type
	end if
	if( pentry(id)%type /= value_type ) then
	  call list_error(id,type_error)
	end if
	pentry(id)%top = pentry(id)%top + 1
	pentry(id)%array(pentry(id)%top) = value
	end subroutine list_add_d

!--------------------

	logical function list_get_i(id,value)
	integer id
	integer value
	double precision dvalue
	list_get_i = list_get_d(id,dvalue)
	end function list_get_i

	logical function list_get_r(id,value)
	integer id
	real value
	double precision dvalue
	list_get_r = list_get_d(id,dvalue)
	end function list_get_r

	logical function list_get_d(id,value)
	integer id
	double precision value
	integer n
	list_get_d = .false.
	value = 0.
	n = pentry(id)%top
	if( n == 0 ) return
	list_get_d = .true.
	value = pentry(id)%array(n)
	pentry(id)%top = pentry(id)%top - 1
	end function list_get_d

!--------------------

	logical function list_remove_i(id,value)
	integer id
	integer value
	double precision dvalue
	dvalue = dble(value)
	list_remove_i = list_remove_d(id,dvalue)
	end function list_remove_i

	logical function list_remove_r(id,value)
	integer id
	real value
	double precision dvalue
	dvalue = dble(value)
	list_remove_r = list_remove_d(id,dvalue)
	end function list_remove_r

	logical function list_remove_d(id,value)
	integer id
	double precision value
	integer i,n
	list_remove_d = .false.
	n = pentry(id)%top
	if( n == 0 ) return
	i = findloc(pentry(id)%array(1:n),value,1)
	!write(6,*) 'i: ',i
	if( i == 0 ) return
	list_remove_d = .true.
	pentry(id)%array(i) = pentry(id)%array(pentry(id)%top)
	pentry(id)%top = pentry(id)%top - 1
	end function list_remove_d

!--------------------

	subroutine list_get_entries_i(id,n,values)
	integer id
	integer n
	integer values(:)
	integer nsize
	double precision, allocatable :: dvalues(:)
	nsize = size(values)
	n = pentry(id)%top
	if( n > nsize ) call list_error(id,size_error)
	allocate(dvalues(n))
	call list_get_entries_d(id,n,dvalues)
	values(1:n) = nint(dvalues(1:n))
	end subroutine list_get_entries_i

	subroutine list_get_entries_r(id,n,values)
	integer id
	integer n
	real values(:)
	integer nsize
	double precision, allocatable :: dvalues(:)
	nsize = size(values)
	n = pentry(id)%top
	if( n > nsize ) call list_error(id,size_error)
	allocate(dvalues(n))
	call list_get_entries_d(id,n,dvalues)
	values(1:n) = real(dvalues(1:n))
	end subroutine list_get_entries_r

	subroutine list_get_entries_d(id,n,values)
	integer id
	integer n
	double precision values(:)
	integer nsize
	nsize = size(values)
	n = pentry(id)%top
	if( n > nsize ) call list_error(id,size_error)
	values(1:n) = pentry(id)%array(1:n)
	end subroutine list_get_entries_d

!--------------------

	logical function list_has_entry_i(id,value)
	integer id
	integer value
	list_has_entry_i = list_has_entry_d(id,dble(value))
	end function list_has_entry_i

	logical function list_has_entry_r(id,value)
	integer id
	real value
	list_has_entry_r = list_has_entry_d(id,dble(value))
	end function list_has_entry_r

	logical function list_has_entry_d(id,value)
	integer id
	double precision value
	integer n
	list_has_entry_d = .false.
	n = pentry(id)%top
	if( n == 0 ) return
	list_has_entry_d = any( pentry(id)%array(1:n) == value )
	end function list_has_entry_d

!--------------------

	integer function list_fill(id)
	integer id
	list_fill = pentry(id)%top
	end function list_fill

!--------------------

	subroutine list_clear(id)
	integer id
	pentry(id)%top = 0
	end subroutine list_clear

!--------------------

	logical function list_is_empty(id)
	integer id
	list_is_empty = ( pentry(id)%top == 0 )
	end function list_is_empty

!--------------------

	subroutine list_info(id)
	integer id
	write(6,*) 'list_info: ',id,pentry(id)%top &
     &			,pentry(id)%max,pentry(id)%type
	end subroutine list_info

!===============================================================
	end module list
!===============================================================

	subroutine list_test

	use list

	implicit none

	integer, parameter :: nloop = 500
	integer, save :: ndim = 1000
	integer, allocatable :: vals(:)
	integer val,value,nl,n,i,id,ind,nop,nfill,ne,nmax,perc
	logical bdebug,bverb,bwrite
	real r

	bdebug = .true.
	bdebug = .false.
	bverb = .true.
	bverb = .false.
	bwrite = .false.
	bwrite = .true.

	call list_init(id)
	allocate(vals(ndim))

	call random_seed
	val = 0
	ind = 0
	nop = 0
	nmax = 0

	do nl=1,nloop
	  call list_rand_int(1,10,n)	!add 1-10 new values
	  if( bdebug ) write(6,*) 'add a total of values: ',n
	  do i = 1,n
	    val = val + 1
	    if( bverb ) write(6,*) 'add: ',val
	    call list_add(id,val)
	    ind = ind + 1
	  end do
	  nmax = nmax + n
	  nop = nop + n
	  n = nmax
	  call list_rand_int(1,val,n)
	  if( bdebug ) write(6,*) 'remove a total of values: ',n
	  do i = 1,n
	    call list_rand_int(1,15,value)
	    if( bverb ) write(6,*) 'remove value: ',value
	    if( list_remove(id,value) ) then
	      if( bverb ) write(6,*) 'value removed: ',value
	      ind = ind - 1
	      nop = nop + 1
	    else
	      if( bverb ) write(6,*) 'value not removed: ',value
	      if( list_has_entry(id,value) ) then
		write(6,*) '*** value not removed: ',value
		stop 'error stop list_test: did not remove value'
	      end if
	    end if
	    nfill = list_fill(id)
	    if( nfill > ndim ) then
	      deallocate(vals)
	      ndim = max(ndim*2,nfill)
	      allocate(vals(ndim))
	    end if
	    call list_get_entries(id,ne,vals)
	    if( bverb ) write(6,'(20i5)') vals(1:ne)
	  end do
	  nfill = list_fill(id)
	  if( nfill /= ind ) then
	    stop 'error stop list_test: nfill/=ind'
	  end if
	  if( bwrite ) then
	    if( mod(nl,nloop/100) == 0 ) then
		perc = (100*nl)/nloop
		write(6,*) 'loop: ',perc,nl,ind,nfill
	    end if
	  end if
	end do

	call list_delete(id)

	write(6,*) 'list test successfully finished: ',nloop,nop,val

	end

!******************************************************************

        subroutine list_assert(bcheck,text,id)
        use list
        implicit none
        logical bcheck
        character*(*) text
        integer id
        if( .not. bcheck ) then
          write(6,*) 'list_assertion: ',trim(text)
          call list_info(id)
          stop 'error stop list_assert: assertion failed'
        end if
        end

	subroutine list_rand_int(min,max,irand)

	implicit none
	integer min,max
	integer irand
	real r

	call random_number(r)
	irand = min + (max-min+1)*r

	end

!******************************************************************

	program list_main
	call list_test
	end program list_main

!******************************************************************
