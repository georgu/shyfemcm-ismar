
!--------------------------------------------------------------------------
!
!    Copyright (C) 2018-2019  Georg Umgiesser
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

! indirect heap implementation
!
! the sorting value is in priority
! the index into the (static) array is in irheap
! the index from the (static) array into the heap is in invheap

!**************************************************************

! revision log :
!
! 24.01.2018	ggu	changed VERS_7_5_41
! 14.02.2019	ggu	changed VERS_7_5_56
! 12.09.2026	ggu	implemented indirect heap

! notes :
!
! highest value is at the top (root)
! top has index 1

! Usage:
!
!	iheap_init(nmx)		just allocates
!	iheap_init(n,pa)	initializes and sets priority
!	iheap_init(nmx,n,pa)	initializes and leaves room for addition
!
!	iheap_sort(n,pa)	sorts array pa (does not need call to init)
!	iheap_make		makes heap structure of priority
!
!	iheap_print(text)	prints heap with text
!	iheap_print_tree(text)	prints heap in tree format with text
!	iheap_check_notext	checks heap structure
!	iheap_check_text(text)	checks heap structure
!
!	iheap_has_value()	returns true if heap has still values
!	iheap_replace(ir,p)	replaces priority at ir with p
!	iheap_add(p)		adds p at end of array
!	iheap_remove(ir)	removes array element at ir
!	iheap_pop(p)		pops head of heap and returns p
!	iheap_pop(ir,p)		pops head of heap and returns ir and p
!
!	nmx			maximum size of heap
!	n			size of heap
!	pa(n)			priority array
!	p			priority value		
!	ir			index into array pa
!	text			text string
!
! Typical usage:
!
!	iheap_sort(n,pa)
!
!	iheap_init(nmx)
!	do
!	  iheap_add(p)
!	  ...
!	end do
!	do while( iheap_has_value() )
!	  iheap_pop(p)
!	  ...
!	end do

!**************************************************************

!==============================================================
	module iheap
!==============================================================

	implicit none

	private

	integer, save :: nmax = 0
	integer, save :: nheap = 0
	integer, save :: narray = 0

	integer, allocatable, save :: irheap(:)		!index heap -> array
	integer, allocatable, save :: invheap(:)	!index array -> heap
	real, allocatable, save :: priority(:)

	public &
     &				  iheap_init		&
     &				, iheap_sort		&
     &				, iheap_make		&
     &				, iheap_print		&
     &				, iheap_check		&
     &				, iheap_has_value	&
     &				, iheap_replace		&
     &				, iheap_add		&
     &				, iheap_remove		&
     &				, iheap_pop		&
     &							&
     &				, iheap_print_tree

!     &				, iheap_retire		&
!     &				, iheap_swap		&
!     &				, iheap_adjust		&
!     &				, iheap_promote		&
!     &				, iheap_demote

        INTERFACE iheap_init
        MODULE PROCEDURE          iheap_init_1 &
     &                          , iheap_init_2 &
     &                          , iheap_init_3
        END INTERFACE

        INTERFACE iheap_check
        MODULE PROCEDURE          iheap_check_notext &
     &                          , iheap_check_text
        END INTERFACE

        INTERFACE iheap_pop
        MODULE PROCEDURE          iheap_pop_pir &
     &                          , iheap_pop_p
        END INTERFACE

!==============================================================
	contains
!==============================================================

	subroutine iheap_init_1(nmx)

	integer nmx

	call iheap_allocate(nmx)

	end subroutine iheap_init_1

!*************************************

	subroutine iheap_init_2(n,p)

	integer n
	real p(n)

	call iheap_init_3(n,n,p)

	end subroutine iheap_init_2

!*************************************

	subroutine iheap_init_3(nmx,n,p)

	integer nmx
	integer n
	real p(n)
	
	integer i

	if( nmx < n ) stop 'error stop iheap_init_3: nmx < n'

	call iheap_allocate(nmx)

	nheap = n
	narray = n

	priority(1:n) = p
	do i=1,n
	  irheap(i) = i
	  invheap(i) = i
	end do

	call iheap_make
	
	end subroutine iheap_init_3

!*************************************

	subroutine iheap_allocate(nmx)

	integer nmx

	if( allocated(irheap) ) deallocate(irheap,priority,invheap)

	nmax = nmx
	nheap = 0
	narray = 0

	allocate(irheap(nmax),priority(nmax))
	allocate(invheap(nmax))

	irheap = 0
	invheap = 0
	priority = 0

	end subroutine iheap_allocate

!*************************************

	subroutine iheap_error(text)

	character*(*) text

	stop 'error stop iheap_error: '//trim(text)

	end subroutine iheap_error

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine iheap_sort(n,p)

! sort an array

	implicit none

	integer n
	real p(n)

	integer i

	call iheap_init(n,p)			! also calles iheap_make
	call iheap_check('check in iheap_sort')
	call iheap_retire

	do i=1,n
	  p(i) =  priority(irheap(i))
	end do

	end

!******************************************************************

	subroutine iheap_make

! given an array makes a heap out of it

	implicit none

	integer i

	if( nheap .lt. 2 ) return

	do i=nheap/2,1,-1
	  call iheap_demote(i)
	end do

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine iheap_print(text)

! prints the heap

	implicit none

	character*(*) text

	integer i,ir

	write(6,*) 'heap print: ',trim(text)

	do i=1,nheap
	  ir = irheap(i)
	  write(6,*) i,ir,invheap(ir),priority(ir)
	end do

	end

!******************************************************************

	subroutine iheap_print_tree(text)

! prints the heap in tree format

	implicit none

	character*(*) text

	integer i
	integer ia(nheap)
	real p

	real, parameter :: fact = 1000

	write(6,*) '======================================='
	write(6,*) 'general heap tree print: ',trim(text)
	write(6,*) '======================================='

	ia = nint(fact*priority(1:nheap))
	call iheap_print_tree_internal('priority unsorted',ia)

	do i=1,nheap
	  p = priority(irheap(i))
	  ia(i) = nint(fact*p)
	end do
	call iheap_print_tree_internal('priority sorted',ia)

	ia = irheap(1:nheap)
	call iheap_print_tree_internal('index',ia)

	end
	
!******************************************************************

	subroutine iheap_print_tree_internal(text,ia)

! prints the heap in tree format (internal)

	implicit none

	character*(*) text
	integer ia(nheap)

	integer i,il,n,ns,ne,nt,nspace,nsp
	integer levels
	real rl
	character*80 form

	if( nheap <= 0 ) return

	rl = nheap
	levels = 1 + log(rl)/log(2.)

	write(6,*) 'heap tree print: ',levels,' levels: ',trim(text)

	nspace = 16
	if( nheap > 7 ) nspace = 32
	nspace = 32

	n = 1
	do il=1,nheap
	  nspace = nspace - 4
	  nsp = nspace
	  if( il == 1 ) nsp = nsp - 2
	  ns = n
	  ne = 2*n - 1
	  if( ns > nheap ) exit
	  if( ne > nheap ) ne = nheap
	  write(form,'(a,i2,a)') '(',nsp,'x,10i4)'
	  write(6,form) ia(ns:ne)
	  n = n * 2
	end do

	end

!******************************************************************

	subroutine iheap_check_notext

! checks heap structure

	implicit none

	call iheap_check_text(' ')

	end

!******************************************************************

	subroutine iheap_check_text(text)

! checks heap structure

	implicit none

	character*(*) text

	integer i,j

	do i=1,nheap/2
	  j = i+i
	  if( priority(irheap(i)) .lt. priority(irheap(j)) ) goto 99
	  if( j .eq. nheap ) cycle
	  j = j + 1
	  if( priority(irheap(i)) .lt. priority(irheap(j)) ) goto 99
	end do

	do i=1,nheap
	  if( invheap(irheap(i)) /= i ) goto 98
	  if( irheap(i) == 0 ) goto 97
	  !if( invheap(i) == 0 ) goto 97
	end do

	return
   97	continue
	write(6,*) 'error iheap_check: '//trim(text)
	write(6,*) nheap,i,irheap(i),invheap(i)
	stop 'error stop iheap_check: heap inconsistency: 0 in index'
   98	continue
	write(6,*) 'error iheap_check: '//trim(text)
	write(6,*) nheap,i,irheap(i),invheap(irheap(i))
	stop 'error stop iheap_check: heap inconsistency'
   99	continue
	write(6,*) 'error iheap_check: '//trim(text)
	write(6,*) nheap,i,j,priority(irheap(i)),priority(irheap(j))
	stop 'error stop iheap_check: heap property violated'
	end

!******************************************************************
!******************************************************************
!******************************************************************

	function iheap_has_value()

	implicit none

	logical iheap_has_value

	iheap_has_value = ( nheap > 0 )

	end

!******************************************************************

	subroutine iheap_replace(ir,p)

! replaces value in heap

	implicit none

	integer ir
	real p

	integer i
	real pold

	i = invheap(ir)		!this is position in heap

	pold = priority(ir)
	priority(ir) = p

	call iheap_adjust(i)

	end

!******************************************************************

	subroutine iheap_add(p)

! inserts value at the end and then propmotes it to the right place

	implicit none

	real p

	nheap = nheap + 1
	if( nheap > nmax ) call iheap_error('cannot insert... heap is full')

	irheap(nheap) = nheap
	invheap(nheap) = nheap
	priority(nheap) = p

	call iheap_promote(nheap)

	end

!******************************************************************

	subroutine iheap_remove(ir)

! removes value at index ir from heap

	implicit none

	integer ir

	integer i

	i = invheap(ir)		!this is position in heap
	call iheap_swap(i,nheap)

	irheap(nheap) = 0
	invheap(ir) = 0
	priority(ir) = 0
	nheap = nheap - 1

	call iheap_adjust(i)

	end

!******************************************************************

	subroutine iheap_pop_p(p)

! removes value with highest priority (index 1) from heap and returns it

	implicit none

	real p

	integer ir

	call iheap_pop_pir(ir,p)

	end

!******************************************************************

	subroutine iheap_pop_pir(ir,p)

! removes value with highest priority (index 1) from heap and returns it

	implicit none

	integer ir
	real p

	integer i

	i = 1
	ir = irheap(i)
	p = priority(ir)
	call iheap_swap(i,nheap)

	irheap(nheap) = 0
	invheap(ir) = 0
	priority(ir) = 0
	nheap = nheap - 1

	call iheap_adjust(i)

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine iheap_retire

! used for sorting -> the sorted array (from min to max) is now in the heap
!
! the heap itself is destroyed

	implicit none

	integer i
	integer n

	n = nheap			!save nheap

	do i=nheap,2,-1
	  call iheap_swap(1,i)
	  nheap = nheap - 1
	  call iheap_demote(1)
	end do

	nheap = n

	end

!******************************************************************
!******************************************************************
!******************************************************************
! internal routines ... not to be called from outside
!******************************************************************
!******************************************************************
!******************************************************************

	subroutine iheap_swap(i1,i2)

! swap two values in the heap

	implicit none

	integer i1,i2

	integer ir1,ir2

	if( i1 == i2 ) return

	ir1 = irheap(i1)
	ir2 = irheap(i2)

	irheap(i1) = irheap(i2)
	irheap(i2) = ir1

	invheap(ir1) = i2
	invheap(ir2) = i1

	end

!******************************************************************

	subroutine iheap_adjust(l)

! adjusts entry l to correct place using both promotion and demotion

	implicit none

	integer l

	call iheap_promote(l)
	call iheap_demote(l)

	end

!******************************************************************

	subroutine iheap_promote(l)

! promotes entry l to right place (from below to top)

	implicit none

	integer l

	integer i,j
	integer ir
	real p

	i = l				!adjust this index
	ir = irheap(i)			!this is the original index
	p = priority(irheap(i))		!value to promote
	j = l/2				!node to compare with

	do while( j .ge. 1 )

	  if( p <= priority(irheap(j)) ) exit

	  call iheap_swap(i,j)
	  i = j
	  j = j/2

	end do

	end

!******************************************************************

	subroutine iheap_demote(l)

! demotes entry l to right place (from top downwards)

	implicit none

	integer l

	integer i,j,jlast
	integer ir,ir2
	real p

	i = l				!adjust this index
	ir = irheap(i)			!this is the original index
	p = priority(irheap(i))		!value to promote
	j = l+l				!node to compare with

	do while( j .le. nheap )

	  if( j .lt. nheap ) then
	    if( priority(irheap(j)) .lt. priority(irheap(j+1)) ) then
	      j = j + 1
	    end if
	  end if

	  if( p >= priority(irheap(j)) ) exit

	  call iheap_swap(i,j)
	  i = j
	  j = j + j

	end do

	end

!==============================================================
	end module iheap
!==============================================================

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine iheap_test(ndim)

	use iheap

	implicit none

	!integer, parameter :: ndim = 9
	integer ndim

	real ra(ndim)
	real ra1(ndim)
	real ra3(ndim)
	real raa
	integer i,n,j
	logical berror
	logical, save :: debug = .false.

	j = 0
	n = ndim
	n = 7
	berror = .false.

	write(6,*) 'running basic test'

	do i=1,n
	  call random_number(raa)
	  ra(i) = raa
	end do
	ra1 = ra
	ra3 = ra

	if( debug ) then
	  write(6,*) 'print original'
	  call print_array(n,ra)
	end if

	write(6,*) 'sorting'
	call iheap_sort(n,ra1)
	if( debug ) then
	  call iheap_print('print after sort (should be from low to high)')
	end if
	call check_ordered(n,ra1)

	call iheap_init(n,ra3)
	if( debug ) then
	  call iheap_print('print after init (should have heap structure)')
	end if
	call iheap_check('check after make')

	do i=2,n
	  if( ra1(i) .lt. ra1(i-1) ) then
	    write(6,*) 'not sorted... ',i,ra1(i),ra1(i-1)
	    berror = .true.
	  end if
	end do

	if( berror ) then
	  write(6,*) 'there have been errors...   n = ',n
	else
	  write(6,*) 'test passed.   n = ',n
	end if

	end
	  
!******************************************************************

	subroutine iheap_test_priority_queue(ndim)

	use iheap

	implicit none

	integer ndim

	integer i,ir,ir1,n
	real ra(ndim)
	real raa,p,pp,pold
	logical, save :: debug = .false.

	n = ndim
	write(6,*) 'testing priority queue'

	do i=1,n
	  call random_number(raa)
	  ra(i) = raa
	end do
	pold = maxval(ra(1:n))

	call iheap_init(n,ra)

	i = 0
	do while( iheap_has_value() )
	  i = i + 1
	  call iheap_pop(ir,p)
	  ra(ir) = 0
	  if( ir < n ) then
	    ir1 = ir + 1
	    pp = 0.8*ra(ir1)
	    if( pp > 0 ) then
	      ra(ir1) = pp
	      call iheap_replace(ir1,pp)
	    end if
	  end if
	  if( ir > 1 ) then
	    ir1 = ir - 1
	    pp = 0.8*ra(ir1)
	    if( pp > 0 ) then
	      ra(ir1) = pp
	      call iheap_replace(ir1,pp)
	    end if
	  end if
	  if( debug ) then
	    write(6,*) i,p
	    write(6,*) '===================================='
	    call iheap_print('priority queue')
	    call iheap_check('priority queue')
	    call iheap_print_tree('print_tree')
	    write(6,*) '===================================='
	  end if
	  if( p > pold ) goto 99
	  pold = p
	end do

	write(6,*) 'priority queue test passed n = ',n

	return
   99	continue
	stop 'error stop iheap_test_priority_queue: priorities out of order'
	end

!******************************************************************

	subroutine print_array(n,ra)

	implicit none

	integer n
	real ra(n)

	integer i

	do i=1,n
	  write(6,*) i,ra(i)
	end do

	end

!******************************************************************

	subroutine check_ordered(n,ra)

	implicit none

	integer n
	real ra(n)

	integer i

	do i=2,n
	  if( ra(i) < ra(i-1) ) then
	    write(6,*) 'elements out of order: ',i,ra(i),ra(i-1)
	    stop 'error stop check_ordered: out of order'
	  end if
	end do

	end

!******************************************************************

	subroutine init_random

    IMPLICIT NONE
    INTEGER :: n_seeds
    INTEGER, ALLOCATABLE :: seed_array(:)
    REAL :: rand_val
    INTEGER :: i

    ! 1. Query the size of the seed array required by this compiler
    CALL RANDOM_SEED(SIZE=n_seeds)
    ALLOCATE(seed_array(n_seeds))

    ! 2. Populate the array with your fixed "old" seed pattern
    ! (Avoid setting all zeros or identical tiny values, as some PRNGs misbehave)
    DO i = 1, n_seeds
        seed_array(i) = 12345 + (i * 987)
    END DO

    ! 3. Force the PRNG into this specific starting state
    CALL RANDOM_SEED(PUT=seed_array)

	end

!******************************************************************
!******************************************************************
!******************************************************************
	subroutine iheap_run_test
	!call init_random	!gives reproduceable results
	call iheap_test(9)
	call iheap_test_priority_queue(20)
	end
!******************************************************************
!******************************************************************
!******************************************************************
	program iheap_main
	call iheap_run_test
	end
!******************************************************************

