
!--------------------------------------------------------------------------
!
!    Copyright (C) 1999,2010,2015,2018-2019  Georg Umgiesser
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

! sort routines
!
! contents :
!
! subroutine sortir(ipoint,value,ndim)	indirect sort of real array
! subroutine sortii(ipoint,ivalue,ndim)	indirect sort of interger array
!
! subroutine sort(n,ra)			heapsort, real, direct
! subroutine sort_int(n,ra)		heapsort, integer, direct
! subroutine isort(n,ra,index)		heapsort, integer, indirect
!
! function locate(n,ixx,index,ix)	locates index for ix in array ixx
!
! revision log :
!
! 15.03.1999	ggu	routines written from scratch (Numerical receipes)
! 23.03.2010	ggu	changed v6.1.1
! 10.07.2015	ggu	changed VERS_7_1_50
! 16.12.2015	ggu	changed VERS_7_3_16
! 09.11.2018	ggu	test routines introduced
! 21.12.2018	ggu	changed VERS_7_5_53
! 16.02.2019	ggu	changed VERS_7_5_60
! 02.07.2019	ggu	new routines for indirect real sort, locater_all()
! 05.07.2019	ggu	new routine create_permutation()
!
!***********************************************************
!
! please do not use the first two routines
! they are slow and buggy
!
!***********************************************************

	subroutine sortir(ipoint,value,ndim)

! indirect sort of real array (direct comparison)
!
! ipoint	pointer to values (must be already defined)
! value		values to be sorted
! ndim		total number of values to be sorted

	implicit none

	integer ndim
	integer ipoint(ndim)
	real value(ndim)

	integer i,k,n,n1

	do i=1,ndim-1
	  do k=ndim,i+1,-1
	    n=ipoint(k)
	    n1=ipoint(k-1)
	    if(value(n).lt.value(n1)) then
		ipoint(k)=n1
		ipoint(k-1)=n
	    end if
	  end do
	end do

	end

!***********************************************************

	subroutine sortii(ipoint,ivalue,ndim)
!
! indirect sort of interger array (direct comparison)
!
! ipoint	pointer to values (must be already defined)
! ivalue	values to be sorted
! ndim		total number of values to be sorted
!
	implicit none

	integer ndim
	integer ipoint(ndim),ivalue(ndim)

	integer i,k,n,n1

	do i=1,ndim-1
	  do k=ndim,i+1,-1
	    n=ipoint(k)
	    n1=ipoint(k-1)
	    if(ivalue(n).lt.ivalue(n1)) then
		ipoint(k)=n1
		ipoint(k-1)=n
	    end if
	  end do
	end do

	end

!***********************************************************
!***********************************************************
!***********************************************************

      subroutine sort_real(n,ra)

! heapsort (real, direct)

      implicit none

      integer n
      real ra(n)
!      integer ra(n)

      integer i,ir,j,l
      real rra
!      integer rra

      if(n.lt.2) return

      l=n/2+1
      ir=n
      do
        if(l.gt.1)then
          l=l-1
          rra=ra(l)
        else
          rra=ra(ir)
          ra(ir)=ra(1)
          ir=ir-1
          if(ir.eq.1)then
            ra(1)=rra
            return
          endif
        endif
        i=l
        j=l+l
	do while(j.le.ir)
          if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
          endif
          if(rra.lt.ra(j))then
            ra(i)=ra(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
	end do
        ra(i)=rra
      end do

      end

!***********************************************************

      subroutine sort_int(n,ra)

! heapsort (integer, direct)

      implicit none

      integer n
!      real ra(n)
      integer ra(n)

      integer i,ir,j,l
!      real rra
      integer rra

      if(n.lt.2) return

      l=n/2+1
      ir=n
      do
        if(l.gt.1)then
          l=l-1
          rra=ra(l)
        else
          rra=ra(ir)
          ra(ir)=ra(1)
          ir=ir-1
          if(ir.eq.1)then
            ra(1)=rra
            return
          endif
        endif
        i=l
        j=l+l
	do while(j.le.ir)
          if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
          endif
          if(rra.lt.ra(j))then
            ra(i)=ra(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
	end do
        ra(i)=rra
      end do

      end

!******************************************************************

      subroutine isort(n,ra,index)

! heapsort (integer, indirect)

      implicit none

      integer n
      integer ra(n),index(n)

      integer i,ir,j,l,irra
      integer rra

      do i=1,n
	index(i)=i
      end do

      if(n.lt.2) return

      l=n/2+1
      ir=n
      do
        if(l.gt.1)then
          l=l-1
	  irra=index(l)
          rra=ra(irra)
        else
	  irra=index(ir)
          rra=ra(irra)
          index(ir)=index(1)
          ir=ir-1
          if(ir.eq.1)then
            index(1)=irra
            return
          endif
        endif
        i=l
        j=l+l
	do while(j.le.ir)
          if(j.lt.ir)then
            if(ra(index(j)).lt.ra(index(j+1)))j=j+1
          endif
          if(rra.lt.ra(index(j)))then
            index(i)=index(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
	end do
        index(i)=irra
      end do

      end

!******************************************************************

      subroutine isortr(n,ra,index)

! heapsort (real, indirect)

      implicit none

      integer n
      real ra(n)
      integer index(n)

      integer i,ir,j,l,irra
      real rra

      do i=1,n
	index(i)=i
      end do

      if(n.lt.2) return

      l=n/2+1
      ir=n
      do
        if(l.gt.1)then
          l=l-1
	  irra=index(l)
          rra=ra(irra)
        else
	  irra=index(ir)
          rra=ra(irra)
          index(ir)=index(1)
          ir=ir-1
          if(ir.eq.1)then
            index(1)=irra
            return
          endif
        endif
        i=l
        j=l+l
	do while(j.le.ir)
          if(j.lt.ir)then
            if(ra(index(j)).lt.ra(index(j+1)))j=j+1
          endif
          if(rra.lt.ra(index(j)))then
            index(i)=index(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
	end do
        index(i)=irra
      end do

      end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine sort_invert(n,ia)

	implicit none

	integer n
	integer ia(n)

	integer i,iaux,ip

	do i=1,n/2
	  ip = n + 1 - i
	  iaux = ia(i)
	  ia(i) = ia(ip)
	  ia(ip) = iaux
	end do

	end

!******************************************************************

      function locate(n,ixx,index,ix)

! locates index for ix in array ixx ( ixx(index(i)) is sorted )

      implicit none

      integer locate
      integer n
      integer ixx(n)
      integer index(n)
      integer ix

      integer jl,ju,jm

      jl=1
      ju=n

      do while(ju.gt.jl)
        jm=(ju+jl)/2
        if(ix.gt.ixx(index(jm)))then
          jl=jm+1
        else
          ju=jm
        endif
      end do

      if(ju.le.0) then                    !empty list
        locate=0
      else if(ix.eq.ixx(index(ju))) then  !found
        locate=index(ju)
      else                                !not found
        locate=0
      end if

      return
      end

!******************************************************************

      function locater(n,rxx,index,rx)

! locates index for rx in array rxx ( rxx(index(i)) is sorted )

      implicit none

      integer locater
      integer n
      real rxx(n)
      integer index(n)
      real rx

      integer jl,ju,jm

      jl=1
      ju=n

      do while(ju.gt.jl)
        jm=(ju+jl)/2
        if(rx.gt.rxx(index(jm)))then
          jl=jm+1
        else
          ju=jm
        endif
      end do

      if(ju.le.0) then                    !empty list
        locater=0
      else if(rx.eq.rxx(index(ju))) then  !found
        locater=index(ju)
      else                                !not found
        locater=0
      end if

      return
      end

!**********************************************************

      subroutine locater_all(n,rxx,index,rx,ju1,ju2)

! locates index for rx in array rxx ( rxx(index(i)) is sorted )
! all indices are located, they lie between ju1 and ju2

      implicit none

      integer n
      real rxx(n)
      integer index(n)
      real rx
      integer ju1,ju2

      integer jl,ju,jm

      jl=1
      ju=n

      do while(ju.gt.jl)
        jm=(ju+jl)/2
        if(rx.gt.rxx(index(jm)))then
          jl=jm+1
        else
          ju=jm
        endif
      end do

      ju1 = 0
      ju2 = 0
      if(ju.le.0) then                    !empty list
        ju1 = 0
        ju2 = 0
	return
      else if(rx.eq.rxx(index(ju))) then  !found
        ju1 = ju
        ju2 = ju
      else                                !not found
        ju1 = 0
        ju2 = 0
	return
      end if

! now look for other possible indices

      do ju=ju-1,1,-1
        if(rx.eq.rxx(index(ju))) ju1 = ju
      end do

      do ju=ju+1,n,+1
        if(rx.eq.rxx(index(ju))) ju2 = ju
      end do

      return
      end

!**********************************************************

	subroutine create_permutation(n,index)

! creates permuted index

	implicit none

	integer n
	integer index(n)	!returned permutation

	integer i,nn,iaux,ierr
	integer ihelp(n)
	real r
	logical, parameter :: btest = .true.
	logical is_permutation

	do i=1,n
	  index(i) = i
	end do

	do i=n,2,-1
	  call random_number(r)
	  nn = min(i,1 + int(i*r))
	  iaux = index(i)
	  index(i) = index(nn)
	  index(nn) = iaux
	end do

	if( .not. btest ) return

! from here test routines

	if( .not. is_permutation(n,index) ) then
	  stop 'error stop create_permutation: index'
	end if

	!do i=1,n
	!  write(6,*) i,index(i)
	!end do

	end

!**********************************************************

	subroutine permute_r(n,index,array)

! permutes array with index
!
! index is an array that contains the permutation

	implicit none

	integer n
	integer index(n)	!permutation
	real array(n)		!original on entry - permuted on return

	real aux(n)

	integer i

	do i=1,n
	  aux(index(i)) = array(i)
	end do

	array = aux

	end

!**********************************************************

	function is_permutation(n,index)

! checks if index is permutation

	implicit none

	logical is_permutation
	integer n
	integer index(n)

	integer i,ip
	integer ihelp(n)

	is_permutation = .false.

	ihelp = 0

	do i=1,n
	  ip = index(i)
	  if( ip < 1 .or. ip > n ) return
	  ihelp(ip) = 1
	end do

	do i=1,n
	  if( ihelp(i) /= 1 ) return
	end do

	is_permutation = .true.

	end

!**********************************************************
!**********************************************************
!**********************************************************
! testing
!**********************************************************
!**********************************************************
!**********************************************************

	subroutine test_permutation(n)

	implicit none

	integer n
	integer index(n)

	call create_permutation(n,index)

	write(6,*) 'ok test_permutation ',n

	end

!**********************************************************

	subroutine test_sort_real_indirect(n)

	implicit none

	integer n

	integer i,n1,n2
	real ra(n)
	integer index(n)
	real r
	logical, parameter :: bwrite = .false.

	do i=1,n
	  call random_number(r)
	  ra(i) = 100 * n * r
	  if( bwrite ) write(6,*) i,ra(i)
	end do

        call isortr(n,ra,index)

	if( bwrite ) write(6,*) 1,ra(index(1))
	do i=2,n
	  n1 = index(i-1)
	  n2 = index(i)
	  if( ra(n2) < ra(n1) ) then
	    write(6,*) i,ra(n1),ra(n2)
	    stop 'error stop test_int_direct: not sorted'
	  end if
	  if( bwrite ) write(6,*) i,ra(n2)
	end do

	write(6,*) 'ok test_sort_int_indirect ',n

	end

!**********************************************************

	subroutine test_sort_int_indirect(n)

	implicit none

	integer n

	integer i,n1,n2
	integer ra(n)
	integer index(n)
	real r
	logical, parameter :: bwrite = .false.

	do i=1,n
	  call random_number(r)
	  ra(i) = 100 * n * r
	  if( bwrite ) write(6,*) i,ra(i)
	end do

        call isort(n,ra,index)

	if( bwrite ) write(6,*) 1,ra(index(1))
	do i=2,n
	  n1 = index(i-1)
	  n2 = index(i)
	  if( ra(n2) < ra(n1) ) then
	    write(6,*) i,ra(n1),ra(n2)
	    stop 'error stop test_int_direct: not sorted'
	  end if
	  if( bwrite ) write(6,*) i,ra(n2)
	end do

	write(6,*) 'ok test_sort_int_indirect ',n

	end

!**********************************************************

	subroutine test_sort_int_direct(n)

	implicit none

	integer n

	integer i
	integer ra(n)
	real r
	logical, parameter :: bwrite = .false.

	do i=1,n
	  call random_number(r)
	  ra(i) = 100 * n * r
	  if( bwrite ) write(6,*) i,ra(i)
	end do

        call sort_int(n,ra)

	if( bwrite ) write(6,*) 1,ra(1)
	do i=2,n
	  if( ra(i) < ra(i-1) ) then
	    write(6,*) i,ra(i-1),ra(i)
	    stop 'error stop test_int_direct: not sorted'
	  end if
	  if( bwrite ) write(6,*) i,ra(i)
	end do

	write(6,*) 'ok test_sort_int_direct ',n

	end

!**********************************************************

	subroutine test_sort_real_direct(n)

	implicit none

	integer n

	integer i
	real ra(n)
	real r
	logical, parameter :: bwrite = .false.

	do i=1,n
	  call random_number(r)
	  ra(i) = r
	  if( bwrite ) write(6,*) i,ra(i)
	end do

        call sort_real(n,ra)

	if( bwrite ) write(6,*) 1,ra(1)
	do i=2,n
	  if( ra(i) < ra(i-1) ) then
	    write(6,*) i,ra(i-1),ra(i)
	    stop 'error stop test_int_direct: not sorted'
	  end if
	  if( bwrite ) write(6,*) i,ra(i)
	end do

	write(6,*) 'ok test_sort_real_direct ',n

	end

!**********************************************************

	subroutine test_sort_all

	call test_sort_int_direct(100)
	call test_sort_real_direct(100)
	call test_sort_int_indirect(100)
	call test_sort_real_indirect(100)

	call test_sort_int_direct(10000)
	call test_sort_real_direct(10000)
	call test_sort_int_indirect(10000)
	call test_sort_real_indirect(10000)

	call test_permutation(100)
	call test_permutation(10000)

	end

!**********************************************************
!	program test_sort_main
!	call test_sort_all
!	end
!**********************************************************

