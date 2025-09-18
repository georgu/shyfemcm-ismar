
!--------------------------------------------------------------------------
!
!    Copyright (C) 2012,2015-2017,2019  Georg Umgiesser
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

! omp administration routines
!
! revision log :
!
! 25.01.2012	ggu	changed VERS_6_1_42
! 30.03.2012	ggu	changed VERS_6_1_51
! 29.08.2012	ggu	changed VERS_6_1_56
! 26.02.2015	ggu	changed VERS_7_1_5
! 05.05.2015	ggu	changed VERS_7_1_10
! 09.09.2016	ggu	changed VERS_7_5_17
! 12.01.2017	ggu	changed VERS_7_5_21
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 24.07.2024	ggu	revisited and omp enabled
! 29.01.2025	ggu	nomp==-1 uses maximal available procs (OMP_NUM_THREADS)
! 18.09.2025	ggu	some debug code, documentation

!***************************************************************

!===============================================================
	module omp_admin
!===============================================================

!$      use omp_lib

	logical, save :: bomp_init = .false.
	integer, save :: nomp_max = 0
	integer, save :: nomp_use = 0

	logical, save :: bomp_debug = .false.

!===============================================================
	end module omp_admin
!===============================================================

	subroutine openmp_init

	use omp_admin

	implicit none

	bomp_init = .true.
	nomp_max = 1
!$      nomp_max = omp_get_max_threads()
	nomp_use = nomp_max

	if( bomp_debug ) then
	  write(6,*) 'openmp_init: nomp_max,nomp_use: ',nomp_max,nomp_use
	end if

	end

!***************************************************************

	subroutine openmp_get_max_threads(n)

	use omp_admin

	implicit none

	integer n

	if( .not. bomp_init ) call openmp_init

	n = nomp_max

	end

!***************************************************************

	subroutine openmp_set_num_threads(n)

! this sets the numer of threads that are used

	use omp_admin

	implicit none

	integer n,m

	if( .not. bomp_init ) call openmp_init

	nomp_use = n
	if( n < 0 ) nomp_use = nomp_max		!use maximal possible value
	nomp_use = min(nomp_use,nomp_max)
	nomp_use = max(nomp_use,1)
	
	if( bomp_debug ) then
	  write(6,*) 'openmp_set_num_threads: nomp_max,nomp_use: ' &
     &				,nomp_max,nomp_use
	end if

!$	call omp_set_num_threads(nomp_use)

	if( bomp_debug ) then
	  m = 0
!$	  m = omp_get_max_threads()
	  write(6,*) '...',n,nomp_use,m
	end if

	end

!***************************************************************

	subroutine openmp_get_num_threads(n)

! this gets the numer of threads that are used

	use omp_admin

	implicit none

	integer n

	if( .not. bomp_init ) call openmp_init

	n = nomp_use

	end

!***************************************************************

	subroutine openmp_get_thread_num(it)

! in a parallel section it returns the actual thread number [0,nomp_use-1]

	use omp_admin

	implicit none

	integer it

	if( .not. bomp_init ) call openmp_init

	it = 0
!$	it = omp_get_thread_num ()

	end

!***************************************************************

        function openmp_is_master()

	use omp_admin

        implicit none
    
        logical openmp_is_master

	if( .not. bomp_init ) call openmp_init

        openmp_is_master = .true.
!$	openmp_is_master = ( omp_get_thread_num() == 0 )

        end

!***************************************************************

	subroutine openmp_parallel_code(text)

	use omp_admin

	implicit none

	character*(*) text

	if( .not. bomp_init ) call openmp_init

	text = 'serial'
	if ( nomp_max > 1 ) text = 'omp'

	end

!***************************************************************

	function openmp_is_parallel()

! true if omp is available

	use omp_admin

	implicit none

	logical openmp_is_parallel

	if( .not. bomp_init ) call openmp_init

	openmp_is_parallel = ( nomp_max > 1 )

	end

!***************************************************************

        function openmp_in_parallel()

! true if in parallel section

	use omp_admin

        implicit none

        logical openmp_in_parallel

	if( .not. bomp_init ) call openmp_init

        openmp_in_parallel = .false.
!$	openmp_in_parallel = omp_in_parallel()

        end

!***************************************************************
!***************************************************************
!***************************************************************

        subroutine get_clock_count(count)

        implicit none

        integer count
        integer count_rate,count_max

        call system_clock(count,count_rate,count_max)

        end

!***************************************************************

        subroutine get_clock_count_diff(count0,dcount)

        implicit none

        integer count0,dcount
        integer count,count_rate,count_max

        call system_clock(count,count_rate,count_max)

        dcount = count - count0

        end

!***************************************************************

        subroutine get_clock_time(time)

        implicit none

        real time

        call cpu_time(time)

        end

!***************************************************************

        subroutine get_clock_time_diff(time0,dtime)

        implicit none

        real time0,dtime

        real time

        call cpu_time(time)

        dtime = time - time0

        end

!***************************************************************

        function openmp_get_wtime()

! gets time

        implicit none

        double precision openmp_get_wtime

        openmp_get_wtime = 0.

        end

!***************************************************************
!***************************************************************
!***************************************************************

        subroutine omp_compute_chunk(imax,nchunk)

	use omp_admin

        implicit none

        integer imax,nchunk

        nchunk = 1

        end

!***************************************************************

        subroutine omp_compute_minmax(nchunk,imax,i,iend)

        implicit none

        integer nchunk,imax,i,iend

        iend = i

        end

!***************************************************************

        subroutine omp_set_minmax(nuse,nmax,ichunk,istart,iend)

	use omp_admin

        implicit none

        integer nuse,nmax,ichunk,istart,iend
	integer nchunk

	if( .not. bomp_init ) call openmp_init

	if( nuse <= 0 ) nuse = nomp_use
	nchunk = nmax / nuse
	istart = 1 + (ichunk-1)*nchunk
	iend = ichunk*nchunk
	if( ichunk == nuse .and. iend /= nmax ) iend = nmax

	end

!***************************************************************
!***************************************************************
!***************************************************************
! test routines
!***************************************************************
!***************************************************************
!***************************************************************

	subroutine multi_omp_shell(nwant,nloop)

	implicit none
	integer nwant,nloop
	integer nres,it,i,nl,nuse
	integer nchunk,ie,is,ntot,nt,nr
	integer count,dcount
	double precision rres,rr

	nuse = nwant
	if( nuse == 0 ) call openmp_get_max_threads(nuse)
	ntot = 0
	rres = 0
	call openmp_set_num_threads(nuse)

	call get_clock_count(count)

!! !$OMP PARALLEL 
!! !$OMP SINGLE
!! !$OMP TASK PRIVATE(i,it,nres) SHARED(nl,nloop,nuse)    DEFAULT(NONE)
!$OMP PARALLEL DO PRIVATE(i,it,is,ie,nt,nr,rr) & 
!$OMP & SHARED(nloop,nl,nuse,nchunk) DEFAULT(NONE) &
!$OMP & REDUCTION(+:rres,ntot)
	do i=1,nuse
          call omp_set_minmax(nuse,nloop,i,is,ie)
	  call multi_omp(is,ie,nt,rr)
	  !write(6,*) i,is,ie,nt,rr
	  ntot = ntot + nt 
	  rres = rres + rr 
	end do
!$OMP END PARALLEL DO
!! !$OMP END TASK
!! !$OMP END SINGLE
!! !$OMP TASKWAIT
!! !$OMP END PARALLEL      

	write(6,*) 'final: ',ntot,rres

        call get_clock_count_diff(count,dcount)
	write(6,*) 'time of execution: ',nuse,dcount

	end

!***************************************************************

	subroutine multi_omp(is,ie,ntot,rres)
	implicit none
	integer is,ie,ntot,i
	double precision rres
	ntot = 0
	rres = 0
	do i=is,ie
	  rres = rres + 2*i
	  ntot = ntot + 1
	end do
	end

!***************************************************************

	subroutine test_omp_admin
	use omp_admin
	implicit none
	integer nmax,nuse
	integer nloop,nres,it,i
	write(6,*) 'testing omp admin routines'
	call openmp_get_max_threads(nmax)
	call openmp_set_num_threads(nmax/2)
	call openmp_get_num_threads(nuse)
	write(6,*) 'nmax = ',nmax
	write(6,*) 'nuse = ',nuse

	nloop = 100000000
	call srand(4637)

	call multi_omp_shell(1,nloop)
	call multi_omp_shell(3,nloop)
	call multi_omp_shell(nmax,nloop)
	call multi_omp_shell(10,nloop)
	call multi_omp_shell(20,nloop)
	call multi_omp_shell(0,nloop)

	end

!***************************************************************

!	program main_omp_admin
!	call test_omp_admin
!	end

!***************************************************************

