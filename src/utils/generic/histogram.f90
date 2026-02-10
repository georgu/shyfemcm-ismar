
!--------------------------------------------------------------------------
!
!    Copyright (C) 2004,2010-2011,2015,2019  Georg Umgiesser
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

! routines dealing with histogram
!
! contents :
!
! subroutine histo_init(nbin,bin0,dbin,rbin)
! subroutine histo_insert_value(value)
! 
! revision log :
!
! 16.03.2004	ggu	routines written from scratch
! 23.03.2010	ggu	changed v6.1.1
! 31.05.2011	ggu	changed VERS_6_1_23
! 09.01.2015	ggu	changed VERS_7_0_12
! 16.02.2019	ggu	changed VERS_7_5_60
! 28.09.2023	ggu	include substituted with module
! 22.09.2024	ggu	revisited - allocate arrays
!
!****************************************************************

!---------------------------------------------------------------------------
!	nbin					bin size
!	bin0					first value of separator
!	dbin					regular bin spacing
!	rbin(nbin)				arbitrary bin definition
!	n					total number of values
!	values(n)				values to be inserted
!	value					single value to be inserted
!
!       call histo_init(nbin,bin0,dbin)
!       call histo_init(nbin,rbin)
!
!	call histo_make(n,values,nbin,dbin,bin0)
!	call histo_make(n,values,nbin,dbin)
!	call histo_make(n,values,nbin)
!	call histo_make(n,values)		# default nbin = 20
!	
!	call histo_insert(value)
!	call histo_insert(n,values)
!
!	call histo_return(nbin,ac,ic)
!	call histo_info
!
! calling sequence:
!
!	call histo_init(nbin,bin0,dbin)
!	call histo_insert(n,values)
!	call histo_return(nbin,ac,ic)
!
! or
!
!       call histo_make(n,values,nbin)		# or similar
!	call histo_return(nbin,ac,ic)
!
!---------------------------------------------------------------------------

!================================================================
	module mod_histo
!================================================================

	implicit none

	private

        integer, save :: ncbin = 0
        integer, allocatable, save :: icount(:)
        real, allocatable, save :: acenter(:)
        real, allocatable, save :: abin(:)

        INTERFACE histo_init
        MODULE PROCEDURE histo_init_auto, histo_init_bins
        END INTERFACE

        INTERFACE histo_make
        MODULE PROCEDURE histo_make_0,histo_make_1,histo_make_2,histo_make_3
        END INTERFACE

        INTERFACE histo_insert
        MODULE PROCEDURE histo_insert_value , histo_insert_values
        END INTERFACE

	public :: histo_init , histo_insert &
     &			, histo_return , histo_make , histo_info

!================================================================
	contains
!================================================================

	subroutine histo_alloc(nbin)

	integer nbin

	if( ncbin /= 0 ) deallocate(icount,abin,acenter)
	allocate(icount(nbin),abin(nbin),acenter(nbin))
	ncbin = nbin
	abin = 0.
	acenter = 0.
	icount = 0

	end

!****************************************************************

	subroutine histo_info

	integer i

	write(6,*) 'histo info:'
	write(6,*) 'ncbin = ',ncbin
	do i=1,ncbin
	  write(6,*) i,abin(i),acenter(i),icount(i)
	end do
	write(6,*) 'end histo info:'

	end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine histo_make_3(n,values,nbin,dbin,bin0)

! sets up and computes histogram automatically

	implicit none

	integer n
	real values(n)
	integer nbin
	real dbin,bin0

	integer i
	
	call histo_init_auto(nbin,bin0,dbin)

	call histo_insert(n,values)

	end

!****************************************************************

	subroutine histo_make_2(n,values,nbin,dbin)

! sets up and computes histogram automatically

	implicit none

	integer n
	real values(n)
	integer nbin
	real dbin

	integer i
	real bin0,vmin
	
	vmin = minval(values)
	bin0 = vmin + dbin

	call histo_make_3(n,values,nbin,dbin,bin0)

	end

!****************************************************************

	subroutine histo_make_1(n,values,nbin)

! sets up and computes histogram automatically

	implicit none

	integer n
	real values(n)

	integer nbin,i
	real vmin,vmax,dv,dbin,bin0
	real rnext
	
	vmin = minval(values)
	vmax = maxval(values)

	dv = vmax - vmin
	dbin = dv / nbin
	!bin0 = vmin + dbin/2.
	bin0 = vmin + dbin

	call histo_make_3(n,values,nbin,dbin,bin0)

	end

!****************************************************************

	subroutine histo_make_0(n,values)

! sets up and computes histogram automatically

	implicit none

	integer n
	real values(n)

	integer nbin
	
	nbin = 20
	call histo_make_1(n,values,nbin)

	end

!****************************************************************
!****************************************************************
!****************************************************************

	subroutine histo_return(nbin,ac,ic)

	implicit none

	integer nbin
	real ac(nbin)
	integer ic(nbin)

	if( nbin == 0 ) then
	  nbin = ncbin
	  return
	end if

	if( nbin < ncbin ) then
	  error stop 'error stop histo_return: nbin<ncbin'
	end if

	nbin = ncbin
	ac = acenter
	ic = icount

	end

!****************************************************************

        subroutine histo_init_auto(nbin,bin0,dbin)

! sets up icount and abin

        implicit none

        integer nbin            !total number of bins
        real bin0               !first bin (limit)
        real dbin               !regular bin size (0 => use rbin)

        integer i

	call histo_alloc(nbin)

        do i=1,nbin
          abin(i) = bin0 + (i-1) * dbin
        end do

	call histo_make_center(nbin,abin,acenter)

        end

!****************************************************************

        subroutine histo_init_bins(nbin,rbin)

! sets up icount and abin

        implicit none

        integer nbin            !total number of bins
        real rbin(nbin)         !bin size limits (upper)

	call histo_alloc(nbin)

	abin = rbin

	call histo_make_center(nbin,abin,acenter)

        end

!****************************************************************

        subroutine histo_make_center(nbin,ab,ac)

! sets up icount and abin

        implicit none

        integer nbin            !total number of bins
	real ab(nbin)
	real ac(nbin)

        integer i

        do i=2,nbin-1
          ac(i) = 0.5*( ab(i) + ab(i-1) )
        end do
	ac(1) = ac(2) - (ab(2)-ab(1))
	ac(nbin) = ac(nbin-1) + (ab(nbin-1)-ab(nbin-2))

        end

!****************************************************************
!****************************************************************
!****************************************************************

        subroutine histo_insert_value(value)

        implicit none

        real value

        integer i

        do i=1,ncbin-1
          if( value .le. abin(i) ) exit
        end do

        icount(i) = icount(i) + 1

        end

!****************************************************************

        subroutine histo_insert_values(n,values)

        implicit none

	integer n
        real values(n)

        integer i

        do i=1,n
	  call histo_insert_value(values(i))
        end do

        end

!****************************************************************

!================================================================
	end module mod_histo
!================================================================

	subroutine histo_test1(text,n,amed,drange,sigma)

	use mod_histo

	implicit none

	character*(*) text
	integer n
	real amed,drange,sigma

	integer i,nbin
	real r,vmin,vmax
	real, allocatable :: ac(:)
	integer, allocatable :: ic(:)
	real, allocatable :: values(:)
	character*80 file

	nbin = 20
	allocate(ac(nbin),ic(nbin))
	allocate(values(n))

	do i=1,n
	  call random_number(r)
	  r = r * drange + amed
	  values(i) = r
	end do
	vmin = minval(values)
	vmax = maxval(values)
	  
	write(6,*) 'n,vmin,vmax: ',n,vmin,vmax

	call histo_make(n,values,nbin)
        call histo_return(nbin,ac,ic)

	file=trim(text)//'.tmp'
	open(1,file=file,status='unknown',form='formatted')

	write(6,*) nbin
	do i=1,nbin
	  write(6,*) i,ac(i),ic(i)
	  write(1,*) ac(i),ic(i)
	end do

	close(1)

	end

!****************************************************************

	subroutine histo_test
	call histo_test1('s00',1000000,300.,100.,0.)
	call histo_test1('s30',1000000,300.,100.,30.)
	end

!****************************************************************
	!program main_histo_test
	!call histo_test
	!end
!****************************************************************

