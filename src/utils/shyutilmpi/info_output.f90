
!--------------------------------------------------------------------------
!
!    Copyright (C) 2024  Georg Umgiesser
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

!--------------------------------------------------------------------------
!
! writes info output to file
!
! revision log :
!
! 01.12.2024    ggu     written from scratch
! 03.12.2024    ggu     allow for optional format
! 05.11.2025    ggu     routine getinfo() transferred here from iostr.f90
!
!--------------------------------------------------------------------------

!================================================================
        module mod_info_output
!================================================================

	implicit none

	private

        logical, parameter :: binfo_do = .true.	!enables info output
        !logical, parameter :: binfo_do = .false.	!enables info output

        logical, parameter :: binfo = binfo_do		!accessible outside
        logical, save :: binfo_init = .false.		!has been initialized?
	integer, save :: iuinfo = 0			!unit for writing info

        INTERFACE info_output
        MODULE PROCEDURE  info_output_scalar &
     &                   ,info_output_array
        END INTERFACE

	public :: binfo
	public :: info_output

!================================================================
        contains
!================================================================

	subroutine info_output_scalar(text,what,scalar,btime,format)

	character*(*) text
	character*(*) what
	real scalar
	logical, optional :: btime
	character*(*), optional :: format

	real array(1)

	array(1) = scalar
	call info_output_array(text,what,1,array,btime,format)

	end subroutine info_output_scalar

!****************************************************************

	subroutine info_output_array(text,what,n,array,btime,format)

	use shympi

	character*(*) text
	character*(*) what
	integer n
	real array(n)
	logical, optional :: btime		!write time string
	character*(*), optional :: format

	logical bt
	integer i,nl
	character*20 aline
	character*80 string

	if( .not. binfo_do ) return

	if( .not. binfo_init ) then
	  if( bmpi_master ) then
	    call getinfo(iuinfo)		!FIXME - if error routine hangs
	  end if
	  binfo_init = .true.
	end if

	bt = .false.
	if( present(btime) ) bt = btime

	if( what == ' ' ) then
	  !nothing
	else if( what == 'none' ) then
	  !nothing
	else if( what == 'max' ) then
	  call shympi_array_reduce('max',array)
	else if( what == 'sum' ) then
	  call shympi_array_reduce('sum',array)
	else
	  write(6,*) 'what = ',trim(what)
	  stop 'error stop info_output: unknown what'
	end if
	
	if( .not. bmpi_master ) return

	if( bt ) then
	  call get_act_timeline(aline)
	  string = trim(text)//': '//aline
	  nl = len_trim(string)
	else
	  string = trim(text)//': '
	  nl = len_trim(string) + 1
	end if

	if( present(format) ) then
	  write(iuinfo,format) string(1:nl),array
	else
	  write(iuinfo,*) string(1:nl),array
	end if

	end subroutine info_output_array

!================================================================
        end module mod_info_output
!================================================================

        subroutine getinfo(iunit)

! gets unit of info file

        implicit none

        integer iunit

        integer ifemop

        integer, save :: iu = 0

        if( iu .le. 0 ) then
          iu = ifemop('.inf','formatted','new')
          if( iu .le. 0 ) then
            write(6,*) 'error in opening info file'
            stop 'error stop getinfo'
          end if
        end if

        iunit = iu

        end

!**********************************************************************

