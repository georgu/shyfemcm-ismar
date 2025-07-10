
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

! error stop routines
!
! revision log :
!
! 12.12.2024	ggu	written from scratch
! 10.07.2025	ggu	added sucess_stop()

!===============================================================
	module mod_error_stop
!===============================================================

	implicit none

	private

	integer, save :: i_code_error = 7
	integer, save :: i_code_success = 0

        INTERFACE error_stop
        MODULE PROCEDURE          error_stop_2 &
     &                          , error_stop_1 &
     &                          , error_stop_0 &
     &                          , error_stop_i0
        END INTERFACE

	public :: error_stop
	public :: success_stop

!===============================================================
	contains
!===============================================================

	subroutine error_stop_2(routine,text)

	character*(*) routine
	character*(*) text

	write(6,*) 'error stop '//trim(routine)//': ',trim(text)
	flush(6)
	call exit(i_code_error)

	end subroutine error_stop_2

!******************************************************************

	subroutine error_stop_1(text)

	character*(*) text

	write(6,*) 'error stop: ',trim(text)
	flush(6)
	call exit(i_code_error)

	end subroutine error_stop_1

!******************************************************************

	subroutine error_stop_0

	write(6,*) 'error stop'
	flush(6)
	call exit(i_code_error)

	end subroutine error_stop_0

!******************************************************************

	subroutine error_stop_i0(ierr)

	integer ierr

	write(6,*) 'error stop'
	flush(6)
	call exit(ierr)

	end subroutine error_stop_i0

!******************************************************************

	subroutine success_stop

	call exit(i_code_success)

	end subroutine success_stop

!===============================================================
	end module mod_error_stop
!===============================================================

