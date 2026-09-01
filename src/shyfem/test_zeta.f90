
!--------------------------------------------------------------------------
!
!    Copyright (C) 2021-2023  Georg Umgiesser, Celia Laurent
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

!  testing zeta values
! 
!  revision log :
!
! 13.03.2021	ggu	First round for PETSC solver inclusion
! 23.04.2021	clr	implemented
! 26.09.2023	clr	conversion of fem3d file contents from f77

!==================================================================
        module mod_test_zeta
!==================================================================

        implicit none

	integer, private, save :: nn_step = 0

!==================================================================
        contains
!==================================================================

	subroutine test_zeta_init
	use mod_rungekutta, only : get_rungekutta_idiag_weights
	use basin
	use shympi

	implicit none

        character*4 my_id_s,n_threads_s
	real ai_ll

#if !defined(test_zeta)
	return
#endif

        write(n_threads_s,'(i4.4)')n_threads
        write(my_id_s,'(i4.4)')my_id+1

	call get_rungekutta_idiag_weights(1,ai_ll)

	if (ai_ll == 0. ) then
           open(unit = 99998,  &
     &          file = "test_zeta."//trim(n_threads_s)// &
     &                 ".EXPLICIT_"//trim(my_id_s))
	else 
           open(unit = 99998,  &
     &          file = "test_zeta."//trim(n_threads_s)// &
     &                 ".IMPLICIT_"//trim(my_id_s))
        endif

        if(nkn<100) then
          nn_step=1
        else
          nn_step=nkn/100
        endif

	end subroutine test_zeta_init

!******************************************************************

	subroutine test_zeta_write

	use basin
	use mod_hydro

	implicit none

	integer nn

#if !defined(test_zeta)
	return
#endif

        write(99998,'(i10)') 0   
        do nn=1,nkn,nn_step 
          write(99998,'(i10,E25.15)') nn,znv(nn)   
        end do 

	end subroutine test_zeta_write

!******************************************************************

	subroutine test_zeta_close

	implicit none

#if !defined(test_zeta)
	return
#endif

	close(99998)

	end subroutine test_zeta_close

!==================================================================
        end module mod_test_zeta
!==================================================================

