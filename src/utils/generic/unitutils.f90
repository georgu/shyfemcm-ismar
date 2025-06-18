
!--------------------------------------------------------------------------
!
!    Copyright (C) 2018-2020  Georg Umgiesser
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

! revision log :
!
! 06.07.2018	ggu	changed VERS_7_5_48
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 18.03.2020	ggu	handle case if no node number is given (j=0)
! 12.06.2025	ggu	if string is available add it to file name

!**************************************************************************

	module shyelab_unit

	logical, save :: b_use_new_format = .false.
	integer, save :: iunit = 100
	character*80, save :: string_save = ' '

	end module shyelab_unit

!***************************************************************
!***************************************************************
!***************************************************************

	subroutine shyelab_unit_set_new_format(bval)

	use shyelab_unit

	implicit none

	logical bval

	b_use_new_format = bval

	end

!***************************************************************

	subroutine shyelab_unit_get_new_format(bval)

	use shyelab_unit

	implicit none

	logical bval

	bval = b_use_new_format

	end
	
!***************************************************************

	subroutine get_new_unit(iu)

	use shyelab_unit

	implicit none

	integer iu
	logical bopen

	do
	  iunit = iunit + 1
	  inquire(unit=iunit,opened=bopen)
	  if( .not. bopen ) exit
	end do

	iu = iunit

	end

!***************************************************************

	subroutine set_iunit_string(string)

! transforms " ()" into underscores

	use shyelab_unit

	character*(*) string

	string_save = string

	if( .not. b_use_new_format ) string_save = ' '

	do i=1,len_trim(string_save)
	  if( string_save(i:i) == ' ' ) string_save(i:i) = '_'
	  if( string_save(i:i) == '(' ) string_save(i:i) = '_'
	  if( string_save(i:i) == ')' ) string_save(i:i) = '_'
	end do

	end
	
!***************************************************************
!***************************************************************
!***************************************************************

	subroutine make_iunit_name(short,modi,dim,j,iu)

	use shyelab_unit

	implicit none

	character*(*) short,modi,dim
	integer j
	integer iu

	logical bopen
	character*80 numb
	character*80 dimen
	character*80 string
	character*80 extension
	character*80 name

	if( j <= 0 ) then	!no node number given
	  numb = '0'
	else
          write(numb,'(i5)') j
          numb = adjustl(numb)
	end if

	string = ' '
	extension = ' '

	if( b_use_new_format ) then
	  if( string_save /= ' ' ) then
	    string = '.' // trim(string_save)
	  end if
	  extension = '.txt'
	end if

	dimen = '.' // trim(dim) // '.'
	name = trim(short)//trim(modi)//trim(dimen) &
     &			//trim(numb)//trim(string)//trim(extension)

	call get_new_unit(iu)

	inquire(file=name,opened=bopen)
	if( bopen ) goto 99
	inquire(unit=iu,opened=bopen)
	if( bopen ) goto 98

        !write(6,*) 'opening file : ',iu,'  ',trim(name)
        open(iu,file=name,form='formatted',status='unknown')

	return
   99	continue
	write(6,*) 'file already open: ',trim(name)
	stop 'error stop make_iunit_name: internal error (1)'
   98	continue
	write(6,*) 'unit already open: ',iu
	stop 'error stop make_iunit_name: internal error (2)'
	end

!***************************************************************

