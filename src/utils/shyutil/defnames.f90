
!--------------------------------------------------------------------------
!
!    Copyright (C) 1997-1998,2000-2002,2007,2010-2012  Georg Umgiesser
!    Copyright (C) 2014,2017-2019  Georg Umgiesser
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

! create default names
!
! contents :
!
! subroutine mkname(dir,name,ext,file)            makes file name
! subroutine defmak(defdir,defnam,ext,file)	  makes file with defaults
! function ideffi(defdir,defnam,ext,form,status)  opens file in default dir
! function ifemop(ext,form,status)		  opens file with default name
! function ifemopa(text,ext,form,status)	  opens file with default name
! function ifem_open_file(ext,status)		  opens unformated file
! function ifem_test_file(ext,status)		  tries to open unformated file
!
! revision log :
!
! 23.05.1997	ggu	$$EXTENS - default extension may be overwritten
! 18.06.1997	ggu	restructured - idefna,idefop,idefts,idefun deleted
! 16.01.1998	ggu	idefop reintroduced -> to avoid link error
! 21.01.1998	ggu	in mkname: give extension with or without dot
! 08.08.2000	ggu	new routine ifemop
! 27.11.2001	ggu	error message rewritten
! 11.10.2002	ggu	new subroutine deffile
! 07.03.2007	ggu	new routine ifem_open_file
! 23.03.2010	ggu	changed v6.1.1
! 29.04.2010	ggu	new routine ifem_open_file
! 03.05.2010	ggu	new routine ifem_choose_file() and add_extension()
! 02.07.2011	ggu	idefna,idefop finally deleted
! 13.07.2011	ggu	cleaned from old structures
! 18.08.2011	ggu	bug fix in idefbas -> use status passed in
! 01.06.2012	ggu	changed VERS_6_1_53
! 18.06.2014	ggu	changed VERS_6_1_77
! 09.05.2017	ggu	add_extension -> to subst_extension, new add_extension
! 05.12.2017	ggu	changed VERS_7_5_39
! 07.12.2017	ggu	changed VERS_7_5_40
! 05.10.2018	ggu	avoid run time error in subst_extension()
! 16.10.2018	ggu	changed VERS_7_5_50
! 16.02.2019	ggu	changed VERS_7_5_60
! 03.05.2019	ggu	new routines to get extension and name
! 21.05.2019	ggu	changed VERS_7_5_62
! 08.04.2024	ggu	new routine make_name_with_number()
! 10.05.2024	ggu	new routine parse_file_name()
! 21.11.2024	ggu	in make_name_with_number() use no dot if ext == ''
!
! notes :
!
! exchange ideffi with ifemop (many apperances)
! eliminate call to getpar/getfnm
!
!**************************************************************

        subroutine mkname(dir,name,ext,file)

! makes file name given its constituents
!
! dir   directory
! name  name
! ext   extension (with or without dot)
! file  created file name (return)

        implicit none

! arguments
        character*(*) dir,name,ext,file
! local
        integer nall,nstart,nend,naux
! function
        integer ichafs,ichanm

        nall=1
        file=' '

        nstart=ichafs(dir)
        nend=ichanm(dir)
        if(nend.gt.0) then
		file(nall:)=dir(nstart:nend)
        	nall=nall+nend-nstart+1
	end if

        nstart=ichafs(name)
        nend=ichanm(name)
        if(nend.gt.0) then
		file(nall:)=name(nstart:nend)
       		nall=nall+nend-nstart+1
	end if

	call subst_extension(file,ext,.false.)

	return

	end

!**************************************************************
!**************************************************************
!**************************************************************

	module default_names

        character*80, save :: def_bas = ' '
        character*80, save :: def_nam = ' '

	end module default_names

!**************************************************************

	subroutine set_default_names(defbas,defnam)

! sets default names to be used with def routines
!
! must be called before any other routine in subdef can be used

	use default_names

	implicit none

        character*(*) defbas,defnam
	
	def_bas = defbas
	def_nam = defnam

	end

!**************************************************************

	subroutine get_default_names(defbas,defnam)

! gets default names to be used with def routines

	use default_names

        implicit none

        character*(*) defbas,defnam
	
	defbas = def_bas
	defnam = def_nam

	end
	
!**************************************************************
!**************************************************************
!**************************************************************

        subroutine def_make(ext,file)

! makes file with defaults supplied
!
! ext   extension (with dot)
! file  created file name (return)

	use default_names

        implicit none

        character*(*) ext,file

        character*80 dir,name

	name = def_nam
	dir = ' '
        !call getfnm('datdir',dir)	! this has to be deleted
        call getfnm('runnam',name)	! this has to be deleted

	call mkname(dir,name,ext,file)

	end

!**************************************************************

        subroutine defmak(defdir,defnam,ext,file)

! FIXME -> take defdir,defnam from common block

! makes file with defaults supplied
!
! defdir   directory
! defnam  name
! ext   extension (with dot)
! file  created file name (return)

        implicit none

        character*(*) defdir,defnam,ext,file
	character*80 dir,name

	dir = ' '
        !if( defdir .ne. ' ' ) call getfnm(defdir,dir)
        call getfnm(defnam,name)

	call mkname(dir,name,ext,file)

	end

!**************************************************************

        function ideffi(defdir,defnam,ext,form,status)

! FIXME -> substitute with ifemop (many appearances)

! opens file in default dir

! defdir   directory
! defnam  name
! ext   extension (with dot)
! form  formatted ?
! status open status

        implicit none

	integer ideffi
        character*(*) defdir,defnam,ext,status,form
	character*80 file
	integer ifileo

	call defmak(defdir,defnam,ext,file)
        call def_make(ext,file)
	ideffi=ifileo(0,file,form,status)

	end

!**************************************************************

	function idefbas(basnam,status)

	implicit none

	integer idefbas
	character*(*) basnam
	character*(*) status

	character*80 name
	integer ifileo

        name = basnam
        call add_extension(name,'.bas')

        idefbas=ifileo(0,name,'unform',status)

	end

!**************************************************************
!**************************************************************
!**************************************************************
! opening of default simulation
!**************************************************************
!**************************************************************
!**************************************************************

        function ifemop(ext,form,status)

! opens file with default name (run) and extension given for fem model
! returns with error code

! ext   extension (with dot)
! form  formatted ?
! status open status

        implicit none

	integer ifemop
        character*(*) ext,status,form

	character*80 file,defdir,defnam
	integer ifileo

        call def_make(ext,file)
	ifemop=ifileo(0,file,form,status)

	end

!**************************************************************

        function ifemopa(text,ext,form,status)

! opens file with default name (run) and extension given for fem model
! in case of error exits with error message text

! text  error message
! ext   extension (with dot)
! form  formatted ?
! status open status

        implicit none

	integer ifemopa
        character*(*) text,ext,status,form

	character*80 file,defdir,defnam
	integer ifileo

        call def_make(ext,file)
	ifemopa=ifileo(0,file,form,status)

	if( ifemopa .le. 0 ) then
	  write(6,*) 'error opening file ',file
	  write(6,*) text
	  stop 'error stop ifemopa: opening file'
	end if

	end

!**************************************************************
!**************************************************************
!**************************************************************
! unformatted opening of default simulation
!**************************************************************
!**************************************************************
!**************************************************************

        function ifem_open_file(ext,status)

! opens unformated file with default name (run) and extension given
! in case of error exits 

! ext		extension (with dot)
! status	open status

        implicit none

	integer ifem_open_file
        character*(*) ext,status

	integer ifem_test_file

	ifem_open_file = ifem_test_file(ext,status)

	if( ifem_open_file .le. 0 ) then
	  stop 'error stop ifem_open_file: opening file'
	end if

	end

!**************************************************************

        function ifem_test_file(ext,status)

! tries to open unformated file with default name (run) and extension given

! ext		extension (with dot)
! status	open status

        implicit none

	integer ifem_test_file
        character*(*) ext,status

	character*80 file,defdir,defnam
	character*80 form
	integer ifileo

	form = 'unform'
        call def_make(ext,file)

	ifem_test_file = ifileo(0,file,form,status)

	if( ifem_test_file .le. 0 ) then
	  write(6,*) 'cannot open file ',file
	end if

	end

!**************************************************************

        function ifem_choose_file(ext,status)

! tries to open unformated file with default name (run) and extension given
! insists on extension -> if name has extension substitute it with ext

! ext		extension (with dot)
! status	open status

        implicit none

	integer ifem_choose_file
        character*(*) ext,status

	character*80 file,defdir,defnam
	character*80 form
	integer ifileo

	form = 'unform'
        call def_make(ext,file)
	call subst_extension(file,ext,.true.)

	ifem_choose_file = ifileo(0,file,form,status)

	if( ifem_choose_file .le. 0 ) then
	  write(6,*) 'cannot open file ',file
	end if

	end

!**************************************************************
!**************************************************************
!**************************************************************

	subroutine parse_file_name(file,dir,name,ext)

	implicit none

	character*(*) file	!file name on entry
	character*(*) dir	!directory on return
	character*(*) name	!name on return
	character*(*) ext	!extension on return

	integer idir,iext,istart,iend
	logical, parameter :: bback = .true.

	idir = index(file,'/',bback)
	if( idir == 0 ) then
	  dir = " "
	  istart = 1
	else
	  dir = file(1:idir)
	  istart = idir + 1
	end if

	iext = index(file,'.',bback)
	if( iext == 0 ) then
	  ext = " "
	  iend = len_trim(file)
	else
	  ext = file(iext+1:)
	  iend = iext - 1
	end if

	name = file(istart:iend)

	end

!**************************************************************
!**************************************************************
!**************************************************************

	subroutine add_extension(name,ext)

! adds extension to file name if not already there
!
! ext should have .

	implicit none

	character*(*) name
	character*(*) ext

	integer nall,next,n
	integer ichanm

	if( ext(1:1) /= '.' ) then
	  write(6,*) 'ext: ',ext
	  write(6,*) 'extension must have .'
	  stop 'error stop add_extension: no dot'
	end if

	nall = ichanm(name)
	next = ichanm(ext)

	n = nall - next + 1
	if( n > 0 ) then
	  if( name(n:nall) == ext ) then	!already there
	    !nothing
	  else
	    name(nall+1:) = ext
	  end if
	else
	  name(nall+1:) = ext
	end if

	end

!**************************************************************

	subroutine delete_extension(name,ext)

! deletes extension from file
!
! ext should have .

	implicit none

	character*(*) name
	character*(*) ext

	integer nall,next,n
	integer ichanm

	if( ext(1:1) /= '.' ) then
	  write(6,*) 'ext: ',ext
	  write(6,*) 'extension must have .'
	  stop 'error stop add_extension: no dot'
	end if

	nall = ichanm(name)
	next = ichanm(ext)

	n = nall - next + 1
	if( n > 0 ) then
	  if( name(n:nall) == ext ) then	!extension is there
	    name(n:) = ' '
	  end if
	end if

	end

!**************************************************************

	subroutine change_extension(name,extold,extnew)

! changes extension with new one
!
! ext should have .

	implicit none

	character*(*) name
	character*(*) extold,extnew

	integer nall,next,n
	integer ichanm

	if( extold(1:1) /= '.' .or. extnew(1:1) /= '.' ) then
	  write(6,*) 'extold: ',extold
	  write(6,*) 'extnew: ',extnew
	  write(6,*) 'extension must have .'
	  stop 'error stop add_extension: no dot'
	end if

	nall = ichanm(name)
	next = ichanm(extold)

	n = nall - next + 1
	if( n > 0 ) then
	  if( name(n:nall) == extold ) then	!extension is there
	    name(n:) = extnew
	  end if
	end if

	end

!**************************************************************

	subroutine subst_extension(name,ext,bforce)

! substitutes extension with given one
!
! if bforce is true substitutes given extension
! otherwise only adds if not already there
!
! extension must have 3 chars
!
! name		file name (with or without extension)
! ext		extension (with or without dot)
! bforce	force substitution of extension, even if already there

	implicit none

	character*(*) name
	character*(*) ext
	logical bforce

	integer nall,n,nstart,nend

        integer ichafs,ichanm

	nall = 1 + ichanm(name)

	n = nall - 4		!here should be the dot
	if( n .gt. 0 ) then	!has extension
	 if( name(n:n) .eq. '.' ) then	!has extension
	  if( bforce ) then	!substitute extension
	    nall = n
	  else
	    return		!leave extension
	  end if
	 end if
	end if

        nstart=ichafs(ext)
        nend=ichanm(ext)

	if( nend .gt. 0 ) then
	   if( ext(nstart:nstart) .ne. '.' ) then !add dot if not in ext
		name(nall:nall) = '.'
		nall = nall + 1
	   end if
	   name(nall:)=ext(nstart:nend)
	end if

	end

!**************************************************************

	subroutine check_extension(file,ext)

! finds extension of file and returns it

	implicit none

	character*(*) file
	character*(*) ext

	integer i

	ext = ' '

	i = index(file,'.',.true.)
	if( i == 0 ) return		!no extension

	ext = file(i+1:)

	end

!**************************************************************

	subroutine check_name_and_extension(file,name,ext)

! finds name and extension of file and returns it

	implicit none

	character*(*) file
	character*(*) name,ext

	integer i

	name = ' '
	ext = ' '

	i = index(file,'.',.true.)
	if( i == 0 ) return		!no extension

	name = file(:i-1)
	ext = file(i+1:)

	end

!**************************************************************

        subroutine make_name_with_number(base,n,ext,name)

! makes name as: base.n.ext and returns it back in name
!
! extension ext should be supplied without "."

        implicit none

        character*(*) base
        integer n
        character*(*) ext
        character*(*) name

        character*80 aux

        write(aux,'(i10)') n
        aux = adjustl(aux)

	if( ext == ' ' ) then
          name = trim(base) // '.' // trim(aux)
	else
          name = trim(base) // '.' // trim(aux) // '.' // trim(ext)
	end if

        end

!**************************************************************
!**************************************************************
!**************************************************************
! test routines
!**************************************************************
!**************************************************************
!**************************************************************

        subroutine test_defnames_check(file)

        implicit none

        character*(*) file

        character*80 dir        !directory on return
        character*80 name       !name on return
        character*80 ext        !extension on return

        call parse_file_name(file,dir,name,ext)

        write(6,*) '=================='
        write(6,*) trim(file)
        write(6,*) trim(dir)
        write(6,*) trim(name)
        write(6,*) trim(ext)
        write(6,*) '=================='

        end

!**************************************************************

	subroutine test_defnames

	implicit none

        call test_defnames_check('/home/georg/test.f90')
        call test_defnames_check('test.f90')
        call test_defnames_check('/home/georg/test')
        call test_defnames_check('/test.f90')

	end

!**************************************************************

!        program main_test_defnames
!        call defnames
!        end

!**************************************************************

