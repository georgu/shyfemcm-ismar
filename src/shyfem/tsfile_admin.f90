
!--------------------------------------------------------------------------
!
!    Copyright (C) 2012-2016,2018-2020  Georg Umgiesser
!    Copyright (C) 2016  Erik Pascolo
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

! reading and interpolation of external files
!
! contents :
!
! revision log :
!
! 29.10.2012	ggu	created from scratch
! 05.11.2012	ggu	changed VERS_6_1_60
! 17.12.2012	ggu	changed VERS_6_1_61a
! 25.01.2013	ggu	changed VERS_6_1_62
! 03.05.2013	ggu	changed VERS_6_1_63
! 13.06.2013	ggu	changed VERS_6_1_65
! 17.06.2013	ggu	do not pass function into subroutine
! 18.06.2014	ggu	changed VERS_6_1_77
! 27.06.2014	ggu	changed VERS_6_1_78
! 02.07.2014	ggu	new framework finished
! 10.07.2014	ggu	only new file format allowed
! 18.07.2014	ggu	changed VERS_7_0_1
! 15.01.2015	ggu	changed VERS_7_1_1
! 26.02.2015	ggu	changed VERS_7_1_5
! 21.05.2015	ggu	changed VERS_7_1_11
! 13.07.2015	ggu	changed VERS_7_1_51
! 20.11.2015	ggu	changed VERS_7_3_15
! 22.02.2016	ggu&erp	new files for generic tracer (nvar>1)
! 06.06.2016	ggu	tracer_file routines changed
! 10.06.2016	ggu	changed VERS_7_5_13
! 09.09.2016	ggu	changed VERS_7_5_17
! 25.02.2018	ggu	file cleaned - time is now double
! 03.04.2018	ggu	changed VERS_7_5_43
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 14.02.2020	ggu	new routine ts_file_exists()
! 04.03.2020	ggu	iunit converted to id
! 05.04.2022	ggu	in tracer_init_file only set existing layers
! 27.10.2022	ggu	tracer_init_file also working for 2d arrays
! 22.09.2026	ggu	new routines vel_* and scalar_*
!
!*******************************************************************	
!*******************************************************************	
!*******************************************************************	
! T/S routines (nvar == 1, on nodes) *******************************
!*******************************************************************	
!*******************************************************************	
!*******************************************************************	

	subroutine ts_file_open(file,dtime,np,nlv,id)

! opens T/S file

	use intp_fem_file

	implicit none

	character*(*), intent(in) :: file	!name of file
	double precision, intent(in) :: dtime	!initial time
	integer, intent(in) :: np		!number of points expected
	integer, intent(in) :: nlv		!vertical dimmension
	integer, intent(out) :: id		!unit number (return)

	integer nvar,nexp,lexp,nintp
	integer nodes(1)
	real vconst(1)

	nvar = 1
	nexp = np
	lexp = nlv
	nintp = 2
	nodes = 0
	vconst = 0.

!$OMP CRITICAL
	call iff_init(dtime,file,nvar,nexp,lexp,nintp &
     &                                  ,nodes,vconst,id)
!$OMP END CRITICAL

	end

!*******************************************************************	

	subroutine ts_next_record(dtime,id,nlvddi,nkn,nlv,value)

        use levels, only : ilhkv
	use intp_fem_file

	implicit none

	double precision, intent(in) :: dtime
	integer, intent(in) :: id
	integer, intent(in) :: nlvddi
	integer, intent(in) :: nkn
	integer, intent(in) :: nlv
	real, intent(out) :: value(nlvddi,nkn)

	integer ldim,ndim,ivar
	integer k,lmax
        real vmin,vmax
	character*80 string

!--------------------------------------------------------------
! read new data
!--------------------------------------------------------------

	ivar = 1
	ndim = nkn
	ldim = nlvddi

	!write(6,*)'reading T/S values: ',dtime

	call iff_read_and_interpolate(id,dtime)
	call iff_time_interpolate(id,dtime,ivar,ndim,ldim,value)

	do k=1,nkn
	  lmax = ilhkv(k)
	  value(lmax+1:nlvddi,k) = 0
	end do

!--------------------------------------------------------------
! some statistics
!--------------------------------------------------------------

        call conmima(nlvddi,value,vmin,vmax)

        !write(6,*) 'min/max: ',vmin,vmax

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	end

!*******************************************************************	

	subroutine ts_file_close(id)

! closes T/S file

	use intp_fem_file

	implicit none

	integer, intent(in) :: id

	call iff_forget_file(id)

	end

!*******************************************************************	

	subroutine ts_file_descrp(id,name)

	use intp_fem_file

	implicit none

	integer, intent(in) :: id
	character*(*), intent(in) :: name

	call iff_set_description(id,0,name)

	end

!*******************************************************************	

	subroutine ts_file_exists(file,bexist)

! checks if file exists (and no read error)

	use intp_fem_file

	implicit none

	character*(*), intent(in) :: file		!name of file
	logical, intent(out) :: bexist

	call iff_file_exists(file,bexist)

	end

!*******************************************************************	
!*******************************************************************	
!*******************************************************************	
!****** tracer routines (nvar can be greater than 1) ***************	
!*******************************************************************	
!*******************************************************************	
!*******************************************************************	

	subroutine tracer_file_open(file,dtime,nvar,np,nlv,val0,id)

! opens tracer file

	use intp_fem_file

	implicit none

	character*(*), intent(in) :: file	!name of file
	double precision, intent(in) :: dtime	!time
	integer, intent(in) :: nvar		!number of (state) variables
	integer, intent(in) :: np		!number of points expected
	integer, intent(in) :: nlv		!number of vertical levels
	real, intent(in) :: val0(nvar)		!default initial condition
	integer, intent(out) :: id		!id of file (return)

	integer nexp,lexp,nintp
	integer nodes(nvar)

	nexp = np
	lexp = nlv
	nintp = 2
	nodes = 0

!$OMP CRITICAL
	call iff_init(dtime,file,nvar,nexp,lexp,nintp &
     &                                  ,nodes,val0,id)
!$OMP END CRITICAL

	if( id <= 0 ) then
	  write(6,*) 'Cannot open file: ',file
	  stop 'error stop tracer_file_open: error file open'
	end if

	end

!*******************************************************************	

	subroutine tracer_file_next_record(dtime,id &
     &					,nvar,nlvddi,nkn,nlv,value)

! reads next record of tracer

	use intp_fem_file

	implicit none

	double precision, intent(in) :: dtime
	integer, intent(in) :: id
	integer, intent(in) :: nvar
	integer, intent(in) :: nlvddi
	integer, intent(in) :: nkn
	integer, intent(in) :: nlv
	real, intent(out) :: value(nlvddi,nkn,nvar)

	integer ldim,ndim,ivar
	character*80 string

!--------------------------------------------------------------
! read new data
!--------------------------------------------------------------

	ndim = nkn
	ldim = nlvddi

	call iff_read_and_interpolate(id,dtime)
	do ivar=1,nvar
	  call iff_time_interpolate(id,dtime,ivar,ndim,ldim &
     &					,value(:,:,ivar))
	end do

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	end

!*******************************************************************	

	subroutine tracer_file_close(id)

! closes tracer file

	use intp_fem_file

	implicit none

	integer, intent(in) :: id

	call iff_forget_file(id)

	end

!*******************************************************************	

	subroutine tracer_file_descrp(id,text)

! sets description for file

	use intp_fem_file

	implicit none

	integer, intent(in) :: id
	character*(*), intent(in) :: text

	call iff_set_description(id,0,text)

	end

!*******************************************************************	

	subroutine tracer_file_exists(file,bexist)

! checks if file exists (and no read error)

	use intp_fem_file

	implicit none

	character*(*), intent(in) :: file		!name of file
	logical, intent(out) :: bexist

	call iff_file_exists(file,bexist)

	end

!*******************************************************************	
!*******************************************************************	
!*******************************************************************	
!****** velocity routines (nvar == 2, on elements) *****************	
!*******************************************************************	
!*******************************************************************	
!*******************************************************************	

	subroutine vel_file_open(file,dtime,np,nlv,id)

! opens vel file

	use intp_fem_file

	implicit none

	character*(*), intent(in) :: file	!name of file
	double precision, intent(in) :: dtime	!time
	integer, intent(in) :: np		!number of points expected
	integer, intent(in) :: nlv		!number of vertical levels
	integer, intent(out) :: id		!id of file (return)

	integer nvar,nexp,lexp,nintp
	integer nodes(1)
	real val0(2)				!default initial condition

	nvar = 2
	nexp = np
	lexp = nlv
	nintp = 2
	nodes = 0
	val0 = 0.

!$OMP CRITICAL
	call iff_init(dtime,file,nvar,nexp,lexp,nintp &
     &                                  ,nodes,val0,id)
!$OMP END CRITICAL

	if( id <= 0 ) then
	  write(6,*) 'Cannot open file: ',file
	  stop 'error stop vel_file_open: error file open'
	end if

	end

!*******************************************************************	

	subroutine vel_file_next_record(dtime,id &
     &					,nvar,nlvddi,nel,nlv,value)

! reads next record of vel

	use intp_fem_file

	implicit none

	double precision, intent(in) :: dtime
	integer, intent(in) :: id
	integer, intent(in) :: nvar
	integer, intent(in) :: nlvddi
	integer, intent(in) :: nel
	integer, intent(in) :: nlv
	real, intent(out) :: value(nlvddi,nel,nvar)

	integer ldim,ndim,ivar
	character*80 string

!--------------------------------------------------------------
! read new data
!--------------------------------------------------------------

	ndim = nel
	ldim = nlvddi

	call iff_read_and_interpolate(id,dtime)
	do ivar=1,nvar
	  call iff_time_interpolate(id,dtime,ivar,ndim,ldim &
     &					,value(:,:,ivar))
	end do

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	end

!*******************************************************************	

	subroutine vel_file_close(id)

! closes vel file

	use intp_fem_file

	implicit none

	integer, intent(in) :: id

	call iff_forget_file(id)

	end

!*******************************************************************	

	subroutine vel_file_descrp(id,text)

! sets description for file

	use intp_fem_file

	implicit none

	integer, intent(in) :: id
	character*(*), intent(in) :: text

	call iff_set_description(id,0,text)

	end

!*******************************************************************	

	subroutine vel_file_exists(file,bexist)

! checks if file exists (and no read error)

	use intp_fem_file

	implicit none

	character*(*), intent(in) :: file		!name of file
	logical, intent(out) :: bexist

	call iff_file_exists(file,bexist)

	end

!*******************************************************************	
!*******************************************************************	
!*******************************************************************	
! scalar routines (nvar == 1, on nodes or elements) ****************
!*******************************************************************	
!*******************************************************************	
!*******************************************************************	

	subroutine scalar_file_open(file,dtime,np,nlv,id)

! opens scalar file

	use intp_fem_file

	implicit none

	character*(*), intent(in) :: file	!name of file
	double precision, intent(in) :: dtime	!initial time
	integer, intent(in) :: np		!number of points expected
	integer, intent(in) :: nlv		!vertical dimmension
	integer, intent(out) :: id		!unit number (return)

	integer nvar,nexp,lexp,nintp
	integer nodes(1)
	real vconst(1)

	nvar = 1
	nexp = np
	lexp = nlv
	nintp = 2
	nodes = 0
	vconst = 0.

!$OMP CRITICAL
	call iff_init(dtime,file,nvar,nexp,lexp,nintp &
     &                                  ,nodes,vconst,id)
!$OMP END CRITICAL

	end

!*******************************************************************	

	subroutine scalar_next_record(dtime,id,nlvddi,np,nlv,value)

	use intp_fem_file

	implicit none

	double precision, intent(in) :: dtime
	integer, intent(in) :: id
	integer, intent(in) :: nlvddi
	integer, intent(in) :: np
	integer, intent(in) :: nlv
	real, intent(out) :: value(nlvddi,np)

	integer ldim,ndim,ivar
	integer k,lmax
        real vmin,vmax
	character*80 string

!--------------------------------------------------------------
! read new data
!--------------------------------------------------------------

	ivar = 1
	ndim = np
	ldim = nlvddi

	call iff_read_and_interpolate(id,dtime)
	call iff_time_interpolate(id,dtime,ivar,ndim,ldim,value)

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	end

!*******************************************************************	
!*******************************************************************	
!*******************************************************************	
! generic routines *************************************************
!*******************************************************************	
!*******************************************************************	
!*******************************************************************	

	subroutine generic_file_close(id)

! closes file

	use intp_fem_file

	implicit none

	integer, intent(in) :: id

	call iff_forget_file(id)

	end

!*******************************************************************	

	subroutine generic_file_descrp(id,text)

! sets description for file

	use intp_fem_file

	implicit none

	integer, intent(in) :: id
	character*(*), intent(in) :: text

	call iff_set_description(id,0,text)

	end

!*******************************************************************	

	subroutine generic_file_exists(file,bexist)

! checks if file exists (and no read error)

	use intp_fem_file

	implicit none

	character*(*), intent(in) :: file		!name of file
	logical, intent(out) :: bexist

	call iff_file_exists(file,bexist)

	end

!*******************************************************************	
!*******************************************************************	
!*******************************************************************	
! special routines *************************************************
!*******************************************************************	
!*******************************************************************	
!*******************************************************************	

	subroutine tracer_file_init(what,file_init,dtime &
     &				,nvar,nlvddi,nlv,nkn,val0,val)

! initialization of tracer from file

        use levels, only : ilhkv

        implicit none

	character*(*) what
	character*(*) file_init
        double precision dtime
        integer nvar
        integer nlvddi
        integer nlv
        integer nkn
        real val0(nvar)			!default for vals if no file is given
        real val(nlvddi,nkn,nvar)

        integer id,iv
	integer k,l,lmax
        character*80 file

        call getfnm(file_init,file)

	do iv=1,nvar
	  do k=1,nkn
            lmax = min(nlvddi,ilhkv(k))
            do l=1,lmax
              val(l,k,iv) = val0(iv)
	    end do
	  end do
	end do

        if( file == ' ' ) return

        write(6,*) 'tracer_init: opening file for ',trim(what)
        write(6,*) '   file name: ',trim(file)
        write(6,*) '   variables: ',nvar

        call tracer_file_open(file,dtime,nvar,nkn,nlv,val0,id)
        call tracer_file_descrp(id,what)
        call tracer_file_next_record(dtime,id,nvar,nlvddi,nkn,nlv,val)
        call tracer_file_close(id)

	write(6,*) 'tracer_init: successful init for ',trim(what)

	end

!*******************************************************************	
!*******************************************************************	
!*******************************************************************	

