
!--------------------------------------------------------------------------
!
!    Copyright (C) 1992-1993,1995,1997-2001,2003,2006  Georg Umgiesser
!    Copyright (C) 2009-2012,2014-2015,2017-2019  Georg Umgiesser
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

! routines to set up constants
!
! contents :
!
! subroutine cstinit            set up constants and parameters
! subroutine cstcheck           checks parameters read
! subroutine cstsetup           sets up modules
! subroutine cstfile(strfile)   reads files (str and bas)
!
! subroutine cktime		check time parameters
! subroutine ckdate		check date parameter
! subroutine ckcori		set coriolis parameter
!
! revision log :
!
! 05.08.1992	ggu	$$ibtyp3 - implementation of ibtyp=3
! 31.08.1992	ggu	$$impli - implicit time step
! 23.11.1992	ggu	$$ibtyp11 - implementation of ibtyp=11,51
! 27.10.1993	ggu	$$roger - implementation of ibtyp=70 (nsea)
! 07.04.1995	ggu	!$$conzfl - conz compared to iflag (bug)
! 07.04.1995	ggu	!$$baroc - impl. of baroclinic salt/temp (21/22)
! 02.06.1997	ggu	$$EXTINW - extinw changed to ipint
! 13.06.1997	ggu	!$$kranf - check if kranf <= krend
! 29.06.1997	ggu	no cstdim in file
! 29.04.1998	ggu	module for semi-implicit time-step in own routine
! 12.08.1998	ggu	new parameter dlat -> specify latitude for coriolis
! 03.09.1998	ggu	call bocche to adjust depth at Venice inlets
! 06.11.1998	ggu	call to huniqu to set up hkv and hev
! 22.01.1999	ggu	oxygen modules introduced
! 04.01.2000	ggu	cstset -> cstcheck, new cstsetup
! 24.10.2001	ggu	Write on use of Coriolis
! 14.08.2003	ggu	set depth values transferred to newini
! 23.03.2006	ggu	use routine set_timeunit() to set time unit
! 11.02.2009	ggu	in cstfile set fixed name for STR file
! 23.03.2010	ggu	changed v6.1.1
! 09.04.2010	ggu	changed v6.1.3
! 28.09.2010	ggu	changed VERS_6_1_11
! 15.07.2011	ggu	new call to read basin
! 04.11.2011	ggu	changed VERS_6_1_35
! 14.02.2012	ggu	changed VERS_6_1_44
! 07.03.2014	ggu	changed VERS_6_1_72
! 18.07.2014	ggu	changed VERS_7_0_1
! 20.10.2014	ggu	new routine ckdate()
! 10.11.2014	ggu	time and date routines transfered to subtime.f
! 26.11.2014	ggu	changed VERS_7_0_7
! 23.12.2014	ggu	changed VERS_7_0_11
! 05.06.2015	ggu	changed VERS_7_1_12
! 10.07.2015	ggu	changed VERS_7_1_50
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 30.07.2015	ggu	read str-file from command line
! 31.07.2015	ggu	changed VERS_7_1_84
! 09.05.2017	ggu	changed VERS_7_5_26
! 05.10.2017	ggu	cstfile called with file name
! 05.12.2017	ggu	changed VERS_7_5_39
! 19.04.2018	ggu	changed VERS_7_5_45
! 16.10.2018	ggu	changed VERS_7_5_50
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 21.05.2019	ggu	changed VERS_7_5_62
! 10.05.2024	ggu	read and process also grd files
! 17.09.2024	lrp&ggu	call sleep(1) after writing basin
! 16.11.2024	ggu	in handle_grid_read() use ibasin to decide what to do
! 03.12.2024    lrp     ww3 restored
!
!************************************************************************

	subroutine cstinit

! set up constants and parameters

	use pkonst
	use mkonst
	use befor_after

	implicit none

! parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call nlsinh
	call check_parameter_values('after nlsinh')

! default names %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call fnminh
	call check_parameter_values('after fnminh')

! /mkonst/ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	eps1=1.0E-5
	eps2=1.0E-6
	pi=3.141592653
	flag=-9988765.0
	high=1.0E+35
        higi=.21474e+10   ! +/- 2147483647

! /pkonst/ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! these parameters are set one more time in cstset (from STR file)

	grav=9.81
	fcor=0.
	dcor=0.
	dirn=0.
	rowass=1025.
	roluft=1.225

! other modules %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!	call inexta
	!call inflxa
	!call invola
	call inarea
	call inbnds

	call inclos
!	call inoxy	!oxygen
!	call inlgr	!float tracking

	call modules(M_INIT)

! other parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	call putpar('flag',flag)

	end

!*******************************************************************

	subroutine cstcheck

! sets up and checks parameters read

	use mod_bnd
	use mod_bound_geom
	use basin, only : nkn,nel,ngr,mbw
	use pkonst
	use mkonst
	use befor_after

	implicit none

	integer ibarcl
        real getpar
        integer nkbnd

! parameters read from STR file

	grav=getpar('grav')
	rowass=getpar('rowass')
	roluft=getpar('roluft')

	call setup_date	!set up and check date parameter
	call setup_time	!set up and check and correct time parameters

!	call ckexta	!extra output points
	!call ckflxa	!flux sections
	!call ckvola	!flux sections
	call ckarea	!chezy values

! re-allocate boundary arrays

        call mod_bound_geom_reinit(nkn,nrb)
	call mod_bnd_reinit(nbc)
	call ckbnds	!boundary conditions

	call ckclos
!	call ckoxy	!oxygen

	ibarcl=nint(getpar('ibarcl'))
	if( ibarcl == 0 ) then
	  call putpar('itemp',0.)
	  call putpar('isalt',0.)
	end if

	call modules(M_CHECK)

	call ckcori	!sets dlat for coriolis parameter

	end

!*******************************************************************

	subroutine cstsetup

! sets up modules

	use befor_after
	use mod_info_output

	implicit none

	logical bw

	!call flxini
	!call volini

	bw = print_not_quiet_once()

	if( bw ) write(6,*) 'cstsetup: setting up modules'
	call modules(M_SETUP)
	if( bw ) write(6,*) 'cstsetup: finished setting up modules'

	end

!********************************************************************

	subroutine cstfile(strfile)

! reads files (str and bas)

	implicit none

	character*(*) strfile

	integer nstr,nbas
	character*80 basnam

	integer idefbas,ifileo

!------------------------------------------------------------
! get STR parameter file name and open it
!------------------------------------------------------------

	if( strfile == ' ' ) then
	  write(6,*) 'Usage: shyfem str-file'
	  stop 'error stop cstfile: internal error (1)'
	else
	  nstr = ifileo(0,strfile,'form','old')
	  if( nstr <= 0 ) then
	    write(6,*) 'error opening STR file: ',trim(strfile)
	    stop 'error stop cstfile: no such STR file'
	  end if
	end if

!------------------------------------------------------------
! read STR parameter file
!------------------------------------------------------------

	call nlsh2d(nstr)
	if( nstr .ne. 5 ) close(nstr)

!------------------------------------------------------------
! get STR parameter file
!------------------------------------------------------------

        call getfnm('basnam',basnam)
	call handle_basin_read(basnam)

!------------------------------------------------------------
! end of routine
!------------------------------------------------------------

	end

!********************************************************************

	subroutine handle_basin_read(basnam)

	use basin
	use mod_info_output

	implicit none

	character*(*) basnam

	logical breadgrd
	logical bwrite
	integer nbas
	integer idefbas
	character*80 dir,name,ext
	character*80 grid_file,basin_file
	logical filex

	bwrite = print_not_quiet_once()

	if( bwrite ) write(6,*) 'reading basin: ',trim(basnam)

	call parse_file_name(basnam,dir,name,ext)

	grid_file = trim(dir) // trim(name) // '.grd'
	basin_file = trim(dir) // trim(name) // '.bas'

	if( ext == ' ' ) then			! no extension given
	  if( filex(basin_file) ) then
	    breadgrd = .false.
	  else if( filex(grid_file) ) then
	    breadgrd = .true.
	  else
	    write(6,*) '*** no basin or grid file found for: ',trim(basnam)
	    write(6,*) 'need either .bas or .grd file'
	    stop 'error stop handle_basin_read: no bas or grd file'
	  end if
	else if( ext == 'bas' ) then
	  if( .not. filex(basnam) ) then
	    write(6,*) '*** no such file: ',trim(basnam)
	    stop 'error stop handle_basin_read: no bas file'
	  end if
	  breadgrd = .false.
	else if( ext == 'grd' ) then
	  if( .not. filex(basnam) ) then
	    write(6,*) '*** no such file: ',trim(basnam)
	    stop 'error stop handle_basin_read: no grd file'
	  end if
	  breadgrd = .true.
	else
	  write(6,*) '*** unknown extension in basin file found: ',trim(ext)
	  write(6,*) 'extension must be .grd, .bas, or no extension'
	  stop 'error stop handle_basin_read: unknown extension'
	end if
	
	if( breadgrd ) then		!read grd and produce bas file
	  if( bwrite ) write(6,*) 'reading grid file: ',trim(grid_file)
	  call handle_grid_read(grid_file)
	  call sleep(1)
	end if

	if( bwrite ) write(6,*) 'preparing to read bas file: ',trim(basin_file)
        nbas = idefbas(basin_file,'old')
	call basin_read(nbas)
	if( bwrite ) write(6,*) 'finished reading bas file: ',trim(basin_file)
	close(nbas)

	!stop 'programmed stop in handle_basin_read'

	end

!********************************************************************

	subroutine handle_grid_read(grid_file)

	use mod_info_output

	implicit none

	character*(*) grid_file

	logical bauto,blinfo
	logical bopti,bnepart,brenumber
	logical bquiet,bsilent
	integer ibasin
	real eps_area
	character*80 grdpart

	real getpar

	ibasin = nint(getpar('ibasin'))

	bopti = .false.
	brenumber = .false.
	if( ibasin >= 0 ) brenumber = .true.
	if( ibasin == 1 ) bopti = .true.
	if( bopti ) brenumber = .true.
	
	eps_area = 0.
	bauto = .true.
	blinfo = .false.
	bnepart = .false.
	grdpart = ' '

	if( print_not_quiet_once() ) then
	  if( bopti ) then
	    write(6,*) 'reading grd file with optimzation: ',trim(grid_file)
	  else
	    write(6,*) 'reading grd file without optimzation: ',trim(grid_file)
	  end if
	end if

	bquiet = is_quiet()
	bsilent = is_silent()

	call shypre_sub(grid_file,bquiet,bsilent,bauto,blinfo &
     &                          ,bopti,brenumber,bnepart,eps_area,grdpart)

	end

!********************************************************************

	subroutine ckcori

! set coriolis parameter

	use pkonst

	implicit none

	real dlat

	real getpar

	call bas_get_geom(dcor,dirn)	! from basin

	dlat=getpar('dlat')

	if( dlat .le. 90. ) dcor = dlat
	call putpar('dlat',dcor)

	end

!********************************************************************

