
!--------------------------------------------------------------------------
!
!    Copyright (C) 1997-2012,2014-2020  Georg Umgiesser
!    Copyright (C) 2008,2010,2014  Christian Ferrarin
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

! time management routines
!
! contents :
!
!
! revision log :
!
! 23.09.1997	ggu	boundn deleted -> no access to data structure
! 20.03.1998	ggu	minor changes to priout
! 29.04.1998	ggu	new module for semi-implicit time-step
! 07.05.1998	ggu	check for error on return of nrdvecr
! 19.06.1998	ggu	version number is character
! 22.01.1999	ggu	oxygen section added
! 26.01.1999	ggu	new comp3d added
! 11.08.1999	ggu	new compatibility array hlhv initialized
! 19.11.1999	ggu	new routines for section vol
! 20.01.2000	ggu	common block /dimdim/ eliminated
! 04.02.2000	ggu	no priout, dobefor/after, pritime, endtime
! 15.05.2000	ggu	hm3v substituted
! 26.05.2000	ggu	copright statement adjourned
! 21.11.2001	ggu	routines to handle advective index (aix)
! 27.11.2001	ggu	routine to handle info file (getinfo)
! 11.10.2002	ggu	aix routines deleted
! 07.02.2003	ggu	routine added: changeimp, getaz; deleted getaza
! 10.08.2003	ggu	call adjust_chezy instead sp135r
! 14.08.2003	ggu	femver transfered to subver, not called in nlsh2d
! 20.08.2003	ggu	tsmed substituted by ts_shell
! 01.09.2003	ggu	call wrousa
! 03.09.2004	ggu	call admrst, comp3d renamed to init_3d (not used)
! 03.09.2004	ggu	nlv, hlv initialized in nlsh2d (FIXME)
! 28.09.2004	ggu	read lagrangian section
! 01.12.2004	ggu	new routine set_timestep for variable time step
! 17.01.2005	ggu	get_stab_index to newcon.f, error stop in set_timestep
! 14.03.2005	ggu	syncronize idt with end of simulation (set_timestep)
! 07.11.2005	ggu	handle new section sedtr for sediments
! 23.03.2006	ggu	changed time step to real
! 23.05.2007	ggu	recall variable time step pars at every time step
! 02.10.2007	ggu	bug fix in set_timestep for very small rindex
! 10.04.2008	ccf	output in netcdf format
! 28.04.2008	ggu	in set_timestep new call to advect_stability()
! 03.09.2008	ggu	in nlsh2d different error message
! 20.11.2008	ggu	init_3d deleted, nlv initialized to 0
! 18.11.2009	ggu	new format in pritime (write also time step)
! 22.02.2010	ggu	new call to hydro_stability to compute time step
! 22.02.2010	ccf	new routine for tidal pot. (tideforc), locaus deleted
! 26.02.2010	ggu	in set_timestep compute and write ri with old dt
! 22.03.2010	ggu	some comments for better readability
! 29.04.2010	ggu	new routine set_output_frequency() ... not finished
! 04.05.2010	ggu	shell to compute energy
! 22.02.2011	ggu	in pritime() new write to terminal
! 20.05.2011	ggu	changes in set_timestep(), element removal, idtmin
! 31.05.2011	ggu	changes for BFM
! 01.06.2011	ggu	idtmin introduced
! 12.07.2011	ggu	new routine next_output(), revised set_output_frequency
! 14.07.2011	ggu	new routines for original time step
! 13.09.2011	ggu	better error check, rdtitl() more robust
! 23.01.2012	ggu	new section "proj"
! 24.01.2012	ggu	new routine setup_parallel()
! 10.02.2012	ggu	new routines to initialize and access time common block
! 05.03.2014	ggu	code prepared to repeat time step (irepeat) - not ready
! 05.03.2014	ggu	new routines get_last/first_time()
! 10.04.2014	ccf	new section "wrt" for water renewal time
! 29.10.2014	ggu	do_() routines transfered from newpri.f
! 10.11.2014	ggu	time management routines transfered to this file
! 26.11.2014	ggu	changed VERS_7_0_7
! 05.12.2014	ggu	changed VERS_7_0_8
! 12.12.2014	ggu	changed VERS_7_0_9
! 19.12.2014	ggu	accept date also as string
! 23.12.2014	ggu	fractional time step introduced
! 07.01.2015	ggu	fractional time step without rounding (itsplt=3)
! 26.02.2015	ggu	changed VERS_7_1_5
! 30.04.2015	ggu	changed VERS_7_1_9
! 23.09.2015	ggu	time step is now working with dt as double
! 10.10.2015	ggu	use bsync as global to check for syncronization
! 23.09.2016	ggu	cleaned set_timestep()
! 30.09.2016	ggu	changed VERS_7_5_18
! 20.10.2017	ggu	new get_absolute_act_time(),get_absolute_ref_time()
! 04.11.2017	ggu	changed VERS_7_5_34
! 05.12.2017	ggu	changed VERS_7_5_39
! 22.02.2018	ggu	changed VERS_7_5_42
! 23.02.2018	ggu	most parts converted from int to double
! 29.03.2018	ggu	bug fix for syncronization step
! 03.04.2018	ggu	changed VERS_7_5_43
! 03.04.2018	ggu	changed VERS_7_5_44
! 13.04.2018	ggu	hydro_stability includes explicit gravity wave
! 13.04.2018	ggu	set_timestep now is working with mpi
! 16.04.2018	ggu	write warning if time step is over recommended one
! 09.02.2019	ggu	bug fix for syncronization of last time step
! 16.02.2019	ggu	changed VERS_7_5_60
! 21.05.2019	ggu	changed VERS_7_5_62
! 15.09.2019	ggu	small changes to account for synchorization time step
! 08.02.2020	ggu	utilities in this new file
! 16.02.2020	ggu	itunit eliminated
! 16.02.2020	ggu	new routines get_time_iterations(), get_ddt()
! 03.04.2020	ggu	new routine get_real_time()
! 18.05.2022	ggu	new routines cpu_time_*()
! 24.10.2025	ggu	new routine is_time_absolute()
!
!**********************************************************************
!**********************************************************************
!**********************************************************************

	module femtime

	implicit none

	integer, save :: itanf,itend,idt,nits,niter,it

        integer, save :: itunit,idtorig

	double precision, save :: t_act,dt_act,dt_orig,atime0,dtanf,dtend
	double precision, parameter :: d1000 = 1000*365*86400

	logical, save :: bsync

	character*20, save :: aline_act

	integer, parameter :: ncpu = 10
	double precision, save :: cputime(ncpu)
	double precision, save :: acutime(ncpu)

	end module femtime

!**********************************************************************

	function is_time_absolute(dtime)

! true if in initialization phase

	use femtime

	implicit none

	logical is_time_absolute
	double precision dtime

	is_time_absolute = dtime > d1000	! > 1000 years

	end

!**********************************************************************

	subroutine is_time_first(bfirst)

! true if in initialization phase

	use femtime

	implicit none

	logical bfirst

	bfirst = t_act .eq. dtanf

	end

!**********************************************************************

	subroutine is_time_last(blast)

! true if in last time step

	use femtime

	implicit none

	logical blast

	blast = t_act .eq. dtend

	end

!**********************************************************************

        subroutine get_act_dtime(dtact)

! returns actual time (double)

	use femtime

        implicit none

	double precision dtact

	dtact = t_act

	end

!**********************************************************************

        subroutine get_timeline(dtime,aline)

! returns time as string

	use femtime

        implicit none

	double precision dtime
	character*(*) aline

	double precision atime

	atime = atime0 + dtime
	call dts_format_abs_time(atime,aline)

	end

!**********************************************************************

        subroutine get_act_timeline(aline)

! returns actual time as string

	use femtime

        implicit none

	character*(*) aline

	aline = aline_act

	end

!**********************************************************************

        subroutine get_absolute_act_time(atime)

! returns actual time

	use femtime

        implicit none

	double precision atime

	atime = t_act + atime0

	end

!**********************************************************************

        subroutine get_absolute_ref_time(atime_ref)

! returns actual time

	use femtime

        implicit none

	double precision atime_ref

	atime_ref = atime0

	end

!**********************************************************************

        subroutine get_passed_dtime(dtime)

! returns time passed since start of simulation

	use femtime

        implicit none

	double precision dtime

	dtime = t_act - dtanf

	end

!**********************************************************************

        subroutine get_first_dtime(dtime)

! returns first (initial) time

	use femtime

        implicit none

	double precision dtime

	dtime = dtanf

	end

!**********************************************************************

        subroutine get_last_dtime(dtime)

! returns end time

	use femtime

        implicit none

	double precision dtime

	dtime = dtend

	end

!**********************************************************************

        subroutine get_ddt(ddt)

! returns time step (in seconds - double version)

	use femtime

        implicit none

	double precision ddt		!time step (return)

	ddt = dt_act

	end

!**********************************************************************

        subroutine get_timestep(dt)

! returns time step (in seconds - real version)

	use femtime

        implicit none

	real dt		!time step (return)

	dt = dt_act

	end

!**********************************************************************

        subroutine get_orig_timestep(dt)

! returns original real time step (in real seconds)

	use femtime

        implicit none

	real dt		!time step (return)

	dt = idtorig

	end

!**********************************************************************

        subroutine get_time_iterations(nit_done,nit_todo)

! returns iterations, already done and still to do

	use femtime

        implicit none

	integer nit_done
	integer nit_todo

	nit_done = niter
	nit_todo = nits

	end

!********************************************************************
!********************************************************************
!********************************************************************

	subroutine convert_to_dtime(aline,dtime)

	implicit none

	character*20 aline
	double precision dtime

	double precision atime,atime0

	call convert_to_atime(aline,atime)

        call get_absolute_ref_time(atime0)
        dtime = atime - atime0

	end subroutine convert_to_dtime

!********************************************************************

	subroutine convert_to_atime(aline,atime)

	use iso8601

	implicit none

	character*20 aline
	double precision atime

	integer date,time,ierr

        call string2date(aline,date,time,ierr)
        if( ierr /= 0 ) stop 'error converting date'
        call dts_to_abs_time(date,time,atime)

	end subroutine convert_to_atime

!********************************************************************
!********************************************************************
!********************************************************************

        subroutine get_real_time(atime,aline)

! returns real time as string and absolute time

        implicit none

        double precision atime
        character*20 aline

        character*20 date,time,zone
        integer values(8)

        call date_and_time(date,time,zone,values)

        date = date(1:4)//'-'//date(5:6)//'-'//date(7:8)
        time = time(1:2)//':'//time(3:4)//':'//time(5:6)

        aline = trim(date)//'::'//trim(time)
	call convert_to_atime(aline,atime)

        end

!********************************************************************
!********************************************************************
!********************************************************************

        subroutine cpu_time_init

	use femtime

        implicit none

	cputime = 0.
	acutime = 0.

        end

!********************************************************************

        subroutine cpu_time_start(itime)

	use femtime

        implicit none

	integer itime

	real time

	if( itime < 1 ) stop 'error stop cpu_time_accum: itime'
	if( itime > ncpu ) stop 'error stop cpu_time_accum: itime'

	call cpu_time(time)
	cputime(itime) = time

	end

!********************************************************************

        subroutine cpu_time_end(itime)

	use femtime

        implicit none

	integer itime

	real time
	double precision dtime

	if( itime < 1 ) stop 'error stop cpu_time_accum: itime'
	if( itime > ncpu ) stop 'error stop cpu_time_accum: itime'

	call cpu_time(time)
	dtime = time - cputime(itime)
	acutime(itime) = acutime(itime) + dtime

        end

!********************************************************************

        subroutine cpu_time_get(itime,time)

	use femtime

        implicit none

	integer itime
	real time

	if( itime < 1 ) stop 'error stop cpu_time_accum: itime'
	if( itime > ncpu ) stop 'error stop cpu_time_accum: itime'

	time = acutime(itime)

        end

!********************************************************************
!********************************************************************
!********************************************************************

