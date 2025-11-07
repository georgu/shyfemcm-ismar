
!--------------------------------------------------------------------------
!
!    Copyright (C) 1994-1995,1997-2020  Georg Umgiesser
!    Copyright (C) 2008,2010,2014,2018  Christian Ferrarin
!    Copyright (C) 2016  William McKiver
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

! utility routines for shyfem main routine
!
! contents :
!
! subroutine prilog                     prints parameters for main
! subroutine priout(mode)               writes output files
! subroutine pritst(id)                 test output of constants and variables
! subroutine init_3d			sets up 3D vertical arrays
! subroutine nlsh2d(iunit)		read STR  parameter file for FE model
! subroutine rdtitl			reads title section
!
! subroutine impini		initializes parameters for semi-implicit time
! function bimpli(it)		checks if semi-implicit time-step is active
! function getimp		gets weight for semi-implicit time-step
! subroutine setimp(dtime,aweigh) sets parameters for semi-implicit time-step
!
! revision log :
!
! 20.01.1994	ggu	$$conz - impl. of concentration in bnd(12,.)
! 07.04.1995	ggu	$$baroc - impl. of baroclinic salt/temp (21/22)
! 01.06.1997	ggu	complete revision
! 23.09.1997	ggu	boundn deleted -> no access to data structure
! 18.03.1998	ggu	use variable section instead name
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
! 09.04.2010	ggu	changed v6.1.3
! 29.04.2010	ggu	new routine set_output_frequency() ... not finished
! 03.05.2010	ggu	changed VERS_6_1_8
! 04.05.2010	ggu	shell to compute energy
! 22.07.2010	ggu	changed VERS_6_1_9
! 28.09.2010	ggu	changed VERS_6_1_11
! 15.12.2010	ggu	changed VERS_6_1_14
! 22.02.2011	ggu	in pritime() new write to terminal
! 01.03.2011	ggu	changed VERS_6_1_20
! 14.04.2011	ggu	changed VERS_6_1_22
! 20.05.2011	ggu	changes in set_timestep(), element removal, idtmin
! 31.05.2011	ggu	changes for BFM
! 01.06.2011	ggu	idtmin introduced
! 12.07.2011	ggu	new routine next_output(), revised set_output_frequency
! 14.07.2011	ggu	new routines for original time step
! 13.09.2011	ggu	better error check, rdtitl() more robust
! 18.10.2011	ggu	changed VERS_6_1_33
! 23.01.2012	ggu	new section "proj"
! 24.01.2012	ggu	new routine setup_parallel()
! 10.02.2012	ggu	new routines to initialize and access time common block
! 21.03.2012	ggu	changed VERS_6_1_50
! 30.03.2012	ggu	changed VERS_6_1_51
! 01.06.2012	ggu	changed VERS_6_1_53
! 21.06.2012	ggu	changed VERS_6_1_54
! 26.06.2012	ggu	changed VERS_6_1_55
! 29.08.2012	ggu	changed VERS_6_1_56
! 03.05.2013	ggu	changed VERS_6_1_63
! 10.05.2013	ggu	changed VERS_6_1_64
! 05.12.2013	ggu	changed VERS_6_1_70
! 28.01.2014	ggu	changed VERS_6_1_71
! 05.03.2014	ggu	code prepared to repeat time step (irepeat) - not ready
! 05.03.2014	ggu	new routines get_last/first_time()
! 10.04.2014	ccf	new section "wrt" for water renewal time
! 05.05.2014	ggu	changed VERS_6_1_74
! 07.07.2014	ggu	changed VERS_6_1_79
! 18.07.2014	ggu	changed VERS_7_0_1
! 21.10.2014	ggu	changed VERS_7_0_3
! 29.10.2014	ggu	do_() routines transfered from newpri.f
! 10.11.2014	ggu	shyfem time management routines to new file subtime.f
! 26.11.2014	ggu	changed VERS_7_0_7
! 01.12.2014	ccf	handle new section waves for wave module
! 19.12.2014	ggu	changed VERS_7_0_10
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 26.02.2015	ggu	changed VERS_7_1_5
! 01.04.2015	ggu	changed VERS_7_1_7
! 05.05.2015	ggu	changed VERS_7_1_10
! 21.05.2015	ggu	changed VERS_7_1_11
! 10.07.2015	ggu	changed VERS_7_1_50
! 13.07.2015	ggu	changed VERS_7_1_51
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 24.09.2015	ggu	call initialization for irv before reading STR file
! 10.10.2015	ggu	changed VERS_7_3_2
! 05.11.2015	ggu	changed VERS_7_3_12
! 26.05.2016	ggu	new check for sections: count_sections()
! 16.06.2016	wmk	added check for section nonhyd 
! 31.03.2017	ggu	changed VERS_7_5_24
! 04.11.2017	ggu	changed VERS_7_5_34
! 14.11.2017	ggu	changed VERS_7_5_36
! 05.12.2017	ggu	changed VERS_7_5_39
! 07.12.2017	ggu	changed VERS_7_5_40
! 22.02.2018	ggu	changed VERS_7_5_42
! 03.04.2018	ggu	changed VERS_7_5_43
! 19.04.2018	ggu	changed VERS_7_5_45
! 26.04.2018	ggu	changed VERS_7_5_46
! 11.05.2018	ggu	semi.h deleted and substituted with module
! 06.07.2018	ggu	changed VERS_7_5_48
! 18.12.2018	ggu	changed VERS_7_5_52
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 21.05.2019	ggu	changed VERS_7_5_62
! 06.11.2019	ggu	femtime eliminated
! 16.02.2020	ggu	femtime finally eliminated
! 18.03.2020	ggu	admrst() substituted with rst_write_restart()
! 13.09.2024    lrp     iatm and coupling with atmospheric model
! 15.11.2024    ggu     double for energy introduced
!
!************************************************************

	subroutine prilog

! writes output to terminal or log file

	use shympi
	use simul
	use befor_after

	implicit none

	character*80 name
        integer nrb,nbc
        integer nkbnd,nbnds
	integer niter,nits
	double precision atime0,dtanf,dtend
	real dt
	character*20 aline

	if( .not. shympi_is_master() ) return

        nrb = nkbnd()
        nbc = nbnds()

	call getfnm('runnam',name)
	write(6,*)
	write(6,*) '     Name of run :'
	write(6,*) trim(name)

	write(6,*)
	write(6,*) '     Description of run :'
	write(6,*) trim(descrp)
	write(6,*)

	call get_absolute_ref_time(atime0)
	call get_first_dtime(dtanf)
	call get_last_dtime(dtend)
	call get_orig_timestep(dt)
	call get_time_iterations(niter,nits)

	call dts_format_abs_time(atime0+dtanf,aline)
	write(6,*) '     start time = ',aline
	call dts_format_abs_time(atime0+dtend,aline)
	write(6,*) '     end time   = ',aline
	write(6,*) '     time step  =   ',dt
	write(6,*) '     Iterations to go :',nits

	call getfnm('basnam',name)
	write(6,*)
	write(6,*) '     Name of basin :'
	write(6,*) trim(name)

	write(6,*)
	write(6,*) '     Description of basin :'
	call bas_info

	write(6,*) '     Description of boundary values :'
	write(6,*)
	write(6,*) '     nbc,nrb     :',nbc,nrb

	write(6,*)
	write(6,*) '     Values from parameter file :'
	write(6,*)

	call pripar(6)
	call check_parameter_values('prilog')

	call prbnds		!prints boundary info

	!call prexta		!prints extra points

	!call pretsa		!prints extra time series points

	call prflxa		!prints flux sections

	!call prvola		!prints flux sections

	call prarea		!prints chezy values

	call prclos		!prints closing sections

!	call proxy		!prints oxygen section

!	call prlgr		!prints float coordinates

	call modules(M_PRINT)

	write(6,*)
	write(6,1030)
	write(6,*)

	return
 1030   format(1x,78('='))
	end

!********************************************************************

	subroutine pritst(id)

! test output of all constants and variables
!
! id    identifier

	use basin
	use simul
	use befor_after

	implicit none

	integer id

	write(6,*)
	write(6,*) '============================================='
	write(6,*) '================ test output ================'
	write(6,*) '============================================='
	write(6,1) '================ id = ',id,' ================'
	write(6,*) '============================================='
	write(6,*)

	call bas_info

	write(6,*) '/descrp/'
	write(6,*) descrp

	!call tsexta

	!call tsetsa

	call tsflxa

	!call tsvola

	call tsarea

	call tsbnds

	call tsclos

!	call tsoxy	!oxygen

!	write(6,*) '/close/'	!deleted on 28.05.97

	call modules(M_TEST)
	
	write(6,*)
	write(6,*) '============================================='
	write(6,*) '============================================='
	write(6,*) '============================================='
	write(6,*)

	return
    1   format(1x,a,i6,a)
	end

!********************************************************************
!********************************************************************
!********************************************************************

	subroutine do_init

! to do before time loop

	use mod_trace_point
	use befor_after

	implicit none

	call trace_point('calling wrboxa')
	call wrboxa
	call trace_point('calling wrousa')
	call wrousa
	call trace_point('finished do_init')

	end

!*******************************************************************

	subroutine nlsh2d(iunit)

! read STR  parameter file for FE model
!
! iunit		unit number of file

	use levels
	use nls
	use befor_after

	implicit none

	integer iunit

!---------------------------------------------------------------
!---------------------------------------------------------------

	character*6 section,extra,last
	logical bdebug
	integer nsc,num,iline
!	integer nrdsec,nrdveci,nrdvecr
	integer nrdsec,nrdvecr
	character*80 vers,aline

	logical hasreadsec

	bdebug = .true.
	bdebug = .false.

        nlv = 0         !is initialized really only in adjust_levels

	nsc = 0
	if(iunit.le.0) goto 63
	last = ' '

	if( bdebug ) write(6,*) 'start reading STR file'

	call nrdini(iunit)
	call mod_irv_initialize

! read loop over sections %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	do while( nrdsec(section,num,extra) .ne. 0 )

		if( bdebug ) write(6,*) 'new section: ',section,num

		call setsec(section,num)		!remember section
		call count_sections(section)

		nsc = nsc + 1

		if(section.eq.'title') then
			call rdtitl
			if( nsc .ne. 1 ) goto 66
		else if(section.eq.'end') then
			goto 69
		else if(section.eq.'para') then
			call nrdins(section)
		else if(section.eq.'proj') then
			call nrdins(section)
		else if(section.eq.'extra') then
			call rdexta
		else if(section.eq.'extts') then
			!call rdetsa
			call section_deleted(section,'use section $extra')
		else if(section.eq.'area') then
			call rdarea
		else if(section.eq.'name') then
			call nrdins(section)
		else if(section.eq.'bound') then
			call rdbnds(num)
		else if(section.eq.'nudge') then
			call rdnudge(num)
		else if(section.eq.'float') then
			!call rdfloa(nfldin)
			call section_deleted(section,'use section $lagrg')
		else if(section.eq.'close') then
			call rdclos(num)
		else if(section.eq.'flux') then
			call rdflxa
		!else if(section.eq.'vol') then
		!	call rdvola
		else if(section.eq.'wrt') then		!water renewal time
                        call nrdins(section)
		else if(section.eq.'wind') then
			call section_deleted(section,'use wind file')
		else if(section.eq.'oxypar') then	!oxygen
			call nrdins(section)
		else if(section.eq.'oxyarr') then	!oxygen
			!call rdoxy
			call section_deleted(section,'use ecological model')
		else if(section.eq.'bfmsc')then        ! BFM ECO TOOL
                        call nrdins(section)
		else if(section.eq.'levels') then
			call read_hlv
                else if(section.eq.'lagrg')then
                        call nrdins(section)
                else if(section.eq.'sedtr')then         !sediment
                        call readsed
                else if(section.eq.'waves')then         !wave
                        call nrdins(section)
                else if(section.eq.'atm')then           !atmosphere
                        call nrdins(section)
                else if(section.eq.'mudsec')then        !fluid mud
                        call readmud			!ARON
                else if(section.eq.'nonhyd')then        !NH model
                        call nrdins(section)	
                else if(section.eq.'connec')then       !connectivity
                        call nrdins(section)	
		else					!try modules
			goto 97
			call modules(M_READ)
			if( .not. hasreadsec() ) then	!sec has been handled?
				goto 97			! -> no
			end if
		end if

		last = section
	end do		!loop over sections

	call count_sections(' ')	!check section count

	if( bdebug ) write(6,*) 'finished reading STR file'

! end of read %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	return
   63	continue
	write(6,*) 'Cannot read STR file on unit : ',iunit
	stop 'error stop : nlsh2d'
   66	continue
	write(6,*) 'section $title must be first section'
	stop 'error stop : nlsh2d'
   69	continue
	write(6,*) 'section $end is not allowed on its own'
	write(6,*) 'last section found was: ',last
	write(6,*) 'no section has been opened'
	call nls_return_line(aline,iline)
	write(6,*) trim(aline),'     line number = ',iline
	stop 'error stop : nlsh2d'
   77	continue
	if( nlv .eq. -1 ) then
	  write(6,*) 'read error in section $levels'
	  write(6,*) 'nlv,nlvdi: ',nlv,nlvdi
	  stop 'error stop nlsh2d: read error'
	else
	  write(6,*) 'dimension error in section $levels'
	  write(6,*) 'nlvdi = ',nlvdi,'   number of data read = ',-nlv
	  stop 'error stop nlsh2d: dimension error'
	end if
   97	continue
	write(6,*) 'Cannot handle section : ',section
	stop 'error stop nlsh2d: no such section'
	end

!************************************************************************

	subroutine section_deleted(section,compat)

	implicit none

	character*(*) section,compat

	write(6,*) 'the following section has been removed: '
	write(6,*) section
	write(6,*) 'please remove section from parameter file'
	write(6,*) 'and substitute with a compatible solution'
        if( compat .ne. ' ' ) then
          write(6,*) 'a compatible solution could be: ',compat
        end if

	stop 'error stop section_deleted: section has been removed'
	end

!************************************************************************

        subroutine read_hlv

	use levels
	use basin, only : nkn,nel,ngr,mbw
        use nls

        implicit none

        integer n,nlvddi

	call nls_init_section

        n = nls_read_vector()
        !call levels_init(nkn,nel,n)
	call levels_hlv_init(n)
	nlv = n
	call levels_get_dimension(nlvddi)
	if( n > nlvddi ) then
	  write(6,*) 'nlv,nlvddi: ',nlv,nlvddi
	  stop 'error stop read_hlv: dimension error'
	end if
        call nls_copy_real_vect(n,hlv)

	call nls_finish_section

        end subroutine read_hlv

!************************************************************************

	subroutine rdtitl

! reads title section

	use simul

	implicit none

	character*80 line,extra

	integer nrdlin,nrdsec
	integer num

	if( nrdlin(line) .eq. 0 ) goto 65
	line = adjustl(line)
	descrp=line
	call putfnm('title',line)

	if( nrdlin(line) .eq. 0 ) goto 65
	line = adjustl(line)
	call putfnm('runnam',line)

	if( nrdlin(line) .eq. 0 ) goto 65
	line = adjustl(line)
	call putfnm('basnam',line)

	!if( nrdlin(line) .gt. 0 ) goto 65
	if( nrdsec(line,num,extra) .eq. 0 ) goto 65
	if( line .ne. 'end' ) goto 64

	return
   64	continue
	write(6,*) 'error in section $title'
	stop 'error stop rdtitl: no end found'
   65	continue
	write(6,*) 'error in section $title'
	stop 'error stop rdtitl: cannot read title section'
	end

!**********************************************************************
!**********************************************************************
!**********************************************************************
!
! routines for handling semi-implicit time-step
!
! the only routines that should be necessary to be called are
! setimp(dtime,weight) and getazam(az,am)
!
! setimp sets the implicit parameter until dtime to weight
! getazam returns az,am with the actual weight
!
! usage: call setimp in a program that would like to change the
!	 weight for a limited time (closing sections etc...)
!	 call getazam when az,am are needed
!	 (getaz if only az is needed)
!
!**********************************************************************
!**********************************************************************
!**********************************************************************

!======================================================================
	module semi_implicit
!======================================================================

	double precision, save :: 	dtimpl
	real, save :: 			weight = 0.5
	logical, save :: 		binit = .false.

!======================================================================
	end module semi_implicit
!======================================================================

!**********************************************************************

	subroutine getaz(azpar)

! returns actual az

	use semi_implicit

	implicit none

	real azpar

	double precision dtime

	real ampar
	real getpar

	call get_act_dtime(dtime)
	azpar=getpar('azpar')
	ampar=0.			!dummy

	call changeimp(dtime,azpar,ampar)

	end

!**********************************************************************

	subroutine getazam(azpar,ampar)

! returns actual az,am

	use semi_implicit

	implicit none

	real azpar
	real ampar

	double precision dtime

	real getpar

	call get_act_dtime(dtime)
	azpar=getpar('azpar')
	ampar=getpar('ampar')

	call changeimp(dtime,azpar,ampar)

	end

!**********************************************************************

	subroutine changeimp(dtime,azpar,ampar)

! changes parameters for semi-implicit time-step if necessary

	use semi_implicit

	implicit none

	double precision dtime
	real azpar,ampar

	if( binit .and. dtime .le. dtimpl ) then
	  azpar = weight
	  ampar = weight
	end if

	end

!**********************************************************************

	subroutine setimp(dtime,aweigh)

! sets parameters for semi-implicit time-step

	use semi_implicit

	implicit none

	double precision dtime
	real aweigh

	dtimpl = dtime
	weight = aweigh

	write(6,*) 'implicit parameters changed: ',dtimpl,weight

	end

!**********************************************************************

	function getimp()

! gets weight for semi-implicit time-step

	use semi_implicit

	implicit none

	real getimp

	getimp = weight

	end

!**********************************************************************
!**********************************************************************
!**********************************************************************

        subroutine getinfo_0(iunit)

! gets unit of info filea 
!
! routine transferred to info_output.f90 in shyutilmpi

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
!**********************************************************************
!**********************************************************************

	subroutine setup_omp_parallel

	use shympi

	implicit none

	logical bm
	integer n,nomp
	real getpar
	logical openmp_is_parallel

	call openmp_get_max_threads(n)
	nomp = nint(getpar('nomp'))
	if( nomp > 0 ) then
	  nomp = min(nomp,n)
	else if( nomp == 0 ) then
	  nomp = 1
	else	!nomp < 0 ... use all threads
	  nomp = n
	end if
	call openmp_set_num_threads(nomp)
	call putpar('nomp',float(nomp))

	if( .not. shympi_is_master() ) return

	write(6,*) 'start of setup of parallel OMP threads'

	if( openmp_is_parallel() ) then
	  write(6,*) 'the program can run in OMP parallel mode'
	else
	  write(6,*) 'the program cannot run in OMP parallel mode'
	end if
	  
	write(6,*) 'maximum available OMP threads: ',n
	write(6,*) 'for simulation used OMP threads: ',nomp
	write(6,*) 'end of setup of parallel OMP threads'

	end

!********************************************************************
!********************************************************************
!********************************************************************

        subroutine total_energy

! writes info on total energy to info file

	use shympi
	use mod_info_output

	implicit none

	!real kenergy,penergy,tenergy,ksurf
	double precision kenergy,penergy,tenergy,ksurf
	character*20 aline
	character*80 format
	logical debug
	real array(4)

	integer, save :: iuinfo = 0

	if( .not. binfo ) return

	debug = .false.

	if( iuinfo .eq. 0 ) then
          call getinfo(iuinfo)  !unit number of info file
	end if

	call energ3d(kenergy,penergy,ksurf,-1)
	!call energ3d(kenergy,penergy,ksurf,0)

	array(1) = kenergy
	array(2) = penergy
	array(3) = ksurf
	array(1:3) = (/real(kenergy),real(penergy),real(ksurf)/)
	call shympi_array_reduce('sum',array(1:3))

	if( debug ) write(6,*) 'penergy: ',my_id,penergy
!	kenergy = shympi_sum(kenergy)
!	penergy = shympi_sum(penergy)
!	ksurf = shympi_sum(ksurf)
	if( debug ) write(6,*) 'penergy total: ',my_id,penergy

	tenergy = kenergy + penergy

!	if(shympi_is_master()) then
!	  call get_act_timeline(aline)
!	  write(iuinfo,1000) ' energy: ',aline &
!     &				,kenergy,penergy,tenergy,ksurf
! 1000	  format(a,a20,4e12.4)
!	end if
	
	array(4) = array(3)
	array(3) = array(1) + array(2)
	format = '(a,4e12.5)'
	call info_output(' energy','none',4,array,.true.,format)

	end

!********************************************************************
!********************************************************************
!********************************************************************

	subroutine count_sections(section)

	implicit none

	character*(*) section

	integer, parameter :: ndim = 100
	character*80, save :: sections(ndim) = ' '
	integer, save :: count(ndim) = 0
	integer, save :: ntot = 0

	integer i,ierr,ic

!----------------------------------------------------------
! check sections
!----------------------------------------------------------

	if( section == ' ' ) then	!check sections
	  ierr = 0
	  do i=1,ntot
	    ic = count(i)
	    if( ic > 1 ) then
	      if( sections(i) /= 'bound' ) cycle
	      if( sections(i) /= 'close' ) cycle
	      write(6,*) 'section ',trim(sections(i)),' - count = ',ic
	      ierr = 1
	    end if
	  end do
	  if( ierr > 0 ) then
	    stop 'error stop count_sections: not unique sections'
	  end if
	  return
	end if

!----------------------------------------------------------
! count section
!----------------------------------------------------------

	do i=1,ntot
	  if( section == sections(i) ) then
	    count(i) = count(i) + 1
	    return
	  end if
	end do

	ntot = ntot + 1
	if( ntot > ndim ) then
	  stop 'error stop count_sections: ndim'
	end if

!----------------------------------------------------------
! insert new section
!----------------------------------------------------------

	sections(ntot) = section
	count(ntot) = 1

!----------------------------------------------------------
! end of routine
!----------------------------------------------------------

	end

!********************************************************************
