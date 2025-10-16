
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998-2006,2008-2020  Georg Umgiesser
!    Copyright (C) 2006-2008,2014-2015,2019  Christian Ferrarin
!    Copyright (C) 2008  Andrea Cucco
!    Copyright (C) 2019  Petras Zemlys
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

! routines system dependent
!
! contents :
!
! subroutine nlsinh             initializes the hp parameter file
! subroutine nlsina             initializes the ap parameter file
! subroutine fnminh             initializes default names
!
! revision log :
!
! 22.01.1998	ggu	no more strdir, apndir -> read apnpar.str from local
! 18.03.1998	ggu	no '[not defined...]' in basdir...
! 20.03.1998	ggu	for dos systems .memory -> _memory
! 22.03.1998	ggu	some changes for Technital integrated
! 30.04.1998	ggu	some changes for flux sections
! 13.05.1998	ggu	added colfil (apncol.str)
! 25.05.1998	ggu	documentation added
! 26.05.1998	ggu	documentation for wind file added
! 19.06.1998	ggu	some useless comments deleted
! 24.06.1998	ggu	qflux (heat flux) added
! 22.07.1998	ggu	more on documentation
! 23.07.1998	ggu	documentation
! 12.08.1998	ggu	new parameter dlat -> specify latitude for coriolis
! 02.09.1998	ggu	hack: change depth in Venice inlets (hlido,...)
! 24.11.1998	ggu	switch for biological reactor introduced
! 12.02.1999	ggu	default of parameter file changed to apnstd.str
! 12.02.1999	ggu	new parameters for plotting
! 13.04.1999	ggu	new parameter itide
! 19.04.1999	ggu	itide changed to rtide
! 27.05.1999	ggu	icust introduced
! 31.05.1999	ggu	dval default changed
! 28.09.1999	ggu	new flag regflx
! 19.11.1999	ggu	new parameters for section vol
! 08.08.2000	ggu	hlvmin is now percentage of last layer thickness
! 03.12.2001	ggu	new parameters (rlhdif)
! 11.10.2002	ggu	rstot new meaning
! 14.10.2002	ggu	atpar,adpar,aapar -> default set to 1.0 (implicit)
! 05.10.2003	ggu	changes in color handling of post routines
! 04.12.2003	ggu	sediment and wave module integrated
! 05.03.2004	ggu	bio variable for initialization
! 03.09.2004	ggu	restart variables
! 22.09.2004	ggu	nlsina_3d(): always use 3d file for plot (level=0)
! 28.09.2004	ggu	lagrangian routines integrated (LAGR)
! 05.10.2004	ggu	some more documentation
! 02.12.2004	ggu	documentation for variable time step
! 06.12.2004	ggu	new subroutine nlsina_legvar
! 17.01.2005	ggu	new parameters for horizontal diffusion
! 15.03.2005	ggu	cleaned horiz. diff., read new section legvar
! 19.05.2005	ggu	added time for custom routine (tcust)
! 04.11.2005	ggu	new parameter itlin (semi-lagrangian), some corrections
! 07.11.2005	ggu	new parameter itvd (TVD)
! 07.11.2005	ggu	parameters deleted: isedi,sedref,sedgrs
! 08.11.2005	ggu	do not set old parameters, some in nlsina_para
! 08.11.2005	ggu	more documentation
! 16.02.2006	ggu	new flag itoxi and file toxi
! 23.03.2006	ggu	new variable itunit for unit of time step
! 12.09.2006	ggu	time and date introduced
! 28.09.2006	ggu	new iprogr for progress of simulation
! 18.10.2006	ccf	new params iwave and dtwave for wave model
! 15.11.2006	ggu	new parameters to construct stretched vert. coordinates
! 16.11.2006	ggu	use itemp,isalt to decide about advection
! 29.11.2006	ggu	rwhpar for horizontal diffusion in lagrangian model
! 27.08.2007	ccf	isphe = 1  for spherical coordinate system
! 07.04.2008	aac	parameters for ersem ecological model call
! 10.04.2008	ccf	netcdf and gotmpa parameters
! 17.04.2008	ggu	new param ievap (to compute evaporation mass flux)
! 28.04.2008	ggu	rstot deleted, contau introduced
! 29.04.2008	ggu&aac	new vars for ERSEM
! 16.06.2008	ggu	old parts deleted, new documentation
! 17.09.2008	ggu	new interpretation for level (-1 = bottom)
! 11.10.2008	ggu	zfranco added
! 10.11.2008	ggu	new variable ilytyp, cleaned
! 19.11.2008	ggu	new variable noslip
! 24.11.2008	ggu	new variable vreps (mass error)
! 06.12.2008	ggu	new variables nbfix and nsigma
! 09.01.2009	ggu	documentation
! 21.01.2009	ggu	new variable vrerr (stop if mass error)
! 29.01.2009	ggu	various changes, better documentation
! 13.03.2009	ggu	bugfix: some parameters had default section not set
! 07.05.2009	ggu	new parameter ityrst
! 15.06.2009	ggu	new parameters for plot: isphe, reggrd, reggry
! 11.09.2009	ggu	new section $sect
! 09.10.2009	ggu	new parameter sclvel
! 13.10.2009	ggu	documentation is $sect, new hvmax, lvmax
! 22.02.2010	ggu	new parameter itdrag
! 23.02.2010	ggu	new parameter regdst
! 23.03.2010	ggu	changed v6.1.1
! 26.03.2010	ggu	new parameters for arrows in section plot
! 28.09.2010	ggu	new value for icor
! 29.09.2010	ggu	new param vmode,rxscal,ryscal
! 15.12.2010	ggu	nsigma renamed to nbsig, nsigma used for sigma layers
! 16.12.2010	ggu	changed VERS_6_1_15
! 21.12.2010	ggu	new parameter rwscal
! 27.01.2011	ggu	changed VERS_6_1_17
! 16.02.2011	ggu	new default for isphe, new routine fnm_aquabc_init()
! 25.02.2011	ggu	new param wsmax to catch errors in wind type
! 01.03.2011	ggu	changed VERS_6_1_20
! 23.03.2011	ggu	new parameter itvdv
! 24.03.2011	ggu	new parameters iheat,hdecay,botabs
! 14.04.2011	ggu	changed VERS_6_1_22
! 01.06.2011	ggu	new parameter idtmin
! 14.07.2011	ggu	changed VERS_6_1_27
! 18.08.2011	ggu	new parameter isoinp (interpolate inside element)
! 18.09.2011	ggu	change default for isphe for output (-1)
! 18.10.2011	ggu	changed VERS_6_1_33
! 03.11.2011	ggu	new parameter hsigma (hybrid)
! 18.11.2011	ggu	new subroutine nlsinh_proj() for projection
! 24.01.2012	ggu	new parameter nomp
! 14.02.2012	ggu	changed VERS_6_1_44
! 19.03.2012	ggu	changed VERS_6_1_49
! 02.05.2012	ggu	new default for ndccol (-> 0)
! 01.06.2012	ggu	changed VERS_6_1_53
! 29.08.2012	ggu	changed VERS_6_1_56
! 12.09.2012	ggu	changed VERS_6_1_57
! 24.10.2012	ggu	new parameter dxmin
! 10.05.2013	ggu	new parameters idtbox,itmbox, more comments
! 10.05.2013	ggu	new parameter inohyd
! 16.05.2013	ggu	file name bound renamed to zinit
! 13.06.2013	ggu	changed VERS_6_1_65
! 12.09.2013	ggu	changed VERS_6_1_67
! 25.10.2013	ggu	changed VERS_6_1_68
! 28.01.2014	ggu	changed VERS_6_1_71
! 07.03.2014	ggu	changed VERS_6_1_72
! 28.03.2014	ggu	some new params for lagrangian
! 10.04.2014	ccf	new section "wrt" for water renewal time
! 05.05.2014	ggu	changed VERS_6_1_74
! 30.05.2014	ggu	new default for dragco, new metpnt
! 18.07.2014	ggu	changed VERS_7_0_1
! 20.10.2014	ggu	new default for date (-1)
! 30.10.2014	ggu	changed VERS_7_0_4
! 05.11.2014	ggu	changed VERS_7_0_5
! 26.11.2014	ggu	changed VERS_7_0_7
! 01.12.2014	ccf	wave parameters moved to subwave.f
! 23.12.2014	ggu	changed VERS_7_0_11
! 09.01.2015	ggu	changed VERS_7_0_12
! 15.01.2015	ggu	changed VERS_7_1_1
! 26.02.2015	ggu	changed VERS_7_1_5
! 01.04.2015	ggu	changed VERS_7_1_7
! 30.04.2015	ggu	changed VERS_7_1_9
! 06.05.2015	ccf	new parameters itmoff and offlin
! 21.05.2015	ggu	changed VERS_7_1_11
! 29.09.2015	ggu	new boundary file surfvel
! 29.09.2015	ggu	new initial file uvinit, new flgrst
! 10.10.2015	ggu	changed VERS_7_3_2
! 19.10.2015	ggu	changed VERS_7_3_6
! 22.10.2015	ggu	changed VERS_7_3_8
! 05.11.2015	ggu	changed VERS_7_3_12
! 16.11.2015	ggu	changed VERS_7_3_14
! 20.11.2015	ggu	changed VERS_7_3_15
! 16.12.2015	ggu	changed VERS_7_3_16
! 01.02.2016	ggu	some plot params shifted to para section (bbgray, etc.)
! 19.02.2016	ggu	changed VERS_7_5_2
! 22.02.2016	ggu	new file for initial conditions bfmini
! 05.04.2016	ggu	new parameter iaicef for ice free areas
! 15.04.2016	ggu	changed VERS_7_5_8
! 25.05.2016	ggu	changed VERS_7_5_10
! 07.06.2016	ggu	changed VERS_7_5_12
! 10.06.2016	ggu	changed VERS_7_5_13
! 14.06.2016	ggu	changed VERS_7_5_14
! 17.06.2016	ggu	changed VERS_7_5_15
! 27.06.2016	ggu	changed VERS_7_5_16
! 09.09.2016	ggu	changed VERS_7_5_17
! 30.09.2016	ggu	changed VERS_7_5_18
! 05.10.2016	ggu	changed VERS_7_5_19
! 12.01.2017	ggu	changed VERS_7_5_21
! 13.02.2017	ggu	changed VERS_7_5_23
! 31.03.2017	ggu	changed VERS_7_5_24
! 13.04.2017	ggu	changed VERS_7_5_25
! 09.05.2017	ggu	changed VERS_7_5_26
! 26.05.2017	ggu	default of ishyff is 1 - no nos or ous files
! 13.06.2017	ggu	changed VERS_7_5_29
! 11.07.2017	ggu	changed VERS_7_5_30
! 26.09.2017	ggu	changed VERS_7_5_32
! 14.11.2017	ggu	changed VERS_7_5_36
! 05.12.2017	ggu	changed VERS_7_5_39
! 07.12.2017	ggu	changed VERS_7_5_40
! 24.01.2018	ggu	changed VERS_7_5_41
! 22.02.2018	ggu	changed VERS_7_5_42
! 03.04.2018	ggu	changed VERS_7_5_43
! 19.04.2018	ggu	changed VERS_7_5_45
! 06.07.2018	ggu	changed VERS_7_5_48
! 16.10.2018	ggu	changed VERS_7_5_50
! 25.10.2018	ggu	changed VERS_7_5_51
! 18.12.2018	ggu	changed VERS_7_5_52
! 27.12.2018	ggu	changed VERS_7_5_54
! 12.02.2019	ccf	introduced ibstrs for computing the bottom stess
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	new variables for plotting basin
! 21.05.2019	ggu	changed VERS_7_5_62
! 04.07.2019	ggu	new description for WW3
! 27.09.2019	pzy	new variables for aquabc
! 22.10.2019	ggu	some update of documentation
! 16.02.2020	ggu	itunit not supported anymore
! 05.03.2020	ggu	documentation upgraded
! 09.03.2020	ggu	more documentation upgraded, spell check
! 27.03.2020	ggu	documentation on ibarcl and nudging
! 11.11.2020	ggu	new parameter idtice and icemod
! 30.03.2021	ggu	new parameters idtdbg,itmdbg
! 23.04.2021	clr	new parameters petsc_zcfg and amgx_zcfg
! 15.02.2022	ggu	new parameter iage
! 15.02.2022	ggu	new parameters for limiting_scalar
! 20.07.2023    lrp     new parameter nzadapt
! 13.09.2024    lrp     new parameter iatm
! 16.11.2024	ggu	new parameter ibasin
! 23.11.2024	ggu	new parameter pqual
! 03.12.2024    lrp     new parameter irain for the coupled model
! 01.04.2025    ggu     new value 5 for ibarcl
! 16.04.2025    ggu     new parameter wsettl
!
!************************************************************************

	subroutine nlsinh

! initializes the parameter file for the main FE model

	implicit none

	call nlsinh_general
	call nlsinh_lagrg
        call nlsinh_wrt
	call nlsinh_bfmsc
	call nlsinh_proj
	call nlsinh_undoc
	call nlsinh_georg
	call nlsinh_unused
	call nlsinh_waves
        call nlsinh_atm
	call nlsinh_nonhydro
	call nlsinh_connect

	end

!************************************************************************

	subroutine nlsinh_general

	use para

	implicit none

! $para section

	call sctpar('para')		!sets default section

! DOCS	START	S_para_h
!
! This section |$para| defines the general behavior of the simulation,
! gives various constants of parameters and determines what
! output files are written. In the following the meaning of
! all possible parameters is given.
! 
! Note that the only compulsory parameters in this section are
! the ones that chose the duration of the simulation and the
! integration time step. All other parameters are optional.
!
! DOCS	COMPULS		Compulsory time parameters
!
! These parameters are compulsory parameters that define the
! period of the simulation. They must be present in all cases.
!
! |itanf|	Start time of simulation. (Default 0)
! |itend|	End time of simulation.
! |idt|		Time step of integration.

	call addpar('itanf',0.)
	call addpar('itend',0.)
	call addpar('idt',0.)

!c------------------------------------------------------------------------

! DOCS	OUTPUT		Output parameters
!
! The following parameters deal with the output frequency
! and start time to external files. The content of the various
! output files should be looked up in the appropriate section.
!
! The default for the time step of output of the files is 0 which
! means that no output file is written. If the time step of the
! output files is equal to the time step of the simulation then
! at every time step the output file is written. The default start time
! of the output is |itanf|, the start of the simulation.
!
! |idtout|, |itmout|	Time step and start time for writing to file OUT,
!			the file containing the general hydrodynamic results.

	call addpar('idtout',0.)
	call addpar('itmout',-1.)

! |idtext|, |itmext|	Time step and start time for writing to file EXT,
!			the file containing hydrodynamic data of extra points.
!			The extra points for which the data is written
!			to this file are given in section |extra| of
!			the parameter file.

	call addpar('idtext',0.)
	call addpar('itmext',-1.)

! |idtrst|		Time step for writing the restart
!			file (extension RST). No restart file is written
!			with |idtrst| equal to 0. A negative value
!			is also possible for the time step. In this case
!			the time step used is |-idtrst|, but the file is
!			overwritten every time. It therefore contains 
!			always only the last written restart record. The
!			special value of |idtrst = -1| will write only the
!			last time step of the simulation in the restart file. 
!			This is useful if you want to start another
!			simulation from the last output. (Default 0)
! |itmrst|		Start time for writing the restart file. If
!			not given it is the beginning of the simulation.
! |itrst|		Time to use for the restart. If a restart
!			is performed, then the file name containing
!			the restart data has to be specified in |restrt|
!			and the time record corresponding to |itrst|
!			is used in this file. A value of -1 is also possible.
!			In this case the last record in the restart file
!			is used for the restart and the simulation starts
!			from this time. Be aware that a value of -1 changes
!			the parameter |itanf| to the time of the last
!			record found in |restrt|.
! |ityrst|		Type of restart. If 0 and the restart file is not
!			found the program will exit with an error. Otherwise
!			the program will simply continue with a cold start.
!			If |ityrst| is 1 and the given time record is not
!			found in the file it will exit with error. If
!			it is 2 it will initialize all values from the
!			first time record after |itrst|. Therefore, the
!			value of 2 will guarantee that the program will not
!			abort and continue running, but it might not
!			be doing what you intended. (Default 0)
! |flgrst|		This variable indicates which variables are read
!			from the restart file. By default all available
!			variables are read and used. If some variables
!			are not wanted (because, e.g., you want to restart
!			from a different T/S field), this fact can be indicated
!			in |flgrst|. 1 indicates restart of hydro values,
!			10 the depth values, 100 T/S values, 1000 the tracer
!			concentration, 10000 vertical velocities and
!			100000 the ecological variables. Therefore, a value
!			of 10111 indicates a restart of everything except
!			the tracer and the ecological values. The default
!			value for |flgrst| is -1, which means 111111.

	call addpar('idtrst',0.)
	call addpar('itmrst',-1.)
	call addpar('itrst',0.)
	call addpar('ityrst',0.)
	call addpar('flgrst',-1.)

! |idtres|, |itmres|	Time step and start time for writing to file RES,
!			the file containing residual hydrodynamic data.

	call addpar('idtres',0.)
	call addpar('itmres',-1.)

! |idtrms|, |itmrms|	Time step and start time for writing to file RMS,
!			the file containing hydrodynamic data of root mean
!			square velocities.

	call addpar('idtrms',0.)
	call addpar('itmrms',-1.)

! |idtflx|, |itmflx|	Time step and start time for writing to file FLX,
!			the file containing discharge data through defined
!			sections.
!			The transects for which the discharges are computed
!			are given in section |flux| of
!			the parameter file.

	call addpar('idtflx',0.)
	call addpar('itmflx',-1.)

! |idtstb|, |itmstb|	Time step and start time for writing a file *.stb.shy
!			with the stability index for debug reasons.
!			The file is only written if a scalar is computed,
!			such as S/T or tracer.

	call addpar('idtstb',0.)
	call addpar('itmstb',-1.)

! |idtsti|, |itmsti|	Time step and start time for writing a file *.sti.shy
!			with the maximum allowed time step for the simulation.

	call addpar('idtsti',0.)
	call addpar('itmsti',-1.)

! |idtdbg|, |itmdbg|	Time step and start time for writing a debug file
!			with values of various variables. In order to
!			use this feature you have to run |shyfem| with
!			the command line option |-debout|. Please be
!			aware that this feature will create extremely
!			big files. Use at you own risk.

	call addpar('idtdbg',0.)
	call addpar('itmdbg',-1.)

! |idtvol|, |itmvol|	Time step and start time for writing to file VOL,
!			the file containing volume information of areas
!			defined by transects.
!			The transects that are used to compute the volumes
!			are given in section |volume| of
!			the parameter file.

	call addpar('idtvol',0.)
	call addpar('itmvol',-1.)

! |netcdf|		This parameter chooses output in NetCDF format 
!			if |netcdf| is 1, else the format is unformatted
!			FORTRAN files. (Default 0)

	call addpar('netcdf',0.)

! |idtoff|	handles offline mode (default 0):
!		\begin{description}
!		\item[0] do nothing (no offline routines called)
!		\item[$>$0] write offline data file (.off) 
!			with time step |idtoff|.
!		\item[$<$0] reads offline data from file |offlin|
!			defined in section |name|. Usage:
!			\begin{description}
!			\item[-1] uses offline hydro results
!			\item[-2] uses offline T/S results
!			\item[-4] uses offline turbulence results
!			\end{description}
!			Combinations are possible. A value of -3 would
!			read hydro and T/S results.
!		\end{description}

	call addpar('idtoff',0.)

! |itmoff|	Start time for writing to file OFF,
!		the file containing data for offline runs.

	call addpar('itmoff',-1.)

! |idtmet|, |itmmet|	Time step and start time for writing meteo
!			variables read from file. 
! |imetout|		This parameters indicates what meteo parameters
!			should be output. For wind set it to 1,
!			for heat data 10,
!			for rain 100, for ice 1000, and for heat fluxes 10000. 
!			Combinations are possible, e.g., 11 writes wind 
!			and heat data, and 11111 writes all available data
!			to the file. (Default 0)

	call addpar('idtmet',0.)
	call addpar('itmmet',-1.)
	call addpar('imetout',0.)

!c------------------------------------------------------------------------

! DOCS	BASIN		Handling of grd files
!
! The model normally reads a file |.bas| that contains the numerical
! mesh in an unformatted format. This file can be produced by running
! |shypre| on a |.grd| file that contains the numerical grid in
! ASCII format. However, |shyfem| can also directly process a |.grd| file
! and transform it into a |.bas| file, eliminating the need for a
! pre-processing of the |.grd| file by |shypre|.
!
! If the basin name specified in the |STR| file has no extension, then
! |shyfem| will first look for a |.bas| file and read it. If not available,
! it will look for a |.grd| file, read and process it. If also a |.grd|
! file is not found, |shyfem| will exit with an error.
!
! If the basin name in the |STR| file has an extension |.bas|, 
! |shyfem| will only look for a |.bas| file, and, if not found, 
! it will exit with an error. If the basin file in the |STR| file has
! the extension |.grd|, |shyfem| will look for a |.grd| file and
! process it to create a |.bas| file, even if such a file already exists.
! If no file |.grd| is found, |shyfem| will exit with an error.
!
! If a |.grd| file has to be read, then the parameter |ibasin| will
! regulate, how the file will be processed. A value of -1 will just
! read the file and not change the numbering of nodes and element. 
! A value of 0 will just renumber the elements, and a value of 1 will
! optimize the bandwidth, renumbering both nodes and elements.
! These values correspond to the command line options of |shypre|
! as |-noopti| for -1, |-renumber| for 0, and |-opti| for 1.

! |ibasin|		Treatment of the |.grd| file. A value of -1
!			will leave the numbering of nodes and elements
!			as in the |.grd| file, a value of 0 just
!			renumbers the elements, and a value of 1
!			optimizes the bandwidth, renumbering both 
!			nodes and elements. (Default -1)

	call addpar('ibasin',-1.)

! The model, if run in MPI mode, also checks the quality of the partition
! that is automatically produced. If the quality is too bad the execution
! of the program stops.
!
! The quality index is an empirical value that gives an indication
! on the quality of the partition. Its value depends on the ratio
! of ghost elements to the total number of elements. Smaller numbers
! of the index are better.
!
! If the value of the index is higher than the value in |pqual| the
! program stops executing. If you are satisfied anyway with the partition
! produced, just set the value of |pqual| higher than the computed
! quality index and |shyfem| will run. i
!
! The partition can be checked with either
! |grid -FT partition.np.node.grd| or also |grid -FTo partition.np.node.grd|
! where np is the number of the desired domains.
!
! It is important to note that the
! quality index does not concern the quality of the results, just
! the speed of execution of the program.

! |pqual|		Maximum quality index allowed. (Default 2.5)

	call addpar('pqual',2.5)

!c------------------------------------------------------------------------

! DOCS	TIME_DATE	General time and date parameters
!
! A time and date can be assigned to the simulation. These values
! refer to the time 0 of the FEM model. The format for the date is
! YYYYMMDD and for the time HHMMSS.
! You can also give a time zone if your time is not referring to
! GMT but to another time zone such as MET.

! |date|                The real date corresponding to time 0. (Default 0)
! |time|                The real time corresponding to time 0. (Default 0)
! |tz|                  The time zone you are in. This is 0 for GMT, 1 for MET
!                       and 2 for MEST (MET summer time). (Default 0)

        call addpar('date',-1.)
        call addpar('time',0.)
        call addpar('tz',0.)

!c------------------------------------------------------------------------

! DOCS	TERMS		Model parameters
!
! The next parameters define the inclusion or exclusion of
! certain terms of the primitive equations.
!
! |ilin|	Linearization of the momentum equations. If |ilin| 
!		is different from 0 the advective terms are not 
!		included in the computation. (Default 1)
! |itlin|	This parameter decides how the advective (non-linear)
!		terms are computed. The value of 0 (default) uses
!		the usual finite element discretization over a single
!		element. The value of 1 chooses a semi-lagrangian
!		approach that is theoretically stable also for
!		Courant numbers higher than 1. It is however recommended
!		that the time step is limited using |itsplt| and
!		|coumax| described below. (Default 0)
! |iclin|	Linearization of the continuity equation. If |iclin|
!		is different from 0 the depth term in the continuity
!		equation is taken to be constant. (Default 0)

	call addpar('ilin',1.)
	call addpar('itlin',0.)
	call addpar('iclin',0.)

!c undocumented: rlin can lower effect of non linear terms

	call addpar('rlin',1.)

! The next parameters allow for a variable time step in the
! hydrodynamic computations. This is especially important for the
! non-linear model (|ilin=0|) because in this case the criterion
! for stability cannot be determined a priori and in any case the
! time integration will not be unconditionally stable.
!
! The variable time steps allows for longer basic time steps
! (here called macro time steps) which have to be set in |idt|.
! It then computes the optimal time step (here micro time step)
! in order to not exceed the given Courant number. However,
! the value for the macro time step will never be exceeded.
!
! Normally time steps are always given in full seconds. This is still
! true when specifying the macro time step |idt|. In older versions also
! the computed micro time steps also had to be full integer values.
! Starting from version 7.1 also fractional time steps are allowed.
! This gives the possibility to have time steps smaller than 1 s.
!
! |itsplt|	Type of variable time step computation. If this value
!		is 0, the time step will be kept constant at its initial
!		value. A value of 1 divides the initial time step into
!		(possibly) equal parts, but makes sure that at the end
!		of the micro time steps one complete macro time
!		step has been executed. The mode |itsplt| = 2
!		does not care about the macro time step, but always
!		uses the biggest time step possible. In this case
!		it is not assured that after some micro time steps
!		a macro time step will be recovered. Please note
!		that the initial macro time step will never be exceeded.
!		In any case, the time step will always be rounded to the
!		next lower integer value. This is not the case with
!		|itsplt| = 3 where the highest possible fractional time step 
!		will be used. (Default 0)
! |coumax|	Normally the time step is computed in order to not
!		exceed the Courant number of 1. However, in some cases
!		the non-linear terms are stable even for a value higher
!		than 1 or there is a need to achieve a lower Courant number.
!		Setting |coumax| to the desired Courant number
!		achieves exactly this effect. (Default 1)
! |idtsyn|	In case of |itsplt| = 2 this parameter makes sure that
!		after a time of |idtsyn| the time step will be synchronized
!		to this time. Therefore, setting |idtsyn| = 3600 means
!		that there will be a time stamp every hour, even if the model
!		has to take one very small time step in order to reach that
!		time. This parameter is useful
!		only for |itsplt| = 2 and its default value of
!		0 does not make any synchronization.
! |idtmin|	This variable defines the smallest time step possible
!		when time step splitting is enabled. Normally the smallest
!		time step is 1 second. Please set |idtmin| to values
!		smaller than 1 in order to allow for fractional time steps.
!		A value of 0.001 allows for time steps of down to
!		1 millisecond. (Default 1)

	call addpar('itsplt',0.)
	call addpar('coumax',1.)
	call addpar('idtsyn',0.)
	call addpar('idtmin',1.)
	call addpar('tfact',0.)		!still to comment FIXME

! These parameters define the weighting of time level in the 
! semi-implicit algorithm. With these parameters the damping
! of gravity or Rossby waves can be controlled. Only modify them if
! you know what you are doing.
!
! |azpar|	Weighting of the new time level of the transport
!		terms in the continuity equation. (Default 0.5)
! |ampar|	Weighting of the new time level of the pressure
!		term in the momentum equations. (Default 0.5)
! |afpar|	Weighting of the new time level of the Coriolis
!		term in the momentum equations. (Default 0.5)
! |avpar|	Weighting of the new time level of the non-linear
!		advective terms in the momentum equations. (Default 0.0)

	call addpar('azpar',0.5)
	call addpar('ampar',0.5)
	call addpar('afpar',0.5)
	call addpar('avpar',0.0)

! The next parameters define the weighting of time level for the
! vertical stress and advection terms. They guarantee the stability
! of the vertical system. For this reason they are normally set to
! 1 which corresponds to a fully implicit discretization. Only
! modify them if you know what you are doing.
!
! |atpar|	Weighting of the new time level of the vertical
!		viscosity in the momentum equation. (Default 1.0)
! |adpar|	Weighting of the new time level of the vertical
!		diffusion in the scalar equations. (Default 1.0)
! |aapar|	Weighting of the new time level of the vertical
!		advection in the scalar equations. (Default 1.0)

	call addpar('atpar',1.0)	!time weighting for vertical viscosity
	call addpar('adpar',1.0)	!time weighting for vertical diffusion
	call addpar('aapar',1.0)	!time weighting for vertical advection

!c------------------------------------------------------------------------

! DOCS	CORIOLIS	Coriolis parameters
!
! The next parameters define the parameters to be used
! with the Coriolis terms.
!
! |icor|	If this parameter is 0, the Coriolis terms are
!		not included in the computation. A value of 1
!		uses a beta-plane approximation with a variable
!		Coriolis parameter $f$, whereas a value of
!		2 uses an f-plane approximation where the
!		Coriolis parameter $f$ is kept constant over the
!		whole domain. (Default 0)
! |dlat|	Average latitude of the basin. This is used to
!		compute the Coriolis parameter $f$. This parameter 
!		is not used if spherical coordinates are used 
!		(|isphe|=1) or if a coordinate 	projection is set 
!		(|iproj| $>$0). (Default 0)
! |isphe|	If 0 a cartesian coordinate system is used,
!		if 1 the coordinates are in the spherical system (lat/lon).
!		Please note that in case of spherical coordinates the
!		Coriolis term is always included in the computation, even
!		with |icor| = 0. If you really do not want to use the
!		Coriolis term, then please set |icor| = -1. The default is
!		-1, which means that the type of coordinate system will 
!		be determined automatically.

	call addpar('icor',0.)
	call addpar('dlat',100.)
	call addpar('isphe',-1.)

!c------------------------------------------------------------------------

! DOCS	DEPTH		Depth parameters
!
! The next parameters deal with handling depth values of the basin.
!
! |href|	Reference depth. If the depth values of the basin and
!		the water levels are referred to mean sea level, |href|
!		should be 0 (default value). Else this value is
!		subtracted from the given depth values. For example,
!		if |href = 0.20| then a depth value in the basin
!		of 1 meter will be reduced to 80 centimeters.

	call addpar('href',0.)

! |hzmin|	Minimum total water depth that will remain in a
!		node if the element becomes dry. (Default 0.01 m)
! |hzoff|	Total water depth at which an element will be
!		taken out of the computation because it becomes dry.
!		(Default 0.05 m)
! |hzon|	Total water depth at which a dry element will be
!		re-inserted into the computation.
!		(Default 0.10 m)

	call addpar('hzmin',0.01)
	call addpar('hzoff',0.05)
	call addpar('hzon',0.10)

! |hmin|	Minimum water depth (most shallow) for the whole
!		basin. All depth values of the basin will be adjusted
!		so that no water depth is shallower than |hmin|.
!		(Default is no adjustment)
! |hmax|	Maximum water depth (deepest) for the whole
!		basin. All depth values of the basin will be adjusted
!		so that no water depth is deeper than |hmax|.
!		(Default is no adjustment)

	call addpar('hmin',-99999.)
	call addpar('hmax',+99999.)

!c------------------------------------------------------------------------
!c friction - description in file subn35.f
!c------------------------------------------------------------------------

! \input{P_friction.tex}

	call addpar('ireib',0.)
	call addpar('czdef',0.)
	call addpar('iczv',1.)
	call addpar('uvmin',0.2) !Minimum current speed for ireib=10

!c------------------------------------------------------------------------

! DOCS	PHYSICAL		Physical parameters
!
! The next parameters describe physical values that can be adjusted
! if needed.
!
! |rowass|	Average density of sea water. (Default 1025 \densityunit)
! |roluft|	Average density of air. (Default 1.225 \densityunit)
! |grav|	Gravitational acceleration. (Default 9.81 \accelunit)

        call addpar('rowass',1025.)
        call addpar('roluft',1.225)
        call addpar('grav',9.81)

!c------------------------------------------------------------------------
!c wind - description in file submeteo2.f
!c------------------------------------------------------------------------

! \input{P_wind.tex}

        call addpar('iwtype',1.)
	call addpar('itdrag',0.)
	call addpar('dragco',2.5e-3)
	call addpar('wsmax',50.)
	call addpar('wslim',-1.)
	call addpar('rfact',1.)

	call addpar('ibwind',0.)	!custom - not yet tested

!c------------------------------------------------------------------------

! DOCS  meteo                      Meteo and heat flux parameters

! The next parameters deal with the heat and meteo forcing.

! |iheat|	The type of heat flux algorithm (Default 1):
!		\begin{description}
!		\item[1] As in the AREG model
!		\item[2] As in the POM model
!		\item[3] Following A. Gill
!		\item[4] Following Dejak
!		\item[5] As in the GOTM model
!		\item[6] Using the COARE3.0 module
!		\item[7] Read heat fluxes directly from file. The columns
!		in the data file must be "time srad qsens qlat qlong".
!		\item[8] MFS heat fluxes as Pettenuzzo et al., 2010
!		\end{description}
!		Except when |iheat| is 7, the time series file has
!		the columns "time srad airt rhum cc".

	call addpar('iheat',1.)		!type of heat flux routine

! |ihtype|	Different ways of how to specify water vapor content
!		are possible. Normally relative humidity has to be
!		given (|ihtype|=1). However, also wet bulb temperature
!		(|ihtype|=2), dew point temperature (|ihtype|=3), or
!		specific humidity (|ihtype|=4) can
!		be given. (Default 1).

	call addpar('ihtype',1.)	!type of water vapor

! |isolp|	The type of solar penetration parameterization
!		by one or more exponential decay curves.
!		|isolp| = 0 sets an e-folding decay of radiation 
!		(one exponential decay curve) as function of depth |hdecay|.
!		|isolp| = 1 sets a profile of solar radiation with two length
!		scale of penetration. Following the Jerlov (Jerlov, N. G., 1968
!		Optical Oceanography, Elsevier, 194pp) classification the type
!		of water is clear water (type I). 
!		(Default 0) 

        call addpar('isolp',0.)         !type of solar penetration   

! |iwtyp|	The water types from clear water (type I) to the most 
!		turbid water (coastal water 9) following the classification of
!		Jerlov (Jerlov, N. G., 1968 Optical Oceanography, 
!		Elsevier, 194pp). The possible values for |iwtyp| are:
!		\begin{description}
!		\item[0] clear water type I 
!		\item[1] type IA
!		\item[2] type IB
!		\item[3] type II
!		\item[4] type III
!		\item[5] type 1
!		\item[6] type 3
!		\item[7] type 5
!		\item[8] type 7
!		\item[9] type 9
!		\end{description}

        call addpar('iwtyp',0.)         !water type (works if |isolp|=1)

! |hdecay|	Depth of e-folding decay of radiation [m]. If |hdecay| = 0 
!		everything is absorbed in first layer (Default 0).

	call addpar('hdecay',0.)	!depth of e-folding decay of radiation

! |botabs|	Heat absorption at bottom [fraction] (Default 0).
!		\begin{description}
!		\item[=0] everything is absorbed in last layer
!		\item[=1] bottom absorbs remaining radiation
!		\end{description}

	call addpar('botabs',0.)	!heat absorption at bottom

! |albedo|	General albedo (Default 0.06).

	call addpar('albedo',0.06)	!general albedo

! |albed4|	Albedo for temp below 4 degrees (Default 0.06).

	call addpar('albed4',0.06)	!albedo for temp below 4 degrees

! |ievap| 	Compute evaporation mass flux (Default 0).
!		This option computes internally the evaporation rate
!		and adds it to the mass balance.

	call addpar('ievap',0.)		!compute evaporation mass flux

! |irain|	Compute precipitation mass flux (Default 1).
!		If no precipitation file is given or if |irain|=0, 
!		the precipitation rate
!		is set to zero. If a file is found and irain is set to
!		the default value, then the precipitation rate is considered
!		into the mass balance. Note that, if some evaporation data is
!		known this can be added into the precipation file as a net mass
!		balance. In this case please set |ievap|=0 to not add
!		evaporation twice.

        call addpar('irain',1.)         !compute precipitation mass flux

!c------------------------------------------------------------------------

! DOCS	3D			Parameters for 3d
!
! The next parameters deal with the layer structure in 3D.

! |dzreg|	Normally the bottom of the various layers are given in
!		section |\$levels|. If only a regular vertical grid is desired
!		then the parameter |dzreg| can be used. It specifies the spacing
!		of the vertical layers in meters. (Default is 0, which means 
!		that the layers are specified explicitly in |\$levels|.

	call addpar('dzreg',0.)		!regular vertical grid

! The last layer (bottom layer) is treated in a special way. Depending on
! the parameter |ilytyp| there are various cases to be considered. A value
! of 0 leaves the last layer as it is, even if the thickness is very small.
! A value of 1 will always eliminate the last layer, if it has not full
! layer thickness. A value of 2 will do the same, but only if the last layer
! is smaller than |hlvmin| (in units of fraction). Finally, a value of
! 3 will add the last layer to the layer above, if its layer thickness
! is smaller than |hlvmin|.
!
! |ilytyp|	Treatment of last (bottom) layer. 0 means no adjustment,
!		1 deletes the last layer, if it is not a full layer,
!		2 only deletes it
!		if the layer thickness is less than |hlvmin|, and 3
!		adds the layer thickness to the layer above if it is smaller
!		than |hlvmin|. Therefore, 1 and 2 might change the
!		total depth and layer structure, while 3 only might
!		change the layer structure. The value of 1 will always
!		give you full layers at the bottom. (Default 3)
! |hlvmin|	Minimum layer thickness for last (bottom) layer used when
!		|ilytyp| is 2 or 3. The unit is fractions of the nominal
!		layer thickness. Therefore, a value of 0.5 indicates that
!		the last layer should be at least half of the full
!		layer. (Default 0.25)

	call addpar('ilytyp',3.00)	!type of depth adjustment
	call addpar('hlvmin',0.25)	!min percentage of last layer thickness

! With $z-$layers the treatment of the free-surface must be addressed.
! What happens if the water level falls below the first $z$-level? 
! A $z-$star type vertical grid deformation can be deployed.
! The next parameter specifies the number of surface layers that are moving.

! |nzadapt|	Parameter that controls the number of surface $z-$layers
!		that are moving. The value $|nzadapt|\le 1$ corresponds
!		to standard $z-$layers (Default). Then, some care is
!		needed to define the first interface sufficiently deep to
!		avoid the well-known "drying" of the first layer. The
!		value of $|nzadapt| = N_{tot}$, with $N_{tot}$ the total 
!		numberof $z$-layers, is the conventional $z-$star 
!		(all layers are moving). Other values of 
!		$1<|nzadapt|<N_{tot}$ corresponds to move, at minimum,
!		the first |nzadapt| surface layers with $z-$star.
!		These feature is still experimental.

        call addpar('nzadapt',0.)      ! z-layers parameter

! The above parameters are dealing with zeta layers, where every layer
! has constant thickness, except the surface layer which is varying with
! the water level. The next parameters deal with sigma layers where all
! layers have varying thickness with the water level.

! |nsigma|	Number of sigma layers for the run. This parameter can
!		be given in alternative to specifying the sigma layers
!		in |\$levels|. Only regularly spaced sigma levels
!		will be created. (Default 0)
! |hsigma|	This is still an experimental feature. It allows to use
!		sigma layers above zeta layers. |hsigma| is the depth where
!		the transition between these two types of layers
!		is occurring. (Default 10000)

	call addpar('nsigma',0.)	!number of sigma layers
	call addpar('hsigma',10000.)	!lower depth of sigma layers (hybrid)

! The next parameters deal with vertical diffusivity and viscosity.

! |difmol|	Vertical molecular diffusivity parameter for
!		temperature, salinity, and tracer.
!		(Default 1.0e-06)
! |diftur|	Vertical turbulent diffusivity parameter for 
!		temperature, salinity, and tracer.
!		(Default 0)
! |vismol|	Vertical molecular viscosity parameter for momentum.
!		(Default 1.0e-06)
! |vistur|	Vertical turbulent viscosity parameter for momentum.
!		(Default 0)

        call addpar('difmol',1.0e-06)	!molecular vertical diffusivity
	call addpar('diftur',0.)	!diffusion parameter (vertical), cvpar
        call addpar('vismol',1.0e-06)	!molecular vertical viscosity
        call addpar('vistur',0.)	!turbulent vertical viscosity (nau)

! Instead of setting fixed values for viscosity and diffusivity, these
! parameters can also be computed by a turbulence closure scheme. 
! The parameter |iturb| defines what turbulence scheme is going to be used.
! A value of 0 uses no turbulence scheme. In this case be sure that 
! |vistur| and |diftur| have been set manually. |iturb=1| uses the GOTM
! routines for turbulence. In order to use this value |SHYFEM| has to be
! compiled with GOTM support. |iturb=2| uses a k-epsilon module, and a value
! of 3 uses the Munk-Anderson model. The recommended value for |iturb| 
! is 1 (GOTM module). (Default 0)

	call addpar('iturb',0.)		!use turbulence closure scheme

!c------------------------------------------------------------------------

! The next parameters deal with horizontal diffusion.
!c horizontal diffusion (Smagorinsky)
!c typical values for ahpar with Smagorinski: 0.2 - 0.4
!c idhtyp: 0=constant  1=weigthed with area  2=Smagorinsky

	call addpar('idhtyp',0.)	!type of horizontal diffusion/viscosity
	call addpar('dhlen',1.) 	!length scale for idhtyp=1
	call addpar('noslip',0.) 	!no slip conditions on boundary

	call addpar('idtype',2.) 	!type of hor diffus (delete after test)

	call addpar('ahpar',0.)		!austausch coefficient (still same??)

! |dhpar|	Horizontal diffusion parameter (general).
!		(Default 0)

	call addpar('dhpar',0.)		!diffusion parameter

! The next parameters deal with the control of the scalar transport 
! and diffusion equation. You have possibility to prescribe the tvd scheme
! desired and to limit the Courant number.
!
! |itvd|	Type of the horizontal advection scheme used for 
!		the transport and diffusion
!		equation. Normally an upwind scheme is used (0), but setting
!		the parameter |itvd| to a value greater than 0 
!		choses a TVD scheme. A value of 1 will use a TVD scheme
!		based on the average gradient, and a value of 2 will use
!		the gradient of the upwind node (recommended).
!		This feature is still experimental, so use with care. 
!		(Default 0)
! |itvdv|	Type of the vertical advection scheme used for 
!		the transport and diffusion
!		equation. Normally an upwind scheme is used (0), but setting
!		the parameter |itvdv| to 1 choses a TVD scheme. This feature
!		is still experimental, so use with care. (Default 0)
! |rstol|	Normally the internal time step for scalar advection is
!		automatically adjusted to produce a Courant number of 1
!		(marginal stability). You can set |rstol| to a smaller value 
!		if you think there are stability problems. (Default 1)

	call addpar('itvd',0.)		!horizontal TVD scheme?
	call addpar('itvdv',0.)		!vertical TVD scheme?
	call addpar('rstol',1.)		!limit time step to this Courant number

!c------------------------------------------------------------------------

! DOCS	VAR		Various parameters
!
! The next parameters describe various parameters not related to
! the above parameters.

! |tauvel|	If you have velocity observations given in file
!		|surfvel| then you can specify the relaxation
!		parameter $\tau$ in the variable |tauvel|. (Default 0,
!		which means no assimilation of velocities)

	call addpar('tauvel',0.)	!horizontal TVD scheme?

! |rtide|	If |rtide| = 1 the model calculates equilibrium tidal 
!		potential and load tides and uses these to force the 
!		free surface (Default 0).

	call addpar('rtide',0.)

! |ltidec|	Calibration factor for calculating the loading tide, 
!		which is computed in function of the total water depth as
!		$\beta=ltidec*H$. Usually it has a value of order 1e-6. 
!		If 0 no loading tide is computed (Default 0).

	call addpar('ltidec',0.)

! |ibstrs|	Call parameter for the routine bstress. If equal to 1 it 
!		computes (in function of currents and waves) and writes 
!		the bottom shear stess into a .bstress.shy file.
!		If 0 no bottom stess is computed (Default 0).

	call addpar('ibstrs',0.)

!c------------------------------------------------------------------------

! DOCS	ST	Temperature and salinity
!
! The next parameters deal with the transport and diffusion
! of temperature and salinity (T/S). 

! |ibarcl|	In order to compute
!		T/S the parameter |ibarcl| must be different from 0. 
!		Different values indicate different ways to compute 
!		the advection and diffusion of the variables T/S 
!		and how their values is used in
!		the momentum equation. (Default 0)
!		\begin{description}
!		\item[0] No computation of T/S.
!		\item[1] Computation of T/S. The T/S field has a feedback
!			 through the baroclinic term on the momentum
!			 equation. This corresponds to a full
!			 baroclinic model.
!		\item[2] The model runs in diagnostic mode. No advection
!			 of the T/S field is done, but the baroclinic term
!			 is computed and used in the momentum equations.
!			 For this to work T/S fields have to be provided
!			 in external files |tempobs| and |saltobs|.
!		\item[3] The model computes the advection and diffusion
!			 of T/S, but their value is not used in the baroclinic
!			 terms in momentum equation. This corresponds to
!			 treating T/S as tracers.
!		\item[4] The model runs in baroclinic mode, but
!			 uses nudging for a basic data assimilation of T/S.
!			 For this to work T/S fields have to be provided
!			 in external files |tempobs| and |saltobs|.
!			 Moreover, a time scale for the nudging has to be
!			 provided, either as a constant using
!			 |temptaup| and |salttaip| or in spatially and
!			 temporarily varying fields given in external
!			 files |temptau| and |salttau|.
!		\item[5] As for case 4, but no baroclinic mode. It therefore
!			 is as in case 3, but with nudging.
!		\end{description}

	call addpar('ibarcl',0.)	!compute baroclinic contributions ?

! In case |ibarcl| is different from 0 the variables T/S will
! be computed. However, they may be selectively turned off setting
! one of the two parameters |itemp| or |isalt| explicitly to 0.
!
! |itemp|	Flag if the computation on the temperature is done.
!		A value different from 0 computes the transport
!		and diffusion of the temperature. (Default 1)
! |isalt|	Flag if the computation on the salinity is done.
!		A value different from 0 computes the transport
!		and diffusion of the salinity. (Default 1)

	call addpar('itemp',1.)		!compute temperature ?
	call addpar('isalt',1.)		!compute salinity ?

! The next parameters set the initial conditions for temperature and salinity.
! Both the average value and and a stratification can be specified.
! Initial conditions can also be given in external files |tempin| and |saltin|.
!
! |temref|	Reference (initial) temperature of the water in
!		centigrade. (Default 0)
! |salref|	Reference (initial) salinity of the water in
!		psu (practical salinity units) or ppt.
!		(Default 0)
! |tstrat|	Initial temperature stratification in units of [C/km].
!		A positive value indicates a stable stratification.
!		(Default 0)
! |sstrat|	Initial salinity stratification in units of [psu/km].
!		A positive value indicates a stable stratification.
!		(Default 0)

	call addpar('temref',0.)	!reference temperature for T/S runs
	call addpar('salref',0.)	!reference salinity for T/S runs

	call addpar('sstrat',0.)	!salt stratification
	call addpar('tstrat',0.)	!temp stratification

! The next parameters deal with horizontal diffusion of temperature
! and salinity. These parameters overwrite the general parameter for
! horizontal diffusion |dhpar|.

! |thpar|	Horizontal diffusion parameter for temperature.
!		(Default 0)
! |shpar|	Horizontal diffusion parameter for salinity.
!		(Default 0)

	call addpar('thpar',-1.)	!horiz. diff. coeff. for temp.
	call addpar('shpar',-1.)	!horiz. diff. coeff. for sal.

! When using nudging (|ibarcl=4|) time scales for nudging have to given.
! They can be given in the parameters |temptaup| and |salttaup| or
! in external files |temptau| and |salttau|.

! |temptaup|	Time scale for temperature nudging in seconds.
!		(Default 0)
! |salttaup|	Time scale for salinity nudging in seconds.
!		(Default 0)

	call addpar('temptaup',-1.)	!tau for temperature
	call addpar('salttaup',-1.)	!tau for salinity

!c------------------------------------------------------------------------

! DOCS	CC	Concentrations
!
! The next parameters deal with the transport and diffusion
! of a conservative substance. The substance is dissolved in
! the water and acts like a tracer.
!
! |iconz|	Flag if the computation on the tracer is done.
!		A value different from 0 computes the transport
!		and diffusion of the substance. If greater than 1
!		|iconz| concentrations are simulated. (Default 0)
!c |iconza|	The concentration are normally written as instantaneous
!c		values. If |iconza| is different from 0, instead of
!c		the instantaneous values the minimum, average and maximm
!c		values are written. (Default 0)
! |conref|	Reference (initial) concentration of the tracer in
!		any unit. (Default 0)
! |taupar|	Decay rate for concentration if different from 0. In
!		this case |taupar| is the decay rate (e-folding time) in days.
!		This parameter is also used for multi-concentration runs.
!		In this case either one value has to be given that is used
!		for all concentrations, or |iconz| values have to be given,
!		one for each concentration.
!		(Default 0)
! |wsettl|	Settling rate for concentration if different from 0. In
!		this case |wsettl| is the sinking velocity in [m/s].
!		This parameter is also used for multi-concentration runs.
!		In this case either one value has to be given that is used
!		for all concentrations, or |iconz| values have to be given,
!		one for each concentration.
!		(Default 0)
! |idecay|	Type of decay used. If 0 no decay is used.
!		A value of 1 uses the value of |taupar| as exponential decay.
!		A value of 2 uses a formulation of Chapra, where the
!		decay rate depends on T,S,light and settling. In this
!		case the value of |taupar| is ignored.
!		(Default 0)
! |iage|	Age concentration has been enabled. In the variable conz
!		the value is referred to the age of the water body,
!		not the concentration. Boundary values for conz should
!		be set to 0. (Default 0)

	call addpar('iconz',0.)		!compute concentration ?
!c	call addpar('iconza',0.)	!average values ? - not working
	call addpar('conref',0.)	!reference concentration
	call addpar('idecay',0.)	!type of decay
	call addpar('iage',0.)		!age concentration

	call para_deprecate('contau','taupar')	!no contau anymore -> use taupar
	call para_add_array_value('taupar',0.)	!decay rate [days]
	call para_add_array_value('wsettl',0.)	!sinking velocities [m/s]

! |chpar|	Horizontal diffusion parameter for the tracer.
!		This value overwrites the general parameter for
!		horizontal diffusion |dhpar|. (Default 0)

	call addpar('chpar',-1.)	!diffusion parameter

!c------------------------------------------------------------------------

! DOCS	CC	Concentrations
!
! The next parameters deal with limiting the scalars of temperature,
! salinity, and concentration to reasonable values. They should only
! be used if some strange values are computed (due to evaporation in
! salt marshes etc..). Using these routines will not conserve the
! total quantity of these scalars.

! Limiting values can be setup separately for the scalars and for both
! limiting the minimum and maximum values. If both min and max have been
! setup, the max limiter must be>= min. If no values are given for the
! parameters below, no limiter will be implemented.

! |tlimit0|, |tlimit1|	The min and max limiting values for temperature.
! |slimit0|, |slimit1|	The min and max limiting values for salinity.
! |climit0|, |climit1|	The min and max limiting values for concentration.

	call addpar('tlimit0',-999.)
	call addpar('tlimit1',-999.)
	call addpar('slimit0',-999.)
	call addpar('slimit1',-999.)
	call addpar('climit0',-999.)
	call addpar('climit1',-999.)

!c------------------------------------------------------------------------

! DOCS	STOUTPUT	Output for scalars
!
! The next parameters define the output frequency of the
! computed scalars (temperature, salinity, generic concentration) to file.
!
! |idtcon|, |itmcon|	Time step and start time for writing to file
!			|.conz.shy| (concentration) and
!			|.ts.shy| (temperature and salinity).
! |irho|		Flag to indicate if the density is also written 
!			together with T/S. A value different from 0 
!			writes the density to file. (Default 0)
! |iskin|		Flag to indicate if the skin temperature is written
!			to file |.tskin.shy|.
!			A value different from 0 writes the 
!			skin temperature to file. (Default 0)

	call addpar('idtcon',0.)	!time step for output
	call addpar('itmcon',-1.)	!minimum time for output
	call addpar('irho',0.)		!write rho
	call addpar('iskin',0.)		!write skin temperature

! DOCS	END

!c------------------------------------------------------------------------
!c------------------------------------------------------------------------
!c------------------------------------------------------------------------
!c------------------------------------------------------------------------
!c------------------------------------------------------------------------

!c------------------------------------------------------------------------
!c still to be commented below here
!c------------------------------------------------------------------------

!c------------------------------------------------------------------------

	call addpar('ihydro',1.)	!compute hydrodynamics (for debug)

!c------------------------------------------------------------------------

	call addpar('ibio',0.)		!run biological reactor
        call addpar('imerc',0.)		!run mercury reactor
        call addpar('issedi',0.)	!run simple sediments
	call addpar('itoxi',0.)		!run toxicological routines

!c------------------------------------------------------------------------

!c call routines in flux.f

	call addpar('regflx',0.)	!compute regular fluxes

!c------------------------------------------------------------------------

!c custom call

	call addpar('icust',0.)		!call custom routine
	call addpar('tcust',0.)		!time for custom routine

	call addpar('ishyff',1.)	!SHYFEM file format
					!0=old 1=new 2=both

	call addpar('ipoiss',0.)	!solve poisson equation

!c rain

	call addpar('zdist',0.)		!distributed water level

	call addpar('idtbox',0.)	!for boxes
	call addpar('itmbox',-1.)

	call addpar('ihwadv',0.)        !vert advection of horiz momentum

!c ice

	call addpar('iaicef',-99.)	!area code for ice free condition
	call addpar('icemod',0.)	!use ice model?
	call addpar('idtice',0.)	!time step for ice model

       ! whether to output - instantaneous salt temp rho fields (icon_av=0),  
       !                   - time-averaged fields (icon_av=1),  
       !                   - vertically averged and time-averaged fields (icon_av=1)
	call addpar('icon_av',0.) 
       ! whether to output - instantaneous hydrodynamic fields (iout_av=0),  
       !                   - time-averaged fields (iout_av=1),  
	call addpar('iout_av',0.) 
       ! whether to output - instantaneous ecological fields (ieco_av=0),  
       !                   - time-averaged fields (ieco_av=1),  
       !                   - vertically averaged and time-averaged fields (ieco_av=2)
	call addpar('ieco_av',0.) 
       ! whether to output - instantaneous layer thickness fields (ilayers=0),  
       !                   - time-averaged thickness fields (ilayers=1),  
       !                   - vertically averaged and time-averaged fields (ilayers=2)
	call addpar('ilayers',-1.) 

	end

!************************************************************************

	subroutine nlsinh_lagrg

	implicit none

! $lagrg section

! DOCS	START	P_lagrg
!
! Section |$lagrg| describes the use of the Lagrangian Particle Module.
!
! The lagrangian particles can be released:
! \begin{itemize}
! \item inside the given areas (filename |lgrlin|). If this file is not 
!      specified they are released over the whole domain. The amount of
!      particles released and the time step are specified by |nbdy| and 
!      |idtl|.
! \item at selected times and location, e.g. along a drifter track
!      (filename |lgrtrj|). |nbdy| particles are released at the times
!      and location specified in the file.
! \item as initial particle distribution (filename |lgrini|) at time
!      |itlgin|. This file has the same format as the lagrangian output.
! \item at the open boundaries either as particles per second or per
!      volume flux (parameter |lgrpps|).
! \end{itemize}
!
! |lgrlin| and |lgrtrj| are mutually exclusive. 
!
! The lagrangian module runs between the times |itlanf| and |itlend|. If one
! or both are missing, the simulation extremes are substituted. Inside
! the lagrangian simulation window, the release of particles inside a given
! area is controlled by the parameters |idtl|, |itranf| and |itrend|. 
! |itranf| gives the time of the first release, |itrend| the time for 
! the last release. If not given they are set equal to the extremes of 
! the lagrangian simulation. |idtl| is giving the time step of release.
!
! The output frequency of the results can be controlled by 
! |idtlgr| and |itmlgr|.
!
! Please find all details here below.

        call sctpar('lagrg')             !sets default section
        call sctfnm('lagrg')

! |ilagr|	Switch that indicates that the lagrangian module
!		should be run (default 0):
!		\begin{description}
!		\item[0] do nothing
!		\item[1] surface lagrangian
!		\item[2] 2d lagrangian
!		\item[3] 3d lagrangian
!		\end{description}

	call addpar('ilagr',0.)         !LAGR

! |nbdymax|	Maximum numbers of particles that can be in the domain.
!		This should be the maximum number of particles
!		that can be created and inserted. Use 0 to not limit
!		the number of particles (on your own risk). This
!		parameter must be set and has no default.

        call addpar('nbdymax',-1.)

! |nbdy|	Total numbers of particles to be released in the domain each
!		time a release of particles takes place. 
!		(Default 0)

        call addpar('nbdy',0.)

! |rwhpar|	A horizontal diffusion can be defined for the lagrangian model.
!		Its value can be specified in |rwhpar| and the units are 
!		[m**2/s]. If |rwhpar<0| the diffusion parameter depends on
!		the local diffusivity (see |idhtyp|) (Default 0)

	call addpar('rwhpar',0.)	!diffusion for lagrangian model

! |itlanf, itlend|	The start and end time for the lagrangian module.
!			If not given, the module runs for the whole simulation.

        call addpar('itlanf',-1.)
        call addpar('itlend',-1.)

! |itmlgr, idtlgr|	Initial time and time step for the output to file
!			of the particles. if |idtlgr| is 0, 
!			no output is written. (Default 0)

        call addpar('itmlgr',-1.)
        call addpar('idtlgr',0.)

! |idtl|	The time step used for the release of particles. If this
!		is 0 particles are released only once at the beginning
!		of the lagrangian simulation. No particles are released 
!		for a value of less than 0. (Default 0)

        call addpar('idtl',0.)

! |itranf, itrend|	Initial and final time for the release of particles.
!			If not specified the particles are released over the
!			whole lagrangian simulation period.

        call addpar('itranf',-1.)
        call addpar('itrend',-1.)

! |ipvert|	Set the vertical distribution of particles:
!		\begin{description}
!		\item[0] releases particles only in surface layer
!		\item[$>$ 0] release n particles regularly
!		\item[$<$ 0] release n particles randomly
!		\end{description}

        call addpar('ipvert',0.)

! |linbot|	Set the bottom layer for vertical releases (Default -1, bottom layer)

        call addpar('linbot',-1.)

! |lintop|	Set the top layer for vertical releases (Default 1, surface layer)
        call addpar('lintop',1.)

! |stkpar|	Calibration parameter for parameterizing the stokes drift 
!		induced by waves (and wind). Only affect particle of the sea
!		surface (layer = 1). The wind file is needed even in offline 
!		mode (Default 0). 

        call addpar('stkpar',0.)

! |dripar|	Parameter to account for drifter inertia by multiplying
!		the advective transports. Usually it assumes values between 
!		0.9 and 1.2 (Default 1). 

        call addpar('dripar',1.)

! |lbeach|	Parameter to account for particles beaching on the shore. It 
!		assumes values between 0 (no beaching) and 1 (Default 0).

        call addpar('lbeach',0.)

! |lgrlin|	File name that contains closed lines of the area where
!		the particles have to be released. If not given, the particles
!		are released over the whole domain.

        call addfnm('lgrlin',' ')

! |lgrtrj|	File name that contains a drifter trajectory with time 
!		(yyyy-mm-dd::HH:MM:SS) and position (x and y) of release 
!		of |nbdy| particles.

        call addfnm('lgrtrj',' ')

! |lgrini|	File name that contains initial particle distribution.
!		It has the same format as the lagrangian output.

        call addfnm('lgrini',' ')

! |itlgin|	Time to use for the initialization of the particle 
!		distribution from file (|lgrini|). 

        call addpar('itlgin',-1.)

!c To compute the transit time of the particles to leave
!c a specified area. |artype| is the flag to detect the
!c element outside the defined area. 
!c Default = -1, i.e., no transit times are computed

        call addpar('artype',-1.)

!c still to be commented

        call addpar('lcust',0.)
        call addpar('ldecay',0.)
        call addpar('ioil',0.)
        call addpar('ilarv',0.)
        call addpar('ised',0.)


! DOCS	END

	end

!************************************************************************

        subroutine nlsinh_wrt

        implicit none

! $wrt section

! DOCS  START   P_wrt
!
! Section |$wrt| contains parameters for computing water renewal time.
! During runtime if writes a .jas file with time series of total tracer
! concentration in the basin and WRT computed according to different methods.
! Nodal values of computed WRT are written in the .wrt file.
! Frequency distributions of WRTs are written in the .frq file.
!
! Please find all details here below.

        call sctpar('wrt')             !sets default section
        call sctfnm('wrt')

! |idtwrt|	Time step to reset concentration to c0. Use 0 if no reset
!		is desired. Use -1 if no renewal time computation is desired
!		(Default -1).

        call addpar('idtwrt',-1.)

! |itmin|	Time from when to compute renewal time (-1 for start of sim)
!		(Default -1)

        call addpar('itmin',-1.)

! |itmax|	Time up to when to compute renewal time (-1 for end of sim)
!		(Default -1).

        call addpar('itmax',-1.)

! |c0|		Initial concentration of tracer (Default 1).

        call addpar('c0',1.)

! |iaout|	Area code of elements out of lagoon (used for init and retflow).
!		Use -1 to if no outside areas exist. (Default -1).

        call addpar('iaout',-1.)

! |percmin|	Percentage to reach after which the computation is stopped.
!		Use 0 if no premature end is desired (Default 0).

        call addpar('percmin',0.)

! |iret|	Equal to 1 if return flow is used. If equal to 0 the 
!		concentrations outside are explicitly set to 0 (Default 1).

        call addpar('iret',1.)

! |istir|	If equal to 1 simulates completely stirred tank 
!		(replaces at every time step conz with average conz)
!		(Default 0).

        call addpar('istir',0.)

! |iadj|	Adjust renewal time for tail of distribution (Default 1).

        call addpar('iadj',1.)

! |ilog|	Use logarithmic regression to compute renewal time (Default 0).

        call addpar('ilog',0.)

! |ctop|	Maximum to be used for frequency curve (Default 0).

        call addpar('ctop',0.)

! |ccut|	Cut renewal time at this level (for res time computation)
!		(Default 0).

        call addpar('ccut',0.)

! |wrtrst|	If reset times are not regularly distributed (e.g., 1 month)
!		it is possible to give the exact times when a reset should
!		take place. |wrtrst| is a file name where these reset times
!		are specified, one for each line. For every line two integers
!		indicating date and time for the reset must be specified.
!		If only one value is given, time is taken to be 0. The format
!		of date is "YYYYMMDD" and for time "hhmmss". If the file
!		wrtrst is given |idtwrt| should be 0.

        call addfnm('wrtrst',' ')

! DOCS  END

        end

!************************************************************************

	subroutine nlsinh_bfmsc

! Parameters for BFM ECOLOGICAL module        !BFMSC

	implicit none

        call sctpar('bfmsc')             !sets default section
        call sctfnm('bfmsc')

! BFM ECOLOGICAL MODEL CALL 

	call addpar('ibfm',0.)
	call addpar('ibtanf',0.)
	call addpar('ibtend',0.)
	call addpar('itmbfm',-1.)
	call addpar('idtbfm',0.)
	call addpar('bligth',1.) !light flag=1 max/min light W/m**2 in nml file

	end

!************************************************************************

	subroutine nlsinh_proj

	implicit none

! $proj section

! DOCS  START   P_proj
!
! Section |$proj| handles the projection from cartesian to
! geographical coordinate system. If |iproj| $>$0 the projected geographical 
! coordinates can be used for computing spatially variable Coriolis parameter 
! and tidal potential even if the basin is in cartesian coordinate system 
! (|isphe| = 0) .
!
! Please find all details here below.

        call sctpar('proj')             !sets default section
        call sctfnm('proj')

! |iproj|	Switch that indicates the type of projection
!		(default 0):
!		\begin{description}
!		\item[0] do nothing
!		\item[1] Gauss-Boaga (GB)
!		\item[2] Universal Transverse Mercator (UTM)
!		\item[3] Equidistant cylindrical (CPP)
!		\item[4] UTM non standard
!		\end{description}

	call addpar('iproj',0.)

! |c\_fuse|	Fuse for GB (1 or 2, default 0)

	call addpar('c_fuse',0.)

! |c\_zone|	Zone for UTM (1-60, default 0)

	call addpar('c_zone',0.)

! |c\_lamb|	Central meridian for non-std UTM (default 0)

	call addpar('c_lamb',0.)

! |c\_x0|	x0 for GB and UTM (default 0)

	call addpar('c_x0',0.)

! |c\_y0|	y0 for GB and UTM (default 0)

	call addpar('c_y0',0.)

! |c\_skal|	Scale factor for non-std UTM (default 0.9996)

	call addpar('c_skal',0.9996)

! |c\_phi|	Central parallel for CPP (default 0.9996)

	call addpar('c_phi',0.)

! |c\_lon0|	Longitude origin for CPP (default 0)
	call addpar('c_lon0',0.)

! |c\_lat0|	Latitude origin for CPP (default 0)
	call addpar('c_lat0',0.)

! DOCS  END

	end

!************************************************************************

	subroutine nlsinh_undoc

! not documented parameters 

	implicit none

	call sctpar('para')		!sets default section

!c undocumented parameters

	call addpar('iclose',0.)
	call addpar('itsmed',0.)	!averages for T/S

!c internally used parameters

	call addpar('flag',0.)
	call addpar('zconst',-999.)	!set to flag
	call addpar('volmin',1.)	!minimum volume to remain in el.

!c debug

	call addpar('vreps',5.e-4)
	call addpar('vrerr',1.e-1)
	call addpar('levdbg',0.)	!debug level (0 = no, 9 = max)

!c distance for advective terms

	call addpar('nadist',0.)

!c new for scaling time step

	call addpar('itunit',1.)	!this is not supported anymore

!c experimental stuff

        !call addpar('nbsig',0.)         !sigma layers to read in for OBC

        call addpar('sedim',0.)         !sedimentation for theseus
        call addfnm('hsedim',' ')         !sedimentation hev file for theseus

        call addpar('nomp',0.)          !number of threads to use (-1 for all)

	end

!************************************************************************

	subroutine nlsinh_georg

! parameters used for private projects

	implicit none

	call sctpar('para')		!sets default section

	call addpar('hlido',999.)	!maximum depth at lido
	call addpar('hmala',999.)	!maximum depth at lido
	call addpar('hchio',999.)	!maximum depth at lido

	call addpar('zrise',0.)		!sea level rise
	call addpar('zsalv',999.)	!save-guarding level
	call addpar('zfranc',0.)	!extra security for forecast

	end

!************************************************************************

	subroutine nlsinh_unused

! parameters not used anymore -> to be deleted

	implicit none

	call sctpar('para')		!sets default section

	call addpar('iprogr',0.)

	end

! ********************************************************************

! This subroutine defines the simulation wave parameter

        subroutine nlsinh_waves

        implicit none

! $waves section

! DOCS  START   P_waves
!
! Parameters in section |$waves| activate the wind wave module and define
! which kind of wind wave model has to be used. These parameters must
! be in section |waves|.

        call sctpar('waves')             !sets waves section
        call sctfnm('waves')

! |iwave|	Type of wind wave model and coupling procedure (default 0):
!		\begin{description}
!		\item[0] No wind wave model called 
!		\item[1] The parametric wind wave model is called 
!		(see file subwave.f)
!		\item[$>$1] The spectral wind wave model WWMIII is called
!		\item[2] ... wind from SHYFEM, radiation stress formulation
!		\item[3] ... wind from SHYFEM, vortex force formulation 
!		\item[4] ... wind from WWMIII, radiation stress formulation
!		\item[5] ... wind from WWMIII, vortex force formulation 
!		\item[11] The spectral wind wave model WaveWatch WW3 is called
!		\end{description}
!		When the vortex force formulation is chosen the wave-supported
!		surface stress is subtracted from the wind stress, in order to
!		avoid double counting of the wind forcing in the flow model.
!		Moreover, the use of the wave-depended wind drag coefficient 
!		could be adopted setting |itdrag| = 3.

        call addpar('iwave',0.)

!
! |dtwave|	Time step for coupling with WWMIII. Needed only for
!		|iwave| $>$ 1 (default 0).

        call addpar('dtwave',0.)

! |idtwav|, |itmwav|	Time step and start time for writing to file wav,
!			the files containing output wave variables (significant 
!			wave height, wave period, mean wave direction).

        call addpar('idtwav',0.)
        call addpar('itmwav',-1.)

!
! DOCS  END
!

        end

!************************************************************************

! This subroutine defines the ocean-atmosphere coupling

        subroutine nlsinh_atm

! $atm section

! DOCS  START   P_atm
!
! Parameters in section |$atm| activate the ocean-atmosphere coupling.

        implicit none

        call sctpar('atm')           !sets default section
        call sctfnm('atm')

! |iatm|        Ocean-atmosphere coupling (default 0):
!               \begin{description}
!               \item[0] No coupling with the atmosphere.
!               \item[1] SHYFEM runs coupled with an atmospheric model within the ESMF framwork.
!                        For now the atmospheric model available is WRF.
!                        You should have a look the the specific section to see how to run
!                        a coupled SHYFEM-WRF simulation.
!               \end{description}

        call addpar('iatm',0.)        !0=SHYFEM standalone, 1=coupled with atmospheric model

! |idtatm|      Time step (only in seconds for now) for writing the atmospheric
!               fields computed by the atmosphere component and
!               remapped onto the SHYFEM grid, to a vtk file. This is useful
!               especially for debugging. For operational runs it is
!               better not to specify it (no atmospheric fields are printed).

        call addpar('idtatm',2592000.)!time step for output

!
! DOCS  END
!

        end

!************************************************************************

	subroutine nlsinh_nonhydro

! parameters for non hydrostatic model (experimental)

	implicit none

        call sctpar('nonhyd')           !sets default section

	call addpar('inohyd',0.)	! 0=hydrostatic  1=non-hydrostatic

        call addpar('aqpar',0.5)	!implicit parameter

	call addpar('islope',0.)	!type of grid for poisson equation

        call addpar('ivwadv',0.)        !vert advection of vert momentum
        call addpar('inhflx',0.)        !flux upwind for horiz advect of w
	call addpar('inhadj',0.)        !choice for correction of U,V,eta
        call addpar('inhwrt',0.)        !output every inhwrt time steps
        call addpar('inhbnd',0.)        !exclude NH dynamics for boundaries
        call addpar('iwvel',0.)         !write vertical velocity
        call addpar('iqpnv',0.)         !write NH pressure !DWNH
        call addpar('nqdist',0.)        !distance for NH pressure terms

	end

!************************************************************************

	subroutine nlsinh_connect

! parameters for connectivity (experimental)

	implicit none

        call sctpar('connec')          !sets default section
        call sctfnm('connec')

	call addpar('icnn',0.)	        !number of stations for connectivity

        call addpar('radcnn',0.)	!radius of release area
        call addpar('ppscnn',0.)	!particles per second to be released
        call addpar('pldcnn',0.)	!pelagic larval duration
        call addpar('idtcnn',0.)	!output every idtcnn

        call addfnm('statcnn',' ')      !connectivity file with stations

	end

!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************

	subroutine fnminh

! initializes default names

	use para

	implicit none

	logical haspar

! $name section

	call sctfnm('name')		!sets default section

	call para_deprecate('basdir','')
	call para_deprecate('datdir','')
	call para_deprecate('tmpdir','')
	call para_deprecate('defdir','')

! DOCS	START	S_name
!
! In section |$name| the names of input files can be
! given. All directories default to the current directory,
! whereas all file names are empty, i.e., no input files are
! given.
!
! DOCS	FILENAME	File names
!
! Strings in section |$name| enable the specification of files
! that account for initial conditions or forcing.
!
! |zinit|	Name of file containing initial conditions for water level
! |uvinit|	Name of file containing initial conditions for velocity

	call addfnm('zinit',' ')
        call addfnm('uvinit',' ')

! |wind|	File with wind data. The file may be either
!		formatted or unformatted. For the format of the unformatted
!		file please see the section where the WIN
!		file is discussed.
!		The format of formatted ASCII file
!		is in standard time-series format, with the
!		first column containing the time in seconds and
!		the next two columns containing the wind data.
!		The meaning of the two values depend on the
!		value of the parameter |iwtype| in the |para| section.
! |qflux|	File with heat flux data. This file must be in
!		a special format to account for the various parameters
!		that are needed by the heat flux module to run. Please
!		refer to the information on the file |qflux|.
! |rain|	File with rain data. This file is a standard time series
!		with the time in seconds and the rain values
!		in mm/day. The values may include also evaporation. Therefore,
!		also negative values (for evaporation) are permitted.
! |ice|		File with ice cover. The values range from 0 (no ice cover)
!		to 1 (complete ice cover).

	call addfnm('wind',' ')
        call addfnm('rain',' ')
        call addfnm('qflux',' ')
        call addfnm('ice',' ')

! |surfvel|	File with surface velocities from observation. These
!		data can be used for assimilation into the model.
! |restrt|	Name of the file if a restart is to be performed. The
!		file has to be produced by a previous run
!		with the parameter |idtrst| different
!		from 0. The data record to be used in the file for the
!		restart must be given by time |itrst|.
! |gotmpa|	Name of file containing the parameters for the
!		GOTM turbulence model (iturb = 1).

        call addfnm('surfvel',' ')
	call addfnm('restrt',' ')
	call addfnm('gotmpa',' ')

! |tempin|	Name of file containing initial conditions for temperature
! |saltin|	Name of file containing initial conditions for salinity
! |conzin|	Name of file containing initial conditions for concentration

        call addfnm('tempin',' ')
        call addfnm('saltin',' ')
        call addfnm('conzin',' ')

! |tempobs|	Name of file containing observations for temperature
! |saltobs|	Name of file containing observations for salinity

        call addfnm('tempobs',' ')
        call addfnm('saltobs',' ')

! |temptau|	Name of file containing the time scale for nudging 
!		of temperature
! |salttau|	Name of file containing the time scale for nudging 
!		of salinity

        call addfnm('temptau',' ')
        call addfnm('salttau',' ')

! |bfmini|	Name of file containing initial conditions for bfm
! |offlin|	Name of the file if a offline is to be performed. The
!		file has to be produced by a previous run
!		with the parameter |idtoff| greater than 0.

        call addfnm('bfmini',' ')
	call addfnm('offlin',' ')

! |fricin|	Name of file containing initial values for friction. The values
!		must be compatible with the friction type given by |ireib|.

        call addfnm('fricin',' ')

! DOCS	END

!c non-documented -> try first	HACK	-> initial conditions

        call addfnm('bioin',' ')
        call addfnm('biosin',' ')
        call addfnm('toxi',' ')
        call addfnm('mercin',' ')	!mercury


! |petsc_zcfg|	Name of file containing the configuration of 
!		PETSc solver for zeta
! |amgx_zcfg|	Name of file containing the configuration of 
!		AmgX solver for zeta

	call addfnm('petsc_zcfg',' ')
	call addfnm('amgx_zcfg',' ')

!c ACQUBC

	call fnm_aquabc_init

!c DOCS	DELWAQ		delwaq model

	call addfnm('voldwq',' ')
	call addfnm('flowdwq',' ')
	call addfnm('areadwq',' ')
	call addfnm('femdwq',' ')
	call addfnm('pntfem',' ')
	call addfnm('ftodwq',' ')

!c DOCS	INTERNAL	internal administration - blank section

	call sctfnm(' ')		!sets default section

	call addfnm('title',' ')

	if( .not. haspar('basnam') ) call addfnm('basnam',' ')
	if( .not. haspar('runnam') ) call addfnm('runnam',' ')
	if( .not. haspar('apnnam') ) call addfnm('apnnam',' ')

	call addfnm('apnfil','apnstd.str')
	call addfnm('memfil','.memory')

	call addfnm('pltfil','plot')

	end

!************************************************************************

        subroutine fnm_aquabc_init

        implicit none

!c for model aquabc (curonian)

        call addfnm('biocon',' ')
        call addfnm('bioscon',' ')
        call addfnm('biolight',' ')
        call addfnm('bioaow',' ')  !ascii output for WC
        call addfnm('bioaos',' ')  !ascii output for BS

        !call addfnm('bioph',' ')
        !call addfnm('biotemp',' ')

        call addfnm('bioload',' ') ! point source loads, not tested yet
        call addfnm('bbs_lev',' ') ! BS depth levels
        call addfnm('settl',' ')   ! settling velocities
        call addfnm('ldbs_m',' ')  ! nutrient load from BS A/D forcing, mud
        call addfnm('ldbs_s',' ')  ! nutrient load from BS A/D forcing, sand
        call addfnm('ldrdx_m',' ') ! nutrient load from WC redox forcing, mud

        end

!************************************************************************

