
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
! subroutine nlsina             initializes the ap parameter file (post)
! subroutine fnmina             initializes default names (post)
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
! 21.10.2023    ggu     only post processing parameters
! 18.09.2024    ggu     new parameters lgrcol, lgrtyp, rfaccol
! 13.11.2025    ggu     new parameter rcircle
!
!************************************************************************

	subroutine nlsina

! initializes the parameter file for the post processing routines

	implicit none

	call nlsina_general
	call nlsina_para
	call nlsina_color
	call nlsina_arrow
	call nlsina_legend
	call nlsina_legvar
	call nlsina_sect

	end

!************************************************************************

	subroutine nlsina_general

! $para section with general variables

	implicit none

! $para section for plotting

! DOCS	START	S_para_general
!
! The next parameters decide what to plot. They can be set directly in
! the STR file here. However, the preferred way to set them is through
! the back end |plots| that automatically sets these parameters.

	call sctpar('para')		!sets default section
	call sctfnm('para')		!sets default section

! The parameter |iwhat| is compulsory and defines the variable to be plotted.
! Please note again that this parameter will be set if invoking the plot
! program through |plots|.
!
! |iwhat|		Flag that determines what to plot. If 0 then
!			the program asks interactively for it. (Default 0)
!			\begin{description}
!			\item[1] Plot basin (grid and isolines of depth)
!			\item[2] Plot velocities
!			\item[3] Plot transports
!			\item[4] Plot water levels
!			\item[5] Plot concentration
!			\item[6] Plot temperature
!			\item[7] Plot salinity
!			\item[8] Plot rms-velocity
!			\item[9] (not used)
!			\item[10] Plot generic scalar (see |ivar|)
!			\item[11] Plot wind vectors
!			\item[12] Plot lagrangian particles
!			\item[13] Plot wave data
!			\end{description}

	call addpar('iwhat',0.)		!what to plot

! The next parameters define other aspects of the plot. In particular their
! meaning is
!
! |iauto|		Normally the simulation name and basin are
!			asked interactively when running the plotting
!			program. However, if |iauto| is different from 0
!			then the file |.memory| is read and the information
!			contained in it is used. The file memory can be
!			set through the program |memory| and it is
!			changed when other parameters are inputted
!			interactively.
! |level|		For 3D applications it indicates the vertical level
!			for which the plot is desired. 1 indicates the 
!			surface layer, 2 the second layer from the surface etc.
!			0 gives integrated quantities, and a value of -1
!			indicates the bottom layer (which refers to different
!			layers for every element). (Default 0)
! |ivar|		Variable to be plotted for scalar quantities. In the
!			file NOS more then one scalars can be contained. To
!			choose the desired scalar, the id of the scalar has 
!			to be given in |ivar|. Note that if only one variable
!			type is contained in the file, then it is not
!			necessary to set |ivar|.

	call addpar('iauto',0.)		!silent mode
	call addpar('level',0.) 	!level (-1 -> bottom   0 -> integr.)
	call addpar('ivar',0.)		!what variable to plot

!c still to document

	call addfnm('varnam',' ')	!variable name to be plotted
	call addpar('isect',0.)		!vertical section

! The next variables define the time of the plots. Even if the names of two
! of the variables are identical to variables used in the finite element
! model, the meaning is different. In the simulation model, |itanf, itend|
! define the initial and final time of simulation, whereas here these
! variables define the initial and final time of the plot of results.
! if none of the time variables are set, then all time records are plotted.

! |itanf|		Initial time for plotting. (Default is 
!			first data record)
! |itend|		Final time for plotting. (Default is 
!			last data record)
! |nout|		Frequency for plotting. A value of 0 or 1
!			plots every record, 2 plots every other record etc.
!			A negative value works as a positive one,
!			but starts plotting from the first record. So
!			-2 plots records 1,3,5,etc.
!			(Default 1) 

	call addpar('itanf',-1.)	!time start
	call addpar('itend',-1.)	!time end
	call addpar('nout',1.)		!time frequency

	call addpar('atanf',-1.)	!time start (absolute)
	call addpar('atend',-1.)	!time end (absolute)

! DOCS	END

	end

!************************************************************************

	subroutine nlsina_para

! $para section

	implicit none

! DOCS	START	S_para_a
!
! The parameters in section |$para| set generic values for the plot.

	call sctpar('para')		!sets default section
	call sctfnm('para')		!sets default section

! Some of the parameters set coordinates in the plot. For example, the
! values |x0, y0| and |x1, y1| indicate the actual plotting area, which can
! be bigger or smaller than the extension of the numerical grid.
!
! Normally, values have to be in meters (the same as the coordinates in the
! numerical grid). However, also relative coordinates can be used. If all
! values given are in the range between -1 and +2, then these values
! are interpreted as relative coordinates. Therefore, $x$ coordinates of
! 0 indicate the left border, and 1 the right border. The upper left quarter
! of the domain can be chosen with (|x0, y0|) = (0,0.5) and
! (|x1, y1|) = (0.5,1).

! |x0, y0|		Lower left corner of the plotting area.
!			(Default is whole area)
! |x1, y1|		Upper right corner of the plotting area.
!			(Default is whole area)

	call addpar('x0',0.)		!dimension of plot
	call addpar('y0',0.)		!dimension of plot
	call addpar('x1',0.)		!dimension of plot
	call addpar('y1',0.)		!dimension of plot

! The next values give the position, where the legend (scale bar and
! true north) is plotted. This legend will only be plotted if
! the coordinates are not geographical (lat/lon) but cartesian.

! |x0leg, y0leg|	Lower left corner of the area
!			where the legend is plotted.
! |x1leg, y1leg|	Upper right corner of the area.
!			where the legend (north and scale) is plotted.

	call addpar('x0leg',0.)		!dimension of legend
	call addpar('y0leg',0.)		!dimension of legend
	call addpar('x1leg',0.)		!dimension of legend
	call addpar('y1leg',0.)		!dimension of legend

! |lblank|		The legend is plotted over a white rectangle.
!			Sometimes this blanking is not desirable.
!			If you do not want to have a white box below the
!			legend set |lblank| to 0. (Default 1)

	call addpar('lblank',1.)

! |cislnd|		It is possible to plot all islands in gray color.
!			Setting |cislnd| to a value between 0 (black) and 
!			1 (white) will achieve this. A negative value
!			will not fill islands with gray color.
!			(Default -1)

	call addpar('cislnd',-1.)	!plot also outer island: 2 <= c <= 3

! |dgray|		It is possible to plot all dry areas in gray color.
!			Setting |dgray| to a value between 0 (black) and 
!			1 (white) will achieve this. A negative value
!			will not fill dry areas with gray color.
!			(Default -1)

! |hgray|		Whereas |dgray| is normally only coloring
!			elements that are dry, you can also color elements
!			shallower than a given depth |hgray|. E.g., a value
!			for |hgray| of -0.5 will plot in gray all
!			elements with depth lower than -0.5 m (salt
!			marshes). (Default -10000)

	call addpar('dgray',-1.)
	call addpar('hgray',-10000.)	!gray all elems with h < hgray

! |dxygrd|		Grid size if the results are interpolated on
!			a regular grid. A value of 0 does
!			not use a regular grid but the original
!			finite element grid for plotting. (Default 0)
! |typls|		Typical length scale to be used when scaling
!			velocity or transport arrows. If |dxygrd| is
!			given this length is used and |typls| is not used.
!			If not given it is computed from the basin
!			parameters. (Default 0)
! |typlsf|		Additional factor to be used with typls to
!			determine the length of the maximum or
!			reference vector. This is the easiest way
!			to scale the velocity arrows 
!			with an overall factor. (Default 1)
! |velref|		Reference value to be used when scaling arrows.
!			If given, a vector with this value will have a length
!			of |typls|*|typlsf| on the map, or, in case
!			|dxygrd| is given, |dxygrd|*|typlsf|. If not set
!			the maximum value of the velocity/transport
!			will be used as |velref|. (Default 0)
! |velmin|		Minimum value for which an arrow will be plotted.
!			With this value you can eliminate small arrows
!			in low dynamic areas. (Default 0)
! |velmax|		Maximum value for which an arrow will be plotted.
!			With this value you can eliminate arrows that are
!			too big. This is useful if you would like to study an
!			area with low current speed but adjacent area have
!			high current speeds that would overplot the area.
!			(Default -1, no limitation)

	call addpar('dxygrd',0.)	!grid size for regular grid
	call addpar('typls',0.)		!typical length scale for arrow plot
	call addpar('typlsf',1.)	!factor for typical length scale
	call addpar('velref',0.)	!reference velocity for length scale
	call addpar('velmin',0.)	!minimum velocity to be plotted
	call addpar('velmax',-1.)	!maximum velocity to be plotted

! |isphe|		If 0 a cartesian coordinate system is used,
!			If 1 the coordinates are in the spherical 
!			system (lat/lon). Among other, this
!			indicates that the $x$-coordinates will be multiplied
!			by a factor that accounts for the visual deformation
!			using lat/lon coordinates.
!			The default is -1, which means that the 
!			type of coordinate system will 
!			be determined automatically. (Default -1)
! |reggrd|		If different from 0 it plots a regular grid over
!			the plot for geographical reference. The value of
!			|reggrd| gives the spacing of the regular grid lines.
!			The units must be according to the units used for
!			the coordinates. With value of -1 the regular grid is
!			determined automatically. (Default -1)
! |regdst|		This value gives the number of intervals
!			that are used to sub-divide the grid given by
!			|reggrd| with a black and white scale around
!			the plot. If 0 it tries to determine automatically
!			the sub-intervals (2 or 4). A value of -1 does
!			not plot the subgrid scale. (Default 0)
! |reggry|		If plotting the regular overlay grid this gives
!			the gray value used for the grid. 0 is black, and
!			1 is white. A value of 1 does not plot the
!			overlay grid, but still writes the labels. 
!			(Default 1)

	call addpar('isphe',-1.)	!spherical coordinate system
	call addpar('reggrd',-1.)	!regular grid spacing
	call addpar('regdst',0.)	!regular micro grid spacing
	call addpar('reggry',1.)	!gray value

! |bndlin|		Name of file that gives the boundary line
!			that is not part of the finite element domain.
!			The file must be in GRD format. An older BND
!			format is also accepted, but deprecated.
!			(Default is no file)

	call addfnm('bndlin'," ")	!name of boundary line file

! |ioverl|		Create overlay of velocity vectors on scalar value.
!			With the value of 0 no overlay is created, 1
!			creates an overlay with the velocity speed.
!			The value of 2 overlays vertical velocities
!			3 water levels and 4 overlays bathymetry.(Default 0)
! |inorm|		Normally the horizontal velocities are plotted
!			in scale. The value of |inorm| can change this
!			behavior. A value of 1 normalizes velocity vectors
!			(all vectors are the same length), whereas 2
!			scales from a given minimum velocity |velmin|.
!			Finally, the value of 3 uses a logarithmic scale.
!			(Default 0)

	call addpar('ioverl',0.)	!overlay in color
	call addpar('inorm',0.)		!vertical velocity as overlay

! The next parameters give the choice to selectively avoid to plot areas
! of the basin and to apply different gray tones for the boundary and
! net lines when plotting the basin.
! Please remember that when working with gray tones the value should
! be between 0 (black) and 1 (white).
!
! |ianopl|	Area code for which no plot has to be produced. Normally 
!		the whole basin is plotted, but with this parameter some
!		areas can be excluded. 
!		Only one of |ianopl| and |iadopl| can be used.
!		(Default -1)
! |iadopl|	Area code for which a plot has to be produced. Normally 
!		the whole basin is plotted, but with this parameter only
!		elements with this area code will be plotted.
!		Only one of |ianopl| and |iadopl| can be used.
!		(Default -1)
! |bgray|	Gray value used for the finite element grid when plotting
!		the bathymetry. (Default 0.8)
! |bbgray|	Gray value used for the boundary of the finite element grid.
!		(Default 0)
! |bsgray|	Gray value used to plot the finite element grid over
!		a scalar or velocity plot. This is basically useful
!		for debugging reasons. The default is to not plot
!		the grid (Default -1.0)

        call addpar('ianopl',-1.)      !do not plot these areas
        call addpar('iadopl',-1.)      !plot only these areas
	call addpar('bgray',0.8)       !gray value for bathymetry
	call addpar('bbgray',0.0)      !gray value for boundary
	call addpar('bsgray',-1.0)     !gray value for plotting maps

! The next two parameters handle the plotting of the lagrangian particles.
!
! |lgrtrj|	If equal 1 plot trajectories
!		instead of particle position. (Default 0)

	call addpar('lgrtrj',0.)	!lagrangian trajectories or position

! |lgmean|	Plot mean positions/trajectories.
!		With the value of 0 no mean pos/traj are created, 1
!		plot mean pos/traj together with single values, 
!		2 plot only mean pos/trajs, 3 as 2 but the first
!		trajectory is plot in thick line. (Default 0)

	call addpar('lgmean',0.)	!lagrangian mean pos/trajs

! The next parameters handle the plotting of the basin.

! |ibox|	If set to 1 plots box information if the area code is
!		used to indicate the box number. This parameter
!		is only useful for the plotting of the box model. (Default 0)
! |inumber|	If set to 1 plots node and element numbers on top
!		of the grid. This is only useful in debug mode. (Default 0)
! |icolmin|	If set to 1 uses minimum number of colors possible.
!		(Default 0)

        call addpar('ibox',0.)
        call addpar('inumber',0.)
        call addpar('icolmin',0.)

! DOCS	END

!c not documented

	call addpar('ifreg',0.)		!plot regular grid from fem file
	call addpar('iexreg',1.)	!plot half valid boxes in reg grid

	call addfnm('obspnt'," ")	!name of file with obs points
	call addfnm('metpnt'," ")	!name of file with meteo points
	call addfnm('spcvel'," ")	!name of file for velocity points

!c only for compatibility ... are not used anymore

	call addpar('traref',0.)	!reference transport for length scale
	call addpar('tramin',0.)	!minimum transport to be plotted

!c internally needed

	call addpar('dirn',0.)
	call addpar('href',0.)
	call addpar('hzmin',0.01)
	!call addpar('hzoff',0.05)
	!call addpar('hlvmin',0.25)

	end

!************************************************************************

	subroutine nlsina_color

! $color section

	use para

	implicit none

! DOCS	START	S_color
!
! Section |$color| deals with the definition of the colors
! to be used in the plot. A color bar is plotted too.

	call sctpar('color')		!sets default section
	call sctfnm('color')		!sets default section

! |icolor|	Flag that determines the type of color table
!		to be used. 0 stands for gray scale, 1 for
!		HSB color table. Other possible values are
!		2 (from white to blue), 3 (from white to red),
!		4 (from blue over white to red) and 
!		5 (from blue over black to red).
!		Values 6 and 7 indicate non-linear HSB color tables.
!		(Default 0)

	call addpar('icolor',0.)	!use color (1.) or not (0.)

! |colfil|	A color table can also be read from file. An example
!		of the format can be found in directory |femplot/color|
!		in the file |colormap.dat|. The variable |colfil|
!		indicates the file where the color table is being
!		read from. The default is not to read a color table file.
! |coltab|	If a color table file has been read then the variable
!		|coltab| indicates the name of the color table that
!		is going to be used. The default is to not use any
!		of the color tables if no name is specified.

        call addfnm('colfil',' ')
        call addfnm('coltab',' ')

! |isoval|	Array that defines the values for the isolines
!		and colors that are to be plotted. Values given must
!		be in the unit of the variable that will be plotted,
!		i.e., meters for water levels etc. 
! |color|	Array that gives the color indices for the
!		plotting color to be used. Ranges are from
!		0 to 1. The type of the color depends on the 
!		variable |icolor|. For the gray scale table
!		0 represents black and 1 white. Values in between
!		correspond to tones of gray. For the HSB color table
!		going from 0 to 1 gives the color of the rainbow.
!		There must be one more value in |color| than in |isoval|.
!		The first color in |color| refers to values less
!		than |isoval(1)|, the second color in |color| to
!		values between |isoval(1)| and |isoval(2)|. The last
!		color in |color| refers to values greater than the last
!		value in |isoval|.

	!call addpar('isoval',0.)
	!call addpar('color',0.)
	call para_add_array_value('isoval',0.)
	call para_add_array_value('color',0.)

! |x0col, y0col|	Lower left corner of the area where the
!			color bar is plotted.
! |x1col, y1col|	Upper right corner of the area where the
!			color bar is plotted.

	call addpar('x0col',0.)		!dimension of color bar
	call addpar('y0col',0.)		!dimension of color bar
	call addpar('x1col',0.)		!dimension of color bar
	call addpar('y1col',0.)		!dimension of color bar

! |cblank|		The color bar is plotted over a white rectangle.
!			Sometimes this blanking is not desirable.
!			If you do not want to have a white box below the
!			legend set |cblank| to 0. (Default 1)

	call addpar('cblank',1.)

! |faccol|	Factor for the values that are written to the 
!		color bar legend. This enables you, e.g., to give water level
!		results in mm (|faccol = 1000|). (Default 1)
! |rfaccol|	Same as |faccol| but inverse factor. This allows you to
!		lower the values easily, e.g., giving times in days
!		instead of seconds (|rfaccol = 86400|). (Default 1)
! |ndccol|	Decimals after the decimal point for the values
!		written to the color bar legend. Use the value |-1|
!		to not write the decimal point. A value of 0 automatically
!		computes the number of decimals needed. (Default 0)
! |legcol|	Text for the description of the color bar. This text
!		is written above the color bar.

	call addpar('faccol',1.)	!factor for color bar
	call addpar('rfaccol',1.)	!inverse factor for color bar
	call addpar('ndccol',-1.)	!decimals after point
	call addfnm('legcol'," ")	!legend for colorbar

! It is not necessary to give all values for isolines and colors above.
! A faster way is to give only the minimum and maximum values and fix
! the number of isovalues to be used.

! |niso|		Total number of isolines to use. (Default is |nisodf|)
! |nisodf|		Default number of isolines to use. (Default 5)
! |colmin, colmax|	Minimum and maximum color index used. Defaults are
!			0.1 and 0.9 respectively. The value of |colmax| can
!			be smaller than |colmin| which inverts the color
!			index used.
! |valmin, valmax|	Minimum and maximum value for isovalues to be used.
!			There is no default.
! |rfiso|		Defines function to be used to compute intermediate
!			values between |valmin| and |valmax|. If 0 or 1
!			the values are linearly interpolated. Else they
!			are computed by $y=x^n$ where $n$ is |rfiso|
!			and $x=\frac{v-v_{min}}{v_{max}-v{min}}$. Values
!			for |rfiso| greater than 0 capture higher detail in
!			the lower values, whereas values less than 1 do
!			the opposite.
!			(Default 0)
! |ipllog|		Indicates the usage of a logarithmic color scale.
!			The possible values are 0-3. The value of 0
!			indicates not to use a logarithmic scale.
!			If 1, the values of
!			the scale are 1,10,100,etc., if 2 the values
!			1,2,10,20,100,etc. are used, and for 3 the values
!			are 1,2,5,10,20,50,100,etc. (Default 0)
! |dval|		Difference of values between isolines. If this
!			value is greater then 0 the values for isolines 
!			and the total number of isolines are computed 
!			automatically using also |valmin| and |valmax|. 
!			(Default 0)

        call addpar('niso',0.)         !total number of isovalues
        call addpar('nisodf',5.)       !default number of isovalues to use
        call addpar('colmin',0.1)      !min color [0..1]
        call addpar('colmax',0.9)      !max color [0..1]
        call addpar('valmin',0.)       !min isovalue
        call addpar('valmax',0.)       !max isovalue
	call addpar('rfiso',0.)	       !function for intermediate values
	call addpar('ipllog',0.)       !logarithmic scale
	call addpar('dval',0.)	       !increment for autom. color sep.

! Since there is a great choice of combinations between the parameters,
! we give here the following rules how the values for colors and isolines
! are determined.
!
! If colors are given in array |color|, they are used, else |colmin| and
! |colmax| or their respective defaults are used to determine the color bar.
! If |isoval| is given it is used, else |valmin| and |valmax| are used.
! If |valmin| and |valmax| are not given they are computed every time
! for each plot and the minimum and maximum value in the basin are used.
! In any case, if |isoval| is specified the total number of isovalues
! is known and |niso| is ignored. However, if |isoval| is not given
! then first |dval| is used to decide how many isovalues to plot, and
! if |dval| is 0 then the |niso| and finally |nisodf| is used.
!
! Other parameters that can be changed are the following.

! |nisomx|	Maximum for |niso| allowed. This is especially useful
!		when the value for |niso| is determined automatically.
!		It avoids you to plot 1000 isolines due to wrong settings
!		of |dval|. However, if you want to use 50 isovalues
!		then just set |niso| and |nisomx| to 50. (Default 20)
! |nctick|	Number of values to be written in color bar. If |niso| is high
!		the labels on the color bar become unreadable. Therefore
!		you can use |nctick| to write only some of the values to
!		the color bar. For example, if |valmin| is 0 and |valmax| is
!		5 and you use many isolines, then setting |nctick| to 6 would
!		give you labels at values 0,1,2,3,4,5. If |nctick| is 0
!		then all labels are written. (Default 0)
! |isolin|	Normally the isolines are not drawn on the plot, just
!		the colors are used to show the value in the different
!		parts of the plot. A value different from 0 plots also 
!		the isolines. In this case |isolin| gives the number of
!		isolines to be plotted. A good choice is to make this
!		equal to |nctick|, so that the isolines correspond to the 
!		values	written on the colorbar. For compatibility, a value of
!		1 plots all isolines. (Default 0)
! |isoinp|	Normally inside elements the values are interpolated
!		(|isoinp=1|).
!		Sometimes it is useful to just plot the value of the
!		node without interpolation inside the element. This can
!		be accomplished by setting |isoinp=0|. Setting instead
!		|isoinp| to a value of 2 plots a constant value in
!		the element. Finally, for |isoinp=3| a circle with radius
!		|rcircle| is plotted around the nodes that contain a value.
!		This is usefull for emphasizing patterns that normally 
!		are too small to be seen.
!		(Default 1)
! |rcircle|	Radius for |isoinp=3|. (Default 0)

        call addpar('nisomx',20.)      !maximum number of isovalues allowed
        call addpar('nctick',0.)       !default number of ticks to use
        call addpar('isolin',0.)       !plot isolines with color ?
        call addpar('isoinp',1.)       !interpolate inside elements
        call addpar('rcircle',0.)      !radius of circle

! Next parameters are for the lagrangian model and the way
! how to plot the particles.

! |lgrtyp|	Type of plot desired. The value of 0 indicates just
!		to plot particles with the same color. A value of 1
!		uses the color information to plot the time that the
!		particles is in the basin. (Default 0)
! |lgrcol|	Color [0-1] to be used to plot the particles when |lgrtyp|
!		is 0. The actual color depends on the color table chosen.
!		(Default 0)

        call addpar('lgrtyp',0.)       !type of plot
        call addpar('lgrcol',0.)       !default color

! DOCS	END

	end

!************************************************************************

	subroutine nlsina_arrow

! $arrow section

	implicit none

! DOCS	START	S_arrow
!
! Parameters in section |$arrow| deal with the reference arrow that is plotted
! in a legend. The arrow regards the plots where the velocity or
! the transport is plotted.

	call sctpar('arrow')		!sets default section
	call sctfnm('arrow')		!sets default section

! |x0arr, y0arr|	Lower left corner of the area where the
!			reference arrow is plotted.
! |x1arr, y1arr|	Upper right corner of the area where the
!			reference arrow is plotted.

	call addpar('x0arr',0.)		!dimension of color bar
	call addpar('y0arr',0.)		!dimension of color bar
	call addpar('x1arr',0.)		!dimension of color bar
	call addpar('y1arr',0.)		!dimension of color bar

! |ablank|		The arrow legend is plotted over a white rectangle.
!			Sometimes this blanking is not desirable.
!			If you do not want to have a white box below the
!			legend set |ablank| to 0. (Default 1)

	call addpar('ablank',1.)

! |facvel|		Factor for the value that is written to the 
!			arrow legend for the velocity.
!			This enables you, e.g., to give 
!			velocities in mm/s (|facvel = 1000|). (Default 1)
! |ndcvel|		Decimals after the decimal point for the values
!			written to the arrow legend. 
!			Use the value |-1|
!			to not write the decimal point. (Default 2)
! |legvel|		Text for the description of the arrow legend.
!			This text is written above the arrow.
! |arrvel|		Length of arrow in legend (in velocity
!			units). If not given the arrow length will be computed
!			automatically. (Default 0)
! |sclvel|		Additional factor to be used for the arrow
!			in the legend. When the arrow length will be
!			computed automatically, this parameter gives
!			the possibility to change the length of the
!			reference vector. This is an easy way
!			to scale the velocity arrow
!			with an overall factor. Not used if
!			|arrvel| is given. (Default 1)

	call addpar('facvel',1.)	!factor for velocity
	call addpar('ndcvel',2.)	!decimals after point (velocity)
	call addfnm('legvel'," ")	!legend for velocity
	call addpar('arrvel',0.)	!length of arrow
	call addpar('sclvel',1.)	!factor for arrow

! DOCS	END

	end

!************************************************************************

	subroutine nlsina_legend

! $legend section

	implicit none

! DOCS	START	S_legend

! In section |$legend| annotations in the plots can be given. The
! section consists of a series of lines that must contain the 
! following information:
!
! The first value is a keyword that specifies what has to be plotted.
! Possible values are |text|, |line|, |vect|, |rect|, |circ|  and also
! |wid|, |col|, and |ctxt|. These correspond to different types of information
! that is inserted into the plot such as text, line, vector, rectangle
! or circle (filled or just outline). Moreover, the color and
! line width of the pen can be controlled by with |wid| and |col|,
! and the centering of text can be set with |ctxt hc vc|, where the values
! -1 mean flush left/bottom (default), 0 mean center, 
! and +1 mean flush right/top.
!
! In case of |text| the starting position (lower left corner) is given,
! then the point size of the font and the text that is inserted. |line|
! needs the starting and end point of the line. The same with |vect|,
! but in this case also the relative tip size must be given as a final
! parameter. |rect| needs the coordinates of the lower left corner and 
! upper right corner of the rectangle. It also needs the color used for
! the filling of the rectangle (0-1) or the flag -1 which only draws the
! outline of the rectangle without filling it. |circ| needs the center
! point and the radius and a fill color (see rectangle). Finally |wid| needs 
! the relative width of the line and |col| the stroke color used when plotting
! lines.
!
! A small example of an annotation that explains the above parameters
! would be:
!
! \begin{verbatim}
! $legend
! text  30500 11800     15  'Chioggia'   #text, 15pt
! line  30500 11800 35000 15000          #line
! vect  30500 11800 35000 15000 0.1      #arrow, tip size 0.1
! rect  30500 11800 35000 15000 0.1      #rectangle, fill color 0.1
! rect  30500 11800 35000 15000 -1       #rectangle (outline, no fill)
! circ  30500 11800 5000 -1              #circle (outline, no fill)
! wid   5                                #set line width to 5
! col   0.5                              #set color to 0.5
! $end
! \end{verbatim}
!
! There is also an old way to specify the legend that does not use
! keywords. However, this way is deprecated and unsupported and is therefore
! not described anymore in this manual.

! DOCS	END

	end

!************************************************************************

	subroutine nlsina_legvar

! $legvar section (variable legend)

	implicit none

! DOCS	START	S_legvar
!
! In section |$legvar| variable fields like the date and wind vectors
! may be inserted into the plot.

	call sctpar('legvar')		!sets default section
	call sctfnm('legvar')		!sets default section

! A time and date can be assigned to the simulation results. These values
! refer to the time 0 of the FEM model. The format for the date is
! YYYYMMDD and for the time HHMMSS. 
! Please note that the date should not be
! given as YYYYMMDD because due to precision problems this will not work. 
! You can also give a time zone if your time is not referring to 
! GMT but to another time zone such as MET. Please note that you have to give
! this information only if the simulation does not contain it already.
! Normally, this information is already assigned during the simulation runs.

! |date|		The real date corresponding to time 0. (Default 0)
! |time|		The real time corresponding to time 0. (Default 0)
! |tz|			The time zone you are in. This is 0 for GMT, 1 for MET
!			and 2 for MEST (MET summer time). (Default 0)
! |tzshow|		The time zone you want to show your results. If your
!			time zone is GMT (0) and you want to show the results
!			referred to MET (+1) set this to +1. Please note that
!			you have to set this variable only if you want to
!			show results in a different time zone than the one
!			given in |tz|. (Default 0)

	call addpar('date',-1.)
	call addpar('time',0.)
	call addpar('tz',0.)
	call addpar('tzshow',0.)

! The information of date and time may be written to the plot. This
! is done with the following parameters.

! |xdate, ydate|	Starting point for the date text (lower left corner).
! |sdate|		Point size of the text. (Default 18)
! |idate|		Output mode. If 0 no date is written to the
!			plot, else the date and time is written. (Default 0)

	call addpar('xdate',0.)
	call addpar('ydate',0.)
	call addpar('sdate',18.)         !size
	call addpar('idate',0.)         !mode

! Wind data can be used to insert a wind vector into the figure.
! This is useful because in the case of variable wind 
! the direction and speed of the wind that was blowing
! in the moment of the plot is shown.
!
! Since only one wind vector can be plotted, the wind data must consist
! of one value for each time. The same ASCII file that is used
! in the STR file can be used.

! |xwind, ywind|	Starting point where the wind arrow is plotted.
! |iwtype|		Type of wind data. The same as the one in the
!			STR file. If this parameter is 0 then no
!			wind vector is plotted. (Default 0)
! |lwwind|		Line width of the wind vector. (Default 0.1)
! |scwind|		Scaling parameter of the wind vector. This depends
!			on the size of your plot. If your wind is 10 m/s
!			and you want the vector to stretch over a distance
!			of 5 km on the plot then you have to choose
!			the value of 500 (10*500=5000) for |scwind|.
!			(Default 1)
! |wfile|		Name of the file containing the wind data. This 
!			may be the same file than the one used in the
!			STR file to run the program.

	call addpar('xwind',0.)
	call addpar('ywind',0.)
	call addpar('iwtype',1.)         !mode
	call addpar('lwwind',0.1)        !line width
	call addpar('scwind',1.)         !scale
	call addfnm('wfile',' ')	 !name of wind file

! The wind vector is also given a text legend with the speed of the
! wind written out. The next parameters decide where and how this information
! is put.

! |xtwind, ytwind|	Starting point for the legend text (lower left corner).
! |stwind|		Point size of the text. (Default 18)
! |wtext|		Text used for the legend (Default 'Wind speed')
! |wunit|		Unit for the wind speed (Default 'm/s')

	call addpar('xtwind',0.)
	call addpar('ytwind',0.)
	call addpar('stwind',18.)         !size
	call addfnm('wtext','Wind speed') !legend for wind
	call addfnm('wunit','m/s')	  !unit for wind

! DOCS	END

	end

!*******************************************************************

	subroutine nlsina_sect

! $sect section (vertical section)

	implicit none

! DOCS	START	S_sect
!
! In this section a vertical section plot can be set up. This will
! allow to plot 3D variables not horizontally for each layer, but
! along a given section in the vertical

	call sctpar('sect')		!sets default section
	call sctfnm('sect')		!sets default section

! |vsect|		Name of the file containing the node list defining
!			the section. The nodes must be adjacent in the
!			numerical grid used. There should be one node
!			on each line. The program |make_line_nodes.sh|
!			can be used to produce such a node list from a
!			given line created with |grid|. If this file
!			is not given, no vertical section will be plotted.

	call addfnm('vsect','')		!name of line defining vertical section

! The next parameters decide about how the plot is scaled in the vertical
! and if a grid overlay is created.
!
! |ivert|		This parameter decides what is plotted on the
!			vertical axis. If 0, depth is plotted on the axis.
!			If 1, vertical layers are plotted. If 2, depth is
!			plotted, but it is scaled with a logarithm, giving
!			more importance to the layers close to the surface.
!			(Default 0)
! |ivgrid|		If 1 a grid is plotted as overlay over the
!			vertical section plot. This is mainly needed
!			for a better orientation where the nodes and 
!			depth values or layers are situated. (Default 0)
! |hvmax, lvmax|	Normally the whole vertical extension is plotted
!			for the section. Sometimes only the upper part
!			of the water column may be needed. With these
!			two parameters the maximum depth (|hvmax|) or
!			the maximum layer (|lvmax|) to be plotted can
!			be imposed. The values may be also higher than
!			the maximum depth/layer. In this case extra space
!			is added below the plotting area. Only one
!			of the two parameters can be set.
! |bsmt|		Set to 1 to have smooth bottom profile even in 
!			case of z-layer structure. (Default 0)

	call addpar('ivert',0.)		!0: depth  1: layers  2: log depth
	call addpar('ivgrid',0.)	!plot grid layout over plot
	call addpar('hvmax',0.)		!maximum depth to be plotted
	call addpar('lvmax',0.)		!maximum layer to be plotted
	call addpar('bsmt',0.)		!1: smooth bottom

! The vertical section plot also creates a color bar. This color bar is
! normally put on the right side of the plot. If there is space inside the
! plot it might be placed on top of the section plot. In this case the
! following parameters can be used. Please note that in this case it is
! the relative position (from 0 to 1) that has to be specified.
!
! |x0sect, y0sect|	Lower left corner of the area where the
!			color bar is plotted.
! |x1sect, y1sect|	Upper right corner of the area where the
!			color bar is plotted.

	call addpar('x0sect',0.)
	call addpar('y0sect',0.)
	call addpar('x1sect',0.)
	call addpar('y1sect',0.)

! The following parameters set titles for the plot, the axis, and for
! an extra description of the starting and end point of the plot.
! They have some (intelligent) default values.
!
! |vtitle|		Title of the plot.
! |xtitle|		Title for the x-axis.
! |ytitle|		Title for the y-axis.
! |ltitle|		Title for start point of the line. (No default)
! |rtitle|		Title for end point of the line. (No default)

	call addfnm('vtitle','Section Plot')	!title for plot
	call addfnm('xtitle','Distance [m]')	!title for x-axis
	call addfnm('ytitle','Depth [m]')	!title for y-axis
	call addfnm('ltitle','')		!title for left point
	call addfnm('rtitle','')		!title for right point

! When plotting velocities you can decide if using for the color the normal 
! velocity across the section or the tangent velocity.
! In any case, in order to visualize also the velocity
! tangent to the section arrwos are used. The next parameters deal with
! the scaling of these arrows.

! |vmode|	When plotting velocities as default the normal velocity 
!		across the section is used for the color plot (|vmode|=0).
!		If you want to use the tangential velocity, please set
!		|vmode|=1. (Default 0)
! |avscal|	This parameter defines the horizontal scale for the 
!		velocity vector. It defines the length scale 
!		in units of x-coordinates of the vertical section so
!		that a velocity of 1 m/s fits comfortably
!		into the plot. If 0 the scale is computed automatically.
!		Please note that in this case the velocities will be
!		scaled differently for every plot. Setting |avscal| 
!		therefore guarantees that the velocity arrows will
!		be comparable between plots. If |avscal| is negative,
!		instead of using x-coordinates units, the legth scale is
!		in [cm]. Therefore, a value of -1 will represent a
!		velocity of 1 m/s with 1 cm on the plot. This unit for
!		the length scale is sometimes easier to understand. 
!		(Default 0)
! |rvscal|	Extra factor that multiplies the scale factor. If your
!		automatic scale gives you vectors which are too small, you
!		can use |rvscal| to increase them. (Default 1)
! |rwscal|	Extra factor for the vertical scale. Normally the
!		vertical scale is computed automatically, If you dont
!		like the size of the vertical vectors you can control
!		it with this parameter. A value of 2 will give you
!		a vertical vector twice as big a the default. (Default 1)
! |dxmin|	Sometimes there are two many arrows plotted horizontally.
!		The value of |dxmin| gives the minimum distance arrows have
!		to be apart to be plotted. A value of 0 plots all arrows.
!		(Default 0)
! |svtip|	The (relative) tip size of the arrow. It specifies how
!		big the arrow will be drawn. A value of 0 only draws the
!		arrow line without tip, and a negative value inhibits
!		drawing arrows at all. (Default 0.3)
! |rxscal,ryscal|	In case arrows are plotted, also a reference
!			vector is plotted. The size of this reference
!			vector s computed automatically, but can be 
!			controlled additionally by the parameters
!			|rxscal,ryscal|, which are in relative units
!			with respect to the reference box plotted.
!			(Default 0.6)

	call addpar('vmode',0.)
	call addpar('avscal',0.)
	call addpar('rvscal',1.)
	call addpar('rwscal',1.)
	call addpar('dxmin',0.)
	call addpar('svtip',0.3)
	call addpar('rxscal',0.6)
	call addpar('ryscal',0.6)

!c still to do: plot vector legend

! DOCS	END

	end

!*******************************************************************
!*******************************************************************
!*******************************************************************
!*******************************************************************
!*******************************************************************

	subroutine fnmina

! initializes default names - not sure if needed

	use para

	implicit none

	logical haspar

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

