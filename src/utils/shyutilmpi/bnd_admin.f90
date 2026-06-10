
!--------------------------------------------------------------------------
!
!    Copyright (C) 1997-2001,2003-2006,2008-2019  Georg Umgiesser
!    Copyright (C) 2005,2014  Christian Ferrarin
!    Copyright (C) 2008  Andrea Cucco
!    Copyright (C) 2012  Aaron Roland
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

! bnd administration routines
!
! contents :
!
! subroutine inbnds
! subroutine rdbnds(ibc)
! subroutine ckbnds
! subroutine prbnds
! subroutine tsbnds
!
! function zvbnds(ibc)			 returns value of open boundary
! subroutine stybnd(ibc,ibtyp)		 sets type of open boundary
! subroutine infobnd(ibc,ibtype,...)	 returns important info on boundary ibc
! function itybnd(ibc)			 returns type of open boundary
! function nbnds			 returns total number of open b.
! function nkbnds(ibc)			 returns total number of nodes of ibc
! function kbnd(i)			 returns i th node of all b. nodes
! subroutine kanfend(ibc,kranf,krend)    returns index of first and last bnode
! function kbnds(ibc,i)			 returns i th node of boundary ibc
! subroutine irbnds(ibc,ndim,idim,nodes) returns nodes of boundary ibc
! subroutine setget_bnd_par(ibc,ientry,value,bset) sets/gets value at ientry
! subroutine setbnd(ibc,value,barray)	 sets boundary ibc to value in barray
!
! subroutine setbc(value,array,flag)	 sets all open boundaries to value
! subroutine chkibc(ibc,errtext)	 checks if ibc is in bounds
!
! notes :
!
! for scalars: -999. uses ambient value
!
! what to do when adding a new parameter to boundary section:
!	add parameter in rdbnds() with description
!	add parameter to array in mod_bnd.f
! what to do when adding a new file name to boundary section:
!	add file name in rdbnds() with description
!	add file name to array in mod_bnd.f
!	if needed add file name in get_boundary_file()
!
! revision log :
!
! 01.08.1997	ggu	$$1stnode - first boundary node was not registered
! 23.09.1997	ggu	introduced conzn -> name of file for conz values
! 30.09.1997	ggu	write error in prbnds
! 03.12.1997	ggu	introduced tempn,saltn (see conzn)
! 27.03.1998	ggu	new utility routines
! 25.05.1998	ggu	documentation (DOCS)
! 20.06.1998	ggu	documentation for momentum input
! 13.07.1998	ggu	implemented ibtyp = 4
! 23.07.1998	ggu	documentation
! 24.08.1998	ggu	BC for concentration is bnd(20,..)
! 24.08.1998	ggu	BC for maximum input level is bnd(12,..) -> levmax
! 27.08.1998	ggu	accessor function levbnd for levmax
! 22.01.1999	ggu	new subroutine setbnd
! 08.07.1999	ggu	bug with iqual -> ibtyp not respected
! 20.01.2000	ggu	call to rdbnds without dimension -> use getdim
! 07.04.2000	ggu	new subroutine setbc
! 07.05.2001	ggu	introduced new variable zfact
! 25.09.2001	ggu	introduced bio2dn
! 07.08.2003	ggu	check for nrb in rdbnds
! 15.10.2004	ggu	new boundary types and sedin
! 02.03.2005	ggu	new nbdim for 3D boundary values
! 02.03.2005	ggu	some new helper functions
! 07.11.2005	ccf	introduced sed2dn
! 16.02.2006	ggu	introduced tox3dn
! 07.04.2008	aac	introduced bfm1bc bfm2bc bfm3bc OB condition for ERSEM
! 17.04.2008	ggu	deleted infobnd(), levbnd()
! 28.04.2008	ggu	call to nrdpar in double precision
! 29.04.2008	ggu&aac	new boundary arrays for ERSEM
! 30.05.2008	ggu	eliminated numbers for parameters
! 03.06.2008	ggu	new parameters levmin, kref
! 06.06.2008	ggu	completely restructured
! 02.04.2009	ggu	intpol default is 0, some unused routines deleted
! 20.04.2009	ggu	new variable ztilt, ndim substituted with nbvdim
! 23.03.2010	ggu	changed v6.1.1
! 28.09.2010	ggu	changed VERS_6_1_11
! 17.02.2011	ggu	changed VERS_6_1_18
! 23.02.2011	ggu	new parameters tramp and levflx implemented
! 01.03.2011	ggu	changed VERS_6_1_20
! 21.06.2012	ggu&aar	new file names for mud module
! 26.06.2012	ggu	changed VERS_6_1_55
! 29.11.2013	ggu	allow for non continous boundary numbering
! 05.12.2013	ggu	changed VERS_6_1_70
! 28.01.2014	ggu	changed VERS_6_1_71
! 28.03.2014	ggu	new parameter lgrpps
! 05.05.2014	ggu	changed VERS_6_1_74
! 16.06.2014	ggu	new include file bnd.h (with nbndim)
! 25.06.2014	ggu	new routine exists_bnd_name()
! 21.10.2014	ggu	changed VERS_7_0_3
! 29.10.2014	ccf	include vel3dn boundary file
! 03.11.2014	ggu	nbdim deleted
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 23.06.2015	ggu	setbc() deleted, nrz,nrq eliminated
! 10.07.2015	ggu	changed VERS_7_1_50
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 30.07.2015	ggu	changed VERS_7_1_83
! 05.11.2015	ggu	changed VERS_7_3_12
! 18.12.2015	ggu	changed VERS_7_3_17
! 14.01.2016	ggu	check of boundaries considers mpi subdomains
! 15.02.2016	ggu	check if boundary is given twice
! 19.02.2016	ggu	changed VERS_7_5_3
! 22.02.2016	ggu	new files bfmbcn integrated
! 01.04.2016	ggu	restructured - arrays transfered to mod_bnd.f
! 15.04.2016	ggu	changed VERS_7_5_8
! 07.06.2016	ggu	changed VERS_7_5_12
! 30.09.2016	ggu	changed VERS_7_5_18
! 31.03.2017	ggu	changed VERS_7_5_24
! 13.04.2017	ggu	use array feature of para for kbound
! 05.12.2017	ggu	changed VERS_7_5_39
! 03.04.2018	ggu	changed VERS_7_5_43
! 19.04.2018	ggu	changed VERS_7_5_45
! 25.10.2018	ggu	changed VERS_7_5_51
! 18.12.2018	ggu	changed VERS_7_5_52
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 20.03.2023	ggu	relax condition of no holes for ibtyp == 3
! 24.05.2023	ggu	in ckbnds() more debug, new boundary_debug()
! 03.12.2024	ggu	do not complain for holes in zeta boundary
! 25.02.2026	ggu	updated documentation for nad
!
!************************************************************************

	subroutine inbnds

! initializes boundary parameters

	use mod_bnd

	implicit none

	nbc = 0
	nrb = 0

	end

!************************************************************************

	subroutine rdbnds(ibc)

! reads boundary info from STR file

	use mod_bnd
	use mod_bound_geom
	use nls
	use para

	implicit none

	integer ibc

	character*80 name,text,file
	double precision dvalue
	real value
	integer n
	integer i,kranf,krend,kref
	integer iweich,id,nbnd
	integer nrdpar
	logical ball

	real getpar

!------------------------------------------------------
! add new elements to boundary arrays
!------------------------------------------------------

	n = max(nbc,ibc)	!allow for non contiguous numbering
	call mod_bnd_adjust(n)
	do i=nbc+1,n
	  call init_bnd_par(i)
	  call init_bnd_file(i)
	end do
	nbc = n

!------------------------------------------------------
! check if boundary section is unique
!------------------------------------------------------

        call get_bnd_ipar(ibc,'kranf',kranf)
	if( kranf > 0 ) then
	  write(6,*) 'ibc = ',ibc
	  stop 'error stop rdbnds: boundary defined twice'
	end if

!------------------------------------------------------
! start initializing boundary section
!------------------------------------------------------

	call sctpar('bound')
	call sctfnm('bound')

! DOCS	START	S_bound
!
! These parameters determine the open boundary nodes and the type
! of the boundary: level or flux boundary. At the first the water levels
! are imposed, on the second the fluxes are prescribed.
!
! There may be multiple sections |bound| in one parameter input file,
! describing all open boundary conditions necessary. Every section
! must therefore be supplied with a boundary number. The numbering
! of the open boundaries must
! be increasing. The number of the boundary must be specified
! directly after the keyword |bound|, such as |bound1| or |bound 1|.
!
! |kbound|	Array containing the node numbers that are part of the
!		open boundary. The node numbers must form one contiguous
!		line with the domain (elements) to the left. This
!		corresponds to an anti-clockwise sense. The type of
!		boundary depends on the	value of |ibtyp|. In case this value
!		is 1 or 2 at least two nodes must be given.

	call para_add_array_value('kbound',0.)

! |ibtyp|	Type of open boundary. 
!		\begin{description}
!		\item[0] No boundary values specified
!		\item[1] Level boundary. At this open boundary
!			 the water level is imposed and the prescribed
!			 values are interpreted as water levels in meters.
!			 If no value for |ibtyp| is specified this
!			 is the default.
!		\item[2] Flux boundary. Here the discharge in \dischargeunit
!			 has to be prescribed.
!		\item[3] Internal flux boundary. As with |ibtyp = 2| a
!			 discharge has to be imposed, but the node where
!			 discharge is imposed can be an internal node
!			 and need not be on the outer boundary of
!			 the domain. For every node in |kbound| the
!			 volume rate specified will be added to the
!			 existing water volume. This behavior is different
!			 from the |ibtyp = 2| where the whole boundary
!			 received the discharge specified.
!		\item[4] Momentum input. The node or nodes may be internal.
!			 This feature can be used to describe local 
!			 acceleration of the water column. 
!			 The unit is force / density [\maccelunit].
!			 In other words it is the rate of volume 
!			 [\dischargeunit] times the velocity [m/s] 
!			 to which the water is accelerated.
!		\end{description}
! |iqual|	If the boundary conditions for this open boundary
!		are equal to the ones of boundary |i|, then
!		setting |iqual = i| copies all the values of
!		boundary |i| to the actual boundary. Note that the
!		value of |iqual| must be smaller than the number
!		of the actual boundary, i.e., boundary |i| must have
!		been defined before. (This feature is temporarily
!		not working; please do not use.)

	call addpar('ibtyp',1.)
	call addpar('iqual',0.)

! The next parameters give a possibility to specify the file name
! of the various input files that are to be read by the model.
! Values for the boundary condition can be given at any time step.
! The model interpolates in between given time steps if needed. The
! grade of interpolation can be given by |intpol|.
!
! All files are in ASCII and share a common format.
! The file must contain two columns, the first giving the
! time of simulation in seconds that refers to the value
! given in the second column. The value in the second
! column must be in the unit of the variable that is given.
! The time values must be in increasing order.
! There must be values for the whole simulation,
! i.e., the time value of the first line must be smaller
! or equal than the start of the simulation, and the time
! value of the last line must be greater or equal than the
! end of the simulation.

! |boundn|	File name that contains values for the boundary condition.
!		The value of the variable given in the second column
!		must be in the unit determined by |ibtyp|, i.e.,
!		in meters for a level boundary, in \dischargeunit for
!		a flux boundary and in \maccelunit for a momentum
!		input.
! |zfact|	Factor with which the values from |boundn|
!		are multiplied to form the final value of the
!		boundary condition. E.g., this value can be used to
!		set up a quick sensitivity run by multiplying
!		all discharges by a factor without generating
!		a new file. (Default 1)

	call addfnm('boundn',' ')
	call addpar('zfact',1.)

! |levmin, levmax|	A point discharge normally distributes its discharge
!		over the whole water column. If it is important that in
!		a 3D simulation the water mass discharge is concentrated 
!		only in some levels, the parameters |levmin| and |levmax|
!		can be used. They indicate the lowest and deepest level over
!		which the discharge is distributed. Default values are 0, which
!		indicate that the discharge is distributed over the
!		whole water column. Setting only |levmax| distributes from
!		the surface to this level, and setting only |levmin|
!		distributes from the bottom to this level.

	call addpar('levmin',0.)
	call addpar('levmax',0.)

! |conzn, tempn, saltn|	File names that contain values for the respective
!			boundary condition, i.e., for concentration,
!			temperature and salinity. The format is the same
!			as for file |boundn|. The unit of the values
!			given in the second column must the ones of the
!			variable, i.e., arbitrary unit for concentration,
!			centigrade for temperature and psu (per mille)
!			for salinity.

	call addfnm('conzn',' ')
	call addfnm('tempn',' ')
	call addfnm('saltn',' ')

! |vel3dn|	File name that contains current velocity values for the 
!		boundary condition.  The format is the same as for file 
!		|tempn| but it has two variables: 
!	 	current velocity in x and current velocity in y.
!		Velocity can be nudged or imposed depending on the value 
!		of |tnudge| (mandatory). The unit is [m/s].

	call addfnm('vel3dn',' ')

! |tnudge|	Relaxation time for nudging of boundary velocity.
!		For |tnudge| = 0. velocities are imposed, for
!		|tnudge| > 0. velocities are nudged. The
!		default is -1 which means do nothing. Unit is [s].
!		(Default -1)

        call addpar('tnudge',-1.)

! The next variables specify the name of the boundary value file
! for different modules. Please refer to the documentation of the
! single modules for the units of the variables.

! |bio2dn|	File name that contains values for the ecological
!		module (EUTRO-WASP).
! |sed2dn|	File name that contains values for the sediment
!		transport module.
!		The unit of the values given
!		in the second and following columns (equal to the 
!		number of defined grainsize in parameter |sedgrs|).

! |mud2dn|	File name that contains values for the fluid mud
!		module.
! |lam2dn|	File name that contains values for the fluid mud
!		module (boundary condition for the structural parameter, 
!		to be implemented).
! |dmf2dn|	File name that contains values for the fluid mud
!		module (boundary conditions for the advection of flocsizes,
!		to be implemented).
! |tox3dn|	File name that contains values for the toxicological
!		module.

	call addfnm('bio2dn',' ')
	call addfnm('sed2dn',' ')
	call addfnm('mud2dn',' ')
	call addfnm('lam2dn',' ')
	call addfnm('dmf2dn',' ')
	call addfnm('tox3dn',' ')

!c File name for OB condition in ERSEM MODULE - undocumented
!c ... will be removed sooner or later ...

	call addfnm('bfm1bc',' ')
	call addfnm('bfm2bc',' ')
	call addfnm('bfm3bc',' ')

! |bfmbcn|	File name that contains values for the bfm module.

	call addfnm('bfmbcn',' ')

! |mercn|	File name that contains values for the mercury module.

	call addfnm('mercn',' ')

! |s4mern|      File name that contains values for the mercury module.

        call addfnm('s4mern',' ')


! |intpol|	Order of interpolation for the boundary values read
!		in files. Use for 1 for stepwise (no) interpolation,
!		2 for linear and 4 for cubic interpolation. 
!		The default is linear interpolation, except for
!		water level boundaries (|ibtyp=1|) where cubic
!		interpolation is used.

	call addpar('intpol',0.)

! The next parameters can be used to impose a sinusoidal water level
! (tide) or flux at the open boundary. These values are used if no
! boundary file |boundn| has been given. The values must be in the unit
! of the intended variable determined by |ibtyp|.

! |ampli|	Amplitude of the sinus function imposed. (Default 0)
! |period|	Period of the sinus function. (Default 43200, 12 hours)
! |phase|	Phase shift of the sinus function imposed. A positive value
!		of one quarter of the period reproduces a cosine
!		function. (Default 0)
! |zref|	Reference level of the sinus function imposed. If only
!		|zref| is specified (|ampli = 0|) a constant value
!		of |zref| is imposed on the open boundary.

	call addpar('ampli',0.)
	call addpar('period',43200.)
	call addpar('phase',0.)
	call addpar('zref',0.)

! With the next parameters a constant value can be imposed for the 
! variables of concentration, temperature and salinity. In this case
! no file with boundary values has to be supplied. The default for all
! values is 0, i.e., if no file with boundary values is supplied and
! no constant is set the value of 0 is imposed on the open boundary.
! A special value of -999 is also allowed. In this case the value
! imposed is the ambient value of the parameter close to the boundary.

! |conz, temp, salt|	Constant boundary values for concentration,
!			temperature and salinity respectively. If these
!			values are set no boundary file has to be supplied.
!			(Default 0)

	call addpar('conz',0.)			!$$conz
	call addpar('temp',0.)			!$$baroc
	call addpar('salt',0.)			!$$baroc

! The next two values are used for constant momentum input. 
! This feature can be used to describe local acceleration of the
! water column. The values give the input of momentum 
! in x and y direction. The unit is force / density (\maccelunit).
! In other words it is the rate of volume (\dischargeunit) times
! the velocity (m/s) to which the water is accelerated.
!
! These values are used if 
! boundary condition |ibtyp = 4| has been chosen and
! no boundary input file has been given.
! If the momentum input is varying then it may be specified with
! the file |boundn|. In this case the file |boundn| must contain
! three columns, the first for the time, and the other two for
! the momentum input in $x,y$ direction.
!
! Please note that this feature is temporarily not available.
!
! |umom, vmom|		Constant values for momentum input. (Default 0)

	call addpar('umom',0.)
	call addpar('vmom',0.)

! The next two values can be used
! to achieve the tilting of the open boundary if only one water level value
! is given. If only |ktilt| is given then the boundary values
! are tilted to be in equilibrium with the Coriolis force. This may avoid
! artificial currents along the boundary. |ktilt| must be a boundary node
! on the boundary.
!
! If |ztilt| is given the tilting of the boundary is explicitly set
! to this value. The tilting of the first node of the boundary is set 
! to $-|ztilt|$
! and the last one to $+|ztilt|$. The total amount of tilting is
! therefore is $2 \cdot |ztilt|$. If |ktilt| is not specified
! then a linear interpolation between the first and the last boundary
! node will be carried out. If also |ktilt| is specified then
! the boundary values are arranged that the water levels are 
! tilted around |ktilt|, e.g., $-|ztilt|$ at the first boundary node,
! 0 at |ktilt|, and $+|ztilt|$ at the last boundary node.
!
! |ktilt|		Node of boundary around which tilting should
!			take place. (Default 0, i.e., no tilting)
! |ztilt|		Explicit value for tilting (unit meters).
!			(Default 0)

	call addpar('ktilt',0.)
	call addpar('ztilt',0.)

! Other parameters:

! |igrad0|		If different from 0 a zero gradient boundary
!			condition will be implemented. This is already the
!			case for scalars under outflowing conditions. However,
!			with |igrad0| different from 0 this conditions
!			will be used also for inflow conditions. (Default 0)

	call addpar('igrad0',0.)	!use 0 gradient for scalars

! |tramp|		Use this value to start smoothly a discharge
!			boundary condition. If set it indicates the
!			time (seconds) that will be used to increase
!			a discharge from 0 to the desired value (Default 0)

	call addpar('tramp',0.)		!start smoothly for discharge

! |levflx|		If discharge is depending on the water level
!			(e.g., lake outflow) then this parameter indicates to
!			use one of the possible outflow curves. Please
!			note that the flow dependence on the water level
!			must be programmed in the routine $|level\_flux()|$.
!			(Default 0)

	call addpar('levflx',0.)	!use level-discharge relationship

! |nad|			On the open boundaries it is sometimes convenient
!			to not compute the non-linear terms in the momentum
!			equation because instabilities may occur. Setting 
!			the parameter |nad| to a value different greater than 0
!			indicates that in the first |nad| nodes from the
!			boundary the non linear terms are switched off. For
!			|nad=0| all terms are computed for the boundary. For
!			the default value of -1 the value of |nad| is taken
!			from the global value of |nadist|.
!			(Default -1)

	call addpar('nad',-1.)		!no advective terms for this boundary

! |lgrpps|		Indicates the number of particles released at
!			the boundary for the lagrangian module. If positive
!			it is the number of particles per second released
!			along the boundary. If negative its absolute
!			value indicates the particles per volume flux
!			(unit \dischargeunit) released along the boundary.
!			(Default 0)

	call addpar('lgrpps',0.)	!particles per second for lagrange
					!if negative parts per volume flux

! DOCS	END

!c undocumented
	call addpar('kref',0.)		!not working...
! here add dummy variables
	call addpar('zval',0.)
	call addpar('kmanf',0.)
	call addpar('kmend',0.)

!------------------------------------------------------
! start reading loop
!------------------------------------------------------

	kranf = nrb + 1
	call set_bnd_ipar(ibc,'kranf',kranf)	!position of starting node
	call addpar('kranf',float(kranf))

	iweich=1
	do while(iweich.ne.0)
	    iweich=nls_insert_variable('bound',name,dvalue,text)
	    value = dvalue
	    if( iweich .lt. 0 ) goto 92
	    if( iweich .eq. 2 .and. name .ne. 'kbound' ) goto 93
	    if( iweich .eq. 4 ) goto 93
	    if( name .eq. 'kbound' ) then	!$$1stnode
		!nrb=nrb+1
		!call mod_irv_init(nrb)
		!irv(nrb)=nint(value)
	    else if( iweich .eq. 3 ) then	!file name
                call set_bnd_file(ibc,name,text)
	    !else if( iweich .ne. 0 ) then
	    else if( iweich .eq. 1 ) then
                call set_bnd_par(ibc,name,value)
	    end if
	end do

	call para_get_array_size('kbound',n)
	call mod_irv_init(nrb+n)
	call para_get_array_value('kbound',n,n,irv(nrb+1:))
	nrb=nrb+n

	krend = nrb
	call set_bnd_ipar(ibc,'krend',krend)	!position of end node
	call addpar('krend',float(krend))

!------------------------------------------------------
! copy all values and files to private boundary arrays
!------------------------------------------------------

	call copy_bnd_par(ibc)
	call copy_bnd_file(ibc)

!------------------------------------------------------
! write out all values of boundary section (only for debug)
!------------------------------------------------------

	ball = .true.
        !call check_bnd_par_entries(ibc,ball)
        !call check_bnd_file_entries(ibc,ball)

!------------------------------------------------------
! delete public section
!------------------------------------------------------

	call check_parameter_values('before deleting section')
	call delete_section('bound')
	call check_parameter_values('after deleting section')

!------------------------------------------------------
! end of routine
!------------------------------------------------------

	return
   92	continue
	write(6,*) 'Error in name list read : ',iweich
	stop 'error stop : rdbnds'
   93	continue
	write(6,*) 'Variable not allowed in vector context : ',name
	stop 'error stop : rdbnds'
   94   continue
        write(6,*) 'Boundary conditions out of order : ',ibc
        stop 'error stop : rdbnds'
	end

!********************************************************************

	subroutine ckbnds

! checks boundary information read from STR

	use mod_bnd
	use mod_bound_geom
	use shympi

	implicit none

	logical bstop
	integer istop
	integer i,k,ibc
	integer iqual,ibtyp,kranf,krend,ktilt,knode,kref
	integer levmax,levmin
	integer intpol
	integer kmanf,kmend,krtot,kmtot
	real period
	real ztilt
	integer ipint
	character*80 file

	bstop = .false.

	call shympi_syncronize

	do i=1,nbc

	 ibc = i

         call get_bnd_ipar(ibc,'iqual',iqual)
	 if(iqual.ge.i) then
	   write(6,'(a,i2,a)') 'Section BOUND ',i,' :'
	   write(6,*) '   Value for iqual must be'
	   write(6,*) '   smaller than actual boundary number'
	   write(6,*) '   iqual = ',iqual
	   bstop=.true.
	 else if(iqual.lt.0) then
	   write(6,'(a,i2,a)') 'Section BOUND ',i,' :'
	   write(6,*) '   Value for iqual must be'
	   write(6,*) '   greater than 0'
	   write(6,*) '   iqual = ',iqual
	   bstop=.true.
	 else if( iqual .gt. 0 ) then
           call get_bnd_ipar(iqual,'ibtyp',ibtyp)
           call set_bnd_ipar(ibc,'ibtyp',ibtyp)
           call get_bnd_ipar(iqual,'levmax',levmax)
           call set_bnd_ipar(ibc,'levmax',levmax)
           call get_bnd_ipar(iqual,'levmin',levmin)
           call set_bnd_ipar(ibc,'levmin',levmin)
	 end if

         call get_bnd_ipar(ibc,'ibtyp',ibtyp)
         if(ibtyp.ge.0.and.ibtyp.le.3) then !$$ibtyp3
	 else if(ibtyp.eq.4) then		 !$$ibtyp11
	 else if(ibtyp.eq.5) then		 !$$ibtyp11
	 else if(ibtyp.eq.11) then		 !$$ibtyp11
	 else if(ibtyp.ge.30.and.ibtyp.le.33) then	 !$$ibtyp11
	 else if(ibtyp.ge.50.and.ibtyp.le.53) then	 !$$ibtyp11
	 else if(ibtyp.eq.70) then	 !$$roger
	 else
	   write(6,'(a,i2,a)') 'Section BOUND ',i,' :'
	   write(6,*) '   Value not allowed for ibtyp'
	   write(6,*) '   ibtyp = ',ibtyp
	   bstop=.true.
	 end if

	 if( ibtyp == 1 ) then
	   call get_boundary_file(ibc,'zeta',file)
           call get_bnd_par(ibc,'period',period)
	   if( period .le. 0. .and. file .eq. ' ' ) then
		write(6,'(a,i2,a)') 'section BOUND ',i,' :'
		write(6,*) '   Period must be > 0'
		write(6,*) '   period = ',period
		bstop=.true.
	   end if
	 end if

         call get_bnd_ipar(ibc,'ktilt',ktilt)
	 if( ktilt .gt. 0 ) then
	   knode=ipint(ktilt)		!$$EXTINW
	   if(knode.le.0) then
		write(6,'(a,i2,a)') 'Section BOUND ',i,' :'
		write(6,*) '   tilt node not found ',ktilt
		bstop=.true.
	   end if
           call set_bnd_ipar(ibc,'ktilt',knode)
	 end if

         call get_bnd_ipar(ibc,'kref',kref)
	 if( kref .gt. 0 ) then
	   knode=ipint(kref)		!$$EXTINW
	   if(knode.le.0) then
		write(6,'(a,i2,a)') 'Section BOUND ',i,' :'
		write(6,*) '   reference node not found ',kref
		bstop=.true.
	   end if
           call set_bnd_ipar(ibc,'kref',knode)
	 end if

         call get_bnd_ipar(ibc,'intpol',intpol)
	 if(intpol.lt.0.or.intpol.gt.4) then
		write(6,'(a,i2,a)') 'Section BOUND ',i,' :'
		write(6,*) '   erroneous value for INTPOL : ',intpol
		bstop=.true.
	 end if

         call get_bnd_ipar(ibc,'kranf',kranf)
         call get_bnd_ipar(ibc,'krend',krend)
	 if( kranf .gt. krend ) then	!$$kranf
	   write(6,'(a,i2,a)') 'section BOUND ',i,' :'
	   write(6,*) '   No nodes given for boundary'
	   bstop=.true.
	 end if

	 kmanf = 0
	 kmend = 0
	 krtot = krend-kranf+1

	 istop = 0
	 do k=kranf,krend
	   if( k == 0 ) cycle
	   knode=ipint(irv(k))		!$$EXTINW
	   if(knode.le.0) then
             if( .not. bmpi ) then
               write(6,'(a,i2,a)') 'Section BOUND ',i,' :'
               write(6,*) '   boundary node not found ',irv(k)
	     end if
           else if( .not. shympi_is_inner_node(knode) ) then
             !knode = 0		!we keep ghost node
	   end if
	   if( knode > 0 ) then
	     if( kmanf == 0 ) kmanf = k
	     kmend = k
	   else
	     istop = istop + 1
	   end if
	   irv(k)=knode
	   !write(6,*) k,knode,kmanf,kmend,istop
	 end do

	 kmtot = kmend-kmanf+1
         call set_bnd_ipar(ibc,'kmanf',kmanf)
         call set_bnd_ipar(ibc,'kmend',kmend)
!         write(6,'(a,10i5)') 'boundary: ',my_id,ibc,istop
!      &			,kranf,krend,kmanf,kmend

         if( istop > 0 ) then
           write(6,'(a,10i5)') 'boundary: ',my_id,ibc,istop          &
      &			,kranf,krend,kmanf,kmend
           if( shympi_is_parallel() ) then
             if( istop == krtot ) then
               !write(6,*) 'boundary not in this domain... ok'
               call set_bnd_ipar(ibc,'ibtyp',0)
             else if( istop == krtot-kmtot ) then
               write(6,*) 'boundary in more than one domain... ok'
             else
	       !if( ibtyp == 1 .or. ibtyp == 2 ) then	!no holes allowed
	       if( ibtyp == 2 ) then	!no holes allowed
	         write(6,*) my_id,krtot,kmtot,istop,krtot-kmtot
	         write(6,*) '*** there are holes in the node list'
	         call boundary_debug(ibc)
		 call error_stop('ckbnds','internal error (1)')
	       end if
             end if
           else
	     call error_stop('ckbnds','no MPI and missing nodes')
           end if
         end if

	end do

	call shympi_syncronize
	if( bstop ) call error_stop('ckbnds','general error')

	end

!********************************************************************

	subroutine boundary_debug(ibc)

	use mod_bnd
	use mod_bound_geom
	use shympi

	implicit none

	integer ibc

	integer k,kranf,krend,kmanf,kmend,knode
	integer id,iu
	integer ipint

	iu = 600 + my_id

	write(iu,*) '==============================='
	write(iu,*) 'boundary_debug: ',my_id
	write(iu,*) '==============================='

         call get_bnd_ipar(ibc,'kranf',kranf)
         call get_bnd_ipar(ibc,'krend',krend)
         call get_bnd_ipar(ibc,'kmanf',kmanf)
         call get_bnd_ipar(ibc,'kmend',kmend)

	write(iu,*) 'boundary_debug: ',kranf,krend
	write(iu,*) 'boundary_debug: ',kmanf,kmend

	do k=kranf,krend
	   if( k == 0 ) cycle
	   knode=irv(k)		!$$EXTINW
	   id = -1
	   if( knode > 0 ) id = id_node(knode)
	   write(iu,*) k,knode,id
	end do

	write(iu,*) '==============================='
	write(iu,*) 'end of boundary_debug: ',my_id
	write(iu,*) '==============================='

	flush(6)

	end

!********************************************************************

	subroutine prbnds

	use mod_bnd
	use mod_bound_geom
	use mod_info_output

	implicit none

	integer i,ibc
	integer ibtyp,kranf,krend
	integer intpol
	integer ipext
	logical, parameter :: ball = .false.	!write all file names

	if( nbc .le. 0 ) return

	if( print_verbose() ) then

	write(6,*)
	write(6,*) '====== info on open boundaries ========='
	write(6,*)

	do ibc=1,nbc
          call get_bnd_ipar(ibc,'kranf',kranf)
          call get_bnd_ipar(ibc,'krend',krend)

	  write(6,'(a,i9)') 'boundary: ',ibc
	  if( kranf > 0 .and. krend > 0 ) then
	    write(6,'(a,i9)') ' boundary nodes: ',krend-kranf+1
	    write(6,'(8i9)') (ipext(irv(i)),i=kranf,krend)
	  end if

	  call check_bnd_par_entries(ibc,ball)
	  call check_bnd_file_entries(ibc,ball)
	end do

	write(6,*)
	write(6,*) '==== end info on open boundaries ======='
	write(6,*)

	end if

	return
	end

!********************************************************************

	subroutine tsbnds

	use mod_bnd
	use mod_bound_geom

	implicit none

	integer j,ibc
	logical, parameter :: ball = .true.	!write all file names

	write(6,*) '/bnd/'
	write(6,*) nbc
	do ibc=1,nbc
	  call check_bnd_par_entries(ibc,ball)
	  call check_bnd_file_entries(ibc,ball)
	end do

	write(6,*) '/irv/'
	write(6,*) nrb
	do j=1,nrb
	  write(6,*) irv(j)
	end do

	end

!********************************************************************
!********************************************************************
!********************************************************************
!      utility routines
!********************************************************************
!********************************************************************
!********************************************************************

	function zvbnds(ibc)

! returns value of open boundary

	implicit none

	real zvbnds
	integer ibc

	real zval

        call chkibc(ibc,'zvbnds:')

        call get_bnd_par(ibc,'zval',zval)
	zvbnds = zval

	end

!********************************************************************

	subroutine setbnds(ibc,zval)

! sets value of open boundary

	implicit none

	integer ibc
	real zval

        call chkibc(ibc,'setbnds:')

        call set_bnd_par(ibc,'zval',zval)

	end

!********************************************************************

	subroutine stybnd(ibc,ibtyp)

! sets type of open boundary

	implicit none

	integer ibc
	integer ibtyp

        call chkibc(ibc,'stybnd:')

        call set_bnd_ipar(ibc,'ibtyp',ibtyp)

	end

!********************************************************************

	function itybnd(ibc)

! returns type of open boundary

	implicit none

	integer itybnd
	integer ibc

	integer ibtyp

        call chkibc(ibc,'itybnd:')

        call get_bnd_ipar(ibc,'ibtyp',ibtyp)
	itybnd = ibtyp

	end

!********************************************************************

	function nbnds()

! returns total number of open boundaries

	use mod_bnd

	implicit none

	integer nbnds


	nbnds = nbc

	end

!********************************************************************

	function nkbnd()

! returns total number of open boundary nodes

	use mod_bnd

	implicit none

	integer nkbnd


	nkbnd = nrb

	end

!********************************************************************

	function nkbnds(ibc)

! returns total number of nodes of boundary ibc

	implicit none

	integer nkbnds
	integer ibc

	integer kranf,krend

        call chkibc(ibc,'nkbnds:')

        call get_bnd_ipar(ibc,'kranf',kranf)
        call get_bnd_ipar(ibc,'krend',krend)

	nkbnds = krend - kranf + 1
	if( kranf == 0 .or. krend == 0 ) nkbnds = 0

	end

!********************************************************************

	function kbnd(i)

! returns i th node of all boundary nodes

	use mod_bnd
	use mod_bound_geom

	implicit none

	integer kbnd
	integer i

        if( i .gt. nrb ) then
            write(6,*) 'i, nrb, nbc : ',i,nrb,nbc
            stop 'error stop kbnd: i out of bounds'
        end if

	kbnd = irv(i)

	end

!********************************************************************

	subroutine kanfend(ibc,kranf,krend)

! returns index of first and last boundary node of boundary ibc

	implicit none

	integer ibc
        integer kranf,krend

        call chkibc(ibc,'kanfend:')

        call get_bnd_ipar(ibc,'kranf',kranf)
        call get_bnd_ipar(ibc,'krend',krend)

	end

!********************************************************************

	subroutine kmanfend(ibc,kmanf,kmend)

! returns index of first and last boundary node of boundary ibc
! mpi version

	implicit none

	integer ibc
        integer kmanf,kmend

        call chkibc(ibc,'kmanfend:')

        call get_bnd_ipar(ibc,'kmanf',kmanf)
        call get_bnd_ipar(ibc,'kmend',kmend)

	end

!********************************************************************

	function kbndind(ibc,i)

! returns global index of i th node of boundary ibc

	implicit none

	integer kbndind
	integer ibc
	integer i

        integer kranf,krend,n

        call chkibc(ibc,'kanfend:')

        call get_bnd_ipar(ibc,'kranf',kranf)
        call get_bnd_ipar(ibc,'krend',krend)

	n = krend - kranf + 1

	if( i .lt. 1 .or. i .gt. n ) then
	    write(6,*) 'i, imax, ibc : ',i,n,ibc
	    stop 'error stop kbndind: i out of bounds'
	end if

	kbndind = kranf + i - 1

	end

!********************************************************************

	function kbnds(ibc,i)

! returns i th node of boundary ibc

	use mod_bound_geom

	implicit none

	integer kbnds
	integer ibc
	integer i

	integer kranf,krend,n

        call chkibc(ibc,'kbnds:')

        call get_bnd_ipar(ibc,'kranf',kranf)
        call get_bnd_ipar(ibc,'krend',krend)

	n = krend - kranf + 1

	if( i .lt. 1 .or. i .gt. n ) then
	    write(6,*) 'i, imax, ibc : ',i,n,ibc
	    stop 'error stop kbnds: i out of bounds'
	end if

	kbnds = irv(kranf+i-1)

	end

!********************************************************************

	subroutine irbnds(ibc,ndim,idim,nodes)

! returns nodes of boundary ibc (maximum ndim)

	use mod_bound_geom

	implicit none

	integer ibc			!number of open boundary  (in)
	integer ndim			!dimension of nodes()     (in)
	integer idim			!total number of nodes    (out)
	integer nodes(ndim)		!boundary nodes           (out)

	integer i,imaxi
	integer kranf,krend

        call chkibc(ibc,'irbnds:')

        call get_bnd_ipar(ibc,'kranf',kranf)
        call get_bnd_ipar(ibc,'krend',krend)

	idim = krend - kranf + 1

	imaxi = min(ndim,idim)

	do i=1,imaxi
	  nodes(i) = irv(i+kranf-1)
	end do

	end

!********************************************************************

        subroutine ksinget(ibc,ampli,period,phase,zref)
        call chkibc(ibc,'ksinget:')
        call get_bnd_par(ibc,'ampli',ampli)
        call get_bnd_par(ibc,'period',period)
        call get_bnd_par(ibc,'phase',phase)
        call get_bnd_par(ibc,'zref',zref)
        end

!********************************************************************

	subroutine setbnd(ibc,value,barray)

! sets boundary ibc to value in barray (apparently not used)

	use mod_bound_geom

	implicit none

	integer ibc
	real value
	real barray(1)

	integer kranf,krend,k,kn

        call chkibc(ibc,'setbnd:')

        call get_bnd_ipar(ibc,'kranf',kranf)
        call get_bnd_ipar(ibc,'krend',krend)

	do k=kranf,krend
	  kn = irv(k)
	  barray(kn) = value
	end do

	end

!********************************************************************

        subroutine chkibc(ibc,errtext)
 
! checks if ibc is in bounds
 
	use mod_bnd

        implicit none
 
        integer ibc
	character*(*) errtext
 
        if( ibc .lt. 1 .or. ibc .gt. nbc ) then
	    write(6,*) errtext
            write(6,*) 'ibc, nbc : ',ibc,nbc
            stop 'error stop bnd: ibc out of bounds'
        end if

	end

!********************************************************************

	subroutine is_closed(ibc,nbc,icl)

	implicit none

	integer ibc,nbc,icl

	integer i
        integer nbnds,itybnd

        icl=0
        nbc = nbnds()

	if( ibc .gt. 0 ) then
	  if( itybnd(ibc) .lt. 0 ) icl = icl + 1
	else
          do i=1,nbc
            if( itybnd(i) .lt. 0 ) icl = icl + 1
          end do
	end if

	end

!********************************************************************

	subroutine get_oscil(ibc,dtime,zvalue)

	implicit none

	integer ibc
	double precision dtime
	real zvalue

        real pi
        parameter( pi=3.141592653 )

	real ampli,period,phase,zref

	call ksinget(ibc,ampli,period,phase,zref)

        zvalue = zref+ampli*sin(2.*pi*(dtime+phase)/period)

	end

!********************************************************************
!********************************************************************
!********************************************************************
!
! routines dealing with parameter values
!
!********************************************************************
!********************************************************************
!********************************************************************

	subroutine init_bnd_par(ibc)

! initializes boundary ibc

	use mod_bnd

	implicit none

	integer ibc

	bnd(:,ibc) = 0.

	end

!********************************************************************

        subroutine setget_bnd_par(ibc,ientry,name,value,bset)

! sets/gets value at entry ientry

	use mod_bnd

        implicit none

        integer ibc
        integer ientry
	character*(*) name
        real value
        logical bset

        call chkibc(ibc,'setget_bnd_par:')

        if( ientry < 1 .or. ientry > nbvdim ) then
          write(6,*) 'no such entry: ',ientry,ibc,name
	  stop 'error stop setget_bnd_par: no such entry'
        end if

        if( bset ) then
          bnd(ientry,ibc) = value
        else
          value = bnd(ientry,ibc)
        end if

        end

!********************************************************************
 
        subroutine set_bnd_par(ibc,name,value)
        character*(*) name
        id = iget_bnd_par_id(name,.true.)
        call setget_bnd_par(ibc,id,name,value,.true.)
        end

        subroutine get_bnd_par(ibc,name,value)
        character*(*) name
        id = iget_bnd_par_id(name,.true.)
        call setget_bnd_par(ibc,id,name,value,.false.)
        end

        subroutine set_bnd_ipar(ibc,name,ivalue)
        character*(*) name
        id = iget_bnd_par_id(name,.true.)
	value = ivalue
        call setget_bnd_par(ibc,id,name,value,.true.)
        end

        subroutine get_bnd_ipar(ibc,name,ivalue)
        character*(*) name
        id = iget_bnd_par_id(name,.true.)
        call setget_bnd_par(ibc,id,name,value,.false.)
	ivalue = nint(value)
        end

!********************************************************************

	subroutine copy_bnd_par(ibc)

! copies parameter values

	use mod_bnd

	implicit none

	integer ibc

	integer id
	real value
	character*6 name

	real getpar

	do id=1,nbvdim
	  call get_bnd_par_name(id,name)
	  value = getpar(name)
          call set_bnd_par(ibc,name,value)
	end do

	end

!********************************************************************

	subroutine get_bnd_npar(nbnd)

	use mod_bnd

        implicit none

	integer nbnd

	nbnd = nbvdim

	end

!********************************************************************

        function exists_bnd_par(name)

! tests if parameter name exists

        implicit none

        logical exists_bnd_par
        character*(*) name

	integer iget_bnd_par_id

	exists_bnd_par = iget_bnd_par_id(name,.false.) > 0

	end

!********************************************************************

        function iget_bnd_par_id(name,berror)

! gets id given a parameter name

	use mod_bnd

        implicit none

        integer iget_bnd_par_id
        character*(*) name
	logical berror		!raises error if name not existing

        integer id
        character*6 bname

        do id=1,nbvdim
	  call get_bnd_par_name(id,bname)
          if( bname .eq. name ) then
            iget_bnd_par_id = id
            return
          end if
        end do

	if( berror ) then
          write(6,*) 'unknown parameter name for boundary: ',name
          stop 'error stop iget_bnd_par_id: name'
	end if

	iget_bnd_par_id = 0

        end

!********************************************************************

	subroutine get_bnd_par_name(id,name)

! gets parameter name given id

	use mod_bnd

	implicit none

	integer id
        character*(*) name

	if( id < 1 .or. id > nbvdim ) then
	  name = ' '
	else
	  name = bnd_par_names(id)
	end if

	end

!********************************************************************

	subroutine check_bnd_par_entries(ibc,ball)

! writes parameter values for given boundary

	use mod_bnd

	implicit none

	integer ibc
	logical ball	!write all parameter names

        integer i
	real value
	character*6 name

	!write(6,*) 'check_bnd_par_entries: ',ibc
        do i=1,nbvdim
	  call get_bnd_par_name(i,name)
          call setget_bnd_par(ibc,i,' ',value,.false.)
	  if( ball .or. value /= 0. ) then
	    write(6,*) name,value
	  end if
	end do

        end

!********************************************************************
!********************************************************************
!********************************************************************
!
! routines dealing with file names
!
!********************************************************************
!********************************************************************
!********************************************************************

	subroutine init_bnd_file(ibc)

! initializes boundary ibc

	use mod_bnd

	implicit none

	integer ibc

	bnd_file(:,ibc) = ' '

	end

!********************************************************************

        subroutine setget_bnd_file(ibc,ientry,name,file,bset)

! sets/gets file at entry ientry

	use mod_bnd

        implicit none

        integer ibc
        integer ientry
	character*(*) name
	character*(*) file
        logical bset

        call chkibc(ibc,'setget_bnd_file:')

        if( ientry < 1 .or. ientry > nbfdim ) then
          write(6,*) 'no such entry: ',ientry,ibc,name
	  stop 'error stop setget_bnd_file: no such entry'
        end if

        if( bset ) then
          bnd_file(ientry,ibc) = file
        else
          file = bnd_file(ientry,ibc)
        end if

        end

!********************************************************************
 
        subroutine set_bnd_file(ibc,name,file)
        character*(*) name,file
        id = iget_bnd_file_id(name,.true.)
        call setget_bnd_file(ibc,id,name,file,.true.)
        end

        subroutine get_bnd_file(ibc,name,file)
        character*(*) name,file
        id = iget_bnd_file_id(name,.true.)
        call setget_bnd_file(ibc,id,name,file,.false.)
        end

!********************************************************************

	subroutine copy_bnd_file(ibc)

! copies file names

	use mod_bnd

	implicit none

	integer ibc

	integer id
	character*80 file
	character*6 name

	do id=1,nbfdim
	  call get_bnd_file_name(id,name)
	  call getfnm(name,file)
          call set_bnd_file(ibc,name,file)
	end do

	end

!********************************************************************

	subroutine get_bnd_nfile(nbnd)

	use mod_bnd

        implicit none

	integer nbnd

	nbnd = nbfdim

	end

!********************************************************************

        function exists_bnd_file(name)

! tests if file name exists

        implicit none

        logical exists_bnd_file
        character*(*) name

	integer iget_bnd_file_id

	exists_bnd_file = iget_bnd_file_id(name,.false.) > 0

	end

!********************************************************************

        function iget_bnd_file_id(name,berror)

! gets id given a file name

	use mod_bnd

        implicit none

        integer iget_bnd_file_id
        character*(*) name
	logical berror		!raises error if name not existing

        integer id
        character*6 bname

        do id=1,nbfdim
	  call get_bnd_file_name(id,bname)
          if( bname .eq. name ) then
            iget_bnd_file_id = id
            return
          end if
        end do

	if( berror ) then
          write(6,*) 'unknown file name for boundary: ',name
          stop 'error stop iget_bnd_file_id: name'
	end if

	iget_bnd_file_id = 0

        end

!********************************************************************

	subroutine get_bnd_file_name(id,name)

! gets file name given id

	use mod_bnd

	implicit none

	integer id
        character*(*) name

	if( id < 1 .or. id > nbfdim ) then
	  name = ' '
	else
	  name = bnd_file_names(id)
	end if

	end

!********************************************************************

	subroutine check_bnd_file_entries(ibc,ball)

! writes file names for given boundary

	use mod_bnd

	implicit none

	integer ibc
	logical ball	!write all file names

        integer i
	character*6 name
	character*80 file

	!write(6,*) 'check_bnd_file_entries: ',ibc
        do i=1,nbfdim
	  call get_bnd_file_name(i,name)
          call setget_bnd_file(ibc,i,' ',file,.false.)
	  if( ball .or. file /= ' ' ) then
	    write(6,*) name,'  ',trim(file)
	  end if
	end do

        end

!********************************************************************
!********************************************************************
!********************************************************************

	subroutine get_boundary_file(ibc,what,file)

	implicit none

	integer ibc
	character*(*) what
	character*(*) file

	character*6 name

	logical exists_bnd_file

        if( what .eq. 'zeta' ) then
	  name = 'boundn'
        else if( what .eq. 'conz' ) then
	  name = 'conzn'
        else if( what .eq. 'temp' ) then
	  name = 'tempn'
        else if( what .eq. 'salt' ) then
	  name = 'saltn'
        else if( what .eq. 'vel' ) then
	  name = 'vel3dn'
        else if( what .eq. 'sedt' ) then
	  name = 'sed2dn'
        else if( what .eq. 'lagvebio' ) then
	  name = 'bio2dn'
        else if( what .eq. 'toxi' ) then
	  name = 'tox3dn'
        else if( what .eq. 'bfm1' ) then
	  name = 'bfm1bc'
        else if( what .eq. 'bfm2' ) then
	  name = 'bfm2bc'
        else if( what .eq. 'bfm3' ) then
	  name = 'bfm3bc'
        else if( what .eq. 'bfm' ) then
	  name = 'bfmbcn'
        else if( what .eq. 'mercury' ) then
	  name = 'mercn'
        else if( what .eq. 's4mercury' ) then
          name = 's4mern'
        else
          if( exists_bnd_file(what) ) then	!use name given in what
	    name = what
	  else
            write(6,*) 'keyword not recognized: ',what
            write(6,*) 'boundary: ',ibc
            stop 'error stop get_boundary_file'
	  end if
        end if

        call get_bnd_file(ibc,name,file)

	end

!********************************************************************
!********************************************************************
!********************************************************************

