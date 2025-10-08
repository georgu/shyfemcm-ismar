
!--------------------------------------------------------------------------
!
!    Copyright (C) 2002-2012,2014-2016,2018-2019  Georg Umgiesser
!    Copyright (C) 2014-2016  Christian Ferrarin
!    Copyright (C) 2016  Ivan Federico
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

! heat flux module (file administration)
!
! contents :
!
! subroutine qflux_init
!		initializes routines
! subroutine qflux3d(it,dt,nkn,nlvdi,temp,dq)
!		computes new temperature (forced by heat flux) - 3d version
!
! revision log :
!
! 01.02.2002	ggu	new as administrative routines
! 11.04.2003	ggu	albedo introduced in heat2t
! 16.08.2004	ggu	some comments on usage, heat2t copied to subqfxu4.f
! 22.02.2005	ggu	subroutine qflux2d deleted
! 23.03.2006	ggu	changed time step to real
! 20.08.2007	ggu	prepared for evaporation - salinity feedback
! 17.04.2008	ggu	save evaporation in evapv for later use
! 17.11.2008	ggu	new call to qfcheck_file() before opening
! 10.03.2009	ggu	new qflux_read(), new call to meteo_get_values()
! 11.03.2009	ggu	new routine qflux_compute()
! 27.08.2009	ggu	new call to heatgill, new routine heatareg
! 11.11.2009	ggu	handle abosrption of short rad in more than one layer
! 23.03.2010	ggu	changed v6.1.1
! 04.03.2011	ggu	new routine heatgotm
! 23.03.2011	ggu	new routine check_heat() to check for Nan, new iheat
! 25.03.2011	ggu	new parameters iheat,hdecay,botabs
! 14.04.2011	ggu	changed VERS_6_1_22
! 14.02.2012	ggu	changed VERS_6_1_44
! 30.03.2012	ggu	changed VERS_6_1_51
! 29.04.2014	ccf	read qsens, qlat, long from file
! 16.06.2014	ccf	new routine heatcoare, which also update wind stress
! 20.06.2014	ccf	new routine for computing sea surface skin temperature
! 18.07.2014	ggu	changed VERS_7_0_1
! 05.11.2014	ggu	changed VERS_7_0_5
! 26.11.2014	ggu	changed VERS_7_0_7
! 05.12.2014	ggu	changed VERS_7_0_8
! 19.01.2015	ggu	changed VERS_7_1_3
! 26.02.2015	ggu	changed VERS_7_1_5
! 30.04.2015	ggu	changed VERS_7_1_9
! 05.06.2015	ggu	changed VERS_7_1_12
! 10.07.2015	ggu	changed VERS_7_1_50
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 18.09.2015	ccf	do not compute heat fluxes in dry nodes
! 18.09.2015	ccf	checks only for levdbg > 2
! 26.10.2015	ggu	critical omp sections introduced (eliminated data race)
! 18.12.2015	ggu	changed VERS_7_3_17
! 07.04.2016	ggu	compute total evaporation
! 15.04.2016	ggu	changed VERS_7_5_8
! 04.05.2016	ccf	do not pass albedo into heat2t
! 04.05.2016	ggu	include effect of ice cover
! 25.05.2016	ggu	changed VERS_7_5_10
! 21.07.2016	ivn	isolp = 1, 2 length scale solar penetration (Jerlov) 
! 09.09.2016	ggu	changed VERS_7_5_17
! 29.09.2016	ivn	bug fix for isolp = 1
! 24.01.2018	ggu	changed VERS_7_5_41
! 03.04.2018	ggu	changed VERS_7_5_43
! 16.10.2018	ggu	changed VERS_7_5_50
! 14.02.2019	ggu	changed VERS_7_5_56
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 10.12.2019	ggu	ice cover introduced, handle ice-water sensible heat
! 11.11.2020	ggu	new ice model integrated, cleaned, documented
! 13.11.2020	ggu	bug fix: correct rain from [m/s] to [mm/d]
! 14.11.2020	ggu	allow for ice and other heat fluxes to coexist
! 29.04.2022	ggu	check impossible heat flux values
! 09.07.2022	ggu	kspecial for special output
! 05.09.2022	ggu	debug output and use ustar in ice computation
! 09.05.2023    lrp     introduce top layer index variable
! 06.09.2024    lrp     nuopc-compliant
! 13.09.2024    lrp     iatm and coupling with atmospheric model
! 21.09.2024    ggu     new handle_skin() for handling skin temperature output
! 08.02.2025    ggu     some new debug code
!
! notes :
!
! qflux_init is called in subn11.f
! qflux_read is called in subn11.f
! qflux3d is called in newbcl.f
! qfget is called in qflux3d (here)
!
! to use heat flux module insert in section "name" of STR file:
!
! $name
!     qflux = 'qflux.dat'
! $end
!
! if such a file is given the heat flux module is used and works
!   without any other intervention (ibarcl>0)
!
! other notes in subqfxf.f
!
!*****************************************************************************

	subroutine qflux_compute(byes)

! returnes flag -> compute heat flux or not

	implicit none

	logical byes

	character*80 file

        call getfnm('qflux',file)

	byes = ( file /= ' ' )

	end

!*****************************************************************************

	subroutine qflux_init

! initializes routines

	implicit none

	character*80 file

        call getfnm('qflux',file)

	if( file .ne. ' ' ) then
	  call qfcheck_file(file)
	  call qfinit(file)
	end if

	end

!*****************************************************************************

	subroutine qflux_read(it)

! reads new meteo data

	implicit none

	integer it

	call qfmake(it)

	end

!*****************************************************************************

	subroutine qflux3d(dtime,dt,nkn,nlvddi,temp,dq)

! computes new temperature (forced by heat flux) - 3d version

	use mod_meteo
	use mod_ts
	use levels
        use basin, only : xgv,ygv,ipv
        use mod_hydro_print  
	use heat_const
	use meteo_forcing_module, only : iatm

	implicit none

	double precision dtime
	real dt
	integer nkn
	integer nlvddi
        real temp(nlvddi,nkn) 
	double precision dq	!total energy introduced [(W/m**2)*dt*area = J]

        logical byes
        logical buseice,bicecover
        integer levdbg,iud,ksd,ksd1,iudbg
	integer k,kext,kssalt
	integer l,lmax,lmin,kspec,ilist
	integer mode
	integer days,im,ih
	integer ys(8)
	real hm,sm,tm,tnew
	real salt,tfreeze
	real albedo
	real adecay,qsbottom
	real qtot
        real qs,ta,tb,uw,cc,ur,p,e,r,q
        real ddlon,ddlat  
        real dp,uuw,vvw  
	real qsens,qlat,qlong,evap,qrad,qsdown,qsens_orig
        real qswa,qice,qsurface,qsolar
	real cice,fice_free,fice_cover
	real ev,eeff
	real area
	real evaver
        real uub,vvb        
	real ik1,ik2
	real tice,u,v,uv

	real, parameter :: epsice = 1.e-5

	double precision ddq
	double precision atime
	character*20 aline

	real, save, allocatable :: dtw(:)	!Warm layer temp. diff
	real, save, allocatable :: tws(:)	!Skin temperature (deg C)

	real hb			!depth of modelled T in first layer [m]
	real usw		!surface friction velocity [m/s]
	real qss		!Net shortwave flux 
	real cd			!wind drag coefficient

!     coefficient after Paul e Simon (1977) up to type IV
!     for coastal water (1-3-5-7-9) coefficient fitting data with Jerlov(1968)
      real,parameter, dimension(10) :: Rt =              &
           &          (/ 0.58, 0.62, 0.67, 0.77, 0.78,   &
           &         0.605, 0.61, 0.37, 0.338, 0.277 /)
      real,parameter, dimension(10) :: k1 =              &
           &         (/ 0.35, 0.6,  1.0,  1.5,  1.4,     &
           &         0.387, 0.353, 2.395, 1.8, 1.505 /)
      real,parameter, dimension(10) :: k2 =              &
           &       (/ 23.0, 20.0, 17.0, 14.0, 7.9,       &
           &          5.05, 3.55, 0.335, 0.329, 0.325 /)

	integer itdrag

! functions
	real depnode,areanode,getpar
	integer ifemopa,ipext
	logical is_dry_node
! save
	integer, save :: icall = 0
	integer, save :: n93 = 0

	integer, save :: iheat
	integer, save :: isolp  
	integer, save :: iwtyp
	real, save :: hdecay
	real, save :: botabs

	real, save :: aice = 1.
	real, save :: fice_pen = 0.3		!penetration through ice

	integer, save :: kdebug = 0		!debug for this node
	integer, save :: kspecial = 0		!special output for this node
	logical, save :: bdebug = .false.
	logical, save :: bfluxdebug = .true.
	logical, save :: baverevap = .false.
	logical, save :: bwind = .false.
	logical, save :: bice = .false.
	logical, save :: bheat = .false.
	logical, save :: batm = .false.
	logical, save :: bqflux = .false.

!---------------------------------------------------------
! start of routine
!---------------------------------------------------------

	dq = 0.
	if( icall < 0 ) return

!---------------------------------------------------------
! iheat		1=areg  2=pom  3=gill  4=dejak  5=gotm  
!		6=COARE 3.0
!               7=read flux from file
!		8=MFS routines (Pettenuzzo et al., 2010)
! hdecay	depth of e-folding decay of radiation
!		0. ->	everything is absorbed in first layer
! botabs	1. ->	bottom absorbs remaining radiation
!		0. ->	everything is absorbed in last layer
!---------------------------------------------------------
! format of heat file containing time series (4 data columns):
!    time srad airt rhum cc
! in case of iheat==7 the columns are:
!    time srad qsens qlat qlong
!---------------------------------------------------------

!---------------------------------------------------------
! initialize on first call
!---------------------------------------------------------

	if( icall .eq. 0 ) then
	  iheat = nint(getpar('iheat'))
          isolp = nint(getpar('isolp'))   
	  iwtyp  = nint(getpar('iwtyp'))+1
	  hdecay = getpar('hdecay')
	  botabs = getpar('botabs')

	  call set_iheat(iheat)

	  call shyice_init
	  call shyice_is_active(bice)		!ice model is active
	  bheat = iheat > 0			!we have to compute heat fluxes
	  batm = iatm > 0			!have ocean-atmosphere coupling
	  call qflux_compute(bqflux)		!has qflux file?

	  if( .not. bqflux .and. .not. batm ) icall = -1
	  if( .not. bheat .and. .not. bice ) icall = -1

	  !---------------------------------------------------------
	  ! exit if no files, no atm-oce coupling, no heat fluxes, no ice model 
	  !---------------------------------------------------------

	  if( icall < 0 ) return

!$OMP CRITICAL
	  write(6,*) 'qflux3d routines are active'
	  write(6,*) 'qflux3d: bqflux,bheat,bice: ',bqflux,bheat,bice
	  write(6,*) 'qflux3d: iheat,hdecay,botabs: ',iheat,hdecay,botabs  
!$OMP END CRITICAL

	  allocate(dtw(nkn))
	  allocate(tws(nkn))
	  dtw = 0.
	  do k=1,nkn
	    lmin = jlhkv(k)
	    tws(k) = temp(lmin,k)
	  end do

          itdrag = nint(getpar('itdrag'))
	  bwind = itdrag .eq. 4
          if ( bwind ) then   
            if ( iheat .eq. 6 .or. iheat .eq. 8 ) then   
              write(6,*) 'itdrag = ',itdrag
              write(6,*) 'iheat = ',iheat
            else
              write(6,*) 'Erroneous value for itdrag = ',itdrag
              write(6,*) 'Use itdrag = 4 only with iheat = 6 or 8'
              stop 'error stop qflux3d: itdrag'
            end if
          end if

	  call shyice_init_output
	  call meteo_get_ice_usage(aice)	!use ice (1) or not (0)
	  if( bice ) fice_pen = 0.

          levdbg = nint(getpar('levdbg'))
	  bdebug = levdbg .ge. 2 
	end if

	icall = icall + 1

!---------------------------------------------------------
! set other parameters
!---------------------------------------------------------

	mode = +1	!use new time step for depth

	adecay = 0.
	if( hdecay .gt. 0. ) adecay = 1. / hdecay

	ik1 = 0.
	ik2 = 0.
	if ( isolp .eq. 1 ) then
          if (iwtyp > 10 ) then
            write(6,*) 'not valid value for iwtyp = ',iwtyp-1
            write(6,*) 'use values from 0 (Jrv t-I) to 10 (Jrv t-9)'
            stop 'error stop qflux3d: iwtyp'
          end if
          ik1 = 1. / k1(iwtyp)
          ik2 = 1. / k2(iwtyp)
	end if

!---------------------------------------------------------
! set date parameters
!---------------------------------------------------------

	call get_act_timeline(aline)
	call get_absolute_act_time(atime)
	call dts_from_abs_time_to_ys(atime,ys)
	call dts_from_abs_time_to_days_in_year(atime,days)
	days = days - 1
	im = ys(2)
	ih = ys(4)

!---------------------------------------------------------
! loop over nodes
!---------------------------------------------------------

        ddq = 0.
	ilist = 0
	kssalt = 0

	do k=1,nkn

	  kext = ipv(k)
	  lmax = ilhkv(k)
	  lmin = jlhkv(k)

	  if( kssalt > 0 .and. saltv(1,k) > 50. ) then
	    if( ilist == 0 ) write(167,*) aline
	    ilist = ilist + 1
	    write(167,*) ilist,k,kext,is_dry_node(k),saltv(1,k)
	  end if
	  if( kext == kssalt .and. saltv(1,k) > 50. ) then
	    iudbg = 170
	    call check_set_unit(iudbg)
	    call check_node(k)
	    call check_elems_around_node(k)
	  end if

	  if (is_dry_node(k)) then	!do not compute if node is dry
	    dtw(k)   = 0.
	    tws(k)   = temp(lmin,k)
	    evapv(k) = 0.
	    cycle
	  end if

	  tm = temp(lmin,k)
	  salt = saltv(lmin,k)
	  call get_ice_cover(k,cice)
	  !if( cice < 0. .or. cice > 1. ) goto 98
	  if( cice < -epsice .or. cice > 1.+epsice ) goto 98
	  fice_cover = aice*cice
	  fice_free  = 1. - fice_cover
	  bicecover = fice_cover > 0.
	  buseice = bice .and. bicecover	!we use ice model
	  tfreeze = -0.0575*salt
	  qsens = 0.
	  qlat = 0.
	  qlong = 0.
	  evap = 0.
	  area = areanode(lmin,k)
          if( isolp .eq. 0 .and. hdecay .le. 0. ) lmax = 1   

	  !------------------------------------------------
	  ! solar radiation - positive from air to sea
	  !------------------------------------------------

          call meteo_get_heat_values(k,qs,ta,ur,tb,uw,cc,p)
	  call make_albedo(tm,albedo)
	  qsdown = qs * (1. - albedo)
	  qss = fice_free*qsdown + fice_cover*fice_pen*qsdown
	  if( .not. bheat ) qss = 0.

	  !------------------------------------------------
	  ! heat exchange - qlong, qlat, qsens positive from sea to air
	  !------------------------------------------------

	  if( iheat .eq. 0 ) then
	    buseice = bice			!we use ice model heat fluxes
	  else if( iheat .eq. 1 ) then
	    call heatareg (ta,p,uw,ur,cc,tm,qsens,qlat,qlong,evap)
	  else if( iheat .eq. 2 ) then
	    call heatpom  (ta,p,uw,ur,cc,tm,qsens,qlat,qlong,evap)
	  else if( iheat .eq. 3 ) then
	    call heatgill (ta,p,uw,ur,cc,tm,qsens,qlat,qlong,evap)
	  else if( iheat .eq. 4 ) then
	    !call rh2wb(ta,p,ur,tb)	!FIXME
	    call heatlucia(ta,p,uw,tb,cc,tm,qsens,qlat,qlong,evap)
	  else if( iheat .eq. 5 ) then
	    call heatgotm (ta,p,uw,ur,cc,tm,qsens,qlat,qlong,evap)
	  else if( iheat .eq. 6 ) then
	    call get_pe_values(k,r,ev,eeff)
	    call heatcoare(ta,p,uw,ur,cc,tws(k),r,qss,qsens,qlat,qlong,evap,cd)
	    if ( bwind ) windcd(k) = cd
	  else if( iheat .eq. 7 ) then
	    qsens = ta
	    qlat  = ur
	    !qlong = -cc  !change sign of long wave radiation given by ISAC
	    qlong = abs(cc)
	    evap  = qlat / (2.5008e6 - 2.3e3 * tm)	!pom, gill, gotm
          else if( iheat .eq. 8 ) then
            ddlon = xgv(k)    
            ddlat = ygv(k)   
            uub = uprv(lmin,k)  
            vvb = vprv(lmin,k)  
	    call meteo_get_heat_extra(k,dp,uuw,vvw)
       call heatmfsbulk(days,im,ih,ddlon,ddlat,ta,p,uuw,vvw,dp,          &
           &                   cc,tm,uub,vvb,qsens,qlat,qlong,evap,qswa,cd)   
            qss= fice_free*qswa  !albedo (monthly) already in qshort1 -> qswa  
            if ( bwind ) windcd(k) = cd   
          else
            write(6,*) 'iheat = ',iheat
            stop 'error stop qflux3d: value for iheat not allowed'
          end if

	  if (bdebug) call check_heat(k,tm,qsens,qlat,qlong,evap)

	  !------------------------------------------------
	  ! heat transfer from water to ice - qice positive from water to ice
	  !------------------------------------------------

	  qice = 0.
	  if( bicecover .and. .not. bice ) then	!ice on node but no ice model
	    tice = 0.
	    call getuv(lmin,k,u,v)
	    uv = sqrt(u*u+v*v)		!current speed
	    call ice_water_exchange(tm,tice,uv,qice)
	    if( abs(qice) > 100000. ) goto 99
	    if( .not. bheat ) qice = 0.
	  else
	    tice = 0.
	    uv = 0.
	  end if

	  !------------------------------------------------
	  ! final heat exchange 
	  ! qlong, qlat, qsens are positive from sea to air
	  ! qss (solar radiation) is positive from air to sea
	  ! qrad, qss, qtot are positive from air to sea
	  !------------------------------------------------

	  qlong = fice_free * qlong
	  qlat = fice_free * qlat
	  qsens_orig = qsens
	  qsens = fice_free * qsens + fice_cover * qice
          qrad =  - ( qlong + qlat + qsens )
	  qtot = qss + qrad
	  qsolar = qss

          !------------------------------------------------
          ! debug output
          !------------------------------------------------

	  ksd1 = 877
	  ksd = 1461
	  ksd1 = 0
	  ksd = 0
          if( kext == ksd .or. kext == ksd1 ) then
            iud = 667
            write(iud,*) '-------------------------'
            write(iud,*) atime
            write(iud,*) aline
            write(iud,*) k,kext
            write(iud,*) fice_cover,qice,aice,cice
            write(iud,*) fice_free,fice_pen
            write(iud,*) qlong,qlat,qsens,qss
            write(iud,*) tm,tice,uv,qice
            write(iud,*) tm
            write(iud,*) '-------------------------'
          end if

	  !------------------------------------------------
	  ! distribute short wave radiation to layers and compute heat budget
	  !------------------------------------------------

	  qsurface = qrad

          do l=lmin,lmax
	    if( .not. bheat ) cycle	!no heat flux computed
	    if( buseice ) cycle		!we use heat flux from ice model
            hm = depnode(l,k,mode)
            tm = temp(l,k)
            if (isolp == 0) then  
              if( hdecay .le. 0. ) then
                qsbottom = 0.
              else
                qsbottom = qss * exp( -adecay*hm )
                if( l .eq. lmax ) qsbottom = botabs * qsbottom
              end if
            else if (isolp == 1) then
              !Jerlov classification (1968) from Type I to 9
               qsbottom = qss * (Rt(iwtyp) * exp(- ik1 * hm )                 &
           &                            +  (1-Rt(iwtyp)) * exp(- ik2 * hm))
              if( l .eq. lmax ) qsbottom = botabs * qsbottom
            else
              write(6,*) 'Erroneous value for isolp = ',isolp
              write(6,*) 'Use isolp = 0 (hdecay) or 1 (Jerlov t-I)'
              stop 'error stop qflux3d: isolp'
            end if
	    if( abs(qsurface) > 100000. ) goto 99
            call heat2t(dt,hm,qss-qsbottom,qsurface,tm,tnew)
            if (bdebug) call check_heat2(k,l,qss,qsbottom,qsurface,     &
           &                                   albedo,tm,tnew)
            tnew = max(tnew,tfreeze)
            temp(l,k) = tnew
            albedo = 0.
            qsurface = 0.	!different from 0 only in first layer
            qss = qsbottom
          end do

	  if( k == kspecial ) then
	    write(177,*) aline,qsolar,qrad,tnew
	  end if

!         ---------------------------------------------------------
!         compute sea surface skin temperature
!         ---------------------------------------------------------

	  tm   = temp(lmin,k)
	  hb   = depnode(lmin,k,mode) * 0.5
          usw  = max(1.e-5, sqrt(sqrt(tauxnv(k)**2 + tauynv(k)**2)))
	  call tw_skin(qss,qrad,tm,hb,usw,dt,dtw(k),tws(k))

!         ---------------------------------------------------------
!         handle ice model
!         ---------------------------------------------------------

	  if( k == kdebug ) write(444,*) icall,bice,bicecover,buseice
	  if( bice ) then
            hm = depnode(lmin,k,mode)
            tm = temp(lmin,k)
            sm = saltv(lmin,k)
	    call get_pe_values(k,r,e,eeff)
	    r = r * 1000 * 86400	!r is in [m/s] -> convert to [mm/d]
	    call shyice_run(k,qs,ta,ur,uw,cc,p,r,hm,tm,sm,dt)
	    if( buseice ) then	!ice on water or bheat == .false.
              temp(lmin,k) = tm
              saltv(lmin,k) = sm
	      evap = 0.
	      dtw(k) = 0.
	      call shyice_get_tsurf(k,tws(k))
	    end if
	  end if

!	  ---------------------------------------------------------
!	  evap is in [kg/(m**2 s)] -> convert it to [m/s]
!	  evap is normally positive -> we are loosing mass
!	  if batm is true -> evapv is diagnostic (already in rain)
!	  ---------------------------------------------------------

	  evap = evap / rhow			!evaporation in m/s
	  evapv(k) = evap			!positive if loosing mass
	  qsensv(k) = qsens
	  qlatv(k) = qlat
	  qlongv(k) = qlong

          ddq = ddq + qtot * dt * area

	end do

	if( bfluxdebug ) call debug_fluxes(aline)

	dq = ddq

!---------------------------------------------------------
! compute total evaporation
!---------------------------------------------------------

	if( baverevap ) then
	  call aver_nodal(evapv,evaver)	!in evaver is average of evaporation m/s
	  write(678,*) dtime,evaver
	end if


!---------------------------------------------------------
! skin temperature and special output
!---------------------------------------------------------

	call handle_skin(dtime,tws)

	call shyice_write_output

	k = min(nkn,1000)
	k = -1
	if( k .gt. 0 ) then
	  if( n93 .eq. 0 ) then
	    n93 = ifemopa('opening file 93','.93','form','new')
	  end if
	  write(n93,*) 'qflux3d: ',dtime,temp(1,k)
	end if

!---------------------------------------------------------
! end of routine
!---------------------------------------------------------

	return
   98	continue
	write(6,*) 'impossible value for cice:'
	write(6,*) k,ipv(k),cice
	stop 'error stop qflux3d: impossible value for cice'
   99	continue
	write(6,*) 'error in heat flux'
	write(6,*) 'k',k,ipext(k),iheat
	write(6,*) 'l',l,lmin,lmax,hm
	write(6,*) 'qurf',qsurface,qrad
	write(6,*) 'qlong',qlong,qlat,qsens,qss
	write(6,*) 'qice',qsens_orig,qice
	write(6,*) 'fice',fice_free,fice_cover
	write(6,*) 'ta',ta,p,uw,ur,cc
	write(6,*) 'tm',tm,tws(k),r
	flush(6)
	stop 'error stop qflux3d: impossible heat flux'
	end

!*****************************************************************************

	subroutine ice_water_exchange(tw,ti,uv,qice)

! Li at al.,
! Heat transfer at ice-water interface under conditions of low flow velocities
! Journal of Hydrodynamics,2016,28(4):603-609
! DOI: 10.1016/S1001-6058(16)60664-9

! use ustar instead of uv in formula
! The role of the heat exchange coefficient at the ice/ocean interface
! in Bohai Sea ice modeling
! Bin Jia et al., Continental Shelf Research 241 (2022) 104735
! (attention: formula 8 is wrong, needs velocity difference squared)

! qice is heat flux [W/m**2] from water to ice (positive)

	implicit none

	real tw			!temperature of water
	real ti			!temperature of ice
	real uv			!current speed
	real qice		!ice heat flux (return)

	real ustar

	real, parameter :: rhow = 1000.	!density of water
	real, parameter :: cw = 4179.6	!specific heat of water
	real, parameter :: ch = 1.1E-3	!heat transfer coefficient
	real, parameter :: cdw = 0.0055	!ice ocean drag coefficient

	ustar = sqrt(cdw)*uv
	qice = rhow*cw*ch*ustar*(tw-ti)

	if( abs(qice) > 100000. ) then
	  write(6,*) tw,ti,uv,qice
	  write(6,*) rhow,cw,ch,cdw
	  write(6,*) ustar,tw-ti,sqrt(cdw)
	  stop 'error stop ice_water_exchange: qice'
	end if

	end

!*****************************************************************************

	subroutine check_heat(k,tm,qsens,qlat,qlong,evap)

	use mod_debug

	implicit none

	integer k
	real tm,qsens,qlat,qlong,evap

	if( is_nan(tm) ) goto 98
	if( is_nan(qsens) ) goto 98
	if( is_nan(qlat) ) goto 98
	if( is_nan(qlong) ) goto 98
	if( is_nan(evap) ) goto 98

	return
   98	continue
	write(6,*) k,tm,qsens,qlat,qlong,evap
	stop 'error stop check_heat: Nan found'
	end

!*****************************************************************************

	subroutine check_heat2(k,l,qs,qsbottom,qrad,albedo,tm,tnew)

	use mod_debug

	implicit none

	integer k,l
	real qs,qsbottom,qrad,albedo,tm,tnew

	if( is_nan(qs) ) goto 98
	if( is_nan(qsbottom) ) goto 98
	if( is_nan(qrad) ) goto 98
	if( is_nan(albedo) ) goto 98
	if( is_nan(tm) ) goto 98
	if( is_nan(tnew) ) goto 98

	if( abs(tm) .gt. 100. ) goto 97
	if( abs(tnew) .gt. 100. ) goto 97

	return
   97	continue
	write(6,*) k,l,tm,tnew
	stop 'error stop check_heat2: t out of range'
   98	continue
	write(6,*) k,l,qs,qsbottom,qrad,albedo,tm,tnew
	stop 'error stop check_heat2: Nan found'
	end

!*****************************************************************************

	subroutine handle_skin(dtime,tws)

! handles output of skin temperature

	use basin

	implicit none

	double precision dtime
	real tws(nkn)

	integer, parameter :: nvar = 1
	logical, parameter :: b2d = .true.
	double precision, save :: da_out(4)
	integer, save :: icall = 0
	integer id,iskin

	logical has_output_d,next_output_d
	real getpar

	if( icall < 0 ) return

	if( icall == 0 ) then
	  icall = -1
	  iskin = nint(getpar('iskin'))
	  if( iskin == 0 ) return			!no skin output wanted
          call init_output_d('itmcon','idtcon',da_out)
          if( .not. has_output_d(da_out) ) return	!has no skin output
          call shyfem_init_scalar_file('tskin',nvar,b2d,id)
          da_out(4) = id
	  icall = 1
	end if

        if( next_output_d(da_out) ) then
          id = nint(da_out(4))
          call shy_write_scalar_record(id,dtime,14,1,tws)
          call shy_sync(id)
        end if

	end

!*****************************************************************************

	subroutine debug_fluxes(aline)

	use mod_meteo

	implicit none

	character*20 aline

	integer, save :: icall = 0
	integer, save :: k1,k2
	integer k

	integer ipint

	return

	if( icall == 0 ) then
          k = 10337
          k = ipint(k)
          if( k <= 0 ) stop 'error k1'
          k1 = k

	  k = 8636
          k = ipint(k)
          if( k <= 0 ) stop 'error k2'
          k2 = k
        end if

	icall = icall + 1

	k = k1
	write(661,1000) aline,evapv(k),qsensv(k),qlatv(k),qlongv(k)
	k = k2
	write(662,1000) aline,evapv(k),qsensv(k),qlatv(k),qlongv(k)
	flush(661)
	flush(662)

	return
 1000	format(a20,6f10.2)
	end

!*****************************************************************************

