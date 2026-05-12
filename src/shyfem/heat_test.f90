
!**************************************************************

	program heat_test

! tests heat modules for one single node

	implicit none

	integer i,k,ierr
	integer iheat
	integer imax
	real ta,p,uw,ur,cc,tm,qsens,qlat,qlong,evap,tb,tws,zeta
	real cd,r,ev,eeff
	real fice_free,fice_cover,fice_pen
	real hm,hm_base,hm_min
	real qice
	real qss,qrad,qsens_orig,qsolar,qtot
	real albedo,qs,qsdown
	real qsbottom,qsurface
	real dt,tnew
	double precision dtime

	iheat = 1
	tm = 10.
	hm_base = 0.10
	hm_min = 0.02
	dt = 300
	imax = 1000
	dtime = 0

        fice_cover = 0
        fice_free  = 1. - fice_cover
	fice_pen = 0.3
	qice = 0.
	r = 0			! rain

	do

	  dtime = dtime + dt

          call get_heat_values(dtime,qs,ta,ur,tb,uw,cc,p,zeta,ierr)
	  if( ierr < 0 ) exit
	  if( ierr > 0 ) goto 99

          !call make_albedo(tm,albedo)
	  albedo = 0.06
          qsdown = qs * (1. - albedo)
          qss = fice_free*qsdown + fice_cover*fice_pen*qsdown

          if( iheat .eq. 0 ) then
            !buseice = bice                      !we use ice model heat fluxes
          else if( iheat .eq. 1 ) then
            call heatareg (ta,p,uw,ur,cc,tm,qsens,qlat,qlong,evap)
          else if( iheat .eq. 2 ) then
            call heatpom  (ta,p,uw,ur,cc,tm,qsens,qlat,qlong,evap)
          else if( iheat .eq. 3 ) then
            call heatgill (ta,p,uw,ur,cc,tm,qsens,qlat,qlong,evap)
          else if( iheat .eq. 4 ) then
            !call rh2wb(ta,p,ur,tb)     !FIXME
            call heatlucia(ta,p,uw,tb,cc,tm,qsens,qlat,qlong,evap)
          else if( iheat .eq. 5 ) then
            call heatgotm (ta,p,uw,ur,cc,tm,qsens,qlat,qlong,evap)
          else if( iheat .eq. 6 ) then
            !call get_pe_values(k,r,ev,eeff)
            call heatcoare(ta,p,uw,ur,cc,tws(k),r,qss,qsens,qlat,qlong,evap,cd)
           ! if ( bwind ) windcd(k) = cd
          else if( iheat .eq. 7 ) then
            qsens = ta
            qlat  = ur
            !qlong = -cc  !change sign of long wave radiation given by ISAC
            qlong = abs(cc)
            evap  = qlat / (2.5008e6 - 2.3e3 * tm)      !pom, gill, gotm
          else if( iheat .eq. 8 ) then
            !ddlon = xgv(k)
            !ddlat = ygv(k)
            !uub = uprv(lmin,k)
            !vvb = vprv(lmin,k)
            !call meteo_get_heat_extra(k,dp,uuw,vvw)
            !call heatmfsbulk(days,im,ih,ddlon,ddlat,ta,p,uuw,vvw,dp, &
           !&                   cc,tm,uub,vvb,qsens,qlat,qlong,evap,qswa,cd)
            !qss= fice_free*qswa  !albedo (monthly) already in qshort1 -> qswa  
            !if ( bwind ) windcd(k) = cd
          else
            write(6,*) 'iheat = ',iheat
            stop 'error stop qflux3d: value for iheat not allowed'
          end if

          qlong = fice_free * qlong
          qlat = fice_free * qlat
          qsens_orig = qsens
          qsens = fice_free * qsens + fice_cover * qice
          qrad =  - ( qlong + qlat + qsens )
          qtot = qss + qrad
          qsolar = qss

	  qsurface = qrad
	  qsbottom = 0.

	  ! only for one layer

	  hm = hm_base + zeta
	  !hm = 1.
	  hm = 0.001
	  !if( hm < hm_min) cycle
	  call heat2t(dt,hm,qss-qsbottom,qsurface,tm,tnew)
	  tm = tnew
	  write(6,'(f12.0,7f9.2)') dtime,qss,qsens,qlong,qlat,ta,hm,tm

	end do

	return
   99	continue
	stop 'error stop heat_test: read error'
	end

!**************************************************************

	subroutine get_heat_values(dtime,qs,ta,ur,twb,uw,cc,p,zeta,ierr)

	implicit none

	double precision dtime
	real qs,ta,ur,twb,uw,cc,p,zeta
	integer ierr

	logical bwind,bspeed
	integer iqunit,iwunit
	integer nvar,nintp
	integer datetime(2)
	integer iwtype
	real f(4),values(4)
	character*80 varline
	character*80 dir,qfile,wfile,zfile
	character*128 file

	integer, save :: id_q,id_w,id_z
	integer, save :: icall = 0
	logical iff_ts_has_data

	nintp = 2
	iwtype = 3
	ierr = 0

	bwind = .true.		!want wind
	bspeed = .false.	!has speed
	if( iwtype == 2 .or. iwtype == 4 ) then
	  write(6,*) 'cannot handle yet iwtype ',iwtype
	  stop 'error stop get_heat_values: iwtype'
	else if( iwtype == 0 ) then
	  bwind = .false.
	else if( iwtype == 3 ) then
	  bspeed = .true.
	end if

	if( icall == 0 ) then
	  icall = 1
	  dir = '/home/georg/work/ogs/bugs/heatflux/run_2005_heatflux/input'
	  qfile = 'qflux_test_2005_orario.csv'
	  wfile = 'wind05_06.txt'
	  zfile = 'dsl2005_2006gmt.txt'

	  file = trim(dir) // '/' // trim(qfile)
	  call ts_get_file_info(file,nvar)
	  if( nvar /= 4 ) stop 'error stop qflux: nvar/=4'
	  call iff_ts_init(dtime,file,nintp,nvar,id_q)

	  if( bwind ) then
	    file = trim(dir) // '/' // trim(wfile)
	    call ts_get_file_info(file,nvar)
	    if( nvar /= 3 ) stop 'error stop wind: nvar/=3'
	    call iff_ts_init(dtime,file,nintp,nvar,id_w)
	  end if

	  file = trim(dir) // '/' // trim(zfile)
	  call ts_get_file_info(file,nvar)
	  if( nvar /= 1 ) stop 'error stop wind: nvar/=1'
	  call iff_ts_init(dtime,file,nintp,nvar,id_z)

	end if

	ierr = -1
	if( .not. iff_ts_has_data(id_q,dtime) ) return
	if( bwind ) then
	  if( .not. iff_ts_has_data(id_w,dtime) ) return
	end if
	if( .not. iff_ts_has_data(id_z,dtime) ) return
	ierr = 0

	call iff_ts_intp(id_q,dtime,values)
	qs = values(1)
	ta = values(2)
	ur = values(3)
	cc = values(4)
	if( bwind ) then
	  call iff_ts_intp(id_w,dtime,values)
	  if( bspeed ) then
	    uw = values(1)
	  else
	    uw = sqrt( values(1)**2 + values(2)**2 )
	  end if
	  p = values(3)
	else
	  uw = 0.
	  p = 1013.25
	end if
	call iff_ts_intp(id_z,dtime,values)
	zeta = values(1)

        if( p > 10000 ) p = 0.01 * p   !Pascal to mb
        call rh2wb(ta,p,ur,twb)

	!write(6,'(f12.0,7f9.2)') dtime,qs,ta,ur,cc,uw,p,zeta

	end

!**************************************************************

