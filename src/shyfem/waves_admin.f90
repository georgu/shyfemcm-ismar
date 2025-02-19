
!--------------------------------------------------------------------------
!
!    Copyright (C) 2006,2008-2009,2013-2014,2019  Christian Ferrarin
!    Copyright (C) 2008  Andrea Cucco
!    Copyright (C) 2010-2011,2014-2019  Georg Umgiesser
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

! waves subroutines
!
! revision log :
!
! 18.10.2006	ccf	integrated into main tree
! 19.06.2008	aac&ccf	udate to 3D version
! 16.04.2009	ccf	update to new WWMII-2 version, both 2D and 3D
! 23.03.2010	ggu	changed v6.1.1
! 08.10.2010	ggu	changed VERS_6_1_13
! 17.02.2011	ggu	changed VERS_6_1_18
! 18.02.2011	ggu	compiler warnings/errors adjusted
! 25.10.2013	ccf	upgrade compatibility with WWMIII
! 04.11.2014	ccf	rewritten
! 05.12.2014	ggu	changed VERS_7_0_8
! 19.12.2014	ggu	changed VERS_7_0_10
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 21.01.2015	ggu	computing fetch for geographical coordinates (bug fix)
! 10.02.2015	ggu	randomixe change of coordinates if on node (bug fix)
! 26.02.2015	ggu	changed VERS_7_1_5
! 26.02.2015	ggu	changed VERS_7_1_6
! 05.06.2015	ggu	changed VERS_7_1_12
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 24.07.2015	ggu	changed VERS_7_1_82
! 23.09.2015	ggu	changed VERS_7_2_4
! 10.10.2015	ggu	changed VERS_7_3_2
! 12.10.2015	ggu	changed VERS_7_3_3
! 12.10.2015	ggu	changed VERS_7_3_4a
! 10.03.2016	ggu	in parametric wave module fix segfault (allocatable)
! 15.04.2016	ggu	parametric wave module cleaned
! 28.04.2016	ggu	changed VERS_7_5_9
! 30.05.2016	ggu	changed VERS_7_5_11
! 10.06.2016	ggu	changed VERS_7_5_13
! 12.04.2017	ggu	routines integrated to compute bottom stress, new module
! 09.05.2017	ggu	changed VERS_7_5_26
! 13.06.2017	ggu	changed VERS_7_5_29
! 05.12.2017	ggu	changed VERS_7_5_39
! 03.04.2018	ggu	changed VERS_7_5_43
! 01.02.2019	ggu	bug fix in parwaves: do not compute H/P for wind == 0
! 10.02.2019	ggu	bug fix for FPE (GGUZ0)
! 12.02.2019	ccf	stress computed in substress.f
! 11.11.2020	ggu	get_ice_all() renamed to get_ice_cover_all()
! 20.03.2022	ggu	converted convert_time -> convert_time_d
! 30.03.2022	ggu	bug: in write_wwm nlev was not set before call
! 09.05.2023    lrp     introduce top layer index variable
! 25.01.2024    ggu     in parwaves() introduced OMP parallelization
! 12.02.2025	ccf	removed wwm
!
!**************************************************************
! DOCS  START   S_wave_par
!
! This empirical wave module is used to calculate the wave height 
! and period from wind speed, fetch and depth using the EMPIRICAL 
! PREDICTION EQUATIONS FOR SHALLOW WATER \cite{shoreprot:84}.
!
! Activate setting iwave = 1 in the .str file
!
! The wave module writes in the WAV file the following output:
! \begin{itemize}
! \item significant wave height [m], variable 231
! \item mean wave period [s], variable 232
! \item mean wave direction [deg], variable 233
! \end{itemize}
!
! The time step and start time for writing to file WAV 
! are defined by the parameters |idtwav| and |itmwav| in the |waves|
! section.  If |idtwav| is not
! defined, then the wave module does not write any results. The wave 
! results can be plotted using |plots -wav|.
!
! DOCS  END
!**************************************************************

        subroutine init_wave

! initialize arrays 

	use mod_waves

	implicit none

	real getpar		!get parameter function

        iwave = nint(getpar('iwave'))
	if ( iwave .eq. 1 ) then
	  call parwaves
	end if

	end subroutine init_wave

!******************************************************************

        subroutine get_wave_values(k,wh,wmp,wpp,wd)

! returns significant wave heigh, wave periods (mean and peak) and 
! mean wave direction

	use mod_waves

        implicit none

        integer k
        real wh                 !sign. wave height [m]
        real wmp                !mean wave period [s]
        real wpp                !peak wave period [s]
        real wd                 !mean wave direction [deg]

        wh  = waveh(k)
        wmp = wavep(k)
        wpp = wavepp(k)
        wd  = waved(k)

        end subroutine get_wave_values

!*********************************************************************

        function has_waves() 

! gives indication if waves are computed

	use mod_waves

        implicit none

        logical has_waves

        if( iwave .le. 0 ) then
          has_waves = .false.
	else
          has_waves = .true.
	end if

        end function has_waves

!******************************************************************
!******************************************************************
!******************************************************************
!******************************************************************
!******************************************************************
!
! This routine is used to calculate the wave height and period
! from wind speed, fetch and depth using the EMPIRICAL PREDICTION
! EQUATIONS FOR SHALLOW WATER (Shore Protection Manual, 1984).
! It considers a homogeneous wind field all over the domain.
! It works only with cartesian coordinate system.
!
! notes :
!
! Hs            significant wave height
! Tm            mean period
! Hrms          rms height (Hrms = Hs/sqrt(2))
! Tp            peak period (Tm = 0.781 Tp)
!
! Tm = 11 sqrt(Hs/g)    if no info on Tm is available
!
! for sediment transport use T=Tp and Uw = sqrt(2) Urms
!
! Uw            orbital velocity
! Urms          std. dev. of orbital velocity
!
! Dispersion relation:
!
! o             omega, frequency, o = 2 pi / T
! k             wave number, k = 2 pi / L
! L             wave length
!
! o**2 = gk tanh(kh)    dispersion relation
!
! zeta = o**2 h / g
! eta = k h
!
! ==> zeta = eta tanh(eta)
!
! easy computation (better than 1 %)
!
! zeta < 1      eta = sqrt(zeta) * (1 + 0.2 * zeta)
! zeta > 1      eta = zeta( 1 + 0.2 * exp(2 - 2 * zeta)
!
! Orbital velocity:
!
! Uw = pi H / (T sinh(kh))
!
!**************************************************************

!==============================================================
	module mod_parwaves
!==============================================================

	implicit none

        real, parameter :: z0 = 5.e-4
        real, parameter :: awice = 1.	!use wave reduction due to ice cover

! --- input variable

        real, save, allocatable :: winds(:) !wind speed at 10m [m/s]
        real, save, allocatable :: windd(:) !wind direction [degree north]
        real, save, allocatable :: fet(:)   !wind fetch length [m]
        real, save, allocatable :: daf(:)   !averaged depth along the fetch [m]

! --- output variable

        real, save, allocatable :: waeh(:)	!wave height [m]
        real, save, allocatable :: waep(:)	!wave period [s]
        real, save, allocatable :: waed(:)	!wave direction (same as wind)

!==============================================================
	end module mod_parwaves
!==============================================================

        subroutine parwaves

! called for iwave == 1

	use mod_meteo
	use mod_waves
	use basin
	use mod_parwaves

        implicit none

! --- aux variable

        real, save, allocatable :: v1v(:)	!aux variable
        real, save, allocatable :: icecover(:)	!ice cover

! --- local variable

	logical debug
        !integer ie,icount,ii,k
	integer id,nvar
	double precision dtime

        real getpar
        real depele             !element depth function [m]

        logical has_output_d,next_output_d

        integer, save :: icall = 0		!initialization parameter

	debug = .true.
	debug = .false.

! ----------------------------------------------------------
! Initialization
! ----------------------------------------------------------

        if( icall .le. -1 ) return

        if( icall .eq. 0 ) then

!         --------------------------------------------------
!         Initialize state variables
!         --------------------------------------------------

          iwave = nint(getpar('iwave'))
          if( iwave .le. 0 ) icall = -1
          if( iwave .gt. 1 ) icall = -1
          if( icall .le. -1 ) return

	  allocate(waeh(nel))
	  allocate(waep(nel))
	  allocate(waed(nel))
          waeh = 0.
          waep = 0.
          waed = 0.

	  allocate(winds(nel),windd(nel))
	  allocate(fet(nel),daf(nel))
	  allocate(icecover(nkn),v1v(nkn))

!         --------------------------------------------------
!         Initialize output
!         --------------------------------------------------

	  nvar = 3
          call init_output_d('itmwav','idtwav',da_wav)
          if( has_output_d(da_wav) ) then
	    call shyfem_init_scalar_file('wave',nvar,.true.,id)
	    da_wav(4) = id
          end if

          write(6,*) 'parametric wave model initialized...'
          icall = 1
        endif

! -------------------------------------------------------------------
! normal call
! -------------------------------------------------------------------

	call get_ice_cover_all(icecover)

!       -------------------------------------------------------------
!	get wind speed and direction
!       -------------------------------------------------------------

	call compute_wind_values(icecover)

!       -------------------------------------------------------------
!	get wind fetch (fet is fetch on elements)
!       -------------------------------------------------------------

        call fetch(windd,fet,daf)

!       -------------------------------------------------------------------
!       compute wave values
!       -------------------------------------------------------------------

        call compute_wave_parameters

!       -------------------------------------------------------------------
!       copy to global values
!       -------------------------------------------------------------------

        call e2n2d(waeh,waveh,v1v)
        call e2n2d(waep,wavep,v1v)
        call e2n2d(waed,waved,v1v)

	wavepp = wavep                  !peak period is not computed

!       -------------------------------------------------------------------
!       write results to file (WAV)
!       -------------------------------------------------------------------

        if( next_output_d(da_wav) ) then
	  id = nint(da_wav(4))
	  call get_act_dtime(dtime)
	  call shy_write_scalar_record2d(id,dtime,231,waveh)
	  call shy_write_scalar_record2d(id,dtime,232,wavep)
	  call shy_write_scalar_record2d(id,dtime,233,waved)
	  call shy_sync(id)
	end if

!       -------------------------------------------------------------------
!       end of routine
!       -------------------------------------------------------------------

        end

!**************************************************************
!**************************************************************
!**************************************************************

        subroutine compute_wave_parameters

        use basin

        implicit none

        integer ie
        integer nchunk

!$OMP PARALLEL 
!$OMP SINGLE

        call omp_compute_chunk(nel,nchunk)

        do ie=1,nel
!$OMP TASK FIRSTPRIVATE(ie) &
!$OMP&     SHARED(nchunk) &
!$OMP&     DEFAULT(NONE)
          call compute_waves_on_element(ie)
!$OMP END TASK
        end do

!$OMP END SINGLE
!$OMP TASKWAIT  
!$OMP END PARALLEL      

        end

!**************************************************************

	subroutine compute_waves_on_element(ie)

	use basin
	use mod_meteo
	use mod_parwaves

	implicit none

	integer ie

        real ah1,ah2,ah3,eh1,eh2,eh3,eh4
        real at1,at2,at3,et1,et2,et3,et4
!------------------------------------------------------ Hurdle and Stive
!        parameter(ah1=0.25,ah2=0.6,eh1=0.75)
!        parameter(eh2=0.5,ah3=4.3e-05,eh3=1.,eh4=2.)
!        parameter(at1=8.3,at2=0.76,et1=0.375)
!        parameter(et2=1./3.,at3=4.1e-05,et3=1.,et4=3.)
!------------------------------------------------------ SPM
        parameter(ah1=0.283,ah2=0.53,eh1=3./4.)
        parameter(eh2=1.,ah3=0.00565,eh3=1./2.,eh4=1.)
        parameter(at1=7.54,at2=0.833,et1=3./8.)
        parameter(et2=1.,at3=0.0379,et3=1./3.,et4=1.)

        real, parameter :: g = 9.81		!gravity acceleration [m2/s]

	integer icount
	real dep,depe
	real wis,wid
	real gh,gx,hg
	real auxh,auxh1,auxt,auxt1
	real hbr			!limiting wave height [m]

	real depele

          icount = 1

!         -----------------------------------------------------------------
!	  get averaged depth along the fetch
!         -----------------------------------------------------------------

          dep = daf(ie)
          depe = depele(ie,+1)
          dep = depele(ie,+1)
10        continue

!         -----------------------------------------------------------------
!	  calculate wave height, period and direction
!         -----------------------------------------------------------------

	  wis = winds(ie)
	  wid = windd(ie)

!         -----------------------------------------------------------------
!	  method of SPM
!         -----------------------------------------------------------------

	  if( wis > 0. ) then
            gh = (g*dep)/(wis**2.)
            gx = (g*fet(ie))/(wis**2.)
            hg = dep / (g*wis**2.)
            auxh = ah2*gh**eh1
            auxh1 = ah2*hg**eh1
            auxt = at2*gh**et1
            auxt1 = ah2*gx**eh1

            waeh(ie) = (tanh(auxh))**eh4
            waeh(ie) = (ah3*gx**eh3) / waeh(ie)
            waeh(ie) = (tanh(waeh(ie)))**eh2
            waeh(ie) = ah1*tanh(auxh)*waeh(ie)
            waeh(ie) = waeh(ie) * wis**2 / g
          
            waep(ie) = (tanh(auxt))**et4
            waep(ie) = (at3*gx**et3) / waep(ie)
            waep(ie) = (tanh(waep(ie)))**et2
            waep(ie) = at1*tanh(auxt)*waep(ie)
            waep(ie) = waep(ie) * wis / g
	  else
	    waeh(ie) = 0.
	    waep(ie) = 0.
	  end if

!          waeh(ie) = 0.283 * tanh(0.530*(gh**(3./4.)))*
!     %            tanh((0.00565*(gx**0.5))/
!     %            (tanh(0.530*(gh**(3./4.)))))*((wis**2)/g)
!
!          waep(ie) = 7.54*tanh(0.833*(gh**(3./8.)))*
!     %            tanh((0.0379*(gx**(1./3.)))/
!     %            (tanh(0.833*(gh**(3./8.)))))*(wis/g)
!

!         -----------------------------------------------------------------
!	  method of hurdle and stive
!         -----------------------------------------------------------------

!          waeh(ie) = (tanh(auxh))**eh4
!          waeh(ie) = (ah3*gx**eh3) / waeh(ie)
!          waeh(ie) = (tanh(waeh(ie)))**eh2
!          waeh(ie) = ah1*tanh(auxh)*waeh(ie)
!          waeh(ie) = waeh(ie) * wis**2 / g
         
!          waep(ie) = (tanh(auxt1))**et4
!          waep(ie) = (at3*gx**et3) / waep(ie)
!          waep(ie) = (tanh(waep(ie)))**et2
!          waep(ie) = at1*tanh(auxt)*waep(ie)
!          waep(ie) = waep(ie) * wis / g

          waed(ie) = wid

!         -----------------------------------------------------------------
!	  limiting wave height
!         -----------------------------------------------------------------

          hbr = 0.50*depe
          if( waeh(ie).gt.hbr ) then
            if (icount.gt.1) then
              waeh(ie) = hbr
              go to 20
            end if
            dep = depe
            icount = 2
            goto 10
          end if 
20        continue

        end

!**************************************************************

	subroutine compute_wind_values(icecover)

	use basin
	use mod_meteo
	use mod_parwaves

	implicit none

	real icecover(nkn)
	
	integer ie,ii,k
        integer nchunk
	real wx,wy,fice
	
!$OMP PARALLEL 
!$OMP SINGLE

        call omp_compute_chunk(nel,nchunk)

	do ie=1,nel

!$OMP TASK FIRSTPRIVATE(ie) &
!$OMP&     PRIVATE(wx,wy,ii,k,fice) &
!$OMP&     SHARED(nel,nchunk,icecover,nen3v,wxv,wyv,winds,windd) &
!$OMP&     DEFAULT(NONE)

	  wx = 0.
	  wy = 0.
	  do ii=1,3
	    k = nen3v(ii,ie)
	    fice = 1. - awice*icecover(k)
	    wx = wx + fice * wxv(k)
	    wy = wy + fice * wyv(k)
	  end do
	  wx = wx / 3.
	  wy = wy / 3.
          call c2p(wx,wy,winds(ie),windd(ie))

!$OMP END TASK

	end do

!$OMP END SINGLE
!$OMP TASKWAIT  
!$OMP END PARALLEL      

	end

!**************************************************************

        subroutine fetch(windd,fet,daf)

! This subroutine computes the wind fetch for each element of the
! grid given the wind direction.

	use basin, only : nkn,nel,ngr,mbw
	use coordinates

        implicit none
  
        real windd(nel)		!wind direction [degree north]
        real fet(nel)		!wind fetch length [m] (return)
        real daf(nel)		!averaged depth along the fetch [m] (return)

        real xe,ye		!element point coordinates [m]
	real fff,ddd
        real rad,wdir,wid
        integer ie,ierr,nchunk

        rad = 45. / atan (1.)

! --- loop over elements

!$OMP PARALLEL 
!$OMP SINGLE

        call omp_compute_chunk(nel,nchunk)

        do ie = 1,nel
          
!$OMP TASK FIRSTPRIVATE(ie,rad) &
!$OMP&     PRIVATE(xe,ye,wid,wdir,ierr,fff,ddd) &
!$OMP&     SHARED(nel,nchunk,windd,fet,daf) &
!$OMP&     DEFAULT(NONE)
        
          call baric_cart(ie,xe,ye)

	  wid = windd(ie)
          wdir = wid / rad		!from deg to rad

	  ierr = 0
	  call fetch_element(ie,xe,ye,wdir,fff,ddd,ierr)

	  if( ierr .ge. 1000 ) then
	    write(6,*) 'warning: iteration exceeded: ',ie
            stop 'error stop fetch: iterations'
	  end if

          fet(ie) = fff
          daf(ie) = ddd
	  
!$OMP END TASK

        end do

!$OMP END SINGLE
!$OMP TASKWAIT  
!$OMP END PARALLEL      

	end

!**************************************************************

        subroutine fetch_element(ie,xein,yein,wdir,fff,ddd,ierr)

! This subroutine computes the wind fetch for element ie
! given the wind direction.

	use basin, only : nkn,nel,ngr,mbw

        implicit none
  
	integer ie		!element
	real xein,yein		!initial point
	real wdir		!wind direction
	real fff		!fetch computed (return)
	real ddd 		!average depth computed (return)
	integer ierr		!error code
				!on entry if >0 write debug information
				!on return number of iterations executed

        real xe,ye		!element point coordinates [m]
        real xnew,ynew		!new coordinates [m]
        real de			!distance between points [m]
        real depele             !element depth function [m]
        real dep		!element depth [m]
        integer iie,ii,ienew,icount,ieold
	logical bdebug

        iie = ie
        ieold = ie
	bdebug = ierr > 0
	fff = 0.
	ddd = 0.
	xe = xein
	ye = yein

! --- calculate fetch and averaged depth along the fetch

	if( bdebug ) then
	  write(156,*) '=========================='
	  write(156,*) iie
	end if

        icount = 0
	ienew = 1	!just to enter the while loop

	do while( ienew > 0 .and. icount < 1000 )
          call intersect(iie,xe,ye,wdir,ienew,xnew,ynew,ieold,bdebug)
          dep = depele(iie,+1)
          de = ((xnew-xe)**2 + (ynew-ye)**2)**0.5
          fff = fff + de
          ddd = ddd + dep*de
          icount = icount + 1
	  if( bdebug ) then
	    write(156,*) '-------------------'
	    write(156,*)  icount
	    write(156,*)  iie,ienew,ieold
	    write(156,*)  xe,ye,xnew,ynew
	    write(156,*)  de,dep,fff,ddd
	    write(156,*) '-------------------'
	  end if
          ieold = iie
          iie = ienew
          xe = xnew
          ye = ynew
	end do

        if( fff > 0. ) ddd = ddd/fff
        if(ienew.lt.0) fff = fff + 50000.	!open boundary

	if( bdebug ) then
	  write(156,*) icount,fff,ddd
	  write(156,*) '=========================='
	end if
 
	ierr = icount

	return
        end
           
!******************************************************************

        subroutine intersect(iie,x,y,wdir,ien,xn,yn,ieold,bdebug)

! this routine computes the coordinate of the intersection beetwen the
! line and one of the border line of the element

	use mod_geom
	use basin
	use coordinates

        implicit none

        integer iie		!element number
        real x,y		!start point cooridnates [m]
        real wdir		!direction to search [radians]
        integer ien		!next element number (return)
        real xn,yn		!intersection coordinates [m] (return)
	integer ieold
	logical bdebug

        real x0(3),y0(3)        !element vertices coordinates [m]
        real xg0(3),yg0(3)      !element vertices coordinates [degrees]
        real x3,y3,x4,y4	!node points coordiantes [m]
        real xi,yi		!intersection point coordinate [m]
	real xf,yf		!far away point
	real d			!distance
	real rx,ry
	double precision a(3),b(3),c(3)
        integer iint,i,ii,iii
        integer ienew

        integer segsegint	!intersection function

	d = 5000000.
        xf = x + d*sin(wdir)
        yf = y + d*cos(wdir)

        ien = 0
	xn = 0.
	yn = 0.
        call getexy_cart(iie,x0,y0)
	if( bdebug ) then
	  write(156,*) '.............'
	  write(156,*) iie
	  write(156,*) x0
	  write(156,*) y0
	  write(156,*) x,y,xf,yf
	end if
 
        do i = 1,3

          ii=mod(i,3)+1
          iii=mod(ii,3)+1

          x3=x0(ii)
          y3=y0(ii)
          x4=x0(iii)
          y4=y0(iii)

2         continue
          iint = segsegint(x,y,xf,yf,x3,y3,x4,y4,xi,yi)

          ienew = ieltv(i,iie)
	
	  if( bdebug ) then
	    write(156,*) i,iint,iie,ienew,ieold
	  end if

          if(iint.gt.0.and.ienew.ne.ieold)then	!intersection
            if(iint.eq.3)then	 		!intersection with node
              !x = x + 1.
              !y = y + 1.
	      call random_number(rx)
	      call random_number(ry)
              x = x + 10.*(rx-0.5)
              y = y + 10.*(ry-0.5)
	      if( bdebug ) then
	  	write(156,*) 9,0,0,0,0
		write(156,*) 'warning: node intersection: ',i,x,y
	      end if
              go to 2
            else
	      if( bdebug .and. ien .gt. 0 ) then
	  	write(156,*) 9,0,0,0,0
		write(156,*) 'warning: ien already set: ',ien,xn,yn
	      end if
              xn = xi
              yn = yi
              ien = ienew
            end if
          end if

        end do

	if( bdebug ) then
	  if( ien .eq. 0 ) then
	    write(156,*) 9,0,0,0,0
	    write(156,*) 'warning: ien not set: ',ien,xn,yn
	  end if
	  write(156,*) 0,0,0,0,0
	  write(156,*) ien,xn,yn
	  write(156,*) '.............'
	end if

        end

!******************************************************************
!******************************************************************
