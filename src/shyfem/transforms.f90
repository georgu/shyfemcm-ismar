
!--------------------------------------------------------------------------
!
!    Copyright (C) 2003-2005,2008,2010-2017,2019  Georg Umgiesser
!    Copyright (C) 2004,2015  Christian Ferrarin
!    Copyright (C) 2012  Debora Bellafiore
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

! routines for various transformations
!
! contents :
!
! subroutine ttov			transforms transports to velocities
! subroutine vtot			transforms velocities to transports
! subroutine uvtop0			transforms bar. transp. to nod. veloc.
! subroutine uvtopr			transforms velocities to nodal values
! subroutine prtouv			transforms nodal values to element vel.
! subroutine uvint			computation of barotropic transports
! subroutine austau(vv)			computes aux vectors for austausch term
! subroutine baro2l			distribute barotropic velocities 
! subroutine setxv			sets obsolete data structure xv
! subroutine getuv(l,k,u,v)             accessor routine to get velocities u/v
! subroutine copy_uvz                   copies u/v/z to old time step
! subroutine make_prvel                 makes print velocities and xv
! subroutine init_uv                    initializes uvz values
! subroutine e2n2d(elv,nov,aux)         transforms element to nodal values
! subroutine n2e2d(nov,elv)             transforms nodal to element values 2D
! subroutine n2e3d(nlvdi,nov,elv)       transforms nodal to element values 3D
!
! revision log :
!
! 10.08.2003	ggu	new routines copy_uvz, make_prvel, init_uvz
! 18.09.2003	ggu	new routine e2n2d
! 04.12.2003	ggu	new routine n2e2d
! 30.03.2004	ccf	new routine n2e2d and n2e3d
! 15.10.2004	ggu	new routine smagorinsky started
! 14.01.2005	ggu	bug (nlvdi) in n2e3d fixed
! 24.02.2005	ggu	smagorinsky into subdif.f
! 15.03.2005	ggu	austv eliminated, austau() into subdif.f
! 27.06.2005	ggu	bug in vtot corrected
! 10.04.2008	ggu	copy velocities at nodes in copy_uvz()
! 01.03.2010	ggu	new version of n2e3d()
! 11.03.2010	ggu	new routine check_volume(); init w only if no restart
! 23.03.2010	ggu	changed v6.1.1
! 15.12.2010	ggu	changed VERS_6_1_14
! 16.02.2011	ggu	new routine e2n3d() and e2n3d_minmax()
! 27.01.2012	dbf&ggu	routines adapted for sigma levels
! 05.12.2013	ggu	changed VERS_6_1_70
! 18.06.2014	ggu	changed VERS_6_1_77
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 05.06.2015	ggu	changed VERS_7_1_12
! 13.07.2015	ggu	changed VERS_7_1_51
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 24.07.2015	ggu	changed VERS_7_1_82
! 18.09.2015	ggu	changed VERS_7_2_3
! 23.09.2015	ggu	changed VERS_7_2_4
! 05.11.2015	ggu	changed VERS_7_3_12
! 03.12.2015	ccf&ggu	code optimized
! 16.12.2015	ggu	changed VERS_7_3_16
! 18.12.2015	ggu	changed VERS_7_3_17
! 07.04.2016	ggu	new routine aver_nodal()
! 15.04.2016	ggu	changed VERS_7_5_8
! 19.05.2016	ggu	use where construct where possible
! 17.06.2016	ggu	adjust code to reflect that wprv now starts from 1
! 09.09.2016	ggu	changed VERS_7_5_17
! 05.12.2017	ggu	changed VERS_7_5_39
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 30.03.2021	ggu	new routine compute_velocities()
! 07.04.2022	ggu	ie_mpi introduced computing print velocities
! 10.04.2022	ggu	ie_mpi and double in uvint (compiler issue with INTEL)
! 09.05.2023    lrp     introduce top layer index variable
! 05.06.2023    lrp     introduce z-star
! 20.03.2026    ggu     new routines compute_derived_variables() and uvint_old()
!
!****************************************************************************

	subroutine ttov

! transforms transports to velocities

	use mod_layer_thickness
	use mod_hydro_vel
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer ie

	if( .not. mod_layer_thickness_is_initialized() ) then
	  write(6,*) 'layer thickness is not initialized'
	  stop 'error stop ttov: no layerthickness'
	end if

	ulnv = 0.
	vlnv = 0.
	where( hdenv > 0. )
	  ulnv = utlnv / hdenv
	  vlnv = vtlnv / hdenv
	end where

	end

!******************************************************************

	subroutine vtot

! transforms velocities to transports

	use mod_layer_thickness
	use mod_hydro_vel
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	utlnv = ulnv * hdenv
	vtlnv = vlnv * hdenv

	end

!******************************************************************

	subroutine uvtop0

! transforms barotropic transports to nodal velocities

	use mod_geom_dynamic
	use mod_depth
	use mod_hydro_baro
	use mod_hydro_print
	use mod_hydro
	use evgeom
	use basin
	use shympi

	implicit none

! local
	logical bcolin
	integer ie,k,ii,ie_mpi
	real aj,zm,hm
	real vv(nkn)
! function
	real getpar
	integer iround

	bcolin=iround(getpar('iclin')).ne.0

	up0v = 0.
	vp0v = 0.
	vv   = 0.

	do ie_mpi=1,nel
	  ie = ip_sort_elem(ie_mpi)
	  if( iwegv(ie) /= 0 ) cycle
	  aj=ev(10,ie)
	  zm=0.
	  do ii=1,3
	    zm=zm+zenv(ii,ie)
	  end do
	  zm=zm/3.
	  hm=hev(ie)
	  if(.not.bcolin) hm=hm+zm
	  do ii=1,3
	    k=nen3v(ii,ie)
	    vv(k)=vv(k)+aj
	    up0v(k)=up0v(k)+aj*unv(ie)/hm
	    vp0v(k)=vp0v(k)+aj*vnv(ie)/hm
	  end do
	end do

       !call shympi_comment('shympi_elem: exchange up0v, vp0v')
        call shympi_exchange_and_sum_2d_nodes(up0v)
        call shympi_exchange_and_sum_2d_nodes(vp0v)
        call shympi_exchange_and_sum_2d_nodes(vv)

	where ( vv > 0. ) 
          up0v = up0v / vv
          vp0v = vp0v / vv
	end where

        call shympi_comment('uvtop0: exchange up0v, vp0v')
	call shympi_exchange_2d_node(up0v)
	call shympi_exchange_2d_node(vp0v)

	return
	end

!******************************************************************

	subroutine uvtopr

! transforms velocities to nodal values

	use mod_geom_dynamic
	use mod_hydro_print
	use mod_hydro_vel
	use evgeom
	use levels
	use basin
	use shympi

	implicit none

	integer ie,l,k,ii,ie_mpi
	integer lmax,lmin
	real aj
	!real vv(nlvdi,nkn)
	real, allocatable :: vv(:,:)

	allocate(vv(nlvdi,nkn))
	uprv = 0.
	vprv = 0.
	vv   = 0.

! baroclinic part

	do ie_mpi=1,nel
	  ie = ip_sort_elem(ie_mpi)
	  if ( iwegv(ie) /= 0 ) cycle
          lmax = ilhv(ie)
	  lmin = jlhv(ie)
	  aj=ev(10,ie)
	  do l=lmin,lmax
	    do ii=1,3
	      k=nen3v(ii,ie)
	      vv(l,k)=vv(l,k)+aj
	      uprv(l,k)=uprv(l,k)+aj*ulnv(l,ie)
	      vprv(l,k)=vprv(l,k)+aj*vlnv(l,ie)
	    end do
	  end do
	end do

        !call shympi_comment('shympi_elem: exchange uprv, vprv')
        call shympi_exchange_and_sum_3d_nodes(uprv)
        call shympi_exchange_and_sum_3d_nodes(vprv)
        call shympi_exchange_and_sum_3d_nodes(vv)

	where ( vv > 0. ) 
	  uprv = uprv / vv
	  vprv = vprv / vv
	end where

        call shympi_comment('uvtopr: exchange uprv, vprv')
	call shympi_exchange_3d_node(uprv)
	call shympi_exchange_3d_node(vprv)

! vertical velocities -> we compute average over one layer

	do l=1,nlv
	  wprv(l,:)=0.5*(wlnv(l,:)+wlnv(l-1,:))
	end do

	deallocate(vv)

	end

!******************************************************************

	subroutine prtouv

! transforms nodal values to element values (velocities)

	use mod_geom_dynamic
	use mod_hydro_print
	use mod_hydro_vel
	use levels
	use basin

	implicit none

	integer ie,l,k,ii
	integer lmin,lmax
	real u,v

! baroclinic part

	do ie=1,nel
	 if( iwegv(ie) .eq. 0 ) then
	  lmax = ilhv(ie)
	  lmin = jlhv(ie)
	  do l=lmin,lmax
	    u=0.
	    v=0.
	    do ii=1,3
	      k=nen3v(ii,ie)
	      u=u+uprv(l,k)
	      v=v+vprv(l,k)
	    end do
	    ulnv(l,ie)=u/3.
	    vlnv(l,ie)=v/3.
	  end do
	 else
	    ulnv(:,ie)=0.
	    vlnv(:,ie)=0.
	 end if
	end do

! vertical velocities -> from layer average to interface values

	wlnv(nlv,:) = 0.
	do l=nlv-1,0,-1
	  wlnv(l,:)=2.*wprv(l+1,:)-wlnv(l+1,:)
	end do
	wlnv(0,:) = 0.

	return
	end

!*****************************************************************

	subroutine uvint

! computation of barotropic part of transports

	use mod_hydro_baro
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw
	use shympi

	implicit none

	integer ie,l,ie_mpi
	double precision u,v	!needed for bit2bit compatibility with INTEL

	do ie_mpi=1,nel
	  ie = ip_sort_elem(ie_mpi)
	  u=0.
	  v=0.
	  do l=jlhv(ie),ilhv(ie)
	    u=u+utlnv(l,ie)
	    v=v+vtlnv(l,ie)
	  end do
	  unv(ie)=u
	  vnv(ie)=v
	end do

	return
	end

!*****************************************************************

	subroutine uvint_old

! computation of barotropic part of transports

	use mod_hydro_baro
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw
	use shympi

	implicit none

	integer ie,l,ie_mpi
	double precision u,v	!needed for bit2bit compatibility with INTEL

	do ie_mpi=1,nel
	  ie = ip_sort_elem(ie_mpi)
	  u=0.
	  v=0.
	  do l=jlhv(ie),ilhv(ie)
	    u=u+utlov(l,ie)
	    v=v+vtlov(l,ie)
	  end do
	  uov(ie)=u
	  vov(ie)=v
	end do

	return
	end

!*****************************************************************

	subroutine check_volume

! checks for negative volume (depth)
!
! only first layer has to be checked

	use mod_layer_thickness
	use mod_hydro
	use levels, only : nlvdi,nlv,jlhv,jlhkv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	logical bstop,bsigma
	integer nsigma
	integer k,ie,ke,iee,ii
	real h,z
	real hsigma

	integer ipext,ieext

	bstop = .false.

	call get_sigma(nsigma,hsigma)
	bsigma = nsigma .gt. 0

	if( .not. bsigma ) then		!not sure we need this - FIXME
	 do k=1,nkn
	  h = hdknv(jlhkv(k),k)
	  if( h .le. 0. ) then
	    z = znv(k)
	    ke = ipext(k)
	    write(6,*) 'negative depth in node (layer 1): '
	    write(6,*) '   node (int/ext): ',k,ke
	    write(6,*) '   depth,zeta    : ',h,z
	    bstop = .true.
	  end if
	 end do
	end if

	do ie=1,nel
	  h = hdenv(jlhv(ie),ie)
	  if( h .le. 0. ) then
	    iee = ieext(ie)
	    write(6,*) 'negative depth in elem (layer 1): '
	    write(6,*) '   element (int/ext): ',ie,iee
	    write(6,*) '   depth: ',h
	    write(6,*) '   water levels: ',(zenv(ii,ie),ii=1,3)
	    bstop = .true.
	  end if
	end do

	if( bstop ) then
	  write(6,*) 'Negative volumes in first layer found.'
	  write(6,*) 'Drying of an element with more than one layer'
	  write(6,*) 'is not allowed.'
	  write(6,*) 'This can be due to either a computed water level'
	  write(6,*) 'which is too low or because of instability.'
	  write(6,*) 'In the first case check the boundary data and/or'
	  write(6,*) 'increase the layer thickness of the first layer.'
	  write(6,*) 'In the second case look for the causes of the'
	  write(6,*) 'instability and maybe reduce the time step.'
	  stop 'error stop check_volume: negative volume'
	end if

	end

!*****************************************************************

	subroutine baro2l

! distribute barotropic velocities onto layers (only in dry elements)

	use mod_geom_dynamic
	use mod_hydro_baro
	use mod_hydro_vel
	use mod_layer_thickness
	use mod_hydro
	use levels
	use basin

	implicit none

	logical bsigma
	integer nsigma
	integer ie,ilevel,jlevel,ii,l
	real hsigma,weight,htot
  
! functions
        integer ieext

	call get_sigma(nsigma,hsigma)
	bsigma = nsigma .gt. 0

	if( bsigma ) return

	do ie=1,nel
	  if( iwegv(ie) .gt. 0 ) then	!dry
	    ilevel=ilhv(ie)
	    jlevel=jlhv(ie)
	    !wettying with more then one layer
	    !allowed: jlevel.ne.ilevel case
	    !we put a uniform velocity across layers
	    htot = 0.0
            do l=jlevel,ilevel
	      htot = htot + hdenv(l,ie)  
	    end do
	    do l=jlevel,ilevel		
	      weight = hdenv(l,ie)/htot
	      utlnv(l,ie) = unv(ie) * weight
	      vtlnv(l,ie) = vnv(ie) * weight
	    end do
	  end if
	end do

	end

!******************************************************************

	subroutine setxv

! sets obsolete data structure xv

	use mod_hydro_print
	use mod_hydro
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	xv(1,:) = up0v(:)
	xv(2,:) = vp0v(:)
	xv(3,:) = znv(:)

	end

!******************************************************************

	subroutine getuv(l,k,u,v)

! accessor routine to get velocities u/v

	use mod_hydro_print

	implicit none

        integer l,k
        real u,v

        u = uprv(l,k)
        v = vprv(l,k)

        end

!******************************************************************

	subroutine copy_uvz

! copies u/v/z to old time step

	use mod_hydro_baro
	use mod_hydro_print
	use mod_hydro_vel
	use mod_hydro

	implicit none

	zov   = znv
	zeov  = zenv

	utlov = utlnv
	vtlov = vtlnv
	wlov  = wlnv

	uov   = unv			!$$UVBARO
	vov   = vnv
	ulov  = ulnv
	vlov  = vlnv

        upro  = uprv
        vpro  = vprv

	end

!******************************************************************

	subroutine make_prvel

! makes print velocities (and xv) from new level arrays

	use basin
	use shympi

	implicit none

	call uvtopr	!computes uprv,vprv,wprv
	call uvtop0	!computes up0v,vp0v
	call setxv	!sets xv from up0v,vp0v,znv

	end

!******************************************************************

	subroutine compute_velocities

! computes horizontal velocities from zenv, utlnv, vtlnv, hdenv

	implicit none

	call ttov			!velocities ulnv/vlnv
	call make_prvel			!nodal values uprv/vprv/wprv/up0v/vp0v

	end

!******************************************************************

	subroutine compute_derived_variables

! computes a lot of things for old time level

	use basin
	use mod_hydro
	use mod_hydro_baro

	implicit none

	call ttov			!velocities ulnv/vlnv
	call make_prvel			!nodal values uprv/vprv/wprv/up0v/vp0v
	call uvint			!barotropic transports unv/vnv

	end

!******************************************************************

	subroutine compute_old_derived_variables

! computes a lot of things for old time level

	use basin
	use mod_hydro
	use mod_hydro_vel
	use mod_layer_thickness

	implicit none

        ulov = 0.
        vlov = 0.
        where( hdeov > 0. )
          ulov = utlov / hdeov
          vlov = vtlov / hdeov
        end where

	call uvint_old			!old barotropic transports uov/vov

	end

!******************************************************************

	subroutine init_uv

! initializes uvz values from zenv, utlnv, vtlnv, hdenv

	use basin

	implicit none

	logical rst_use_restart
	real dzeta(nkn)

	if( .not. rst_use_restart(5) ) then
	  call hydro_vertical(dzeta)	!vertical velocities
	end if

	call compute_velocities

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine e2n2d(elv,nov,aux)

! transforms element values to nodal values (weights are area)
!
! (2D version)

	use evgeom
	use basin

	implicit none

! arguments
        real elv(nel)     !array with element values (in)
        real nov(nkn)     !array with nodal values (out)
        real aux(nkn)     !aux array (nkn)

! local
        integer k,ie,ii
        real area,value

!-----------------------------------------------------------
! initialize arrays
!-----------------------------------------------------------

        nov = 0.
        aux = 0.

!-----------------------------------------------------------
! accumulate values
!-----------------------------------------------------------

        do ie=1,nel
          area = 4.*ev(10,ie)
          value = elv(ie)
          do ii=1,3
            k = nen3v(ii,ie)
            nov(k) = nov(k) + area*value
            aux(k) = aux(k) + area
          end do
        end do

!-----------------------------------------------------------
! compute final value
!-----------------------------------------------------------

        nov = nov / aux

!-----------------------------------------------------------
! end of routine
!-----------------------------------------------------------

        end

!******************************************************************

	subroutine e2n3d(nlvddi,elv,nov,aux)

! transforms element values to nodal values (weights are area)
!
! (3D version)

	use evgeom
	use levels
	use basin

	implicit none

! arguments
	integer nlvddi
        real elv(nlvddi,nel)     !array with element values (in)
        real nov(nlvddi,nkn)     !array with nodal values (out)
        real aux(nlvddi,nkn)     !aux array (nkn)

! local
        integer k,ie,ii,l,lmax,lmin
        real area,value

!-----------------------------------------------------------
! initialize arrays
!-----------------------------------------------------------

        nov = 0.
        aux = 0.

!-----------------------------------------------------------
! accumulate values
!-----------------------------------------------------------

        do ie=1,nel
          area = 4.*ev(10,ie)
	  lmax = ilhv(ie)
	  lmin = jlhv(ie)
	  do l=lmin,lmax
            value = elv(l,ie)
            do ii=1,3
              k = nen3v(ii,ie)
              nov(l,k) = nov(l,k) + area*value
              aux(l,k) = aux(l,k) + area
	    end do
          end do
        end do

!-----------------------------------------------------------
! compute final value
!-----------------------------------------------------------

	where ( aux > 0. ) 
	  nov = nov / aux
	end where

!-----------------------------------------------------------
! end of routine
!-----------------------------------------------------------

        end

!******************************************************************

	subroutine e2n3d_minmax(mode,nlvddi,elv,nov)

! transforms element values to nodal values (no weights - use min/max)
!
! (3D version)

	use levels
	use basin

	implicit none

! arguments
	integer mode		!min (-1) or max (+1)
	integer nlvddi		!vertical dimension
        real elv(nlvddi,nel)      !array with element values (in)
        real nov(nlvddi,nkn)      !array with nodal values (out)

! local
        integer k,ie,ii,l,lmax,lmin
        real rinit,value

!-----------------------------------------------------------
! initialize arrays
!-----------------------------------------------------------

	if( mode .eq. -1) then
	  rinit = 1.e+30
	else if( mode .eq. 1 ) then
	  rinit = -1.e+30
	else
	  stop 'error stop e2n3d_minmax: unknown mode'
	end if

        nov = rinit

!-----------------------------------------------------------
! accumulate values
!-----------------------------------------------------------

        do ie=1,nel
	  lmax = ilhv(ie)
	  lmin = jlhv(ie)
	  do l=lmin,lmax
            value = elv(l,ie)
            do ii=1,3
              k = nen3v(ii,ie)
	      if( mode .eq. 1 ) then
                nov(l,k) = max(nov(l,k),value)
	      else
                nov(l,k) = min(nov(l,k),value)
	      end if
	    end do
          end do
        end do

!-----------------------------------------------------------
! end of routine
!-----------------------------------------------------------

        end

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine n2e2d(nov,elv)

! transforms nodal values to element values
!
! (2D version)

	use basin

        implicit none

        real nov(nkn)     !array with nodal values (in)
        real elv(nel)     !array with element values (out)

        integer k,ie,ii
        real acu,value

!-----------------------------------------------------------
! convert values
!-----------------------------------------------------------

        do ie=1,nel
          acu = 0.
          do ii=1,3
            k = nen3v(ii,ie)
            value = nov(k)
            acu = acu + value
          end do
          elv(ie) = acu / 3.
        end do

!-----------------------------------------------------------
! end of routine
!-----------------------------------------------------------

        end

!******************************************************************

        subroutine n2e3d(nlvddi,nov,elv)

! transforms nodal values to element values
!
! (3D version)

	use levels
	use basin

        implicit none

! arguments
        integer nlvddi		!vertical dimension of arrays
        real nov(nlvddi,nkn)	!array with nodal values (in)
        real elv(nlvddi,nel)	!array with element values (out)

! local
        integer k,ie,ii,l,lmax,lmin
        real acu,value

!-----------------------------------------------------------
! convert values
!-----------------------------------------------------------

        do ie=1,nel
          lmax = ilhv(ie)
	  lmin = jlhv(ie)
          do l = lmin,lmax
	    acu = 0.
            do ii=1,3
              k = nen3v(ii,ie)
              value = nov(l,k)
              acu = acu + value
	    end do
            elv(l,ie) = acu / 3.
          end do
	end do

!-----------------------------------------------------------
! end of routine
!-----------------------------------------------------------

        end

!******************************************************************

	subroutine aver_nodal(val,aver)

! computes average of val (defined on nodes) over total basin

	use evgeom
	use basin

	implicit none

        real val(nkn)     !array with element values (in)
	real aver	  !average (out)

        integer k,ie,ii
        double precision area,value
        double precision accum,area_tot

!-----------------------------------------------------------
! initialize arrays
!-----------------------------------------------------------

	accum = 0.
	area_tot = 0.

!-----------------------------------------------------------
! accumulate values
!-----------------------------------------------------------

        do ie=1,nel
          area = 4.*ev(10,ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
            value = val(k)
	    accum = accum + value * area
	    area_tot = area_tot + area
          end do
        end do

!-----------------------------------------------------------
! compute final value
!-----------------------------------------------------------

	if ( area_tot > 0. ) accum = accum / area_tot

	aver = accum

!-----------------------------------------------------------
! end of routine
!-----------------------------------------------------------

        end

!******************************************************************

