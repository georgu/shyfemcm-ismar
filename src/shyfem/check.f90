
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998-2019  Georg Umgiesser
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

! routines for various checks
!
! contents :
!
! subroutine test3d(iunit,nn)           test output for new variables
! subroutine check_all			checks arrays for sanity (shell)
! subroutine check_fem			checks arrays for sanity
! subroutine check_values		checks important variables
! subroutine tsmass(ts,z,nlvdi,tstot)   computes mass of T/S or any conc. ts
! subroutine debug_dry			writes debug information on dry areas
! subroutine debug_node(k)		writes debug information on node k
! subroutine mimafem(string)		writes some min/max values to stdout
! subroutine mass_conserve		checks mass conservation
!
! subroutine check_set_unit(iu)		sets unit for debug output
! subroutine check_get_unit(iu)		gets unit for debug output
! subroutine check_node(k)		debug info on node k
! subroutine check_elem(ie)		debug info on element ie
! subroutine check_nodes_in_elem(ie)	debug info on nodes in element ie
! subroutine check_elems_around_node(k) debug info on elements around node k
! subroutine check_nodes_around_node(k) debug info on nodes around node k
!
! revision log :
!
! 24.08.1998	ggu	levdbg used for debug
! 26.08.1998	ggu	subroutine tsmass transferred from newbcl0
! 26.08.1998	ggu	subroutine convol, tstvol transferred from newcon1
! 22.10.1999	ggu	igrv, ngrv eliminated (subst by ilinkv, lenkv)
! 07.03.2000	ggu	arrays for vertical turbulent diffusion coefficient
! 05.12.2001	ggu	always execute tstvol, more debug info
! 11.10.2002	ggu	call to setaix deleted
! 09.01.2003	ggu	some variables saved in contst
! 27.03.2003	ggu	new routine value_check
! 31.07.2003	ggu	eliminate compiler warnings
! 31.07.2003	ggu	eliminate useless variables
! 10.08.2003	ggu	new routine check_fem
! 03.09.2003	ggu	routines check and sanity_check deleted
! 03.09.2003	ggu	renamed value_check to check_values, new check_all
! 13.03.2004	ggu	write total volume to inf file
! 15.03.2005	ggu	call to check_austausch() eliminated
! 23.03.2006	ggu	changed time step to real
! 23.08.2007	ggu	test for boundary nodes using routines in testbndo.h
! 27.09.2007	ggu	deleted tstvol,tstvol1,contst, new mass_conserve
! 24.06.2008	ggu	bpresv deleted
! 06.12.2008	ggu	read vreps from STR file
! 21.01.2009	ggu	new var vrerr to stop if mass error is too high
! 23.03.2009	ggu	more debug for vrerr, new routine check_node()
! 02.04.2009	ggu	new routine check_elem()
! 06.04.2009	ggu	new check_elems_around_node, check_nodes_in_elem
! 26.02.2010	ggu	in test3d() write also meteo data
! 23.03.2010	ggu	changed v6.1.1
! 08.04.2010	ggu	more info in checks (depth and area)
! 20.12.2010	ggu	changed VERS_6_1_16
! 17.02.2011	ggu	changed VERS_6_1_18
! 17.05.2011	ggu	new routine check_set_unit() to set output unit
! 31.05.2011	ggu	changed VERS_6_1_23
! 12.07.2011	ggu	loop only over actual nodes/elements, not dimensions
! 15.07.2011	ggu	new routines for CRC computation
! 19.08.2011	ggu	changed VERS_6_1_29
! 25.10.2011	ggu	hlhv eliminated
! 04.11.2011	ggu	changed VERS_6_1_35
! 30.03.2012	ggu	changed VERS_6_1_51
! 05.12.2013	ggu	changed VERS_6_1_70
! 15.05.2014	ggu	write mass error only for levdbg >= 3
! 19.12.2014	ggu	changed VERS_7_0_10
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 10.07.2015	ggu	changed VERS_7_1_50
! 13.07.2015	ggu	changed VERS_7_1_51
! 17.07.2015	ggu	changed VERS_7_1_53
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 24.07.2015	ggu	changed VERS_7_1_82
! 17.09.2015	ggu	in mass_conserve aux variables are local
! 29.09.2015	ggu	changed VERS_7_2_5
! 16.12.2015	ggu	changed VERS_7_3_16
! 28.04.2016	ggu	changed VERS_7_5_9
! 14.06.2016	ggu	changed VERS_7_5_14
! 17.06.2016	ggu	changed VERS_7_5_15
! 05.12.2017	ggu	changed VERS_7_5_39
! 03.04.2018	ggu	changed VERS_7_5_43
! 19.04.2018	ggu	changed VERS_7_5_45
! 14.02.2019	ggu	changed VERS_7_5_56
! 16.02.2019	ggu	changed VERS_7_5_60
! 29.08.2020	ggu	added new routine check_nodes_around_node()
! 27.03.2021	ggu	some femtime.h eliminated (not all), cleanup
! 31.05.2021	ggu	write time line in node/elem debug
! 02.04.2023	ggu	only master writes to iuinfo
! 15.03.2024	ggu	new routine debug_write_var()
! 10.09.2024	ggu&lrp	bug fix: do not compute mass balance on domain boundary
! 10.09.2024    lrp     relax check_values
! 08.02.2025    ggu     more info in check_* routines
!
!*************************************************************

	subroutine test3d(iunit,nn)

! test output for new variables
!
! nn	number of first array elements to be printed

	use mod_meteo
	use mod_internal
	use mod_geom_dynamic
	use mod_depth
	use mod_ts
	use mod_diff_visc_fric
	use mod_hydro_vel
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer iunit
	integer nn

	logical bmeteo
	integer i,l,nk,ne
	integer iu,ii
	double precision dtime
	character*20 aline

	bmeteo = .false.

	iu = iunit
	if( iu .le. 0 ) iu = 6

	call get_act_dtime(dtime)
	call get_act_timeline(aline)

	write(iu,*) 'time:',dtime,'  ',aline
	write(iu,*) 'nn  :',nn
	write(iu,*) 'nkn :',nkn
	write(iu,*) 'nel :',nel
	write(iu,*) 'nlvdi,nlv :',nlvdi,nlv	

	write(iu,*) 'hlv :'
	write(iu,*) (hlv(l),l=1,nlv)
	write(iu,*) 'hldv :'
	write(iu,*) (hldv(l),l=1,nlv)

	if(nn.eq.0) then
		nk=nkn
		ne=nel
	else
		nk=min(nn,nkn)		!$$cmplerr
		ne=min(nn,nel)
	end if

	write(iu,*) 'ilhv :'
	write(iu,*) (ilhv(i),i=1,ne)
	write(iu,*) 'fcorv :'
	write(iu,*) (fcorv(i),i=1,ne)
	write(iu,*) 'hev :'
	write(iu,*) (hev(i),i=1,ne)
	write(iu,*) 'iwegv :'
	write(iu,*) (iwegv(i),i=1,ne)

	write(iu,*) 'zov :'
	write(iu,*) (zov(i),i=1,nk)
	write(iu,*) 'znv :'
	write(iu,*) (znv(i),i=1,nk)
	write(iu,*) 'zeov :'
	write(iu,*) ((zeov(ii,i),ii=1,3),i=1,ne)
	write(iu,*) 'zenv :'
	write(iu,*) ((zenv(ii,i),ii=1,3),i=1,ne)

	if( bmeteo ) then
	write(iu,*) 'ppv :'
	write(iu,*) (ppv(i),i=1,nk)
	write(iu,*) 'wxv :'
	write(iu,*) (wxv(i),i=1,nk)
	write(iu,*) 'wyv :'
	write(iu,*) (wyv(i),i=1,nk)
	write(iu,*) 'tauxnv :'
	write(iu,*) (tauxnv(i),i=1,nk)
	write(iu,*) 'tauynv :'
	write(iu,*) (tauynv(i),i=1,nk)
	end if

	do l=1,nlv
	write(iu,*)
	write(iu,*) 'level :',l
	write(iu,*) 'ulov :'
	write(iu,*) (ulov(l,i),i=1,ne)
	write(iu,*) 'vlov :'
	write(iu,*) (vlov(l,i),i=1,ne)
	write(iu,*) 'wlov :'
	write(iu,*) (wlov(l-1,i),i=1,nk)
	write(iu,*) 'ulnv :'
	write(iu,*) (ulnv(l,i),i=1,ne)
	write(iu,*) 'vlnv :'
	write(iu,*) (vlnv(l,i),i=1,ne)
	write(iu,*) 'wlnv :'
	write(iu,*) (wlnv(l-1,i),i=1,nk)
	write(iu,*) 'utlov :'
	write(iu,*) (utlov(l,i),i=1,ne)
	write(iu,*) 'vtlov :'
	write(iu,*) (vtlov(l,i),i=1,ne)
	write(iu,*) 'utlnv :'
	write(iu,*) (utlnv(l,i),i=1,ne)
	write(iu,*) 'vtlnv :'
	write(iu,*) (vtlnv(l,i),i=1,ne)
	write(iu,*) 'visv :'
	write(iu,*) (visv(l,i),i=1,nk)
	write(iu,*) 'difv :'
	write(iu,*) (difv(l,i),i=1,nk)
	write(iu,*) 'tempv :'
	write(iu,*) (tempv(l,i),i=1,nk)
	write(iu,*) 'saltv :'
	write(iu,*) (saltv(l,i),i=1,nk)
	write(iu,*) 'difhv :'
	write(iu,*) (difhv(l,i),i=1,ne)
	end do

	end

!******************************************************************

	subroutine check_all

! checks arrays for sanity

	implicit none

        integer levdbg
        real getpar

        levdbg = nint(getpar('levdbg'))

        if( levdbg .ge. 5 ) call check_fem
        if( levdbg .ge. 3 ) call check_values

	!call mimafem('panic')

	end

!******************************************************************

	subroutine check_fem

! checks arrays for sanity

	implicit none

!-------------------------------------------------------
! check geom arrays
!-------------------------------------------------------

	call check_ev
	call check_geom

!-------------------------------------------------------
! check vertical structure
!-------------------------------------------------------

	call check_vertical

!-------------------------------------------------------
! check various arrays
!-------------------------------------------------------

	call check_eddy
	call check_coriolis
	call check_chezy

!-------------------------------------------------------
! end of routine
!-------------------------------------------------------

	end

!******************************************************************

	subroutine check_values

! checks important variables

	use mod_layer_thickness
	use mod_ts
	use mod_hydro_baro
	use mod_hydro_vel
	use mod_hydro
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*16 text
	
	integer nocheck
	real, parameter :: zero = 0.
	real, parameter :: zmax = 10.
	real, parameter :: vmax = 10.
	real, parameter :: umax = 100000.
	real, parameter :: hmax = 100000.
	real, parameter :: smin = -1.
	real, parameter :: smax = 70.
	real, parameter :: tmin = -30.
	real, parameter :: tmax = 70.

	real getpar

        nocheck = nint(getpar('nocheck'))
	if( nocheck > 0 ) return

	text = '*** check_values'

	!call check1Dr(nkn,zov,-zmax,zmax,text,'zov')
	!call check1Dr(nkn,znv,-zmax,zmax,text,'znv')

	!call check2Dr(3,3,nel,zeov,-zmax,zmax,text,'zeov')
	!call check2Dr(3,3,nel,zenv,-zmax,zmax,text,'zenv')

	!call check1Dr(nel,unv,-umax,umax,text,'unv')
	!call check1Dr(nel,vnv,-umax,umax,text,'vnv')

	!call check2Dr(nlvdi,nlv,nel,utlnv,-umax,umax,text,'utlnv')
	!call check2Dr(nlvdi,nlv,nel,vtlnv,-umax,umax,text,'vtlnv')

	!call check2Dr(nlvdi,nlv,nel,ulnv,-vmax,vmax,text,'ulnv')
	!call check2Dr(nlvdi,nlv,nel,vlnv,-vmax,vmax,text,'vlnv')

	call check2Dr(nlvdi,nlv,nkn,tempv,tmin,tmax,text,'tempv')
	call check2Dr(nlvdi,nlv,nkn,saltv,smin,smax,text,'saltv')

	call check2Dr(nlvdi,nlv,nkn,hdknv,zero,hmax,text,'hdknv')
	call check2Dr(nlvdi,nlv,nkn,hdkov,zero,hmax,text,'hdkov')

	call check2Dr(nlvdi,nlv,nel,hdenv,zero,hmax,text,'hdenv')
	call check2Dr(nlvdi,nlv,nel,hdeov,zero,hmax,text,'hdeov')

	end

!**********************************************************************

        subroutine tsmass(ts,mode,nlvdi,tstot)

! computes mass of T/S or any concentration ts

        implicit none

        integer nlvdi          !dimension of levels
        real ts(nlvdi,1)       !concentration on nodes
!        real z(3,1)             !water level
	integer mode
        real tstot              !total computed mass of ts

	double precision scalcont

	if( mode .ne. 1 .and. mode .ne. -1 ) then
	  write(6,*) 'mode = ',mode
	  stop 'error stop tsmass: wrong value for mode'
	end if

	tstot = scalcont(mode,ts)

        end

!************************************************************

	subroutine debug_dry

! writes debug information on dry areas

	use mod_geom_dynamic
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer ie,iweg
	double precision dtime
	character*20 aline

	iweg = 0
	do ie=1,nel
	  if( iwegv(ie) .ne. 0 ) iweg = iweg + 1
	end do

	call get_act_dtime(dtime)
	call get_act_timeline(aline)

	write(6,*) 'drydry... ',dtime,'  ',aline

	end

!*************************************************************

	subroutine debug_node(k)

! writes debug information on final volume around node k (internal)

	use mod_geom_dynamic
	use mod_depth
	use mod_hydro_baro
	use mod_hydro_print
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin
	use mkonst
	use femtime

	implicit none

	integer k

! common

	integer ie,ii,kk,l,i
	integer ilevel
	integer iweg
	real flux,dzvol,avvol
	real diff,rdiff
	real aj,uv0,uv1
	real b,c
	real dt,az,azt,azpar

	real getpar
	integer ipext,ieext

	integer, save :: netot
	integer, save, allocatable :: kinf(:,:)

	integer kmem
	save kmem
	data kmem / 0 /

	if( k == 0 ) then
	  allocate(kinf(2,ngr))
	  kinf = 0
	end if

	if( k .ne. kmem ) then
	  netot = 0
          do ie=1,nel
            do ii=1,3
	      kk = nen3v(ii,ie)
	      if( kk .eq. k ) then
	        netot = netot + 1
	        if( netot .gt. ngr ) then
		  stop 'error stop debug_node: ngr'
	        end if
	        kinf(1,netot) = ie
	        kinf(2,netot) = ii
	      end if
	    end do
	  end do
	  kmem = k
	  write(6,*) 'new node for debug...'
	  write(6,*) k,ipext(k),netot
	  do i=1,netot
	    ie = kinf(1,i)
	    ii = kinf(2,i)
	    write(6,*) ie,ieext(ie),ii
	  end do
	end if

! compute inflow into column and volume of column

	call get_timestep(dt)
	call getaz(azpar)
	az = azpar
	azt = 1. - az

	flux = 0.	! flux into water column
	dzvol = 0.	! volume change due to water level change
	avvol = 0.	! average volume of water column
	iweg = 0

	do i=1,netot
	  ie = kinf(1,i)
	  ii = kinf(2,i)
          aj=ev(10,ie)
          ilevel=ilhv(ie)
          kk=nen3v(ii,ie)
	  if( kk .ne. k ) stop 'error stop debug_node: internal error'
          b=ev(ii+3,ie)
          c=ev(ii+6,ie)
          uv0=0.
          uv1=0.
          do l=ilevel,1,-1
            uv1=uv1+utlnv(l,ie)*b+vtlnv(l,ie)*c
            uv0=uv0+utlov(l,ie)*b+vtlov(l,ie)*c
          end do
          uv1=unv(ie)*b+vnv(ie)*c
          uv0=uov(ie)*b+vov(ie)*c
          flux  = flux  + dt*12.*aj*(uv0*azt+uv1*az)
          dzvol = dzvol + 4.*aj*(zenv(ii,ie)-zeov(ii,ie))
          avvol = avvol + 4.*aj*hev(ie)
	  if( iwegv(ie) .ne. 0 ) iweg = iweg + 1
        end do

	diff = abs(flux-dzvol)
	rdiff = 0.
	if( avvol .gt. 0. ) rdiff = diff/avvol

	write(6,*) 'debug... ',it,diff,rdiff,iweg

	end

!*************************************************************

        subroutine mimafem(string)

! writes some min/max values to stdout

	use mod_layer_thickness
	use mod_ts
	use mod_hydro_print
	use mod_hydro_vel
	use mod_hydro
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        character*(*) string

	integer ie,l,k
	real u,v,w,z,s,t,c
        real high
        real zmin,zmax,umin,umax,vmin,vmax
        real hknmax,hkomax,henmax,heomax
        real utomax,utnmax,vtomax,vtnmax
        real hlvmax,h1vmax
        real bprmax

	integer ipext,ieext

	double precision dtime
	character*20 aline

!-----------------------------------------------------
! initial check and write
!-----------------------------------------------------

        !return  !FIXME

	call get_act_dtime(dtime)
	call get_act_timeline(aline)

        write(6,*) '------------------ ',trim(string)
        write(6,*) '   time: ',dtime,'  ',aline

!-----------------------------------------------------
! check water levels and barotropic velocities
!-----------------------------------------------------

        high = 1.e+30

        zmin =  high
        zmax = -high
        umin =  high
        umax = -high
        vmin =  high
        vmax = -high

	do k=1,nkn
	  z = znv(k)
	  u = up0v(k)
	  v = vp0v(k)
          zmin = min(zmin,z)
          zmax = max(zmax,z)
          umin = min(umin,u)
          umax = max(umax,u)
          vmin = min(vmin,v)
          vmax = max(vmax,v)
	end do

        write(6,*) zmin,zmax,umin,umax,vmin,vmax

!-----------------------------------------------------
! check of layer thickness
!-----------------------------------------------------

        hknmax = -high
        hkomax = -high
        do k=1,nkn
          do l=1,nlv
            hknmax = max(hknmax,hdknv(l,k))
            hkomax = max(hkomax,hdkov(l,k))
          end do
        end do

        henmax = -high
        heomax = -high
        do ie=1,nel
          do l=1,nlv
            henmax = max(henmax,hdenv(l,ie))
            heomax = max(heomax,hdeov(l,ie))
          end do
        end do

        write(6,*) hknmax,hkomax,henmax,heomax

!-----------------------------------------------------
! check of transports
!-----------------------------------------------------

        utomax = -high
        utnmax = -high
        vtomax = -high
        vtnmax = -high
        do ie=1,nel
          do l=1,nlv
            utomax = max(utomax,abs(utlov(l,ie)))
            utnmax = max(utnmax,abs(utlnv(l,ie)))
            vtomax = max(vtomax,abs(vtlov(l,ie)))
            vtnmax = max(vtnmax,abs(vtlnv(l,ie)))
          end do
        end do

        write(6,*) utomax,utnmax,vtomax,vtnmax

!-----------------------------------------------------
! end of routine
!-----------------------------------------------------

        write(6,*) '------------------'

        end

!*************************************************************

	subroutine vol_mass(mode)

! computes and writes total water volume

	use shympi
	use mod_info_output

        implicit none

	integer mode

        real mtot              !total computed mass of ts
	real mtot_orig
	double precision dmtot
	double precision masscont
	character*20 aline

	integer, save :: iuinfo = 0

	if( .not. binfo ) return

	if( iuinfo == 0 ) then
	  iuinfo = -1
          if(shympi_is_master()) call getinfo(iuinfo)
	end if

	if( mode .ne. 1 .and. mode .ne. -1 ) then
	  write(6,*) 'mode = ',mode
	  stop 'error stop vol_mass: wrong value for mode'
	end if

	dmtot = masscont(mode)
	mtot_orig = real(dmtot)

	!dmtot = shympi_sum(dmtot)
	!mtot = real(dmtot)
        !if( iuinfo > 0 ) then
	!  call get_act_timeline(aline)
	!  write(iuinfo,*) 'total_volume: ',aline,mtot
	!end if
	call info_output('total_volume','sum',mtot_orig,.true.)

        end

!*************************************************************

	subroutine mass_conserve

! checks mass conservation of single boxes (finite volumes)

	use mod_bound_geom
	use mod_bound_dynamic
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin
	use shympi
	use mkonst
	use mod_info_output

	implicit none

	logical berror,bdebug
	integer ie,l,ii,k,lmin,lmax,mode,ks,kss,ie_mpi
	integer levdbg
	real am,az,azt,dt,azpar,ampar
	real areafv,b,c
	real ffn,ffo,ff
	real vmax,vrmax,vdiv,vdiff,vrdiff
	real abot,atop
	real volo,voln
	real ubar,vbar
	real vbmax,vlmax,vrbmax,vrlmax
	real verrvol(4)
	real vrwarn,vrerr
	real qinput
	character*20 aline
	double precision vtotmax,vvv,vvm
	real, allocatable :: vf(:,:)
	real, allocatable :: va(:,:)

	real volnode,areanode,getpar

	integer, save :: iuinfo = 0

!----------------------------------------------------------------
! initialize
!----------------------------------------------------------------

	if( iuinfo == 0 ) then
	  iuinfo = -1
          if(shympi_is_master()) call getinfo(iuinfo)
	end if

	vrwarn = getpar('vreps')
	vrerr = getpar('vrerr')
	levdbg = nint(getpar('levdbg'))

	if( levdbg .le. 1 ) return

	mode = +1
        call getazam(azpar,ampar)
	az = azpar
	am = ampar
        azt = 1. - az
	call get_timestep(dt)

	call get_act_timeline(aline)

	allocate(vf(nlvdi,nkn),va(nlvdi,nkn))
	vf = 0.
	va = 0.

!----------------------------------------------------------------
! compute horizontal divergence
!----------------------------------------------------------------

        do ie_mpi=1,nel
          ie = ip_sort_elem(ie_mpi)
          areafv = 4. * ev(10,ie)               !area of triangle / 3
	  lmin = jlhv(ie)
          lmax = ilhv(ie)
          do l=1,lmax
            do ii=1,3
                k=nen3v(ii,ie)
                b = ev(ii+3,ie)
                c = ev(ii+6,ie)
                ffn = utlnv(l,ie)*b + vtlnv(l,ie)*c
                ffo = utlov(l,ie)*b + vtlov(l,ie)*c
                ff = ffn * az + ffo * azt
                vf(l,k) = vf(l,k) + 3. * areafv * ff
                va(l,k) = va(l,k) + areafv
            end do
          end do
        end do

!----------------------------------------------------------------
! include vertical divergence
!----------------------------------------------------------------

	ks = 1000
	ks = 5071
	ks = 0
	if( ks .gt. 0 ) then
	  k = ks
	  lmin = jlhkv(k)
	  lmax = ilhkv(k)
	  write(77,*) '------------- mass_conserve'
	  write(77,*) k,lmin,lmax
	  write(77,*) (vf(l,k),l=lmin,lmax)
	  write(77,*) (wlnv(l,k),l=lmin,lmax)
	  vtotmax = 0.
	  do l=lmin,lmax
	    vtotmax = vtotmax + vf(l,k)
	  end do
	  write(77,*) 'from box: ',vtotmax
	end if

	vtotmax = 0.
	do k=1,nkn_inner
	  lmin = jlhkv(k)
          lmax = ilhkv(k)
	  abot = 0.
	  vvv = 0.
	  vvm = 0.
	  do l=lmax,1,-1
	    atop = va(l,k)
	    vdiv = wlnv(l,k)*abot - wlnv(l-1,k)*atop
	    vf(l,k) = vf(l,k) + vdiv + mfluxv(l,k)
	    abot = atop
	    vvv = vvv + vdiv
	    vvm = vvm + mfluxv(l,k)
	    if( k .eq. ks ) write(77,*) 'vdiv: ',l,vf(l,k),vdiv,vvv
	  end do
	  vtotmax = max(vtotmax,abs(vvv))
	  if( k .eq. ks ) write(77,*) 'vvv: ',vvv,vvm
	end do

!----------------------------------------------------------------
! check mass balance in boxes
!----------------------------------------------------------------

	kss = 5226
	kss = 0
	berror = .false.
	vrmax = 0.
	vmax = 0.
	do k=1,nkn_inner
	  if( is_zeta_boundary(k) ) cycle
	  if( is_external_boundary(k) ) cycle
	  bdebug = k .eq. kss
	  if( bdebug ) write(78,*) '============================='
	  berror = .false.
	  lmin = jlhkv(k)
          lmax = ilhkv(k)
	  do l=lmin,lmax
	    voln = volnode(l,k,+1)
	    volo = volnode(l,k,-1)
	    vdiv = vf(l,k)
	    vdiff = voln - volo - vdiv * dt
	    vdiff = abs(vdiff)
	    vrdiff = vdiff / volo
	    vmax = max(vmax,vdiff)
	    vrmax = max(vrmax,vrdiff)
	    if( bdebug ) then
	        write(78,*) l,k
	        write(78,*) volo,voln,vdiff,vrdiff
	        write(78,*) vdiv,vdiv*dt
	    end if
	    if( vrdiff .gt. vrerr ) then
		berror = .true.
	        write(6,*) 'mass_conserve: ',l,k
	        write(6,*) volo,voln,vdiff,vrdiff
	        write(6,*) vdiv,vdiv*dt
	    end if
	    if( k == kss .and. l == 1 ) then
	      write(182,*) aline,vdiff,vrdiff
	    end if
	  end do
	  if( berror ) call check_node(k)
	  if( bdebug ) then
		call check_set_unit(78)
		call check_node(k)
		write(78,*) '============================'
	  end if
	end do

	vlmax = vmax		!absolute error for each box
	vrlmax = vrmax		!relative error for each box

!----------------------------------------------------------------
! barotropic
!----------------------------------------------------------------

	do k=1,nkn
	    vf(1,k) = 0.
	    va(1,k) = 0.
	end do

        do ie_mpi=1,nel
          ie = ip_sort_elem(ie_mpi)

          areafv = 4. * ev(10,ie)               !area of triangle / 3

	  ubar = 0.
	  vbar = 0.
          lmax = ilhv(ie)
          do l=1,lmax
	    ubar = ubar + az * utlnv(l,ie) + azt * utlov(l,ie)
	    vbar = vbar + az * vtlnv(l,ie) + azt * vtlov(l,ie)
	  end do

          do ii=1,3
                k=nen3v(ii,ie)
                b = ev(ii+3,ie)
                c = ev(ii+6,ie)
		ff = ubar * b + vbar * c
                vf(1,k) = vf(1,k) + 3. * areafv * ff
                va(1,k) = va(1,k) + areafv
          end do
        end do

	if( ks .gt. 0 ) write(77,*) 'from baro: ',vf(1,ks)

	vrmax = 0.
	vmax = 0.
	do k=1,nkn_inner
	 !if( is_inner(k) ) then
	 if( .not. is_external_boundary(k) ) then
           lmax = ilhkv(k)
	   voln = 0.
	   volo = 0.
	   qinput = 0.
	   do l=1,lmax
	     voln = voln + volnode(l,k,+1)
	     volo = volo + volnode(l,k,-1)
	     qinput = qinput + mfluxv(l,k)
	   end do
	   vdiv = vf(1,k) + rqv(k)
	   vdiv = vf(1,k) + qinput	!should be the same
	   vdiff = voln - volo - vdiv * dt
	   if( k .eq. ks ) write(77,*) 'vdiff: ',vdiff
	   !if( vdiff .gt. 0.1 ) write(6,*) 'baro error: ',k,vdiff
	   vdiff = abs(vdiff)
	   vrdiff = vdiff / volo
	   vmax = max(vmax,vdiff)
	   vrmax = max(vrmax,vrdiff)
	 end if
	end do

	vbmax = vmax		!absolute error for water column
	vrbmax = vrmax		!relative error for water column

	verrvol = (/vlmax,vbmax,vrbmax,vrlmax/)

!----------------------------------------------------------------
! write diagnostic output
!----------------------------------------------------------------

!	vbmax 		!absolute error for water column
!	vrbmax 		!relative error for water column
!	vlmax 		!absolute error for each box
!	vrlmax 		!relative error for each box

	if( vrlmax .gt. vrwarn ) then
	  if( levdbg .ge. 3 ) then
	    write(6,*) 'mass error: ',vbmax,vlmax,vrbmax,vrlmax
	  end if
	  if( vrlmax .gt. vrerr ) then
	    write(6,*) 'mass error of matrix solution is very high'
	    write(6,*) 'the relative mass error is = ',vrlmax
	    write(6,*) 'the limit of the mass error is vrerr = ',vrerr
	    write(6,*) 'Probably there is some problem with the solution'
	    write(6,*) 'of the system matrix. However, if you think'
	    write(6,*) 'you can live with this mass error, then please'
	    write(6,*) 'increase the value of vrerr in the STR file.'
	    stop 'error stop mass_conserve: mass error too high'
	  end if
	end if

	verrvol = (/vbmax,vlmax,vrbmax,vrlmax/)
	call info_output('mass_balance','max',4,verrvol)

	vlmax = shympi_max(vlmax)
	vbmax = shympi_max(vbmax)
	vrbmax = shympi_max(vrbmax)
	vrlmax = shympi_max(vrlmax)

	!if( iuinfo > 0 ) then
	!  write(iuinfo,*) 'mass_balance: ',vbmax,vlmax,vrbmax,vrlmax
	!end if

	deallocate(vf,va)

!----------------------------------------------------------------
! end of routine
!----------------------------------------------------------------

	end

!*************************************************************
!*************************************************************
!************ CRC computation ********************************
!*************************************************************
!*************************************************************

	subroutine check_crc

	use mod_internal
	use mod_depth
	use mod_layer_thickness
	use mod_bound_dynamic
	use mod_area
	use mod_ts
	use mod_diff_visc_fric
	use mod_hydro_print
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin, only : nkn,nel,ngr,mbw
	use femtime

	implicit none

	integer icrc,iucrc
	save iucrc
	data iucrc /0/		! unit

	integer ifemop

	icrc = 0		! level of output [0-10]

	if( icrc .le. 0 ) return

	if( iucrc .le. 0 ) then
	  iucrc = ifemop('.crc','form','new')
	  if( iucrc .le. 0 ) stop 'error stop check_crc: open file'
	end if

	write(iucrc,*) '====================================='
	write(iucrc,*) ' crc check at time = ',it
	write(iucrc,*) '====================================='

	call check_crc_1d(iucrc,'znv',nkn,znv)
	call check_crc_1d(iucrc,'zenv',3*nel,zenv)
	call check_crc_2d(iucrc,'utlnv',nlvdi,nel,ilhv,utlnv)
	call check_crc_2d(iucrc,'vtlnv',nlvdi,nel,ilhv,vtlnv)
	call check_crc_2d(iucrc,'saltv',nlvdi,nkn,ilhkv,saltv)
	call check_crc_2d(iucrc,'tempv',nlvdi,nkn,ilhkv,tempv)
	call check_crc_2d(iucrc,'rhov',nlvdi,nkn,ilhkv,rhov)

	if( icrc .le. 1 ) return

	!call check_crc_1d(iucrc,'ev',evdim*nel,ev)	!FIXME - double
	call check_crc_1d(iucrc,'hev',nel,hev)
	call check_crc_1d(iucrc,'fcorv',nel,fcorv)
	call check_crc_2d(iucrc,'visv',nlvdi,nkn,ilhkv,visv)
	call check_crc_2d(iucrc,'difv',nlvdi,nkn,ilhkv,difv)
	call check_crc_2d(iucrc,'hdknv',nlvdi,nkn,ilhkv,hdknv)
	call check_crc_2d(iucrc,'hdenv',nlvdi,nel,ilhv,hdenv)

	if( icrc .le. 2 ) return

	call check_crc_2d(iucrc,'ulnv',nlvdi,nel,ilhv,ulnv)
	call check_crc_2d(iucrc,'vlnv',nlvdi,nel,ilhv,vlnv)
	call check_crc_2d(iucrc,'mfluxv',nlvdi,nkn,ilhkv,mfluxv)
	call check_crc_2d(iucrc,'areakv',nlvdi,nkn,ilhkv,areakv)
	call check_crc_2d(iucrc,'wlnv',nlvdi+1,nkn,ilhkv,wlnv)
	call check_crc_2d(iucrc,'wprv',nlvdi,nkn,ilhkv,wprv)

	if( icrc .le. 3 ) return

	end

!*************************************************************

	subroutine check_crc_2d(iu,text,nlvdi,n,levels,array)

	use mod_debug

	implicit none

	integer iu
	character*(*) text
	integer nlvdi
	integer n
	integer levels(n)
	real array(nlvdi,n)

	integer crc,nlv

	nlv = 0		! use levels

	call checksum_2d(nlvdi,n,nlv,levels,array,crc)
	write(iu,*) text,'   ',n,nlvdi,nlv,crc

	end

!*************************************************************

	subroutine check_crc_1d(iu,text,n,array)

	use mod_debug

	implicit none

	integer iu
	character*(*) text
	integer n
	real array(n)

	integer crc

	call checksum_1d(n,array,crc)
	write(iu,*) text,'   ',n,crc

	end

!*************************************************************
!*************************************************************
!*************************************************************
!*************************************************************
!*************************************************************

	module check_unit

	integer, save :: iucheck = 6

	end module check_unit

!*************************************************************

	subroutine check_set_unit(iu)

	use check_unit
	use shympi

	implicit none

	integer iu

	iucheck = iu
	if( iu /= 6 ) iucheck = iucheck + my_id

	end

!*************************************************************

	subroutine check_get_unit(iu)

	use check_unit

	implicit none

	integer iu

	iu = iucheck

	end

!*************************************************************
!*************************************************************
!*************************************************************
!*************************************************************
!*************************************************************

	subroutine check_node(k)

! writes debug information on node k

	use mod_geom_dynamic
	use mod_depth
	use mod_layer_thickness
	use mod_bound_dynamic
	use mod_area
	use mod_ts
	use mod_diff_visc_fric
	use mod_hydro_vel
	use mod_hydro
	use mod_meteo
	use levels
	use basin
	use mod_hydro_print
        use mod_nohyd
	use femtime
	use shympi

	implicit none

	integer k

	logical binner,buniq
	integer iu
	integer l,lmax,lmin,kk
	character*20 aline

	integer ipext
	real volnode

	call check_get_unit(iu)
        lmin = jlhkv(k)
	lmax = ilhkv(k)
	call get_act_timeline(aline)
	binner = shympi_is_inner_node(k)
	buniq = shympi_is_unique_node(k)

	write(iu,*) '-------------------------------- check_node'
	write(iu,*) 'time:            ',aline
	write(iu,*) 'my,id,inner,uniq:',my_id,id_node(k),binner,buniq
	write(iu,*) 'it,idt,k,kext:   ',it,idt,k,ipext(k)
	write(iu,*) 'lmin,lmax,inodv: ',lmin,lmax,inodv(k)
	write(iu,*) 'xgv,ygv:         ',xgv(k),ygv(k)
	write(iu,*) 'zov,znv:         ',zov(k),znv(k)
	write(iu,*) 'hkv,hkv+znv:     ',hkv(k),hkv(k)+znv(k)
	write(iu,*) 'evapv:           ',evapv(k)
	write(iu,*) 'hdkov:           ',(hdkov(l,k),l=1,lmax)
	write(iu,*) 'hdknv:           ',(hdknv(l,k),l=1,lmax)
	write(iu,*) 'areakv:          ',(areakv(l,k),l=1,lmax)
	write(iu,*) 'volold:          ',(volnode(l,k,-1),l=1,lmax)
	write(iu,*) 'volnew:          ',(volnode(l,k,+1),l=1,lmax)
	write(iu,*) 'wlnv:            ',(wlnv(l,k),l=0,lmax)
	write(iu,*) 'mfluxv:          ',(mfluxv(l,k),l=1,lmax)
	write(iu,*) 'tempv:           ',(tempv(l,k),l=1,lmax)
	write(iu,*) 'saltv:           ',(saltv(l,k),l=1,lmax)
	write(iu,*) 'visv:            ',(visv(l,k),l=0,lmax)
	write(iu,*) 'difv:            ',(difv(l,k),l=0,lmax)
	write(iu,*) 'qpnv:            ',(qpnv(l,k),l=1,lmax)
	write(iu,*) 'uprv:            ',(uprv(l,k),l=1,lmax)
	write(iu,*) 'vprv:            ',(vprv(l,k),l=1,lmax)
	write(iu,*) '-------------------------------------------'

	end

!*************************************************************

	subroutine check_elem(ie)

! writes debug information on element ie

	use mod_geom_dynamic
	use mod_depth
	use mod_layer_thickness
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin
	use femtime
	use shympi

	implicit none

	integer ie

	logical binner,buniq
	integer iu
	integer l,lmin,lmax,ii
	real zmed
	character*20 aline

	integer ieext

	call check_get_unit(iu)
        lmin = jlhv(ie)
	lmax = ilhv(ie)
	zmed = sum(zenv(:,ie))/3.
	call get_act_timeline(aline)
	binner = shympi_is_inner_elem(ie)
	buniq = shympi_is_unique_elem(ie)

	write(iu,*) '-------------------------------- check_elem'
	write(iu,*) 'time:             ',aline
	write(iu,*) 'my_id,inner,uniq: ',my_id,binner,buniq
	write(iu,*) 'id_elem:          ',id_elem(:,ie)
	write(iu,*) 'it,idt,ie,ieext:  ',it,idt,ie,ieext(ie)
	write(iu,*) 'lmin,lmax:        ',lmin,lmax
        write(iu,*) 'iwegv,iwetv:      ',iwegv(ie),iwetv(ie)
	write(iu,*) 'area:             ',ev(10,ie)*12.
	write(iu,*) 'nen3v  :          ',(nen3v(ii,ie),ii=1,3)
	write(iu,*) 'hev,hev+zenv:     ',hev(ie),hev(ie)+zmed
	write(iu,*) 'hm3v:             ',(hm3v(ii,ie),ii=1,3)
	write(iu,*) 'zeov:             ',(zeov(ii,ie),ii=1,3)
	write(iu,*) 'zenv:             ',(zenv(ii,ie),ii=1,3)
	write(iu,*) 'zov:              ',(zov(nen3v(ii,ie)),ii=1,3)
	write(iu,*) 'znv:              ',(znv(nen3v(ii,ie)),ii=1,3)
	write(iu,*) 'hdeov:            ',(hdeov(l,ie),l=1,lmax)
	write(iu,*) 'hdenv:            ',(hdenv(l,ie),l=1,lmax)
	write(iu,*) 'utlov:            ',(utlov(l,ie),l=1,lmax)
	write(iu,*) 'vtlov:            ',(vtlov(l,ie),l=1,lmax)
	write(iu,*) 'utlnv:            ',(utlnv(l,ie),l=1,lmax)
	write(iu,*) 'vtlnv:            ',(vtlnv(l,ie),l=1,lmax)
	write(iu,*) 'ulnv:             ',(ulnv(l,ie),l=1,lmax)
	write(iu,*) 'vlnv:             ',(vlnv(l,ie),l=1,lmax)
	write(iu,*) '-------------------------------------------'

	end

!*************************************************************

	subroutine check_nodes_in_elem(ie)

! writes debug information on nodes in element ie

	use basin

	implicit none

	integer ie

	integer ii,k,iu
	integer ieext

	call check_get_unit(iu)

	write(iu,*) '-------------------------------------------'
	write(iu,*) 'checking nodes in element: ',ie,ieext(ie)
	write(iu,*) '-------------------------------------------'

	do ii=1,3
	  k = nen3v(ii,ie)
	  call check_node(k)
	end do

	end

!*************************************************************

	subroutine check_elems_around_node(k)

! writes debug information on elements around node k

	use basin

	implicit none

	integer k

	integer ie,ii,kk,iu
	logical bdebug

	integer ipext

	call check_get_unit(iu)

	write(iu,*) '-------------------------------------------'
	write(iu,*) 'checking elements around node: ',k,ipext(k)
	write(iu,*) '-------------------------------------------'

	do ie=1,nel
	  bdebug = .false.
	  do ii=1,3
	    kk = nen3v(ii,ie)
	    if( kk .eq. k ) bdebug = .true.
	  end do
	  if( bdebug ) call check_elem(ie)
	end do

	end

!*************************************************************

	subroutine check_nodes_around_node(k)

! writes debug information on nodes around node k

	use basin

	implicit none

	integer k

	integer n,i,kk,iu
	integer nodes(ngr)

	integer ipext

	call check_get_unit(iu)

	write(iu,*) '-------------------------------------------'
	write(iu,*) 'checking nodes around node: ',k,ipext(k)
	write(iu,*) '-------------------------------------------'

	call get_nodes_around(k,ngr,n,nodes)

	do i=1,n
	  kk = nodes(i)
	  call check_node(kk)
	end do

	end

!*************************************************************

        subroutine debug_write_var

! this writes specific values on nodes and elements to file for comparison

	use mod_hydro
	use mod_hydro_vel
	use mod_ts
	use mod_diff_visc_fric
        use basin
        use levels
        use shympi

        implicit none

        integer nelg,nkng
	integer nktot,netot
        integer ie_int,ie_ext
        integer ik_int,ik_ext
        integer i,ie,ik,lmax
	integer iu,iubase,iumax
	integer ihour
        double precision dtime
        character*80 name

        integer, save, allocatable :: ies(:)
        integer, save, allocatable :: iks(:)
        integer, save, allocatable :: iue(:)
        integer, save, allocatable :: iuk(:)
        integer, save :: nie = 0
        integer, save :: nik = 0
        integer, save :: icall_local = 0

        integer ieint, ipint
        integer ieext, ipext

	return

        iubase = 700	! base unit
	iumax = 0
	nktot = 20	! how many nodes to write (0 for none)
	netot = 20	! how many elems to write (0 for none)

        nelg = nel_global
        nkng = nkn_global

        call get_act_dtime(dtime)

        if( icall_local == 0 ) then

          allocate(ies(nel))
          allocate(iks(nkn))
          allocate(iue(nel))
          allocate(iuk(nkn))

	  if( netot > 0 ) then
            do i=1,nelg,nelg/netot
	      ie_ext = ip_ext_elem(i)
              ie_int = ieint(ie_ext)
              if( ie_int > 0 .and. shympi_is_unique_elem(ie_int) ) then
                nie = nie + 1
                ies(nie) = ie_int
                iu = iubase + ie_ext
		iumax = max(iu,iumax)
                iue(nie) = iu
                call make_name_with_number('e',ie_ext,'aux',name)
                open(iu,file=name,status='unknown',form='formatted')
              end if
            end do
	  end if

	  if( nktot > 0 ) then
            do i=1,nkng,nkng/nktot
	      ik_ext = ip_ext_node(i)
              ik_int = ipint(ik_ext)
              if( ik_int > 0 .and. shympi_is_unique_node(ik_int) ) then
                nik = nik + 1
                iks(nik) = ik_int
                iu = iumax + ik_ext
                iuk(nik) = iu
                call make_name_with_number('n',ik_ext,'aux',name)
                open(iu,file=name,status='unknown',form='formatted')
              end if
            end do
	  end if

          icall_local = 1
        end if

	!write(6,*) 'writing debug to files: ',nie,iue(1),nik,iuk(1)

	ihour = dtime/3600
	if( 3600*ihour /= dtime ) return

        do i=1,nie
          iu = iue(i)
          ie = ies(i)
          ie_ext = ieext(ie)
          lmax = ilhv(ie)
          write(iu,*) dtime
          write(iu,*) ie_ext,lmax
          write(iu,*) 'ze: ',zenv(:,ie)
          write(iu,*) 'ut: ',utlnv(1:lmax,ie)
          write(iu,*) 'vt: ',vtlnv(1:lmax,ie)
	  flush(iu)
        end do

        do i=1,nik
          iu = iuk(i)
          ik = iks(i)
          ik_ext = ipext(ik)
          lmax = ilhkv(ik)
          write(iu,*) dtime
          write(iu,*) ik_ext,lmax
          write(iu,*) 'zn: ',znv(ik)
          write(iu,*) 'wn: ',wlnv(0:lmax,ik)
          write(iu,*) 'salt: ',saltv(1:lmax,ik)
          write(iu,*) 'temp: ',tempv(1:lmax,ik)
          write(iu,*) 'vis: ',visv(0:lmax,ik)
          write(iu,*) 'dif: ',difv(0:lmax,ik)
	  flush(iu)
        end do

	call shympi_barrier

        end

!*************************************************************

