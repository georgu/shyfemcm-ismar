
!--------------------------------------------------------------------------
!
!    Copyright (C) 2010-2012,2014-2019  Georg Umgiesser
!    Copyright (C) 2012,2016  Christian Ferrarin
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

! routines for stability computations
!
! revision log :
!
! 19.02.2010	ggu	new file to contain stability computations
! 26.02.2010	ggu	internal_stability restructured (on element)
! 08.03.2010	ggu	run only down to avail layers in stability (bug fix)
! 23.03.2010	ggu	changed v6.1.1
! 26.01.2011	ggu	robs for nudging implemented
! 16.02.2011	ggu	pass robs to subroutines, write stb-ind to nos file
! 20.05.2011	ggu	allow for elimination of elems due to high rstab
! 31.05.2011	ggu	changed VERS_6_1_23
! 01.06.2011	ggu	wsink for stability integrated
! 12.07.2011	ggu	new routine output_stability()
! 14.07.2011	ggu	new routine output_stability_node()
! 21.06.2012	ggu&ccf	variable vertical sinking velocity integrated
! 07.03.2014	ggu	changed VERS_6_1_72
! 08.04.2014	ggu	use rlin to determine advective stability
! 05.05.2014	ggu	changed VERS_6_1_74
! 26.11.2014	ggu	changed VERS_7_0_7
! 23.12.2014	ggu	changed VERS_7_0_11
! 09.01.2015	ggu	changed VERS_7_0_12
! 19.01.2015	ggu	changed VERS_7_1_3
! 20.05.2015	ggu	always compute stability, call to conzstab changed
! 05.06.2015	ggu	changed VERS_7_1_12
! 13.07.2015	ggu	changed VERS_7_1_51
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 23.09.2015	ggu	changed VERS_7_2_4
! 20.10.2015	ggu	in output_stability() bug that icall was not adjouned
! 23.10.2015	ggu	changed VERS_7_3_9
! 09.11.2015	ggu	changed VERS_7_3_13
! 25.05.2016	ggu	changed VERS_7_5_10
! 20.10.2016	ccf	pass rtauv for differential nudging
! 12.01.2017	ggu	changed VERS_7_5_21
! 05.12.2017	ggu	changed VERS_7_5_39
! 13.04.2018	ggu	re-structured, included gravity wave stability
! 16.04.2018	ggu	use also ilin to compute stability
! 16.02.2019	ggu	changed VERS_7_5_60
! 06.11.2019	ggu	femtime eliminated
! 30.03.2021	ggu	better error output
! 20.03.2022	ggu	upgraded to da_out
! 01.06.2022	ggu	in gravity_wave_stability() set hz to min 0
! 29.03.2023	ggu	exchange rindex,tindex,gindex, write to info file
! 02.04.2023    ggu     only master writes to iuinfo
! 09.05.2023    lrp     introduce top layer index variable
! 24.05.2023    ggu     debug section introduced in hydro_internal_stability()
! 01.12.2023    ggu     gindex not reduced in here
!
!*****************************************************************
!*****************************************************************
!*****************************************************************
!
! typical usage :
!
! call scalar_stability before advection of similar variables
!
! notes :
!
!	scal3sh							newcon
!		scalar_info_stability				newstab
!		scalar_stability				newstab
!			scalar_compute_stability		newstab
!				conzstab			newcon
!
!	set_timestep						subtime
!		hydro_stability					newstab
!			hydro_internal_stability		newstab
!				momentum_advective_stability	newexpl
!				momentum_viscous_stability	newexpl
!				gravity_wave_stability		newstab
!
! saux is never computed in conzstab - may be a bug		!FIXME
!
!*****************************************************************

!==================================================================
        module stab
!==================================================================

	implicit none

	integer, parameter :: ndim_stab = 50
	integer, save :: nentry = 0
	real, save :: rkind(2,ndim_stab) = 0.

!==================================================================
        end module stab
!==================================================================

	subroutine scalar_compute_stability(robs,rtauv,wsink,wsinkv &
     &					,rkpar,ai_ll,rindex,saux)

! computes stability index

	use mod_rungekutta, only : get_rungekutta_sdiag_weights, &
     &				   get_rungekutta_cdiag_weights, &
     &				   get_rungekutta_c_weights
	use mod_diff_visc_fric
	use levels, only : nlvdi,nlv
	use basin
	use shympi

	implicit none

	real robs
	real rtauv(nlvdi,nkn)
	real wsink
	real wsinkv(0:nlvdi,nkn)
        real rkpar
        real ai_ll
        real rindex
        real saux(nlvdi,nkn)

        real as_ll,ac_ll,c_l
        real difmol
	real ddt
        integer isact,istot

!----------------------------------------------------------------
! set parameters
!----------------------------------------------------------------

	call get_rungekutta_sdiag_weights(1,as_ll)
	call get_rungekutta_cdiag_weights(1,ac_ll)
	call get_rungekutta_c_weights(1,c_l)

        isact = 1
        istot = 1
        difmol = 0.
	ddt = 1.		!always for 1 sec
	rindex = 0.

!----------------------------------------------------------------
! call conzstab
!----------------------------------------------------------------

        call conzstab( &
     &          ddt,robs,rtauv,wsink,wsinkv,rkpar,difhv,difv &
     &		,difmol,ai_ll/c_l,as_ll/c_l,ac_ll/c_l &
     &          ,rindex,istot,isact,nlvdi,nlv)

!----------------------------------------------------------------
! propagate to all domains
!----------------------------------------------------------------

	rindex = shympi_max(rindex)

!----------------------------------------------------------------
! end of routine
!----------------------------------------------------------------

	end

!*****************************************************************

        subroutine scalar_basic_stability(dt,rkpar,rindex)

! computes scalar stability without nudging and sinking

	use mod_rungekutta, only : get_rungekutta_idiag_weights
	use levels, only : nlvdi,nlv
	use basin

        implicit none

	real dt
        real rkpar
        real rindex

        integer istot
	real ai_ll
	real robs
	real wsink
	real, allocatable :: wsinkv(:,:)
	real, allocatable :: rtauv(:,:)
	real, allocatable :: saux(:,:)

	allocate(wsinkv(0:nlvdi,nkn))
	allocate(rtauv(nlvdi,nkn))
	allocate(saux(nlvdi,nkn))

	wsinkv = 0.
	rtauv = 0.
	saux = 0.
	robs = 0.
	wsink = 0.

	call get_rungekutta_idiag_weights(1,ai_ll)
	call scalar_compute_stability(robs,rtauv,wsink,wsinkv,rkpar,ai_ll, &
     &					rindex,saux)

	end

!*****************************************************************

        subroutine scalar_stability(dt,robs,rtauv,wsink,wsinkv &
     &					,rkpar,rindex,istot,saux)

! gets stability index (if necessary computes it)

	use mod_rungekutta, only : get_rungekutta_idiag_weights
	use levels, only : nlvdi,nlv
	use basin

        implicit none

	real dt
	real robs
	real rtauv(nlvdi,nkn)
	real wsink
	real wsinkv(0:nlvdi,nkn)
        real rkpar
        real rindex
        integer istot
	real saux(nlvdi,nkn)

	real ai_ll

!----------------------------------------------------------------
! compute stability index
!----------------------------------------------------------------

	call get_rungekutta_idiag_weights(1,ai_ll)
	call scalar_compute_stability(robs,rtauv,wsink,wsinkv,rkpar,ai_ll, &
     &					rindex,saux)

!----------------------------------------------------------------
! scale to real time step dt
!----------------------------------------------------------------

	rindex = dt * rindex
	istot = 1 + rindex

!----------------------------------------------------------------
! end of routine
!----------------------------------------------------------------

        end

!*****************************************************************

        subroutine scalar_info_stability(dt,robs,rtauv,wsink,wsinkv &
     &					,rkpar,rindex,istot,saux)

! gets stability index (if necessary computes it)
!
! it is assumed that the program will exit after this call

	use mod_rungekutta, only : get_rungekutta_idiag_weights
	use levels, only : nlvdi,nlv
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        real dt
	real robs
	real rtauv(nlvdi,nkn)
	real wsink
	real wsinkv(0:nlvdi,nkn)
        real rkpar
        real rindex
        integer istot
	double precision dtime
	real saux(nlvdi,nkn)

	integer ia,id
	integer l,k
	real aindex
	real ai_ll

!----------------------------------------------------------------
! compute stability index
!----------------------------------------------------------------

	call get_rungekutta_idiag_weights(1,ai_ll)
	call scalar_compute_stability(robs,rtauv,wsink,wsinkv,rkpar,ai_ll, &
     &					rindex,saux)
	rindex = dt * rindex
	istot = 1 + rindex

!----------------------------------------------------------------
! write to terminal
!----------------------------------------------------------------

	ia = 1
	aindex = saux(1,1)

	do k=1,nkn
	  do l=1,nlv
	    if( saux(l,k) .gt. aindex ) then
	      aindex = saux(l,k)
	      ia = k
	    end if
	  end do
	end do

!ggu protect
	write(6,*) 'scalar_info_stability:'
	write(6,*) rkpar,ai_ll,rindex,istot
	write(6,*) ia,aindex,dt*aindex
	id = 0
	call get_act_dtime(dtime)
	call shy_write_scalar(id,'sta',dtime,1,778,nlvdi,saux)
!ggu protect

!----------------------------------------------------------------
! end of routine
!----------------------------------------------------------------

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

        subroutine hydro_stability(dt,rindex)

! computes stability index for hydro timestep - no error output

        implicit none

        real dt
        real rindex

	call hydro_internal_stability(0,dt,rindex)

	end

!**********************************************************************

        subroutine error_stability(dt,rindex)

! computes stability index for hydro timestep - with error output
!
! after this call the program should abort

        implicit none

        real dt
        real rindex

	call hydro_internal_stability(1,dt,rindex)

	end

!**********************************************************************

        subroutine eliminate_stability(rmax)

! eliminates elements with stability index higher than rmax

        implicit none

        real rmax

	integer mode
        real rindex,dt

	mode = 2
	dt = 0.
	rindex = rmax
	call hydro_internal_stability(mode,dt,rindex)

	end

!**********************************************************************

        subroutine hydro_internal_stability(mode,dt,rindex)

! computes stability index for hydro timestep (internal)
!
! mode = 0		normal call, compute stability
! mode = 1		error call, compute stability and write error message
! mode = 2		eliminate elements with r>rindex

	use levels
	use basin, only : nkn,nel,ngr,mbw
	use shympi
	use mod_info_output

        implicit none

	integer mode		!0: normal call  1:error output
        real dt			!time step to be used
        real rindex		!stability index (return)

	logical bdebug
	integer ie,l,lmax,lmin,iweg,ilin,ibarcl,iu
	integer, save :: iuinfo = 0
	integer, save :: icall = 0
        real rkpar,ahpar,rlin
	real dindex,aindex,tindex,sindex,gindex
	real rmax
	real array(4)

	logical openmp_in_parallel

	real, allocatable :: sauxe1(:,:)
	real, allocatable :: sauxe2(:,:)
	real, allocatable :: sauxe3(:)
	real, allocatable :: sauxe(:,:)

	real getpar
	logical is_i_nan
  
	if( iuinfo == 0 ) then
          iuinfo = -1
          if(shympi_is_master()) call getinfo(iuinfo)
        end if

        rkpar = 0.
	ahpar = getpar('ahpar')
	ibarcl = nint(getpar('ibarcl'))
	rlin = getpar('rlin')
	ilin = nint(getpar('ilin'))
	if( ibarcl == 1 ) ilin = 0	!for baroclinic use advection stability
	rlin = rlin * (1-ilin)

	allocate(sauxe1(nlvdi,nel),sauxe2(nlvdi,nel))
	allocate(sauxe(nlvdi,nel),sauxe3(nel))
	sauxe1 = 0.
	sauxe2 = 0.

	rmax = 1.e+30
	if( mode .eq. 2 ) rmax = rindex
	if( mode .eq. 2 ) then
		write(6,*) 'eliminating rmax: ',rmax
	end if

	call momentum_advective_stability(rlin,aindex,sauxe1)
	call momentum_viscous_stability(ahpar,dindex,sauxe2)
	call gravity_wave_stability(gindex,sauxe3)

	do ie=1,nel
	  sauxe(:,ie) = sauxe1(:,ie) + sauxe2(:,ie) + sauxe3(ie)
	end do

	if( .not. openmp_in_parallel() ) then
          call output_stability(dt,sauxe)	!in case write to file
	end if

	tindex = 0.
	do ie=1,nel
	  lmax = ilhv(ie)
	  lmin = jlhv(ie)
	  iweg = 0
	  do l=lmin,lmax
	    sindex = sauxe(l,ie)
	    if( sindex .ge. rmax ) iweg = 1
	    tindex = max(tindex,sindex)
	  end do
	  if( iweg .gt. 0 ) then
	    write(6,*) 'eliminating element for stability: ',ie
	    write(569,*) 'eliminating element for stability: ',ie
	    call check_set_unit(570)
	    call check_elem(ie)
	    call set_element_dry(ie)
	  end if
	end do

	tindex = tindex*dt
	aindex = aindex*dt
	dindex = dindex*dt
	gindex = gindex*dt

	array = (/tindex,aindex,dindex,gindex/)
	call shympi_array_reduce('max',array)
	call info_output('rindex','none',4,array,.false.)

	if( mode .eq. 1 ) then		!error output
	  write(6,*) 'hydro_internal_stability: '
	  write(6,*) 'tindex,aindex,dindex,gindex: '
	  write(6,*) array
	  call output_errout_stability(dt,sauxe)
	end if

	rindex = array(1)

	deallocate(sauxe1,sauxe2,sauxe3,sauxe)

	!if( iuinfo > 0 ) then
	!  write(iuinfo,*) 'rindex: ',array
	!end if

	!array = (/tindex,aindex,dindex,gindex/)
	!call info_output('rindex','max',4,array,.false.)

	icall = icall + 1

        end

!*****************************************************************
!*****************************************************************
!*****************************************************************

        subroutine output_errout_stability(dt,sauxe)

! outputs stability index for hydro timestep (internal) (error handling)

	use levels
	use basin, only : nkn,nel,ngr,mbw

        implicit none

        real dt
	real sauxe(nlvdi,nel)

	logical bnos
	integer ie,l,lmax
	integer iemax,iee
	real tindex
	real sauxn(nlvdi,nkn)

	integer icall,iustab,ifnos
	save icall,iustab,ifnos
	data icall,iustab,ifnos /0,0,0/

	integer ieext

	icall = icall + 1

	iemax = 0
	tindex = 0.

	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
	    if( sauxe(l,ie) .gt. tindex ) then
	      tindex = sauxe(l,ie)
	      iemax = ie
	    end if
	  end do
	end do

	ie = iemax
	iee = ieext(ie)
	write(6,*) 'errout_stability: '
	write(6,*) '  ie-intern   ie-extern' // &
     &			'   stab-index    dt*stab-index'
	write(6,*) ie,iee,tindex,tindex*dt

	end

!*****************************************************************

        subroutine output_stability(dt,sauxe)

! outputs stability index for hydro timestep (internal)

	use levels
	use basin

        implicit none

        real dt
	real sauxe(nlvdi,nel)	!stability index for element

	real, save, allocatable :: smax(:)

	logical bnos
	integer ie,ii,k,l,lmax
	integer ia,id,idc,nvar
	real sindex,smin
	real sx,sn
	real dtorig
	double precision dtime
	logical next_output_d,has_output_d,is_over_output_d

	integer, save :: icall = 0
	double precision, save :: da_out(4)

	if( icall .lt. 0 ) return

	if( icall .eq. 0 ) then
	  call init_output_d('itmsti','idtsti',da_out)
	  call increase_output_d(da_out)
	  if( .not. has_output_d(da_out) ) icall = -1
	  if( icall .lt. 0 ) return
	  nvar = 1
          call shyfem_init_scalar_file('sti',nvar,.true.,id)
          da_out(4) = id
	  allocate(smax(nkn))
	  smax = 0.
	end if

	icall = icall + 1

	if( .not. is_over_output_d(da_out) ) return 

	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
	    sindex = sauxe(l,ie)
	    do ii=1,3
	      k = nen3v(ii,ie)
	      smax(k) = max(smax(k),sindex)
	    end do
	  end do
	end do

	if( next_output_d(da_out) ) then
	  call get_orig_timestep(dtorig)
	  do k=1,nkn				!convert to time step
	    if( smax(k) > 0 ) then
	      smax(k) = 1./smax(k)
	      if( smax(k) > dtorig ) smax(k) = dtorig
	    else
	      smax(k) = dtorig
	    end if 
	  end do
	  call get_act_dtime(dtime)
	  id = nint(da_out(4))
          idc = 95
          call shy_write_scalar_record2d(id,dtime,idc,smax)
	  smax = 0.
	end if

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine gravity_wave_stability(gindex,garray)

! if the time discretization is explicit:
! outputs stability index related to gravity waves (internal).
! if the time discretization is implicit:
! the time discretization is unconditionally stable and you can
! jump this subroutine.
! Note: the selection of explicit or implicit discretization
! is based on the coefficients of the first stage.

	use mod_rungekutta, only : get_rungekutta_idiag_weights, &
     &				   get_rungekutta_c_weights
	use basin
	use mod_hydro
	use shympi
	use pkonst

	implicit none

	real gindex
	real garray(nel)

	integer ie,ii,ii1,k1,k2
	real distmin,d,dx,dy
	real ai_ll,c_l
	real hz,ri

	integer, save :: icall = 0
	real,save,allocatable :: dist(:)
        real, parameter :: high = 1.e+30

	real getpar

	gindex = 0.
	garray = 0.

	if( icall < 0 ) return

	if( icall == 0 ) then
	  call get_rungekutta_idiag_weights(1,ai_ll)
	  call get_rungekutta_c_weights(1,c_l)
	  if ( ai_ll/c_l >= 0.5 ) then
	    icall = -1
	  else
	    !ok
	  end if
	  if( icall < 0 ) return
	  icall = 1

	  allocate(dist(nel))
          do ie=1,nel
           distmin = high
           do ii=1,3
            ii1 = mod(ii,3) + 1
            k1 = nen3v(ii,ie)
            k2 = nen3v(ii1,ie)
            call compute_distance(xgv(k1),ygv(k1),xgv(k2),ygv(k2),dx,dy)
            d = dx**2 + dy**2
            distmin = min(distmin,d)
	   end do
	   dist(ie) = sqrt(distmin)
          end do
	end if

	gindex = 0.
	do ie=1,nel
	  hz = maxval( hm3v(:,ie) + zenv(:,ie) )
	  hz = maxval( hm3v(:,ie) )
	  hz = max(hz,0.)
	  ri = sqrt(grav*hz) / dist(ie)
	  gindex = max(gindex,ri)
	  garray(ie) = ri
	end do

	!gindex = shympi_max(gindex)

	!write(6,*) nel,gindex,1./gindex
	!stop

	return

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine parallel_test

! tests parallel implementation

	use levels, only : nlvdi,nlv
	use basin

	implicit none

	real dt,rkpar,ai_ll,rindex
	real robs,wsink
	real rtauv(nlvdi,nkn)
	real wsinkv(0:nlvdi,nkn)
	real saux(nlvdi,nkn)

	ai_ll = 0.
	rkpar = 0.
	robs = 0.
	wsink = 0.
	rtauv = 0.

	write(6,*) 'parallel test...'
	call scalar_compute_stability(robs,rtauv,wsink,wsinkv,rkpar,ai_ll, &
     &					rindex,saux)
	write(6,*) 'parallel is ok.'

	end

!*****************************************************************

