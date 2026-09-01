
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998-2000,2002-2003,2006-2007,2006-2007  Georg Umgiesser
!    Copyright (C) 2009-2020  Georg Umgiesser
!    Copyright (C) 2017  Christian Ferrarin
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

! subroutines for computing discharge / flux
!
! contents :
!
! subroutine inflxa
! subroutine rdflxa
! subroutine ckflxa
! subroutine prflxa
! subroutine tsflxa
!
! subroutine wrflxa				write of flux data
!
! subroutine flxscs(n,kflux,iflux,ai_ll,fluxes)	flux through sections
! subroutine flxsec(n,kflux,iflux,ai_ll,fluxes)	flux through section
!
! subroutine flxini				initializes flux routines
! subroutine flx_init(kfluxm,kflux,nsect,iflux)	sets up array iflux
! subroutine flxinf(m,kflux,iflux)		sets up one info structure
! function igtnsc(k1,k2)			gets number of internal section
!
! revision log :
!
! 30.04.1998	ggu	newly written routines (subpor deleted)
! 07.05.1998	ggu	check nrdveci on return for error
! 08.05.1998	ggu	restructured with new comodity routines
! 13.09.1999	ggu	type of node computed in own routine flxtype
! 19.11.1999	ggu	iskadj into sublin
! 20.01.2000	ggu	old routines substituted, new routine extrsect
! 20.01.2000	ggu	common block /dimdim/ eliminated
! 20.01.2000	ggu	common block /dimdim/ eliminated
! 01.09.2002	ggu	ggu99 -> bug in flx routines (how to reproduce?)
! 26.05.2003	ggu	in flxnov substituted a,b with b,c
! 26.05.2003	ggu	new routine make_fluxes (for lagrangian)
! 10.08.2003	ggu	do not call setweg, setnod, setkan
! 23.03.2006	ggu	changed time step to real
! 28.09.2007	ggu	use testbndo to determine boundary node in flxtype
! 28.04.2009	ggu	links re-structured
! 23.03.2010	ggu	changed v6.1.1
! 23.02.2011	ggu	new routine call write_node_fluxes() for special output
! 01.03.2011	ggu	changed VERS_6_1_20
! 01.06.2011	ggu	documentation to flxscs() changed
! 14.07.2011	ggu	changed VERS_6_1_27
! 21.09.2011	ggu	some lower-level subroutines copied to subflx.f
! 07.10.2011	ggu	adjusted for 3d flux routines
! 18.10.2011	ggu	changed VERS_6_1_33
! 19.10.2011	ggu	added T/S variables, created fluxes_*() routines
! 19.10.2011	ggu	added conz variables, created fluxes_template()
! 09.12.2011	ggu	changed VERS_6_1_38
! 01.06.2012	ggu	changed VERS_6_1_53
! 10.05.2013	ggu	introduced subflxa.h, common routines to subflxu.f
! 25.10.2013	ggu	changed VERS_6_1_68
! 07.03.2014	ggu	changed VERS_6_1_72
! 26.11.2014	ggu	changed VERS_7_0_7
! 19.12.2014	ggu	changed VERS_7_0_10
! 19.01.2015	ggu	changed VERS_7_1_3
! 20.05.2015	ggu	modules introduced
! 10.07.2015	ggu	changed VERS_7_1_50
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 09.11.2015	ggu	changed VERS_7_3_13
! 12.04.2016	ggu	fluxes_template adjourned
! 15.04.2016	ggu	fluxes_template debugged and finished
! 09.09.2016	ggu	changed VERS_7_5_17
! 31.03.2017	ggu	changed VERS_7_5_24
! 22.09.2017	ccf	added total sediment concentration
! 26.10.2017	ggu	reads itable, chflx and section description
! 04.11.2017	ggu	changed VERS_7_5_34
! 14.11.2017	ggu	changed VERS_7_5_36
! 17.11.2017	ggu	sediment output adapted to new framework
! 17.11.2017	ggu	changed VERS_7_5_38
! 27.03.2018	ggu	new code for generic flux handling (fluxes_generic())
! 03.04.2018	ggu	changed VERS_7_5_43
! 06.07.2018	ggu	changed VERS_7_5_48
! 16.02.2019	ggu	changed VERS_7_5_60
! 06.03.2020	ggu	new flux0d, get_barotropic_flux()
! 27.05.2022	ggu	changed to be used with mpi, fluxes now double
! 30.05.2022	ggu	more changes for mpi
! 18.05.2023	ggu	in flx_write() call flx_collect_3d()
! 22.05.2023	ggu	need fluxes_r for write
! 23.10.2024    ggu     module definition taken out to mod_flux.f90
! 10.12.2024    ggu     only reduce one time, then write
!
! notes :
!
! These routines can also be used internally to compute the flux
! over various sections. The following calling sequence must be respected:
!
! call flx_init(kfluxm,kflux,nsect,iflux)		initializes iflux
!
! call flxscs(kfluxm,kflux,iflux,ai_ll,fluxes) computes fluxes
!
! Initialization can be done anytime.
!
! nsect		total number of sections
! kfluxm	total number of nodes defining sections
! kflux()	node numbers defining sections
!
! calling sequence :
!
! rdflxa
!	flx_read_section
! ckflxa
!	convert_nodes
! prflxa
! tsflxa
!
! flxini
!
! wrflxa	administration
!	init:
!		flx_alloc_arrays
!		flux_initialize
!			flx_init
!				klineck
!				flxinf
!					igtnsc
!			get_nlayers
!			correct_nlayers
!			correct_iflux
!		fluxes_init_d
!		flx_file_open
!			flux_file_open
!				flx_write_header
!				flx_write_header2
!	loop:
!		flxscs
!			flxsec
!				flx3d
!		fluxes_accum_d
!		fluxes_aver_d
!		flx_write
!			flux_write
!				flx_write_record
!		fluxes_init_d
!
!
! there are 3 fluxes given for every section: fluxes(3,ns)
!
!	fluxes(1,is)		total flux
!	fluxes(2,is)		positive flux
!	fluxes(3,is)		negative flux
!
! in three 3D the definition is: fluxes(0:nlv,3,ns)
!
! in fluxes(0,:,:) the verticall integrated flux is stored
!
!******************************************************************
!******************************************************************
!******************************************************************

        subroutine flx_read_section(n,ns)

        use mod_flux
        use nls

        integer n,ns

	call nls_init_section

        !n = nls_read_vector()
	call nls_read_isctable(n,ns)
        kfluxm = n
	nsect = ns

        if( n > 0 ) then
          allocate(kflux(n))
          allocate(kflux_ext(n))
          allocate(iflux(3,n))
          allocate(itable(2,ns))
          allocate(chflx(ns))
	  chflx = ' '
	  call nls_copy_isctable(n,ns,kflux,itable,chflx)
	  kflux_ext = kflux
          !call nls_copy_int_vect(n,kflux)
        end if

	call nls_finish_section

        end subroutine flx_read_section

!******************************************************************

        subroutine flx_alloc_arrays0(nl,ns)

	use mod_flux

	implicit none

	integer nl	!layers
	integer ns	!sections

	if( nl == nl_flux .and. ns == ns_flux ) return

	if( nl > 0 .or. ns > 0 ) then
	  if( nl == 0 .or. ns == 0 ) then
            write(6,*) 'nl,ns: ',nl,ns
            stop 'error stop flx_alloc_arrays: incompatible parameters'
	  end if
	end if

	!write(6,*) 'flx_alloc_arrays: ',nl,ns

	if( ns_flux > 0 ) then
          deallocate(nlayers)
          deallocate(fluxes)
          deallocate(fluxes_r)
          deallocate(flux0d)
          deallocate(masst)
          deallocate(saltt)
          deallocate(tempt)
          deallocate(conzt)
          deallocate(ssctt)
	end if

	nl_flux = nl
	ns_flux = ns

	if( ns == 0 ) return

        allocate(nlayers(ns))
        allocate(fluxes(0:nl,3,ns))
        allocate(fluxes_r(0:nl,3,ns))
        allocate(flux0d(ns))

        allocate(masst(0:nl,3,ns))
        allocate(saltt(0:nl,3,ns))
        allocate(tempt(0:nl,3,ns))
        allocate(conzt(0:nl,3,ns))
        allocate(ssctt(0:nl,3,ns))

	flux0d = 0.

	bflxalloc = .true.

	end

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine mod_flx(mode)
 
	use befor_after

        implicit none
 
        integer mode
 
        if( mode .eq. M_AFTER ) then
           call wrflxa
        else if( mode .eq. M_INIT ) then
           ! nothing
        else if( mode .eq. M_READ ) then
           call rdflxa
        else if( mode .eq. M_CHECK ) then
           call ckflxa
        else if( mode .eq. M_SETUP ) then
           ! nothing
        else if( mode .eq. M_PRINT ) then
           call prflxa
        else if( mode .eq. M_TEST ) then
           call tsflxa
        else if( mode .eq. M_BEFOR ) then
           ! nothing
        else
           write(6,*) 'unknown mode : ', mode
           stop 'error stop mod_flx'
        end if
 
        end

!******************************************************************

        subroutine rdflxa

	use mod_flux
	use mod_info_output

        implicit none

	integer n,ns

        call flx_read_section(n,ns)

	if( print_not_quiet_once() ) then
	write(6,*) 'running rdflxa: ',n,ns
	end if

        if( n .lt. 0 ) then
          write(6,*) 'read error in section $flux'
          stop 'error stop rdflxa'
        end if

        end

!******************************************************************

        subroutine ckflxa

! converts external to internal nodes

	use mod_flux
	use mod_info_output

        implicit none

	if( kfluxm <= 0 ) return

	if( print_not_quiet_once() ) then
	write(6,*) 'running ckflxa: '
	end if

	call convert_nodes(kfluxm,kflux)

        end

!******************************************************************

	subroutine prflxa

	use mod_flux
	use mod_info_output

	implicit none

	integer nnode,ifirst,ilast
	integer ntotal,ns
	integer i,ii

	integer ipext
	logical nextline

	if( kfluxm == 0 ) return

	if( print_not_quiet_once() ) then

	write(6,*)
	write(6,*) 'flux section :'
	write(6,*)
	write(6,*) 'nsect,kfluxm ',nsect,kfluxm
	write(6,*)

	ns = 0
	nnode = 0

	do while( nextline(kflux,kfluxm,nnode,ifirst,ilast) )
	  ns = ns + 1
	  ntotal = ilast - ifirst + 1
	  write(6,*) 'section : ',ns,ntotal,'  ',trim(chflx(ns))
	  write(6,'(12i6)') (ipext(kflux(i)),i=ifirst,ilast)
	end do

	end if

	end

!******************************************************************

	subroutine tsflxa

	use mod_flux

	implicit none

	integer i,ii

	write(6,*) '/kfluxc/'
	write(6,*) nsect,kfluxm
	write(6,*) (kflux(i),i=1,kfluxm)

	write(6,*) '/iflux/'
	write(6,*) ((iflux(ii,i),ii=1,3),i=1,kfluxm)

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine get_barotropic_flux(is,flux0)

	use mod_flux

	implicit none

	integer is
	real flux0

	if( bflxalloc ) then
	  flux0 = flux0d(is)
	else
	  flux0 = 0.
	end if

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine wrflxa

! administers writing of flux data

	use mod_rungekutta, only : get_rungekutta_idiag_weights
	use mod_conz, only : cnv
	use mod_ts
	use mod_sediment, only : tcn
	use levels, only : nlvdi,nlv
	use mod_flux
	use shympi
	use mod_trace_point
	use mod_info_output

	implicit none

	integer j,i,l,lmax,ivar,nvers,nsaux
	integer idtflx,ierr,iv,is,n,ns
	integer ip,ipstart,ipend,np
	integer nl,nlg
	integer nbflx
	integer iunit6
	real ai_ll,dt
	double precision atime0,atime

	integer ifemop,ipext
	real getpar
	double precision dgetpar
	logical has_output_d,next_output_d,is_over_output_d,is_last_output_d

        double precision, allocatable, save :: array(:)
        double precision, allocatable, save :: flux_global(:,:,:)
        double precision, allocatable, save :: fluxes_global(:,:)
        integer, allocatable, save :: ivinfo(:,:)
        integer, allocatable, save :: nlayersg(:)

        double precision, save :: trm,trs,trt,trc,trsc
	logical, save :: btemp,bsalt,bconz,bsedi
	integer, save :: nvar
	integer, save :: icall = 0

!-----------------------------------------------------------------
! start of code
!-----------------------------------------------------------------

        if( icall .eq. -1 ) return

!-----------------------------------------------------------------
! initialization
!-----------------------------------------------------------------

        if( icall .eq. 0 ) then
		call trace_point('initializing flux')

		btemp = ( nint(getpar('itemp')) > 0 )
		bsalt = ( nint(getpar('isalt')) > 0 )
		bconz = ( nint(getpar('iconz')) == 1 )
		bsedi = ( nint(getpar('isedi')) > 0 )

		nvar = 1
		if( btemp ) nvar = nvar + 1
		if( bsalt ) nvar = nvar + 1
		if( bconz ) nvar = nvar + 1
		if( bsedi ) nvar = nvar + 1

		call init_output_d('itmflx','idtflx',da_out)
		call increase_output_d(da_out)
                if( .not. has_output_d(da_out) ) icall = -1

                if( kfluxm .le. 0 ) icall = -1
                if( nsect .le. 0 ) icall = -1
                if( icall .eq. -1 ) return

		if( print_not_quiet_once() ) then
		write(6,*) 'initializing flux sections'
		end if
		bflxinit = .true.

       		call flx_alloc_arrays(nlvdi,nsect)
		call flux_initialize(kfluxm,kflux,iflux,nsect,nlayers,nlmax)

		nlayers_global = nlayers
		call shympi_array_reduce('max',nlayers_global)
		if( any(nlayers_global /= nlayers) ) then
		  write(6,'(10i5)') nlmax,nlayers
		  write(6,'(10i5)') nlmax,nlayers_global
		  stop 'error stop wrflxa: nlayers_global /= nlayers'
		end if

		call fluxes_init_d(nlvdi,nsect,nlayers,trm,masst)
		if( bsalt ) then
		  call fluxes_init_d(nlvdi,nsect,nlayers,trs,saltt)
		end if
		if( btemp ) then
		  call fluxes_init_d(nlvdi,nsect,nlayers,trt,tempt)
		end if
		if( bconz ) then
		  call fluxes_init_d(nlvdi,nsect,nlayers,trc,conzt)
		end if
		if( bsedi ) then
		  call fluxes_init_d(nlvdi,nsect,nlayers,trsc,ssctt)
		end if

		call flx_file_open(nvar)

		allocate(ivinfo(3,nvar))
		allocate(array((nlv_global+1)*3*ns_flux*nvar))
		allocate(flux_global(0:nlv_global,3*ns_flux,nvar))
		allocate(fluxes_global(0:nlv_global,3*ns_flux))
		ivinfo = 0
		array = 0.
		flux_global = 0.
		fluxes_global = 0.

!               here we could also compute and write section in m**2

		icall = 1
		return

        end if

!-----------------------------------------------------------------
! normal call
!-----------------------------------------------------------------

        if( .not. is_over_output_d(da_out) ) return

	call get_timestep(dt)
	call get_rungekutta_idiag_weights(1,ai_ll)

!	-------------------------------------------------------
!	accumulate results
!	-------------------------------------------------------

	ivar = 0
	call flxscs(kfluxm,kflux,iflux,ai_ll,fluxes,ivar,rhov)
	call fluxes_accum_d(nlvdi,nsect,nlayers,dt,trm,masst,fluxes)

	flux0d(:) = fluxes(0,1,:)		!remember barotropic fluxes

	if( bsalt ) then
	  ivar = 11
	  call flxscs(kfluxm,kflux,iflux,ai_ll,fluxes,ivar,saltv)
	  call fluxes_accum_d(nlvdi,nsect,nlayers,dt,trs,saltt,fluxes)
	end if
	if( btemp ) then
	  ivar = 12
	  call flxscs(kfluxm,kflux,iflux,ai_ll,fluxes,ivar,tempv)
	  call fluxes_accum_d(nlvdi,nsect,nlayers,dt,trt,tempt,fluxes)
	end if
	if( bconz ) then
	  ivar = 10
	  call flxscs(kfluxm,kflux,iflux,ai_ll,fluxes,ivar,cnv)
	  call fluxes_accum_d(nlvdi,nsect,nlayers,dt,trc,conzt,fluxes)
	end if
	if( bsedi ) then
	  ivar = 800
	  call flxscs(kfluxm,kflux,iflux,ai_ll,fluxes,ivar,tcn)
	  call fluxes_accum_d(nlvdi,nsect,nlayers,dt,trsc,ssctt,fluxes)
	end if


!	-------------------------------------------------------
!	time for output?
!	-------------------------------------------------------

        if( .not. next_output_d(da_out) ) return

	call trace_point('writing flux')

        nbflx = nint(da_out(4))
        call get_absolute_act_time(atime)

!	-------------------------------------------------------
!	average and write results
!	-------------------------------------------------------

	ierr = 0
	ip = 0
	iv = 0
	nl = nl_flux
	nlg = nlv_global
	ns = ns_flux
	n = 3 * ns_flux

	ivar = 0
	iv = iv + 1
	ipstart = ip
	call fluxes_aver_d(nlvdi,nsect,nlayers,trm,masst,fluxes)
        call add_to_flx_array(ip,nl,nlg,ns,fluxes,array)
	ivinfo(:,iv) = (/ivar,ipstart,ip/)
	!call flx_collect_3d(nl,n,fluxes,flux_global(:,:,iv))

	if( bsalt ) then
	  ivar = 11
	  iv = iv + 1
	  ipstart = ip
	  call fluxes_aver_d(nlvdi,nsect,nlayers,trs,saltt,fluxes)
          call add_to_flx_array(ip,nl,nlg,ns,fluxes,array)
	  ivinfo(:,iv) = (/ivar,ipstart,ip/)
	  !call flx_collect_3d(nl,n,fluxes,flux_global(:,:,iv))
	end if
	if( btemp ) then
	  ivar = 12
	  iv = iv + 1
	  ipstart = ip
	  call fluxes_aver_d(nlvdi,nsect,nlayers,trt,tempt,fluxes)
          call add_to_flx_array(ip,nl,nlg,ns,fluxes,array)
	  ivinfo(:,iv) = (/ivar,ipstart,ip/)
	  !call flx_collect_3d(nl,n,fluxes,flux_global(:,:,iv))
	end if
	if( bconz ) then
	  ivar = 10
	  iv = iv + 1
	  ipstart = ip
	  call fluxes_aver_d(nlvdi,nsect,nlayers,trc,conzt,fluxes)
          call add_to_flx_array(ip,nl,nlg,ns,fluxes,array)
	  ivinfo(:,iv) = (/ivar,ipstart,ip/)
	  !call flx_collect_3d(nl,n,fluxes,flux_global(:,:,iv))
	end if
	if( bsedi ) then
	  ivar = 800
	  iv = iv + 1
	  ipstart = ip
	  call fluxes_aver_d(nlvdi,nsect,nlayers,trsc,ssctt,fluxes)
          call add_to_flx_array(ip,nl,nlg,ns,fluxes,array)
	  ivinfo(:,iv) = (/ivar,ipstart,ip/)
	  !call flx_collect_3d(nl,n,fluxes,flux_global(:,:,iv))
	end if

	if( iv /= nvar ) goto 91

!	-------------------------------------------------------
!	global reduce and write
!	-------------------------------------------------------

        call shympi_operate_all(.false.)
        call shympi_array_reduce('sum',array(1:ip))
        call shympi_operate_all(.true.)

	if( bmpi_master ) then
	do iv=1,nvar
	  !write(6,*) 'writing fluxes for iv,nvar: ',iv,nvar
	  ivar = ivinfo(1,iv)
	  ipstart = ivinfo(2,iv)
	  ipend = ivinfo(3,iv)
	  np = ipend-ipstart
	  n = 3 * ns_flux
          call flx_array_to_vals(np,array(ipstart+1:ipend) &
     &			,nlg,n,fluxes_global)
	  !if( any( fluxes_global /= flux_global(:,:,iv) ) ) goto 92
	  call flx_write_global(atime,ivar,fluxes_global(:,:))
	end do
	end if

!	-------------------------------------------------------
!	reset variables
!	-------------------------------------------------------

	call fluxes_init_d(nlvdi,nsect,nlayers,trm,masst)

	if( bsalt ) then
	  call fluxes_init_d(nlvdi,nsect,nlayers,trs,saltt)
	end if
	if( btemp ) then
	  call fluxes_init_d(nlvdi,nsect,nlayers,trt,tempt)
	end if

	if( bconz ) then
	  call fluxes_init_d(nlvdi,nsect,nlayers,trc,conzt)
	end if

	if( bsedi ) then
	  call fluxes_init_d(nlvdi,nsect,nlayers,trsc,ssctt)
	end if

        if( is_last_output_d(da_out) ) then
	  call flux_file_close(da_out)
	end if

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

	return
   91   continue
        write(6,*) 'iv,nvar: ',iv,nvar
        write(6,*) 'iv is different from nvar'
        stop 'error stop wrflxa: internal error (1)'
   92   continue
        write(6,*) 'arrays are different'
	write(6,*) 'nvar,iv: ',nvar,iv
	write(6,*) 'ns_flux: ',ns_flux
	write(6,*) 'ivinfo: ',ivinfo(:,iv)
	do n=1,ns_flux*3
	  do l=0,nlv_global
	    write(6,*) l,n,fluxes_global(l,n),flux_global(l,n,iv)
	  end do
	end do
        stop 'error stop wrflxa: arrays are different'
   99   continue
        write(6,*) 'Error opening FLX file :'
        stop 'error stop wrflxa: opening flx file'
   98   continue
        write(6,*) 'Error writing header of FLX file'
        write(6,*) 'unit,ierr :',nbflx,ierr
        stop 'error stop wrflxa: writing flx header'
   97   continue
        write(6,*) 'Error writing file FLX'
        write(6,*) 'unit,ierr,iv :',nbflx,ierr,iv
        stop 'error stop wrflxa: writing flx record'

	end

!**********************************************************************
!**********************************************************************
!**********************************************************************

	subroutine fluxes_generic(ext,ivbase,nscal,scal)

! administers writing of flux data
!
! serves as a template for new variables
! please adapt to your needs
!
! this routine must be called after the other fluxes have already been set up
!
! here the number of scalars and the scalar values are passed into the routine
! you can also import them through other means (modules, etc..)
!
	use mod_rungekutta, only : get_rungekutta_idiag_weights
	use levels, only : nlvdi,nlv
	use basin, only : nkn
	use mod_flux
	use simul

	implicit none

	character*(*) ext		!extension of new file (e.g., '.fxw')
	integer ivbase			!base for variable numbers
	integer nscal			!how many tracers to compute/write
	real scal(nlvdi,nkn,nscal)	!tracer

	integer i,ivar,nvers
	integer idtflx
	integer nvar,ierr
	integer kext(kfluxm)
	real ai_ll,dt
	double precision atime,atime0
	character*80 title,femver

        !double precision, save :: da_out(4)
        integer, save :: nbflx = 0

	double precision, save, allocatable :: trs(:)
	double precision, save, allocatable :: scalt(:,:,:,:)	!accumulator

	integer ifemop,ipext
	logical has_output_d,next_output_d,is_over_output_d
	double precision dgetpar

!-----------------------------------------------------------------
! start of code
!-----------------------------------------------------------------

        if( nbflx .eq. -1 ) return

!-----------------------------------------------------------------
! initialization
!-----------------------------------------------------------------

        if( nbflx .eq. 0 ) then

                call init_output_d('itmflx','idtflx',da_out)
                call increase_output_d(da_out)
                if( .not. has_output_d(da_out) ) nbflx = -1

                if( kfluxm .le. 0 ) nbflx = -1
                if( nsect .le. 0 ) nbflx = -1
                if( nbflx .eq. -1 ) return

		if( .not. bflxinit ) goto 94

        	allocate(trs(nscal))
        	allocate(scalt(0:nlvdi,3,nsect,nscal))

        	call flx_alloc_arrays(nlvdi,nsect)
		call get_nlayers(kfluxm,kflux,nsect,nlayers,nlmax)

		do i=1,nscal
		  call fluxes_init_d(nlvdi,nsect,nlayers,trs(i) &
     &				,scalt(0,1,1,i))
		end do

                nbflx = ifemop(ext,'unform','new')
                if( nbflx .le. 0 ) goto 99
		write(6,*) 'flux file opened: ',nbflx,' ',ext
		da_out(4) = nbflx

	        nvers = 0
		nvar = nscal
		idtflx = nint(da_out(1))
                call flx_write_header(nbflx,0,nsect,kfluxm,idtflx &
     &                                  ,nlmax,nvar,ierr)
                if( ierr /= 0 ) goto 98

                title = descrp
                call get_shyfem_version_and_commit(femver)
                call get_absolute_ref_time(atime0)

                do i=1,kfluxm
                  kext(i) = ipext(kflux(i))
                end do

                call flx_write_header2(nbflx,0,nsect,kfluxm &
     &                          ,kext,nlayers &
     &                          ,atime0,title,femver,chflx,ierr)
                if( ierr /= 0 ) goto 98

        end if

!-----------------------------------------------------------------
! normal call
!-----------------------------------------------------------------

        if( .not. is_over_output_d(da_out) ) return

	call get_timestep(dt)
	call get_rungekutta_idiag_weights(1,ai_ll)

!	-------------------------------------------------------
!	accumulate results
!	-------------------------------------------------------

	do i=1,nscal
	  ivar = ivbase + i
	  call flxscs(kfluxm,kflux,iflux,ai_ll,fluxes,ivar,scal(1,1,i))
	  call fluxes_accum_d(nlvdi,nsect,nlayers,dt,trs(i) &
     &			,scalt(0,1,1,i),fluxes)
	end do

!	-------------------------------------------------------
!	time for output?
!	-------------------------------------------------------

        if( .not. next_output_d(da_out) ) return

!	-------------------------------------------------------
!	average and write results
!	-------------------------------------------------------

        call get_absolute_act_time(atime)

	do i=1,nscal
	  ivar = ivbase + i
	  call fluxes_aver_d(nlvdi,nsect,nlayers,trs(i) &
     &			,scalt(0,1,1,i),fluxes)
	  fluxes_r = fluxes
          call flx_write_record(nbflx,nvers,atime,nlvdi,nsect,ivar &
     &                          ,nlayers,fluxes_r,ierr)
          if( ierr /= 0 ) goto 97
	end do

!	-------------------------------------------------------
!	reset variables
!	-------------------------------------------------------

	do i=1,nscal
	  call fluxes_init_d(nlvdi,nsect,nlayers,trs(i) &
     &			,scalt(0,1,1,i))
	end do

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

	return
   94   continue
        write(6,*) 'Flux section has not been initialized'
        stop 'error stop fluxes_template: no initialization'
   97   continue
        write(6,*) 'Error writing data record of FLX file'
        write(6,*) 'unit,ierr :',nbflx,ierr
        stop 'error stop fluxes_template: writing flx record'
   98   continue
        write(6,*) 'Error writing headers of FLX file'
        write(6,*) 'unit,ierr :',nbflx,ierr
        stop 'error stop fluxes_template: writing flx header'
   99	continue
	write(6,*) 'extension: ',ext
        stop 'error stop fluxes_template: cannot open flx file'
	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine flx_file_open(nvar)

	use mod_flux

	implicit none

	integer nvar

	call flux_file_open('.flx',da_out,nvar,nsect,kfluxm &
     &                          ,kflux_ext,nlayers,chflx)

	end

!******************************************************************

	subroutine flx_write(atime,ivar,flux_local)

	use mod_flux
	use shympi

	implicit none

	double precision atime
	integer ivar
	double precision flux_local(0:nl_flux,3,ns_flux)
	double precision flux_global(0:nlv_global,3,ns_flux)

	integer nbflx,nl,ns,n,ng

	nbflx = da_out(4)
	nl = nl_flux
	ng = nlv_global
	ns = ns_flux
	n = 3*ns

	call flx_collect_3d(nl,n,flux_local,flux_global)
        call flux_write(nbflx,atime,ivar,ng,ns &
     &                          ,nlayers,flux_global)

	end

!******************************************************************

	subroutine flx_write_global(atime,ivar,flux_global)

	use mod_flux
	use shympi

	implicit none

	double precision atime
	integer ivar
	double precision flux_global(0:nlv_global,3,ns_flux)

	integer nbflx,nl,ng,ns,n

	nbflx = da_out(4)
	nl = nl_flux
	ng = nlv_global
	ns = ns_flux
	n = 3*ns

        call flux_write(nbflx,atime,ivar,ng,ns &
     &                          ,nlayers,flux_global)

	end

!******************************************************************
!******************************************************************
!******************************************************************

        subroutine add_to_flx_array(ip,nl,nlg,ns,val,array)

        use basin

        implicit none

        integer ip
        integer k
        integer nl,nlg,ns
        double precision val(0:nl,3,ns)
        double precision array(*)

	integer ntot,nend
	integer n,m

	ntot = (nlg+1)*3*ns
	nend = ip + ntot
        array(ip+1:ip+ntot) = 0.

	do n=1,ns
	  do m=1,3
            array(ip+1:ip+nl+1) = val(0:nl,m,n)
            ip = ip + nlg + 1
	  end do
	end do

	if( nend /= ip ) stop 'error stop add_to_flx_array: internal error'

        end

!******************************************************************

        subroutine flx_array_to_vals(np,array,nlg,n,vals)

        implicit none

        integer np
        double precision array(np)
        integer nlg,n
        double precision vals(0:nlg,n)

        integer ip,j,k,i

        ip = 0
	do i=1,n
          vals(0:nlg,i) = array(ip+1:ip+nlg+1)
          ip = ip + nlg + 1
        end do

	if( ip /= np ) stop 'error stop flx_array_to_vals: internal'

        end

!******************************************************************

