
!--------------------------------------------------------------------------
!
!    Copyright (C) 1998,2000,2010,2014-2015,2017-2019  Georg Umgiesser
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

! extra file administration routines
!
! revision log :
!
! 08.05.1998	ggu	changed read of node numbers through nrdveci
! 20.01.2000	ggu	common block /dimdim/ eliminated
! 01.02.2000	ggu	module framework introduced
! 23.03.2010	ggu	changed v6.1.1
! 07.03.2014	ggu	changed VERS_6_1_72
! 18.06.2014	ggu	changed VERS_6_1_77
! 27.06.2014	ggu	changed VERS_6_1_78
! 26.11.2014	ggu	changed VERS_7_0_7
! 19.12.2014	ggu	changed VERS_7_0_10
! 19.01.2015	ggu	changed VERS_7_1_3
! 05.05.2015	ggu	changed VERS_7_1_10
! 20.05.2015	ggu	modules introduced
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 23.09.2015	ggu	changed VERS_7_2_4
! 18.12.2015	ggu	changed VERS_7_3_17
! 20.10.2017	ggu	new framework - read also table with strings
! 04.11.2017	ggu	changed VERS_7_5_34
! 14.11.2017	ggu	changed VERS_7_5_36
! 22.11.2017	ccf	write also waves and sediment concentration
! 05.12.2017	ggu	changed VERS_7_5_39
! 24.01.2018	ggu	changed VERS_7_5_41
! 03.04.2018	ggu	changed VERS_7_5_43
! 11.04.2018	ggu	for mpi use global vertical levels
! 19.04.2018	ggu	changed VERS_7_5_45
! 11.05.2018	ggu	changed VERS_7_5_47
! 06.07.2018	ggu	changed VERS_7_5_48
! 16.02.2019	ggu	changed VERS_7_5_60
! 02.04.2022	ggu	close ext file explicitly
! 18.09.2023	ggu	also write rho to ext file
! 04.12.2024	ggu	only do one reduce call for mpi
! 09.12.2024	ggu	extra finished and tested
! 10.11.2025	ggu	also write concv (multi concentration) to file
!
!******************************************************************
!******************************************************************
!******************************************************************

!==================================================================
        module extra
!==================================================================

        implicit none

        integer, save :: knausm = 0
        integer, save, allocatable :: knaus(:)
        integer, save, allocatable :: knext(:)
        character*80, save, allocatable :: chext(:)

!==================================================================
        contains
!==================================================================

!==================================================================
        end module extra
!==================================================================

	subroutine extra_read_section(n)

	use extra
	use nls

	integer n

	call nls_init_section

	!n = nls_read_vector()
	n = nls_read_ictable()
	knausm = n

	if( n > 0 ) then
	  allocate(knaus(n))
	  allocate(knext(n))
	  allocate(chext(n))
	  !call nls_copy_int_vect(n,knaus)
	  call nls_copy_ictable(n,knaus,chext)
	  knext = knaus
	end if

	call nls_finish_section

	end subroutine extra_read_section

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine mod_ext(mode)

	use befor_after

	implicit none

	integer mode

	if( mode .eq. M_AFTER ) then
	   call wrexta
	else if( mode .eq. M_INIT ) then
	   call inexta
	else if( mode .eq. M_READ ) then
	   call rdexta
	else if( mode .eq. M_CHECK ) then
	   call ckexta
	else if( mode .eq. M_SETUP ) then
	   call wrexta				!ggu 7/5/2001 -> write it=0
	else if( mode .eq. M_PRINT ) then
	   call prexta
	else if( mode .eq. M_TEST ) then
	   call tsexta
	else if( mode .eq. M_BEFOR ) then
!	   nothing
	else
	   write(6,*) 'unknown mode : ', mode
	   stop 'error stop mod_ext'
	end if

	end

!******************************************************************

	subroutine inexta

	implicit none

	end

!******************************************************************

	subroutine rdexta

	use extra

	implicit none

	integer n
	logical handlesec

	if( .not. handlesec('extra') ) return

	call extra_read_section(n)

	if( n .lt. 0 ) then
	  write(6,*) 'read error in section $extra'
	  stop 'error stop rdexta'
	end if

	end

!******************************************************************

	subroutine ckexta

	use extra
	use shympi

	implicit none

	integer i,ke,ki
	integer ipint
	logical bstop

	bstop = .false.

        do i=1,knausm
	  ke = knaus(i)
          ki = ipint(ke)                !$$EXTINW
          if( ki <= 0 ) then
	    if( .not. shympi_exist_node(ke) ) then
              write(6,*) 'section EXTRA : node not found ',ke
              bstop = .true.
	    end if
	  else if( .not. shympi_is_unique_node(ki) ) then
	    ki = 0
          end if
          knaus(i) = ki
        end do

	if( bstop ) stop 'error stop: ckexta'

	end

!******************************************************************

	subroutine prexta

	use extra

	implicit none

	integer i,k

        if(knausm.le.0) return

        write(6,*)
        write(6,*) 'extra section : ',knausm
	do i=1,knausm
	  k = knext(i)
          write(6,*) i,k,'  ',trim(chext(i))
	end do
        write(6,*)

	end

!******************************************************************

	subroutine tsexta

	use extra

	implicit none

	integer i

        write(6,*) '/knausc/'
        write(6,*) knausm
	do i=1,knausm
          write(6,*) i,knext(i),'  ',trim(chext(i))
	end do

	end

!******************************************************************
!******************************************************************
!******************************************************************

	subroutine wrexta

! writes and administers ext file

	use mod_hydro
	use mod_hydro_print
	use mod_ts
	use mod_conz, only : cnv,conzv
        use mod_waves, only: waveh,wavep,waved
        use mod_sediment, only : tcn
	use mod_depth
	use basin
	use levels
	use extra
	use shympi
	use simul
	use mod_trace_point

	implicit none

	integer, parameter :: mmax = 3
	logical blast
	integer nbext,ierr
	integer ivar,i,m,j,k,iv,nlv2d,nlv3d,nlv_d
	integer lmax,l,nl,nlg,ip,ipstart,ipend,iptot,ic
	integer iu
	real href,hzmin,nzadapt
	double precision atime,atime0
	character*80 femver,title
	integer kext(knausm)
	real hdep(knausm)
	real x(knausm)
	real y(knausm)
	real vals(nlv_global,knausm,mmax)
	real vals_aux(nlv_global,knausm,mmax)
	real, allocatable, save :: array(:)
	integer, save :: nvar
	integer, save :: iconz
	logical, save :: btemp,bsalt,brho,bconz,bwave,bsedi,bcmulti
	integer, save, allocatable :: il(:)
	integer, save, allocatable :: kind(:,:)
	integer, save, allocatable :: ivinfo(:,:)
	character*80 strings(knausm)

	integer ideffi,ipext
	real getpar
	logical has_output_d,next_output_d,is_last_output_d

	double precision, save :: da_out(4) = 0
	integer, save :: icall = 0

	if( icall .eq. -1 ) return

!--------------------------------------------------------------
! initialization
!--------------------------------------------------------------

	if( icall .eq. 0 ) then
	  call trace_point('initializing extra')
          call init_output_d('itmext','idtext',da_out)
	  call assure_initial_output_d(da_out)
          if( .not. has_output_d(da_out) ) icall = -1
	  if( knausm .le. 0 ) icall = -1
	  if( icall .eq. -1 ) return

          btemp = ( nint(getpar('itemp')) > 0 )
          bsalt = ( nint(getpar('isalt')) > 0 )
          brho  = ( nint(getpar('irho'))  > 0 )
          bconz = ( nint(getpar('iconz')) > 0 )
          bsedi = ( nint(getpar('isedi')) > 0 )
          bwave = ( nint(getpar('iwave')) > 0 )

          iconz = nint(getpar('iconz'))
          bcmulti = iconz > 1

          nvar = 2				!includes zeta and vel
          if( btemp ) nvar = nvar + 1
          if( bsalt ) nvar = nvar + 1
          if( brho  ) nvar = nvar + 1
          if( bconz ) nvar = nvar + iconz
          if( bsedi ) nvar = nvar + 1
          if( bwave ) nvar = nvar + 1

	  if( shympi_is_master() ) then
	    nbext=ideffi('datdir','runnam','.ext','unform','new')
            if(nbext.le.0) goto 99
	    da_out(4) = nbext
            call ext_write_header(nbext,0,knausm,nlv_global,nvar,ierr)
            if( ierr /= 0 ) goto 98
	  end if

	  allocate(il(knausm))
	  allocate(kind(2,knausm))
	  allocate(array(nlv_global*knausm*mmax*nvar))
	  allocate(ivinfo(6,nvar))
	  array = 0.
	  ivinfo = 0

	  il = 0
	  kext = 0
	  hdep = 0.
	  x = 0.
	  y = 0.
	  strings = ' '
	  title = descrp
	  href = getpar('href')
	  hzmin = getpar('hzmin')
	  nzadapt = getpar('nzadapt')
	  call collect_header(knausm,kext,hdep,x,y,il,strings,kind)
	  call get_shyfem_version_and_commit(femver)
	  call get_absolute_ref_time(atime0)
	  if( shympi_is_master() ) then
            call ext_write_header2(nbext,0,knausm,nlv_global &
     &                          ,atime0 &
     &                          ,href,hzmin,nzadapt,title,femver &
     &                          ,kext,hdep,il,x,y,strings,hlv_global &
     &                          ,ierr)
            if( ierr /= 0 ) goto 96
	  end if
        end if

	icall = icall + 1

!--------------------------------------------------------------
! write file ext
!--------------------------------------------------------------

        if( .not. next_output_d(da_out) ) return

	call trace_point('writing extra')

        nbext = nint(da_out(4))
	call get_absolute_act_time(atime)

	nlv2d = nlv_global
	nlv3d = nlv_global
	!nlv2d = 1	!to be tested
	nlg = nlv_global
	vals = 0.
	iv = 0
	ip = 0

!	-------------------------------------------------------
!	barotropic velocities and water level
!	-------------------------------------------------------

	ivar = 1
	iv = iv + 1
	m = 3
	if( iv > nvar ) goto 88
	if( m > mmax ) goto 88
	nl = 1
	ipstart = ip
	do j=1,knausm
	  k = knaus(j)
	  call add_to_array(ip,k,nl,nlg,up0v,array)
	  call add_to_array(ip,k,nl,nlg,vp0v,array)
	  call add_to_array(ip,k,nl,nlg,znv,array)
	end do
	ivinfo(:,iv) = (/ivar,nl,nlv2d,m,ipstart,ip/)

!	-------------------------------------------------------
!	velocities
!	-------------------------------------------------------

	ivar = 2
	iv = iv + 1
	m = 2
	if( iv > nvar ) goto 88
	if( m > mmax ) goto 88
	nl = nlv
	ipstart = ip
	do j=1,knausm
	  k = knaus(j)
	  call add_to_array(ip,k,nl,nlg,uprv,array)
	  call add_to_array(ip,k,nl,nlg,vprv,array)
	end do
	ivinfo(:,iv) = (/ivar,nl,nlv3d,m,ipstart,ip/)

!	-------------------------------------------------------
!	temperature
!	-------------------------------------------------------

	m = 1
	if( m > mmax ) goto 88
	nl = nlv

	ipstart = ip
	if( btemp ) then
	  ivar = 12
	  iv = iv + 1
	  if( iv > nvar ) goto 88
	  do j=1,knausm
	    k = knaus(j)
	    call add_to_array(ip,k,nl,nlg,tempv,array)
	  end do
	  ivinfo(:,iv) = (/ivar,nl,nlv3d,m,ipstart,ip/)
	end if

!	-------------------------------------------------------
!	salinity
!	-------------------------------------------------------

	ipstart = ip
	if( bsalt ) then
	  ivar = 11
	  iv = iv + 1
	  if( iv > nvar ) goto 88
	  do j=1,knausm
	    k = knaus(j)
	    call add_to_array(ip,k,nl,nlg,saltv,array)
	  end do
	  ivinfo(:,iv) = (/ivar,nl,nlv3d,m,ipstart,ip/)
	end if

!       -------------------------------------------------------
!       density
!       -------------------------------------------------------

	ipstart = ip
        if( brho ) then
          ivar = 13
          iv = iv + 1
	  if( iv > nvar ) goto 88
          do j=1,knausm
            k = knaus(j)
	    call add_to_array(ip,k,nl,nlg,rhov,array)
          end do
	  ivinfo(:,iv) = (/ivar,nl,nlv3d,m,ipstart,ip/)
        end if

!	-------------------------------------------------------
!	concentration
!	-------------------------------------------------------

	if( bconz ) then
          do ic=1,iconz
	    ipstart = ip
            ivar = 300 + ic
            if( .not. bcmulti ) ivar = 10
	    iv = iv + 1
	    if( iv > nvar ) goto 88
	    do j=1,knausm
	      k = knaus(j)
              if( bcmulti ) then
                call add_to_array(ip,k,nl,nlg,conzv(:,:,ic),array)
              else
                call add_to_array(ip,k,nl,nlg,cnv,array)
              end if
            end do
	    ivinfo(:,iv) = (/ivar,nl,nlv3d,m,ipstart,ip/)
	  end do
	end if

!	-------------------------------------------------------
!	total suspended sediment concentration
!	-------------------------------------------------------

	ipstart = ip
	if( bsedi ) then
	  ivar = 800
	  iv = iv + 1
	  if( iv > nvar ) goto 88
	  do j=1,knausm
	    k = knaus(j)
	    call add_to_array(ip,k,nl,nlg,tcn,array)
	  end do
	  ivinfo(:,iv) = (/ivar,nl,nlv3d,m,ipstart,ip/)
	end if

!       -------------------------------------------------------
!       waves
!       -------------------------------------------------------

	ipstart = ip
        if( bwave ) then
          ivar = 230
          iv = iv + 1
	  m = 3
	  if( iv > nvar ) goto 88
	  if( m > mmax ) goto 88
	  nl = 1
          do j=1,knausm
            k = knaus(j)
	    call add_to_array(ip,k,nl,nlg,waveh,array)
	    call add_to_array(ip,k,nl,nlg,wavep,array)
	    call add_to_array(ip,k,nl,nlg,waved,array)
          end do
	  ivinfo(:,iv) = (/ivar,nl,nlv2d,m,ipstart,ip/)
        end if

!--------------------------------------------------------------
! error check
!--------------------------------------------------------------

	if( iv /= nvar ) goto 91

!--------------------------------------------------------------
! reduce and write
!--------------------------------------------------------------

	call shympi_operate_all(.false.)
	call shympi_array_reduce('sum',array(1:ip))
	call shympi_operate_all(.true.)
	ierr = 0

	if( shympi_is_master() ) then
	  !write(6,*) 'nvar = ',nvar
	  do iv=1,nvar
	    ivar = ivinfo(1,iv)
	    nl = ivinfo(2,iv)
	    nlv_d = ivinfo(3,iv)
	    m = ivinfo(4,iv)
	    ipstart = ivinfo(5,iv)
	    ipend = ivinfo(6,iv)
	    vals = 0.
	    call array_to_vals(ipstart,ipend,array,nlg,knausm,m,vals)
	    !write(6,'(a,8i7)') 'ivinfo: ',ivinfo(:,iv),iv,my_id
	    !write(6,*) iv,ivar,vals(1,1,1)
            call ext_write_record(nbext,0,atime,knausm,nlv_d &
     &                                  ,ivar,m,il,vals,ierr)
            if( ierr /= 0 ) goto 97
	  end do
	end if
	call shympi_barrier

!--------------------------------------------------------------
! finalize
!--------------------------------------------------------------

        if( is_last_output_d(da_out) ) then
	  close(nbext)
	  call shympi_barrier
        end if

!--------------------------------------------------------------
! end of routine
!--------------------------------------------------------------

	return
   88   continue
	write(6,*) 'm,mmax: ',m,mmax
	write(6,*) 'iv,nvar: ',iv,nvar
	stop 'error stop wrexta: incompatible parameters'
   91   continue
	write(6,*) 'iv,nvar: ',iv,nvar
	write(6,*) 'iv is different from nvar'
	stop 'error stop wrexta: internal error (1)'
   92   continue
	iu = 750 + my_id
	write(iu,*) 'iv,ivar: ',iv,ivar
	write(iu,*) 'm,knausm,nlv_d: ',m,knausm,nlv_d
	stop 'error stop wrexta: values differ'
   99   continue
	write(6,*) 'Error opening EXT file :'
	stop 'error stop wrexta: opening ext file'
   96   continue
	write(6,*) 'Error writing second header of EXT file'
	write(6,*) 'unit,ierr :',nbext,ierr
	stop 'error stop wrexta: writing second ext header'
   98   continue
	write(6,*) 'Error writing first header of EXT file'
	write(6,*) 'unit,ierr :',nbext,ierr
	stop 'error stop wrexta: writing first ext header'
   97   continue
	write(6,*) 'Error writing record of EXT file'
	write(6,*) 'unit,iv,ierr :',nbext,iv,ierr
	write(6,*) btemp,bsalt,brho,bconz,bwave,bsedi
	stop 'error stop wrexta: writing ext record'
	end

!*********************************************************

	subroutine collect_header(n,kext,hdep,x,y,il,strings,kind)

	use basin
	use levels
	use mod_depth
	use extra
	use shympi

	implicit none

	integer n
	integer kext(n)
	real hdep(n)
	real x(n)
	real y(n)
	integer il(n)
	character*(*) strings(n)
	integer kind(2,n)

	integer j,k,ke,ki,id

	if( n /= knausm ) stop 'error stop collect_header: internal error'

	do j=1,knausm
	  ke = knext(j)
	  call shympi_find_node(ke,ki,id)
	  kind(1,j) = ki
	  kind(2,j) = id
	end do

	do j=1,knausm
	  k = knaus(j)
          kext(j) = knext(j)
	  strings(j) = chext(j)

	  call shympi_collect_node_value(k,hkv_max,hdep(j))
	  call shympi_collect_node_value(k,xgv,x(j))
	  call shympi_collect_node_value(k,ygv,y(j))
	  call shympi_collect_node_value(k,ilhkv,il(j))
	end do

	call shympi_barrier
	return

	if( shympi_is_master() ) then
	  do j=1,knausm
	    ke = knext(j)
	    write(6,*) '++ ',ke,x(j),y(j),hdep(j)
	  end do
	end if
	call shympi_finalize
	stop

	end

!*********************************************************

	subroutine add_to_array(ip,k,nl,nlg,val,array)

	use basin

	implicit none

	integer ip
	integer k
	integer nl,nlg
	real val(nl,nkn)
	real array(*)

	array(ip+1:ip+nlg) = 0.
	if( k > 0 ) then
	  array(ip+1:ip+nl) = val(1:nl,k)
	end if
	  
	ip = ip + nlg
	
	end

!*********************************************************

	subroutine array_to_vals(ipstart,np,array,nlg,knausm,m,vals)

	implicit none

	integer ipstart,np
	real array(np)
	integer nlg,knausm,m
	real vals(nlg,knausm,m)

	integer ip,j,k

	ip = ipstart
	do j=1,knausm
	  do k=1,m
	    vals(1:nlg,j,k) = array(ip+1:ip+nlg)
	    ip = ip + nlg
	  end do
	end do

	end

!*********************************************************

