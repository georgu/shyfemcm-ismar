
!--------------------------------------------------------------------------
!
!    Copyright (C) 2026  Georg Umgiesser
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

! routine to combine various nc files into one
!
! revision log :
!
! 12.11.2026    ggu     written from scratch

!===================================================================
	module mod_ncucombine
!===================================================================

	logical, save :: binfo = .false.
	logical, save :: bdebug = .false.
	logical, save :: bverbose = .false.
	logical, save :: bquiet = .false.
	logical, save :: bsilent = .false.
	logical, save :: bwrite = .false.
	character*80, save :: variable
	character*80, save :: attribute
	integer, save :: irec = 0

!===================================================================
	end module mod_ncucombine
!===================================================================

        program ncucombine

	use ncf
	use mod_ncucombine
	use dts
	use clo

        implicit none

        integer ncid,ncid_out,ncid_in
	integer, allocatable :: ncidss(:)

	integer ierr
	integer n
        integer nvars,natts,ngatts,idunlim,nrecs,nrecsmax
        integer nc,ia,varid,attid

	real xmin,ymin,xmax,ymax,dx,dy
	real regpar(7)

	integer nt,nx,ny,nz
	integer nxx,nyy
	integer nit

        integer id
	character*80 ncfile
	character*80 file
	character*80 name
	character*80 varname
	character*20 aline

	type(var_item) :: vitem
	type(att_item) :: aitem
	type(dim_item) :: ditem

	type(var_item) :: vzitem,vtitem

	type(nc_item) :: nitem
	type(nc_item), allocatable :: nitemss(:)

	logical bexist
	integer n3
	integer ic
	integer nkn,nel
	integer varcid,varxid,varyid
	integer vartid,varzid
	integer varttid,varzzid
	integer ndims
	integer dimids(10)
	integer dimlens(10)
	real, save :: flag = -999.
	double precision, save :: dflag = -999.
	double precision w
	double precision dtime,dtime1,dtime2,atime,atimeglobal
	double precision, allocatable :: weight(:)
	double precision, allocatable :: dtimes(:),atimes(:)
	double precision, allocatable :: x(:), y(:), c(:,:), z(:)
	double precision, allocatable :: areg(:,:,:)
	double precision, allocatable :: aaux(:,:)
	logical, allocatable :: bdry(:,:)
	logical, allocatable :: bdryacum(:,:)
	real, allocatable :: raux(:,:)
	real, allocatable :: rx(:), ry(:), rz(:)
	real, allocatable :: rxx(:)
	real, allocatable :: ryy(:)
	character*256 units

!---------------------------------------------------------------------
! open first netcdf file
!---------------------------------------------------------------------

	call ncucombine_init

        nc = clo_number_of_files()

	write(6,*) 'number of files: ',nc
	if( nc <= 0 ) then
	  stop 'error stop ncucombine: internal error (1)'
	end if

	ic = 1
        call clo_get_file(ic,file)
	inquire(file=file,exist=bexist)
	if( .not. bexist ) goto 92
	call ncf_open_read(file,ncid)
	write(6,*) '====================================='
	write(6,*) 'first file: ',trim(file),ncid
	write(6,*) '====================================='
	nitem = ncf_get_nitem(ncid)
	call identify_variable(ncid,'longitude',varxid,nrecs)
	call get_xy_dim(ncid,varxid,nxx)
	call identify_variable(ncid,'latitude',varyid,nrecs)
	call get_xy_dim(ncid,varyid,nyy)
	allocate(rx(nxx),ry(nyy))
	allocate(rxx(nxx),ryy(nyy))
	allocate(aaux(nxx,nyy))
	allocate(raux(nxx,nyy))
	call get_xy(ncid,varxid,nxx,rxx)
	call get_xy(ncid,varyid,nyy,ryy)

        vartid = 0
        name = 'time'
        call identify_variable(ncid,name,vartid,nrecs)
        if( vartid == 0 ) goto 96
	call ncf_init_time(ncid,vartid,nit)
	if( nit /= nrecs ) stop 'error stop ncucombine: nrecs /= nit'
        call ncf_var_inf(ncid,vartid,vitem)
        call ncf_var_get(ncid,vitem)
	call ncf_att_name_string(ncid,vartid,'units',units)
	call ncf_get_reference_time(atime0)
	call dts_format_abs_time(atime0,aline)
	write(6,*)  'atime0 ',atime0,aline
	write(6,*)  'units ',trim(units)
	atimeglobal = atime0
	allocate(dtimes(nit),atimes(nit))
	allocate(weight(nit))
	allocate(areg(nxx,nyy,nit))
	allocate(bdry(nxx,nyy))
	allocate(bdryacum(nxx,nyy))
	bdryacum = .false.
	areg = 0.
	weight = 0.
        vtitem = vitem
	call ncf_get_data(ncid,vartid,dtimes)
	atimes(:) = dtimes(:) + atime0

        varzid = 0
        name = 'water_surface_height_above_reference_datum'
        call identify_variable(ncid,name,varzid,nrecs)
        if( varzid == 0 ) goto 96
        call ncf_var_inf(ncid,varzid,vitem)
        !call ncf_print_variable(ncid,varzid)
        vzitem = vitem

	call ncf_close(ncid)

!---------------------------------------------------------------------
! loop over all files and accumulate
!---------------------------------------------------------------------

	nrecsmax = 0

	do ic=1,nc
          call clo_get_file(ic,file)
	  inquire(file=file,exist=bexist)
	  if( .not. bexist ) goto 92
	  call ncf_open_read(file,ncid)
	  write(6,*) '====================================='
	  write(6,*) 'new file: ',trim(file),ncid
	  write(6,*) '====================================='
	  nitem = ncf_get_nitem(ncid)
	  call identify_variable(ncid,'longitude',varxid,nrecs)
	  call get_xy_dim(ncid,varxid,nx)
	  call identify_variable(ncid,'latitude',varyid,nrecs)
	  call get_xy_dim(ncid,varyid,ny)
	  if( nx /= nxx .or. ny /= nyy ) then
	    write(6,*) 'nx: ',nx,nxx
	    write(6,*) 'ny: ',ny,nyy
	    stop 'error stop ncucombine: dimensions not compatible'
	  end if
	  call get_xy(ncid,varxid,nx,rx)
	  call get_xy(ncid,varyid,ny,ry)
	  if( any(rxx/=rx) ) stop 'error stop ncucombine: rx not compatible'
	  if( any(ryy/=ry) ) stop 'error stop ncucombine: ry not compatible'

          vartid = 0
          name = 'time'
          call identify_variable(ncid,name,vartid,nrecs)
	  if( nrecs > nrecsmax ) nrecsmax = nrecs
          if( vartid == 0 ) goto 96
	  call ncf_init_time(ncid,vartid,nit)
	  if( nit /= nrecs ) stop 'error stop ncucombine: nrecs /= nit'
          call ncf_var_inf(ncid,vartid,vitem)
          call ncf_var_get(ncid,vitem)
	  call ncf_get_reference_time(atime0)
	  call dts_format_abs_time(atime0,aline)
	  write(6,*)  'atime0 ',atime0,aline
          vtitem = vitem

          varzid = 0
          name = 'water_surface_height_above_reference_datum'
          call identify_variable(ncid,name,varzid,nrecs)
          if( varzid == 0 ) goto 96
          call ncf_var_inf(ncid,varzid,vitem)
          vzitem = vitem

	  do irec=1,nrecs
	    call ncf_get_time(ncid,irec,atime)
	    if( atime /= atimes(irec) ) goto 91
	    call dts_format_abs_time(atime,aline)
	    call ncf_get_record(ncid,varzid,irec,aaux)
	    w = 1.
	    bdry = .false.
	    where( aaux == dflag ) bdry = .true.
	    bdryacum = bdryacum .or. bdry
	    areg(:,:,irec) = areg(:,:,irec) + w * aaux(:,:)
	    weight(irec) = weight(irec) + w
	    if( irec == nrecs) write(6,*) irec,atime,aline
	    !write(6,*) irec,atime,aline
	  end do

	  call ncf_close(ncid)
	end do

!---------------------------------------------------------------------
! topology is in rx,ry,nen3v - now write grid file
!---------------------------------------------------------------------

	do irec = 1,nrecsmax
	  w = weight(irec)
	  areg(:,:,irec) = areg(:,:,irec) / w
	  where( bdryacum ) areg(:,:,irec) = dflag
	  !write(6,*) irec,weight(irec)
	end do

!---------------------------------------------------------------------
! write combined output file
!---------------------------------------------------------------------

	ncid_in = ncid
	call prepare_output_file(ncid_in,nx,ny,units,nrecs,ncid_out)

	call ncf_start_data_mode(ncid_out)

	call ncf_var_id(ncid_out,'lon',varid)
	if( varid == 0 ) stop 'error stop ncucombine: cannot find lon'
	call ncf_var_put(ncid_out,varid,rxx)
	call ncf_var_id(ncid_out,'lat',varid)
	if( varid == 0 ) stop 'error stop ncucombine: cannot find lat'
	call ncf_var_put(ncid_out,varid,ryy)
	
	varname = 'time'
        call ncf_var_id(ncid_out,varname,varttid)
	varname = 'eta'
        call ncf_var_id(ncid_out,varname,varzzid)

	do irec=1,nrecsmax
	  atime = atimes(irec)
	  dtime = atime - atimeglobal
	  if( bverbose ) write(6,*) 'irec = ',irec,dtime
	  call ncf_var_put(ncid_out,varttid,irec,dtime)
	  raux(:,:) = areg(:,:,irec)
	  call ncf_var_put(ncid_out,varzzid,irec,raux)
	end do

	if( .not. bsilent ) then
	  write(6,*) nrecsmax,' records have been handled'
	end if

!---------------------------------------------------------------------
! close netcdf file
!---------------------------------------------------------------------

	call ncf_close(ncid_out)

!---------------------------------------------------------------------
! end of routine
!---------------------------------------------------------------------

	return
   91	continue
	write(6,*) 'times are different: ',atime,atimes(irec)
	stop 'error stop ncucombine: different atimes'
   92	continue
	write(6,*) 'no such file: ',trim(file)
	stop 'error stop ncucombine: no such file'
   96	continue
	write(6,*) 'cannot find variable ',trim(name)
	stop 'error stop ncucombine: no such variable'
   97	continue
	write(6,*) 'error in dimensions: ',nx,ny,n3
	stop 'error stop ncucombine: size of dimensions'
   98	continue
	write(6,*) 'wrong number of dimensions: ',vitem%ndims
	stop 'error stop ncucombine: number of dimensions'
	end program

!*********************************************************************

	subroutine prepare_output_file(ncid_in,nx,ny,units,nrecs,ncid)

	use ncf

	implicit none

	integer ncid_in
	integer nx,ny,nrecs
	integer ncid
	character*(*) units

	integer, parameter :: ndims = 3

	integer vartid,varxid,varyid,varzid
	integer len(ndims)
	character*9 name(ndims)

	call ncf_open_write('out.nc',ncid)
	len = (/nx,ny,0/)
	name = (/'longitude','latitude ','time     '/)
	call ncf_make_dim(ncid,ndims,len,name)

	call ncf_make_var(ncid,'time',nfc_double_id,1,(/3/),vartid)
	call ncf_make_var(ncid,'lon',nfc_real_id,1,(/1/),varxid)
	call ncf_make_var(ncid,'lat',nfc_real_id,1,(/2/),varyid)
	call ncf_make_var(ncid,'eta',nfc_real_id,3,(/1,2,3/),varzid)

	call make_attributes_time(ncid,vartid,units)
	call make_attributes_longitude(ncid,varxid)
	call make_attributes_latitude(ncid,varyid)
	call make_attributes_waterlevel(ncid,varzid)

	end

!*********************************************************************

	subroutine write_grd(nkn,nel,rx,ry,ic,flag,rz)

! writes grid from information contained in ugrid file

	implicit none

	integer nkn,nel
	real rx(nkn)
	real ry(nkn)
	integer ic(3,nel)
	real flag
	real rz(nkn)

	integer k,ie
	integer iu
	integer res(2)
	logical, save :: bconnect = .true.
	logical, save :: bcorrect = .false.
	real z

	iu = 20
	open(iu,file='ugrid.grd',status='unknown',form='formatted')
	write(iu,*)

	do k=1,nkn
	  z = rz(k)
	  if( z == flag ) then
	    write(iu,1000) 1,k,0,rx(k),ry(k)
	  else
	    write(iu,1000) 1,k,0,rx(k),ry(k),z
	  end if
 1000	  format(i1,2i10,3f14.6)
	end do

	write(iu,*)

	res = findloc(ic,0)
	if( any(res/=0) ) then
	  bcorrect = .true.
	  write(6,*) 'zeros found in index... must correct res = ',res
	end if

	if( bcorrect ) then
	  write(6,*) 'warning: correcting element index'
	  ic = ic + 1
	end if

	if( bconnect ) then
	do ie=1,nel
	  write(iu,2000) 2,ie,0,3
	  write(iu,*) '   ',ic(:,ie)
 2000	  format(i1,3i10)
	end do
	end if

	write(iu,*)
	close(iu)

	end subroutine

!*********************************************************************

	subroutine identify_variable(ncid,name,varid,nrecs)

! identifies variable in nc file

	use ncf
	use mod_ncucombine

	implicit none

	integer ncid
	character*(*) name
	integer varid
	integer nrecs

	logical :: bprint
	integer nvars,natts,attid,idunlim,itime,ndims
        type(nc_item) :: nitem
        type(var_item) :: vitem
        type(att_item) :: aitem
        type(dim_item) :: ditem
	character*80 vstring,astring,attname,varname

	bprint = .not. bquiet

	nitem = ncf_get_nitem(ncid)
        call ncf_set_var_name_length(ncid,nitem)
	nvars = nitem%nvars

        do varid=1,nvars
          call ncf_var_inf(ncid,varid,vitem)
          call ncf_natts(vitem,natts)
	  attname = 'standard_name'
	!write(6,*) 'eeeeeeeeeeeeeee ',varid
          call ncf_att_id(ncid,varid,attname,attid)
          if( attid > 0 ) then
            call ncf_att_inf(ncid,varid,attid,aitem)
            vstring = vitem%name
            astring = aitem%string
            if( bprint ) write(6,*) varid,attid,'  ' &
     &                          ,trim(vstring),'  ',trim(astring)
            if( astring == name ) exit
          end if
        end do

	!write(6,*) 'eeeeeeeeeeeeeee ',varid
	if( varid > nvars ) varid = 0
	if( varid == 0 ) return

	idunlim = nitem%idunlim
        !call ncf_var_inf(ncid,idunlim,vitem)
	!call ncf_print_variable(ncid,vitem)
	ditem = nitem%ditem
	nrecs = ditem%len(idunlim)

        call ncf_var_inf(ncid,varid,vitem)
	varname = vitem%name

	if( bdebug ) then
	  write(6,*) '-------------------- debug identify_variable start'
	  write(6,*) nitem%idunlim
	  write(6,*) vitem%ndims
	  write(6,*) vitem%dimids
	  write(6,*) vitem%dims
	  write(6,*) vitem%idtime
	  write(6,*) vitem%len,vitem%rlen,vitem%tlen
	  write(6,*) '-------------------- debug identify_variable end'
	end if

!---------------------------------------------------------------------
! check if variable has time dimension
! if it has time dimension it must be the last dimension
!---------------------------------------------------------------------

	ndims = vitem%ndims
	idunlim = nitem%idunlim
	itime = findloc(vitem%dimids,idunlim,1)
	if( bverbose ) write(6,*) 'itime = ',itime
	if( itime == 0 ) then
	  if( bprint ) then
	    write(6,*) 'variable ',trim(varname),' has no time dimension'
	  end if
	else if( itime == ndims ) then
	  if( bprint ) then
	    write(6,*) 'variable ',trim(varname),' has time dimension'
	  end if
	  ndims = ndims - 1	!pop time dimension
	else
	  write(6,*) 'time dimension is not last dimension... aborting'
	  stop 'error stop identify_variable: garbled dimensions'
	end if

	if( bprint ) then
	  write(6,*) 'variable ',trim(varname),' has ',ndims &
     &			,' dimensions (excluding time)'
	end if

!---------------------------------------------------------------------
! basic checks - cannot handle more than 2 dimensions
!---------------------------------------------------------------------

	if( ndims > 2 ) then
	  write(6,*) 'cannot handle more than 2 dimensions (l,k)'
	  stop 'error stop identify_variable: ndims > 2'
	end if

	end

!*********************************************************************

	subroutine get_variable(ncid,varid,varname,ir,flag,n,rz)

! extracts variable from nc file
!
! if variable is 3D, only surface layer is returned
! if variable has time dimension, record ir is returnd
! variable can be identified either by varid or varname

	use ncf
	use mod_ncucombine

	implicit none

	integer ncid
	character*(*) varname
	integer ir
	real flag
	integer n
	real rz(n)

	logical, save :: bprint = .true.
	integer itime
	integer varid
	integer idunlim,ndims
	integer nn,nv
	type(nc_item) :: nitem
	type(var_item) :: vitem
	type(dim_item) :: ditem
	double precision, allocatable :: dz1(:)
	double precision, allocatable :: dz2(:,:)

!---------------------------------------------------------------------
! if no variable desired set rz to flag and return
!---------------------------------------------------------------------

	rz = flag
	if( varid == 0 ) then
	  if( varname == ' ' ) return
	end if

!---------------------------------------------------------------------
! get info on file and variable
!---------------------------------------------------------------------

	nitem = ncf_get_nitem(ncid)
        call ncf_set_var_name_length(ncid,nitem)

	if( varid == 0 ) call ncf_var_id(ncid,varname,varid)
	call ncf_var_inf(ncid,varid,vitem)

!---------------------------------------------------------------------
! check if variable has time dimension
! if it has time dimension it must be the last dimension
!---------------------------------------------------------------------

	ndims = vitem%ndims
	idunlim = nitem%idunlim
	itime = findloc(vitem%dimids,idunlim,1)
	if( itime == 0 ) then	!no time dimentsion
	  !ok
	else if( itime == ndims ) then
	  ndims = ndims - 1	!pop time dimension
	else
	  write(6,*) 'time dimension is not last dimension... aborting'
	  stop 'error stop get_variable: garbled dimensions'
	end if

!---------------------------------------------------------------------
! basic checks - cannot handle more than 2 dimensions
!---------------------------------------------------------------------

	if( ndims > 2 ) then
	  write(6,*) 'cannot handle more than 2 dimensions (l,k)'
	  stop 'error stop get_variable: ndims > 2'
	end if

	nn = vitem%dims(ndims)
	if( nn /= n ) then
	  write(6,*) 'variable not defined on nodes: ',n,nn
	  stop 'error stop get_variable: variable not on nodes'
	end if

!---------------------------------------------------------------------
! extract variable - if 3D only return surface layer
!---------------------------------------------------------------------

	if( ndims == 1 ) then		! 2d variable
	  allocate(dz1(nn))
	  if( itime == 0 ) then
	    call ncf_get_data(ncid,varid,dz1)
	  else
	    call ncf_get_record(ncid,varid,ir,dz1)
	  end if
	  rz = real(dz1)
	else				! 3d variable
	  nv = vitem%dims(1)
	  allocate(dz2(nv,nn))
	  if( itime == 0 ) then
	    call ncf_get_data(ncid,varid,dz2)
	  else
	    call ncf_get_record(ncid,varid,ir,dz2)
	  end if
	  rz = real(dz2(1,:))		!only surface layer
	end if

!---------------------------------------------------------------------
! end of routine
!---------------------------------------------------------------------

	end subroutine

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine ncucombine_init

! initializes command line options

	use clo
	use mod_ncucombine

	implicit none

        call clo_init('ncucombine','nc-file','1.0')

        call clo_add_info('writes grid of ugrid nc-file')

        call clo_add_sep('general options')
        call clo_add_option('debug',.false.,'write also debug messages')
        call clo_add_option('verbose',.false.,'be more verbose')
        call clo_add_option('quiet',.false.,'be quiet in execution')
        call clo_add_option('silent',.false.,'do not write anything')

        call clo_add_sep('what to do')
        call clo_add_option('var var',' ','inserts value of variable var ' &
     &				// 'into grid file')
        call clo_add_option('record ir',1,'uses record ir of variable var' &
     &				// '(default is 1)')

        call clo_add_com('exit status 0 is success')

        call clo_parse_options

        call clo_get_option('debug',bdebug)
        call clo_get_option('verbose',bverbose)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('silent',bsilent)
        call clo_get_option('var',variable)
        call clo_get_option('record',irec)

        call clo_check_files(1)
        !call clo_get_file(1,ncfile)

	if( bsilent ) bquiet = .true.
	if( bdebug ) bverbose = .true.

	call set_ncu_params(bdebug,bverbose,bquiet,bsilent)

	end subroutine 

!*********************************************************************
!*********************************************************************
!*********************************************************************

