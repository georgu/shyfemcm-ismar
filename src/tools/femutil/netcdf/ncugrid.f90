
!--------------------------------------------------------------------------
!
!    Copyright (C) 2025  Georg Umgiesser
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

! routine to write grd file from information in ugrid
!
! revision log :
!
! 15.06.2025    ggu     written from scratch

!===================================================================
	module mod_ncugrid
!===================================================================

	logical, save :: binfo = .false.
	logical, save :: bverbose = .false.
	logical, save :: bquiet = .false.
	logical, save :: bsilent = .false.
	logical, save :: bwrite = .false.
	character*80, save :: variable
	character*80, save :: attribute
	integer, save :: irec = 0

!===================================================================
	end module mod_ncugrid
!===================================================================

        program ncugrid

	use ncf
	use mod_ncugrid

        implicit none
        !include 'netcdf.inc'

        integer ncid

	integer ierr
	integer n
        integer nvars, natts, ngatts, idunlim, nrecs
        integer nc,ia,varid,attid

	integer nt,nx,ny,nz

        integer id
	character*80 ncfile
	character*80 name

	type(var_item) :: vitem
	type(att_item) :: aitem
	type(dim_item) :: ditem

	type(nc_item) :: nitem

	integer n3
	integer nkn,nel
	integer varcid,varxid,varyid
	real, save :: flag = -999.
	double precision, allocatable :: x(:), y(:), c(:,:), z(:)
	real, allocatable :: rx(:), ry(:), rz(:)
	integer, allocatable :: ic(:,:)
	integer, allocatable :: nen3v(:,:)

!---------------------------------------------------------------------
! open netcdf file
!---------------------------------------------------------------------

	call ncugrid_init(ncfile)

	bwrite = .not. bquiet

	call ncf_open_read(ncfile,ncid)

	nitem = ncf_get_nitem(ncid)

!---------------------------------------------------------------------
! find mesh topology of ugrid file
!---------------------------------------------------------------------

	call find_mesh_topology(ncid,varcid,varxid,varyid,ierr)
	if( ierr /= 0 ) then
	  call find_custom_topology(ncid,varcid,varxid,varyid,ierr)
	  if( ierr /= 0 ) then
	    stop 'error stop ncugrid: no mesh topology found'
	  end if
	end if

!---------------------------------------------------------------------
! extract variables that contain mesh topology
!---------------------------------------------------------------------

	call get_connectivity_dim(ncid,varcid,n)
	nel = n
	allocate(nen3v(3,nel))
	call get_connectivity(ncid,varcid,nel,nen3v)

	call get_xy_dim(ncid,varxid,nx)
	call get_xy_dim(ncid,varyid,ny)
	if( nx /= ny ) goto 97
	nkn = nx
	allocate(rx(nkn),ry(nkn))
	call get_xy(ncid,varxid,nkn,rx)
	call get_xy(ncid,varyid,nkn,ry)

!---------------------------------------------------------------------
! topology is in rx,ry,ic - now write grid file
!---------------------------------------------------------------------

!---------------------------------------------------------------------
! write grid file
!---------------------------------------------------------------------

	allocate(rz(nkn))
	rz = flag
	call write_grd(nkn,nel,rx,ry,nen3v,flag,rz)
	if( .not. bsilent ) write(6,*) 'grid file ugrid.grd has been written'

!---------------------------------------------------------------------
! loop over data
!---------------------------------------------------------------------

	varid = 0
	name = 'water_surface_height_above_reference_datum'
	call identify_variable(ncid,name,varid,nrecs)
	if( varid == 0 ) goto 96
	call ncf_var_inf(ncid,varid,vitem)
	call ncf_print_variable(ncid,vitem)

	do irec=1,nrecs
	  if( bverbose ) write(6,*) 'irec = ',irec
	  variable = ' '
	  call get_variable(ncid,varid,variable,irec,flag,nkn,rz)
	end do

	write(6,*) nrecs,' records have been handled'

!---------------------------------------------------------------------
! close netcdf file
!---------------------------------------------------------------------

	call ncf_close(ncid)

!---------------------------------------------------------------------
! end of routine
!---------------------------------------------------------------------

	return
   96	continue
	write(6,*) 'cannot find variable ',trim(name)
	stop 'error stop ncugrid: no such variable'
   97	continue
	write(6,*) 'error in dimensions: ',nx,ny,n3
	stop 'error stop ncugrid: size of dimensions'
   98	continue
	write(6,*) 'wrong number of dimensions: ',vitem%ndims
	stop 'error stop ncugrid: number of dimensions'
	end program

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
	use mod_ncugrid

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
	character*80 vstring,astring,attname

	bprint = bverbose

	nitem = ncf_get_nitem(ncid)
        call ncf_set_var_name_length(ncid,nitem)
	nvars = nitem%nvars

        do varid=1,nvars
          call ncf_var_inf(ncid,varid,vitem)
          call ncf_natts(vitem,natts)
	  attname = 'standard_name'
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

	if( varid > nvars ) varid = 0
	if( varid == 0 ) return

	idunlim = nitem%idunlim
        !call ncf_var_inf(ncid,idunlim,vitem)
	!call ncf_print_variable(ncid,vitem)
	ditem = nitem%ditem
	nrecs = ditem%len(idunlim)

        call ncf_var_inf(ncid,varid,vitem)

	if( bprint ) then
	  write(6,*) nitem%idunlim
	  write(6,*) vitem%ndims
	  write(6,*) vitem%dimids
	  write(6,*) vitem%dims
	  write(6,*) vitem%idtime
	  write(6,*) vitem%len,vitem%rlen,vitem%tlen
	end if

!---------------------------------------------------------------------
! check if variable has time dimension
! if it has time dimension it must be the last dimension
!---------------------------------------------------------------------

	ndims = vitem%ndims
	idunlim = nitem%idunlim
	itime = findloc(vitem%dimids,idunlim,1)
	if( bprint ) write(6,*) 'itime = ',itime
	if( itime == 0 ) then
	  write(6,*) 'variable has no time dimension'
	else if( itime == ndims ) then
	  write(6,*) 'variable has as last dimension time dimension'
	  ndims = ndims - 1	!pop time dimension
	else
	  write(6,*) 'time dimension is not last dimension... aborting'
	  stop 'error stop identify_variable: garbled dimensions'
	end if

	write(6,*) 'variable has ',ndims,' dimensions (excluding time)'

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
	use mod_ncugrid

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

	subroutine get_ncugrid_vars(bverbose_l,bquiet_l,bsilent_l)

	use mod_ncugrid

	implicit none

	logical bverbose_l,bquiet_l,bsilent_l

	bverbose_l = bverbose
	bquiet_l = bquiet
	bsilent_l = bsilent

	end

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine ncugrid_init(ncfile)

! initializes command line options

	use clo
	use mod_ncugrid

	implicit none

	character*(*) ncfile

        call clo_init('ncugrid','nc-file','1.0')

        call clo_add_info('writes grid of ugrid nc-file')

        call clo_add_sep('general options')
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

        call clo_get_option('verbose',bverbose)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('silent',bsilent)
        call clo_get_option('var',variable)
        call clo_get_option('record',irec)

        call clo_check_files(1)
        call clo_get_file(1,ncfile)

	if( bsilent ) bquiet = .true.

	end subroutine 

!*********************************************************************
!*********************************************************************
!*********************************************************************

