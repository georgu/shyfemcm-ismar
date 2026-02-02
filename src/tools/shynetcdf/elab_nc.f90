
!--------------------------------------------------------------------------
!
!    Copyright (C) 2016-2017,2019-2020  Georg Umgiesser
!    Copyright (C) 2018  Christian Ferrarin
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

! netcdf routines for shyelab: elab_nc
!
! revision log :
!
! 07.06.2016	ggu	changed VERS_7_5_12
! 11.07.2017	ggu	changed VERS_7_5_30
! 05.12.2017	ggu	changed VERS_7_5_39
! 25.10.2018	ccf	write field already on regular grid
! 16.02.2019	ggu	changed VERS_7_5_60
! 14.05.2019	ggu	changes in nc_output_record()
! 16.05.2019	ggu	new aux variable to write nc_records (rncaux)
! 16.05.2019	ggu	also hydro output fixed
! 29.01.2020	ggu	always write variables as 3d
! 05.08.2024	ggu	new variable sncglobal
! 03.10.2024	ggu	introduced dim and dir
! 30.01.2026	ggu	change in call to get_lmax_reg() to avoid out of bounds
!
!************************************************************

!===============================================================
        module netcdf_out
!===============================================================

        implicit none

!	next vars have to be passed from outside

        logical, save :: breg = .false.
        integer, save :: nxreg = 0
        integer, save :: nyreg = 0

        real, save, allocatable :: hcoord(:)
        real, save, allocatable :: fmreg(:,:)
        real, save, allocatable :: xlon(:)
        real, save, allocatable :: ylat(:)

!	next vars are set internally

        integer, save :: iwrite_nc = 0
        integer, save :: lmaxreg = 0

        integer, save, allocatable :: var_ids(:)
        integer, save, allocatable :: var_dims(:)
        real, save, allocatable :: var3d(:)

        real, save, allocatable :: value2d(:,:)
        real, save, allocatable :: value3d(:,:,:)
        real, save, allocatable :: vnc3d(:,:,:)

!===============================================================
        end module netcdf_out
!===============================================================

	subroutine nc_output_init(ncid,title,nvar,ivars,b2d,sncglobal)

	use basin
	use levels
	use mod_depth
	use netcdf_out

	implicit none

	integer ncid			!id of file (return)
	character*(*) title		!name of simulation
	integer nvar			!total number of variables to be written
	integer ivars(nvar)		!variable id of SHYFEM
	logical b2d			!should write 2d fields
	character*(*) sncglobal		!file name of extra global values

	logical bugrid			!write unstructured mesh in ugrid format

	logical bhydro
	integer date0,time0
	integer lmax,i,iztype,idim,ivar,ncnlv,idir
	real, save :: ncflag = -999.

	logical has_direction_ivar
	logical is_2d

	call dts_get_date(date0,time0)
	call compute_iztype(iztype)

        ncnlv = nlv
	if ( b2d ) ncnlv = 1
	bugrid = .false.

	!write(6,*) 'ivars: ',ivars

	if( breg ) then
	  allocate(value2d(nxreg,nyreg))
	  allocate(value3d(ncnlv,nxreg,nyreg))
	  allocate(vnc3d(nxreg,nyreg,ncnlv))
	  call get_lmax_reg(nxreg,nyreg,fmreg,nel,ilhv,lmax)
	  lmax = max(1,lmax)	!at least one layer
	  if ( b2d ) lmax = 1
	  lmaxreg = lmax
	  call nc_open_reg(ncid,nxreg,nyreg,lmax &
     &				,ncflag,date0,time0,iztype)
	else
	  allocate(var3d(ncnlv*nkn))
	  bugrid = .true.
	  !call nc_open_fem(ncid,nkn,nel,ncnlv,date0,time0,iztype)
	  call nc_open_ugrid(ncid,nkn,nel,ncnlv,date0,time0,iztype)
	end if

	call nc_global(ncid,title,bugrid,sncglobal)

	bhydro = ivars(1) == 1		!this is a hydro file
	allocate(var_ids(nvar))
	allocate(var_dims(nvar))
	var_ids = 0

	idir = 0
        do i=1,nvar
	  ivar = ivars(i)
	  idim = 3
	  if( bhydro ) then
	    if( i == 2 ) cycle			!do not write second level
	    if( i >= 3 ) ivar = 2				!velocity
	  end if
	  if( has_direction_ivar(ivar) ) idir = idir + 1	!directional
	  if( is_2d(ivar) ) idim = 2
	  !write(6,*) 'nc_output_init: ivar,idim: ',i,ivar,idim,idir
	  var_dims(i) = idim
          call nc_init_variable(ncid,breg,idim,ivar,idir,ncflag,var_ids(i))
        end do

        call nc_end_define(ncid)
	!write(6,*) 'end of define...'

        if( breg ) then
          call nc_write_coords_reg(ncid,nxreg,nyreg,ncnlv &
     &					,xlon,ylat,hcoord,hlv)
        else
          call nc_write_coords_fem(ncid,nkn,nel,ncnlv,xgv,ygv &
     &					,hkv,nen3v,hlv)
        end if
	!write(6,*) 'end of nc_output_init...'

	end

!********************************************************************

	subroutine nc_output_record(ncid,var_id,var_dim,np,sv)

	use basin
	use levels
	use netcdf_out

	implicit none

	integer ncid
	integer var_id
	integer var_dim
	integer np
	!real cv3(nlvdi,nkn)
	real sv(nlvdi,np)

	integer lmax,nx,ny,iwrite
	real, allocatable :: rncaux(:,:,:)
	real, allocatable :: val2d(:)

	if( var_dim == 2 .and. nlvdi /= 1 ) then
	  write(6,*) 'nc_output_record: ',var_dim,nlvdi,nlv
	  stop 'error stop nc_output_record: not a 2d variable'
	end if

	iwrite = iwrite_nc

        if( breg ) then
	  lmax = lmaxreg
	  nx = nxreg
	  ny = nyreg
	  if( nx*ny /= np ) stop 'error stop nc_output_record: breg'
	  if( var_dim == 3 ) then
	    allocate(rncaux(nx,ny,lmax))
            call nc_rewrite_3d_reg(nlvdi,lmax,nx,ny,sv,rncaux)
            call nc_write_data_3d_reg(ncid,var_id,iwrite &
     &					,lmax,nx,ny,rncaux)
	  else
	    allocate(val2d(np))
	    val2d(:) = sv(1,:)
            !call fm2am2d(val2d,nx,ny,fmreg,value2d)
            call nc_write_data_2d_reg(ncid,var_id,iwrite,nx,ny,val2d)
	  end if
	else
	  if( nkn /= np ) stop 'error stop nc_output_record: nkn'
	  if( var_dim == 3 ) then
            call nc_compact_3d(nlv,nlv,nkn,sv,var3d)
            call nc_write_data_3d(ncid,var_id,iwrite,nlv,nkn,var3d)
	  else
	    allocate(val2d(np))
	    val2d(:) = sv(1,:)
            call nc_write_data_2d(ncid,var_id,iwrite,nkn,val2d)
	  end if
        end if

	end

!********************************************************************

	subroutine nc_output_hydro(ncid,znv,uprv,vprv)

	use basin
	use levels
	use netcdf_out

	implicit none

	integer ncid
	real znv(nkn)
	real uprv(nlvdi,nkn)
	real vprv(nlvdi,nkn)

	integer lmax,nx,ny,var_id,iwrite
	real, allocatable :: rncaux(:,:,:)

	iwrite = iwrite_nc

        if( breg ) then
	  lmax = lmaxreg
	  nx = nxreg
	  ny = nyreg
	  allocate(rncaux(nx,ny,lmax))
	  !write(6,*) 'nc_output_hydro: ',nx,ny,lmax,nlv,nlvdi

	  var_id = var_ids(1)
          call fm2am2d(znv,nx,ny,fmreg,value2d)
          call nc_write_data_2d_reg(ncid,var_id,iwrite,nx,ny,value2d)

	  var_id = var_ids(3)
          call fm2am3d(nlv,ilhv,uprv,nlvdi,nx,ny,fmreg,value3d)
          call nc_rewrite_3d_reg(nlvdi,lmax,nx,ny,value3d,rncaux)
          call nc_write_data_3d_reg(ncid,var_id,iwrite &
     &				,lmax,nx,ny,rncaux)

	  var_id = var_ids(4)
          call fm2am3d(nlv,ilhv,vprv,nlvdi,nx,ny,fmreg,value3d)
          call nc_rewrite_3d_reg(nlvdi,lmax,nx,ny,value3d,rncaux)
          call nc_write_data_3d_reg(ncid,var_id,iwrite &
     &				,lmax,nx,ny,rncaux)
        else
	  var_id = var_ids(1)
          call nc_write_data_2d(ncid,var_id,iwrite,nkn,znv)

	  var_id = var_ids(3)
          call nc_compact_3d(nlv,nlv,nkn,uprv,var3d)
          call nc_write_data_3d(ncid,var_id,iwrite,nlv,nkn,var3d)

	  var_id = var_ids(4)
          call nc_compact_3d(nlv,nlv,nkn,vprv,var3d)
          call nc_write_data_3d(ncid,var_id,iwrite,nlv,nkn,var3d)
        end if

	end

!********************************************************************

	subroutine nc_output_time(ncid,dtime)

	use netcdf_out

	implicit none

	integer ncid
	double precision dtime

	integer, save :: icall = 0
	double precision, save :: dtime_old = 0

	if( icall == 0 ) dtime_old = dtime - 1.
	icall = icall + 1
	if( dtime_old == dtime ) return			!already written

	iwrite_nc = iwrite_nc + 1

	call nc_write_dtime(ncid,iwrite_nc,dtime)

	end

!********************************************************************

	subroutine nc_output_final(ncid)

	implicit none

	integer ncid

	call nc_close(ncid)

	end

!********************************************************************

	subroutine nc_output_set_vars(br,nx,ny,hc,fm,x,y)

	use netcdf_out

	implicit none

        logical :: br
        integer :: nx
        integer :: ny
        real :: hc(nx*ny)
        real :: fm(4,nx*ny)
        real :: x(nx)
        real :: y(ny)

	integer nxy

	breg = br
	nxreg = nx
	nyreg = ny

	nxy = nx*ny
        allocate(hcoord(nxy))
        allocate(fmreg(4,nxy))
        allocate(xlon(nx),ylat(ny))

	hcoord = hc
	fmreg = fm
	xlon = x
	ylat = y

	end

!********************************************************************

	subroutine nc_output_get_var_id(iv,var_id)

	use netcdf_out

	implicit none

	integer :: iv
	integer :: var_id

	var_id = var_ids(iv)

	end

!********************************************************************

	subroutine nc_output_get_var_dim(iv,var_dim)

	use netcdf_out

	implicit none

	integer :: iv
	integer :: var_dim

	var_dim = var_dims(iv)

	end

!********************************************************************

