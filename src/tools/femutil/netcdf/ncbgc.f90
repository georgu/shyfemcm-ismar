
!--------------------------------------------------------------------------
!
!    Copyright (C) 2017  Georg Umgiesser
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

! revision log :
!
! 16.05.2017    ggu     written from scratch

!*********************************************************************

!===================================================================
	module mod_ncbgc
!===================================================================

	logical, save :: binfo = .false.
	logical, save :: bverbose = .false.
	logical, save :: bquiet = .false.
	logical, save :: bsilent = .false.
	character*80, save :: variable
	character*80, save :: attribute

!===================================================================
	end module mod_ncbgc
!===================================================================

        program ncbgc

	use ncf
	use mod_ncbgc

        implicit none
        !include 'netcdf.inc'

        integer ncid

        integer nvars, ngatts, idunlim
        integer nc,ia,varid,attid

	integer nt,nx,ny,nz
	character*80 tcoord,xcoord,ycoord,zcoord

	logical bwrite
        integer id,nbox
	character*80 ncfile,file

	type(var_item) :: vitem
	type(att_item) :: aitem
	type(dim_item) :: ditem

	type(nc_item) :: nitem

        double precision, allocatable :: xx(:,:),yy(:,:)
        double precision, allocatable :: bathy(:,:)
        double precision, allocatable :: mask(:,:,:)
        double precision, allocatable :: smask(:,:)
        double precision, allocatable :: daver(:)
        double precision, allocatable :: dbottom(:)
        double precision, allocatable :: dlayer(:)
        integer, allocatable :: layers(:,:)
        integer, allocatable :: boxes(:,:)

!---------------------------------------------------------------------
! get dimensions
!---------------------------------------------------------------------

	file='mask.nc'
	!file='mask_tgrid.nc'
	call ncf_open_read(file,ncid)
	nitem = ncf_get_nitem(ncid)
	!call ncf_file_info(nitem,.true.)

	call get_dims(nitem,nx,ny,nz)
	write(6,*) 'dimesions found: ',nx,ny,nz

	allocate(xx(nx,ny))
	allocate(yy(nx,ny))
	allocate(bathy(nx,ny))
	allocate(mask(nx,ny,nz))
	allocate(smask(nx,ny))
	allocate(layers(nx,ny))
	allocate(boxes(nx,ny))
	allocate(daver(nz))
	allocate(dbottom(nz))
	allocate(dlayer(nz))

	call read_mask(nitem,nx,ny,nz,xx,yy,mask,smask,daver)
	call ncf_close(ncid)

	file='bathy.nc'
	call ncf_open_read(file,ncid)
	nitem = ncf_get_nitem(ncid)
	!call ncf_file_info(nitem,.true.)
	call read_bathy(nitem,nx,ny,xx,yy,smask,bathy)
	call ncf_close(ncid)

	call make_layers(nx,ny,nz,xx,yy,mask,layers)
	call make_depth(nz,daver,dbottom,dlayer)

	file='line_boxes.txt'
	call make_boxes(file,nx,ny,xx,yy,boxes)
	nbox = maxval(boxes)

	write(6,*) 'finished preparing extraction'

!---------------------------------------------------------------------
! open netcdf file
!---------------------------------------------------------------------

	call ncbgc_init(ncfile)

	bwrite = .not. bquiet

!---------------------------------------------------------------------
! print info on file, dimensions and global attributes
!---------------------------------------------------------------------

	call ncf_open_read(ncfile,ncid)

	nitem = ncf_get_nitem(ncid)

	call handle_intp(nitem,nbox,nx,ny,nz,boxes,layers,mask,dlayer)

!---------------------------------------------------------------------
! close netcdf file
!---------------------------------------------------------------------

	call ncf_close(ncid)

!---------------------------------------------------------------------
! end of routine
!---------------------------------------------------------------------

	end

!*********************************************************************

	subroutine ncbgc_init(ncfile)

! initializes command line options

	use clo
	use mod_ncbgc

	implicit none

	character*(*) ncfile

        call clo_init('ncbgc','nc-file','1.0')

        call clo_add_info('information on nc-file')

        call clo_add_sep('general options')
        call clo_add_option('verbose',.false.,'be more verbose')
        call clo_add_option('quiet',.false.,'be quiet in execution')
        call clo_add_option('silent',.false.,'do not write anything')

        call clo_add_sep('what to do')
        call clo_add_option('info',.false.,'give info on nc-file')
        call clo_add_option('var var',' ','info on variable')
        call clo_add_option('att att',' ','info on attribute')

        call clo_add_com('exit status 0 is success')

        call clo_parse_options

        call clo_get_option('verbose',bverbose)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('silent',bsilent)
        call clo_get_option('info',binfo)
        call clo_get_option('var',variable)
        call clo_get_option('att',attribute)

        call clo_check_files(1)
        call clo_get_file(1,ncfile)

	if( bsilent ) bquiet = .true.

	end

!*********************************************************************

