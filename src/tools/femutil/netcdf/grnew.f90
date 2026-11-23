
!--------------------------------------------------------------------------
!
!    Copyright (C) 2022-2026  Georg Umgiesser
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
! 16.05.2017    ggu     started from gr.f90

!*********************************************************************

      program grnew

! this programs reads an existing nc file, 
! changes some infortmation,
! and writes a new file with the changed information

	use ncf

      implicit none

      integer nvars_in, ngatts_in, idunlim
      integer nc,ia,varid
      integer id, ncid, ncid_out
	character*80 file
	logical, save :: bvarelab = .false.
	logical, save :: battelab = .false.

	type(var_item) :: vitem
	type(att_item) :: aitem
	type(dim_item) :: ditem

	type(nc_item), pointer :: nitem_in,nitem_out

!---------------------------------------------------------------------
! open input file
!---------------------------------------------------------------------

        nc = command_argument_count()
	if( nc /= 1 ) then
	  write(6,*) 'Usage: gr file'
	  stop 'error stop: need one file'
	end if
        call get_command_argument(1,file)

!---------------------------------------------------------------------
! read input file and write general information
!---------------------------------------------------------------------

	call ncf_open_read(file,ncid)

	nitem_in => ncf_get_nitem(ncid)

	ngatts_in = nitem_in%ngatts
	idunlim = nitem_in%idunlim
	write(6,*) 'global attributes: ',ngatts_in
	write(6,*) 'unlimited dimension : ',idunlim
	write(6,*) 'dimensions: ',nitem_in%ditem%ndims
	call ncf_print_dimensions(nitem_in%ditem)

	write(6,*) 'global attributes: ',ngatts_in
	call ncf_print_attribute_header(ncid,NF_GLOBAL)
	do id=1,ngatts_in
	  aitem = nitem_in%gitems(id)
	  call ncf_print_attribute(aitem)
	end do

!---------------------------------------------------------------------
! open output file and copy dimensions
!---------------------------------------------------------------------

	call ncf_open_write('out.nc',ncid_out)

	nitem_out => ncf_get_nitem(ncid_out)

	ditem = nitem_in%ditem
	call ncf_make_dim(ncid_out,ditem%ndims,ditem%len,ditem%name)

!---------------------------------------------------------------------
! preparing attributes for variables
!---------------------------------------------------------------------

	write(6,*) 'preparing attributes for variables'
	nvars_in = nitem_in%nvars
	do varid=1,nvars_in
	  call ncf_var_inf(ncid,varid,vitem)
	  call ncf_var_make(ncid_out,vitem)
	  call ncf_print_variable(ncid,vitem)
	  do ia=1,vitem%natts
	    call ncf_att_inf(ncid,varid,ia,aitem)
	    if( battelab) call att_elab(varid,aitem)
	    call ncf_att_put(ncid_out,varid,aitem)
	    call ncf_print_attribute(aitem)
	  end do
	end do

!---------------------------------------------------------------------
! writing output file
!---------------------------------------------------------------------

	write(6,*) 'writing global attributes'
	varid = NF_GLOBAL
	do id=1,ngatts_in
	  aitem = nitem_in%gitems(id)
	  call ncf_att_put(ncid_out,varid,aitem)
	end do

	call ncf_start_data_mode(ncid_out)

	write(6,*) 'writing variables'
	do varid=1,nvars_in
	  call ncf_var_inf(ncid,varid,vitem)
	  call ncf_var_get(ncid,vitem)
	  if( bvarelab) call var_elab(varid,vitem)
	  call ncf_var_put(ncid_out,vitem)
	end do

!---------------------------------------------------------------------
! closing open files
!---------------------------------------------------------------------

	write(6,*) 'closing files'
	call ncf_close(ncid)
	call ncf_close(ncid_out)

!---------------------------------------------------------------------
! end of routine
!---------------------------------------------------------------------

	end

!*********************************************************************

	subroutine att_elab(varid,aitem)

! changes attributes for special varid

	use ncf

	implicit none

	integer varid
	type(att_item) :: aitem

	character*80 name

	if( varid /= 1 ) return

	name = aitem%name

	if( name == 'La1' .or. name == 'La2' ) then
	  write(6,*) '+++ changing attribute: ',varid,'  ',trim(name)
	  aitem%value = aitem%value + 10.
	end if

	end

!*********************************************************************

	subroutine var_elab(varid,vitem)

! changes variable for special varid

	use ncf

	implicit none

	integer varid
	type(var_item) :: vitem

	character*80 name

	if( varid /= 1 ) return

	write(6,*) '+++ changing variable: ',varid,'  ',trim(vitem%name)
	vitem%value = vitem%value + 10.

	end

!*********************************************************************

