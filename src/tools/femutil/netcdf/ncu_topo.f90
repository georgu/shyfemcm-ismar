
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

! routines to handle topology for ugrid files
!
! revision log :
!
! 15.11.2026    ggu     written from scratch

!*********************************************************************

	subroutine find_mesh_topology(ncid,varcid,varxid,varyid,ierr)

! finds mesh topology - returns varcid,varxid,varyid

	use ncf

	implicit none

	integer, intent(in) :: ncid
	integer, intent(out) :: varcid,varxid,varyid
	integer, intent(out) :: ierr

	type(nc_item) :: nitem
	type(var_item) :: vitem
	type(att_item) :: aitem
	type(dim_item) :: ditem
	integer nvars,varid
	integer natts,attid
	integer mtvid,mtaid
	integer ianz
	logical bdebug,bverbose,bquiet,bsilent,bwrite
	character*80 vstring,astring,attname
	character*80 varname,varcname,varxname,varyname
	character*80 vars(2)

	integer iscans

!---------------------------------------------------------------------
! looks up variable containing cf_role = "mesh_topology"
!---------------------------------------------------------------------

	call get_ncu_params(bdebug,bverbose,bquiet,bsilent)
	bwrite = .not. bquiet
	bwrite = bverbose

	attname = 'cf_role'

	mtvid = 0
	mtaid = 0

	nitem = ncf_get_nitem(ncid)
        call ncf_set_var_name_length(ncid,nitem)
	nvars = nitem%nvars
	if( bwrite ) write(6,*) 'total number of variables in file = ',nvars

	do varid=1,nvars
	  call ncf_var_inf(ncid,varid,vitem)
	  call ncf_natts(vitem,natts)
	  call ncf_att_id(ncid,varid,attname,attid)
	  if( attid > 0 ) then
	    call ncf_att_inf(ncid,varid,attid,aitem)
	    vstring = vitem%name
	    astring = aitem%string
	    if( bwrite ) write(6,*) varid,attid,'  ' &
     &				,trim(vstring),'  ',trim(astring)
	    if( astring == 'mesh_topology' ) then
	      mtvid = varid
	      mtaid = attid
	    end if
	  end if
	end do

	if( mtvid > 0 ) then
	  if( .not. bsilent ) then
	    write(6,*) 'mesh topology found: ', trim(vstring),' ',trim(astring)
	  end if
	else
	  if( .not. bsilent ) then
	    write(6,*) 'cannot find ugrid mesh topology...'
	  end if
	  ierr = 1
	  return
	  !stop 'error stop find_mesh_topology: no mesh topology found'
	end if

!---------------------------------------------------------------------
! mesh topology is in mtvid
!---------------------------------------------------------------------

	varid = mtvid
	call ncf_var_inf(ncid,varid,vitem)
	call ncf_natts(vitem,natts)

!---------------------------------------------------------------------
! look for variables containg info on x,y, and connectivity
!---------------------------------------------------------------------

	attname='face_node_connectivity'
	call ncf_att_id(ncid,varid,attname,attid)
	if( attid == 0 ) goto 99
	call ncf_att_inf(ncid,varid,attid,aitem)
	varcname = aitem%string
	if( bwrite ) write(6,*) trim(attname),'  :  ',trim(varcname)
	
	attname='node_coordinates'
	call ncf_att_id(ncid,varid,attname,attid)
	if( attid == 0 ) goto 99
	call ncf_att_inf(ncid,varid,attid,aitem)
	varname = aitem%string
	if( bwrite ) write(6,*) trim(attname),'  :  ',trim(varname)
	
	ianz = iscans(varname,vars,2)
	if( ianz /= 2 ) goto 99
	varxname = vars(1)
	varyname = vars(2)

!---------------------------------------------------------------------
! get varids for variables
!---------------------------------------------------------------------

	call ncf_var_id(ncid,varcname,varcid)
	call ncf_var_id(ncid,varxname,varxid)
	call ncf_var_id(ncid,varyname,varyid)

	if( bwrite ) then
	  write(6,*) 'c  ',varcid,'  :  ',trim(varcname)
	  write(6,*) 'x  ',varxid,'  :  ',trim(varxname)
	  write(6,*) 'y  ',varyid,'  :  ',trim(varyname)
	end if

!---------------------------------------------------------------------
! end of routine
!---------------------------------------------------------------------

	return
   99	continue
	write(6,*) 'attribute not found: ',trim(attname)
	stop 'error stop find_mesh_topology: attribute not found'
	end subroutine

!*********************************************************************
	
	subroutine find_custom_topology(ncid,varcid,varxid,varyid,ierr)

! finds custom mesh topology - returns varcid,varxid,varyid

	use ncf

	implicit none

	integer, intent(in) :: ncid
	integer, intent(out) :: varcid,varxid,varyid
	integer, intent(out) :: ierr

	type(nc_item) :: nitem
	type(var_item) :: vitem
	type(att_item) :: aitem
	type(dim_item) :: ditem
	integer nvars,varid
	integer natts,attid
	integer mtvid,mtaid
	integer ianz
	logical bfound
	logical bdebug,bverbose,bquiet,bsilent,bwrite
	character*80 vstring,astring,attname
	character*80 varname,varcname,varxname,varyname
	character*80 vars(2)

	integer iscans

!---------------------------------------------------------------------
! looks up variable containing cf_role = "mesh_topology"
!---------------------------------------------------------------------

	call get_ncu_params(bdebug,bverbose,bquiet,bsilent)
	bwrite = .not. bquiet
	bwrite = bverbose

	varcid = 0
	varxid = 0
	varyid = 0

	nitem = ncf_get_nitem(ncid)
        call ncf_set_var_name_length(ncid,nitem)
	nvars = nitem%nvars
	if( bwrite ) write(6,*) 'total number of variables in file = ',nvars

	do varid=1,nvars
	  call search_var(ncid,'cf_role','face_node_connectivity',varid,bfound)
	  if( bfound ) varcid = varid
	  call search_var(ncid,'standard_name','latitude',varid,bfound)
	  if( bfound ) varyid = varid
	  call search_var(ncid,'standard_name','longitude',varid,bfound)
	  if( bfound ) varxid = varid
	end do

	if( varcid*varyid*varxid > 0 ) then
	  if( .not. bsilent ) then
	    write(6,*) 'mesh topology found: custom'
	  end if
	  ierr = 0
	  bfound = .true.
	else
	  if( .not. bsilent ) then
	    write(6,*) 'cannot find mesh topologiy...'
	  end if
	  ierr = 3
	  bfound = .false.
	end if

	if( bwrite ) then
	  if( varcid > 0 ) then
	    varid = varcid
	    call ncf_var_inf(ncid,varid,vitem)
	    vstring = vitem%name
	    write(6,*) 'c  ',varcid,'  :  ',trim(vstring)
	  end if
	  if( varxid > 0 ) then
	    varid = varxid
	    call ncf_var_inf(ncid,varid,vitem)
	    vstring = vitem%name
	    write(6,*) 'x  ',varcid,'  :  ',trim(vstring)
	  end if
	  if( varyid > 0 ) then
	    varid = varyid
	    call ncf_var_inf(ncid,varid,vitem)
	    vstring = vitem%name
	    write(6,*) 'y  ',varcid,'  :  ',trim(vstring)
	  end if
	end if

	end

!*********************************************************************

	subroutine search_var(ncid,attname,attvalue,varid,bfound)

! looks for variable with attribute name and value

	use ncf

	implicit none

	integer ncid
	character*(*) attname,attvalue
	integer varid
	logical bfound

	integer attid
	type(att_item) :: aitem
	character*80 astring

	bfound = .false.

	!call ncf_var_inf(ncid,varid,vitem)
	!call ncf_natts(vitem,natts)
	call ncf_att_id(ncid,varid,attname,attid)
	if( attid > 0 ) then
	  call ncf_att_inf(ncid,varid,attid,aitem)
	  !vstring = vitem%name
	  astring = aitem%string
	  if( astring == attvalue ) bfound = .true.
	end if

	end

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine get_connectivity_dim(ncid,varid,nel)

	use ncf

	implicit none

	integer ncid,varid
	integer nel

        logical :: bdebug,bverbose,bquiet,bsilent,bprint
	type(nc_item) :: nitem
	type(var_item) :: vitem

	call get_ncu_params(bdebug,bverbose,bquiet,bsilent)
	bprint = bdebug

        nitem = ncf_get_nitem(ncid)
        call ncf_set_var_name_length(ncid,nitem)

	call ncf_var_inf(ncid,varid,vitem)

        if( bprint ) then
	  write(6,*) '-------------------- debug get_connectivity_dim start'
          write(6,*) nitem%idunlim
          write(6,*) vitem%ndims
          write(6,*) vitem%dimids
          write(6,*) vitem%dims
          write(6,*) vitem%idtime
          write(6,*) vitem%len,vitem%rlen,vitem%tlen
	  write(6,*) '-------------------- debug get_connectivity_dim end'
        end if

	if( vitem%ndims /= 2 ) goto 98
	if( vitem%dims(1) /= 3 ) goto 98

	nel = vitem%dims(2)

	return
   98	continue
	write(6,*) 'error in onnectivity: ',vitem%ndims,vitem%dims(1)
	stop 'error stop get_connectivity_dim: error in parameters'
	end

!*********************************************************************

	subroutine get_connectivity(ncid,varid,n,nen3v)

	use ncf

	implicit none

	integer ncid,varid
	integer n
	integer nen3v(3,n)

        logical, save :: bprint = .true.
	type(nc_item) :: nitem
	type(var_item) :: vitem
	double precision, allocatable :: daux(:,:)

        nitem = ncf_get_nitem(ncid)
        call ncf_set_var_name_length(ncid,nitem)
	call ncf_var_inf(ncid,varid,vitem)

	allocate(daux(3,n))
	call ncf_get_data(ncid,varid,daux)
	nen3v = nint(daux)

	end

!*********************************************************************

	subroutine get_xy_dim(ncid,varid,nkn)

	use ncf

	implicit none

	integer ncid,varid
	integer nkn

        logical :: bdebug,bverbose,bquiet,bsilent,bprint
	type(nc_item) :: nitem
	type(var_item) :: vitem

	call get_ncu_params(bdebug,bverbose,bquiet,bsilent)
	bprint = bdebug

        nitem = ncf_get_nitem(ncid)
        call ncf_set_var_name_length(ncid,nitem)

	call ncf_var_inf(ncid,varid,vitem)

        if( bprint ) then
	  write(6,*) '-------------------- debug get_xy_dim start'
          write(6,*) nitem%idunlim
          write(6,*) vitem%ndims
          write(6,*) vitem%dimids
          write(6,*) vitem%dims
          write(6,*) vitem%idtime
          write(6,*) vitem%len,vitem%rlen,vitem%tlen
	  write(6,*) '-------------------- debug get_xy_dim end'
        end if

	if( vitem%ndims /= 1 ) goto 98

	nkn = vitem%dims(1)

	return
   98	continue
	write(6,*) 'error in xy topology: ',vitem%ndims
	stop 'error stop get_xy_dim: error in parameters'
	end

!*********************************************************************

	subroutine get_xy(ncid,varid,n,xy)

	use ncf

	implicit none

	integer ncid,varid
	integer n
	real xy(n)

	type(nc_item) :: nitem
	type(var_item) :: vitem
	double precision, allocatable :: daux(:)

        nitem = ncf_get_nitem(ncid)
        call ncf_set_var_name_length(ncid,nitem)
	call ncf_var_inf(ncid,varid,vitem)

	allocate(daux(n))
	call ncf_get_data(ncid,varid,daux)
	xy = real(daux)

	end

!*********************************************************************

