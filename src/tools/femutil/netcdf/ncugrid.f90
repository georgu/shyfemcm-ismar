
! **************************************************************************
! 
!      This is part of the netCDF package.
!      Copyright 2006 University Corporation for Atmospheric Research/Unidata.
!      See COPYRIGHT file for conditions of use.
! 
!      This is an example which reads some surface pressure and
!      temperatures. The data file read by this program is produced
!      comapnion program sfc_pres_temp_wr.f. It is intended to illustrate
!      the use of the netCDF fortran 77 API.
! 
!      This program is part of the netCDF tutorial:
!      http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial
! 
!      Full documentation of the netCDF Fortran 77 API can be found at:
!      http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f77
! 
!      $Id: sfc_pres_temp_rd.f,v 1.8 2007/01/24 19:45:09 russ Exp $
! 
! **************************************************************************

! routine to write information on nc files
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

        integer nvars, natts, ngatts, idunlim
        integer nc,ia,varid,attid

	integer nt,nx,ny,nz
	character*80 tcoord,xcoord,ycoord,zcoord

        integer id
	character*80 ncfile

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

	call find_mesh_topology(ncid,varcid,varxid,varyid)

!---------------------------------------------------------------------
! extract variables that contain mesh topology
!---------------------------------------------------------------------

	nx = 0
	ny = 0
	nz = 0
	n3 = 0

	call ncf_var_inf(ncid,varxid,vitem)
	if( bwrite ) write(6,*) trim(vitem%name),vitem%ndims,vitem%dims
	if( vitem%ndims /= 1 ) goto 98
	nx = vitem%dims(1)
	allocate(x(nx))
	allocate(rx(nx))
	call ncf_get_data(ncid,varxid,x)
	rx = x
	
	call ncf_var_inf(ncid,varyid,vitem)
	if( bwrite ) write(6,*) trim(vitem%name),vitem%ndims,vitem%dims
	if( vitem%ndims /= 1 ) goto 98
	ny = vitem%dims(1)
	allocate(y(ny))
	allocate(ry(nx))
	call ncf_get_data(ncid,varyid,y)
	ry = y
	
	if( nx /= ny ) goto 97
	nkn = nx

	call ncf_var_inf(ncid,varcid,vitem)
	if( bwrite ) write(6,*) trim(vitem%name),vitem%ndims,vitem%dims
	if( vitem%ndims /= 2 ) goto 98
	n3 = vitem%dims(1)
	nel = vitem%dims(2)
	allocate(c(n3,nel))
	allocate(ic(n3,nel))
	call ncf_get_data(ncid,varcid,c)
	ic = nint(c)

	if( n3 /= 3 ) goto 97
	
!---------------------------------------------------------------------
! topology is in rx,ry,ic - now write grid file
!---------------------------------------------------------------------

!---------------------------------------------------------------------
! if needed, extract variable for writing
!---------------------------------------------------------------------

	allocate(rz(nkn))
	call get_variable(ncid,variable,irec,flag,nkn,rz)

!---------------------------------------------------------------------
! write grid file
!---------------------------------------------------------------------

	call write_grd(nkn,nel,rx,ry,ic,flag,rz)
	if( .not. bsilent ) write(6,*) 'grid file ugrid.grd has been written'

!---------------------------------------------------------------------
! close netcdf file
!---------------------------------------------------------------------

	call ncf_close(ncid)

!---------------------------------------------------------------------
! end of routine
!---------------------------------------------------------------------

	return
   96	continue
	write(6,*) 'cannot find variable ',trim(variable)
	stop 'error stop ncugrid: no such variable'
   97	continue
	write(6,*) 'error in dimensions: ',nx,ny,n3
	stop 'error stop ncugrid: size of dimensions'
   98	continue
	write(6,*) 'wrong number of dimensions: ',vitem%ndims
	stop 'error stop ncugrid: number of dimensions'
	end program

!*********************************************************************

	subroutine find_mesh_topology(ncid,varcid,varxid,varyid)

! finds mesh topology - returns varcid,varxid,varyid

	use ncf
	use mod_ncugrid

	implicit none

	integer ncid
	integer varcid,varxid,varyid

	type(nc_item) :: nitem
	type(var_item) :: vitem
	type(att_item) :: aitem
	type(dim_item) :: ditem
	integer nvars,varid
	integer natts,attid
	integer mtvid,mtaid
	integer ianz
	character*80 vstring,astring,attname
	character*80 varname,varcname,varxname,varyname
	character*80 vars(2)

	integer iscans

!---------------------------------------------------------------------
! looks up variable containing cf_role = "mesh_topology"
!---------------------------------------------------------------------

	attname = 'cf_role'

	mtvid = 0
	mtaid = 0

	nitem = ncf_get_nitem(ncid)
        call ncf_set_var_name_length(ncid,nitem)
	nvars = nitem%nvars
	if( bwrite ) write(6,*) 'nvars = ',nvars

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

	if( mtvid > 0 .and. .not. bsilent ) then
	  write(6,*) 'mesh topology found: ', trim(vstring),'  ',trim(astring)
	else
	  write(6,*) 'cannot find mesh topology... aborting'
	  stop 'error stop find_mesh_topology: no mesh topology found'
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

	subroutine get_variable(ncid,varname,ir,flag,n,rz)

! extracts variable from nc file
!
! if variable is 3D, only surface layer is returned
! if variable has time dimension, record ir is returnd

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
	if( varname == ' ' ) return

!---------------------------------------------------------------------
! get info on file and variable
!---------------------------------------------------------------------

	nitem = ncf_get_nitem(ncid)
        call ncf_set_var_name_length(ncid,nitem)

	call ncf_var_id(ncid,varname,varid)
	call ncf_var_inf(ncid,varid,vitem)

!---------------------------------------------------------------------
! print info on file and variable
!---------------------------------------------------------------------

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
	write(6,*) 'itime = ',itime
	if( itime == 0 ) then
	  write(6,*) 'variable has no time dimension'
	else if( itime == ndims ) then
	  write(6,*) 'variable has as last dimension time dimension'
	  ndims = ndims - 1	!pop time dimension
	else
	  write(6,*) 'time dimension is not last dimension... aborting'
	  stop 'error stop get_variable: garbled dimensions'
	end if

	write(6,*) 'variable has ',ndims,' dimensions (excluding time)'

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

