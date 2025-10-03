
!--------------------------------------------------------------------------
!
!    Copyright (C) 2011-2013,2015-2020  Georg Umgiesser
!    Copyright (C) 2011,2013  Debora Bellafiore
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

! netcdf utility routines
!
! revision log :
!
! 05.12.2011	ggu&dbf	written from scratch
! 26.03.2012	ggu	standardized implicit routines, compiler warnings
! 20.09.2012	ggu	new routines for regular output
! 25.10.2012	ggu	changed VERS_6_1_59
! 21.01.2013	ggu	routines for handling scalar variables
! 25.01.2013	ggu	new part for nos variable initialization
! 28.01.2013	dbf	different types of vertical coordinates
! 03.05.2013	ggu	changed VERS_6_1_63
! 25.09.2013	ggu	new routines for writing time series
! 25.10.2013	ggu	changed VERS_6_1_68
! 19.01.2015	ggu	changed VERS_7_1_3
! 10.07.2015	ggu	changed VERS_7_1_50
! 17.07.2015	ggu	changed VERS_7_1_53
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 15.04.2016	ggu	changed VERS_7_5_8
! 31.05.2016	ggu	changed time variable to double precision
! 07.06.2016	ggu	changed VERS_7_5_12
! 09.05.2017	ggu	changed VERS_7_5_26
! 16.05.2017	ggu	changed VERS_7_5_27
! 13.06.2017	ggu	changed VERS_7_5_29
! 11.07.2017	ggu	changed VERS_7_5_30
! 02.09.2017	ggu	changed VERS_7_5_31
! 17.11.2017	ggu	changed VERS_7_5_37
! 05.12.2017	ggu	changed VERS_7_5_39
! 24.01.2018	ggu	changed VERS_7_5_41
! 22.02.2018	ggu	changed VERS_7_5_42
! 25.10.2018	ggu	changed VERS_7_5_51
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 14.05.2019	ggu	wrong definition of dimensions in nc_write_data_3d_reg
! 16.05.2019	ggu	new version of nc_rewrite_3d_reg(), new nc_set_quiet()
! 08.01.2020	ggu	allow for double user data
! 29.01.2020	ggu	insert extra information for error message
! 22.04.2021	ggu	accept short data type (NF_SHORT)
! 20.03.2022	ggu	accept int64 data type for time variable
! 20.06.2022	ggu	maximum dimensions updated
! 09.05.2023    lrp     write vertical velocity and density
! 02.03.2024    ccf     new entries for sediments
! 27.06.2024	ggu	better info on variables
! 31.07.2024	ggu	integrated changes for ugrid, new nc_open_ugrid()
! 02.08.2024	ggu	eliminated cf_role for CF compliance
! 05.08.2024	ggu	add error message for xtype==12 formats
! 05.08.2024	ggu	new routine nc_user_defined_global()
! 03.10.2024	ggu	introduced directional variable
! 03.10.2024	ggu	new variables added
! 01.04.2025	ggu	handle variables with no CL description
! 07.04.2025	ggu	new routine handle_time_string()
! 17.06.2025	ggu	deal with description == 'unknown'
! 30.09.2025	ggu	handle variables with no CL description
! 01.10.2025	ggu	accept different reference date
!
! notes :
!
! information on unstructured grids:
! https://publicwiki.deltares.nl/display/NETCDF/Unstructured+grids
!
! ugrid conventions:
! https://ugrid-conventions.github.io/ugrid-conventions/
!
! for non dimensional vertical coordinates (sigma etc) see:
! http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.1/cf-conventions.html
! in appendic D.6 - D.9
!
! information on f77 interface description:
! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f77/
!
! this file implements CF 1.6 compliance
!
! compliance checker:
!	CF:	https://cfchecker.ncas.ac.uk/
!	ugrid:	https://github.com/pp-mo/ugrid-checks
!
! other CF compliance checker:
! 	http://puma.nerc.ac.uk/cgi-bin/cf-checker.pl
! 	http://titania.badc.rl.ac.uk/cgi-bin/cf-checker.pl
!
! still to be implemented:
!	sigma/hybrid coordinates
!	get file name
!
!******************************************************************
!******************************************************************
! module
!******************************************************************
!******************************************************************

!==================================================================
	module netcdf_params
!==================================================================

	!use netcdf

	implicit none

        integer, save :: dimids_2d(5)		!dimensions for 2D case
        integer, save :: dimids_3d(5)		!dimensions for 3D case

        integer, save :: rec_varid		!id for time
        integer, save :: coord_varid(9)		!ids for coordinates

	character*80, save :: time_dim = ' '	!dimension name of time
	character*80, save :: time_var = ' '	!variable name of time

	double precision, save :: dtime0 = 0.	!extra for different reference
	character*80, save :: tref = ' '	!reference date
	!character*80, save :: tref = '1970-01-01::00:00:00'	!reference date

	integer, save :: ntime_recs = 0		!how many time records
	logical, save :: btime_is_char = .false. !time is a string ... special

	logical, save :: bdebug_nc = .false.
	logical, save :: bquiet_nc = .false.
	logical, save :: bverbose_nc = .false.

! if b_use_cf_role == .true. then cf_convention must be CF-1.4
! CF-1.6 does not understand cf_role (-> error)
! CF-1.11 should resolve this problem, but there is no CF checker yet

	logical, save :: b_use_cf_role = .true.
	character*80, save :: cf_convention = 'CF-1.4'
	character*80, save :: ugrid_convention = 'UGRID-1.0'

	INTERFACE
	subroutine nc_handle_err(errcode,string)
	integer errcode
	character*(*), optional :: string
	end subroutine
	END INTERFACE
	INTERFACE
	function nc_has_err(errcode)
	logical nc_has_err
	integer errcode
	end function
	END INTERFACE

	include 'netcdf.inc'

!==================================================================
	end module netcdf_params
!==================================================================

!******************************************************************
!******************************************************************
! open routines
!******************************************************************
!******************************************************************

	subroutine nc_init

	use netcdf_params

	time_dim = ' '
	time_var = ' '

	end

!******************************************************************

	subroutine nc_open_reg(ncid,nx,ny,nlv,flag,date0,time0,iztype)

! opens nc file for writing regular grid

	use netcdf_params

	implicit none

        integer ncid            !identifier (return)
	integer nx,ny,nlv	!size of arrays
	real flag		!flag for no data
        integer date0,time0     !date and time of time 0
	integer iztype          !type of vertical coordinates

	integer lmax
	integer lat_varid,lon_varid,lvl_varid,dep_varid
	integer varid
	integer lvl_dimid,nx_dimid,ny_dimid,rec_dimid
	integer ltext
	integer retval
	integer matrix_dimid(2)

	character*80 file_name
	character*80 text
	character*80 what
	character*80 date

	integer nc_ichanm

!-----------------------------------------
! initialize parameters
!-----------------------------------------

	file_name = 'netcdf_reg.nc'
	file_name = 'netcdf.nc'
	file_name = 'out.nc'

!-----------------------------------------
! Create the file.
!-----------------------------------------

	retval = nf_create(FILE_NAME, nf_clobber, ncid)
	call nc_handle_err(retval,'open_reg create')

!-----------------------------------------
! Define the dimensions. The record dimension is defined to have
! unlimited length - it can grow as needed. In this example it is
! the time dimension.
!-----------------------------------------

	lmax = max(1,nlv)	!be sure to have at least one layer

	retval = nf_def_dim(ncid, 'lon', nx, nx_dimid)
	call nc_handle_err(retval,'open_reg x')
	retval = nf_def_dim(ncid, 'lat', ny, ny_dimid)
	call nc_handle_err(retval,'open_reg y')
	retval = nf_def_dim(ncid, 'level', lmax, lvl_dimid)
	call nc_handle_err(retval,'open_reg level')
	retval = nf_def_dim(ncid, 'time', NF_UNLIMITED, rec_dimid)
	call nc_handle_err(retval,'open_reg time')

	matrix_dimid(1) = nx_dimid
	matrix_dimid(2) = ny_dimid

!-----------------------------------------
! Define the coordinate variables
!-----------------------------------------

!-----------------------------------------
! Assign units attributes to coordinate variables.
!-----------------------------------------

        retval = nf_def_var(ncid, 'lon', NF_REAL, 1, nx_dimid           &
     &                          ,lon_varid)
	call nc_handle_err(retval,'open_reg')
	varid = lon_varid

	what = 'units'
	text = 'degrees_east'
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'longitude'
	call nc_define_attr(ncid,what,text,varid)

	what = 'axis'
	text = 'X'
	call nc_define_attr(ncid,what,text,varid)

!---------------------

        retval = nf_def_var(ncid, 'lat', NF_REAL, 1, ny_dimid           &
     &                          ,lat_varid)
	call nc_handle_err(retval,'open_reg')
	varid = lat_varid

	what = 'units'
	text = 'degrees_north'
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'latitude'
	call nc_define_attr(ncid,what,text,varid)

	what = 'axis'
	text = 'Y'
	call nc_define_attr(ncid,what,text,varid)

!---------------------
! for non dimensional vertical coordinates (sigma etc) see:
! http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.1/cf-conventions.html
! in appendic D.6 - D.9
!---------------------

        retval = nf_def_var(ncid, 'level', NF_REAL, 1, lvl_dimid        &
     &                          ,lvl_varid)
	call nc_handle_err(retval,'open_reg')
	varid = lvl_varid

	what = 'units'
	text = 'm'
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'depth'
	call make_vertical_coordinate(iztype,what,text)
	call nc_define_attr(ncid,what,text,varid)

	what = 'description'
	text = 'bottom of vertical layers'
	call nc_define_attr(ncid,what,text,varid)

	what = 'axis'
	text = 'Z'
	call nc_define_attr(ncid,what,text,varid)

	what = 'positive'
	text = 'down'
	call nc_define_attr(ncid,what,text,varid)

!---------------------

        retval = nf_def_var(ncid, 'total_depth', NF_REAL                &
     &                          , 2, matrix_dimid                       &
     &                          ,dep_varid)
	call nc_handle_err(retval,'open_reg')
	varid = dep_varid

	what = 'units'
	text = 'm'
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'sea_floor_depth_below_sea_surface'
	call nc_define_attr(ncid,what,text,varid)

	what = 'description'
	text = 'total depth at data grid nodes'
	call nc_define_attr(ncid,what,text,varid)

	call nc_define_range(ncid,-100.0,+10000.0,flag,varid)

!---------------------

!	retval = nf_def_var(ncid, 'time', NF_INT, 1, rec_dimid
!     +				,rec_varid)

        retval = nf_def_var(ncid, 'time', NF_DOUBLE, 1, rec_dimid       &
     &                          ,rec_varid)
	call nc_handle_err(retval,'open_reg')
	varid = rec_varid

	what = 'units'
	call nc_convert_date(date0,time0,date)
	text = 'seconds since '//trim(date)
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'time'
	call nc_define_attr(ncid,what,text,varid)

	what = 'calendar'
	text = 'standard'
	call nc_define_attr(ncid,what,text,varid)

	what = 'axis'
	text = 'T'
	call nc_define_attr(ncid,what,text,varid)

!-----------------------------------------
! define dimensions to remember
!-----------------------------------------

	dimids_2d(1) = nx_dimid
	dimids_2d(2) = ny_dimid
	dimids_2d(3) = rec_dimid

	dimids_3d(1) = nx_dimid
	dimids_3d(2) = ny_dimid
	dimids_3d(3) = lvl_dimid
	dimids_3d(4) = rec_dimid

	coord_varid(1) = lon_varid
	coord_varid(2) = lat_varid
	coord_varid(3) = lvl_varid
	coord_varid(4) = dep_varid

!-----------------------------------------
! end of routine
!-----------------------------------------

	end

!******************************************************************

	subroutine nc_open_fem(ncid,nkn,nel,nlv,date0,time0,iztype)

! opens nc file for writing fem grid

	use netcdf_params

	implicit none

	integer ncid		!identifier (return)
	integer nkn,nel,nlv	!size of arrays
	integer date0,time0	!date and time of time 0
	integer iztype		!type of vertical coordinates

	integer lmax
	integer lat_varid,lon_varid,lvl_varid,dep_varid
	integer eix_varid,top_varid
	integer varid
	integer lvl_dimid,node_dimid,elem_dimid,vertex_dimid,rec_dimid
	integer ltext
	integer retval
	integer eix_dimid(2)

	character*80 file_name
	character*80 text
	character*80 what
	character*80 date

	integer nc_ichanm

!-----------------------------------------
! initialize parameters
!-----------------------------------------

	file_name = 'netcdf.nc'
	file_name = 'out.nc'

!-----------------------------------------
! Create the file.
!-----------------------------------------

	retval = nf_create(FILE_NAME, nf_clobber, ncid)
	call nc_handle_err(retval,'open_fem')

!-----------------------------------------
! Define the dimensions. The record dimension is defined to have
! unlimited length - it can grow as needed. In this example it is
! the time dimension.
!-----------------------------------------

	lmax = max(1,nlv)	!be sure to have at least one layer

	retval = nf_def_dim(ncid, 'level', lmax, lvl_dimid)
	call nc_handle_err(retval,'open_fem')
	retval = nf_def_dim(ncid, 'node', nkn, node_dimid)
	call nc_handle_err(retval,'open_fem')
	retval = nf_def_dim(ncid, 'element', nel, elem_dimid)
	call nc_handle_err(retval,'open_fem')
	retval = nf_def_dim(ncid, 'vertex', 3, vertex_dimid)
	call nc_handle_err(retval,'open_fem')
	retval = nf_def_dim(ncid, 'time', NF_UNLIMITED, rec_dimid)
	call nc_handle_err(retval,'open_fem')

	eix_dimid(1) = vertex_dimid
	eix_dimid(2) = elem_dimid

!-----------------------------------------
! Define the coordinate variables
!-----------------------------------------

!-----------------------------------------
! Assign units attributes to coordinate variables.
!-----------------------------------------

        retval = nf_def_var(ncid, 'longitude', NF_REAL, 1, node_dimid   &
     &                          ,lon_varid)
	call nc_handle_err(retval,'open_fem')
	varid = lon_varid

	what = 'units'
	text = 'degrees_east'
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'longitude'
	call nc_define_attr(ncid,what,text,varid)

!---------------------

        retval = nf_def_var(ncid, 'latitude', NF_REAL, 1, node_dimid    &
     &                          ,lat_varid)
	call nc_handle_err(retval,'open_fem')
	varid = lat_varid

	what = 'units'
	text = 'degrees_north'
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'latitude'
	call nc_define_attr(ncid,what,text,varid)

!---------------------
! for non dimensional vertical coordinates (sigma etc) see:
! http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.1/cf-conventions.html
! in appendic D.6 - D.9
!---------------------

        retval = nf_def_var(ncid, 'level', NF_REAL, 1, lvl_dimid        &
     &                          ,lvl_varid)
	call nc_handle_err(retval,'open_fem')
	varid = lvl_varid

	what = 'units'
	text = 'm'
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'depth'
	call make_vertical_coordinate(iztype,what,text)
	call nc_define_attr(ncid,what,text,varid)

	what = 'description'
	text = 'bottom of vertical layers'
	call nc_define_attr(ncid,what,text,varid)

	what = 'axis'
	text = 'Z'
	call nc_define_attr(ncid,what,text,varid)

	what = 'positive'
	text = 'down'
	call nc_define_attr(ncid,what,text,varid)

!---------------------

        retval = nf_def_var(ncid, 'total_depth', NF_REAL, 1, node_dimid &
     &                          ,dep_varid)
	call nc_handle_err(retval,'open_fem')
	varid = dep_varid

	what = 'units'
	text = 'm'
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'sea_floor_depth_below_sea_surface'
	call nc_define_attr(ncid,what,text,varid)

	what = 'description'
	text = 'total depth at nodes'
	call nc_define_attr(ncid,what,text,varid)

!---------------------

        retval = nf_def_var(ncid, 'element_index', NF_INT, 2, eix_dimid &
     &                          ,eix_varid)
	call nc_handle_err(retval,'open_fem')
	varid = eix_varid

	what = 'units'
	text = '1'
	call nc_define_attr(ncid,what,text,varid)

	what = 'long_name'
	text = 'element index of nodes'
	call nc_define_attr(ncid,what,text,varid)

	what = 'description'
	text = 'maps every element to its three vertices'
	call nc_define_attr(ncid,what,text,varid)

        retval = nf_def_var(ncid, 'topology', NF_INT, 0, 0              &
     &                          ,top_varid)
	call nc_handle_err(retval,'open_fem')
	varid = top_varid

	what = 'units'
	text = '1'
	call nc_define_attr(ncid,what,text,varid)

	what = 'long_name'
	text = 'topology data of 2D unstructured mesh'
	call nc_define_attr(ncid,what,text,varid)

	what = 'dimensionality'
	text = '2'
	call nc_define_attr(ncid,what,text,varid)

!---------------------

!	retval = nf_def_var(ncid, 'time', NF_INT, 1, rec_dimid
!     +				,rec_varid)

        retval = nf_def_var(ncid, 'time', NF_DOUBLE, 1, rec_dimid       &
     &                          ,rec_varid)
	call nc_handle_err(retval,'open_fem')
	varid = rec_varid

	what = 'units'
	call nc_convert_date(date0,time0,date)
	text = 'seconds since '//trim(date)
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'time'
	call nc_define_attr(ncid,what,text,varid)

	what = 'calendar'
	text = 'standard'
	call nc_define_attr(ncid,what,text,varid)

	what = 'axis'
	text = 'T'
	call nc_define_attr(ncid,what,text,varid)

!-----------------------------------------
! define dimensions to remember
!-----------------------------------------

	dimids_2d(1) = node_dimid
	dimids_2d(2) = rec_dimid

	dimids_3d(1) = lvl_dimid
	dimids_3d(2) = node_dimid
	dimids_3d(3) = rec_dimid

	coord_varid(1) = lon_varid
	coord_varid(2) = lat_varid
	coord_varid(3) = lvl_varid
	coord_varid(4) = dep_varid
	coord_varid(5) = eix_varid
	coord_varid(6) = top_varid

!-----------------------------------------
! end of routine
!-----------------------------------------

	end

!******************************************************************

	subroutine nc_open_ugrid(ncid,nkn,nel,nlv,date0,time0,iztype)

! opens nc file for writing fem grid with ugrid standard

	use netcdf_params

	implicit none

	integer ncid		!identifier (return)
	integer nkn,nel,nlv	!size of arrays
	integer date0,time0	!date and time of time 0
	integer iztype		!type of vertical coordinates

	integer lmax
	integer lat_varid,lon_varid,lvl_varid,dep_varid
	integer eix_varid,top_varid
	integer varid
	integer lvl_dimid,node_dimid,elem_dimid,vertex_dimid,rec_dimid
	integer ltext
	integer ivalue
	integer retval
	integer eix_dimid(2)
	integer top_dimid(0)

	character*80 file_name
	character*80 text
	character*80 what
	character*80 date

	integer nc_ichanm

!-----------------------------------------
! initialize parameters
!-----------------------------------------

	file_name = 'netcdf.nc'
	file_name = 'out.nc'

!-----------------------------------------
! Create the file.
!-----------------------------------------

	retval = nf_create(FILE_NAME, nf_clobber, ncid)
	call nc_handle_err(retval,'open_fem')

!-----------------------------------------
! Define the dimensions. The record dimension is defined to have
! unlimited length - it can grow as needed. In this example it is
! the time dimension.
!-----------------------------------------

	lmax = max(1,nlv)	!be sure to have at least one layer

	retval = nf_def_dim(ncid, 'nMesh2_node', nkn, node_dimid)
	call nc_handle_err(retval,'open_fem')
	retval = nf_def_dim(ncid, 'nMesh2_face', nel, elem_dimid)
	call nc_handle_err(retval,'open_fem')
	retval = nf_def_dim(ncid, 'Three', 3, vertex_dimid)
	call nc_handle_err(retval,'open_fem')
	retval = nf_def_dim(ncid, 'level', lmax, lvl_dimid)
	call nc_handle_err(retval,'open_fem')
	retval = nf_def_dim(ncid, 'time', NF_UNLIMITED, rec_dimid)
	call nc_handle_err(retval,'open_fem')

	eix_dimid(1) = vertex_dimid
	eix_dimid(2) = elem_dimid

!-----------------------------------------
! Define the Topology
!-----------------------------------------

        retval = nf_def_var(ncid, 'Mesh2', NF_INT, 0, top_dimid &
     &                          ,top_varid)
	call nc_handle_err(retval,'open_fem')
	varid = top_varid

	what = 'cf_role'
	text = 'mesh_topology'
	if( b_use_cf_role ) call nc_define_attr(ncid,what,text,varid)

	what = 'long_name'
	text = 'topology data of 2D unstructured mesh'
	call nc_define_attr(ncid,what,text,varid)

	what = 'topology_dimension'
	ivalue = 2
	call nc_define_attr_int(ncid,what,ivalue,varid)

	what  = 'node_coordinates'
	text = 'Mesh2_node_x Mesh2_node_y'
	call nc_define_attr(ncid,what,text,varid)

	what  = 'face_node_connectivity'
	text = 'Mesh2_face_nodes'
	call nc_define_attr(ncid,what,text,varid)

!---------------------

        retval = nf_def_var(ncid, 'Mesh2_face_nodes', NF_INT, 2, eix_dimid &
     &                          ,eix_varid)
	call nc_handle_err(retval,'open_fem')
	varid = eix_varid

	what = 'cf_role'
	text = 'face_node_connectivity'
	if( b_use_cf_role ) call nc_define_attr(ncid,what,text,varid)

	what = 'long_name'
	text = 'Maps every triangular face to its three corner nodes'
	call nc_define_attr(ncid,what,text,varid)

	what = 'start_index'
	ivalue = 1
	call nc_define_attr_int(ncid,what,ivalue,varid)

	what = 'units'
	text = '1'
	call nc_define_attr(ncid,what,text,varid)

!-----------------------------------------
! Define the coordinate variables
!-----------------------------------------

        retval = nf_def_var(ncid, 'Mesh2_node_x', NF_REAL, 1, node_dimid   &
     &                          ,lon_varid)
	call nc_handle_err(retval,'open_fem')
	varid = lon_varid

	what = 'standard_name'
	text = 'longitude'
	call nc_define_attr(ncid,what,text,varid)

	what = 'long_name'
	text = 'Longitude of 2D mesh nodes'
	call nc_define_attr(ncid,what,text,varid)

	what = 'units'
	text = 'degrees_east'
	call nc_define_attr(ncid,what,text,varid)

!---------------------

        retval = nf_def_var(ncid, 'Mesh2_node_y', NF_REAL, 1, node_dimid   &
     &                          ,lat_varid)
	call nc_handle_err(retval,'open_fem')
	varid = lat_varid

	what = 'standard_name'
	text = 'latitude'
	call nc_define_attr(ncid,what,text,varid)

	what = 'long_name'
	text = 'Latitude of 2D mesh nodes'
	call nc_define_attr(ncid,what,text,varid)

	what = 'units'
	text = 'degrees_north'
	call nc_define_attr(ncid,what,text,varid)

!-----------------------------------------
! Define vertical structure (layer thickness and node total depth)
!-----------------------------------------

!---------------------
! for non dimensional vertical coordinates (sigma etc) see:
! http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.1/cf-conventions.html
! in appendic D.6 - D.9
!---------------------

        retval = nf_def_var(ncid, 'level', NF_REAL, 1, lvl_dimid        &
     &                          ,lvl_varid)
	call nc_handle_err(retval,'open_fem')
	varid = lvl_varid

	what = 'units'
	text = 'm'
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'depth'
	call make_vertical_coordinate(iztype,what,text)
	call nc_define_attr(ncid,what,text,varid)

	what = 'description'
	text = 'bottom of vertical layers'
	call nc_define_attr(ncid,what,text,varid)

	what = 'axis'
	text = 'Z'
	call nc_define_attr(ncid,what,text,varid)

	what = 'positive'
	text = 'down'
	call nc_define_attr(ncid,what,text,varid)

!---------------------

        retval = nf_def_var(ncid, 'total_depth', NF_REAL, 1, node_dimid &
     &                          ,dep_varid)
	call nc_handle_err(retval,'open_fem')
	varid = dep_varid

	what = 'units'
	text = 'm'
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'sea_floor_depth_below_sea_surface'
	call nc_define_attr(ncid,what,text,varid)

	what = 'description'
	text = 'total depth at nodes'
	call nc_define_attr(ncid,what,text,varid)

!-----------------------------------------
! Define time
!-----------------------------------------

        retval = nf_def_var(ncid, 'time', NF_DOUBLE, 1, rec_dimid       &
     &                          ,rec_varid)
	call nc_handle_err(retval,'open_fem')
	varid = rec_varid

	what = 'units'
	call nc_convert_date(date0,time0,date)
	text = 'seconds since '//trim(date)
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'time'
	call nc_define_attr(ncid,what,text,varid)

	what = 'calendar'
	text = 'standard'
	call nc_define_attr(ncid,what,text,varid)

	what = 'axis'
	text = 'T'
	call nc_define_attr(ncid,what,text,varid)

!-----------------------------------------
! define dimensions to remember
!-----------------------------------------

	dimids_2d(1) = node_dimid
	dimids_2d(2) = rec_dimid

	dimids_3d(1) = lvl_dimid
	dimids_3d(2) = node_dimid
	dimids_3d(3) = rec_dimid

	coord_varid(1) = lon_varid
	coord_varid(2) = lat_varid
	coord_varid(3) = lvl_varid
	coord_varid(4) = dep_varid
	coord_varid(5) = eix_varid
	coord_varid(6) = top_varid

!-----------------------------------------
! end of routine
!-----------------------------------------

	end

!******************************************************************

	subroutine nc_open_ts(ncid,node,date0,time0)

! opens nc file for time series write

	use netcdf_params

	implicit none

	integer ncid
	integer node		!total number of nodes for ts
	integer date0,time0

	integer lat_varid,lon_varid
	integer varid
	integer node_dimid,rec_dimid
	integer ltext
	integer retval

	character*80 file_name
	character*80 text
	character*80 what
	character*80 date

	integer nc_ichanm

!-----------------------------------------
! initialize parameters
!-----------------------------------------

	file_name = 'netcdfts.nc'
	file_name = 'netcdf.nc'

!-----------------------------------------
! Create the file.
!-----------------------------------------

	retval = nf_create(FILE_NAME, nf_clobber, ncid)
	call nc_handle_err(retval,'open_ts')

!-----------------------------------------
! Define the dimensions. The record dimension is defined to have
! unlimited length - it can grow as needed. In this example it is
! the time dimension.
!-----------------------------------------

	retval = nf_def_dim(ncid, 'node', 1, node_dimid)
	call nc_handle_err(retval,'open_ts')
	retval = nf_def_dim(ncid, 'time', NF_UNLIMITED, rec_dimid)
	call nc_handle_err(retval,'open_ts')

!-----------------------------------------
! Define the coordinate variables
!-----------------------------------------

!-----------------------------------------
! Assign units attributes to coordinate variables.
!-----------------------------------------

        retval = nf_def_var(ncid, 'longitude', NF_REAL, 1,node_dimid    &
     &                          ,lon_varid)
	write(6,*)'lon_varid',lon_varid
	call nc_handle_err(retval,'open_ts')
	varid = lon_varid

	what = 'units'
	text = 'degrees_east'
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'longitude'
	call nc_define_attr(ncid,what,text,varid)

!---------------------

        retval = nf_def_var(ncid, 'latitude', NF_REAL, 1, node_dimid    &
     &                          ,lat_varid)
	write(6,*)'lat_varid',lat_varid
	call nc_handle_err(retval,'open_ts')
	varid = lat_varid

	what = 'units'
	text = 'degrees_north'
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'latitude'
	call nc_define_attr(ncid,what,text,varid)

!---------------------

!	retval = nf_def_var(ncid, 'time', NF_INT, 1, rec_dimid
!     +				,rec_varid)

        retval = nf_def_var(ncid, 'time', NF_DOUBLE, 1, rec_dimid       &
     &                          ,rec_varid)
	call nc_handle_err(retval,'open_ts')
	varid = rec_varid

	what = 'units'
	call nc_convert_date(date0,time0,date)
	text = 'seconds since '//trim(date)
	call nc_define_attr(ncid,what,text,varid)

	what = 'standard_name'
	text = 'time'
	call nc_define_attr(ncid,what,text,varid)

	what = 'calendar'
	text = 'standard'
	call nc_define_attr(ncid,what,text,varid)

	what = 'axis'
	text = 'T'
	call nc_define_attr(ncid,what,text,varid)

!-----------------------------------------
! define dimensions to remember
!-----------------------------------------

	dimids_2d(1) = node_dimid
	dimids_2d(2) = rec_dimid

	coord_varid(1) = lon_varid
	coord_varid(2) = lat_varid

!-----------------------------------------
! end of routine
!-----------------------------------------

	end

!*****************************************************************


!*****************************************************************

        subroutine make_vertical_coordinate(iztype,what,text)

! defines definition for vertical coordinate

        implicit none

        integer iztype
        character*(*) what
        character*(*) text

        what = 'standard_name'

        if( iztype .eq. 1 ) then        ! z-coordinates
          text = 'depth'
        else if( iztype .eq. 2 ) then   ! sigma-coordinates
          text = 'ocean_sigma_coordinate'
        else if( iztype .eq. 3 ) then   ! hybrid-coordinates
          text = 'ocean_sigma_z_coordinate'
        else
          write(6,*) 'iztype = ',iztype
          stop 'error stop make_vertical_coordinate: unknown iztype'
        end if

        end

!*****************************************************************
!*****************************************************************
! read dimensions
!*****************************************************************
!*****************************************************************

	subroutine nc_open_read(ncid,file)

! opens nc file for reading fem file

	use netcdf_params

	implicit none

	integer ncid		!identifier (return)
	character*(*) file
	integer retval

	!write(6,*) 'opening nc file for read: ',trim(file)

        retval = nf_open(file, NF_NOWRITE, ncid)
	call nc_handle_err(retval,'open_read')

	!retval = nf_inq(ncid, ndims, nvars, ngatts, unlim)
	!call nc_handle_err(retval,'open_read')

	end

!*****************************************************************

	subroutine nc_dims_info(ncid)

	use netcdf_params

	implicit none

	integer ncid

	integer ndims,dim_id,length,i
	character*30 name
	integer retval

	retval = nf_inq_ndims(ncid,ndims)
	call nc_handle_err(retval,'dims_info')

	write(6,*) 'dimensions: '
	do i=1,ndims
	  dim_id = i
	  retval = nf_inq_dim(ncid,dim_id,name,length)
	  call nc_handle_err(retval,'dims_info')
	  write(6,*) dim_id,length,name
	end do

	end

!*****************************************************************

        subroutine nc_get_dim_totnum(ncid,ndims)

	use netcdf_params

        implicit none

        integer ncid
	integer ndims

        integer retval

	retval = nf_inq_ndims(ncid,ndims)
	call nc_handle_err(retval,'get_dim_totnum')

	end

!*****************************************************************

        subroutine nc_get_dim_name(ncid,dim_id,name)

	use netcdf_params

        implicit none

        integer ncid
        integer dim_id
        character*(*) name

        integer retval

        retval = nf_inq_dimname(ncid,dim_id,name)
        call nc_handle_err(retval,'get_dim_name')

        end

!*****************************************************************

	subroutine nc_get_dim_id(ncid,name,dim_id)

	use netcdf_params

	implicit none

	integer ncid
	character*(*) name
	integer dim_id

	integer retval

	retval = nf_inq_dimid(ncid,name,dim_id)
	call nc_handle_err(retval,'get_dim_id')

	end

!*****************************************************************

        subroutine nc_has_dim_name(ncid,name,dim_id)

	use netcdf_params

        implicit none

        integer ncid
        character*(*) name
	integer dim_id

        integer retval

	retval = nf_inq_dimid(ncid,name,dim_id)
	if( nc_has_err(retval) ) dim_id = 0

        end

!*****************************************************************

	subroutine nc_get_dim_len(ncid,dim_id,dim_len)

	use netcdf_params

	implicit none

	integer ncid
	integer dim_id
	integer dim_len

	integer retval

	retval = nf_inq_dimlen(ncid,dim_id,dim_len)
	call nc_handle_err(retval,'has_dim_len')

	end

!*****************************************************************
!*****************************************************************
! read time records
!*****************************************************************
!*****************************************************************

	subroutine nc_get_time_rec(ncid,irec,t)

	use netcdf_params

	implicit none

	integer ncid
	integer irec			!number of record
	double precision t		!time found (return)

	integer nvars,var_id,time_id,i
	character*80 name
	character*80 time,time_d,time_v
	integer istart,icount
	integer xtype
	integer retval
	integer itime
	real rtime
	double precision dtime
	logical bdebug

	bdebug = .true.
	bdebug = .false.

	retval = nf_inq_nvars(ncid,nvars)
	call nc_handle_err(retval,'get_time_rec')

	call nc_get_time_name(time_d,time_v)
	time = time_v

	call nc_get_var_id(ncid,time,time_id)
	if( time_id .eq. 0 ) then
	  stop 'error stop nc_get_time_rec: cannot find time variable'
	end if

	retval = nf_inq_vartype(ncid,time_id,xtype)
	call nc_handle_err(retval,'get_time_rec')
	!write(6,*) 'time_id: ',time_id,xtype

	!write(6,*) 'types: ',NF_INT,NF_FLOAT,NF_DOUBLE,NF_CHAR

	istart = irec
	icount = 1
	if( xtype .eq. NF_INT ) then
	  !write(6,*) 'time is int.........'
	  retval = nf_get_vara_int(ncid,time_id,istart,icount,itime)
	  call nc_handle_err(retval,'get_time_rec')
	  t = itime
	else if( xtype .eq. NF_FLOAT ) then
	  !write(6,*) 'time is real.........'
	  retval = nf_get_vara_real(ncid,time_id,istart,icount,rtime)
	  call nc_handle_err(retval,'get_time_rec')
	  t = rtime
	else if( xtype .eq. NF_DOUBLE ) then
	  !write(6,*) 'time is double.........'
	  retval = nf_get_vara_double(ncid,time_id,istart,icount,dtime)
	  call nc_handle_err(retval,'get_time_rec')
	  t = dtime
	else if( xtype .eq. NF_CHAR ) then
	  !write(6,*) 'time is char.........'
	  call handle_time_string(ncid,irec,t)
	else if( xtype .eq. 10 ) then	!has no name
	  !write(6,*) 'time is int64.........'
	  retval = nf_get_vara_double(ncid,time_id,istart,icount,dtime)
	  call nc_handle_err(retval,'get_time_rec')
	  t = dtime
	else
	  write(6,*) 'unknown type for variable time'
	  write(6,*) 'xtype = ',xtype
	  stop 'error stop nc_get_time_rec: cannot read time'
	end if

	end

!*****************************************************************

	subroutine handle_time_string(ncid,irec,t)

	use netcdf_params

	implicit none

	integer ncid,irec
	double precision t

	logical bdebug
	integer retval
	character*80 time,time_d,time_v
	character*80 name
	integer i,time_id
	integer xtype
	integer itime
	integer ndims,natts
	integer dimid,len,dimlen
	integer iu,it,ierr
	integer, save :: clen,crecs
	integer, allocatable, save :: dimids(:),dimlens(:)
	integer, save :: icall = 0
	double precision atime
	character(len=:), allocatable, save :: ctimes(:)
	character*20 string

	bdebug = .true.
	bdebug = .false.

	call nc_get_time_name(time_d,time_v)
	time = time_v

	call nc_get_var_id(ncid,time,time_id)
	if( time_id .eq. 0 ) then
	  stop 'error stop handle_time_string: cannot find time variable'
	end if

	retval = nf_inq_vartype(ncid,time_id,xtype)
	call nc_handle_err(retval,'get_time_rec')
	if( xtype /= NF_CHAR ) then
	  write(6,*) 'can only handle strings for time: ',xtype
	  stop 'error stop handle_time_string: only time strings'
	end if

	retval =  NF_INQ_VARNDIMS(ncid,time_id,ndims)
	call nc_handle_err(retval,'get_time_rec')
	if( ndims /= 2 ) then
	  write(6,*) 'can only handle time strings with 2 dimensions:', ndims
	  stop 'error stop handle_time_string: ndims /= 2'
	end if

	if( icall == 0 ) then
	  allocate(dimids(ndims),dimlens(ndims))
	  retval = NF_INQ_VAR(ncid,time_id,name,xtype,ndims,dimids,natts)
	  call nc_handle_err(retval,'get_time_rec')
	  do i=1,ndims
	    dimid = dimids(i)
	    retval =  NF_INQ_DIMLEN(NCID,DIMID,len)
	    dimlens(i) = len
	  end do
	  clen = dimlens(1)
	  crecs = dimlens(2)
	  allocate(character(LEN=clen) ::  ctimes(crecs))

	  retval = nf_get_var_text(ncid,time_id,ctimes)
	  call nc_handle_err(retval,'get_time_rec')

	  ntime_recs = crecs
	  btime_is_char = .true.

	end if
	icall = icall + 1

	if( bdebug ) then
	  retval = NF_INQ_VAR(ncid,time_id,name,xtype,ndims,dimids,natts)
	  call nc_handle_err(retval,'get_time_rec')
	  write(6,*) 'time info: '
	  write(6,*) time_id
	  write(6,*) trim(name)
	  write(6,*) xtype
	  write(6,*) ndims
	  write(6,*) dimids
	  do i=1,ndims
	    dimid = dimids(i)
	    dimlen = dimlens(i)
	    write(6,*) 'dimension: ',i,dimid,dimlen
	  end do
	  write(6,*) natts
	  write(6,*) 'time records found: ',crecs
	  do i=1,crecs
	    write(6,*) 'ctimes: ',i,ctimes(i)
	  end do
	end if

	t = 0
	if( irec == 0 ) return

	string = ctimes(irec)
	iu = index(string,'_')
	it = index(string,'T')
	if( iu > 0 ) then
	  string = string(:iu-1) // '::' // string(iu+1:)
	else if( it > 0 ) then
	  string = string(:it-1) // '::' // string(it+1:)
	end if

	if( bdebug ) then
	  write(6,*) 'time string: ',string
	end if
	
	call dts_string2time(string,atime,ierr)
	if( ierr /= 0 ) then
	  write(6,*) 'cannot parse string: ',string
	  stop 'error stop handle_time_string: cannot parse string'
	end if

	t = atime

	!stop 'forced stop in handle_time_string'
	end

!*****************************************************************

	subroutine nc_get_time_recs(ncid,trecs)

	use netcdf_params

	implicit none

	integer ncid
	integer trecs

	integer dim_id,len
	character*80 name
	character*80 time,time_d,time_v
	integer retval
	double precision t

	trecs = 0

	if( ntime_recs > 0 ) then	!we already know
	  trecs = ntime_recs
	  return
	end if

	if( btime_is_char ) then
	  call handle_time_string(ncid,0,t)	!only to set up data structure
	  trecs = ntime_recs
	  return
	end if

	call nc_get_time_name(time_d,time_v)
	if( time_d == ' ' .or. time_v == ' ' ) return

	call nc_get_dim_id(ncid,time_d,dim_id)

	if( dim_id .gt. 0 ) then
	  retval = nf_inq_dim(ncid,dim_id,time_d,len)
	  call nc_handle_err(retval,'get_time_recs')
	  ntime_recs = len
	  trecs = len
	end if

	end

!*****************************************************************
!*****************************************************************
! read variables
!*****************************************************************
!*****************************************************************

	subroutine nc_var_info(ncid,var_id,bverb)

	use netcdf_params

	implicit none

	integer ncid
	integer var_id
	logical bverb

	integer nvars,i,j,ia
	integer type,ndims,natts,xtype,length
	integer dimids(1024)
	double precision avalue
	character*256 name,aname,atext
	character*80 description
	integer retval

	logical, save :: blong = .false.

	blong = bverb

	  retval = nf_inq_var(ncid,var_id,name,type,ndims,dimids,natts)
	  if( bverb ) write(6,*) ncid,var_id,type,ndims,natts,trim(name)
	  call nc_handle_err(retval,'get_var_info')
	  if( ndims .gt. 10 ) stop 'error stop nc_var_info: ndims'
	  call nc_get_var_description(ncid,var_id,description)
	  write(6,1010) var_id,natts,ndims,'   ' &
     &			,trim(name),' / ',trim(description)
 1010     format(3i6,4a)

	  if( blong ) then
	    do ia=1,natts
	      retval = nf_inq_attname(ncid,var_id,ia,aname)
	      call nc_handle_err(retval,'get_var_info')
	      retval = nf_inq_att(ncid,var_id,aname,xtype,length)
	      if( nc_has_err(retval) ) cycle	!no such attribute name
	      atext = ' '
	      call nc_get_var_attrib(ncid,var_id,aname,atext,avalue)
	      if( xtype .ne. NF_CHAR ) then	!attribute is not a string
		write(atext,*) avalue
	      end if
	      write(6,1000) '    ',ia,'  ',trim(aname),'  ',trim(atext)
 1000	      format(a,i5,a,a,a,a)
	    end do
	  end if

	end

!*****************************************************************

	subroutine nc_vars_info(ncid,bverb)

	use netcdf_params

	implicit none

	integer ncid
	logical bverb

	integer nvars,var_id,i
	integer retval

	retval = nf_inq_nvars(ncid,nvars)
	call nc_handle_err(retval,'vars_info')

	write(6,*) 'variables: '
	write(6,*) '   id natts  ndim   name / description'
	do i=1,nvars
	  var_id = i
	  call nc_var_info(ncid,var_id,bverb)
	end do

	end

!*****************************************************************

	subroutine nc_get_var_id(ncid,name,var_id)

! returns var_id = 0 if not found (no error)

	use netcdf_params

	implicit none

	integer ncid
	character*(*) name
	integer var_id

	integer retval

	retval = nf_inq_varid(ncid,name,var_id)
	if( retval .ne. 0 ) then
	  var_id = 0
	  return
	end if
	call nc_handle_err(retval,'get_var_id')

	end

!*****************************************************************

	subroutine nc_get_var_description(ncid,var_id,description)

	use netcdf_params

	implicit none

	integer ncid
	integer var_id
	character*(*) description

	integer retval
	double precision avalue
	integer xtype,length
	character*256 aname,atext

	description = ' '

	aname = 'standard_name'
	retval = nf_inq_att(ncid,var_id,aname,xtype,length)
	if( .not. nc_has_err(retval) ) then
	  call nc_get_var_attrib(ncid,var_id,aname,atext,avalue)
	  description = atext
	end if

	if( description /= ' ' .and. description /= 'unknown' ) return

	aname = 'long_name'
	retval = nf_inq_att(ncid,var_id,aname,xtype,length)
	if( .not. nc_has_err(retval) ) then
	  call nc_get_var_attrib(ncid,var_id,aname,atext,avalue)
	  description = atext
	end if

	end

!*****************************************************************

	subroutine nc_get_var_name(ncid,var_id,name)

	use netcdf_params

	implicit none

	integer ncid
	integer var_id
	character*(*) name

	integer retval

	retval = nf_inq_varname(ncid,var_id,name)
	call nc_handle_err(retval,'get_var_name')

	end

!*****************************************************************

	subroutine nc_get_var_totnum(ncid,nvars)

	use netcdf_params

	implicit none

	integer ncid
	integer nvars

	integer retval

	retval = nf_inq_nvars(ncid,nvars)
	call nc_handle_err(retval,'get_var_totnum')

	end

!*****************************************************************

	subroutine nc_get_var_ndims(ncid,var_id,ndims,dimids)

! if ndims == 0 -> just compute total number of dimensions
!
! if ndims > 0 and smaller than total number of dimensions -> negative return value

	use netcdf_params

	implicit none

	integer ncid
	integer var_id
	integer ndims	!on entry dimension, on return number of dimensions
	integer dimids(ndims)

	integer retval
	integer ndim

	ndim = ndims

	retval = nf_inq_varndims(ncid,var_id,ndims)
	call nc_handle_err(retval,'get_var_ndims')

	if( ndim == 0 ) return
	if( ndims > ndim ) then
	  ndims = -ndims
	  return
	end if

	retval = nf_inq_vardimid(ncid,var_id,dimids)
	call nc_handle_err(retval,'get_var_ndims')

	end

!*****************************************************************

	subroutine nc_check_var_type(ncid,var_id,type)

	use netcdf_params

	implicit none

	integer ncid
	integer var_id
	character*(*) type

	integer xtype
	integer retval

	retval = nf_inq_vartype(ncid,var_id,xtype)
	call nc_handle_err(retval,'get_var_type')

	if( type .eq. 'integer' ) then
	  if( xtype .ne. NF_INT ) goto 99
	else if( type .eq. 'real' ) then
	  if( xtype .ne. NF_REAL ) goto 99
	else if( type .eq. 'double' ) then
	  if( xtype .ne. NF_DOUBLE ) goto 99
	else
	  write(6,*) 'type: ',type
	  stop 'error stop nc_check_type: cannot check type'
	end if

	return
   99	continue
	write(6,*) 'type: ',type
	write(6,*) 'xtype: ',xtype
	stop 'error stop nc_check_type: type mismatch'
	end

!*****************************************************************

	subroutine nc_get_global_attr(ncid,aname,atext)

	use netcdf_params

	implicit none

	integer ncid
	character*(*) aname
	character*(*) atext

	call nc_get_var_attr(ncid,NF_GLOBAL,aname,atext)

	end

!*****************************************************************

	subroutine nc_get_var_attr(ncid,var_id,aname,atext)

! in this routine is a gfortran bug !!!FIXME

	use netcdf_params

	implicit none

	integer ncid
	integer var_id
	character*(*) aname
	character*(*) atext

	logical debug
	integer retval
	integer xtype,ll,cl,sl
	integer, parameter :: cmax = 1000
	character*(cmax) aux
	!character(len=:), allocatable :: aux

	!debug = .true.		!GGU
	!debug = .false.	!GGU
	debug = bdebug_nc

	atext = ' '
	retval = nf_inq_att(ncid,var_id,aname,xtype,ll)
	if( nc_has_err(retval) ) return	!no such attribute name
	if( xtype .ne. NF_CHAR ) return		!attribute is not a string

	if( ll > cmax ) then
	  write(6,*) 'll = ',ll,'  cmax = ',cmax
	  write(6,*) 'variable is too short to receive string'
	  stop 'error stop nc_get_var_attr: ll > cmax'
	end if

	cl = ll + 1
	!allocate(character(len=cl) :: aux)
	aux = ' '
	sl = len(aux)

	if( debug ) write(6,*) 'cl,sl: ',cl,sl,len(aux)

	if( sl < cl ) then
	  write(6,*) 'sl,cl: ',sl,cl
	  write(6,*) 'compiler cannot allocate variable string length'
	  write(6,*) 'please use constant string variable'
	  stop 'error stop nc_get_var_attr: len(aux)<cl'
	end if

	if( debug ) then
	write(6,*) 'cl: ',sl,len(aux),cl
	write(6,*) 'nnn: ',ll,var_id,trim(aname),len(aux)
	end if

	retval = nf_get_att_text(ncid,var_id,aname,aux)

	sl = len(aux)
	if( ichar(aux(ll:ll)) == 0 ) aux(ll:ll) = ' '
	if( ichar(aux(sl:sl)) == 0 ) aux(sl:sl) = ' '
	if( debug ) then
	  write(6,*) 'auxaux1: ',trim(aux)
	  write(6,*) 'auxaux2: ',len_trim(aux),ll,sl
	  write(6,*) 'auxaux3: ',ichar(aux(ll:ll))
	  write(6,*) 'nnn: end of call'
	end if
	call nc_handle_err(retval,'get_var_attr')
	atext = aux

	!deallocate(aux)

	end

!*****************************************************************

	function nc_has_var_attrib(ncid,var_id,aname)

	use netcdf_params

	implicit none

	logical nc_has_var_attrib
	integer ncid
	integer var_id
	character*(*) aname

	integer retval
	integer xtype,len

	nc_has_var_attrib = .false.

	retval = nf_inq_att(ncid,var_id,aname,xtype,len)
	if( nc_has_err(retval) ) return	!no such attribute name

	nc_has_var_attrib = .true.

	end function nc_has_var_attrib

!*****************************************************************

	subroutine nc_get_var_attrib(ncid,var_id,aname,atext,avalue)

	use netcdf_params

	implicit none

	integer ncid
	integer var_id
	character*(*) aname
	character*(*) atext
	double precision avalue

	integer retval
	integer xtype,len
	double precision, allocatable :: daux(:)
	character(len=:), allocatable :: saux		!NEMUNAS_FIX_NEW
	!character*1000 :: saux				!NEMUNAS_FIX_OLD

	atext = ' '
	avalue = 0.

	retval = nf_inq_att(ncid,var_id,aname,xtype,len)
	if( nc_has_err(retval) ) return	!no such attribute name

	if( xtype .eq. 12 ) then
	  write(6,*) 'this is a netcdf-4 which contains string values'
	  write(6,*) 'this file cannot be read by the present routines'
	  write(6,*) 'please convert this file to netcdf-3 with ncks:'
	  write(6,*) '   ncks -3 nc4-file.nc nc3-file.nc'
	  write(6,*) 'you might have to install the nco package'
	  stop 'error stop nc_get_var_attrib: string format unsupported'
	end if

	if( xtype .eq. NF_CHAR ) then
	  allocate(character(len=len) :: saux)		!NEMUNAS_FIX_NEW
	  !if( len > 1000 ) goto 99			!NEMUNAS_FIX_OLD
	  retval = nf_get_att_text(ncid,var_id,aname,saux)
	  if( len == 1 .and. ichar(saux(1:1)) == 0 ) saux=' '	!FIX
	  atext = saux
	else
	  if( len > 1 ) then
	    write(6,*) '*** attributes with multiple values'
	    write(6,*) '*** only returning first of array'
	    write(6,*) '   name,len: ',trim(aname),len
	    allocate(daux(len))
	    retval = nf_get_att_double(ncid,var_id,aname,daux)
	    avalue = daux(1)
	  else
	    retval = nf_get_att_double(ncid,var_id,aname,avalue)
	  end if
	end if

	call nc_handle_err(retval,'get_var_attrib')

	return
   99	continue
	write(6,*) 'len = ',len,'  max possible is 1000'
	stop 'error stop nc_get_var_attrib: len > 1000'
	end

!*****************************************************************

	subroutine nc_get_var_int(ncid,var_id,data)

	use netcdf_params

	implicit none

	integer ncid
	integer var_id
	integer data(*)

	integer retval

	retval = nf_get_var_int(ncid,var_id,data)
	call nc_handle_err(retval,'get_var_int')

	end

!*****************************************************************

	subroutine nc_get_var_real(ncid,var_id,data)

	use netcdf_params

	implicit none

	integer ncid
	integer var_id
	real data(*)

	integer retval

	retval = nf_get_var_real(ncid,var_id,data)
	call nc_handle_err(retval,'get_var_real')

	end

!*****************************************************************

        subroutine nc_get_var_data(ncid,name,trec,ndim,ndimens          &
     &                                  ,dims,data)

! reads time record trec of variable name

	use netcdf_params

	implicit none

	integer ncid
	character*(*) name	!name of variable to read
	integer trec		!number of time record to read
	integer ndim		!dimension of data array
	integer ndimens		!expected dimensionality of data to read
				!(this should exclude time dimension)
	integer dims(ndimens+1)	!length of dimensions (return)
	real data(ndim)		!data (return)

	integer retval
	integer i,dim_id,dim_len
	integer var_id,itime
	integer ndims,nlength
	integer xtype,size
	integer, allocatable :: icount(:)
	integer, allocatable :: istart(:)
	integer, allocatable :: dimids(:)
	character*80, allocatable :: dimn(:)
	character*80 dimname
	character*80 auxname
	character*80 time,time_d,time_v

	double precision, allocatable :: ddata(:)

	call nc_get_time_name(time_d,time_v)
	time = time_d
	!write(6,*) 'name of time dim: ',time

	ndims = 0
	call nc_get_var_id(ncid,name,var_id)
	if( var_id .le. 0 ) then
	  write(6,*) 'Cannot find variable: ',name
	  return
	end if

	retval = nf_inq_varndims(ncid,var_id,ndims)
	call nc_handle_err(retval,'get_var_data')
	if( ndims > ndimens+1 ) goto 95
	allocate(icount(ndims),istart(ndims),dimids(ndims),dimn(ndims))

	retval = nf_inq_vardimid(ncid,var_id,dimids)
	call nc_handle_err(retval,'get_var_data')

	itime = 0
	nlength = 1
	do i=1,ndims
	  dim_id = dimids(i)
	  call nc_get_dim_name(ncid,dim_id,dimname)
	  call nc_get_dim_len(ncid,dim_id,dim_len)
	  dims(i) = dim_len
	  dimn(i) = dimname
	  !write(6,*) 'name of dim: ',i,dimname(1:30)
	  if( dimname .eq. time ) then
	    itime = i
	  else
	    istart(i) = 1
	    icount(i) = dim_len
	    nlength = nlength * dim_len
	  end if
	  !write(6,*) i,dim_id,dim_len,dimname(1:20)
	end do

	if( itime .eq. 0 ) then
	  !write(6,*) 'variable without time... ',itime
	else
	  !write(6,*) 'variable has time... ',itime
	  if( itime .ne. ndims ) goto 99
	  ndims = ndims - 1		! time is always last
	  istart(itime) = trec
	  icount(itime) = 1
	end if

	if( ndims .ne. ndimens ) then
	  write(6,*) 'expected dimension of data: ',ndimens
	  write(6,*) 'real dimension of data: ',ndims
	  do i=1,ndims
	    write(6,*) i,dims(i),'  ',trim(dimn(i))
	  end do
	  stop 'error stop nc_get_var_data: ndimens'
	end if

	if( nlength .gt. ndim ) then
	  write(6,*) 'size of data to be read: ',nlength
	  write(6,*) 'dimension of data array: ',ndim
	  stop 'error stop nc_get_var_data: ndim'
	end if

	retval = nf_inq_vartype(ncid,var_id,xtype)
	call nc_handle_err(retval,'get_var_data')

	if( xtype .eq. NF_FLOAT ) then
	  retval = nf_get_vara_real(ncid,var_id,istart,icount,data)
	  call nc_handle_err(retval,'get_var_data')
	else if( xtype .eq. NF_SHORT ) then
	  retval = nf_get_vara_real(ncid,var_id,istart,icount,data)
	  call nc_handle_err(retval,'get_var_data')
	else if( xtype .eq. NF_DOUBLE ) then
	  allocate(ddata(nlength))
	  retval = nf_get_vara_double(ncid,var_id,istart,icount,ddata)
	  call nc_handle_err(retval,'get_var_data')
	  data = ddata
	  deallocate(ddata)
	else
	  write(6,*) 'cannot handle data type...',xtype
	  stop 'error stop nc_get_var_data: data type'
	end if

	return
   95	continue
	write(6,*) 'dimensions to small: ',ndims,ndimens+1
	stop 'error stop nc_get_var_data: ndimens'
   99	continue
	write(6,*) 'time variable is not last: ',itime,ndims
	stop 'error stop nc_get_var_data: itime'
	end

!*****************************************************************
!*****************************************************************
! define variables for write
!*****************************************************************
!*****************************************************************

	subroutine nc_define_2d_reg(ncid,what,var_id)

	use netcdf_params

	implicit none

	integer ncid
	character*(*) what
	integer var_id				!return

	integer retval

        retval = nf_def_var(ncid, what, NF_REAL, 3, dimids_2d           &
     &                          ,var_id)
	call nc_handle_err(retval,'define_2d_reg')

	end

!*****************************************************************

	subroutine nc_define_3d_reg(ncid,what,var_id)

	use netcdf_params

	implicit none

	integer ncid
	character*(*) what
	integer var_id				!return

	integer retval

        retval = nf_def_var(ncid, what, NF_REAL, 4, dimids_3d           &
     &                          ,var_id)
	call nc_handle_err(retval,'define_3d_reg')

	end

!*****************************************************************

	subroutine nc_define_2d(ncid,what,var_id)

	use netcdf_params

	implicit none

	integer ncid
	character*(*) what
	integer var_id				!return

	integer retval

        retval = nf_def_var(ncid, what, NF_REAL, 2, dimids_2d           &
     &                          ,var_id)
	call nc_handle_err(retval,'define_2d')

	end

!*****************************************************************

	subroutine nc_define_3d(ncid,what,var_id)

	use netcdf_params

	implicit none

	integer ncid
	character*(*) what
	integer var_id				!return

	integer retval

        retval = nf_def_var(ncid, what, NF_REAL, 3, dimids_3d           &
     &                          ,var_id)
	call nc_handle_err(retval,'define_3d')

	end

!*****************************************************************
!*****************************************************************
! define attributes for write
!*****************************************************************
!*****************************************************************

	subroutine nc_define_attr(ncid,what,def,var_id)

	use netcdf_params

	implicit none

	integer ncid
	character*(*) what,def
	integer var_id

	integer ldef
	integer retval

	integer nc_ichanm

        ldef = nc_ichanm(def)
        retval = nf_put_att_text(ncid, var_id, what, ldef               &
     &                          ,def)
	call nc_handle_err(retval,'define_attr')

	end

!*****************************************************************

	subroutine nc_define_attr_real(ncid,what,value,var_id)

	use netcdf_params

	implicit none

	integer ncid
	character*(*) what
	real value
	integer var_id

	integer ldef
	integer retval

        ldef = 1
        retval = nf_put_att_real(ncid, var_id, what, NF_FLOAT, ldef     &
     &                          ,value)
	call nc_handle_err(retval,'define_attr_real')

	end

!*****************************************************************

	subroutine nc_define_attr_int(ncid,what,value,var_id)

	use netcdf_params

	implicit none

	integer ncid
	character*(*) what
	integer value
	integer var_id

	integer ldef
	integer retval

        ldef = 1
        retval = nf_put_att_int(ncid, var_id, what, NF_INT, ldef     &
     &                          ,value)
	call nc_handle_err(retval,'define_attr_int')

	end

!*****************************************************************

	subroutine nc_define_range(ncid,rmin,rmax,flag,var_id)

	use netcdf_params

	implicit none

	integer ncid
	real rmin,rmax,flag
	integer var_id

	real rminmax(2)
	integer retval

	rminmax(1) = rmin
	rminmax(2) = rmax

	if( rmin < rmax ) then

          retval = nf_put_att_real(ncid, var_id, 'valid_min', NF_REAL     &
     &                          ,1,rmin)
	  call nc_handle_err(retval,'define_range')

          retval = nf_put_att_real(ncid, var_id, 'valid_max', NF_REAL     &
     &                          ,1,rmax)
	  call nc_handle_err(retval,'define_range')

	end if

        retval = nf_put_att_real(ncid, var_id, '_FillValue', NF_REAL    &
     &                          ,1,flag)
	call nc_handle_err(retval,'define_range')

	end

!*****************************************************************

	subroutine nc_end_define(ncid)

	use netcdf_params

	implicit none

	integer ncid

	integer retval

	retval = nf_enddef(ncid)
	call nc_handle_err(retval,'end_define')

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine nc_has_vertical_dimension(ncid,name,bvert)

! checks if variable is 3d

	use netcdf_params

	implicit none

	integer ncid
	character*(*) name
	logical bvert

	logical btime
	integer var_id,ndims
	integer retval
	integer, allocatable :: dim_ids(:)

	call nc_get_var_id(ncid,name,var_id)

        retval = nf_inq_varndims(ncid,var_id,ndims)
        call nc_handle_err(retval,'has_vertical_dimension')
	allocate(dim_ids(ndims))

	call nc_get_var_ndims(ncid,var_id,ndims,dim_ids)

	call nc_has_time_dimension(ncid,name,btime)
	if( btime ) ndims = ndims - 1

	bvert = ( ndims == 3 )

	end

!*****************************************************************

	subroutine nc_has_time_dimension(ncid,name,btime)

! checks if variable name has time dimension

	use netcdf_params

	implicit none

	integer ncid
	character*(*) name
	logical btime

	integer var_id,ndims,time_id
	integer retval
	integer, allocatable :: dim_ids(:)
	character*80 tname,time_d,time_v

	call nc_get_var_id(ncid,name,var_id)
        retval = nf_inq_varndims(ncid,var_id,ndims)
        call nc_handle_err(retval,'has_time_dimension')
        allocate(dim_ids(ndims))

	call nc_get_var_ndims(ncid,var_id,ndims,dim_ids)

	time_id = dim_ids(ndims)			!time is always last
        call nc_get_dim_name(ncid,time_id,tname)
	call nc_get_time_name(time_d,time_v)

	btime = tname .eq. time_d

	end

!*****************************************************************

	subroutine nc_set_time_name(time_d,time_v)

	use netcdf_params

	implicit none

	character*(*) time_d,time_v

	if( time_d .ne. ' ' )  time_dim  = time_d	!dimension name of time
	if( time_v .ne. ' ' )  time_var  = time_v	!variable name of time

	end

!*****************************************************************

	subroutine nc_get_time_name(time_d,time_v)

	use netcdf_params

	implicit none

	character*(*) time_d,time_v

	time_d = time_dim
	time_v = time_var

	end

!*****************************************************************
!*****************************************************************
! write coordinate information
!*****************************************************************
!*****************************************************************

        subroutine nc_write_coords_reg(ncid,nx,ny,nlv                   &
     &                                  ,xlon,ylat,depth,hlv)

	use netcdf_params

	implicit none

	integer ncid
	integer nx,ny
	integer nlv
	real xlon(nx)
	real ylat(ny)
	real depth(nx,ny)
	real hlv(nlv)

	integer lon_varid,lat_varid,lvl_varid,dep_varid
	integer eix_varid,top_varid
	integer retval

	lon_varid = coord_varid(1)
	lat_varid = coord_varid(2)
	lvl_varid = coord_varid(3)
	dep_varid = coord_varid(4)

!-----------------------------------------
! write coordinate data
!-----------------------------------------

	retval = nf_put_var_real(ncid, lon_varid, xlon)
	call nc_handle_err(retval,'write_coords_reg')

	retval = nf_put_var_real(ncid, lat_varid, ylat)
	call nc_handle_err(retval,'write_coords_reg')

	retval = nf_put_var_real(ncid, lvl_varid, hlv)
	call nc_handle_err(retval,'write_coords_reg')

	retval = nf_put_var_real(ncid, dep_varid, depth)
	call nc_handle_err(retval,'write_coords_reg')

!-----------------------------------------
! end of routine
!-----------------------------------------

	end

!*****************************************************************

        subroutine nc_write_coords_fem(ncid,nkn,nel,nlv                 &
     &                          ,xgv,ygv,hkv,nen3v,hlv)

	!use mod_depth
	!use levels
	!use basin
	use netcdf_params

	implicit none

	integer ncid
	integer nkn,nel,nlv
	real xgv(nkn),ygv(nkn)
	real hkv(nkn)
	integer nen3v(3,nel)
	real hlv(nlv)

	integer lon_varid,lat_varid,lvl_varid,dep_varid
	integer eix_varid,top_varid
	integer retval

	lon_varid = coord_varid(1)
	lat_varid = coord_varid(2)
	lvl_varid = coord_varid(3)
	dep_varid = coord_varid(4)
	eix_varid = coord_varid(5)
	top_varid = coord_varid(6)

!-----------------------------------------
! write coordinate data
!-----------------------------------------

	retval = nf_put_var_real(ncid, lon_varid, xgv)
	call nc_handle_err(retval,'write_coords_fem')

	retval = nf_put_var_real(ncid, lat_varid, ygv)
	call nc_handle_err(retval,'write_coords_fem')

	retval = nf_put_var_real(ncid, lvl_varid, hlv)
	call nc_handle_err(retval,'write_coords_fem')

	retval = nf_put_var_real(ncid, dep_varid, hkv)
	call nc_handle_err(retval,'write_coords_fem')

	retval = nf_put_var_int(ncid, eix_varid, nen3v)
	call nc_handle_err(retval,'write_coords_fem')

	retval = nf_put_var_int(ncid, top_varid, 2)
	call nc_handle_err(retval,'write_coords_fem')

!-----------------------------------------
! end of routine
!-----------------------------------------

	end

!*****************************************************************

	subroutine nc_write_coords_ts(ncid,lon,lat)

	use netcdf_params

	implicit none

	integer ncid
        real lon(1), lat(1)

	integer lon_varid,lat_varid
	integer retval

	lon_varid = coord_varid(1)
	lat_varid = coord_varid(2)

!-----------------------------------------
! write coordinate data
!-----------------------------------------

	retval = nf_put_var_real(ncid, lon_varid, lon)
	call nc_handle_err(retval,'write_coords_ts')

	retval = nf_put_var_real(ncid, lat_varid, lat)
	call nc_handle_err(retval,'write_coords_ts')

!-----------------------------------------
! end of routine
!-----------------------------------------

	end

!*****************************************************************
!*****************************************************************
! write time record
!*****************************************************************
!*****************************************************************

	subroutine nc_write_time(ncid,irec,it)

	use netcdf_params

	implicit none

	integer ncid
	integer irec
	integer it

	integer retval
	double precision dtime

	dtime = it
	call nc_write_dtime(ncid,irec,dtime)

	end

!*****************************************************************

	subroutine nc_write_dtime(ncid,irec,dtime)

	use netcdf_params

	implicit none

	integer ncid
	integer irec
	double precision dtime

	integer retval
	double precision ctime

	ctime = dtime + dtime0		!convert to other reference
	retval = nf_put_vara_double(ncid, rec_varid, irec, 1, ctime)
	call nc_handle_err(retval,'write_dtime')

	end

!*****************************************************************
!*****************************************************************
! write data
!*****************************************************************
!*****************************************************************

	subroutine nc_write_data_2d_reg(ncid,var_id,irec,nx,ny,var2d)

	use netcdf_params

	implicit none

	integer ncid
	integer var_id
	integer irec
	integer nx,ny
	real var2d(nx,ny)

	integer retval
	integer count(3)
	integer start(3)

	count(1) = nx
	count(2) = ny
	count(3) = 1
	start(1) = 1
	start(2) = 1
	start(3) = irec

	retval = nf_put_vara_real(ncid, var_id, start, count, var2d)
	call nc_handle_err(retval,'write_data_2d_reg')

	end

!*****************************************************************

	subroutine nc_write_data_3d_reg(ncid,var_id,irec,nlv,nx,ny,var3d)

	use netcdf_params

	implicit none

	integer ncid
	integer var_id
	integer irec
	integer nlv
	integer nx,ny
	!real var3d(nlv,nx,ny)
	real var3d(nx,ny,nlv)

	integer retval
	integer count(4)
	integer start(4)

	count(1) = nx
	count(2) = ny
	count(3) = nlv
	count(4) = 1
	start(1) = 1
	start(2) = 1
	start(3) = 1
	start(4) = irec

	retval = nf_put_vara_real(ncid, var_id, start, count, var3d)
	call nc_handle_err(retval,'write_data_3d_reg')

	end

!*****************************************************************

	subroutine nc_rewrite_3d_reg(nlvddi,lmax,nx,ny,var3d,vnc3d)

! re-writes a 3d array to be CF compliant

	use netcdf_params

	implicit none

	integer nlvddi,lmax
	integer nx,ny
	real var3d(nlvddi,nx,ny)
	real vnc3d(nx,ny,lmax)

	integer i,j,l

	do j=1,ny
	  do i=1,nx
	    do l=1,lmax
	      vnc3d(i,j,l) = var3d(l,i,j)
	    end do
	  end do
	end do

	end

!*****************************************************************

	subroutine nc_write_data_2d(ncid,var_id,irec,nkn,var2d)

	use netcdf_params

	implicit none

	integer ncid
	integer var_id
	integer irec
	integer nkn
	real var2d(nkn)

	integer retval
	integer count(2)
	integer start(2)

	count(1) = nkn
	count(2) = 1
	start(1) = 1
	start(2) = irec

	retval = nf_put_vara_real(ncid, var_id, start, count, var2d)
	call nc_handle_err(retval,'write_data_2d')

	end

!*****************************************************************

	subroutine nc_write_data_3d(ncid,var_id,irec,nlv,nkn,var3d)

	use netcdf_params

	implicit none

	integer ncid
	integer var_id
	integer irec
	integer nlv
	integer nkn
	real var3d(nlv,nkn)

	integer retval
	integer count(3)
	integer start(3)

	count(1) = nlv
	count(2) = nkn
	count(3) = 1
	start(1) = 1
	start(2) = 1
	start(3) = irec

	retval = nf_put_vara_real(ncid, var_id, start, count, var3d)
	call nc_handle_err(retval,'write_data_3d')

	end

!*****************************************************************
!*****************************************************************
! compact data for write
!*****************************************************************
!*****************************************************************

	subroutine nc_compact_3d_reg(nlvddi,nlv,nx,ny,var_in,var_out)

	implicit none

	integer nlvddi
	integer nlv,nx,ny
	real var_in(nlvddi,nx,ny)
	real var_out(nlv,nx,ny)

	integer l,ix,iy

	do iy=1,ny
	 do ix=1,nx
	  do l=1,nlv
	    var_out(l,ix,iy) = var_in(l,ix,iy)
	  end do
	 end do
	end do

	end

!*****************************************************************

	subroutine nc_compact_3d(nlvddi,nlv,nkn,var_in,var_out)

	implicit none

	integer nlvddi
	integer nlv,nkn
	real var_in(nlvddi,nkn)
	real var_out(nlv,nkn)

	integer k,l

	do k=1,nkn
	  do l=1,nlv
	    var_out(l,k) = var_in(l,k)
	  end do
	end do

	end

!*****************************************************************
!*****************************************************************
! various routines
!*****************************************************************
!*****************************************************************

	subroutine nc_close(ncid)

	use netcdf_params

	implicit none

	integer ncid

	integer retval

	retval = nf_close(ncid)
	call nc_handle_err(retval,'close')

	end

!*****************************************************************

	subroutine nc_global(ncid,title,bugrid,sncglobal)

! writes global conventions

	use netcdf_params

	implicit none

	integer ncid
	character*(*) title
	logical bugrid
	character*(*) sncglobal

	integer ltext,retval,varid
	character*80 text
	character*80 what
	character*80 cdate

	integer nc_ichanm

	varid = NF_GLOBAL

	what = 'Conventions'
	text = cf_convention
	if( bugrid ) text = trim(cf_convention)//','//trim(ugrid_convention)
	call nc_define_attr(ncid,what,text,varid)

	what = 'title'
	text = title
	call nc_strip(text)
	call nc_define_attr(ncid,what,text,varid)

	call nc_current_time(cdate)
	what = 'history'
	text = 'created on ' // trim(cdate)
	call nc_strip(text)
	call nc_define_attr(ncid,what,text,varid)

	what = 'institution'
	text = 'ISMAR-CNR, Venice, Italy'
	call nc_define_attr(ncid,what,text,varid)

	what = 'source'
	text = 'Model data produced by SHYFEM at ISMAR-CNR'
	call nc_define_attr(ncid,what,text,varid)

	what = 'references'
	text = 'Model info: http://www.ismar.cnr.it/shyfem'
	call nc_define_attr(ncid,what,text,varid)

	what = 'contact'
	text = 'email: georg.umgiesser@ismar.cnr.it'
	call nc_define_attr(ncid,what,text,varid)

	what = 'comment'
	text = 'Data restriction: for academic research use only'
	call nc_define_attr(ncid,what,text,varid)

	call nc_user_defined_global(ncid,sncglobal)

	end

!*****************************************************************

	subroutine nc_user_defined_global(ncid,sncglobal)

! user defined global comments

	use netcdf_params

	implicit none

	integer ncid
	character*(*) sncglobal

	integer iunit,ios,i,varid
	character*80 myfile
	character*256 string
	character*256 key,text

	iunit = 3
	myfile = sncglobal
	if( myfile == ' ' ) return

	varid = NF_GLOBAL

	write(6,*) 'ncglobal file: ',trim(myfile)

	open(iunit,file=myfile,form='formatted',status='old',iostat=ios)

	if( ios /= 0 ) then
	  write(6,*) 'no such file: ',trim(myfile)
	  stop 'error stop nc_user_defined_global: no such file'
	end if

	do
	  read(iunit,'(a)',iostat=ios) string
	  if( ios /= 0 ) exit
	  string =  adjustl(string)
	  if( string == ' ' ) cycle		! empty line
	  i = index(string,':')
	  if( i == 0 ) goto 99			! no ":" found
	  key = string(1:i-1)
	  key = trim(key)
	  text = string(i+1:)
	  text =  adjustl(text)	
	  text = trim(text)
	  if( text == ' ' ) cycle		! no text given
	  !write(6,*) 'string: ',trim(string)
	  !write(6,*) 'key:    ',trim(key)
	  !write(6,*) 'text:   ',trim(text)
	  call nc_define_attr(ncid,key,text,varid)
	end do

	close(iunit)

	return
   99	continue
	write(6,*) 'cannot parse string: ',trim(string)
	write(6,*) 'expecting key: text'
	stop 'error stop nc_user_defined_global: cannot parse'
	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine nc_handle_err(errcode,string)

	use netcdf_params, only : nf_noerr,nf_strerror

	implicit none

	integer errcode
	character*(*), optional :: string

	if( errcode .eq. nf_noerr ) return

	if( present(string) ) write(6,*) trim(string)
	write(6,*) 'Error: ', nf_strerror(errcode)

	stop 'error stop nc_handle_err'
	end

!*****************************************************************

	function nc_has_err(errcode)

	use netcdf_params, only : nf_noerr

	implicit none

	logical nc_has_err
	integer errcode

	nc_has_err = ( errcode .ne. nf_noerr )

	end

!*****************************************************************
!*****************************************************************
! string utility routines
!*****************************************************************
!*****************************************************************

        function nc_ichanm(line)

! computes length of line without trailing blanks
!
! line          line of text
! ichanm        length of line (return value)
!               ... 0 : line is all blank

	implicit none

	integer nc_ichanm
        character*(*) line

	integer i,ndim
        character*1 blank,tab,char
        data blank /' '/

	tab = char(9)

        ndim=len(line)

        do i=ndim,1,-1
          if(line(i:i).ne.blank.and.line(i:i).ne.tab) goto 1
        end do

    1   continue
        nc_ichanm=i

        return
        end

!*****************************************************************

        subroutine nc_strip(line)

! strip blank lines
!
! line          line of text
! ichanm        length of line (return value)
!               ... 0 : line is all blank

	implicit none

        character*(*) line

	integer i,ndim
        character*1 blank,tab,char
        data blank /' '/

	tab = char(9)

        ndim=len(line)

        do i=1,ndim
          if(line(i:i).ne.blank.and.line(i:i).ne.tab) goto 1
        end do

    1   continue
	if( i .gt. ndim ) return
	if( i .eq. 1 ) return

	line(1:) = line(i:)

        return
        end

!*****************************************************************

	subroutine nc_subst_char(line,orig,subst)

! substitutes in line character orig with subst

	implicit none

        character*(*) line
        character*1 orig,subst

	integer i,n

	n = len(line)

	do i=1,n
	  if( line(i:i) .eq. orig ) line(i:i) = subst
	end do

	end

!*****************************************************************
!*****************************************************************
! date and time routines
!*****************************************************************
!*****************************************************************

	subroutine nc_unpack_date(ipack,i1,i2,i3)

! unpacks date in ipack to integers

	implicit none

	integer ipack,i1,i2,i3

	integer iaux

	i1 = 0
	i2 = 0
	i3 = 0

	iaux = ipack
	i1 = iaux/10000
	iaux = iaux - 10000*i1
	i2 = iaux/100
	iaux = iaux - 100*i2
	i3 = iaux

	end

!*****************************************************************

	subroutine nc_format_date(date,year,month,day,hour,min,sec,zone)

! formats date string

	implicit none

	character*(*) date
	integer year,month,day,hour,min,sec
	character*(*) zone	!time zone - can be blank (UTC,MET,etc.)

	date = ' '
	write(date,'(i4,a1,i2,a1,i2)') year,'-',month,'-',day
	call nc_subst_char(date(1:10),' ','0')
	write(date(12:19),'(i2,a1,i2,a1,i2)') hour,':',min,':',sec
	call nc_subst_char(date(12:19),' ','0')
	date(21:) = zone

	end

!*****************************************************************

	subroutine nc_convert_date(date0,time0,date)

! converts date from integer to character string

	use netcdf_params, only : dtime0,tref

	implicit none

	integer date0,time0
	character*(*) date

	integer aux,year,month,day
	integer hour,min,sec
	integer i
	integer ierr
	double precision atime0,atime_ref
	character*80 saux

!-----------------------------------------------
! convert date and time - dtime is relative to this date
!-----------------------------------------------

	aux = date0
	if( aux .lt. 10000 ) aux = 10000 * aux
	call nc_unpack_date(aux,year,month,day)
	if( month .le. 0 ) month = 1
	if( day .le. 0 ) day = 1
	call nc_unpack_date(time0,hour,min,sec)

	call dts_to_abs_time(aux,time0,atime0)

	!write(6,*) 'date0: ',date0,year,month,day
	!write(6,*) 'time0: ',time0,hour,min,sec

!-----------------------------------------------
! prepare date
!-----------------------------------------------

	call nc_format_date(date,year,month,day,hour,min,sec,'UTC')

!-----------------------------------------------
! see if we have to change reference date
!-----------------------------------------------

	if( tref == ' ' ) return	!no, reference date is in date0,time0

	call dts_string2time(tref,atime_ref,ierr)
	if( ierr /= 0 ) goto 99

	dtime0 = atime0 - atime_ref

!-----------------------------------------------
! convert tref to nc format
!-----------------------------------------------

	saux = tref
	i = index(saux,'::')
	if( i == 0 ) goto 99
	saux = saux(:i-1) // ' ' // saux(i+2:)
	saux = trim(saux) // ' UTC'
	date = saux

!-----------------------------------------------
! end of routine
!-----------------------------------------------

	return
   99	continue
	write(6,*) 'error converting date string:'
	write(6,*) 'date: ',trim(date)
	write(6,*) 'tref: ',trim(tref)
	stop 'error stop nc_convert_date: converting date string'
	end

!*****************************************************************

	subroutine nc_current_time(cdate)

	use netcdf_params

	implicit none

	character*(*) cdate

	integer year,month,day,hour,min,sec
	character*10 date,time,zone
	integer value(8)

	call date_and_time(date,time,zone,value)

	year = value(1)
	month = value(2)
	day = value(3)
	hour = value(5)
	min = value(6)
	sec = value(7)

	!if( year .lt. 1000 ) year = 2000 + year

	call nc_format_date(cdate,year,month,day,hour,min,sec,'MET')

	if( .not. bquiet_nc ) then
	  write(6,*) 'cdate: ',trim(cdate),'  ',date,'  ',time
	end if

	end
	
!*****************************************************************

	subroutine nc_set_ref_date(date0)

	use netcdf_params

	implicit none

	character*(*) date0
	integer ierr

	tref = date0

	end

!*****************************************************************
!*****************************************************************
! variable initialization
!*****************************************************************
!*****************************************************************

	subroutine nc_init_variable(ncid,breg,dim,ivar,idir,flag,var_id)

	use netcdf_params

	implicit none

	integer ncid
	logical breg
	integer dim
	integer ivar
	integer idir
	real flag
	integer var_id		! id to be used for other operations (return)

	character*80 name,what,std,units
	character*80 string,short
	real cmin,cmax

	if( ivar .eq. 1 ) then		! water level
	  name = 'water_level'
	  what = 'standard_name'
          std = 'water_surface_height_above_reference_datum'
          units = 'm'
	  cmin = -10.
	  cmax = +10.
	else if( ivar .eq. 2 ) then	! x-velocity
	  if( idir == 1 ) then
	    name = 'u_velocity'
	    std = 'eastward_sea_water_velocity'
	  else if( idir == 2 ) then
	    name = 'v_velocity'
	    std = 'northward_sea_water_velocity'
	  else
	    goto 99
	  end if
	  what = 'standard_name'
	  units = 'm s-1'
	  cmin = -10.
	  cmax = +10.
	else if( ivar .eq. 367 ) then	! generic tracer
	  name = 'tracer'
	  what = 'long_name'
	  std = 'river_influence_index'
	  units = '1'
	  cmin = 0.
	  cmax = 20.
	else if( ivar .eq. 10 ) then	! generic tracer
	  name = 'tracer'
	  what = 'long_name'
	  std = 'generic_tracer'
	  units = '1'
	  cmin = 0.
	  cmax = 110.
	else if( ivar .eq. 11 ) then	! salinity
	  name = 'salinity'
	  what = 'standard_name'
	  std = 'sea_water_salinity'
	  units = '1e-3'
	  cmin = 0.
	  cmax = 200.
	else if( ivar .eq. 12 ) then	! temperature
	  name = 'temperature'
	  what = 'standard_name'
	  std = 'sea_water_temperature'
	  units = 'degC'
	  cmin = -10.
	  cmax = 100.
        else if( ivar .eq. 13 ) then    ! density anomaly
          name = 'density_anomaly'
          what = 'standard_name'
          std = 'sea_water_density_anomaly'
          units = 'kg m-3'
          cmin = -1000.
          cmax = +1000.		
	else if( ivar .eq. 14 ) then	! skin temperature
	  name = 'skin_temperature'
	  what = 'standard_name'
	  std = 'sea_surface_skin_temperature'
	  units = 'degC'
	  cmin = -10.
	  cmax = 100.
        else if( ivar .eq. 17 ) then    ! vertical mass-transfer among layers
          name = 'w_velocity'
          what = 'standard_name'
          std = 'vertical_sea_water_mass_tranfer'
          units = 'm s-1'
          cmin = -10.
          cmax = +10.
	else if( ivar .eq. 19 ) then	! vorticity
	  name = 'vorticity'
	  what = 'standard_name'
	  std = 'ocean_relative_vorticity'
	  units = 's-1'
	  cmin = -100.
	  cmax =  100.
	else if( ivar .eq. 20 ) then	! atmospheric pressure
	  name = 'atmospheric_pressure'
	  what = 'standard_name'
	  std = 'air_pressure_at_mean_sea_level'
	  units = 'Pa'
	  cmin =  90000.
	  cmax = 110000.
	else if( ivar .eq. 21 ) then	! wind velocity
	  if( idir == 1 ) then
	    name = 'x_wind_velocity_at_10m'
	    std = 'eastward_wind'
	  else if( idir == 2 ) then
	    name = 'y_wind_velocity_at_10m'
	    std = 'northward_wind'
	  else
	    goto 99
	  end if
	  what = 'standard_name'
	  units = 'm s-1'
	  cmin =  -100.
	  cmax =  +100.
	else if( ivar .eq. 22 ) then	! solar radiation
	  name = 'solar_radiation'
	  what = 'standard_name'
	  std = 'solar_irradiance'
	  units = 'W m-2'
	  cmin = 0.
	  cmax = 1300.
	else if( ivar .eq. 23 ) then	! air temperature
	  name = 'air_temperature'
	  what = 'standard_name'
	  std = 'air_temperature'
	  units = 'degree Celsius'
	  cmin = -100.
	  cmax = +100.
	else if( ivar .eq. 24 ) then	! relative humidity
	  name = 'relative_humidity'
	  what = 'standard_name'
	  std = 'relative_humidity'
	  units = '%'
	  cmin = 0.
	  cmax = 100.
	else if( ivar .eq. 25 ) then	! cloud cover
	  name = 'cloud_cover'
	  what = 'standard_name'
	  std = 'cloud_area_fraction'
	  units = ''
	  cmin = 0.
	  cmax = 1.
	else if( ivar .eq. 26 ) then	! rain
	  name = 'rain'
	  what = 'standard_name'
	  std = 'rainfall_rate'
	  units = 'mm d-1'
	  cmin = 0.
	  cmax = +100.
	else if( ivar .eq. 28 ) then	! wind speed
	  name = 'wind_speed'
	  what = 'standard_name'
	  std = 'wind_speed'
	  units = 'm s-1'
	  cmin = -100.
	  cmax = +100.
	else if( ivar .eq. 99 ) then	! wrt
	  name = 'water_renewal_time'
	  what = 'long_name'
	  std = 'water_renewal_time'
	  units = '1'
	  cmin = 0.
	  cmax = 500.
	else if( ivar .eq. 81 ) then	! particle age
	  name = 'particle_age'
	  what = 'long_name'
	  std = 'lagrangian_particle_age'
	  units = 'day'
	  cmin = 0.
	  cmax = 50000.
	else if( ivar .eq. 84 ) then	! particle density
	  name = 'particle_density'
	  what = 'long_name'
	  std = 'lagrangian_particle_concentration'
	  units = 'num km-2'
	  cmin = 0.
	  cmax = 50000.
	else if( ivar .eq. 85 ) then	! ice cover
	  name = 'ice_cover'
	  what = 'standard_name'
	  std = 'sea_ice_area_fraction'
	  units = ''
	  cmin = 0.
	  cmax = 1.
	else if( ivar .eq. 86 ) then	! ice thickness
	  name = 'sea_ice_thickness'
	  what = 'long_name'
	  std = 'sea_ice_thickness'
	  units = 'm'
	  cmin = 0.
	  cmax = 100.
	else if( ivar .eq. 231 ) then	! significant wave height
	  name = 'significant_wave_height'
	  what = 'long_name'
	  std = 'sea_surface_wind_wave_significant_height'
	  units = 'm'
	  cmin = 0.
	  cmax = 100.
	else if( ivar .eq. 232 ) then	! significant wave period
	  name = 'significant_wave_period'
	  what = 'long_name'
	  std = 'sea_surface_wave_significant_period'
	  units = 's'
	  cmin = 0.
	  cmax = 100.
	else if( ivar .eq. 233 ) then	! wave direction
	  name = 'wave_direction'
	  what = 'long_name'
	  std = 'sea_surface_wave_to_direction'
	  units = 'degree'
	  cmin = 0.
	  cmax = 360.
	else if( ivar > 300 .and. ivar < 400 ) then	! concentration
	  name = 'concentration'
	  what = 'long_name'
	  std = 'concentration'		! no CF convention name available
	  units = ''
	  cmin = 0.
	  cmax = 0.
	else if( ivar .eq. 800 ) then	! suspended sediment concentration
	  name = 'ssc'
	  what = 'long_name'
	  std = 'suspended_sediment_concentration'
	  units = 'kg m-3'
	  cmin = 0.
	  cmax = 10000.
	else if( ivar .eq. 891 ) then	! sediment erosion-deposition
	  name = 'erosion-deposition'
	  what = 'long_name'
	  std = 'erosion-deposition'
	  units = 'm'
	  cmin = 0.
	  cmax = 0.
	else if( ivar .eq. 892 ) then	! sediment grainsize
	  name = 'sediment_grainsize'
	  what = 'long_name'
	  std = 'average sediment grainsize'
	  units = 'm'
	  cmin = 0.
	  cmax = 0.
	else if( ivar .eq. 893 ) then	! bottom_shear_stress
	  name = 'bottom_shear_stress'
	  what = 'long_name'
	  std = 'bottom_shear_stress'
	  units = 'm'
	  cmin = 0.
	  cmax = 0.
	else if( ivar .eq. 894 ) then	! mud_fraction
	  name = 'mud_fraction'
	  what = 'long_name'
	  std = 'mud_fraction'
	  units = ''
	  cmin = 0.
	  cmax = 1.
	else if( ivar .eq. 895 ) then	! bedload_transport
	  name = 'bedload_transport'
	  what = 'long_name'
	  std = 'bedload_transport'
	  units = 'kg/ms'
	  cmin = 0.
	  cmax = 0.
	else if( ivar .eq. 0 ) then	! no variable description
	  write(6,*) 'unknown variable: ',ivar
	  stop 'error stop nc_init_variable'
	else
	  call make_unknown_variable_name(ivar,short,string,units)
	  if( bverbose_nc ) then
	    write(6,*) 'this variable has no CL description: ',ivar
	    write(6,*) 'string: ',trim(string)
	    write(6,*) 'short: ',trim(short)
	    write(6,*) 'units: ',trim(units)
	  end if
	  name = string
	  name = short
	  what = 'long_name'
	  std = string
	  cmin = 0.
	  cmax = 0.
	end if

	if( bverbose_nc ) then
	  write(6,*) 'nc variable ',ivar,trim(name)
	end if

	if( dim .eq. 2 ) then
	  if( breg ) then
	    call nc_define_2d_reg(ncid,name,var_id)
	  else
	    call nc_define_2d(ncid,name,var_id)
	  end if
	else if( dim .eq. 3 ) then
	  if( breg ) then
	    call nc_define_3d_reg(ncid,name,var_id)
	  else
	    call nc_define_3d(ncid,name,var_id)
	  end if
	end if

	if( units /= ' ' ) call nc_define_attr(ncid,'units',units,var_id)
	call nc_define_attr(ncid,what,std,var_id)
	call nc_define_range(ncid,cmin,cmax,flag,var_id)

	!write(6,*) 'var_id: ',var_id

	return
   99	continue
	write(6,*) 'cannot parse directional flag: ',ivar,idir
	stop 'error stop nc_init_variable: internal error (1)'
	end

!*****************************************************************
!*****************************************************************
! auxiliary routines
!*****************************************************************
!*****************************************************************

	subroutine make_unknown_variable_name(ivar,short,string,units)

	use shyfem_strings

	implicit none

	integer ivar
	character*(*) short
	character*(*) string
	character*(*) units

	integer isub,i
	character*80 saux

	call strings_get_short_name(ivar,short)
        call strings_get_full_name(ivar,saux,isub)
	
	!--------------------------------
	! get rid of parenthesis
	!--------------------------------

	i = index(saux,'(')
	if( i > 0 ) saux(i:) = saux(i+1:)
	i = index(saux,')')
	if( i > 0 ) saux(i:) = saux(i+1:)

	!--------------------------------
	! extract units
	!--------------------------------

	units = ' '
	i = index(saux,'[')
	if( i > 0 ) then
	  units = saux(i+1:)
	  saux = saux(:i-1)
	  i = index(units,']')
	  if( i == 0 ) then
	    stop 'error stop make_unknown_variable_name: cannot parse unit'
	  end if
	  units = units(:i-1)
	end if
	
	!--------------------------------
	! convert blanks to underscores
	!--------------------------------

	do i=1,len_trim(saux)
	  if( saux(i:i) == ' ' ) saux(i:i) = '_'
	end do
        string = saux

	!--------------------------------
	! add substring
	!--------------------------------

        if( isub > 0 ) then
          write(saux,'(i4)') isub
	  saux = adjustl(saux)
          string = trim(string) // '_' // saux
        end if

	end

!*****************************************************************

	subroutine nc_set_verbose(bverbose)

	use netcdf_params

	implicit none

	logical bverbose

	bverbose_nc = bverbose

	end

!*****************************************************************

	subroutine nc_set_quiet(bquiet)

	use netcdf_params

	implicit none

	logical bquiet

	bquiet_nc = bquiet

	end

!*****************************************************************
!*****************************************************************
! test routines
!*****************************************************************
!*****************************************************************

	subroutine wrnetcdf
	implicit none
	end

!*****************************************************************

	function next_record(it,znv)

	implicit none

	logical next_record
	integer it
	real znv(1)

	next_record = .true.

	end

!*****************************************************************

	subroutine test_nc

! this is not working - we first have to read a basin...

!	use basin
!	use levels
!	use mod_depth

	implicit none

	integer ncid
	integer dimids_2d(2)
	integer coord_varid(5)
	integer rec_varid
	integer level_id

	integer it
	integer irec
	integer date0,time0
	integer iztype		!type of vertical coordinates
	real znv(1)
	logical next_record
	character*1 units

	date0 = 20120101
	time0 = 0
	units = 'm'
	iztype = 1

!	call nc_open_fem(ncid,nkn,nel,nlv,date0,time0,iztype)
!	call nc_define_2d(ncid,'water_level',level_id)
!	call nc_define_attr(ncid,'units',units,level_id)
!	call nc_end_define(ncid)
!	call nc_write_coords_fem(ncid,nkn,nel,nlv
!     +			,xgv,ygv,hkv,nen3v,hlv)

!	irec = 0
!	do while( next_record(it,znv) )
!
!	  irec = irec + 1
!	  call nc_write_time(ncid,irec,it)
!	  call nc_write_data_2d(ncid,level_id,irec,nkn,znv)
!
!	end do
!
!	call nc_close(ncid)

	end

!******************************************************************

!	program test_nc_main
!	implicit none
!	call test_nc
!	end

!******************************************************************

