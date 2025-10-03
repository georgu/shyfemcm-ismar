
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
! 03.03.2023    ggu     created
! 15.01.2025    ggu     adjourned
! 01.10.2025    ggu     new routine print_time() and option -time

!===================================================================
	module mod_ncinfo
!===================================================================

	logical, save :: binfo = .false.
	logical, save :: bverbose = .false.
	logical, save :: bquiet = .false.
	logical, save :: bsilent = .false.
	logical, save :: btime = .false.
	character*80, save :: variable
	character*80, save :: attribute

!===================================================================
	end module mod_ncinfo
!===================================================================

        program ncinfo

	use ncf
	use mod_ncinfo

        implicit none
        !include 'netcdf.inc'

        integer ncid

        integer nvars, natts, ngatts, idunlim
        integer nc,ia,varid,attid

	integer nt,nx,ny,nz
	character*80 tcoord,xcoord,ycoord,zcoord

	logical bwrite
        integer id
	character*80 ncfile

	type(var_item) :: vitem
	type(att_item) :: aitem
	type(dim_item) :: ditem

	type(nc_item) :: nitem

!---------------------------------------------------------------------
! open netcdf file
!---------------------------------------------------------------------

	call ncinfo_init(ncfile)

	bwrite = .not. bquiet

!---------------------------------------------------------------------
! print info on file, dimensions and global attributes
!---------------------------------------------------------------------

	call ncf_open_read(ncfile,ncid)

	nitem = ncf_get_nitem(ncid)

	ngatts = nitem%ngatts
	idunlim = nitem%idunlim

	if( binfo ) then
	write(6,*) 'global attributes: ',ngatts
	write(6,*) 'unlimited dimension : ',idunlim
	write(6,*) 'dimensions: ',nitem%ditem%ndims
	call ncf_print_dimensions(nitem%ditem)

	write(6,*) 'global attributes: ',ngatts
	call ncf_print_attribute_header(ncid,NF_GLOBAL)
	do id=1,ngatts
	  aitem = nitem%gitems(id)
	  call ncf_print_attribute(aitem)
	end do
	end if

!---------------------------------------------------------------------
! dimensions and coordinates
!---------------------------------------------------------------------

!	call ncnames_init
!	call ncnames_get_dims_and_coords(ncid,bverbose
!     +				,nt,nx,ny,nz
!     +				,tcoord,xcoord,ycoord,zcoord)

!---------------------------------------------------------------------
! print info on variables and attributes
!---------------------------------------------------------------------

	if( binfo ) then
	  nvars = nitem%nvars
	  write(6,*) 'variables: ',nvars
	  call ncf_print_variable_header(ncid,nitem)
	  do varid=1,nvars
	    call ncf_var_inf(ncid,varid,vitem)
	    call ncf_print_variable(ncid,vitem)
	    if( bverbose ) then
	      call ncf_natts(vitem,natts)
	      if( natts > 0 ) then
	        call ncf_print_attribute_header(ncid,varid)
	        call ncf_print_attributes(ncid,vitem)
	      end if
	    end if
	  end do
	end if

!---------------------------------------------------------------------
! print info on single variable and attribute
!---------------------------------------------------------------------

	if( variable /= ' ' ) then
	  call ncf_var_id(ncid,variable,varid)
	  if( varid > 0 ) then
	    if( attribute /= ' ' ) then
	      call ncf_att_id(ncid,varid,attribute,attid)
	      if( bwrite ) then
	        write(6,*) 'varid = ',varid,'  attid = ',attid
	      end if
	      if( attid > 0 ) then
	        if( bwrite ) then
	          call ncf_print_variable_header(ncid,nitem)
	          call ncf_var_inf(ncid,varid,vitem)
	          call ncf_print_variable(ncid,vitem)
	          call ncf_att_inf(ncid,varid,attid,aitem)
	          call ncf_print_attribute_header(ncid,varid)
	          call ncf_print_attribute(aitem)
	        end if
	        call exit(0)
	      else
		write(6,*) '*** cannot find attribute: ',trim(attribute)
	        call exit(3)
	      end if
	    else
	      if( bwrite ) then
	        write(6,*) 'varid = ',varid
	        call ncf_print_variable_header(ncid,nitem)
	        call ncf_var_inf(ncid,varid,vitem)
	        call ncf_print_variable(ncid,vitem)
	        call ncf_print_attribute_header(ncid,varid)
	        call ncf_print_attributes(ncid,vitem)
	      end if
	      call exit(0)
	    end if
	  else
	    write(6,*) '*** cannot find variable: ',trim(variable)
	    call exit(1)
	  end if
	else if( attribute /= ' ' ) then
	  write(6,*) '*** for attribute also need variable'
	  call exit(5)
	end if

	if( btime ) then
	  call print_time(ncid)
	end if

!---------------------------------------------------------------------
! close netcdf file
!---------------------------------------------------------------------

	call ncf_close(ncid)

!---------------------------------------------------------------------
! end of routine
!---------------------------------------------------------------------

	end

!*********************************************************************

	subroutine print_time(ncid)

! print time records - must still be simplified
!
! time coordinate is hard coded

	use ncf

	implicit none

	integer ncid

	integer i,ierr
	integer varid,attid
	character*80 varname,attname,string,aline
	type(var_item) :: vitem
	double precision atime,atime0,dtime

	varname = 'time'
	attname = 'units'
	call ncf_var_id(ncid,varname,varid)
	call ncf_var_inf(ncid,varid,vitem)
	call ncf_var_get(ncid,vitem)
	!call ncf_print_variable(ncid,vitem)
	call ncf_att_id(ncid,varid,attname,attid)
	call ncf_att_name_string(ncid,varid,attname,string)
	!write(6,*) trim(string)
	i = index(string,'since')
	string = string(i+6:)
	i = index(string,'UTC')
	string = string(:i-2)
	!write(6,*) trim(string)
	i = index(string,' ')
	string = string(:i-1) // '::' // string(i+1:)
	!write(6,*) trim(string)
	call dts_string2time(string,atime0,ierr)
	if( ierr /= 0 ) stop 'error stop print_time: converting time'

	do i=1,vitem%len
	  dtime = vitem%value(i)
	  atime = atime0 + dtime
	  call dts_format_abs_time(atime,aline)
	  write(6,*) vitem%value(i),atime,'  ',trim(aline)
	end do

	end

!*********************************************************************

	subroutine ncinfo_init(ncfile)

! initializes command line options

	use clo
	use mod_ncinfo

	implicit none

	character*(*) ncfile

        call clo_init('ncinfo','nc-file','1.0')

        call clo_add_info('information on nc-file')

        call clo_add_sep('general options')
        call clo_add_option('verbose',.false.,'be more verbose')
        call clo_add_option('quiet',.false.,'be quiet in execution')
        call clo_add_option('silent',.false.,'do not write anything')

        call clo_add_sep('what to do')
        call clo_add_option('info',.false.,'give info on nc-file')
        call clo_add_option('var var',' ','info on variable')
        call clo_add_option('att att',' ','info on attribute')
        call clo_add_option('time',.false.,'writes time records')

        call clo_add_com('exit status 0 is success')

        call clo_parse_options

        call clo_get_option('verbose',bverbose)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('silent',bsilent)
        call clo_get_option('info',binfo)
        call clo_get_option('var',variable)
        call clo_get_option('att',attribute)
        call clo_get_option('time',btime)

        call clo_check_files(1)
        call clo_get_file(1,ncfile)

	if( bsilent ) bquiet = .true.

	end

!*********************************************************************

