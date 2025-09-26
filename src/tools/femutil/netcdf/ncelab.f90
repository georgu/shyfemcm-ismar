
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

!===================================================================
	module mod_ncelab
!===================================================================

	logical, save :: babout = .false.
	logical, save :: binfo = .false.
	logical, save :: bverbose = .false.
	logical, save :: bquiet = .false.
	logical, save :: bsilent = .false.
	logical, save :: bdebug = .false.

	logical, save :: btime = .false.
	logical, save :: bcoords = .false.

!===================================================================
	end module mod_ncelab
!===================================================================

        program ncelab

	use ncf
	use mod_ncelab

        implicit none

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

	call ncelab_init(ncfile)

	call ncf_open_read(ncfile,ncid)

	bwrite = .not. bquiet

!---------------------------------------------------------------------
! print info on file, dimensions and global attributes
!---------------------------------------------------------------------

	if( binfo ) call ncelab_info(ncid)

        call get_dims_and_coords(ncid,bwrite                    &
     &                  ,nt,nx,ny,nz                                    &
     &                  ,tcoord,xcoord,ycoord,zcoord)

        if( btime ) then
          if( bwrite ) call ncf_print_all_time_records(ncid)
          call ncf_print_minmax_time_records(ncid)
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

	subroutine ncelab_info(ncid)

! writes general info on nc-file

	use ncf
	use mod_ncelab

	implicit none
	
	integer ncid

	integer ngatts,idunlim,id
	integer natts,nvars,varid
	type(nc_item) :: nitem
	type(att_item) :: aitem
	type(var_item) :: vitem

	nitem = ncf_get_nitem(ncid)

	ngatts = nitem%ngatts
	idunlim = nitem%idunlim

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

!	call ncnames_init
!	call ncnames_get_dims_and_coords(ncid,bverbose
!     +				,nt,nx,ny,nz
!     +				,tcoord,xcoord,ycoord,zcoord)

        call ncf_set_var_name_length(ncid,nitem)
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

	end

!*********************************************************************

        subroutine ncelab_about

        implicit none

        write(6,*) 'converts nc (netcdf) file to fem file'
        write(6,*)
        write(6,*) 'The file created is either a regular fem file'
        write(6,*) 'or it can be single points to be used for'
        write(6,*) 'boundary conditions given with the option -single'
        write(6,*) 'The domain can be adjusted with -domain'
        write(6,*) 'The variables to be written are given with -vars'
        write(6,*)
        write(6,*) 'The program should recognize most of the'
        write(6,*) 'dimensions and coordinates used in nc files'
        write(6,*) 'If some of these are not recognized you can'
        write(6,*) 'insert them at the end of file nc_dim_coords.f'
        write(6,*) 'and then recompile with "make fem"'
        write(6,*) 'The same is true for the description of the'
        write(6,*) 'variables written to file'

        end

!*********************************************************************
!*********************************************************************
!*********************************************************************

	subroutine ncelab_init(ncfile)

! initializes command line options

	use clo
	use mod_ncelab

	implicit none

	character*(*) ncfile

        call clo_init('ncelab','nc-file','1.0')

        call clo_add_info('elaborates nc-file and creates fem file')

        call clo_add_sep('general options')
        call clo_add_option('about',.false.,'about this program')
        call clo_add_option('info',.false.,'general info on nc-file')
        call clo_add_option('verbose',.false.,'be more verbose')
        call clo_add_option('quiet',.false.,'be quiet in execution')
        call clo_add_option('silent',.false.,'do not write anything')
        call clo_add_option('debug',.false.,'writes debug information')

        call clo_add_sep('what to do')
        call clo_add_option('time',.false. &
     &			,'write available time records to terminal')
        call clo_add_option('coords',.false.,'write coordinate file')

        call clo_add_com('exit status 0 is success')

        call clo_parse_options

        call clo_get_option('about',babout)
        call clo_get_option('info',binfo)
        call clo_get_option('verbose',bverbose)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('silent',bsilent)
        call clo_get_option('debug',bdebug)

        call clo_get_option('time',btime)
        call clo_get_option('coords',bcoords)

	if( babout ) call ncelab_about

        call clo_check_files(1)
        call clo_get_file(1,ncfile)

	if( bsilent ) bquiet = .true.

	end

!*********************************************************************

