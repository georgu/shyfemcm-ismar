
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

! revision log :
!
! 05.10.2026    ggu     written from scratch

!*********************************************************************

	subroutine make_attributes_time(ncid,varid,units)

	use ncf

	implicit none

	integer ncid,varid
	character*(*) units

	call ncf_make_att(ncid,varid,'standard_name','time')
	call ncf_make_att(ncid,varid,'units',units)
	call ncf_make_att(ncid,varid,'calendar','standard')
	call ncf_make_att(ncid,varid,'axis','T')

	end

!*********************************************************************

	subroutine make_attributes_latitude(ncid,varid)

	use ncf

	implicit none

	integer ncid,varid

	call ncf_make_att(ncid,varid,'standard_name','latitude')
	call ncf_make_att(ncid,varid,'long_name','Latitude of nodes')
	call ncf_make_att(ncid,varid,'units','degrees_north')

	end

!*********************************************************************

	subroutine make_attributes_longitude(ncid,varid)

	use ncf

	implicit none

	integer ncid,varid

	call ncf_make_att(ncid,varid,'standard_name','longitude')
	call ncf_make_att(ncid,varid,'long_name','Longitude of nodes')
	call ncf_make_att(ncid,varid,'units','degrees_east')

	end

!*********************************************************************

	subroutine make_attributes_waterlevel(ncid,varid)

	use ncf

	implicit none

	integer ncid,varid

	call ncf_make_att(ncid,varid,'standard_name','water_surface_height_above_reference_datum')
	call ncf_make_att(ncid,varid,'long_name','Sea Surface Elevation')
	call ncf_make_att(ncid,varid,'units','m')
	call ncf_make_att(ncid,varid,'valid_min',-10.)
	call ncf_make_att(ncid,varid,'valid_max',10.)
	call ncf_make_att(ncid,varid,'_FillValue',-999.)
	call ncf_make_att(ncid,varid,'missing_value',-999.)

	end

!*********************************************************************

