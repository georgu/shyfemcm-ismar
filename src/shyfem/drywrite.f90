
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

! routines for writing dry information
!
! revision log :
!
! 05.11.2026    ggu     (from scratch)

!*****************************************************************

	subroutine shy_write_dry_info(id,dtime)

	use basin
	use mod_hydro
	use mod_geom_dynamic

	implicit none

	integer id
	double precision dtime
	integer iwngv(nkn)
	real rwngv(nkn)

	call drywet_encrypt(iwegv,iwngv)
	rwngv = iwngv

	call shy_write_scalar_record(id,dtime,-1,1,znv)
	call shy_write_scalar_record(id,dtime,-74,1,rwngv)

	end

!*****************************************************************

