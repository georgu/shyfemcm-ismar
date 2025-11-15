
!--------------------------------------------------------------------------
!
!    Copyright (C) 2004,2006,2008-2012,2014-2015,2014-2015  Georg Umgiesser
!    Copyright (C) 2017-2020  Georg Umgiesser
!    Copyright (C) 2012  Andrea Cucco
!    Copyright (C) 2014  Christian Ferrarin
!    Copyright (C) 2018  Marco Bajo
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

! modules for restart routines
!
! revision log :
!
! 08.03.2025    ggu     extract module for restart routines
! 19.06.2025    ggu     updated
!
! notes :
!
!  |--------|---------|------------|--------------|--------------|
!  | ityrst | no file | file empty ! itrec/=itrst ! itrec==itrst !
!  |--------|---------|------------|--------------|--------------|
!  |   0    |  error  |   error    |     error    |    warm      |
!  |--------|---------|------------|--------------|--------------|
!  |   1    |  cold   |   error    |     error    |    warm      |
!  |--------|---------|------------|--------------|--------------|
!  |   2    |  cold   |   cold     |     warm     |    warm      |
!  |--------|---------|------------|--------------|--------------|
!
! versions :
!
! <3	not supported
! 3	write hydro
! 4	write hm3v
! 5	write saltv,tempv,rhov
! 6	write conv (tracer)
! 7	write wlnv (vertical velocity)
! 8	write ecological variables
! 9	write date,time
! 10	write hlv
! 11	write atime, write idfile (regular header)
! 12	mercury restart
! 13	write ilhv and ilhkv
! 14	write vertical only for nlv > 1
! 15	write gotm arrays
! 16	adapted for mpi
! 17	write bfm restart
!
!*********************************************************************

!=====================================================================
	module mod_restart
!=====================================================================

	implicit none

	logical, save :: bok_rst = .false.	!restart file has been read
	integer, save :: nvers_rst = 0
	integer, save :: iflag_want_rst  = -1
	integer, save :: iflag_avail_rst = -1

	integer, save :: idfrst = 749652	!id for restart file

	integer, save :: nvmax = 17		!last version of file
	integer, parameter :: nidmax = 9

	integer, save :: id_hydro_rst = 1	!1		hydro
	integer, save :: id_depth_rst = 2	!10		depth
	integer, save :: id_barcl_rst = 3	!100		t/s/rho
	integer, save :: id_conz_rst  = 4	!1000		tracer
	integer, save :: id_wvert_rst = 5	!10000		vertical vel.
	integer, save :: id_eco_rst   = 6	!100000		ecology
	integer, save :: id_merc_rst  = 7	!1000000	mercury
	integer, save :: id_gotm_rst  = 8	!10000000	gotm
	integer, save :: id_bfm_rst   = 9	!100000000	bfm

	integer, save :: ibarcl_rst = 0
	integer, save :: iconz_rst  = 0
	integer, save :: iwvert_rst = 0
	integer, save :: ieco_rst   = 0
	integer, save :: imerc_rst  = 0
	integer, save :: iturb_rst  = 0

	character*20, save :: descript_rst(nidmax) = (/    &
     &		 'hydrodynamics       '    &
     &		,'depth               '    &
     &		,'T/S/rho             '    &
     &		,'tracer concentration'    &
     &		,'vertical velocities '    &
     &		,'ecological model    '    &
     &		,'mercury model       '    &
     &		,'gotm turb model     '    &
     &		,'bfm model           '    &
     &						/)

	real, save, allocatable :: hlvrst(:)
	integer, save, allocatable :: ilhrst(:)
	integer, save, allocatable :: ilhkrst(:)

	logical, save :: bmaster = .false.

        INTERFACE restart_write_value
        MODULE PROCEDURE   restart_write_value_scalar_i  &
     &                    ,restart_write_value_scalars_i &
     &                    ,restart_write_value_2d_i      &
     &                    ,restart_write_value_2d_r      &
     &                    ,restart_write_value_3d_r      &
     &                    ,restart_write_value_3d_d      &
     &                    ,restart_write_value_fix_r
        END INTERFACE

        INTERFACE restart_read_value
        MODULE PROCEDURE   restart_read_value_2d_i    &
     &                    ,restart_read_value_2d_r    &
     &                    ,restart_read_value_3d_r    &
     &                    ,restart_read_value_3d_d    &
     &                    ,restart_read_value_fix_r
        END INTERFACE

!=====================================================================
	contains
!=====================================================================

	function get_nv_global(nv_local)

	use shympi

	integer get_nv_global
	integer nv_local

	integer nv_global

	nv_global = nlv_global
        if( nv_local /= nlv_local ) then
          if( nv_local == nlv_local+1 ) then
            nv_global = nv_global + 1
          else
	    write(6,*) nv_local,nlv_local
            stop 'error stop get_nv_global: nv incompatible'
          end if
        end if

	get_nv_global = nv_global

	end function
	
!--------------------------------

	function get_nn_global(belem,nn)

	use shympi

	integer get_nn_global
	logical belem
	integer nn

	integer nn_global

        if( belem ) then
	  nn_global = nel_global
          if( nn /= nel_local ) goto 99
        else
	  nn_global = nkn_global
          if( nn /= nkn_local ) goto 99
        end if

	get_nn_global = nn_global

        return
   99   continue
        write(6,*) belem,nn,nel_local,nkn_local
        stop 'error stop get_nn_global: internal error'
	end function
	
!--------------------------------

	subroutine restart_write_value_scalar_i(iu,ival)

	integer iu
	integer ival

	if( bmaster ) write(iu) ival

	end subroutine

!--------------------------------

	subroutine restart_write_value_scalars_i(iu,n,ivals)

	integer iu
	integer n
	integer ivals(n)

	if( bmaster ) write(iu) ivals

	end subroutine

!--------------------------------

	subroutine restart_write_value_2d_i(iu,belem,array)

	use shympi

	integer iu
	logical belem
	integer array(:)

	integer nn_global,nn
	integer, allocatable :: garray(:)

	if( bmpi ) then
	  nn = size(array)
	  nn_global = get_nn_global(belem,nn)
	  allocate(garray(nn_global))
	  call shympi_l2g_array(array,garray)
	  if( bmaster ) write(iu) garray
	else
	  write(iu) array
	end if

	end subroutine

!--------------------------------

	subroutine restart_write_value_2d_r(iu,belem,array)

	use shympi

	integer iu
	logical belem
	real array(:)

	integer nn_global,nn
	real, allocatable :: garray(:)

	if( bmpi ) then
	  nn = size(array)
	  nn_global = get_nn_global(belem,nn)
	  allocate(garray(nn_global))
	  call shympi_l2g_array(array,garray)
	  if( bmaster ) write(iu) garray
	else
	  write(iu) array
	end if

	end subroutine

!--------------------------------

	subroutine restart_write_value_3d_r(iu,belem,array)

	use shympi

	integer iu
	logical belem
	real array(:,:)

	integer nn_global,nv_global,nn
	real, allocatable :: garray(:,:)

	if( bmpi ) then
	  nn = size(array,2)
	  nv_global = get_nv_global(size(array,1))
	  nn_global = get_nn_global(belem,nn)
	  allocate(garray(nv_global,nn_global))
	  call shympi_l2g_array(array,garray)
	  if( bmaster ) write(iu) garray
	else
	  write(iu) array
	end if

	end subroutine

!--------------------------------

	subroutine restart_write_value_3d_d(iu,belem,array)

	use shympi

	integer iu
	logical belem
	double precision array(:,:)

	integer nn_global,nv_global,nn
	double precision, allocatable :: garray(:,:)

	if( bmpi ) then
	  nn = size(array,2)
	  nv_global = get_nv_global(size(array,1))
	  nn_global = get_nn_global(belem,nn)
	  allocate(garray(nv_global,nn_global))
	  call shympi_l2g_array(array,garray)
	  if( bmaster ) write(iu) garray
	else
	  write(iu) array
	end if

	end subroutine

!--------------------------------

	subroutine restart_write_value_fix_r(iu,belem,nfix,array)

	use shympi

	integer iu
	logical belem
	integer nfix
	real array(:,:)

	integer nn_global,nv_global,nn
	real, allocatable :: garray(:,:)

	if( bmpi ) then
	  if( size(array,1) /= nfix ) then
	    stop 'error stop restart_write_value_fix: nv incompatible'
	  end if
	  nn = size(array,2)
	  nn_global = get_nn_global(belem,nn)
	  allocate(garray(nfix,nn_global))
	  call shympi_l2g_array(nfix,array,garray)
	  if( bmaster ) write(iu) garray
	else
	  write(iu) array
	end if

	end subroutine

!--------------------------------

	subroutine restart_read_value_2d_i(iu,belem,array)

	use shympi

	integer iu
	logical belem
	integer array(:)

	integer nn_global,nn
	integer, allocatable :: garray(:)

	if( bmpi ) then
	  nn = size(array)
	  nn_global = get_nn_global(belem,nn)
	  allocate(garray(nn_global))
	  read(iu) garray
	  call shympi_g2l_array(garray,array)
	else
	  read(iu) array
	end if

	end subroutine

!--------------------------------

	subroutine restart_read_value_2d_r(iu,belem,array)

	use shympi

	integer iu
	logical belem
	real array(:)

	integer nn_global,nn
	real, allocatable :: garray(:)

	if( bmpi ) then
	  nn = size(array)
	  nn_global = get_nn_global(belem,nn)
	  allocate(garray(nn_global))
	  read(iu) garray
	  call shympi_g2l_array(garray,array)
	else
	  read(iu) array
	end if

	end subroutine

!--------------------------------

	subroutine restart_read_value_3d_r(iu,belem,array)

	use shympi

	integer iu
	logical belem
	real array(:,:)

	integer nn_global,nv_global,nn
	real, allocatable :: garray(:,:)

	if( bmpi ) then
	  nn = size(array,2)
	  nv_global = get_nv_global(size(array,1))
	  nn_global = get_nn_global(belem,nn)
	  allocate(garray(nv_global,nn_global))
	  read(iu) garray
	  call shympi_g2l_array(garray,array)
	else
	  read(iu) array
	end if

	end subroutine

!--------------------------------

	subroutine restart_read_value_3d_d(iu,belem,array)

	use shympi

	integer iu
	logical belem
	double precision array(:,:)

	integer nn_global,nv_global,nn
	double precision, allocatable :: garray(:,:)

	if( bmpi ) then
	  nn = size(array,2)
	  nv_global = get_nv_global(size(array,1))
	  nn_global = get_nn_global(belem,nn)
	  allocate(garray(nv_global,nn_global))
	  read(iu) garray
	  call shympi_g2l_array(garray,array)
	else
	  read(iu) array
	end if

	end subroutine

!--------------------------------

	subroutine restart_read_value_fix_r(iu,belem,nfix,array)

	use shympi

	integer iu
	logical belem
	integer nfix
	real array(:,:)

	integer nn_global,nv_global,nn
	real, allocatable :: garray(:,:)

	if( bmpi ) then
	  if( size(array,1) /= nfix ) then
	    stop 'error stop restart_read_value_fix: nv incompatible'
	  end if
	  nn = size(array,2)
	  nn_global = get_nn_global(belem,nn)
	  allocate(garray(nfix,nn_global))
	  read(iu) garray
	  call shympi_g2l_array(garray,array)
	else
	  read(iu) array
	end if

	end subroutine

!=====================================================================
	end module mod_restart
!=====================================================================

