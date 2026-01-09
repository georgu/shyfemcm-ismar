
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

! routines to cover basin with regular grid
!
! revision log :
!
! 09.01.2026    ggu     written from scratch

!****************************************************************

        subroutine handle_cover_basin(file,nlvddi,np,lmax,il,regpar,data)

	use basin
	use elabutil

	implicit none

	character*(*) file
	integer nlvddi,np,lmax
	integer il(np)
	real regpar(7)
	real data(nlvddi,np)

	logical bw

	bw = .not. bquiet

	write(6,*) 'handling cover option: ',trim(file)

        call read_cover_basin(file,bw)

	if( ncover < 0 ) then
	  call compute_ncover(ncover,nlvddi,np,lmax,il,regpar,data)
	end if

	end

!****************************************************************

        subroutine read_cover_basin(file,bw)

	use basin

	implicit none

	character*(*) file
	logical bw

	logical filex, is_grd_file

        call grd_set_write(bw)

        if( .not. filex(file) ) then
          write(6,*) 'file not existing: ',trim(file)
          stop 'error stop read_cover_basin: no such file'
        else if( basin_is_basin(file) ) then
          if( bw ) write(6,*) 'reading BAS file: ',trim(file)
          call basin_read(file)
          !breadbas = .true.
        else if( is_grd_file(file) ) then
          if( bw ) write(6,*) 'reading GRD file: ',trim(file)
          call grd_read(file)
          call grd_to_basin
          call bas_check_spherical
          call estimate_ngr(ngr)
          !breadbas = .false.
        else
          write(6,*) 'Cannot read this file: ',trim(file)
          stop 'error stop read_cover_basin: unknown format'
        end if

	end

!****************************************************************

	subroutine compute_ncover(ncover,nlvddi,np,lmax,il,regpar,data)

	use basin

	implicit none

	integer ncover
	integer nlvddi,np,lmax
	integer il(np)
	real regpar(7)
	real data(nlvddi,np)

	integer nx,ny
	integer regexpand,nloop
	integer ierr
	real x0,y0,dx,dy,flag,x1,y1
	real regval(np)
	real femval(nkn)

	nloop = 0
	regexpand = 1
        nx = nint(regpar(1))
        ny = nint(regpar(2))
        x0 = regpar(3)
        y0 = regpar(4)
        dx = regpar(5)
        dy = regpar(6)
        flag = regpar(7)

	x1 = x0 + (nx-1)*dx
	y1 = y0 + (nx-1)*dy

        if( nx*ny /= np ) then
          write(6,*) 'np,nx,ny,nx*ny: ',np,nx,ny,nx*ny
          stop 'error stop compute_ncover: incompatible params'
        end if

	if( nlvddi /= lmax ) goto 98
	if( x0 > minval(xgv) ) goto 99
	if( y0 > minval(ygv) ) goto 99
	if( x1 < maxval(xgv) ) goto 99
	if( y1 < maxval(ygv) ) goto 99

	regval = data(1,:)
        call intp_reg( nx, ny, x0, y0, dx, dy, flag &
     &                  , regval, nkn, xgv, ygv, femval, ierr )
	write(6,*) 'flag values found: ',0,ierr
	if( ierr == 0 ) return

	do
	  nloop = nloop + 1
	  call reg_set_flag(nlvddi,np,il,regpar,data)
	  call reg_expand_3d(nlvddi,nx,ny,lmax,regexpand,flag,data)

	  regval = data(1,:)
          call intp_reg( nx, ny, x0, y0, dx, dy, flag &
     &                  , regval, nkn, xgv, ygv, femval, ierr )
	  write(6,*) 'flag values found: ',nloop,ierr
	  if( ierr == 0 ) exit
	end do

	return
   98	continue
	write(6,*) 'nlvddi,lmax: ',nlvddi,lmax
	stop 'error stop compute_ncover: error in vertical dimension'
   99	continue
	write(6,*) 'x0, y0, x1, y1 of regular domain'
	write(6,*) x0, y0, x1, y1
	write(6,*) 'x/y min/max of fem domain'
	write(6,*) minval(xgv),minval(ygv),maxval(xgv),maxval(ygv)
	stop 'error stop compute_ncover: regular domain too small'
	end

!****************************************************************

