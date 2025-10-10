
!--------------------------------------------------------------------------
!
!    Copyright (C) 2025  Georg Umgiesser
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
! 09.10.2025    ggu     created from scratch for smoothing regular field
!
!***************************************************************

	subroutine smooth_regular(nlvdi,nx,ny,cv3,flag,salpha,sloop)

	implicit none

	integer nlvdi,nx,ny
	real cv3(nlvdi,nx,ny)
	real flag		!flag for no data
	real salpha		!smoothing parameter [0-1]
	integer sloop		!how many smoothing loops

	integer iloop,l
	integer ix,iy,ixn,iyn
	integer, allocatable :: icount(:,:)
	double precision cmed
	double precision, allocatable :: cold(:,:),cnew(:,:)
 
	allocate(cnew(nx,ny))
	allocate(cold(nx,ny))
	allocate(icount(nx,ny))

	do iloop=1,sloop

	do l=1,nlvdi
	  cold = cv3(l,:,:)
	  cnew = 0.
	  icount = 0
	  do ix=1,nx
	    do iy=1,ny
	      if( cold(ix,iy) == flag ) cycle
	      ixn = ix - 1
	      if( ix > 1 ) then
		cnew(ixn,iy) = cnew(ixn,iy) + cold(ix,iy)
		icount(ixn,iy) = icount(ixn,iy) + 1
	      end if
	      ixn = ix + 1
	      if( ix < nx ) then
		cnew(ixn,iy) = cnew(ixn,iy) + cold(ix,iy)
		icount(ixn,iy) = icount(ixn,iy) + 1
	      end if
	      iyn = iy - 1
	      if( iy > 1 ) then
		cnew(ix,iyn) = cnew(ix,iyn) + cold(ix,iy)
		icount(ix,iyn) = icount(ix,iyn) + 1
	      end if
	      iyn = iy + 1
	      if( iy < ny ) then
		cnew(ix,iyn) = cnew(ix,iyn) + cold(ix,iy)
		icount(ix,iyn) = icount(ix,iyn) + 1
	      end if
	    end do
	  end do

	  do ix=1,nx
	    do iy=1,ny
	      if( icount(ix,iy) == 0 ) cycle
	      cmed = cnew(ix,iy) / icount(ix,iy)
              cnew(ix,iy) = (1.-salpha)*cold(ix,iy) + salpha*cmed
	    end do
          end do

	  cv3(l,:,:) = real(cnew)

	end do

	end do

	end

!***************************************************************

