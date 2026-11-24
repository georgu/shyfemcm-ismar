
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

! routine to handle pre-optimization
!
! revision log :
!
! 15.06.2025    ggu     written from scratch

!*******************************************************************

	subroutine handle_shypre_optimization(bwrite,bww,bdebug,bauto,bopti &
     &				,nkn,nel,nen3v,ng,ipv &
     &				,kphv,iarnv &
     &				,xgv,ygv,hkv,mbw,ngr)

	implicit none

	logical bwrite,bww,bdebug,bauto,bopti
	integer nkn,nel
	integer nen3v(3,nel)
	integer ng(nkn)
	integer ipv(nkn)
	integer kphv(nkn),iarnv(nkn)
	real xgv(nkn),ygv(nkn),hkv(nkn)
	integer mbw,ngr

	logical bstop
	integer ngr1
	integer nelddi
	integer, allocatable :: iphv(:)
	integer, allocatable :: iknot(:,:)
	integer, allocatable :: kvert(:,:)

	bstop = .false.
	nelddi = nel

	allocate(iphv(nkn))
	allocate(kvert(2,nkn))

	call ininum(nkn,iphv,kphv)

!--------------------------------------------------------
! setup side index and compute ngr and mbw
!--------------------------------------------------------

        if( bwrite ) write(6,*) ' ...setting up side index'

        call estimate_grade(nkn,nel,nen3v,ng,ngr1)
        allocate(iknot(ngr1,nkn))

        call sidei(nkn,nel,nen3v,ng,iknot,ngr1,ngr,bstop)
        if(bstop) goto 99934
        call check_sidei(nkn,nel,nen3v,ipv,ng,ngr1,bstop)
        if(bstop) goto 99935

        call gtest(bdebug,'bandwidth',nelddi,nkn,nel,nen3v)
        call bandw(nel,nen3v,mbw)
        if( bww ) write(6,*) 'Maximum grade of nodes is ',ngr
        if( bww ) write(6,*) 'Bandwidth is ',mbw
        if( bww ) write(6,*) 'Bandwidth is not optimized'

	return
99934   write(6,*)' (34) error: ngr1 is wrong (internal error'
        stop 'error stop handle_optimization: (34)'
99935   write(6,*)' (35) error: irregular connections'
        stop 'error stop handle_optimization: (35)'
	end

!**********************************************************

        subroutine handle_shypre_partition(bwrite,bnepart,grdpart,nn,ne &
     &				,knrank,ierank)

        use grd
        use basin

        implicit none

	logical bwrite
	logical bnepart
	character*80 grdpart
        integer nn,ne
        integer knrank(nn)
        integer ierank(ne)

        logical bgrd
        integer nnpart,nepart
        integer, allocatable :: area_node(:)
        integer, allocatable :: area_elem(:)

	return

        end

!**********************************************************

