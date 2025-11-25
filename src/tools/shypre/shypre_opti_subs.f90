
!--------------------------------------------------------------------------
!
!    Copyright (C) 2024  Georg Umgiesser
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

! pre-processing routine (intermediate version) - can be deleted
!
! revision log :
!
! 08.11.2024    ggu     last intermediate version
 
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
        call knscr(nkn,ngr,ngr1,iknot)

        call gtest(bdebug,'bandwidth',nelddi,nkn,nel,nen3v)
        call bandw(nel,nen3v,mbw)
        if( bww ) write(6,*) 'Maximum grade of nodes is ',ngr
        if( bww ) write(6,*) 'Bandwidth is ',mbw

!--------------------------------------------------------
! bandwidth optimization
!--------------------------------------------------------

        if( .not. bopti ) then
	  if( bww ) write(6,*) 'Bandwidth is not optimized'
	  return
	end if

        call gtest(bdebug,'bandwidth 1',nelddi,nkn,nel,nen3v)
        if( bwrite ) write(6,*) ' ...optimizing of band width'
        call bandop(nkn,ngr,ipv,iphv,kphv,ng,iknot,kvert &
     &                  ,bopti,bauto,bww)
        call gtest(bdebug,'bandwidth 2',nelddi,nkn,nel,nen3v)

        if(bopti) then
          call bandex(nkn,nel,nen3v,kphv,ipv,iarnv,xgv,ygv,hkv) !kphv is n-rank
          call gtest(bdebug,'bandwidth 3',nelddi,nkn,nel,nen3v)
          call bandw(nel,nen3v,mbw)
          call gtest(bdebug,'bandwidth 4',nelddi,nkn,nel,nen3v)
          if( bww ) write(6,*) 'Optimized bandwidth is ',mbw
        end if

	return
99934   write(6,*)' (34) error: ngr1 is wrong (internal error'
        stop 'error stop handle_optimization: (34)'
99935   write(6,*)' (35) error: irregular connections'
        stop 'error stop handle_optimization: (35)'
	end

!**********************************************************

        subroutine bandop(nkn,ngr1,ipv,iphv,kphv,ng,iknot,kvert &
     &                          ,bopti,bauto,bwrite)

! optimize band width

        implicit none

        integer nkn,ngr1
        logical bopti,bauto
        integer ipv(nkn)
        integer iphv(nkn),kphv(nkn)
        integer ng(nkn)
        integer iknot(ngr1,nkn)
        integer kvert(2,nkn)
        logical bwrite,bloop

        integer iantw

        if( bauto ) then
          call ininum(nkn,iphv,kphv)
          call optest('before optimization: ',nkn,ipv,iphv,kphv)
          call cmgrade(nkn,ngr1,ipv,iphv,kphv,ng,iknot,1,4,bwrite)
          call optest('after Cuthill McKee: ',nkn,ipv,iphv,kphv)
          call revnum(nkn,iphv,kphv,bwrite)
          call optest('after reversing nodes: ',nkn,ipv,iphv,kphv)
          call rosen(nkn,ngr1,iphv,kphv,ng,iknot,kvert,bwrite)
          call optest('after Rosen: ',nkn,ipv,iphv,kphv)
          return
        end if

        bloop = .true.

        do while( bloop )
          call ininum(nkn,iphv,kphv)

!         call anneal(nkn,ngr1,kphv,ng,iknot,iphv,kvert)

          if( iantw(' Cuthill McKee algorithm ?') .gt. 0 ) then
            call cmv(nkn,ngr1,ipv,iphv,kphv,ng,iknot,bwrite)
          end if

          if( iantw(' Reverse numbering of nodes ?') .gt. 0 ) then
            call revnum(nkn,iphv,kphv,bwrite)
          end if

          if( iantw(' Rosen algorithm ?') .gt. 0 ) then
            call rosen(nkn,ngr1,iphv,kphv,ng,iknot,kvert,bwrite)
          end if

          bloop = iantw(' Repeat optimization of bandwidth ?') .gt. 0
        end do

        bopti = .true.

        end

!**********************************************************

        subroutine revnum(nkn,iphv,kphv,bwrite)

! reverses numbering of nodes

        implicit none

        integer nkn
        integer iphv(nkn), kphv(nkn)
        logical bwrite

        integer i

        if( bwrite ) write(6,*) 'Applying reverse algorithm...'

        do i=1,nkn
          kphv(i) = nkn+1-kphv(i)
          iphv(kphv(i)) = i
        end do

        end

!**********************************************************


        subroutine zernum(nkn,iphv,kphv)

! zeros numbering of nodes

        implicit none

        integer nkn
        integer iphv(nkn), kphv(nkn)

        integer i

        do i=1,nkn
          iphv(i) = 0
          kphv(i) = 0
        end do

        end

!**********************************************************

        subroutine bandex(nkn,nel,nen3v,kphv &
     &                      ,ipv,iarnv,xgv,ygv,hkv)

! exchange nodes after optimization

        implicit none

        integer nkn,nel
        integer nen3v(3,nel)
        integer kphv(nkn)
        integer ipv(nkn)
        integer iarnv(nkn)
        real xgv(nkn),ygv(nkn)
        real hkv(nkn)

        integer ie,ii
        integer neaux(3,nel)

        do ie=1,nel
          do ii=1,3
            neaux(ii,ie)=nen3v(ii,ie)
          end do
        end do

        do ie=1,nel
          do ii=1,3
            nen3v(ii,ie)=kphv(neaux(ii,ie))
          end do
        end do

!       copy arrays with kphv as rank table

        call icopy(nkn,ipv,kphv)
        call icopy(nkn,iarnv,kphv)
        call rcopy(nkn,xgv,kphv)
        call rcopy(nkn,ygv,kphv)
        call rcopy(nkn,hkv,kphv)

        return
        end

!**********************************************************

        subroutine optest(text,nkn,ipv,iphv,kphv)

        implicit none

        character*(*) text
        integer nkn
        integer ipv(nkn)
        integer iphv(nkn)
        integer kphv(nkn)

        integer i

        do i=1,nkn
          if( iphv(i) .gt. nkn ) goto 99
          if( kphv(i) .gt. nkn ) goto 99
          if( iphv(kphv(i)) .ne. i ) goto 99
          if( kphv(iphv(i)) .ne. i ) goto 99
        end do

        return
   99   continue
        write(6,*) '*** Error in optest: '
        write(6,*) text
        write(6,*) 'error in pointers...'
        write(6,*) 'problem is close to following node...'
        write(6,*) 'node (intern/extern): ',i,ipv(i)
        write(6,*) i,nkn
        write(6,*) iphv(i),kphv(i)
        write(6,*) iphv(kphv(i)),kphv(iphv(i))
        write(6,*) 'maybe the domain is not connected...'
        stop 'error stop optest: pointers'
        end

!**********************************************************

        subroutine knscr(nkn,ngr,ngr1,iknot)

! de-scrambles side index for cmv and rosen

        implicit none

        integer nkn,ngr,ngr1
        integer iknot(ngr1*nkn)

        integer k,j

        do k=1,nkn
          do j=1,ngr
            iknot((k-1)*ngr+j)=iknot((k-1)*ngr1+j)
          end do
        end do

        return
        end

!**********************************************************
!**********************************************************
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

        integer i

        if( bwrite ) write(6,*) 'handle_partition: ',nn,ne

        bgrd = ( grdpart /= ' ' )
        if( .not. bnepart .and. .not. bgrd ) return

        if( bnepart .and. bgrd ) then
          write(6,*) 'only one of -partition and -nepart can be given'
          stop 'error stop handle_partition: options'
        end if

        allocate(area_node(nn))
        allocate(area_elem(ne))

        if( bgrd ) then
          write(6,*) 'reading partitioning file ',trim(grdpart)
          call grd_read(grdpart)
          if( nk_grd /= nn ) goto 99
          if( ne_grd /= ne ) goto 99
          area_node = ianv
          area_elem = iaev
        else
          write(6,*) 'using area code as partitioning information'
          area_node = iarnv
          area_elem = iarv
        end if

        call renumber_partition(nn,area_node,nnpart)
        call renumber_partition(ne,area_elem,nepart)

        if( bgrd ) then
          call icopy(nn,area_node,knrank)
          call icopy(ne,area_elem,ierank)
        end if

        call basin_set_partition(nn,ne,nnpart,nepart,area_node,area_elem)

        return
   99   continue
        write(6,*) nk_grd,nn
        write(6,*) ne_grd,ne
        stop 'error stop handle_partition: incompatibility'
        end

!**********************************************************

        subroutine renumber_partition(n,area,npart)

        implicit none

        integer n
        integer area(n)
        integer npart

        integer i,ia,imax
        integer nmin,nmax
        integer, allocatable :: table_in(:),table_out(:)

        nmin = minval(area)
        nmax = maxval(area)

        if( nmin == nmax ) then         !no partition
          npart = 0
          area = 0
          return
        end if

        if( nmin < 0 .or. nmax < 0 ) then
          write(6,*) nmin,nmax
          stop 'error stop renumber_partition: nmin,nmax'
        end if

        allocate(table_in(0:nmax),table_out(0:nmax))

        table_in = 0
        table_out = 0

        do i=1,n
          ia = area(i)
          table_in(ia) = table_in(ia) + 1
        end do

        imax = -1
        do ia=0,nmax
          if( table_in(ia) > 0 ) then
            imax = imax + 1
            table_out(ia) = imax
          end if
        end do
        npart = imax

        do i=1,n
          ia = area(i)
          area(i) = table_out(ia)
        end do

        end

!**********************************************************

