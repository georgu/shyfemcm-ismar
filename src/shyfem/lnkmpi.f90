
!--------------------------------------------------------------------------
!
!    Copyright (C) 2003,2008-2015,2019  Georg Umgiesser
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

! topological set up routines - correct for mpi domains
!
! revision log :
!
! 13.04.2022	ggu	newly written
! 26.04.2022	ggu	write grd file if needed
! 08.06.2022	ggu	new routine write_grd_domain()
! 15.10.2022    ggu     shympi_exchange_array substituted with shympi_l2g_array
! 12.12.2022    ggu     new routine write_grd_general()
! 18.03.2023    ggu     adapted to new id_elem(0:3,:)
! 20.03.2023    ggu     write_grd_domain() is now called in shympi_setup()
! 29.03.2023    ggu     refactoring
! 07.11.2024    ggu     updated to new connection framework
! 13.11.2024    ggu     indicate other domain (not yet finished)
! 04.02.2025    ggu     new routine check_2id_elements() (not yet finished)
! 22.04.2026    ggu&clc bug fix adjusting kantv (kk<1)
!
!*****************************************************************

	subroutine set_geom_mpi

! sets up geometrical arrays
!
! geom structures have already been setup on local domain
! with this call they are corrected with global boundary node information

	use mod_geom
	use mod_connect
	use basin
	use shympi

        implicit none

        integer i,n
	integer k,kk,kn,kb,k1,k2
	integer ie,ii
	integer ie_ext,ien
	integer id,id_neigh
	integer ide,idk,nid
	integer iunit
	integer, allocatable :: ibound(:)
	integer, allocatable :: ecgv(:,:)	! global ieltv
        integer kerr,ierr
	integer ival_in(3),ival_out(3)

	integer knext,kbhnd

	!if( .not. bmpi ) return

	iunit = 660 + my_id

!-------------------------------------------------------------
! make global ibound array
!-------------------------------------------------------------

	allocate(ibound(nkn_global))
	allocate(ecgv(3,nel_global))

        call connect_init(nkn_global,nel_global,nen3v_global,ierr)
        call connect_get_bound(nkn_global,ibound)
        call connect_get_ecv(nel_global,ecgv)

!-------------------------------------------------------------
! reduce to local
!-------------------------------------------------------------

! ibound are global boundary nodes
! iboundv are local boundary nodes

	call shympi_barrier
	call shympi_g2l_array(ibound,iboundv)

!-------------------------------------------------------------
! now we adjust indices on local domain
!-------------------------------------------------------------

!-------------------------------------------------------------
! adjust kantv
!-------------------------------------------------------------

	do k=1,nkn
	  if( iboundv(k) == 0 ) then			!not a boundary node
	    kantv(:,k) = 0
	  else if( .not. shympi_is_inner_node(k) ) then	!in other domain
	    do i=1,2
	      kk = kantv(i,k)
	      id = -1
	      if( kk > 0 ) id = id_node(kk)
	      if( id /= my_id ) kantv(i,k) = 0		!node outside domain
	    end do
	  end if
	end do

!-------------------------------------------------------------
! adjust ieltv (-1 is reserved for OB)
! > 0  indicates local neighbor of element
! == 0 indicates material boundary
! < 0  indicates open boundary
! <= 1000  indicates other domain
!-------------------------------------------------------------

	do ie=1,nel
	  ie_ext = ipev(ie)
	  nid = id_elem(0,ie)			! how many domains
	  ide = id_elem(1,ie)			! main id of element
	  if( nid == 1 ) cycle			! internal element
	  if( nid == 2 .and. ide == my_id ) cycle ! 2 domains but this is main
	  ! now look for oposite side of node in this domain
	  do ii=1,3
	    if( ieltv(ii,ie) /= 0 ) cycle	! has local neighbor
	    k = nen3v(ii,ie)
	    idk = id_node(k)
	    if( idk /= my_id ) cycle
	    kn = knext(k,ie)
	    kb = kbhnd(k,ie)
	    ! dont do this if both nodes are on boundary
	    if( iboundv(kn) > 0 .and. iboundv(kb) > 0 ) cycle
	    if( nid == 3 ) then			! tripple point
	      do i=1,3
		if( id_elem(i,ie) /= my_id ) exit	! choose any not my_id
	      end do
	      ide = id_elem(i,ie)
	    end if
	    ieltv(ii,ie) = -1000 - ide	! internal (not real) border
	  end do
	end do

!-------------------------------------------------------------
! do sanity checks
!-------------------------------------------------------------

	!call check_2id_elements(ecgv)

	return

	do ie=1,nel
	  do ii=1,3
	    ien = ieltv(ii,ie)
	    if( ien > 0 ) then
	      if( count(ieltv(:,ien)==ie) /= 1 ) goto 98
	    else if( ien <= 1000 ) then
!	      ide = -ien - 1000
!	      ie_ext = ipev(ie)
!	      iel = gfindloc(ip_ext_elem,ie_ext)
!	      if( iel <= 0 ) goto 97
!	      call shympi_receive(my_id,ide,3,ival_in,ival_out)
!	      ien = ecgv(ii,iel)	!look up in global index
!	      ien_ext = ip_ext_elem(ien)
!	      iel1 = gfindloc(ip_ext_elem,ie_ext)
	    end if
	  end do
	end do

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	return
   98	continue
	stop 'error stop set_geom_mpi: internal error (2)'
   99	continue
	stop 'error stop set_geom_mpi: internal error (1)'
	end

!*****************************************************************

	subroutine check_2id_elements(ecgv)

	use mod_geom
	use mod_connect
	use basin
	use shympi

	implicit none

	integer ecgv(3,nel_global)

	integer, parameter :: ndim = 6
	integer ie,ii,ien,ies
	integer ie_ext
	integer nid,ide,ide2
	integer ip,itot,iint,nit,n
	integer iu
	integer buffer_in(nel*ndim)
	integer buffer_out(nel*ndim)
	integer ie_tripple(3,n_threads)
	integer ie_aux(3*n_threads)
	integer tripple_info(5,n_threads)

	integer, allocatable :: ielt_aux1(:,:)
	integer, allocatable :: ielt_aux2(:,:)
	integer, allocatable :: ielt_aux3(:,:)
	integer, allocatable :: ieltv_global(:,:)
	integer, allocatable :: aux1(:)
	integer, allocatable :: aux2(:)

	allocate(ielt_aux1(3,nel))
	allocate(ielt_aux2(3,nel))
	allocate(ielt_aux3(3,nel))
	allocate(ieltv_global(3,nel_global))
	allocate(aux1(3*nel))
	allocate(aux2(3*nel))

	ip = 0
	itot = 0
	iint = 0
	iu = 400 + my_id
	ies = 11311
	ies = 8667
	ies = 51799

	ielt_aux1 = 0
	do ie=1,nel
	  do ii=1,3
	    ien = ieltv(ii,ie)
	    if( ien > 0 ) ien = ipev(ien)
	    ielt_aux1(ii,ie) = ien
	  end do
	end do
	ielt_aux2 = ielt_aux1

	call shympi_exchange_fix_elem(3,ielt_aux1)

	nit = 0
	ie_tripple = 0
	do ie=1,nel
	  nid = id_elem(0,ie)			! how many domains
	  ide = id_elem(1,ie)			! main id of element
	  if( nid == 3 ) then			! adjust for tripple points
	    nit = nit + 1			! next tripple point
	    ie_ext = ipev(ie)
	    ie_tripple(:,nit) = ielt_aux2(:,ie)
	    tripple_info(1,nit) = ie_ext
	    tripple_info(2,nit) = ie
	    tripple_info(3,nit) = ii
	    write(iu,*) 'tripple: ',nit,my_id
	    write(iu,*) ie_tripple(:,nit)
	  end if
	end do

	ie_aux = reshape(ie_tripple,(/3*n_threads/))
	call shympi_array_reduce('max',ie_aux)		! only 1d array
	ie_tripple = reshape(ie_aux,(/3,n_threads/))

	do n=1,nit
	  write(iu,*) ie_tripple(:,n)
	  ie = tripple_info(2,n)
	  ielt_aux1(:,ie) = ie_tripple(:,n)
	end do

	do ie=1,nel
	  ie_ext = ipev(ie)
	  nid = id_elem(0,ie)			! how many domains
	  ide = id_elem(1,ie)			! main id of element
	  if( nid == 3 ) then			! adjust for tripple points
	    do ii=1,3
	      if( ielt_aux1(ii,ie) <= 1000 ) then
	        !if( ielt_aux2(ii,ie) <= 1000 ) goto 96
		ielt_aux1(ii,ie) = ielt_aux2(ii,ie)
	      end if
	    end do
	  end if
	if( ie_ext == ies ) then
	  write(iu,*) 'ie: ',ie_ext,ie,my_id
	  write(iu,*) 'aux1: ',ielt_aux1(:,ie)
	  write(iu,*) 'aux2: ',ielt_aux2(:,ie)
	end if
	end do

	call shympi_l2g_array(ielt_aux1,ieltv_global)
	do ie=1,nel_global
	  write(567,*) ie,ieltv_global(:,ie)
	end do

	iu = 400 + my_id
	write(iu,*) 'nel... ',nkn,nel,nel_unique
	do ie=1,nel
	  ie_ext = ipev(ie)
	  do ii=1,3
	    if( ielt_aux1(ii,ie) /= ielt_aux2(ii,ie) ) then
	      write(iu,*) ie_ext,ie,ii,ielt_aux1(ii,ie),ielt_aux2(ii,ie)
	    end if
	  end do
	end do

	do ie=1,nel
	  nid = id_elem(0,ie)			! how many domains
	  ide = id_elem(1,ie)			! main id of element
	  if( nid == 3 ) then
	    write(6,*) 'tripple point: ',ipev(ie),my_id,ide
	    write(iu,*) 'tripple point: ',ipev(ie),my_id,ide
	    cycle
	  end if
	  if( nid /= 2 ) cycle
	  ide2 = id_elem(2,ie)			! secondary id of element
	  do ii=1,3
	    ien = ieltv(ii,ie)
	    ie_ext = 0
	    if( ien > 0 ) ie_ext = ipev(ien)
	    buffer_in(ip+ii) = ie_ext
	  end do
	  buffer_in(ip+4) = ipev(ie)
	  buffer_in(ip+5) = ide
	  buffer_in(ip+6) = ide2
	  ip = ip + ndim
	  itot = itot + 1
	  if( my_id == ide ) iint = iint + 1
	end do

	call shympi_barrier
	flush(6)
	stop 'programmed stop in check_2id_elements'

	return
   96	continue
	write(6,*) ie,my_id
	write(6,*) ielt_aux1(:,ie)
	write(6,*) ielt_aux2(:,ie)
	stop 'error stop check_2id_elements: internal error (2)'
   98	continue
	iu = 400 + my_id
	write(iu,*) nkn,nel
	do ie=1,nel
	  ie_ext = ipev(ie)
	  do ii=1,3
	    if( ielt_aux1(ii,ie) /= ielt_aux3(ii,ie) ) then
	      write(iu,*) ie_ext,ie,ii,ielt_aux1(ii,ie),ielt_aux3(ii,ie)
	    end if
	  end do
	end do
	stop 'error stop check_2id_elements: internal error (1)'
   99	continue
	stop 'error stop check_2id_elements: cannot do tripple points yet'
	end

!*****************************************************************

	subroutine write_grd_domain

	use shympi
	use basin

	implicit none

	logical, parameter :: blocal = .false.		!write local domains
	integer nout
	integer k,ie,kext,itype,n,i
	real x,y,depth
	real, parameter :: flag = -999.
	real, allocatable :: xg(:),yg(:)
	integer, allocatable :: intype(:),ietype(:)
	integer, allocatable :: inext(:),ieext(:)
	integer, allocatable :: ieaux(:)
	integer, allocatable :: index(:,:)
	character*80 file,text
	character*5 cid

	if( .not. bmpi_debug ) return

	write(6,*) 'write_grd_domain:',my_id,nkn_global,size(id_node)

        allocate(xg(nkn_global))
        allocate(yg(nkn_global))
        allocate(index(3,nel_global))	!element index
        allocate(intype(nkn_global))	!node type
        allocate(ietype(nel_global))	!elem type
        allocate(inext(nkn_global))	!node external number
        allocate(ieext(nel_global))	!elem external number
        allocate(ieaux(nel))

	index = nen3v_global
	ieaux(:) = id_elem(1,:)		!this is main element id

	call shympi_l2g_array(xgv,xg)
	call shympi_l2g_array(ygv,yg)
	call shympi_l2g_array(id_node,intype)
	call shympi_l2g_array(ieaux,ietype)

	do k=1,nkn_global
	  inext(k) = ip_ext_node(k)
	end do
	do ie=1,nel_global
	  ieext(ie) = ip_ext_elem(ie)
	end do

!---------------------------------------------------------------
! write global grid
!---------------------------------------------------------------

	if( shympi_is_master() ) then
	  file = 'domain1.grd'
	  text = 'mpi domains'
	  call write_grd_file(file,text,nkn_global,nel_global,xg,yg &
     &				,index,inext,ieext,intype,ietype)
	end if

	where( id_elem(0,:) > 1 ) ieaux(:) = -1	!two/three domain elem
	call shympi_l2g_array(ieaux,ietype)

	if( shympi_is_master() ) then
	  file = 'domain2.grd'
	  text = 'mpi domains 2'
	  call write_grd_file(file,text,nkn_global,nel_global,xg,yg &
     &				,index,inext,ieext,intype,ietype)
	end if

	call shympi_barrier

!---------------------------------------------------------------
! write local grid
!---------------------------------------------------------------

	if( .not. blocal ) return

	write(cid,'(i5)') my_id
	do i=1,5
	  if( cid(i:i) == ' ' ) cid(i:i) = '0'
	end do

	do k=1,nkn
	  inext(k) = ipv(k)
	end do
	do ie=1,nel
	  ieext(ie) = ipev(ie)
	end do
	intype(1:nkn) = id_node(1:nkn)
	ieaux(:) = 0					!main element id
	where( id_elem(0,:) == 2 ) ieaux(:) = 2		!two domain elem
	where( id_elem(0,:) == 3 ) ieaux(:) = 3		!three domain elem 
	ietype(1:nel) = ieaux(1:nel)
	index(:,1:nel) = nen3v(:,1:nel)

	file = 'domain_' // cid // '.grd'
	!write(6,*) file
	text = 'local mpi domain'

	call write_grd_file(file,text,nkn,nel,xgv,ygv &
     &				,index,inext,ieext,intype,ietype)

!---------------------------------------------------------------
! end of routine
!---------------------------------------------------------------

	call shympi_barrier

	end

!*****************************************************************

	subroutine write_grd_general(file,text &
     &				,intype,ietype,rndepth,redepth)

! writes node and element arrays as type values to grd

	use shympi
	use basin

	implicit none

	character*(*) file,text
	integer :: intype(nkn),ietype(nel)
	real :: rndepth(nkn),redepth(nel)

	integer nout
	integer k,ie,kext,itype,n,i
	real x,y,depth
	real, parameter :: flag = -999.
	real, allocatable :: xg(:),yg(:)
	integer, allocatable :: inext(:),ieext(:)
	integer, allocatable :: index(:,:)
	integer, allocatable :: ingtype(:),iegtype(:)
	real, allocatable :: rngdepth(:),regdepth(:)

	write(6,*) 'write_grd_general:',my_id,nkn_global,size(id_node)

        allocate(xg(nkn_global))
        allocate(yg(nkn_global))
        allocate(index(3,nel_global))	!element index
        allocate(inext(nkn_global))	!node external number
        allocate(ieext(nel_global))	!elem external number

        allocate(ingtype(nkn_global))	!node type
        allocate(iegtype(nel_global))	!elem type
        allocate(rngdepth(nkn_global))	!node depth
        allocate(regdepth(nel_global))	!elem depth

	index = nen3v_global

	call shympi_l2g_array(xgv,xg)
	call shympi_l2g_array(ygv,yg)
	call shympi_l2g_array(intype,ingtype)
	call shympi_l2g_array(ietype,iegtype)
	call shympi_l2g_array(rndepth,rngdepth)
	call shympi_l2g_array(redepth,regdepth)

	do k=1,nkn_global
	  inext(k) = ip_ext_node(k)
	end do
	do ie=1,nel_global
	  ieext(ie) = ip_ext_elem(ie)
	end do

!---------------------------------------------------------------
! write global grid
!---------------------------------------------------------------

	if( shympi_is_master() ) then
	  call write_grd_file_with_depth(file,text &
     &				,nkn_global,nel_global,xg,yg &
     &				,index,inext,ieext,ingtype,iegtype &
     &				,rndepth,regdepth)
	end if

!---------------------------------------------------------------
! end of routine
!---------------------------------------------------------------

	call shympi_barrier

	end

!*****************************************************************
