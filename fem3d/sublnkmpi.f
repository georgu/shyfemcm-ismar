
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
!
!*****************************************************************

	subroutine set_geom_mpi

! sets up geometrical arrays
!
! geom structures have already been setup on local domain
! with this call they are corrected with gloval boundary node information

	use mod_geom
	use basin
	use shympi

        implicit none

        integer i,n
	integer k,kk,kn,kb,k1,k2
	integer ie,ii
	integer id,id_neigh
	integer iunit
	integer, allocatable :: ibound(:)
        integer kerr

	integer knext,kbhnd

	if( .not. bmpi ) return

	iunit = 660 + my_id

!-------------------------------------------------------------
! make global ibound array
!-------------------------------------------------------------

	call write_grd_domain

	!return
	!write(iunit,*) '-------- set_geom_mpi ----------'
	!write(iunit,*) 'global: ',nkn_global,nel_global
	!write(iunit,*) 'local: ',nkn,nel
	allocate(ibound(nkn_global))
	!call link_set_wmpi(.true.)
	call make_links(nkn_global,nel_global,nen3v_global,ibound,kerr)
	!call link_set_wmpi(.false.)
	!stop

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
	  if( iboundv(k) == 0 ) then
	    kantv(:,k) = 0
	  else if( .not. shympi_is_inner_node(k) ) then
	    do i=1,2
	      kk = kantv(i,k)
	      id = id_node(kk)
	      if( id /= my_id ) kantv(i,k) = 0
	    end do
	  end if
	end do

!-------------------------------------------------------------
! adjust ieltv (-1 is reserved for OB)
!-------------------------------------------------------------

	do ie=1,nel
	  do ii=1,3
	    if( ieltv(ii,ie) /= 0 ) cycle	! has neighbor
	    k = nen3v(ii,ie)
	    kn = knext(k,ie)
	    kb = kbhnd(k,ie)
	    if( iboundv(kn) == 0 .or. iboundv(kb) == 0 ) then
	      id_neigh = id_elem(1,ie)
	      if( id_neigh == -1 ) goto 99
	      ieltv(ii,ie) = -1000 - id_neigh	!internal (not real) border
	    end if
	  end do
	end do

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	return
   99	continue
	stop 'error stop set_geom_mpi: internal error (1)'
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
	ieaux(:) = id_elem(0,:)		!this is main element id

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
	  call write_grd_file(file,text,nkn_global,nel_global,xg,yg
     +				,index,inext,ieext,intype,ietype)
	end if

	where( id_elem(1,:) /= -1 ) ieaux(:) = -1	!two/three domain elem
	call shympi_l2g_array(ieaux,ietype)

	if( shympi_is_master() ) then
	  file = 'domain2.grd'
	  text = 'mpi domains 2'
	  call write_grd_file(file,text,nkn_global,nel_global,xg,yg
     +				,index,inext,ieext,intype,ietype)
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
	where( id_elem(1,:) /= -1 ) ieaux(:) = 2	!two/three domain elem
	where( id_elem(2,:) /= -1 
     +		.and. id_elem(1,:) /= id_elem(2,:) ) 
     +				ieaux(:) = 3 		!three domain elem
	ietype(1:nel) = ieaux(1:nel)
	index(:,1:nel) = nen3v(:,1:nel)

	file = 'domain_' // cid // '.grd'
	!write(6,*) file
	text = 'local mpi domain'

	call write_grd_file(file,text,nkn,nel,xgv,ygv
     +				,index,inext,ieext,intype,ietype)

!---------------------------------------------------------------
! end of routine
!---------------------------------------------------------------

	call shympi_barrier

	end

!*****************************************************************

