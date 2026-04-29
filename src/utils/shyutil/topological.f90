

!--------------------------------------------------------------------------
!
!    Copyright (C) 2003,2008-2009,2011,2015-2019  Georg Umgiesser
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

! topological set up routines
!
! contents :
!
! subroutine setnod
!			 sets array inodv
!
! revision log :
!
! 01.08.2003	ggu	created from sublnk.f
! 13.08.2003	ggu	in update_geom do not call setweg and setnod
! 06.11.2008	ggu	better error handling
! 06.04.2009	ggu	nlidim -> nlkdim
! 23.03.2010	ggu	changed v6.1.1
! 02.12.2011	ggu	print_bound_nodes() for writing out boundary nodes
! 09.12.2011	ggu	changed VERS_6_1_38
! 30.03.2012	ggu	changed VERS_6_1_51
! 12.09.2013	ggu	changed VERS_6_1_67
! 18.06.2014	ggu	changed VERS_6_1_77
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 04.05.2015	ggu	remove equivalence - use winkv as local array
! 20.05.2015	ggu	new call to mklenkii to set up lenkiiv
! 10.07.2015	ggu	changed VERS_7_1_50
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 24.07.2015	ggu	changed VERS_7_1_82
! 10.10.2015	ggu	changed VERS_7_3_2
! 16.12.2015	ggu	changed VERS_7_3_16
! 28.04.2016	ggu	changed VERS_7_5_9
! 26.09.2017	ggu	changed VERS_7_5_32
! 05.12.2017	ggu	changed VERS_7_5_39
! 24.01.2018	ggu	changed VERS_7_5_41
! 19.04.2018	ggu	changed VERS_7_5_45
! 16.02.2019	ggu	changed VERS_7_5_60
! 20.05.2020	ggu	new way to compute link structure (still experimental)
! 28.05.2020	ggu	some more changes in constructing link structure
! 13.04.2022	ggu	new call to make_links (ibound)
! 07.11.2024	ggu	new framework for connection
! 22.04.2026	ggu	make clear that errors can be ignored
!
!*****************************************************************

	subroutine set_geom

! sets up geometrical arrays

	use basin
	use mod_geom
	use mod_connect

        implicit none

	logical bverbose
        integer i,n
        integer nli,nbn,nin,nis,ngrd,ngrd1,ngrd2
	integer ibound(nkn)
        integer kerr

	logical belem
	integer ibase,ipf,ipl
	integer k,ndim,kext
	integer ne,ne1,ne2
	integer nk,nk1,nk2
	integer ngr_local
	integer ierr
	integer, allocatable :: nnlist(:,:)
	integer, allocatable :: eelist(:,:)
	integer, allocatable :: ecv(:,:)
	integer, allocatable :: bound(:)
	integer, allocatable :: kant(:,:)
	integer, allocatable :: neaux(:)
	integer, allocatable :: kerror(:)

	integer ipint

	bverbose = .true.
	bverbose = .false.

!-------------------------------------------------------------
! make static arrays - new way
!-------------------------------------------------------------

        call connect_init(nkn,nel,nen3v,ierr)	!this computes internal arrays

	if( ierr > 0 ) then
	  if( .not. allocated(kerror) ) allocate(kerror(ierr))
	  call connect_errors(ierr,kerror)
	  write(6,*) '(errors can be ignored)'
	end if
 
        call connect_get_grade(ngr_local)
	if( ngr /= ngr_local ) ngr = ngr_local	!insert new value into basin
	allocate(nnlist(0:ngr,nkn))
	allocate(eelist(0:ngr,nkn))
	allocate(ecv(3,nel))
	allocate(bound(nkn))
	allocate(kant(2,nkn))
        call connect_get_lists(nkn,ngr,nnlist,eelist)
        call connect_get_ecv(nel,ecv)
        call connect_get_bound(nkn,bound)
        call connect_get_kant(nkn,kant)

	call connect_release

!-------------------------------------------------------------
! allocate arrays
!-------------------------------------------------------------

	call mod_geom_init(nkn,nel,ngr)

	nlist = nnlist
	elist = eelist

!-------------------------------------------------------------
! check maxlnk
!-------------------------------------------------------------

	if( ngr .gt. maxlnk ) goto 98

!-------------------------------------------------------------
! make static arrays - old way
!-------------------------------------------------------------

	ieltv = ecv
	iboundv = bound
	kantv = kant

	if( buselink ) then
	  call make_links_old(nkn,nel,nen3v)
	  call make_links(nkn,nel,nen3v,ibound,kerr)
	end if

!-------------------------------------------------------------
! check results
!-------------------------------------------------------------

	ndim = ngr
	allocate(neaux(ndim))

	do k=1,nkn
	  nk = nlist(0,k)
	  ne = elist(0,k)
	  belem = .false.
	  call get_nodes_around(k,ndim,nk2,neaux)
	  nk1 = nk2
	  if( buselink ) then
	    call get_node_linkp(k,ipf,ipl)
	    call get_node_links(k,nk1,ibase)
	    if( nk /= ipl-ipf+1 ) goto 89
	  end if
	  if( nk /= nk1 .or. nk /= nk2 ) goto 89
	  if( any( nlist(1:nk,k) /= neaux(1:nk) ) ) goto 88
	  belem = .true.
	  call get_elems_around(k,ndim,ne2,neaux)
	  ne1 = ne2
	  if( buselink ) then
	    call get_elem_linkp(k,ipf,ipl)
	    call get_elem_links(k,ne1,ibase)
	    if( ne /= ipl-ipf+1 ) goto 89
	  end if
	  if( ne /= ne1 .or. ne /= ne2 ) goto 89
	  if( any( elist(1:ne,k) /= neaux(1:ne) ) ) goto 88
	end do

	if( any( ecv /= ieltv ) ) goto 87
	if( any( bound /= iboundv ) ) goto 86
	if( any( kant /= kantv ) ) goto 85

	if( bverbose ) then
	  write(6,*) 'set_geom: all connection compatibility checks' // &
     &				' successfully completed'
	end if

!-------------------------------------------------------------
! write some statistics
!-------------------------------------------------------------

	if( bverbose ) then
	  write(6,*) 'set_geom: statistics'
	end if

	nbn = 0
	ngrd = 0
	if( buselink ) then
          ngrd=ilinkv(nkn+1)
          n=0
          do i=1,ngrd
            if(lenkv(i).eq.0) n=n+1
          end do
	  nbn = n
	end if

	n = 0
	do k=1,nkn
	  if( elist(0,k) /= nlist(0,k) ) n = n + 1
	end do
	if( .not. buselink ) nbn = n
	if( n /= nbn ) then
	  write(6,*) n,nbn
	  stop 'error stop set_geom: n/=nbn'
	end if

	nin = nkn - nbn
	nis = (nel-2*nkn+nbn+2)/2

	if( bverbose ) then
          write(6,*) 'nel                      : ',nel
          write(6,*) 'nkn                      : ',nkn
          write(6,*) 'nbn                      : ',nbn
          write(6,*) 'nin                      : ',nin
          write(6,*) 'nis                      : ',nis
	end if

	ngrd1 = 4*nbn+6*nin+6*(nis-1)
        ngrd2 = 3*nel+nbn
	if( .not. buselink ) ngrd = ngrd1
	nli = ngrd/2

	if( bverbose ) then
          write(6,*) 'ngr                      : ',ngr
          write(6,*) 'grades                   : ',ngrd
          write(6,*) 'formula (nis)            : ',ngrd1
          write(6,*) 'formula (nel)            : ',ngrd2
          write(6,*) 'links                    : ',nli
	end if

	if( ngrd .ne. ngrd1 ) goto 99
	if( ierr == 0 .and. ngrd .ne. ngrd2 ) goto 99

	if( bverbose ) then
	  write(6,*) 'set_geom: all statistics checks successfully completed'
	end if

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	return
   85	continue
	stop 'error stop set_geom: kant'
   86	continue
	stop 'error stop set_geom: bound'
   87	continue
	stop 'error stop set_geom: ecv'
   88	continue
	if( belem ) then
	  write(6,*) elist(1:ne,k)
	  write(6,*) neaux(1:ne)
	  stop 'error stop set_geom: elist'
	else
	  write(6,*) nlist(1:ne,k)
	  write(6,*) neaux(1:ne)
	  stop 'error stop set_geom: nlist'
	end if
   89	continue
	if( belem ) then
	  write(6,*) ne,ne1,ne2
	  stop 'error stop set_geom: ne'
	else
	  write(6,*) nk,nk1,nk2
	  stop 'error stop set_geom: nk'
	end if
   98	continue
	write(6,*) 'ngr,maxlnk: ',ngr,maxlnk
	stop 'error stop set_geom: maxlnk'
   99	continue
	write(6,*) ngrd,ngrd1,ngrd2,ierr
	stop 'error stop set_geom: ngrade'
	end

!*****************************************************************

	subroutine check_geom

! checks geometrical arrays

	use mod_geom
	use basin

        implicit none

!-------------------------------------------------------------
! check static arrays
!-------------------------------------------------------------

	if( buselink ) then
          call checklenk(nkn,ilinkv,lenkv,nen3v)
          call checklink(nkn,ilinkv,linkv)
	end if

        call checkkant(nkn,kantv)
        call checkielt(nel,ieltv,nen3v)

	call check_subs

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	end

!*****************************************************************

        subroutine check_subs

! checks various subroutines

	use mod_geom
	use basin, only : nkn,nel,ngr,mbw

        implicit none

	integer k,n,ne,ipf,ipl,ibase
	integer nl,i,ip,ie1,ie2,nn,k1,k2
	integer nodes(maxlnk)
	integer elems(maxlnk)

	if( .not. buselink ) return

!-------------------------------------------------------------
! check element index
!-------------------------------------------------------------

	do k=1,nkn
	  call get_elem_linkp(k,ipf,ipl)
	  call get_elem_links(k,n,ibase)
	  call get_elems_around(k,maxlnk,ne,elems)
	  nl = ipl-ipf+1
	  if( n .ne. ne .or. nl .ne. ne ) goto 99
	  if( ipf .ne. ibase + 1 ) goto 96
	  if( ipl .ne. ibase + n ) goto 96
	  do i=1,n
	    ip = ipf + i - 1
	    if( ibase + i .ne. ip ) goto 98
	    ie1 = lenkv(ibase+i)
	    ie2 = elems(i)
	    if( ie1 .ne. ie2 ) goto 97
	  end do
	end do

!-------------------------------------------------------------
! check node index
!-------------------------------------------------------------

	do k=1,nkn
	  call get_node_linkp(k,ipf,ipl)
	  call get_node_links(k,n,ibase)
	  call get_nodes_around(k,maxlnk,nn,nodes)
	  nl = ipl-ipf+1
	  if( n .ne. nn .or. nl .ne. nn ) goto 89
	  if( ipf .ne. ibase + 1 ) goto 86
	  if( ipl .ne. ibase + n ) goto 86
	  do i=1,n
	    ip = ipf + i - 1
	    if( ibase + i .ne. ip ) goto 88
	    k1 = linkv(ibase+i)
	    k2 = nodes(i)
	    if( k1 .ne. k2 ) goto 87
	  end do
	end do

!-------------------------------------------------------------
! end of routine
!-------------------------------------------------------------

	return
   86	continue
	write(6,*) k
	write(6,*) 'ipf,ipl,ibase,n: ',ipf,ipl,ibase,n
	stop 'error stop check_subs: extreme base pointers'
   87	continue
	write(6,*) k
	write(6,*) 'k1,k2: ',k1,k2
	stop 'error stop check_subs: node numbers'
   88	continue
	write(6,*) k
	write(6,*) 'i: ',i,ip,ibase + i
	stop 'error stop check_subs: base pointer'
   89	continue
	write(6,*) k
	write(6,*) 'n,nn,nl: ',n,nn,nl
	stop 'error stop check_subs: total number of nodes'
   96	continue
	write(6,*) k
	write(6,*) 'ipf,ipl,ibase,n: ',ipf,ipl,ibase,n
	stop 'error stop check_subs: extreme base pointers'
   97	continue
	write(6,*) k
	write(6,*) 'ie1,ie2: ',ie1,ie2
	stop 'error stop check_subs: element numbers'
   98	continue
	write(6,*) k
	write(6,*) 'i: ',i,ip,ibase + i
	stop 'error stop check_subs: base pointer'
   99	continue
	write(6,*) k
	write(6,*) 'n,ne,nl: ',n,ne,nl
	stop 'error stop check_subs: total number of elements'
	end

!****************************************************************
!****************************************************************
!****************************************************************

        subroutine print_bound_nodes(nkn,aux)

	use mod_geom

        implicit none

        integer nkn
        real aux(nkn)

        integer ib,k,kn,kstart

        integer ipext

        ib = 0
        do k=1,nkn
          aux(k) = 0.
        end do

        open(1,file='bnd_nodes.dat',form='formatted',status='unknown')
        do k=1,nkn
          if( aux(k) .eq. 0. ) then
            if( kantv(1,k) .gt. 0 ) then
              kstart = k
              kn = kantv(1,kstart)
              ib = ib + 1
              write(1,*) 'bound ',ib
              write(1,*) ipext(kstart)
              do while( kn .ne. kstart )
                aux(kn) = 1.
                write(1,*) ipext(kn)
                kn = kantv(1,kn)
              end do
            end if
          end if
          aux(k) = 1.
        end do

        close(1)

        end 

!****************************************************************

