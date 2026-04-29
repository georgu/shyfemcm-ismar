
!--------------------------------------------------------------------------
!
!    Copyright (C) 1999,2003-2004,2008-2009,2011,2013  Georg Umgiesser
!    Copyright (C) 2018-2019  Georg Umgiesser
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
! 06.04.1999	ggu	completely restructured
! 04.06.1999	ggu	new statistics are computed
! 08.09.2003	ggu	mode 5 -> write depth values from elements
! 23.09.2004	ggu	interpolq() changed for bathy interpolation
! 02.10.2004	ggu	interpole() for exponential interpolation
! 01.11.2004	ggu	whole program simplyfied
! 06.12.2008	ggu	smoothing introduced
! 06.04.2009	ggu	read param.h
! 29.05.2009	ggu	does only depth limiting and smoothing
! 20.11.2009	ggu	possibility to smooth only on specific areas
! 30.03.2011	ggu	new routines to delete elements
! 13.06.2013	ggu	copy_depth() renamed to transfer_depth()
! 16.04.2018	ggu	partitioning finished
! 20.04.2018	ggu	code added to check partition integrity
! 16.02.2019	ggu	changed VERS_7_5_60
! 28.05.2020	ggu	new checks for connection
! 30.03.2022	ggu	bug fix: nkn_save,nel_save were not initialized
! 13.04.2022    ggu     new call to make_links (ibound)
! 04.04.2023    ggu     minor changes
! 07.11.2024    ggu     ignore connection errors
! 21.11.2024    ggu     call check routines with bdebug, less error messages
! 23.11.2024    ggu     ignore errors in connectivity
! 30.11.2024    ggu     renamed check_color() to check_part_color()
! 22.04.2026    ggu     documentation integrated
!
!****************************************************************

!================================================================
	module mod_save_index
!================================================================

	implicit none

	integer, save, private :: nkn_save = 0
	integer, save, private :: nel_save = 0
	integer, save, private :: ngr_save = 0

	integer, allocatable :: save_ipv(:)
	integer, allocatable :: save_ipev(:)
	integer, allocatable :: save_nen3v(:,:)
	real, allocatable    :: save_xgv(:)
	real, allocatable    :: save_ygv(:)
	integer, allocatable :: save_iarv(:)
	integer, allocatable :: save_iarnv(:)

        INTERFACE copy_array
        MODULE PROCEDURE  copy_array_i &
     &                   ,copy_array_r
        END INTERFACE

!================================================================
	contains
!================================================================

	subroutine mod_save_index_init(nkn,nel,ngr)

	implicit none

        integer nkn, nel, ngr

        if( nkn == nkn_save .and. nel == nel_save .and. ngr == ngr_save ) return

        if( nel > 0 .or. nkn > 0 .or. ngr > 0 ) then
          if( nel == 0 .or. nkn == 0 .or. ngr == 0 ) then
            write(6,*) 'nel,nkn,ngr: ',nel,nkn,ngr
            stop 'error stop mod_save_index: incompatible parameters'
          end if
        end if

	call release_save_index

	nkn_save = nkn
	nel_save = nel
	ngr_save = ngr

	if( nkn == 0 ) return

	allocate(save_ipv(nkn))
	allocate(save_ipev(nel))
	allocate(save_nen3v(3,nel))
	allocate(save_xgv(nkn))
	allocate(save_ygv(nkn))
	allocate(save_iarv(nel))
	allocate(save_iarnv(nkn))

	end subroutine mod_save_index_init

!****************************************************************

	subroutine release_save_index

	implicit none

        if( nkn_save > 0 ) then
	  deallocate(save_ipv)
	  deallocate(save_ipev)
	  deallocate(save_nen3v)
	  deallocate(save_xgv)
	  deallocate(save_ygv)
	  deallocate(save_iarv)
	  deallocate(save_iarnv)
	end if

	nkn_save = 0
	nel_save = 0
	ngr_save = 0

	end subroutine

!****************************************************************

	subroutine make_new_basin(nk,ne,nenv,nodep,elemp)

! uses nenv as new element index and copies arrays to new basin

	use basin

	implicit none

	integer nk,ne
	integer nenv(3,nel)
	integer nodep(nk)
	integer elemp(ne)

	integer k,ie,ng,nn

	if( nk > nkn .or. ne > nel ) then
	  write(6,*) 'can only make smaller index: ',nk,nkn,ne,nel
	  stop 'error stop make_new_basin: nk/ne > nkn/nel'
	end if

	call mod_save_index_init(nkn,nel,ngr)

	save_ipv = ipv
	save_ipev = ipev
	save_nen3v = nen3v
	save_xgv = xgv
	save_ygv = ygv

	nen3v(:,1:ne) = nenv(:,1:ne)

	call copy_array(nk,ipv,nodep)
	call copy_array(ne,ipev,elemp)

	call copy_array(nk,xgv,nodep)
	call copy_array(nk,ygv,nodep)

	call copy_array(ne,iarv,elemp)
	call copy_array(nk,iarnv,nodep)

	!do k=1,nk
	!  aux_ipv(k) = ipv(nodep(k))
	!end do
	!ipv(1:nk) = aux_ipv(1:nk)

	!do ie=1,ne
	!  aux_ipev(ie) = ipev(elemp(ie))
	!end do
	!ipev(1:ne) = aux_ipev(1:ne)

	call estimate_grade(nk,ne,nenv,ng)

	nkn = nk
	nel = ne
	ngr = ng

	end subroutine

!****************************************************************

	subroutine transfer_array(n,array,ipointer,array_new)

	implicit none

	integer n
	integer array(:)
	integer ipointer(n)
	integer array_new(n)

	integer i
	integer aux(n)

	do i=1,n
	  aux(i) = array(ipointer(i))
	end do
	array_new(1:n) = aux(1:n)

	end

!****************************************************************

	subroutine copy_array_i(n,array,ipointer)

	implicit none

	integer n
	integer array(:)
	integer ipointer(n)

	integer i
	integer aux(n)

	do i=1,n
	  aux(i) = array(ipointer(i))
	end do
	array(1:n) = aux(1:n)

	end

!****************************************************************

	subroutine copy_array_r(n,array,ipointer)

	implicit none

	integer n
	real array(:)
	integer ipointer(n)

	integer i
	real aux(n)

	do i=1,n
	  aux(i) = array(ipointer(i))
	end do
	array(1:n) = aux(1:n)

	end

!****************************************************************

	subroutine restore_old_basin

	use basin

	implicit none

	nkn = nkn_save
	nel = nel_save
	ngr = ngr_save

	nen3v = save_nen3v
	ipv = save_ipv
	ipev = save_ipev
	xgv = save_xgv
	ygv = save_ygv

	end subroutine

!****************************************************************

	subroutine estimate_grade(nk,ne,nenv,ng)

	implicit none

	integer nk,ne
	integer nenv(3,ne)
	integer ng

	integer ie,ii,k
	integer, allocatable :: node(:)

	allocate(node(nk))
	node = 0

	do ie=1,ne
	  do ii=1,3
	    k = nenv(ii,ie)
	    node(k) = node(k) + 1
	  end do
	end do

	ng = 1 + maxval(node)

	end subroutine

!================================================================
	end module mod_save_index
!================================================================

        subroutine bas_partition(llfile)

! performs partition on basin - is called through shybas, not inside shyfem

	use mod_geom
	use mod_depth
	use basin
	use grd
	use shympi

	implicit none

	character*(*) llfile

	logical bdebug
	integer ierr
	integer k,i,nl,il,n,ib,in,node,ilext,np,ic
	integer, allocatable :: nc(:)
	real x,y,perc
	real xx(nkn)
	real yy(nkn)

	logical inconvex,inpoly,filex

!-----------------------------------------------------------------
! open and read file containing lines
!-----------------------------------------------------------------

	bdebug = .false.

	if( .not. filex(llfile) ) then
	  write(6,*) 'Cannot open file ',trim(llfile)
	  stop 'error stop bas_partition: no such file'
	end if

	call grd_read(llfile)

	iarnv = 0
	nl = nl_grd

!-----------------------------------------------------------------
! loop over lines
!-----------------------------------------------------------------

        do i=1,nl
          il = i
          ilext = ipplv(il)
          n = ipntlv(il) - ipntlv(il-1)
          ib = ipntlv(il-1)
	  !write(6,*) 'line: ',i,ilext,n
	  do in=1,n
	    node = inodlv(ib+in)
	    x = xv(node)
	    y = yv(node)
	    xx(in) = x
	    yy(in) = y
	    !write(6,*) in,node,x,y
	  end do
	  if( xx(1) == xx(n) .and. yy(1) == yy(n) ) n = n - 1
	  np = 0
	  do k=1,nkn
	    x = xgv(k)
	    y = ygv(k)
	    !write(6,*) 'looking for... ',k,x,y
	    if( inpoly(n,xx,yy,x,y) ) then
	      iarnv(k) = il
	      np = np + 1
	      !write(6,*) 'inside... ',il,n,k,x,y
	    end if
	  end do
	  perc = (100.*np)/nkn
	  !write(6,*) 'inside points found... ',i,il,n,np,perc
        end do

	if( count( iarnv == 0 ) > 0 ) then	!not handled nodes
	  nl = nl + 1
	  where( iarnv == 0 ) iarnv = nl
	end if

!-----------------------------------------------------------------
! write information to terminal
!-----------------------------------------------------------------

	allocate(nc(0:nl))
	nc = 0
	do k=1,nkn
	  ic = iarnv(k)
	  if( ic < 1 .or. ic > nl ) then
	    write(6,*) 'ic,nl: ',ic,nl
	    stop 'error stop bas_partition: internal error (1)'
	  end if
	  nc(ic) = nc(ic) + 1
	end do
	write(6,*) 'Information on domains: '
	write(6,*) '   domain      area     nodes   percent'
	do ic=1,nl
	  write(6,'(3i10,f10.2)') ic-1,ic,nc(ic),(100.*nc(ic))/nkn
	end do

!-----------------------------------------------------------------
! write file
!-----------------------------------------------------------------

        call basin_to_grd

	if( shympi_is_master() ) then
          call grd_write('bas_partition.grd')
          write(6,*) 'The partition has been written to' // &
     &				' bas_partition.grd'
	end if

	call shympi_syncronize
	call check_connectivity(bdebug,ierr)
	call shympi_syncronize
	call check_connections(bdebug,ierr)
	call shympi_syncronize

!-----------------------------------------------------------------
! end of routine
!-----------------------------------------------------------------

	return
	end

!*******************************************************************

	subroutine check_connectivity(bdebug,ierr)

	use basin

	implicit none

	logical bdebug
	integer ierr

	integer ier
	integer ic,nc
	integer icolor(nkn)

	ierr = 0
	icolor = iarnv
	nc = maxval(icolor)

	do ic=0,nc
	  ier = 0
	  call check_part_color(ic,nkn,icolor,ier)
	  ierr = ierr + ier
	end do

	end

!*******************************************************************

	subroutine check_part_color(ic,n,icolor,ierr)

! checks number of doamins with color ic
!
! writes out error if not contigous
! if ierr == -1 on entry does not write error
! on return in ierr is number of sub-domains

	use basin

	implicit none

	integer ic
	integer n
	integer icolor(n)
	integer ierr

	integer cc,i,nfound
	logical bwrite
	logical, save :: bfirst = .true.
	character*80 header

	bwrite = ( ierr /= -1 )
	ierr = 0
	cc = count( icolor == ic )

	do
	  do i=1,n
	    if( icolor(i) == ic ) exit		!start from this node
	  end do
	  if( i > n ) exit
	  call flood_fill_bas(i,n,icolor,nfound)
	  if( nfound /= cc ) then
	    if( bfirst .and. bwrite ) then
	      header = '                               area' &
     &			// '     nthis    ntotal      node'
	      bfirst = .false.
	      write(6,'(a)') header
	    end if
	    if( bwrite ) then
	      write(6,1000) '  area is not connected: ',ic,nfound,cc,ipv(i)
	    end if
	    ierr = ierr + 1
	  end if
	end do

	return
 1000	format(a,4i10)
	end

!*******************************************************************

	subroutine flood_fill_bas(i,n,icolor,nfound)

	use basin

	implicit none

	integer i
	integer n
	integer icolor(n)
	integer nfound

	integer ic,nf,ie,ii,k
	logical bcol,bdone

	ic = icolor(i)
	icolor(i) = -1
	nfound = 1

	do
	  nf = 0
	  do ie=1,nel
	    bcol = .false.
	    bdone = .false.
	    do ii=1,3
	      k = nen3v(ii,ie)
	      if( icolor(k) == ic ) bcol = .true.
	      if( icolor(k) == -1 ) bdone = .true.
	    end do
	    if( bcol .and. bdone ) then
	      do ii=1,3
	        k = nen3v(ii,ie)
	        if( icolor(k) == ic ) then
	          nf = nf + 1
		  icolor(k) = -1
		end if
	      end do
	    end if
	  end do
	  if( nf == 0 ) exit
	  nfound = nfound + nf
	end do
	
	where( icolor == -1 ) icolor = -2

	end

!*******************************************************************

	subroutine check_connections(bdebug,kerr)

! this checks connections with link/lenk data structure

	use basin
	use mod_save_index
	use shympi

	implicit none

	logical bdebug
	integer kerr

	logical bloop
	logical bwrite
	logical bwritegrd
	logical berror
	logical bhandle_errors
	integer nloop
	integer ic,nc,ncol,kext,iu,i,k
	integer nk,ne
	integer nenv(3,nel)
	integer icolor(nkn)
	integer icolor_new(nkn)
	integer icol(nkn)
	integer nodep(nkn)
	integer elemp(nel)
	character*80 grdfile
	integer, allocatable :: kerror(:)

	integer ipint,ipext

	bwrite = .false.
	bwrite = .true.
	berror = .false.
	bwritegrd = .false.
	bwritegrd = .true.
	bhandle_errors = .false.
	bloop = .true.
	nloop = 0

	bwritegrd = bdebug	!only write error grid if in debug mode
	bwrite = bdebug

	allocate(kerror(nkn))

	icolor = iarnv
	nc = maxval(icolor)

	iu = 0
	!if( shympi_is_master() ) iu = 777

	if( shympi_is_master() ) then
	write(6,*) '========================================'
	write(6,*) 'checking total domain... total domains = ',nc
	write(6,*) '========================================'
	end if

	kerr = nkn
	call check_elem_index(nkn,nel,nen3v,kerr,kerror)
	if( kerr > 0 ) then
	  write(6,*) 'global domain has connection errors...'
	  stop 'error stop check_connections: irregular domain'
	end if

	!if( shympi_is_master() ) then
	!  if( bwritegrd ) then
	!    call write_partition_grd(icolor,0)
	!  end if
	!end if

	call shympi_barrier

	!return

!---------------------------------------------
! loop on domains
!---------------------------------------------

	if( shympi_is_master() ) then

	do while( bloop )

	nloop = nloop + 1
	if( nloop > 10 ) exit
	if( nloop > 1 ) exit

	do ic=1,nc
	  ncol = count( icolor == ic )
	  if( ncol == 0 ) cycle

	  if( bwrite ) write(6,*) 'checking domain ',ic,ncol

	  call make_elem_index(.true.,ic,icolor &
     &				,nk,ne,nenv,nodep,elemp)
	  call make_new_basin(nk,ne,nenv,nodep,elemp)
	  call transfer_array(nk,icolor,nodep,icolor_new)
	  call check_elem_index(nk,ne,nenv,kerr,kerror)

	  if( shympi_is_master() ) then
	    if( kerr > 0 ) then
	      write(6,*) 'errors found: ',kerr
	      write(6,'(a)') '    ierr    kint    kext    icol' // &
     &			 '                   x                   y'
	    end if
	    do i=1,kerr
	      k = kerror(i)
	      kext = ipv(k)
	      write(6,1000) i,k,kext,icolor_new(k),xgv(k),ygv(k)
 1000	      format(4i8,2f20.8)
	      if( bwritegrd ) then
	        call write_partition_grd(icolor_new,kext)
	      end if
	    end do
	    if( kerr > 0 ) write(6,*) '(errors can be ignored)'
	  end if
	  !if( kerr > 0 ) berror = .true.

	  !if( bwrite ) write(6,*) 'finished checking domain ',ic

	  !-------------------------------------------
	  ! handle errors
	  !-------------------------------------------

	  kext = 0
	  if( bhandle_errors .and. kerr /= 0 ) then
	    kext = ipext(kerr)
	    !write(6,*) 'adjusting node kerr = ',kerr,kext,ic
	    call restore_old_basin
	    call adjust_domain(ic,nkn,icolor_new,kext)
	    exit
	  end if

	  call restore_old_basin
	end do

	if( bhandle_errors .and. kerr > 0 ) then
	  if( iu > 0 ) write(iu,*) 'check_connections: ',nloop,kerr,kext
	  call elim_isolated_element(icolor,kerr)
	end if

	if( berror ) stop 'error stop check_connections: connections'
	if( ic > nc ) bloop = .false.
	end do	!do while(bloop)

	end if

	!kerr = 0

!---------------------------------------------
! end of loop on domains
!---------------------------------------------

	if( kerr == 0 ) then
	  if( nloop > 1 ) then
	    write(6,*) 'all domains have been corrected...',nloop
	  else
	    !write(6,*) 'no problems found in domains'
	  end if
	else
	  !write(6,*) 'could not correct error in connections...',nloop
	  !write(6,*) 'error occurred in domain ',my_id
	  !call write_partition_grd(icolor,0)
	end if

	iarnv = icolor

	call shympi_syncronize

!---------------------------------------------
! end of routine
!---------------------------------------------

	write(6,*) '========================================'
	write(6,*) 'finished checking domain...'
	write(6,*) '========================================'

	end

!*******************************************************************

	subroutine make_elem_index(bghost,ic,icolor &
     &				,nk,ne,nenv,nodep,elemp)

! makes element index for domain with color ic

	use basin

	implicit none

	logical, intent(in) ::  bghost	!include ghost nodes (and elements)
	integer, intent(in) ::  ic
	integer, intent(in) ::  icolor(nkn)
	integer, intent(out) ::  nk,ne
	integer, intent(out) ::  nenv(3,nel)
	integer, intent(out) ::  nodep(nkn)
	integer, intent(out) ::  elemp(nel)

	integer ie,ii,k,nic
	integer node_aux(nkn)

	ne = 0
	node_aux = 0
	nodep = 0
	elemp = 0

	do ie=1,nel
	  nic = 0
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( icolor(k) == ic ) nic = nic + 1
	  end do
	  if( .not. bghost ) nic = nic - 2	!only good if nic == 3
	  if( nic > 0 ) then	!color found
	    ne = ne + 1
	    nenv(:,ne) = nen3v(:,ie)
	    elemp(ne) = ie
	  end if
	end do
	    
	do ie=1,ne
	  do ii=1,3
	    k = nenv(ii,ie)
	    node_aux(k) = 1		! these nodes are needed
	  end do
	end do

	nk = 0

	do k=1,nkn
	  if( node_aux(k) > 0 ) then
	    nk = nk + 1
	    node_aux(k) = nk
	    nodep(nk) = k
	  end if
	end do

	do ie=1,ne
	  do ii=1,3
	    k = nenv(ii,ie)
	    nenv(ii,ie) = node_aux(k)
	  end do
	end do

	end

!*******************************************************************

	subroutine check_elem_index(nk,ne,nenv,kerr,kerror)

! checks the element structure
!
! in order to get correct error feed back, the original arrays
! nen3v, ipv, ipev have to be altered
! they are saved at the beginning and then restored at the end

	use basin
	use mod_geom
	use mod_connect

	implicit none

	integer nk,ne
	integer nenv(3,ne)
	integer kerr
	integer kerror(kerr)

	integer k,ie,ngrm

	call connect_init(nk,ne,nenv,kerr)
	call connect_errors(kerr,kerror)
	call connect_release

	end

!*******************************************************************

	subroutine translate_color(nkn,nk,nodep,icolor,icol)

	implicit none

	integer nkn,nk
	integer nodep(nkn)
	integer icolor(nkn)
	integer icol(nk)

	integer k

	do k=1,nk
	  icol(k) = icolor(nodep(k))
	end do

	end

!*******************************************************************

	subroutine adjust_domain(ic,nkn,icolor,kerr)

	implicit none

	integer ic
	integer nkn
	integer icolor(nkn)
	integer kerr

	integer k,nc,icc,kint
	integer, allocatable :: count(:)

	integer ipext,ipint

	kint = ipint(kerr)

	nc = maxval(icolor)
	allocate(count(0:nc))
	count = 0

	do k=1,nkn
	  icc = icolor(k)
	  count(icc) = count(ic) + 1
	end do

	!write(6,*) 'colors for domain ',ic,icolor(kint)
	!do icc=1,nc
	!  write(6,*) icc,count(icc)
	!end do

	write(6,*) 're-coloring: ',kerr,kint,icolor(kint),ic

	icolor(kint) = ic

	end

!*******************************************************************

	subroutine count_elements(nkn,nel,nen3v,ic,icolor,netot,neint)

	implicit none

	integer nkn,nel
	integer nen3v(3,nel)
	integer ic
	integer icolor(nkn)
	integer netot,neint

	integer k,ii,ie,icc,icount

	netot = 0
	neint = 0

	do ie=1,nel
	  icount = 0
	  do ii=1,3
	    k = nen3v(ii,ie)
	    icc = icolor(k)
	    if( icc == ic ) icount = icount + 1
	  end do
	  if( icount > 0 )  netot = netot + 1
	  if( icount == 3 ) neint = neint + 1
	end do

	end

!*******************************************************************

	subroutine write_partition_grd(icolor,kext)  

	use basin
	use shympi

	implicit none

	integer icolor(nkn)
	integer kext			!external node for detail - 0 if none

	logical bdetail
	integer iecolor1(nel)
	integer iecolor2(nel)
	integer inext(nkn)
	integer ieext(nel)
	integer k,ie,ii,ic,n3c
	integer nc,ncmin,ncmax
	real window
	integer icc(3)
	character*80 grdfile
	character*80 text,tcolor

	if( .not. shympi_is_master() ) return

	bdetail = kext > 0
	window = 10.

!---------------------------------------------------------------
! prepare arrays
!---------------------------------------------------------------

	do k=1,nkn
	  inext(k) = ipv(k)
	end do
	do ie=1,nel
	  ieext(ie) = ipev(ie)
	end do

!---------------------------------------------------------------
! color elements
!---------------------------------------------------------------

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    icc(ii) = icolor(k)
	  end do
	  iecolor1(ie) = -1		! means no coloring
	  iecolor2(ie) = -1
	  if( all(icc==icc(1)) ) iecolor1(ie) = icc(1)
	  if( icc(1) == icc(2) ) iecolor2(ie) = icc(1)
	  if( icc(2) == icc(3) ) iecolor2(ie) = icc(2)
	  if( icc(3) == icc(1) ) iecolor2(ie) = icc(3)
	end do

	text = 'error domain partitioning - only proper elements colored'
	if( bdetail ) then
	  call grd_name_part1('part_error1',kext,grdfile)
	  write(6,*) 'writing to file ',trim(grdfile)
          call write_grd_file_with_detail(grdfile,text,nkn,nel,xgv,ygv &
     &                ,nen3v,inext,ieext,icolor,iecolor1,kext,window)
	else
	  grdfile = 'part_error1.grd'
	  write(6,*) 'writing to file ',trim(grdfile)
          call write_grd_file(grdfile,text,nkn,nel,xgv,ygv &
     &                ,nen3v,inext,ieext,icolor,iecolor1)
	end if

	text = 'error domain partitioning - all elements colored'
	if( bdetail ) then
	  call grd_name_part1('part_error2',kext,grdfile)
	  write(6,*) 'writing to file ',trim(grdfile)
          call write_grd_file_with_detail(grdfile,text,nkn,nel,xgv,ygv &
     &                ,nen3v,inext,ieext,icolor,iecolor2,kext,window)
	else
	  grdfile = 'part_error2.grd'
	  write(6,*) 'writing to file ',trim(grdfile)
          call write_grd_file(grdfile,text,nkn,nel,xgv,ygv &
     &                ,nen3v,inext,ieext,icolor,iecolor2)
	end if

!---------------------------------------------------------------
! write single domains
!---------------------------------------------------------------

	ncmin = minval(icolor)
	ncmax = maxval(icolor)

	do ic=ncmin,ncmax
	  iecolor1 = -1
	  text = 'domain error for color ' // trim(tcolor)
	  call grd_name_part2('domain_error',ic,kext,grdfile)
	  write(6,*) 'writing to file ',trim(grdfile)
	  do ie=1,nel
	    do ii=1,3
	      k = nen3v(ii,ie)
	      icc(ii) = icolor(k)
	    end do
	    n3c = count( icc == ic )
	    if( n3c == 3 ) iecolor1(ie) = 3
	    if( n3c == 2 ) iecolor1(ie) = 2
	    if( n3c == 1 ) iecolor1(ie) = 1
	  end do
	  if( bdetail ) then
            call write_grd_file_with_detail(grdfile,text,nkn,nel,xgv,ygv &
     &                ,nen3v,inext,ieext,icolor,iecolor1,kext,window)
	  else
            call write_grd_file(grdfile,text,nkn,nel,xgv,ygv &
     &                ,nen3v,inext,ieext,icolor,iecolor1)
	  end if
	end do

!---------------------------------------------------------------
! end of routine
!---------------------------------------------------------------

	end

!*******************************************************************

	subroutine grd_name_part1(text,k,name)

	implicit none

	character*(*) text
	integer k
	character*(*) name

	character*80 ttext

	write(ttext,'(i10)') k
	ttext = '_' // adjustl(ttext(1:10))
	name = trim(text) // trim(ttext) // '.grd'

	end

!*******************************************************************

	subroutine grd_name_part2(text,ic,k,name)

	implicit none

	character*(*) text
	integer ic,k
	character*(*) name

	character*80 ttext1,ttext2

	write(ttext1,'(i10)') ic
	ttext1 = '_' // adjustl(ttext1(1:10))
	write(ttext2,'(i10)') k
	ttext2 = '_' // adjustl(ttext2(1:10))
	name = trim(text) // trim(ttext1) // trim(ttext2) // '.grd'

	end

!*******************************************************************

	subroutine elim_isolated_element(icolor,kerr)

	use basin

	implicit none

	integer icolor(nkn)
	integer kerr

	end

!*******************************************************************

