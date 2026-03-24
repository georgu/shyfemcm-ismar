
!--------------------------------------------------------------------------
!
!    Copyright (C) 2003-2004,2007-2009,2012,2015,2015  Georg Umgiesser
!    Copyright (C) 2017-2019  Georg Umgiesser
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

!  description :
! 
!  main for mesh adjustment
! 
!  contents :
! 
!  shyadj			main
! 
!  revision log :
! 
!  01.01.2003	ggu	written
!  19.05.2003	ggu	some more checks
!  02.12.2004	ggu	some more debug checks (version 1.41)
!  02.12.2004	ggu	copyright statement added
!  20.03.2007	ggu	fake version
!  16.04.2008	ggu	new Makefile structure
!  09.12.2008	ggu	small changes to Makefile
!  26.01.2009	ggu	makedepend, compiler warnings
!  06.04.2009	ggu	new param.h, include basin.h, better error handling
!  24.04.2009	ggu	new call to rdgrd()
!  30.03.2012	ggu	new call to node_info() (was bug)
!  12.10.2015	ggu	changed VERS_7_3_3
!  25.05.2017	ggu	changed VERS_7_5_28
!  16.10.2018	ggu	changed VERS_7_5_50
!  18.12.2018	ggu	changed VERS_7_5_52
!  21.05.2019	ggu	changed VERS_7_5_62
!  14.02.2022	ggu	set all depth values to flag
!  23.02.2026   ggu     checks to avoid negative areas
!  27.02.2026   ggu     write more information to terminal
!  06.03.2026   ggu     completely restructured
! 
!  notes :
! 
!  to use static (non moveable) nodes please consult ls mod_adj_static
! 
! ***************************************************************

	program shyadj

!  adjusts elements after automatic mesh generator

	use mod_adj_grade
	use mod_depth
	use basin

	implicit none

	integer kspecial
	integer nlidim,nlndim
	integer k
	integer nc
	integer nco,nknh,nli
	integer nk,ne,nl,nne,nnl
	real, parameter :: flag = -999.
	logical bstop
	character*80 file

	real, allocatable :: dx(:),dy(:)
	integer, allocatable :: ic(:)
	integer, allocatable :: ianv(:)
	integer, allocatable :: ipaux(:)
	integer, allocatable :: iaux(:)
	real, allocatable :: raux(:)

! ------------------------------------------------------- end declaration

!---------------------------------------------------------------
! set some parameters
!---------------------------------------------------------------

	kspecial = 0		!set to 0 for no debug

	file = ' '
	bstop = .false.

!---------------------------------------------------------------
! parse command line
!---------------------------------------------------------------

        call shyadj_init(file)

!---------------------------------------------------------------
! read grid file with nodes and elements
!---------------------------------------------------------------

	call grd_set_write(.false.)
	call grd_read(file)
	call grd_get_params(nk,ne,nl,nne,nnl)

	if ( .not. bsilent ) write(6,'(a,5i7)') ' grid parameters: ',nk,ne,nl
	if( nk .le. 0 ) stop

	call grd_to_basin

	allocate(dx(nkn),dy(nkn))
	allocate(ic(nkn),ianv(nkn))
	allocate(ipaux(nkn),iaux(nel))
	allocate(raux(nel))

!---------------------------------------------------------------
! save depth information in elements to nodes
!---------------------------------------------------------------

	call mod_depth_init(nkn,nel)
	call grd_get_depth(nk,ne,hkv,hev)

	nknh = 0
	do k=1,nkn
	  if( hkv(k) .ne. flag ) nknh = nknh + 1
	end do

	if( nknh .eq. 0 ) then
	  write(6,*) 'copying element depth to nodes...'
	  call hev2hkv
	end if

!---------------------------------------------------------------
! determine grade and boundary nodes
!---------------------------------------------------------------

	call maxgrd(nkn,nel,nen3v,ngr)	  !determines ngr (on bnd 1 too low)
	if( .not. bquiet ) write(6,*) 'maximum grade: ',ngr
	call mod_adj_grade_init(nkn,ngr)  !allocates global arrays
	call setgrd(nkn,nel,nen3v,ngrade) !sets ngrade (still wrong on bnd)

	call stats('first call')

! make boundary nodes (flag nbound)

	call mkbound(nkn,nel,ngrdi,nen3v,ngrade,nbound,ngri)
        call mkstatic(nkn,ianv,nbound)
        if( bverbose ) write(6,*) 'grading nodes done'
	call stats('boundary nodes')

	call node_info(kspecial)

!---------------------------------------------------------------
! initialize plot
!---------------------------------------------------------------

	if( bplot ) then
	  call qopen
	  call plobas
	end if

        !call smooth_grid(nsmooth,asmooth)	!only for debug

!---------------------------------------------------------------
! first cycle
!---------------------------------------------------------------

!  eliminate 4- grades

	if( .not. bsilent ) then
          write(6,*) '================================='
          write(6,*) 'first cycle...'
          write(6,*) '================================='
	end if

	call chkgrd('first cycle - checking before low')
	call elimlow
	call chkgrd('first cycle - checking after low')

	if( bplot ) call plobas
	call stats('first cycle - 4- grades')
	call node_info(kspecial)

!  eliminate 8+ grades

	call chkgrd('first cycle - checking before 8+')
	call elimhigh(8)
	call chkgrd('first cycle - checking after 8+')

	if( bplot ) call plobas
	call stats('first cycle - 8+ grades')
	call node_info(kspecial)

!  eliminate 7+ grades

	call chkgrd('first cycle - checking before 7+')
	call elimhigh(7)
	call chkgrd('first cycle - checking after 7+')

	if( bplot ) call plobas
	call stats('first cycle - 7+ grades')
	call node_info(kspecial)

!  smoothing

	!call write_grid('new_nosmooth.grd')

	call chkgrd('first cycle - checking before smoothing')
        call smooth_grid(nsmooth,asmooth)
	call chkgrd('first cycle - checking after sdmoothing')

	if( bplot ) call plobas
	call node_info(kspecial)

	!call write_grid('new_smooth1.grd')

!---------------------------------------------------------------
! second cycle
!---------------------------------------------------------------

	if( .not. bsilent ) then
          write(6,*) '================================='
          write(6,*) 'second cycle...'
          write(6,*) '================================='
	end if

	call chkgrd('second cycle - before low/high')
        call elimlow
	call elimhigh(8)
	call elimhigh(7)
	call chkgrd('second cycle - after low/high')

	if( bplot ) call plobas
	call stats('second cycle - after elimhigh')
	call node_info(kspecial)

!  eliminate 5-5 grades

	call chkgrd('second cycle - checking before 5-5')
	call elim_5_5
	call chkgrd('second cycle - checking after 5-5')

	if( bplot ) call plobas
	call stats('second cycle - 5-5 grades')
	call node_info(kspecial)

!  eliminate 5-7-5 grades

	call chkgrd('second cycle - checking before 5-7-5')
	call elim_5_7_5
	call chkgrd('second cycle - checking after 5-7-5')

	if( bplot ) call plobas
	call stats('second cycle - 5-7-5 grades')
	call node_info(kspecial)

!---------------------------------------------------------------
! third cycle
!---------------------------------------------------------------

	if( .not. bsilent ) then
          write(6,*) '================================='
          write(6,*) 'third cycle...'
          write(6,*) '================================='
	end if

	call chkgrd('third cycle - checking before high, 5-5, 5-7-5')
	call elimhigh(8)
        call elimhigh(7)
        call elim_5_5
        call elim_5_7_5
	call chkgrd('third cycle - checking after high, 5-5, 5-7-5')

	if( bplot ) call plobas
        call stats('thrid cycle - end')
	call node_info(kspecial)

!---------------------------------------------------------------
! final smoothing
!---------------------------------------------------------------

	call write_grid('final_before_smoothing.grd')

	if( .not. bsilent ) then
          write(6,*) '================================='
          write(6,*) 'final smoothing...'
          write(6,*) '================================='
	end if

	call chkgrd('final - checking before smoothing')
        call smooth_grid(nsmooth,asmooth)
	call chkgrd('final - checking after smoothing')

	if( bplot ) call plobas
	call node_info(kspecial)

!---------------------------------------------------------------
! write to grd file
!---------------------------------------------------------------

	if( .not. bsilent ) then
          write(6,*) '================================='
          write(6,*) 'writing to grid...'
          write(6,*) '================================='
	end if

	!call plot_nodes_with_grade(4)

	call chkgrd('final check')
        call stats('final solution')
	call node_info(kspecial)
	call show_strange_grades
	hev = flag
	hkv = flag
	hm3v = flag
	call write_grid('adjust_final.grd')

	call qclose	!this is safe to call

	if( .not. bsilent ) then
	  write(6,*) 'Successful completion.'// &
     &			' Output has been written to adjust_final.grd'
	  if( bplot ) write(6,*) 'plot of basin has been written to plot.ps'
	end if

!---------------------------------------------------------------
! end of routine
!---------------------------------------------------------------

	stop
   97	continue
	write(6,*) 'error reading grd file'
	stop 'error stop rdgrd'
   99	continue
	write(6,*) 'error external to internal numbering'
	stop 'error stop ex2in'
	end

! ***********************************************************
! ***********************************************************
! ***********************************************************

	subroutine hev2hkv

!  saves information about depth to nodes

	use mod_depth
	use basin

	implicit none

	integer k,ie,ii
	integer ic(nkn)

	do k=1,nkn
	  hkv(k) = 0.
	  ic(k) = 0
	end do

	do ie=1,nel
	  do ii=1,3
	    k = nen3v(ii,ie)
	    hkv(k) = hkv(k) + hev(ie)
	    ic(k) = ic(k) + 1
	  end do
	end do

	do k=1,nkn
	  if( ic(k) .gt. 0 ) then
	    hkv(k) = hkv(k) / ic(k)
	  end if
	end do

	end

! ***********************************************************

	subroutine stats(text)

	use mod_adj_grade
	use basin

	implicit none

	character*(*) text

	call statgrd(iugrade,text,nkn,ngr,ngrade,nbound)

	end

! ***********************************************************

	subroutine show_strange_grades

	use mod_adj_grade
	use basin

	implicit none

	integer k,kext,n
	integer itot,ntot,nmax
	integer igr(ngrdi)

	igr = 0
        do k=1,nkn
          if( nbound(k) .ne. 0 ) cycle
          n = ngrade(k)
	  igr(n) = igr(n) + 1
        end do

	ntot = igr(3) + igr(4)
	do n=ngrdi,8,-1
	  if( ntot > 50 ) exit
	  ntot = ntot + igr(n)
	end do
	nmax = n

	if( bverbose ) then
	  write(6,*) 'Listing nodes with not fixable grades: ',ntot
	  write(6,*) 'only nodes are shown with grades < 5 and > ',nmax
	end if

	itot = 0
        do k=1,nkn
          if( nbound(k) .ne. 0 ) cycle
          n = ngrade(k)
          if( n < 5 .or. n > 7 ) itot = itot + 1
          if( n < 5 .or. n > nmax ) then
	    call nint2ext(k,kext)
	    if( bverbose ) then
	      write(6,*) 'node ',k,' (extern ',kext,') with grade ',n
	    end if
          end if
        end do

	if( .not. bquiet ) then
	  write(6,*) 'there are nodes with non fixable grades: ',itot
	end if

	end

! ***********************************************************

        subroutine shyadj_init(grdfile)

        use clo
	use mod_adj_grade

        implicit none

        character*(*) grdfile

	integer n
	real f(2)
	character*80 sline

	integer iscanf

        call clo_init('shyadj','grd-file','2.0')

        call clo_add_info('regolarize grd file')

        call clo_add_option('verbose',.false.,'be verbose')
        call clo_add_option('quiet',.false.,'be quiet')
        call clo_add_option('silent',.false.,'be silent')

        call clo_add_option('smooth params',' ','smoothing options')
        call clo_add_option('plot',.false.,'create plot of grades')
        call clo_add_option('check',.false.,'checks consistency of grid')
        call clo_add_option('check_all',.false.,'more checks on consistency')

        call clo_add_sep('additional information')
        call clo_add_com('  params is nsmooth[,asmooth]')
	call clo_add_com('    nsmooth is number of smoothing iterations')
        call clo_add_com('    asmooth is smoothing strength')
        call clo_add_com('    defaults: 50,0.01')

        call clo_parse_options

        call clo_get_option('verbose',bverbose)
        call clo_get_option('quiet',bquiet)
        call clo_get_option('silent',bsilent)
        call clo_get_option('smooth',sline)
        call clo_get_option('plot',bplot)
        call clo_get_option('check',bcheck)
        call clo_get_option('check_all',bcheck_all)

	if( bsilent ) bquiet = .true.
	if( bquiet ) bverbose = .false.

	if( bquiet ) iugrade = 66

	if( bcheck_all ) bcheck = .true.

        call clo_check_files(1)
        call clo_get_file(1,grdfile)

	n = iscanf(sline,f,2)
	if( n < 0 .or. n > 2 ) then
	  write(6,*) 'error in smoothing parameters: ',trim(sline)
	else if( n == 0 ) then
	  !use default
	else if( n == 1 ) then
	  nsmooth = nint(f(1))
	else
	  nsmooth = nint(f(1))
	  asmooth = f(2)
	end if

        call shyfem_set_short_copyright(bquiet)
        if( .not. bsilent ) then
	  call shyfem_copyright('shyadj - regularize finite element grids')
        end if

	if( .not. bsilent ) then
	  write(6,*) 'using smoothing parameters: ',nsmooth,asmooth
	end if

        end

! ***********************************************************

