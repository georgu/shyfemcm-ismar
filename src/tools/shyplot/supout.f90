
!--------------------------------------------------------------------------
!
!    Copyright (C) 2003-2004,2007-2019  Georg Umgiesser
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

!  routines for reading data files
! 
!  revision log :
! 
!  31.10.2003	ggu	new routines velclose(), resetsim()
!  31.10.2003	ggu	new routines for handling wind
!  22.09.2004	ggu	new routines for handling ous file
!  22.09.2004	ggu	use OUS file for 3D, try to use alsways 3D (level=0)
!  05.10.2004	ggu	adjustments -> use always 3D data structure
!  14.03.2007	ggu	new routines for wave plotting
!  17.09.2008	ggu	new routine level_e2k to compute ilhkv from ilhv
!  09.10.2009	ggu	read also pressure from wind file
!  13.10.2009	ggu	set nlv once file is read
!  23.02.2010	ggu	change in reading wind file
!  23.03.2010	ggu	changed v6.1.1
!  20.12.2010	ggu	changed VERS_6_1_16
!  30.03.2011	ggu	new routines to handle fvl files (not yet integrated)
!  31.03.2011	ggu	use fvl routines to exclude areas from plot
!  14.04.2011	ggu	changed VERS_6_1_22
!  12.07.2011	ggu	in prepsim use what is available for dry areas
!  15.07.2011	ggu	changed VERS_6_1_28
!  18.08.2011	ggu	bug fix in nosopen() -> extra comma eliminated 
!  31.08.2011	ggu	new routines for handling EOS files
!  01.09.2011	ggu	changed VERS_6_1_32
!  14.11.2011	ggu	call to init_sigma_info() to setup layer info
!  22.11.2011	ggu	changed VERS_6_1_37
!  09.12.2011	ggu	changed VERS_6_1_38
!  19.12.2011	ggu	new routine level_k2e -> called for nos files
!  24.01.2012	ggu	changed VERS_6_1_41
!  21.06.2012	ggu	changed VERS_6_1_54
!  13.06.2013	ggu	new routines for handling FEM files
!  19.06.2013	ggu	changed VERS_6_1_66
!  03.09.2013	ggu	level_k2e -> level_k2e_sh, level_e2k -> level_e2k_sh
!  12.09.2013	ggu	changed VERS_6_1_67
!  28.01.2014	ggu	changed VERS_6_1_71
!  05.03.2014	ggu	new read for ous and nos files (use date)
!  30.05.2014	ggu	changed VERS_6_1_76
!  07.07.2014	ggu	changed VERS_6_1_79
!  18.07.2014	ggu	changed VERS_7_0_1
!  20.10.2014	ggu	deleted is2d() and out reading routines
!  30.10.2014	ggu	changed VERS_7_0_4
!  26.11.2014	ggu	changed VERS_7_0_7
!  05.12.2014	ggu	changed VERS_7_0_8
!  23.12.2014	ggu	changed VERS_7_0_11
!  19.01.2015	ggu	changed VERS_7_1_2
!  19.01.2015	ggu	changed VERS_7_1_3
!  10.02.2015	ggu	use different file units (more than one can be opened)
!  26.02.2015	ggu	changed VERS_7_1_5
!  01.04.2015	ggu	changed VERS_7_1_7
!  10.07.2015	ggu	changed VERS_7_1_50
!  13.07.2015	ggu	changed VERS_7_1_51
!  17.07.2015	ggu	changed VERS_7_1_52
!  17.07.2015	ggu	changed VERS_7_1_80
!  20.07.2015	ggu	changed VERS_7_1_81
!  14.09.2015	ggu	introduced bwind and bvel (for velocities)
!  05.11.2015	ggu	changed VERS_7_3_12
!  19.02.2016	ggu	changed VERS_7_5_2
!  28.04.2016	ggu	changed VERS_7_5_9
!  25.05.2016	ggu	changed VERS_7_5_10
!  27.06.2016	ggu	changed VERS_7_5_16
!  12.01.2017	ggu	changed VERS_7_5_21
!  13.04.2017	ggu	changed VERS_7_5_25
!  25.05.2017	ggu	changed VERS_7_5_28
!  11.07.2017	ggu	changed VERS_7_5_30
!  14.11.2017	ggu	changed VERS_7_5_36
!  03.04.2018	ggu	changed VERS_7_5_43
!  18.04.2018	ggu	set up bkplot (node to be plotted)
!  16.10.2018	ggu	changed VERS_7_5_50
!  25.10.2018	ggu	changed VERS_7_5_51
!  18.12.2018	ggu	changed VERS_7_5_52
!  21.05.2019	ggu	changed VERS_7_5_62
!  21.10.2023	ggu	many unsued routines deleted
!  02.04.2025	ggu	adjust adjust_no_plot_area() to only plot in some areas
! 
! **********************************************************
! **********************************************************
! **********************************************************

	module supout

	logical bdebug_out
	parameter (bdebug_out = .false.)

	integer, save :: nunit

	integer, save :: iformat,iwave

	integer, save :: nunit_wave,nunit_ous,nunit_nos,nunit_fvl &
     &                  ,nunit_eos,nunit_fem

	real, save :: regp(7)

	end module supout

! **********************************************************

        subroutine velopen

!  opens velocity file (2d and 3d)

        implicit none

        !call ousopen
	stop 'error stop velopen: no more OUS files'

        end

! **********************************************************

        function velnext(it)

!  gets next velocity field (2d and 3d)

        implicit none

        logical velnext
        integer it

        logical outnext, ousnext

        !velnext = ousnext(it)
	velnext = .false.

        end

! **********************************************************

        subroutine velclose

!  closes velocity file

        implicit none

        !call ousclose

        end

! ******************************************************
! ******************************************************
! ******************************************************

        function ous_is_available()
        implicit none
        logical ous_is_available
	ous_is_available = .false.
        end

        subroutine ousclose
        implicit none
        end

        subroutine ousinfo(nvers,nkn,nel,nlv)
        implicit none
        integer nvers,nkn,nel,nlv
        end

        subroutine ousopen
	stop 'error stop velopen: no more OUS files'
	end

        function ousnext(it)
	implicit none
        logical ousnext         !true if record read, flase if EOF
        integer it              !time of record
	ousnext = .false.
	end

! ******************************************************
! ******************************************************
! ******************************************************

        subroutine reset_dry_mask

!  resets mask data structure

	use mod_hydro_plot

	implicit none

        bwater = .true.
        bkwater = .true.
        bplot = .true.
        bkplot = .true.

        end

! ******************************************************

        subroutine adjust_no_plot_area

        use basin
	use mod_hydro_plot

        implicit none

	logical bset
        integer ie,ii,k
        integer ianopl,iadopl,ia
	real getpar

        ianopl = nint(getpar('ianopl'))		!do not plot here
        iadopl = nint(getpar('iadopl'))		!only plot here

	if( ianopl < 0 .and. iadopl < 0 ) return
	if( ianopl >= 0 .and. iadopl >= 0 ) then
	  write(6,*) 'ianopl,iadopl: ',ianopl,iadopl
	  write(6,*) 'only one of ianopl,iadopl can be given'
	  stop 'error stop adjust_no_plot_area: ianopl,iadopl'
	end if

	if( iadopl >= 0 ) then
	  ia = iadopl
	  bset = .true.
	  bwater = .false.
	  bplot = .false.
	  bkplot = .false.
	end if
	if( ianopl >= 0 ) then
	  ia = ianopl
	  bset = .false.
	  bwater = .true.
	  bplot = .true.
	  bkplot = .true.
	end if

        do ie=1,nel
          if( iarv(ie) == ia ) then
	    bwater(ie) = bset
	    bplot(ie) = bset
	    do ii=1,3
	      k = nen3v(ii,ie)
	      bkplot(k) = bset
	    end do
	  end if
        end do

        end

! ******************************************************

	subroutine prepare_dry_mask

!  prepares simulation for use - computes wet and dry areas

	use mod_hydro_plot
	use mod_hydro
	use levels

	implicit none

	logical bshowdry
	integer level
	real href,hzmin,hdry

	logical fvl_is_available
	logical ous_is_available
	integer getlev
	real getpar

! ---------------------------------------------------
!  set up mask of water points
! ---------------------------------------------------

!  set bshowdry = .true.  if you want to show dried out areas
!  set bshowdry = .false. if you want to plot all areas

	hdry = 0.05				!level for drying
	bshowdry = .true.			!show dry areas, else plot all
	!bshowdry = .false.			!show dry areas, else plot all

	if( .not. bshowdry ) hdry = -1.e+30	!no drying

        href = getpar('href')
        hzmin = getpar('hzmin')
	level = getlev()

	call reset_dry_mask
	write(6,*) 'calling prepare_dry_mask...'

	!if( ous_is_available() ) then			!...handle on elements
	if( .false. ) then			!...handle on elements

	  write(6,*) 'using zeta for dry areas'
	  if( bshowdry ) then
            call set_dry_mask(bwater,znv,zenv,href,hzmin) !false if znv/=zenv
	  end if
          call set_level_mask(bwater,ilhv,level)	!element has this level
	  call make_dry_node_mask(bwater,bkwater)	!copy elem to node mask

	!else if( fvl_is_available() ) then		!...handle on nodes
	else if( .false. ) then		!...handle on nodes

	  write(6,*) 'using fvl file for dry areas: ',hdry
	  call set_dry_volume_mask(bkwater,hdry)	!guess if dry using vol
	  call make_dry_elem_mask(bwater,bkwater)	!copy node to elem mask
          call set_level_mask(bwater,ilhv,level)	!element has this level
	  call make_dry_node_mask(bwater,bkwater)	!copy elem to node mask

	else

	  write(6,*) 'no information on dry areas: ',hdry
          call set_level_mask(bwater,ilhv,level)	!element has this level
	  call make_dry_node_mask(bwater,bkwater)	!copy elem to node mask

	end if

        call adjust_no_plot_area
	call make_dry_node_mask(bwater,bkwater)		!copy elem to node mask
        call info_dry_mask(bwater,bkwater)

! ---------------------------------------------------
!  end of routine
! ---------------------------------------------------

	end

! ******************************************************

	subroutine set_dry_volume_mask(bkw,hdry)

	use mod_hydro_plot
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	logical bkw(nkn)
	real hdry

	integer k,idry
	real vol,area,h
	real hhmin,hhmax

	hhmin = 1.e+30
	hhmax = -1.e+30
	idry = 0

	do k=1,nkn
	  vol = fvlv(1,k)
	  area = arfvlv(k)
	  h = vol/area
	  hhmin = min(hhmin,h)
	  hhmax = max(hhmax,h)
	  bkw(k) = h .ge. hdry
	  if( .not. bkw(k) ) idry = idry + 1
	end do

	write(6,*) 'min/max depth: ',hhmin,hhmax,hdry,idry

	end

! ******************************************************
! ******************************************************
! ******************************************************

        subroutine waveini
        implicit none
        end

        subroutine waveclose
        implicit none
        end

        function wavenext(it)
        implicit none
	logical wavenext
	integer it
	wavenext = .false.
        end

	subroutine waveopen
	implicit none
	stop 'error stop waveopen: no more NOS files'
        end

! ******************************************************
! ******************************************************
! ******************************************************

        subroutine nosini
        implicit none
        end

        subroutine nosclose
        implicit none
        end

        function nosnext(it,ivaria,nlvdi,p3)
	use basin
        implicit none
	logical nosnext
	integer it,ivaria,nlvdi
	real p3(nlvdi,nel)
	nosnext = .false.
        end

	subroutine nosopen(type)
	implicit none
	character*(*) type
	stop 'error stop nosopen: no more NOS files'
        end

! ******************************************************
! ******************************************************
! ******************************************************

	subroutine polar2xy00(n,speed,dir,uv,vv)

	implicit none

	integer n
	real speed(n), dir(n)
	real uv(n), vv(n)

	integer i
	real a
	real, parameter :: rad = atan(1.) / 45.

	do i=1,n
	  a = dir(i)
          a = 90. - a + 180.
          do while( a .lt. 0. )
            a = a + 360.
          end do
          a = mod(a,360.)

	  uv(i) = speed(i) * cos( rad * a )
	  vv(i) = speed(i) * sin( rad * a )
	end do

	end

! ******************************************************
! ******************************************************
! ******************************************************


! ******************************************************
! ******************************************************
! ******************************************************
!  routines to read fvl file
! ******************************************************
! ******************************************************
! ******************************************************

	subroutine fvlini

	end

! ******************************************************

	function fvl_is_available()

!  checks if FVL file is opened

	implicit none

	logical fvl_is_available

	fvl_is_available = .false.

	end

! ******************************************************

	subroutine fvlclose

!  closes FVL file

	implicit none

	end

! ******************************************************

	subroutine fvlopen(type)

!  opens FVL file and reads header

	implicit none

	character*(*) type

	stop 'error stop velopen: no more OUS files'

	end

! ******************************************************

	subroutine fvlnext(it,nlvddi,array)

!  reads next FVL record

	use levels

	implicit none

	integer it		!time of record
	integer nlvddi		!dimension of vertical coordinate
	real array(nlvddi,*)	!values for variable

	end

! ******************************************************
! ******************************************************
! ******************************************************
!  element values
! ******************************************************
! ******************************************************
! ******************************************************

	subroutine eosini

!  initializes internal data structure for EOS file

	implicit none

	end

! ******************************************************

	subroutine eosclose

!  closes EOS file

	implicit none

	end

! ******************************************************

	subroutine eosopen(type)

!  opens EOS file and reads header

	use mod_depth
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	character*(*) type

	stop 'error stop velopen: no more EOS files'

	end

! ******************************************************

	function eosnext(it,ivar,nlvddi,array)

!  reads next EOS record - is true if a record has been read, false if EOF

	use levels

	implicit none

	logical eosnext		!true if record read, flase if EOF
	integer it		!time of record
	integer ivar		!type of variable
	integer nlvddi		!dimension of vertical coordinate
	real array(nlvddi,*)	!values for variable

	eosnext = .false.

	end

! ******************************************************
! ******************************************************
! ******************************************************
!  aux routines
! ******************************************************
! ******************************************************
! ******************************************************

	subroutine level_e2k_sh

!  computes max level at nodes from elements

	use levels
	use basin

	implicit none

	call level_e2k(nkn,nel,nen3v,ilhv,ilhkv)

	end

! ******************************************************

	subroutine level_k2e_sh

!  computes level at elems from nodes (not exact)

	use levels
	use basin

	implicit none

	call level_k2e(nkn,nel,nen3v,ilhkv,ilhv)

	end

! ******************************************************

	subroutine array_check(n,a1,a2,text)

!  checks if arrays are equal

	implicit none

	integer n
	real a1(n)
	real a2(n)
	character*(*) text

	integer i

	do i=1,n
	  if( a1(i) .ne. a2(i) ) then
	    write(6,*) text,' : arrays differ'
	    write(6,*) n,i,a1(i),a2(i)
	    stop 'error stop array_check: arrays differ'
	  end if
	end do

	end

! ******************************************************

	subroutine array_i_check(n,a1,a2,text)

!  checks if arrays are equal

	implicit none

	integer n
	integer a1(n)
	integer a2(n)
	character*(*) text

	integer i

	do i=1,n
	  if( a1(i) .ne. a2(i) ) then
	    write(6,*) text,' : arrays differ'
	    write(6,*) n,i,a1(i),a2(i)
	    stop 'error stop array_check: arrays differ'
	  end if
	end do

	end

! ******************************************************
! ******************************************************
! ******************************************************
!  fem files
! ******************************************************
! ******************************************************
! ******************************************************

	subroutine femini

!  initializes internal data structure for FEM files

	use supout

	implicit none

	integer icall
	save icall
	data icall /0/

	if( bdebug_out ) then
	  write(6,*) 'debug_out: femini'
	  write(6,*) icall,nunit_fvl
	end if

	if( icall .ne. 0 ) return

	icall = 1

	nunit_fem = 0
	iformat = 0

	end

! ******************************************************

	subroutine femclose

!  closes FEM file

	use supout

	implicit none

	call femini
	if( nunit_fem .gt. 0 ) close(nunit_fem)
	nunit_fem = 0

	end

! ******************************************************

	subroutine femopen(type)

!  opens FEM file and reads header

	use mod_depth
	use levels
	use basin, only : nkn,nel,ngr,mbw
	use simul
	use supout

	implicit none

	character*(*) type

	character*80 file

	logical bformat,bregdata
	integer nvers,np,it,lmax,ntype
	integer nknaux,nelaux,nlvaux,nvar
	integer nknddi,nelddi,nlvddi
	integer ierr,l
	integer datetime(2)
	real regpar(7)
	double precision dtime
	integer ifemop,fem_file_regular

!  initialize routines

	call femini

!  open file

	np = nkn
	call def_make(type,file)
	call fem_file_read_open(file,np,iformat,nunit)

	if( nunit .le. 0 ) then
                write(6,*) file
		stop 'error stop femopen: cannot open FEM file'
        else
                write(6,*) 'File opened :'
                inquire(nunit,name=file)
                write(6,*) file
                write(6,*) 'Reading FEM file ...'
	end if

!  read first header

        call fem_file_read_params(iformat,nunit,dtime &
     &                          ,nvers,np,lmax,nvar,ntype,datetime,ierr)

	if( ierr .ne. 0 ) then
		write(6,*) 'ierr = ',ierr
		stop 'error stop femopen: error reading header'
	end if

	bregdata = fem_file_regular(ntype) > 0

        write(6,*)
        write(6,*) ' nvers = ', nvers
        write(6,*) '   nkn = ',np,   ' ntype = ',ntype
        write(6,*) '   nlv = ',lmax, '  nvar = ',nvar
        write(6,*)
	if( bregdata ) then
          write(6,*) 'the file is a regular grid'
          write(6,*)
	end if

	if( .not. bregdata .and. nkn .ne. np ) goto 99
	!if( nkn .lt. np ) goto 99
	if( bregdata .and. lmax > 1 ) goto 98

        nlv = lmax
	call allocate_simulation(nel)
	!if( bregdata ) call reallocate_2d_arrays(np) !re-allocate with minimum np
	call get_dimension_post(nknddi,nelddi,nlvddi)
	if( nlvddi .lt. lmax ) goto 99

	it = nint(dtime)

!  read second header

	call fem_file_read_2header(iformat,nunit,ntype,nlv &
     &				,hlv,regpar,ierr)

	if( ierr .ne. 0 ) then
		write(6,*) 'ierr = ',ierr
		stop 'error stop femopen: error reading hlv'
	end if

	call level_k2e_sh
	call init_sigma_info(nlv,hlv)		!sets up hlv

	write(6,*) 'hlv: ',nlv
	write(6,*) (hlv(l),l=1,nlv)

!  rewind for a clean state

	nunit_fem = nunit
	rewind(nunit)

	return
   98	continue
	write(6,*) 'error in parameters : regular - lmax'
	write(6,*) 'the file is regular and 3D'
	write(6,*) 'Cannot handle interpolation yet'
	stop 'error stop femopen'
   99	continue
	write(6,*) 'error in parameters : basin - simulation'
	write(6,*) 'nkn : ',nkn,np
	write(6,*) 'nlv : ',nlvdi,lmax
	write(6,*) 'parameters are different between basin and simulation'
	stop 'error stop femopen'
	end

! ******************************************************

	function femnext(atime,ivar,nlvddi,nkn,array)

!  reads next FEM record - is true if a record has been read, false if EOF

	use levels
	use supout

	implicit none

	logical femnext			!true if record read, flase if EOF
	integer it			!time of record
	double precision atime		!absolute time
	integer ivar			!type of variable
	integer nlvddi			!dimension of vertical coordinate
	integer nkn			!number of points needed
	real array(nlvddi,nkn,1)	!values for variable

	logical bfound,bformat,bregdata,bwind,bvel
	integer ierr
	integer i,iv,ip
	integer nvers,np,lmax,nvar,ntype
	integer datetime(2)
	real regpar(7)
	double precision dtime
	real, allocatable :: p3read(:,:,:)
	real, allocatable :: v1v(:)
	real fact
	real vmin,vmax
	character*80 string

	integer fem_file_regular

	if( nlvddi .ne. nlvdi ) stop 'error stop femnext: nlvddi'

	call femini
	nunit = nunit_fem
	bformat = iformat .eq. 1

	!write(6,*) 'femnext format: ',bformat,iformat

	np = 0
        call fem_file_read_params(iformat,nunit,dtime &
     &                          ,nvers,np,lmax,nvar,ntype,datetime,ierr)
	if( ierr .ne. 0 ) goto 7
	nlv = lmax
	regpar = 0
	call fem_file_read_2header(iformat,nunit,ntype,nlv &
     &				,hlv,regpar,ierr)
	if( ierr .ne. 0 ) goto 7

	bregdata = fem_file_regular(ntype) > 0
	if( .not. bregdata .and. np .ne. nkn ) goto 99
	if( bregdata .and. lmax > 1 ) goto 98

	write(6,*) 'ggggggguuuuu bregdata: ',bregdata,np

	if( bregdata ) then
	  write(6,*) 'plotting regular grid...'
	  !write(6,*) 'not yet ready for regular grid...'
	  !stop
	  allocate(p3read(nlvddi,np,2),v1v(np))
	else
	  allocate(p3read(nlvddi,nkn,2),v1v(nkn))
	end if
	write(6,*) 'p3read allocated: ',nlvddi,nkn,np,2
	p3read = 0.

	!it = nint(dtime)
	call ptime_set_date_time(datetime(1),datetime(2))
	call ptime_set_dtime(dtime)
	call ptime_get_atime(atime)

	ip = 1
	bfound = .false.
	do i=1,nvar
	  if( bfound ) then
            call fem_file_skip_data(iformat,nunit &
     &                          ,nvers,np,lmax &
     &                          ,string,ierr)
	    if( ierr .ne. 0 ) goto 98
	  else
	    write(6,*) 'reading data: ',bregdata,nlvddi,nkn,np,ip
            call fem_file_read_data(iformat,nunit &
     &                          ,nvers,np,lmax &
     &				,string &
     &                          ,ilhkv,v1v &
     &                          ,nlvddi,p3read(1,1,ip),ierr)
	    if( ierr .ne. 0 ) goto 98
	    call string2ivar(string,iv)
	    bfound = iv .eq. ivar
	    bwind = iv .eq. 21
	    bvel = iv .eq. 2
	    if( bfound .and. ( bwind .or. bvel ) ) then
	      ip = ip + 1
	      if( ip .eq. 2 ) bfound = .false.
	    end if
	  end if
	end do

	if( ip == 3 ) then
	  v1v = sqrt( p3read(1,:,1)**2 + p3read(1,:,2)**2 )
	  vmin = minval(v1v)
	  vmax = maxval(v1v)
	  write(6,*) 'speed: ',vmin,vmax
	end if

	if( bregdata ) then
	  write(6,*) 'interpolating from regular grid... ',ip
	  ip = min(2,ip)
	  call fem_interpolate(nlvddi,nkn,np,ip,regpar,ilhkv &
     &					,p3read,array)
	  regp = regpar		!save for later
	else
	  array = p3read
	end if

	deallocate(p3read,v1v)

	if( bfound ) then
	  call level_k2e_sh
	else
	  ivar = 0
	end if

!  set return value

    7	continue

	if( ierr .gt. 0 ) then
		!stop 'error stop femnext: error reading data record'
		write(6,*) '*** femnext: error reading data record'
		femnext = .false.
	else if( ierr .lt. 0 ) then
		femnext = .false.
	else
		femnext = .true.
	end if

!  end

	return
   98	continue
	write(6,*) 'error in parameters : regular - lmax'
	write(6,*) 'the file is regular and 3D'
	write(6,*) 'Cannot handle interpolation yet'
	stop 'error stop femopen'
   99	continue
	write(6,*) 'nkn,np: ',nkn,np
	stop 'error stop femnext: np different from nkn'
	end

! ******************************************************

	subroutine femscale(nlvddi,nkn,fact,array)

	implicit none

	integer nlvddi,nkn
	real fact
	real array(nlvddi,nkn)

	integer k,l

	do k=1,nkn
	  do l=1,nlvddi
	    array(l,k) = fact * array(l,k)
	  end do
	end do

	end

! ******************************************************

	subroutine fem_interpolate(nlvddi,nkn,np,ip,regpar,ilhkv &
     &			,p3reg,array)

!  interpolates from a regular grid (only for 2D)

	implicit none

	integer nlvddi,nkn,np,ip
	real regpar(7)
	integer ilhkv(nkn)
	real p3reg(nlvddi,np,ip)
	real array(nlvddi,nkn,ip)

	integer i,ivar,j
	integer nx,ny
	real x0,y0,dx,dy,flag
	real areg(np)
	real afem(nkn)

	nx = nint(regpar(1))
	ny = nint(regpar(2))
	x0 = regpar(3)
	y0 = regpar(4)
	dx = regpar(5)
	dy = regpar(6)
	flag = regpar(7)
	ilhkv = 1		!works only for 2D

	call setgeo(x0,y0,dx,dy,flag)
	
	do ivar=1,ip
	  do i=1,np
	    areg(i) = p3reg(1,i,ivar)
	  end do
	  !write(6,*) 'reg: ',(areg(j),j=1,np,np/10)

	  call am2av(areg,afem,nx,ny)

	  do i=1,nkn
	    array(1,i,ivar) = afem(i)
	  end do
	  !write(6,*) 'fem: ',(afem(j),j=1,nkn,nkn/10)
	end do

	end

! *****************************************************************
! *****************************************************************
! *****************************************************************

	subroutine allocate_simulation(npd)

	use mod_hydro_plot
	use mod_hydro_print
	use mod_hydro_vel
	use mod_hydro
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	integer npd

	integer np
	real flag

	nlvdi = nlv
	np = max(2*nel,npd)

	call levels_init(nkn,nel,nlvdi)
	call mod_hydro_init(nkn,nel,nlvdi)
	call mod_hydro_vel_init(nkn,nel,nlvdi)
	call mod_hydro_print_init(nkn,nlvdi)
	call mod_hydro_plot_init(nkn,nel,nlvdi,np)

	call mkareafvl			!area of finite volumes

        call get_flag(flag)
	p3 = flag

	write(6,*) 'allocate_simulation: ',nkn,nel,nlvdi,np

	end

! *****************************************************************
! *****************************************************************
! *****************************************************************

	subroutine get_dimension_post(nknddi,nelddi,nlvddi)

	implicit none

	integer nknddi,nelddi,nlvddi

	call get_dimension_post_2d(nknddi,nelddi)
	call get_dimension_post_3d(nlvddi)

	end

! *****************************************************************

	subroutine get_dimension_post_2d(nknddi,nelddi)

	use basin

	implicit none

	integer nknddi,nelddi

	call basin_get_dimension(nknddi,nelddi)

	end

! *****************************************************************

	subroutine get_dimension_post_3d(nlvddi)

	use levels

	implicit none

	integer nlvddi

	call levels_get_dimension(nlvddi)

	end

! *****************************************************************

