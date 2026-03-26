
!--------------------------------------------------------------------------
!
!    Copyright (C) 2009-2017,2019  Georg Umgiesser
!    Copyright (C) 2009,2011  Christian Ferrarin
!    Copyright (C) 2009  Andrea Cucco
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

! tvd routines
!
! contents :
!
! subroutine tvd_init					initializes tvd scheme
! subroutine tvd_grad_3d(cc,gx,gy,aux,nlvdi,nlv)	computes gradients 3D
! subroutine tvd_grad_2d(cc,gx,gy,aux)			computes gradients 2D
! subroutine tvd_get_upwind_c(ie,l,ic,id,cu,cv)		c of upwind node
! subroutine tvd_upwind_init				init x,y of upwind node
! subroutine tvd_fluxes(ie,l,itot,isum,dt,cl,cv,gxv,gyv,f,fl) tvd fluxes
!
! revision log :
!
! 02.02.2009	ggu&aac	all tvd routines into seperate file
! 24.03.2009	ggu	bug fix: isum -> 6; declaration of cl() was missing 0
! 30.03.2009	ggu	bug fix: ilhv was real in tvd_get_upwind()
! 31.03.2009	ggu	bug fix: do not use internal gradient (undershoot)
! 06.04.2009	ggu&ccf	bug fix: in tvd_fluxes() do not test for conc==cond
! 23.03.2010	ggu	changed v6.1.1
! 15.12.2010	ggu	new routines for vertical tvd: vertical_flux_*()
! 27.01.2011	ggu	changed VERS_6_1_17
! 28.01.2011	ggu	bug fix for distance with lat/lon (tvd_fluxes)
! 29.01.2011	ccf	insert ISPHE for lat-long coordinates
! 17.02.2011	ggu	changed VERS_6_1_18
! 01.03.2011	ggu	changed VERS_6_1_20
! 23.03.2011	ccf	get isphe through get_coords_ev()
! 14.07.2011	ggu	changed VERS_6_1_27
! 24.11.2011	ccf	bug in tvd_init -> not resolved...
! 09.12.2011	ggu	changed VERS_6_1_38
! 24.01.2012	ggu	changed VERS_6_1_41
! 05.12.2013	ggu	changed VERS_6_1_70
! 27.06.2014	ggu	changed VERS_6_1_78
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 05.06.2015	ggu	changed VERS_7_1_12
! 13.07.2015	ggu	changed VERS_7_1_51
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 24.07.2015	ggu	changed VERS_7_1_82
! 18.09.2015	ggu	changed VERS_7_2_3
! 31.10.2016	ggu	initialization made faster
! 12.01.2017	ggu	changed VERS_7_5_21
! 16.02.2019	ggu	changed VERS_7_5_60
! 22.09.2020    ggu     correct warnings for PGI compiler
! 03.06.2022    ggu     documentation, adapted for mpi (only itvd==1 is working)
! 09.05.2023    lrp     introduce top layer index variable
! 28.08.2024    ggu     use find_unique_element() to find element
! 03.09.2024    ggu     more on mpi-tvd
! 26.09.2024    ggu     mark old tvd code with TVD_OLD, btvddebug introduced
! 29.09.2024    ggu     use ltvdup to decide if remote data exists
! 05.10.2024    ggu     introduced quad_tree_search()
!
!*****************************************************************
!
! notes :
!
! itvd = 0	no tvd
! itvd = 1	run tvd with gradient information using average
! itvd = 2	run tvd with gradient computed from up/down wind nodes
!
! itvd == 2 is the better scheme
!
! calling sequence :
!
! initialization
!	call tvd_init
!		if( btvd2 ) then
!		  call tvd_upwind_init_shell (for itvd==2)
!			call tvd_upwind_init
!		end if
! time loop (from newcon)
!	call scal3sh
!		if( btvd1 ) call tvd_grad_3d
!		call conz3d_omp
!			call vertical_flux_ie
!			call tvd_fluxes
!				if( btvd2 ) then
!				  call tvd_get_upwind_c(ie,l,ic,id,conu,cv)
!				end if
!
! look out for
!	is_external_boundary and mpi
!	itvd==1 is working in mpi
!	itvd==2 is not working in mpi: must get adjacent element information
!
!*****************************************************************
!*****************************************************************
!*****************************************************************
! tvd initialization
!*****************************************************************
!*****************************************************************
!*****************************************************************

        subroutine tvd_init

! initializes horizontal tvd scheme

	use basin
	use mod_tvd
	use shympi
	use mod_info_output

        implicit none

	integer itvd
	integer, save :: icall = 0
	!logical, save :: bdebug = .false.
	logical, save :: bdebug = .true.
	logical :: btvd2

	real getpar

	if( icall .ne. 0 ) return
	icall = 1

	itvd = nint(getpar('itvd'))

	itvd_type = itvd
	btvd2 = itvd == 2
	!write(6,*) 'tvd type: ',itvd_type

	if( btvd2 ) call mod_tvd_init(nel)

	if( btvd2 ) then
          if( shympi_is_parallel() ) then
	    call tvd_upwind_init_mpi
	    call tvd_mpi_init
	  else
	    call tvd_upwind_init_shell
          end if
	end if

	if( btvd2 ) then
	  if( bdebug ) call write_tvd_debug(nel) !only if btvddebug==.true.
	end if

	if( print_not_quiet_once() ) then
	if( itvd .eq. 0 ) then
	  write(6,*) 'no horizontal TVD scheme used'
	else
	  write(6,*) 'horizontal TVD scheme initialized: ',itvd
	end if
	end if

	end

!*****************************************************************

        subroutine tvd_upwind_init_shell

! initializes position of upwind node (shell) - final version

	use mod_tvd
	use basin

        implicit none

	logical bsphe,binit
	integer isphe
        integer ie,nthreads
	integer it1,idt
	integer, save :: ifreq = 1
	logical, save :: bwrite = .false.
	real proc
	!integer, save :: ifreq = nel/1000

        write(6,*) 'setting up tvd upwind information...'

	call get_coords_ev(isphe)
	bsphe = isphe .eq. 1

	call openmp_get_num_threads(nthreads)

	call get_clock_count(it1)

!$OMP PARALLEL DO PRIVATE(ie,proc) SHARED(nel,bsphe,bwrite,ifreq) DEFAULT(NONE)

          do ie=1,nel
            call tvd_upwind_init_elem(bsphe,ie)
	    if( bwrite .and. mod(ie,ifreq) == 0 ) then
	      proc = 100.*float(ie)/nel
	      write(6,'(2i10,f8.2)') ie,nel,proc
	    end if
	  end do

!$OMP END PARALLEL DO

	call is_init_ev(binit)

	call get_clock_count_diff(it1,idt)
	write(6,*) 'clock count: ',idt,nthreads
	write(6,*) 'binit,nel: ',binit,nel
        write(6,*) '...tvd upwind setup done (itvd=2)'

	end

!*****************************************************************

        subroutine tvd_upwind_init_shell_1

! initializes position of upwind node (shell) - omp not working

	use mod_tvd
	use basin

        implicit none

	logical bsphe
	integer isphe
        integer ie,ies,ieend,nchunk,nthreads
	integer it1,idt

	integer omp_get_num_threads

        write(6,*) 'setting up tvd upwind information...'

	call get_coords_ev(isphe)
	bsphe = isphe .eq. 1

	call openmp_get_num_threads(nthreads)

	call get_clock_count(it1)

!$OMP PARALLEL 
!$OMP SINGLE

	call omp_compute_chunk(nel,nchunk)	!nchunk == 1 if no OMP
	nthreads = 1
!$      nthreads = omp_get_num_threads()
!$	write(6,*) 'using chunk = ',nchunk,nel,nthreads

	do ie=1,nel,nchunk

!$OMP TASK FIRSTPRIVATE(ie)            PRIVATE(ies,ieend) &
!$OMP&     SHARED(nel,nchunk,bsphe)    DEFAULT(NONE)

	  call omp_compute_minmax(nchunk,nel,ie,ieend)

          do ies=ie,ieend
            call tvd_upwind_init_elem(bsphe,ie)
	  end do

!$OMP END TASK

        end do

!$OMP END SINGLE
!$OMP TASKWAIT  
!$OMP END PARALLEL      

	call get_clock_count_diff(it1,idt)
	write(6,*) 'clock count (old): ',idt,nthreads
        write(6,*) '...tvd upwind setup done (itvd=2)'

	end

!*****************************************************************

        subroutine tvd_upwind_init_shell_0

! initializes position of upwind node (shell) - new simplified version

	use mod_tvd
	use basin

        implicit none

	logical bsphe
	integer isphe
        integer ie,ies,ieend,nchunk,nthreads,nmax
	integer it1,idt

	integer omp_get_num_threads,OMP_GET_MAX_THREADS

        write(6,*) 'setting up tvd upwind information...'

	call get_coords_ev(isphe)
	bsphe = isphe .eq. 1

	call openmp_get_max_threads(nmax)
	call openmp_get_num_threads(nthreads)
	write(6,*) 'nthreads,nmax = ',nthreads,nmax

	call get_clock_count(it1)

!$	nchunk = nel/(10*nmax)
!$	write(6,*) 'max threads = ',nmax,nthreads,nchunk

!$OMP PARALLEL SHARED(nel,bsphe,nchunk) PRIVATE(ie)
!$OMP DO SCHEDULE(DYNAMIC,nchunk)

	do ie=1,nel
          call tvd_upwind_init_elem(bsphe,ie)
        end do

!$OMP END PARALLEL      

	call get_clock_count_diff(it1,idt)
	write(6,*) 'clock count (new): ',idt,nthreads
        write(6,*) '...tvd upwind setup done (itvd=2)'

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

        subroutine tvd_upwind_init_mpi

	use basin
	use mod_tvd
	use shympi
	use shympi_tvd

	implicit none

	logical bsphe
	integer isphe
	integer ie

	call get_coords_ev(isphe)
	bsphe = isphe .eq. 1

	do ie=1,nel
          call tvd_upwind_init_elem(bsphe,ie)
	  iatvdup(:,:,ie) = my_id + 1
	end do

	write(6,*) 'tvd_upwind_init_mpi successfully called',nel,my_id

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************

        subroutine tvd_upwind_init_elem(bsphe,ie)

! initializes position of upwind node for one element
!
! sets position and element of upwind node

	use mod_tvd
	use basin
	use levels
	use mod_quad_tree

        implicit none

	logical, intent(in) :: bsphe
	integer, intent(in) :: ie

	logical bdebug
        integer ii,j,k,in
	integer, save :: itot = 0
	integer, save :: itot2 = 0
        integer ienew,ienew2,ieq
	real x,y
	real r

        double precision xc,yc,xd,yd,xu,yu
        double precision dlat0,dlon0                    !center of projection

	logical in_element
	integer ieext

	bdebug = ie == 73
	bdebug = .true.
	bdebug = .false.

          xtvdup(:,:,ie) = 0.
          ytvdup(:,:,ie) = 0.
          ietvdup(:,:,ie) = 0
          ieetvdup(:,:,ie) = 0
          iatvdup(:,:,ie) = 0
          ltvdup(:,:,ie) = 0

          if ( bsphe ) call ev_make_center(ie,dlon0,dlat0)

          do ii=1,3

	    !if( bdebug ) write(6,*) 'looking for elem in tvd: ',ie,ii

            k = nen3v(ii,ie)
            xc = xgv(k)
            yc = ygv(k)
	    if ( bsphe ) call ev_g2c(xc,yc,xc,yc,dlon0,dlat0)

            j = mod(ii,3) + 1
            k = nen3v(j,ie)
            xd = xgv(k)
            yd = ygv(k)
	    if ( bsphe ) call ev_g2c(xd,yd,xd,yd,dlon0,dlat0)

            xu = 2*xc - xd
            yu = 2*yc - yd
	    if ( bsphe ) call ev_c2g(xu,yu,xu,yu,dlon0,dlat0)
	    x = xu
	    y = yu

	    call quad_tree_search(x,y,ieq)
	    ienew = ieq
	    !call find_unique_element(x,y,ienew)
	    if( ienew /= ieq ) stop 'error stop tvd_upwind_init_elem: ie/=ieq'
	    !call find_close_elem(ie,x,y,ienew)		!TVD_OLD

            xtvdup(j,ii,ie) = x
            ytvdup(j,ii,ie) = y
            if( ienew > 0 ) then
              ietvdup(j,ii,ie) = ienew
	      ieetvdup(j,ii,ie) = ipev(ienew)
	      ltvdup(j,ii,ie) = ilhv(ienew)
	    end if

            j = mod(ii+1,3) + 1
            k = nen3v(j,ie)
            xd = xgv(k)
            yd = ygv(k)
	    if ( bsphe ) call ev_g2c(xd,yd,xd,yd,dlon0,dlat0)

            xu = 2*xc - xd
	    yu = 2*yc - yd
	    if ( bsphe ) call ev_c2g(xu,yu,xu,yu,dlon0,dlat0)
	    x = xu
	    y = yu

	    call quad_tree_search(x,y,ieq)
	    ienew = ieq
	    !call find_unique_element(x,y,ienew)
	    if( ienew /= ieq ) stop 'error stop tvd_upwind_init_elem: ie/=ieq'
	    !call find_close_elem(ie,x,y,ienew)		!TVD_OLD

            xtvdup(j,ii,ie) = x
            ytvdup(j,ii,ie) = y
            if( ienew > 0 ) then
              ietvdup(j,ii,ie) = ienew
	      ieetvdup(j,ii,ie) = ipev(ienew)
	      ltvdup(j,ii,ie) = ilhv(ienew)
	    end if

          end do

	  !if( bdebug ) stop
        end

!*****************************************************************
!*****************************************************************
!*****************************************************************
! horizontal tvd schemes
!*****************************************************************
!*****************************************************************
!*****************************************************************

        subroutine tvd_grad_3d(cc,gx,gy,aux,nlvddi)

! computes gradients for scalar cc (average gradient information)
! this is for itvd == 1 only
!
! output is gx,gy

	use evgeom
	use levels
	use basin
	use shympi

        implicit none

	integer nlvddi
	real cc(nlvddi,nkn)	!scalar - input
	real gx(nlvddi,nkn)	!gradient information in x (return)
	real gy(nlvddi,nkn)	!gradient information in y (return)
	real aux(nlvddi,nkn)	!aux array - not used outside
        
        integer k,l,ie,ii,lmax,ie_mpi
	real b,c,area
	real ggx,ggy

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    gx(l,k) = 0.
	    gy(l,k) = 0.
	    aux(l,k) = 0.
	  end do
	end do

        do ie_mpi=1,nel
	  ie = ip_sort_elem(ie_mpi)
          area=ev(10,ie) 
	  lmax = ilhv(ie)
	  do l=1,lmax
            ggx=0
            ggy=0
            do ii=1,3
              k=nen3v(ii,ie)
              b=ev(ii+3,ie)
              c=ev(ii+6,ie)
              ggx=ggx+cc(l,k)*b
              ggy=ggy+cc(l,k)*c
              aux(l,k)=aux(l,k)+area
	    end do
            do ii=1,3
             k=nen3v(ii,ie)
             gx(l,k)=gx(l,k)+ggx*area
             gy(l,k)=gy(l,k)+ggy*area
            end do 
          end do
        end do

        do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax
	    area = aux(l,k)
	    if( area .gt. 0. ) then
	      gx(l,k) = gx(l,k) / area
	      gy(l,k) = gy(l,k) / area
	    end if
	  end do
        end do

	call shympi_exchange_3d_node(gx)
	call shympi_exchange_3d_node(gy)

        end
        
!*****************************************************************

        subroutine tvd_grad_2d(cc,gx,gy,aux)

! computes gradients for scalar cc (only 2D - used in sedi3d)
! this is for itvd == 1 only

	use evgeom
	use basin

        implicit none

	real cc(nkn)
	real gx(nkn)
	real gy(nkn)
	real aux(nkn)

        integer k,ie,ii
	real b,c,area
	real ggx,ggy

	do k=1,nkn
	  gx(k) = 0.
	  gy(k) = 0.
	  aux(k) = 0.
	end do

        do ie=1,nel
          area=ev(10,ie) 
          ggx=0
          ggy=0
          do ii=1,3
              k=nen3v(ii,ie)
              b=ev(ii+3,ie)
              c=ev(ii+6,ie)
              ggx=ggx+cc(k)*b
              ggy=ggy+cc(k)*c
              aux(k)=aux(k)+area
	  end do
          do ii=1,3
             k=nen3v(ii,ie)
             gx(k)=gx(k)+ggx*area
             gy(k)=gy(k)+ggy*area
          end do 
        end do

        do k=1,nkn
	    area = aux(k)
	    if( area .gt. 0. ) then
	      gx(k) = gx(k) / area
	      gy(k) = gy(k) / area
	    end if
        end do

        end
        
!*****************************************************************

        subroutine tvd_get_upwind_c(ie,l,ic,id,cu,cv)

! computes concentration of upwind node (using info on upwind node)

	use mod_tvd
	use levels
	use basin

        implicit none

        integer, intent(in) :: ie,l
	integer, intent(in) :: ic,id
        real, intent(inout) :: cu
        real, intent(in) :: cv(nlvdi,nkn)

	logical bdebug
        integer ienew
        integer ii,k,ipoint,lremote
        real xu,yu
        real c(3)
	real cu2
	double precision xd,yd,xi(3)
	real, parameter :: eps = 1.e-5

	bdebug = ipev(ie) == 2247 .and. l == 1
	bdebug = .false.

	ipoint = itvdup(id,ic,ie)
	lremote = ltvdup(id,ic,ie)

	if( ipoint > 0 ) then			!value computed in other doamin
	  if( lremote .lt. l ) return		!no such level on remote
	  call tvd_get_upwind(l,ipoint,cu)
	  if( bdebug ) write(6,*) 'upwind1: ',ipev(ie),ipev(ienew),ipoint,cu
	  return
	end if

        ienew = ietvdup(id,ic,ie)

        if( ienew .le. 0 ) return
	if( lremote /= ilhv(ienew) ) stop 'error stop tvd_get_upwind_c: lremote'
	if( lremote .lt. l ) return		!no such level in ienew element

        do ii=1,3
          k = nen3v(ii,ienew)
          c(ii) = cv(l,k)
        end do

        xu = xtvdup(id,ic,ie)
        yu = ytvdup(id,ic,ie)

        !call femintp(ienew,c,xu,yu,cu)		!TVD_OLD
        call femintp_xi(ienew,c,xu,yu,cu)

	if( bdebug ) then
	  xd = xu
	  yd = yu
	  call xy2xi(ienew,xd,yd,xi)
	  write(6,*) 'upwind2: ',ipev(ie),ipev(ienew),ipoint,cu
	end if

        end

!*****************************************************************

	subroutine tvd_fluxes(ie,l,itot,isum,dt,cl,cv,gxv,gyv,f,fl)

! computes horizontal tvd fluxes for one element
!
! this is called for itvd == 1 and itvd == 2
! in case itvd == 1 the values gxv,gyv are used to compute grad
! otherwise (itvd==2) grad is computed as grad = cond - conu

	use mod_tvd
	use mod_hydro_vel
	use evgeom
	use levels, only : nlvdi,nlv
	use basin

	implicit none

	integer, intent(in) :: ie,l
	integer, intent(in) :: itot,isum
	double precision, intent(in) :: dt
	double precision, intent(in) :: cl(0:nlvdi+1,3)		!bug fix
	real, intent(in) :: cv(nlvdi,nkn)
        real, intent(in) :: gxv(nlvdi,nkn)
        real, intent(in) :: gyv(nlvdi,nkn)
	double precision, intent(in) :: f(3)
	double precision, intent(out) :: fl(3)

	real eps
	parameter (eps=1.e-8)

        logical btvd2,btvddebug
        logical bdebug
	integer ii,k
        integer ic,kc,id,kd,ip,iop
	integer itot1,itot2
	integer tet1
        real term,fact,grad
        real conc,cond,conf,conu
        real gcx,gcy,dx,dy
        real u,v
        real rf,psi
        real alfa,dis,aj
        real vel
        real gdx,gdy
	real conu_aux(3)

	integer smartdelta

	btvd2 = itvd_type .eq. 2
	btvddebug = .true.
	btvddebug = btvddebug .and. btvd2
	bdebug = .true.
	bdebug = .false.

	if( bdebug ) then
	  write(6,*) 'tvd: ',ie,l,itot,isum,dt
	  write(6,*) 'tvd: ',btvd2,itvd_type
	end if

	fl = 0.
	conu_aux = 0.

	  if( itot .lt. 1 .or. itot .gt. 2 ) return

	  itot2 = itot - 1
	  itot1 = 2 - itot

	  u = ulnv(l,ie)
          v = vlnv(l,ie)
	  aj = 24 * ev(10,ie)

            ip = isum
            !if( itot .eq. 2 ) ip = 6 - ip		!bug fix
	    ip = itot2*(6-ip) + itot1*ip

            do ii=1,3
              if( ii .ne. ip ) then
                !if( itot .eq. 1 ) then			!flux out of one node
                !  ic = ip
                !  id = ii
                !  fact = 1.
                !else					!flux into one node
                !  id = ip
                !  ic = ii
                !  fact = -1.
                !end if
                ic = itot2*ii + itot1*ip
		id = itot2*ip + itot1*ii
		fact = -itot2 + itot1

                kc = nen3v(ic,ie)
                conc = cl(l,ic)
                kd = nen3v(id,ie)
                cond = cl(l,id)

                !dx = xgv(kd) - xgv(kc)
                !dy = ygv(kd) - ygv(kc)
                !dis = sqrt(dx**2 +dy**2)
		! next is bug fix for lat/lon
		iop = 6 - (id+ic)			!opposite node of id,ic
		tet1 = 1+mod(iop,3)
		dx = aj * ev(6+iop,ie)
		!if( tet1 .eq. id ) dx = -dx
		dx = -2*smartdelta(tet1,id) * dx + dx
		dy = aj * ev(3+iop,ie)
		!if( tet1 .eq. ic ) dy = -dy
		dy = -2*smartdelta(tet1,ic) * dy + dy
		dis = ev(16+iop,ie)

                vel = abs( u*dx + v*dy ) / dis          !projected velocity
                alfa = ( dt * vel  ) / dis

                if( btvd2 ) then
                  conu = cond
                  !conu = 2.*conc - cond		!use internal gradient
                  call tvd_get_upwind_c(ie,l,ic,id,conu,cv)
		  conu_aux(ii) = conu
                  grad = cond - conu
                else
                  gcx = gxv(l,kc)
                  gcy = gyv(l,kc)
                  grad = 2. * (gcx*dx + gcy*dy)
                end if

                if( abs(conc-cond) .lt. eps ) then	!BUG -> eps
                  rf = -1.
                else
                  rf = grad / (cond-conc) - 1.
                end if

                psi = max(0.,min(1.,2.*rf),min(2.,rf))  ! superbee
!               psi = ( rf + abs(rf)) / ( 1 + abs(rf))  ! muscl
!               psi = max(0.,min(2.,rf))                ! osher
!               psi = max(0.,min(1.,rf))                ! minmod

                conf = conc + 0.5*psi*(cond-conc)*(1.-alfa)
                term = fact * conf * f(ii)
                fl(ic) = fl(ic) - term
                fl(id) = fl(id) + term
              end if
            end do

	if( btvddebug ) call tvd_debug_accum(ie,l,conu_aux)

	if( bdebug ) then
	  write(6,*) 'tvd: --------------'
	  write(6,*) 'tvd: ',vel,gcx,gcy,grad
	  write(6,*) 'tvd: ',rf,psi,alfa
	  write(6,*) 'tvd: ',conc,cond,conf
	  write(6,*) 'tvd: ',term,fact
	  write(6,*) 'tvd: ',f
	  write(6,*) 'tvd: ',fl
	  write(6,*) 'tvd: ',(cl(l,ii),ii=1,3)
	  write(6,*) 'tvd: ',ic,id,kc,kd
	  write(6,*) 'tvd: --------------'
	end if

	end

!*****************************************************************
!*****************************************************************
!*****************************************************************
! vertical tvd schemes
!*****************************************************************
!*****************************************************************
!*****************************************************************

	subroutine vertical_flux_k(btvdv,k,dt,wsink,cv,vvel,vflux)

! computes vertical fluxes of concentration - nodal version

! do not use this version - use the element version instead !!!!

! ------------------- l-2 -----------------------
!      u              l-1
! ------------------- l-1 -----------------------
!      c               l     ^        d
! ---------------+---  l  ---+-------------------
!      d         v    l+1             c
! ------------------- l+1 -----------------------
!                     l+2             u
! ------------------- l+2 -----------------------

	use mod_layer_thickness
	use mod_hydro_print
	use levels
	use basin, only : nkn,nel,ngr,mbw

	implicit none

	logical btvdv			!use vertical tvd?
	integer k			!node of vertical
	real dt				!time step
	real wsink			!sinking velocity (positive downwards)
	real cv(nlvdi,nkn)		!scalar to be advected
	real vvel(0:nlvdi)		!velocities at interface (return)
	real vflux(0:nlvdi)		!fluxes at interface (return)

        real eps
        parameter (eps=1.e-8)

	integer l,lmax,lu
	real w,fl
	real conc,cond,conu,conf
	real hdis,alfa,rf,psi

	lmax = ilhkv(k)

	do l=1,lmax-1
	  w = wprv(l,k) - wsink

	  if( w .gt. 0. ) then
	    conc = cv(l+1,k)
	  else
	    conc = cv(l,k)
	  end if

	  conf = conc

	  if( btvdv ) then
	    if( w .gt. 0. ) then
	      !conc = cv(l+1,k)
	      cond = cv(l,k)
	      conu = cond
	      lu = l + 2
	      if( lu .le. lmax ) conu = cv(lu,k)
	    else
	      !conc = cv(l,k)
	      cond = cv(l+1,k)
	      conu = cond
	      lu = l - 1
	      if( lu .ge. 1 ) conu = cv(lu,k)
	    end if

	    hdis = 0.5*(hdknv(l,k)+hdknv(l+1,k))
	    alfa = dt * abs(w) / hdis
            if( abs(conc-cond) .lt. eps ) then
              rf = -1.
            else
              rf = (cond-conu) / (cond-conc) - 1.
            end if
            psi = max(0.,min(1.,2.*rf),min(2.,rf))  ! superbee
            conf = conc + 0.5*psi*(cond-conc)*(1.-alfa)
	  end if

	  vvel(l) = w
	  vflux(l) = w * conf
	end do

	vvel(0) = 0.			!surface
	vflux(0) = 0.
	vvel(lmax) = 0.			!bottom
	vflux(lmax) = 0.

	end

!*****************************************************************

	subroutine vertical_flux_ie(btvdv,ie,lmax,lmin,dt,wsink &
     &					,cl,wvel,hold,vflux)

! computes vertical fluxes of concentration - element version

! ------------------- l-2 -----------------------
!      u              l-1
! ------------------- l-1 -----------------------
!      c               l     ^        d
! ---------------+---  l  ---+-------------------
!      d         v    l+1             c
! ------------------- l+1 -----------------------
!                     l+2             u
! ------------------- l+2 -----------------------

	use levels, only : nlvdi,nlv

	implicit none

	logical btvdv				!use vertical tvd?
	integer ie				!element
	integer lmax				!index of bottom layer
	integer lmin				!index of top layer
	double precision dt			!time step
	double precision wsink			!sinking velocity (+ downwards)
	double precision cl(0:nlvdi+1,3)	!scalar to be advected
	double precision hold(0:nlvdi+1,3)	!depth of layers
	double precision wvel(0:nlvdi+1,3)	!velocities at interface
	double precision vflux(0:nlvdi+1,3)	!fluxes at interface (return)

        double precision eps
        parameter (eps=1.e-8)

	integer ii,l,lu
	double precision w,fl
	double precision conc,cond,conu,conf
	double precision hdis,alfa,rf,psi

	double precision, parameter :: zero = 0.
	double precision, parameter :: one = 1.
	double precision, parameter :: two = 2.
	double precision, parameter :: half = 1./2.

	do ii=1,3
	 do l=lmin,lmax-1
	  w = wvel(l,ii) - wsink

	  if( w .gt. 0. ) then
	    conc = cl(l+1,ii)
	  else
	    conc = cl(l,ii)
	  end if

	  conf = conc

	  if( btvdv ) then
	    if( w .gt. 0. ) then
	      !conc = cl(l+1,ii)
	      cond = cl(l,ii)
	      conu = cond
	      lu = l + 2
	      if( lu .le. lmax ) conu = cl(lu,ii)
	    else
	      !conc = cl(l,ii)
	      cond = cl(l+1,ii)
	      conu = cond
	      lu = l - 1
	      if( lu .ge. lmin ) conu = cl(lu,ii)
	    end if

	    hdis = 0.5*(hold(l,ii)+hold(l+1,ii))
	    alfa = dt * abs(w) / hdis
            if( abs(conc-cond) .lt. eps ) then
              rf = -1.
            else
              rf = (cond-conu) / (cond-conc) - 1.
            end if
            psi = max(zero,min(one,two*rf),min(two,rf))  ! superbee
            conf = conc + half*psi*(cond-conc)*(one-alfa)
	  end if

	  vflux(l,ii) = w * conf
	 end do
	 vflux(lmin-1,ii) = 0.
	 vflux(lmax,ii) = 0.
	end do

	end

!*****************************************************************

	function smartdelta(a,b)

	implicit none

	integer smartdelta
        integer, intent(in) :: a,b

        smartdelta=int((float((a+b)-abs(a-b)))/(float((a+b)+abs(a-b))))

	end function smartdelta

!*****************************************************************

	subroutine write_debug_grid(nelems,ielems,x,y)

	use basin

	implicit none

	integer nelems
	integer ielems(nelems)
	real x,y

	integer, save :: ncount = 0
	integer i
	integer ie_ext,ie_int,k,ii,node,iu
	integer ies(3)
	real xe,ye
	character*80 string,file

	integer ieint

	ncount = ncount + 1

	write(string,'(i4)') ncount
	do i=1,4
	  if( string(i:i) == ' ' ) string(i:i) = '0'
	end do
	file='debug_grid_'//trim(string)//'.grd'

	write(6,*) 'writing grid ',trim(file)

	iu = 1
	open(iu,file=file,status='unknown',form='formatted')

	node = 0

	do i=1,nelems
	  ie_ext = ielems(i)
	  ie_int = ieint(ie_ext)
	  ie_int = ie_ext
	  if( ie_int == 0 ) stop 'error stop write_debug_grid: no internal'
	  do ii=1,3
	    node = node + 1
	    ies(ii) = node
	    k = nen3v(ii,ie_int)
	    xe = xgv(k)
	    ye = ygv(k)
	    write(iu,'(i1,2i8,2f16.6)') 1,node,0,xe,ye
	  end do
	  write(iu,'(i1,6i8)') 2,i,0,3,ies(:)
	end do

	node = node + 1
	write(iu,'(i1,2i8,2f16.6)') 1,node,4,x,y

	close(iu)

	end

!*****************************************************************

