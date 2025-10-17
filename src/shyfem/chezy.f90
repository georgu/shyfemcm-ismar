
!--------------------------------------------------------------------------
!
!    Copyright (C) 1988,1990,1997-2000,2003,2005,2005  Georg Umgiesser
!    Copyright (C) 2007-2020  Georg Umgiesser
!    Copyright (C) 2012  Aaron Roland
!    Copyright (C) 2017,2019  Marco Bajo
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

! parameter changing area routines
!
! contents :
!
! function cdf(h,z0)		computes cd from h and z0
!
! subroutine rdarea		reads area section (chezy) from STR file
! subroutine ckarea		checks values for chezy parameters
! subroutine prarea		prints chezy values to log file
! subroutine tsarea		prints test message to terminal
! subroutine inarea		initializes chezy values
!
! revision log :
!
! 31.08.1988	ggu	(writes real chezy on czv)
! 29.11.1988	ggu	(new chezy, iarv array)
! 12.04.1990	ggu	(href)
! 03.06.1990	ggu	(austausch)
! 26.06.1997	ggu	(implicit none, useless parts deleted)
! 25.05.1998	ggu	documentation started
! 21.08.1998	ggu	xv eliminated
! 25.05.1999	ggu	new routine bofric
! 20.01.2000	ggu	common block /dimdim/ eliminated
! 09.08.2003	ggu	bofric now returns array with friction
! 10.08.2003	ggu	completely restructured, counter from 0 to nczdum
! 04.09.2003	ggu	bug fix for missing return in get_chezy_values
! 11.01.2005	ggu	ausv eliminated (was not used anymore)
! 02.04.2007	ggu	in check -> warning only for cz=0 and Chezy/Strickler
! 10.12.2008	ggu	re-organized, deleted sp135r(), use bottom_friction()
! 29.01.2009	ggu	ausdef eliminated (chezy(5,.) is not used anymore)
! 23.03.2010	ggu	changed v6.1.1
! 16.02.2011	ggu	new routines to deal with nodal area code
! 14.04.2011	ggu	changed VERS_6_1_22
! 21.06.2012	ggu&aar	new friction for mud module
! 26.06.2012	ggu	changed VERS_6_1_55
! 25.10.2013	ggu	changed VERS_6_1_68
! 18.06.2014	ggu	changed VERS_6_1_77
! 23.12.2014	ggu	changed VERS_7_0_11
! 19.01.2015	ggu	changed VERS_7_1_3
! 28.04.2015	ggu	czdef is default for all areas not given
! 12.05.2015	ggu	rewritten with modules and allocatable
! 21.05.2015	ggu	changed VERS_7_1_11
! 17.07.2015	ggu	changed VERS_7_1_80
! 20.07.2015	ggu	changed VERS_7_1_81
! 24.07.2015	ggu	changed VERS_7_1_82
! 15.04.2016	ggu	changed VERS_7_5_8
! 30.09.2016	ggu	changed VERS_7_5_18
! 10.04.2017	ggu	compute cd, normalized bottom stress and bottom stress
! 09.05.2017	ggu	bug fix for computing bottom stress
! 03.11.2017	mbj	new documentation for ireib
! 05.12.2017	ggu	changed VERS_7_5_39
! 26.04.2018	ggu	area code adjusted for mpi
! 11.05.2018	ggu	changed VERS_7_5_47
! 14.02.2019	ggu	changed VERS_7_5_56
! 16.02.2019	ggu	changed VERS_7_5_60
! 12.03.2019	mbj	new friction ireib=10
! 05.03.2020	ggu	documentation upgraded
! 02.04.2022	ggu	revisited for mpi, actual chezy now at position 0
! 02.04.2022	ggu	adjust_chezy() adjusted for multi-domain
! 12.04.2022	ggu	global bdebug and iczunit variable for debugging
! 22.03.2023	ggu	relax error conditions for nodes not in same domain
! 03.05.2023	ggu	avoid out of bound error in ckarea()
! 16.10.2025	ggu	new routine init_friction_from_file()
!
!***********************************************************
!***********************************************************
!***********************************************************

!==================================================================
        module chezy
!==================================================================

        implicit none

        integer, save :: nczdum = 0
        real, save, allocatable :: czdum(:,:)

	integer, parameter :: iczact = 0	!index of actual value of chezy
	integer, parameter :: ncztot = 6	!total number of entries

        integer, save :: nz_lines = 0
	character*80, save, allocatable :: cz_lines(:)

	logical, save :: bdebug = .false.
	integer, save :: iczunit = 0

! czdum(j,ia)
! ia is area code [0,namax]
! j is entry describing chezy for area code [0,ncztot]
! 	j == 0	actual chezy value for area ia
! 	j == 1	default chezy value for area ia
! 	j == 2	secondary chezy value for area ia (if k1,k2 are given)
! 	j == 3	node number k1
! 	j == 4	node number k2
! 	j == 5	viscosity for area ia (not used)
! 	j == 6	flag if secondary value is present (1) or not (0)
! nodes k1,k2 describe direction which is used to chose from chezy value

!==================================================================
        contains
!==================================================================

	subroutine chezy_init(n)

	integer n

	nczdum = n
	allocate(czdum(0:ncztot,0:n))

	czdum = -1.

	end subroutine chezy_init

!==================================================================
        end module chezy
!==================================================================

!
!-------------------------------------------------------------------
!
! DOCS  START   P_friction
!
! DOCS  FRICTION		Bottom friction
!
! The friction term in the momentum equations can be written as
! $Ru$ and $Rv$ where $R$ is the variable friction coefficient and
! $u,v$ are the velocities in $x,y$ direction respectively.
! The form of $R$ can be specified in various ways. The value of 
! |ireib| is choosing between the formulations. In the parameter
! input file a value $\lambda$ is specified that is used in 
! the formulas below. In a 2D simulation the Strickler (2) or the Chezy (3)
! formulation is the preferred option, while for a 3D simulation is is
! recommended to use the drag coefficient (5) or the roughness length
! formulation (6).
!
! |ireib|	Type of friction used (default 0):
!		\begin{description}
!		\item[0] No friction used
!		\item[1] $R=\lambda$ is constant
!		\item[2] $\lambda$ is the Strickler coefficient.
!			 In this formulation $R$ is written as
!			 $R = \frac{g}{C^2} \frac{\vert u \vert}{H}$
!			 with $C=k_s H^{1/6}$ and $\lambda=k_s$ is
!			 the Strickler coefficient. In the above
!			 formula $g$ is the gravitational acceleration,
!			 $\vert u \vert$ the modulus of the current velocity
!			 and $H$ the total water depth.
!		\item[3] $\lambda$ is the Chezy coefficient.
!			 In this formulation $R$ is written as
!			 $R = \frac{g}{C^2} \frac{\vert u \vert}{H}$
!			 and $\lambda=C$ is the Chezy coefficient.
!		\item[4] $R=\lambda/H$ with $H$ the total water depth. 
!			 This corresponds to a linear bottom friction.
!		\item[5] $\lambda$ is a constant drag coefficient and $R$ is
!			 computed as $R=\lambda\frac{\vert u \vert}{H}$.
!			 This corresponds to a quadratic bottom friction.
!		\item[6] $\lambda$ is the bottom roughness length and $R$ is
!			 computed through the formula
!			 $R=C\frac{\vert u \vert}{H}$ with 
!			 $C=\big(\frac{0.4}{log(\frac{\lambda+0.5H}
!			 {\lambda})}\big)^2$
!		\item[7] If $\lambda \geq 1$ it specifies the Strickler 
!			 coefficient (|ireib=2|), otherwise it specifies a 
!			 constant drag coefficient (|ireib=5|).
!		\item[8] The bottom roughness length computed with
!			 sedtrans (sediment transport module) is used
!			 to compute the friction (similar to 6).
!		\item[9] Experimental for fluid mud (no documentation).
!		\item[10] Hybrid formulation switching between quadratic (5)
!			  and linear (4) bottom friction. The velocity 
!			  below which linear friction is used has to be 
!			  given in |uvmin|.
!		\end{description}
! |czdef|	The default value for the friction parameter $\lambda$.
!		Depending on the value of |ireib| the coefficient $\lambda$
!		is representing linear friction, a constant drag coefficient,
!		the Chezy or Strickler parameter, or the roughness length.
!		(default 0)
! |iczv|	Normally the bottom friction coefficient 
!		(such as Strickler, Chezy, etc.)
!		is evaluated at every time step (|iczv| = 1).
!		If for some reason this behavior is not desirable,
!		|iczv| = 0 evaluates this value only before the
!		first time step, keeping it constant for the
!		rest of the simulation. Please note that this is
!		only relevant if you have given more than one bottom
!		friction value (inflow/outflow) for an area. The
!		final value of $R$ is computed at every time step
!		anyway. (default 1)
! |uvmin|	Critical velocity for |ireib|=10 below which bottom friction
!		will be used as linear. (Default 0.2)
!
! The value of $\lambda$ may be specified for the whole basin through
! the value of |czdef|. For more control over the friction parameter
! it can be also specified in section |area| where the friction
! parameter depending on the type of the element may be varied. Please
! see the paragraph on section |area| for more information.
!
! DOCS  END
!
! next are experimental settings
!
!		\item[8] As ireib = 6 but using a bottom roughness length 
!			computed by Sedtrans
!		\item[9] As ireib = 6 but using a bottom roughness length 
!			computed by the fluid mud module
!
!-------------------------------------------------------------------
!
!***********************************************************

	subroutine bottom_friction

! computes bottom friction
!
! rfric is given value (in czv)
! rcd is drag coefficient ( tau = rho * rcd * u**2 )
!
! for some formulations no rcd might exists

	use mod_fluidmud
	use mod_layer_thickness
	use mod_roughness
	use mod_diff_visc_fric
	use mod_hydro
	use levels
	use basin
	use pkonst

	implicit none

	real, parameter :: drittl = 1./3.

	integer ie,ii,k,lmax
	integer ireib
	real hzg,alpha
	real hzoff
	real uso,vso,uv,uuvv
	real rfric,rcd,rr,ss,rho0

	real getpar,cdf
	real uvmin

!-------------------------------------------------------------------
! get variables
!-------------------------------------------------------------------

	hzoff = getpar('hzoff')
	ireib = nint(getpar('ireib'))
	uvmin = getpar('uvmin')
	rho0 = rowass

!-------------------------------------------------------------------
! loop on elements
!-------------------------------------------------------------------

	do ie=1,nel

!         ----------------------------------------------------------
!	  get transport in layer
!         ----------------------------------------------------------

	  lmax = ilhv(ie)

          uso = utlov(lmax,ie)
          vso = vtlov(lmax,ie)
	  uv = sqrt(uso*uso+vso*vso)

!         ----------------------------------------------------------
!	  set total depth
!         ----------------------------------------------------------

	  hzg = hdeov(lmax,ie)
          if( hzg .lt. hzoff ) hzg = hzoff

!         ----------------------------------------------------------
!	  get friction parameter
!         ----------------------------------------------------------

	  rfric = czv(ie)

!         ----------------------------------------------------------
!	  compute friction
!         ----------------------------------------------------------

	  if(ireib.eq.0) then
		rr = 0.
		rcd = 0.
          else if(ireib.eq.1) then
                rr = rfric
		rcd = 0.
		if( uv > 0. ) rcd = rr*hzg*hzg/uv
	  else if(ireib.eq.2) then		! Strickler
		rcd = grav/((rfric**2)*(hzg**drittl))
		rr = rcd*uv/(hzg*hzg)
	  else if(ireib.eq.3) then		! Chezy
		rcd = grav/(rfric**2)
		rr = rcd*uv/(hzg*hzg)
          else if(ireib.eq.4) then		! linear friction
                rr = rfric/hzg
		rcd = 0.
		if( uv > 0. ) rcd = rr*hzg*hzg/uv
	  else if(ireib.eq.5) then		! quadratic friction, cd=const
		rcd = rfric
		rr = rcd*uv/(hzg*hzg)
	  else if(ireib.eq.6) then		! rfric is z0
                rcd = cdf(hzg,rfric)
		rr = rcd*uv/(hzg*hzg)
          else if(ireib.eq.7) then		! mixed Strickler / drag
                if( rfric .ge. 1. ) then
		  rcd = grav/((rfric**2)*(hzg**drittl))
                else
		  rcd = rfric
                end if
		rr = rcd*uv/(hzg*hzg)
          else if(ireib.eq.8) then		! use z0 computed by sedtrans
                ss = 0.
                do ii=1,3
                  k = nen3v(ii,ie)
                  ss = ss + z0bk(k)
                end do
                ss = ss / 3.
                rcd = cdf(hzg,ss)
		rr = rcd*uv/(hzg*hzg)
          else if(ireib.eq.9) then		! function of fluid mud (AR:)
                ss = 0.
                do ii=1,3
                  k = nen3v(ii,ie)
                  lmax = ilhkv(k)
                  call set_mud_roughness(k,lmax,alpha) ! (ARON)
                  ss = ss + alpha * rfric ! rfric = ks for this parameterization
                end do
                ss = ss / 3.
                z0bk(k) = ss
                !z0bk(k) = max(z0bkmud(k),ss)
                !ss = rfric	!ARON: do you really need to compute ss above?
                rcd = cdf(hzg,ss)
                rr = rcd*uv/(hzg*hzg)
		!Well not really there are mainls two issues ...
		!1. Rougnes get reduced by mud this is taken into 
		!account by calling the routine above
		!2. We need to apply mixing length for the 1st grid-cell 
		!otherwise turbulence in gotm fully collapse since k-eps 
		!is only valid for isotropic turbulence. 
	  else if(ireib.eq.10) then 	!Hybrid quadratic-linear formulation 
		!for more info see Bajo et al. 2019
                if(uv > uvmin*hzg) then	!quadratic (uvmin=0.2 by default)
                  rcd = rfric
                  rr = rcd*uv/(hzg*hzg)
                else			!linear
                  rr = (rfric*uvmin)/hzg
                  rcd = 0.
                  if( uv > 0. ) rcd = rr*hzg*hzg/uv
                end if
	  else
		write(6,*) 'unknown friction : ',ireib
		stop 'error stop bottom_friction'
	  end if

	  rfricv(ie) = rr
	  rcdv(ie) = rcd
	  uuvv = uv / hzg
	  bnstressv(ie) = rcd * uuvv * uuvv

	end do

!-------------------------------------------------------------------
! end of routine
!-------------------------------------------------------------------

	end

!***********************************************************

        function cdf(h,z0)

! computes cd from h and z0

        implicit none

        real cdf
        real h,z0

        real kappa,cds

        kappa = 0.4

        cds = kappa / log( (z0+0.5*h) / z0 )

        cdf = cds*cds

        end

!***********************************************************

	subroutine init_nodal_area_code

! interpolates area codes from elements to nodes (min or max)

	use basin
	use shympi

	implicit none

	integer init,mode
	integer k,ie,ii,ia

	mode = -1		! -1: use minimum   +1: use maximum

	init = 99999999
	if( mode .gt. 0 ) init = -init

	do k=1,nkn
	  iarnv(k) = init
	end do

	do ie=1,nel
	  ia = iarv(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( mode .eq. -1 ) then
		iarnv(k) = min(iarnv(k),ia)
	    else
		iarnv(k) = max(iarnv(k),ia)
	    end if
	  end do
	end do

	!call shympi_comment('exchanging iarnv')
	call shympi_exchange_2d_node(iarnv)

        if( mode .eq. -1 ) then
          call shympi_exchange_2d_nodes_min(iarnv)
          !call shympi_comment('shympi_elem: exchange iarnv_min')
        else
          call shympi_exchange_2d_nodes_max(iarnv)
          !call shympi_comment('shympi_elem: exchange iarnv_max')
        end if

!       shympi_elem:   exchange min or max iarnv

	end

!***********************************************************

	subroutine get_nodal_area_code(k,ncode)

	use basin

	implicit none

	integer k	!node number
	integer ncode	!nodal area code (return)

	ncode = iarnv(k)

	end

!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************

	subroutine n_chezy_values(nareas)

	use chezy

	implicit none

	integer nareas

	nareas = nczdum

	end

!***********************************************************

	subroutine get_chezy_values(iar,valin,valout)

	use chezy

	implicit none

	integer iar
	real valin,valout

	if( iar .gt. nczdum ) goto 99

	valin = czdum(1,iar)
	valout = czdum(2,iar)

	return
   99	continue
	write(6,*) 'iar,nczdum: ',iar,nczdum
	stop 'error stop get_chezy_values'
	end

!***********************************************************

	subroutine set_chezy_values(iar,valin,valout)

	use chezy

	implicit none

	integer iar
	real valin,valout

	if( iar .gt. nczdum ) goto 99

	czdum(1,iar) = valin
	czdum(2,iar) = valout

	return
   99	continue
	write(6,*) 'iar,nczdum: ',iar,nczdum
	stop 'error stop set_chezy_values'
	end

!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************

	subroutine set_chezy

! initializes chezy arrays

	use mod_diff_visc_fric
	use basin
	use chezy
	use shympi

	implicit none

	integer ie,iar
	real cmin,cmax

	do ie=1,nel
	    iar=iarv(ie)
	    czv(ie)=czdum(iczact,iar)
	end do

	!cmin = minval(czv)
	!cmax = maxval(czv)
	!write(6,*) 'chezy min/max: ',cmin,cmax

	end

!***********************************************************

	subroutine init_chezy

! initializes chezy arrays

	use chezy

	implicit none

	integer i

	do i=0,nczdum
	  czdum(iczact,i)=czdum(1,i)
	end do

	if( bdebug ) then
	  call print_chezy
	  call prarea
	end if

	call set_chezy

	call init_friction_from_file

	end

!***********************************************************

	subroutine adjust_chezy

! adjusts chezy arrays

	use mod_hydro_print
	use basin
	use chezy
	use shympi

	implicit none

	logical bchange,bstop,blocaldebug
	integer i,k1,k2,iflag,id,iu
	real dx,dy,scal
	real cz,czn,czold,cznew
	real, parameter :: flag = -999.
	real, allocatable :: czaux(:)
	real, allocatable :: czaux_all(:,:)

	integer, save :: iczv = 0
	integer, save :: icall = 0

	real getpar

	if( icall == -1 ) return

	if( icall == 0 ) then
	  iczv = nint(getpar('iczv'))
	  if( iczv == 0 ) icall = -1
	  if( icall == -1 ) return
	  icall = 1
	end if

	do i=0,nczdum
	  iflag=nint(czdum(6,i))
	  k1=nint(czdum(3,i))
	  if(iflag.eq.0) then		!do not use direction
	    czdum(iczact,i)=czdum(1,i)
	  else if( k1 == 0 ) then	!cannot determine direction in domain
	    czdum(iczact,i)=flag
	  else
	    k2=nint(czdum(4,i))
	    dx=xgv(k2)-xgv(k1)
	    dy=ygv(k2)-ygv(k1)
	    scal=dx*up0v(k1)+dy*vp0v(k1)
	    if(scal.ge.0.) then
	      cznew = czdum(1,i)
	    else
	      cznew = czdum(2,i)
	    end if
	    czdum(iczact,i) = cznew
	  end if
	end do

	bstop = .false.
	allocate(czaux(0:nczdum))
	allocate(czaux_all(0:nczdum,n_threads))
	czaux = czdum(iczact,:)
	call shympi_gather(czaux,czaux_all)
	do i=0,nczdum
	    cz = czaux_all(i,1)
	    do id=1,n_threads
	      czn = czaux_all(i,id)
	      if( czn == flag ) cycle
	      if( cz == flag ) cz = czn
	      if( cz /= czn ) then
		write(6,*) 'different values for chezy in area',i
		write(6,*) cz,czn
		bstop = .true.
	      end if
	    end do
	    if( cz == flag ) then
	      write(6,*) 'no value for chezy in area',i
	      bstop = .true.
	    end if
	    czaux(i) = cz
	    czdum(iczact,i) = cz
	end do

	blocaldebug = .false.
	if( blocaldebug ) then
	  iu = 700 + my_id
	  write(iu,*) (czdum(0,i),i=0,nczdum)
	end if

	if( bstop ) then
	  write(6,*) 'checzy values of local domain:'
	  call shympi_barrier
	  write(6,*) my_id,czaux
	  flush(6)
	  call shympi_barrier
	  call prarea
	  if( shympi_is_master() ) then
	    write(6,*) 'checzy values of all threads:'
	    do i=0,nczdum
	      write(6,*) i,czaux_all(i,:)
	    end do
	    flush(6)
	  end if
	  call shympi_barrier
	  stop 'error stop adjust_chezy: errors in chezy values'
	end if

	if( bdebug ) call print_chezy

	call set_chezy

	end

!***********************************************************

	subroutine print_chezy

! prints chezy arrays

	use chezy
	use shympi

	implicit none

	integer i
	integer iunit

	iunit = 6
	if( bdebug ) iunit = iczunit

	write(iunit,*) 'Values for chezy by print_chezy (czv) :'
	write(iunit,*) 'table size: ',nczdum
	do i=0,nczdum
	  write(iunit,*) i,czdum(iczact,i)
	end do
	flush(iunit)

	end

!***********************************************************

	subroutine check_chezy

! checks chezy arrays

	use basin
	use chezy

	implicit none

	integer ie,iar
	integer i,j,k
	real cz

	do i=0,nczdum
	  cz = czdum(1,i)
	  if( cz .lt. 0. .or. cz .gt. 1.e+10 ) goto 99
	  cz = czdum(2,i)
	  if( cz .lt. 0. .or. cz .gt. 1.e+10 ) goto 99
	  if( cz .gt. 0. ) then
	    k = nint(czdum(3,i))
	    if( k .lt. 0. .or. k .gt. nkn ) goto 99
	    k = nint(czdum(4,i))
	    if( k .lt. 0. .or. k .gt. nkn ) goto 99
	  end if
	  cz = czdum(iczact,i)
	  if( cz .lt. 0. .or. cz .gt. 1.e+10 ) goto 99
	end do

	do ie=1,nel
	    iar=iarv(ie)
	    if( iar .gt. nczdum ) goto 98
	    cz=czdum(iczact,iar)
	    if( cz .lt. 0. .or. cz .gt. 1.e+10 ) goto 98
	end do

	return
   98	continue
	write(6,*) 'ie,iar,nczdum,cz: ',ie,iar,nczdum,cz
	write(6,*) (czdum(j,iar),j=0,ncztot)
	stop 'error stop check_chezy: error in values (1)'
   99	continue
	write(6,*) 'i,iar,nczdum: ',i,i-1,nczdum
	write(6,*) (czdum(j,i),j=0,ncztot)
	stop 'error stop check_chezy: error in values (2)'
	end

!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************
!***********************************************************

	subroutine rdarea

! reads area section (chezy) from STR file

	use nls
	use chezy

	implicit none

	integer n,i

	n = nls_read_table()

	if( n > 0 ) then
	  nz_lines = n
	  allocate(cz_lines(n))
	  call nls_copy_char_vect(n,cz_lines)
	end if

	if( bdebug ) write(iczunit,*) 'chezy lines: ',n

	end

!***********************************************************

	subroutine parse_area

! parses area section (chezy) from STR file

	use chezy

	implicit none

	character*80 line
	integer ianz,iar,i,n,il
	real f(10)

	integer iscanf

	do il=1,nz_lines
	  line = cz_lines(il)
	  !write(6,*) il,trim(line)
	  ianz = iscanf(line,f,10)
	  if( ianz .gt. 0 ) then
	    iar = nint(f(1))
            if(iar.lt.0) goto 88
	    if( ianz .gt. 7 ) goto 86
            if(iar.gt.nczdum) then
	      write(6,*) 'warning: no such area code... ignoring ',iar
	      cycle
	    end if
	    do i=2,ianz
	      czdum(i-1,iar) = f(i)
	    end do
	  else if( ianz .lt. 0 ) then
			goto 98
	  end if
	end do

	if( nz_lines > 0 ) then
	  nz_lines = 0
	  deallocate(cz_lines)	!we dont need it anymore
	end if

	return
   86   continue
        write(6,*) 'Too much data on line'
        write(6,*) line
        stop 'error stop : rdarea'
   88   continue
        write(6,*) 'Negative area code = ',iar,' not allowed'
        write(6,*) line
        stop 'error stop : rdarea'
   98   continue
        write(6,*) 'Read error in line :'
	write(6,*) line
        stop 'error stop : rdarea'
	end

!***********************************************************

	subroutine ckarea

! checks values for chezy parameters and sets up czdum

	use basin
	use chezy
	use shympi

	implicit none

	integer i,knode,knodeh,ireib,nczmax
	integer ke1,ke2,ki1,ki2
	integer id1,id2
	logical bstop,bpos,busedirection
	real czdef

	integer ipint
	real getpar

	bstop = .false.

! get default values

        ireib=nint(getpar('ireib'))
	!if( ireib .le. 0 ) return

	bpos = ireib .gt. 1 .and. ireib .ne. 5	!must be strictly positive

        czdef=getpar('czdef')

! compute maximum value of area code

	nczmax = 0
	do i=1,nel
          if(iarv(i).gt.nczmax) nczmax=iarv(i)
        end do
	nczmax = shympi_max(nczmax)

! allocate and parse arrays

	call chezy_init(nczmax)
	call parse_area
	call shympi_syncronize

! check read in values

        do i=0,nczdum

         if(czdum(1,i).eq.-1.) czdum(1,i)=czdef

         if(czdum(1,i).lt.0.) then
                write(6,*) 'Friction value cannot be negative:'
                write(6,*) 'area = ',i,'  chezy = ',czdum(1,i)
                bstop=.true.
	 end if

	 if( bpos .and. czdum(1,i).eq.0.) then
                write(6,*) 'Friction value must be positive:'
                write(6,*) 'area = ',i,'  chezy = ',czdum(1,i)
                bstop=.true.
	 end if

         if(czdum(2,i).eq.-1.) then
	   czdum(2,i)=czdum(1,i)
           czdum(5,i)=0.
           czdum(6,i)=0.
	   busedirection = .false.
	 else
           czdum(6,i)=1.
	   busedirection = .true.
	 end if

         if(czdum(3,i).eq.-1. .or. czdum(3,i).eq.0.) then
           czdum(3,i)=0.
	   ke1 = 0
	   ki1 = 0
         else
           ke1=nint(czdum(3,i))
           ki1=ipint(ke1)          !$$EXTINW
           czdum(3,i)=ki1
         end if

         if(czdum(4,i).eq.-1. .or. czdum(4,i).eq.0.) then
           czdum(4,i)=0.
	   ke2 = 0
	   ki2 = 0
         else
           ke2=nint(czdum(4,i))
           ki2=ipint(ke2)          !$$EXTINW
           czdum(4,i)=ki2
         end if

	 if( ke1 == 0 .and. ke2 == 0 ) then
	   !no external nodes given
	   if( busedirection ) then
             write(6,*) 'section AREA : secondary chezy given ' //      &
     &			' but no nodes specified'
             bstop=.true.
	   end if
	 else if( ke1 > 0 .and. ke2 > 0 ) then
	   if( ki1 > 0 .and. ki2 > 0 ) then			!ok
	     !nodes inside domain
	   else if( ki1 == 0 .and. ki2 == 0 ) then		!other domain
	     if( bmpi ) then
               !write(6,*) 'section AREA : nodes not in domain ',ke1,ke2
	     else
               write(6,*) 'section AREA : nodes not found ',ke1,ke2
               bstop=.true.
	     end if
	   else
             write(6,*) 'section AREA : nodes in different domains ',ke1,ke2
	     id1 = -1
	     if( ki1 > 0 ) id1 = id_node(ki1)
	     id2 = -1
	     if( ki2 > 0 ) id2 = id_node(ki2)

	     if( ki1 > 0 .and. id1 /= my_id ) then
               write(6,*) 'section AREA : ids different... ignoring'
	       czdum(3,i) = 0
	     else if( ki2 > 0 .and. id2 /= my_id ) then
               write(6,*) 'section AREA : ids different... ignoring'
	       czdum(4,i) = 0
	     else
               write(6,*) 'section AREA : main domain... cannot ignore'
               write(6,*) 'section AREA : nodes probably not contiguous'
               write(6,*) 'section AREA : czdum: ',my_id,czdum(:,i)
               bstop=.true.
	     end if
	   end if
	 else
           write(6,*) 'section AREA : only one node given ',ke1,ke2
           bstop=.true.
	 end if

         czdum(5,i)=0.
         czdum(iczact,i)=0.

        end do

	call shympi_syncronize

	if( bstop ) stop 'error stop ckarea'

	end

!***********************************************************

	subroutine prarea

! prints chezy values to log file

	use chezy
	use shympi

	implicit none

	integer ianf,i
	integer ipext
	integer iunit

	iunit = 6
	if( bdebug ) iunit = iczunit

	if( .not. bdebug .and. .not. shympi_is_master() ) return
 
        ianf=0
        if(czdum(1,0).eq.0) ianf=1
        write(iunit,*)
        write(iunit,*) 'chezy table written by prarea:',nczdum
        write(iunit,1007)

        do i=ianf,nczdum
            if(czdum(2,i).ne.0.) then			!with two chezy
                write(iunit,1008) i,czdum(0,i)              &
     &				,czdum(1,i),czdum(2,i)                    &
     &                          ,ipext(nint(czdum(3,i)))    &
     &                          ,ipext(nint(czdum(4,i)))    &
     &                          ,nint(czdum(6,i))
            else					!just one chezy
                write(iunit,1008) i,czdum(1,i)
            end if
        end do

	flush(iunit)

	return
 1007   format(' area,cz,cz1,cz2,k1,k2,flag : ')
 1008   format(i5,3e12.4,3i7)
	end

!***********************************************************

	subroutine tsarea

! prints test message to terminal

	use chezy

	implicit none

	integer j,i

        write(6,*) '/chezy/'
        write(6,*) nczdum
        do i=0,nczdum
            write(6,'(1x,6e12.4)') (czdum(j,i),j=0,ncztot)
        end do

	end

!***********************************************************

	subroutine inarea

! initializes chezy values

	use chezy
	use shympi

	implicit none

	iczunit = 500 + my_id

	end

!***********************************************************

        subroutine init_friction_from_file

! initialization of conz from file

	use basin
	use chezy
	use mod_diff_visc_fric

        implicit none

	logical, parameter :: busedef = .false.		!use default in czdef
	logical, parameter :: bshyout = .false.		!writes shy file

        character*80 what
        character*80 file_init
        character*80 file
        integer nvar
        integer nlvddi,nlv
	integer iczv
	integer id
        double precision dtime
	real czdef
        real val0(1)
        real val(nkn)

	real getpar

	what = 'friction init'
	file_init = 'fricin'
        call getfnm(file_init,file)
        if( file == ' ' ) return

	iczv = nint(getpar('iczv'))
	if( iczv > 0 ) then
	  write(6,*) 'initializing friction values from files incompatible'
	  write(6,*) '  with time varying friction coefficient (iczv>0)'
	  stop 'error stop init_friction_from_file: iczv>0'
	end if

        dtime = 0.
	nvar = 1
        nlvddi = 1
        nlv = 1
        val0 = -1.

	czdef = getpar('czdef')
        if( busedef ) val0 = czdef	!set friction to default value

        call tracer_file_init(what,file_init,dtime &
     &                          ,nvar,nlvddi,nlv,nkn,val0,val)

	if( any(val==-1.) ) then
	  write(6,*) 'some nodes are not initialized with friction'
	  write(6,*) 'probably fem file given is too small -> extend'
	  stop 'error stop init_friction_from_file: missing values'
	end if

	call n2e2d(val,czv)	!averages nodal values to elemental values

	if( bshyout ) then
	  call shyfem_init_scalar_file('fric',1,.true.,id)
	  call shy_write_scalar_record2d(id,dtime,61,val)
	  call shy_close_output_file(id)
	  !call shyfem_init_elem_file('frice',1,.true.,id)
	  !call shy_write_elem_record2d(id,dtime,61,val)
	  !call shy_close_output_file(id)
	end if

        end

!***********************************************************

