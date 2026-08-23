
!--------------------------------------------------------------------------
!
!    Copyright (C) 1994,1996,2015-2019  Georg Umgiesser
!    Copyright (C) 2015  Erik Pascolo
!    Copyright (C) 2016  Christian Ferrarin
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
! 09.01.1994	ggu	(from scratch)
! 19.01.1994	ggu	$$flux - flux conserving property
! 20.01.1994	ggu	$$iclin - iclin not used to compute volume
! 20.01.1994	ggu	$$lumpc - evaluate conz. nodewise
! 03.02.1994	ggu	$$itot0 - exception for itot=0 or 3
! 04.02.1994	ggu	$$fact3 - factor 3 missing in transport
! 04.02.1994	ggu	$$azpar - azpar used to compute transport
! 04.02.1994	ggu	$$condry - comute conz also in dry areas
! 07.02.1994	ggu	$$istot - istot for fractional time step
! 01.06.1994	ggu	restructured for 3-d model
! 18.07.1994	ggu	$$htop - use htop instead of htopo for mass cons.
! 09.04.1996	ggu	$$rvadj adjust rv in certain areas
! 20.05.2015	erp	transformed for OMP
! 30.09.2015	ggu	routine cleaned, no reals in conz3d
! 12.10.2015	ggu	changed VERS_7_3_3
! 22.10.2015	ggu	changed VERS_7_3_7
! 23.10.2015	ggu	changed VERS_7_3_9
! 20.11.2015	ggu&erp	chunk size introduced, omp finalized
! 18.12.2015	ggu	changed VERS_7_3_17
! 19.02.2016	ggu	changed VERS_7_5_2
! 20.10.2016	ccf	pass rtauv for differential nudging
! 12.01.2017	ggu	changed VERS_7_5_21
! 05.12.2017	ggu	changed VERS_7_5_39
! 11.05.2018	ggu	compute only unique nodes (needed for zeta layers)
! 06.07.2018	ggu	changed VERS_7_5_48
! 11.10.2018	ggu	code adjusted for sediment deposition (negative loads)
! 01.02.2019	ggu	bug fix for conz==0 with negative loading
! 14.02.2019	ggu	bug fix for conz<0 with negative loading
! 16.02.2019	ggu	changed VERS_7_5_60
! 13.03.2019	ggu	changed VERS_7_5_61
! 09.05.2023    lrp     introduce top layer index variable
! 24.09.2024    ggu     introduced tvd mpi (bmpi)
! 21.04.2026    ggu     use bdry to see if mfluxv should be used
!
!**************************************************************

        subroutine conz3d_omp(curr_stage &
     &			,coeff_erk,coeff_srk,coeff_crk &
     &			,cc,co &
     &			,ddt &
     &                  ,rkpar,difhv,difv &
     &			,difmol,cbound &
     &		 	,itvd,itvdv,gradxv,gradyv &
     &			,cobs,robs,rtauv &
     &			,wsink,wsinkv &
     &			,rload,load &
     &			,istot,isact,nlvddi &
     &                  ,nlev &
     &			,erk_reg,crk_reg,srk_reg)
     
! computes concentration
!
! cc     current stage concentration
! co     old concentration
! cn     new concentration
! caux   aux vector
! clow	 lower diagonal of vertical system
! chig	 upper diagonal of vertical system
! ddt    time step
! rkpar  horizontal turbulent diffusivity
! difhv  horizontal turbulent diffusivity (variable between elements)
! difv   vertical turbulent diffusivity
! difmol vertical molecular diffusivity
! cbound boundary condition (mass flux) [kg/s] -> now concentration [kg/m**3]
! itvd	 type of horizontal transport algorithm used
! itvdv	 type of vertical transport algorithm used
! gradxv,gradyv  gradient vectors for TVD algorithm
! cobs	 observations for nudging
! robs	 use observations for nuding (real)
! rtauv	 variable relaxation coefficient (real)
! wsink	 factor for settling velocity
! wsinkv variable settling velocity [m/s]
! rload	 factor for loading
! load   load (source or sink) [kg/s]
! curr_stage current stage index
! coeff_erk  explicit runge-kutta coefficients at current stage
! coeff_irk  implicit runge-kutta coefficients at current stage
! coeff_srk  stiffly implicit runge-kutta coefficients at current stage
! coeff_crk  special runge-kutta coefficients at current stage for tracer
! istot	 total inter time steps
! isact	 actual inter time step
! nlvddi	 dimension in z direction
! nlv	 actual needed levels
!
! solution of purely diffusional part :
!
! dC/dt = a*laplace(C)    with    c(x,0+)=delta(x)
!
! C(x,t) =  (4*pi*a*t)**(-n/2) * exp( -|x|**2/(4*a*t) )
!
! for n-dimensions and
!
! C(x,t) =  1/sqrt(4*pi*a*t) * exp( -x**2/(4*a*t) )
!
! for 1 dimension
!
! the solution is normalized, i.e.  int(C(x,t)dx) = 1 over the whole area
!
! DPGGU -> introduced double precision to stabilize solution

	use mod_rungekutta, only : n_rkstages
	use mod_bound_geom
	use mod_geom
	use mod_depth
	use mod_diff_aux
	use mod_bound_dynamic
	use mod_area
	use mod_ts
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin
!$	use omp_lib
	use mod_subset
	use shympi

	implicit none

! arguments
	integer, intent(in) :: curr_stage,nlvddi,nlev,itvd,itvdv,istot,isact
	real, intent(in) :: difmol,robs,wsink,rload,ddt,rkpar
	real,dimension(n_rkstages),intent(in) :: coeff_erk
	real,dimension(n_rkstages+1),intent(in) :: coeff_srk
	real,dimension(n_rkstages+1),intent(in) :: coeff_crk
	real,dimension(nlvddi,nkn),intent(inout) :: cc
	real,dimension(nlvddi,nkn),intent(in) :: co,cbound
	real,dimension(nlvddi,nel),intent(in) :: difhv
	real,dimension(nlvddi,nkn),intent(in) :: gradxv,gradyv
	real,dimension(nlvddi,nkn),intent(in) :: cobs,rtauv
	real,dimension(nlvddi,nkn),intent(inout) :: load		!LLL
	real,dimension(0:nlvddi,nkn),intent(in) :: difv,wsinkv
	real,dimension(nlvdi,nkn,n_rkstages-1),intent(inout) :: erk_reg !runge-kutta register array for explicit terms
	real,dimension(nlvdi,nkn,n_rkstages-1),intent(inout) :: crk_reg !runge-kutta register array for explicit terms
	real,dimension(nlvdi,nkn,n_rkstages-1),intent(inout) :: srk_reg !runge-kutta register array for explicit terms
        
	logical :: btvdv,btvd2,is_rk_explicit
	integer :: ie,k,ilevel,ibase,ii,l,n,i,j,x,ies,iend,kl,kend,ntot
	integer :: myid,numthreads,j_init,j_end,knod,k_end,jel
	integer,allocatable,dimension(:) :: subset_l
	real :: time1,time2
	double precision :: dtime1,dtime2
	integer :: nchunk,nthreads,nelems,nnodes
	double precision :: dt
	double precision :: rstot,rso,rsn,rsot,rsnt
	double precision :: timer,timer1,chunk,rest
	
	double precision,dimension(:,:),allocatable :: cn
        double precision,dimension(:,:),allocatable :: cdiag
	double precision,dimension(:,:),allocatable :: clow
	double precision,dimension(:,:),allocatable :: chigh

	! here debug code ------------------------
	logical, save :: bdebug = .false.
	double precision dtime
	integer iudb,ks,lmax
	integer ipint,ipext
	! here debug code ------------------------

        if(nlv.ne.nlev) stop 'error stop conz3d_omp: nlv/=nlev'
	
!----------------------------------------------------------------
! initialize variables and parameters
!----------------------------------------------------------------

!	call cpu_time(time1)
!!$	dtime1 = omp_get_wtime()
	
	ALLOCATE(cn(nlvddi,nkn))
	ALLOCATE(cdiag(nlvddi,nkn))
	ALLOCATE(clow(nlvddi,nkn))
	ALLOCATE(chigh(nlvddi,nkn))

	rstot = istot			!ERIC - what a brown paper bag bug
	rso=(isact-1)/rstot
	rsn=(isact)/rstot
	rsot=1.-rso
	rsnt=1.-rsn

	is_rk_explicit = (coeff_crk(curr_stage+1).eq.0.) &
     &		   .and. (coeff_srk(curr_stage+1).eq.0.)

	dt=ddt/rstot
	
	btvdv = itvdv .gt. 0
	if( btvdv .and. coeff_crk(curr_stage+1) .ne. 0. ) then
	  write(6,*) 'aapar = ',coeff_crk(curr_stage+1)
	  write(6,*) 'itvdv = ',itvdv
	  write(6,*) 'Cannot use vertical TVD scheme'
	  write(6,*) 'together with implicit vertical advection.'
	  write(6,*) 'Please set:'
	  write(6,*) 'itvdv = 0 (no vertical TVD) in the STR file.'
	  write(6,*) 'or use an implicit vertical advection scheme.'
	  stop 'error stop conz3d: vertical tvd scheme'
	end if

	btvd2 = itvd == 2
	if( btvd2 ) call tvd_mpi_run(cc)

        cn=0.
        cdiag=0.
        clow=0.
        chigh=0.

	if (curr_stage .ne. n_rkstages) then 	!if (not last stage)
	  erk_reg(:,:,curr_stage) = 0.		!zeroing-out rk register vectors (they are cumulated)
	  crk_reg(:,:,curr_stage) = 0.
	  srk_reg(:,:,curr_stage) = 0.
	end if

        !call tsdebug(robs,cobs,rtauv)

        nchunk = 1
	nthreads = 1
!!!!$	nthreads = omp_get_num_threads()
	call openmp_get_num_threads(nthreads)
 
	!write(6,*) 'subset_num: ',subset_num,nthreads
	!write(6,*) subset_el
	!write(6,*) sum(subset_el),nkn,nel

      do i=1,subset_num 	! loop over indipendent subset
       
!$     nchunk = subset_el(i) / ( nthreads * 10 )
       nchunk = max(nchunk,1)
       !write(6,*) i,subset_el(i),nchunk,nthreads

!$OMP TASKWAIT 
!!!$OMP TASKGROUP 
       do jel=1,subset_el(i),nchunk

!$OMP TASK &
!$OMP& DEFAULT(NONE) &
!$OMP& FIRSTPRIVATE(jel,i) &
!$OMP& PRIVATE(j,ie) &
!$OMP& SHARED(curr_stage,nlvddi,nlev,itvd,itvdv,istot,isact,nchunk) &
!$OMP& SHARED(difmol,robs,wsink,rload,ddt,rkpar) &
!$OMP& SHARED(rso,rsn,rsot,rsnt,dt,nkn) &
!$OMP& SHARED(cn,cdiag,clow,chigh,subset_el,cc,co) &
!$OMP& SHARED(subset_num,indipendent_subset) &
!$OMP& SHARED(difhv,cbound,gradxv,gradyv,cobs,rtauv,load,difv,wsinkv) &
!$OMP& SHARED(coeff_erk,coeff_srk,coeff_crk)

       do j=jel,jel+nchunk-1 	! loop over elements in subset
		if(j .le. subset_el(i)) then
	        ie = indipendent_subset(j,i)
	        !print *,i,ie
                call conz3d_element( &
     &			 curr_stage &
     &			,coeff_erk,coeff_srk,coeff_crk &
     &                  ,ie,cdiag,clow,chigh,cn,cc,co &
     &			,dt &
     &                  ,rkpar,difhv,difv &
     &			,difmol,cbound &
     &		 	,itvd,itvdv,gradxv,gradyv &
     &			,cobs,robs,rtauv &
     &			,wsink,wsinkv &
     &			,rload,load &
     &			,rso,rsn,rsot,rsnt &
     &			,nlvddi,nlev &
     &			,erk_reg,crk_reg,srk_reg)
		end if
	end do ! end loop over el in subset
!$OMP END TASK
      end do

!!!$OMP END TASKGROUP       
!$OMP TASKWAIT       

       end do ! end loop over subset
       
       if( shympi_partition_on_elements() ) then
         !call shympi_comment('shympi_elem: exchange scalar')
         call shympi_exchange_and_sum_3d_nodes(cn)
         call shympi_exchange_and_sum_3d_nodes(cdiag)
         call shympi_exchange_and_sum_3d_nodes(clow)
         call shympi_exchange_and_sum_3d_nodes(chigh)
       end if

       ntot = nkn
       if( shympi_partition_on_nodes() ) ntot = nkn_unique
!$     nchunk = ntot / ( nthreads * 10 )
       nchunk = max(nchunk,1)

!$OMP TASKWAIT
!!!$OMP TASKGROUP
       do knod=1,ntot,nchunk
!$OMP TASK FIRSTPRIVATE(knod) PRIVATE(k) DEFAULT(NONE)      &
!$OMP& SHARED(cn,cdiag,clow,chigh,cc,cbound,load,nchunk,   &
!$OMP&           rload,dt,nlvddi,ntot)
	 do k=knod,knod+nchunk-1
	 if(k .le. ntot) then
	   call conz3d_nodes(&
     &			curr_stage, &
     &			coeff_erk,coeff_srk,coeff_crk, &
     &			k,cn,cdiag(:,k),clow(:,k),chigh(:,k), &
     &                  cc,cbound,load,rload, &
     &                  is_rk_explicit,dt,nlvddi, &
     &			erk_reg,crk_reg,srk_reg)
         endif
         enddo
!$OMP END TASK 	      
	end do

!!!$OMP END TASKGROUP
!$OMP TASKWAIT       

	cc = real(cn)		!here happens INTEL_BUG
	
	if (bdebug ) then
	  iudb = 990 + my_id
	  ks = 2314
	  k = ipint(ks)
	  call get_act_dtime(dtime)
	  if( k > 0 .and. dtime == 1500. ) then
	    lmax = ilhkv(k)
	    write(iudb,*) 'after: ',dtime,ipext(k)
	    do l=1,lmax
	      write(iudb,*) l,cn(l,k),cc(l,k)
	    end do
	  end if
	end if

	DEALLOCATE(cn)
	DEALLOCATE(cdiag)
	DEALLOCATE(clow)
	DEALLOCATE(chigh)
	
!	call cpu_time(time2)
!!$	dtime2 = omp_get_wtime()
!	write(6,*) time2-time1,dtime2-dtime1

!----------------------------------------------------------------
! end of routine
!----------------------------------------------------------------

	end

!*****************************************************************

       subroutine conz3d_element(  &
     &			curr_stage &
     &			,coeff_erk,coeff_srk,coeff_crk &
     &			,ie &
     &			,cdiag,clow,chigh,cn,cc,co &
     &			,dt &
     &                  ,rkpar,difhv,difv &
     &			,difmol,cbound &
     &		 	,itvd,itvdv,gradxv,gradyv &
     &			,cobs,robs,rtauv &
     &			,wsink,wsinkv &
     &			,rload,load &
     &			,rso,rsn,rsot,rsnt &
     &			,nlvddi,nlev &
     &			,erk_reg,crk_reg,srk_reg)

	use mod_rungekutta, only : n_rkstages, &
     &	                           urk_reg,vrk_reg,a_erk,a_irk
        use mod_bound_geom
	use mod_geom
	use mod_depth
	use mod_diff_aux
	use mod_bound_dynamic
	use mod_area
	use mod_ts
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin
	use mod_layer_thickness
      
      implicit none
      
      integer,intent(in) :: curr_stage,ie,nlvddi,nlev,itvd,itvdv
      real,intent(in) :: difmol,robs,wsink,rload,rkpar
      real,dimension(nlvddi,nkn),intent(in) :: cc,co,cbound
      real,dimension(nlvddi,nel),intent(in) :: difhv
      real,dimension(nlvddi,nkn),intent(in) :: gradxv,gradyv
      real,dimension(nlvddi,nkn),intent(in) :: cobs,rtauv,load
      real,intent(in),dimension(0:nlvddi,nkn) :: wsinkv,difv
      real,dimension(n_rkstages),intent(in) :: coeff_erk
      real,dimension(n_rkstages+1),intent(in) :: coeff_srk
      real,dimension(n_rkstages+1),intent(in) :: coeff_crk
      double precision,intent(in) :: dt
      double precision,intent(in) :: rso,rsn,rsot,rsnt
      double precision,dimension(nlvddi,nkn),intent(inout) :: cdiag
      double precision,dimension(nlvddi,nkn),intent(inout) :: clow
      double precision,dimension(nlvddi,nkn),intent(inout) :: chigh
      double precision,dimension(nlvddi,nkn),intent(inout) :: cn
      real,dimension(nlvdi,nkn,n_rkstages-1),intent(inout) :: erk_reg
      real,dimension(nlvdi,nkn,n_rkstages-1),intent(inout) :: crk_reg
      real,dimension(nlvdi,nkn,n_rkstages-1),intent(inout) :: srk_reg

      logical :: btvdv,btvd
      integer :: k,ii,l,iii,ll,ibase,lstart,ilevel,itot,isum
      integer :: jlevel
      integer :: n,i,iext
      integer :: istage,jstage,mstage
      integer, dimension(3) :: kn
      real :: as_ll,ac_ll,c_l,inv_ae_ll,ai_ll,ai_llm
      double precision :: cexpl,cbm,ccm,waux,loading,wws
      double precision :: aj,rk3,aj4,aj12
      double precision :: hmed,hmbot,hmtop,hmotop,hmobot
      double precision :: hmntop,hmnbot,rvptop,rvpbot,w,aux
      double precision :: flux_tot,flux_tot1,flux_top,flux_bot
      double precision :: rstot,hn,ho,hc,cdummy,alow,adiag,ahigh,rrc
      double precision :: rkmin,rkmax,cconz
      double precision :: rhs_us,rhs_vs,sum_us,sum_vs
      double precision,dimension(curr_stage) :: us,vs
      double precision,dimension(3) :: fw,fd,fl,fnudge_o,fnudge_c
      double precision,dimension(3) :: b,c,f,wdiff
      double precision,dimension(0:nlvddi+1) :: haver,presentl
      double precision,dimension(0:nlvddi+1,3) :: hnew,rtau,cob
      double precision,dimension(0:nlvddi+1,3) :: hold,hcur,vflux,wl
      double precision,dimension(0:nlvddi+1,3) :: cl
      double precision,dimension(0:nlvddi+1,3) :: finu
      double precision,dimension(nlvddi,3) :: clc,clm,clp,cle
	
	if(nlv.ne.nlev) stop 'error stop conz3d_element: nlv/=nlev'

! ----------------------------------------------------------------
!  initialize variables and parameters
! ----------------------------------------------------------------

	btvd = itvd .gt. 0	!flags for tvd scheme
	btvdv = itvdv .gt. 0

!  renaming of special runge-kutta coefficients in Butcher tableaux
        ac_ll = coeff_crk(curr_stage+1)	!diagonal coeff for vertical adv tableau
        as_ll = coeff_srk(curr_stage+1)	!diagonal coeff for stiffly implicit tableau
	c_l = sum(coeff_erk)		!c coeff

! ----------------------------------------------------------------
! global arrays for accumulation of implicit terms
! ----------------------------------------------------------------

! 	 ALLOCATE(fw(3),fd(3),fl(3),fnudge(3),wdiff(3))
! 	 ALLOCATE(b(3),c(3),f(3))
! 	 ALLOCATE(presentl(0:nlvddi+1))
! 	 ALLOCATE(hnew(0:nlvddi+1,3),htnew(0:nlvddi+1,3))
! 	 ALLOCATE(rtau(0:nlvddi+1,3),cob(0:nlvddi+1,3))
! 	 ALLOCATE(vflux(0:nlvddi+1,3),wl(0:nlvddi+1,3))
! 	 ALLOCATE(cl(0:nlvddi+1,3))
! 	 ALLOCATE(clc(nlvddi,3),clm(nlvddi,3))
! 	 ALLOCATE(clp(nlvddi,3),cle(nlvddi,3))

          haver = 0.
	  presentl = 0.		!1. if layer is present
	  hnew = 0.		!as hreal but with zeta_new
	  hold = 0.		!as hreal but with zeta_old
	  hcur = 0.		!as hreal but with zeta_cur
	  cl = 0.		!concentration in layer
	  wl = 0.		!vertical velocity
	  vflux = 0.		!vertical flux
	
!	these are the local arrays for accumulation of implicit terms
!	(maybe we do not need them, but just to be sure...)
!	after accumulation we copy them onto the global arrays

	    cle = 0.
	    clc = 0.
	    clm = 0.
	    clp = 0.
      
	do ii=1,3
          k=nen3v(ii,ie)
	  kn(ii)=k
	  b(ii)=ev(ii+3,ie)
	  c(ii)=ev(ii+6,ie)
	end do

	aj=ev(10,ie)    !area of triangle / 12
	aj4=4.*aj
	aj12=12.*aj
        ilevel=ilhv(ie)
	jlevel=jlhv(ie)

! 	----------------------------------------------------------------
! 	set up vectors for use in assembling contributions
! 	----------------------------------------------------------------
! 
!	note that hdeov is at current stage, hdenv is at the new stage
!	hdkov is at old stage, hdknv is at new stage, hdkcv is at
!	current stage

        do l=jlevel,ilevel
          haver(l) = rso*hdenv(l,ie) + rsot*hdeov(l,ie)
	  presentl(l) = 1.
	  do ii=1,3
	    k=kn(ii)
	    hn = hdknv(l,k)		! there are never more layers in ie
	    ho = hdkov(l,k)		! ... than in k
	    hc = hdkcv(l,k)
	    hold(l,ii) = rso * hn + rsot * ho
	    hnew(l,ii) = rsn * hn + rsnt * ho
	    hcur(l,ii) = rso * hn + rsot * hc
	    cl(l,ii) = cc(l,k)
	    cob(l,ii) = cobs(l,k)	!observations
	    rtau(l,ii) = rtauv(l,k)	!observations
	    wl(l,ii) = wlnv(l,k) - wsink * wsinkv(l,k)
	  end do
	end do

	do l=ilevel+1,nlv
	  presentl(l) = 0.
	end do

! 	----------------------------------------------------------------
! 	set vertical velocities in surface and bottom layer
! 	----------------------------------------------------------------
! 
! 	we do not set wl(0,ii) because otherwise we loose concentration
! 	through surface
! 
! 	we set wl(ilevel,ii) to 0 because we are on the bottom
! 	and there should be no contribution from this element
! 	to the vertical velocity

	do ii=1,3
	  wl(ilevel,ii) = 0.
	end do

! 	----------------------------------------------------------------
! 	compute vertical fluxes (w/o vertical TVD scheme)
! 	----------------------------------------------------------------

	wws = 0.	!sinking already in wl
	call vertical_flux_ie(btvdv,ie,ilevel,jlevel, &
     &			      dt,wws,cl,wl,hcur,vflux)

! ----------------------------------------------------------------
!  loop over levels
! ----------------------------------------------------------------

        do l=jlevel,ilevel

! 	----------------------------------------------------------------
! 	compute advection transport
! 	----------------------------------------------------------------
! 
!	The discharge values us,vs used in the tracer equation are
!	computed by solving the following linear system, whose
!	coefficient matrix is a lower triangular submatrix of the
!	explicit weight matrix A. At each stage l we compute us,vs
!	with forward substitution:

        do istage=1,curr_stage
	  inv_ae_ll = 1./a_erk(istage,istage)
	  ai_ll = a_irk(istage,curr_stage+1)
	  ai_llm = a_irk(istage,curr_stage)
	  rhs_us = ai_ll*utlnv(l,ie) + ai_llm*utlcv(l,ie)
	  rhs_vs = ai_ll*vtlnv(l,ie) + ai_llm*vtlcv(l,ie)
	  do jstage=1,curr_stage-1
	    rhs_us = rhs_us + a_irk(istage,jstage)*urk_reg(l,ie,jstage)
	    rhs_vs = rhs_vs + a_irk(istage,jstage)*vrk_reg(l,ie,jstage)
	  end do
          sum_us = 0.0d0
          sum_vs = 0.0d0
          do mstage=1,istage-1
            sum_us = sum_us + a_erk(istage,mstage) * us(mstage)
            sum_vs = sum_vs + a_erk(istage,mstage) * vs(mstage)
          end do
          us(istage) = (rhs_us - sum_us) * inv_ae_ll
          vs(istage) = (rhs_vs - sum_vs) * inv_ae_ll
        end do


        rk3 = 3. * rkpar * difhv(l,ie)

	cbm=0.
	ccm=0.
	itot=0
	isum=0
	do ii=1,3
	  k=kn(ii)
	  f(ii)=us(curr_stage)*b(ii)+vs(curr_stage)*c(ii)	!$$azpar
	  if(f(ii).lt.0.) then	!flux out of node
	    itot=itot+1
	    isum=isum+ii
	  end if
	  cbm=cbm+b(ii)*cl(l,ii)
	  ccm=ccm+c(ii)*cl(l,ii)

! 	  ----------------------------------------------------------------
! 	  initialization to be sure we are in a clean state
! 	  ----------------------------------------------------------------

	  fw(ii) = 0.
	  !cle(l,ii) = 0.	!ERIC
	  !clc(l,ii) = 0.
	  !clm(l,ii) = 0.
	  !clp(l,ii) = 0.

! 	  ----------------------------------------------------------------
! 	  contributions from horizontal diffusion
! 	  ----------------------------------------------------------------

          waux = 0.
          do iii=1,3
            waux = waux + wdifhv(iii,ii,ie) * cl(l,iii)
          end do
          wdiff(ii) = waux

! 	  ----------------------------------------------------------------
! 	  contributions from vertical diffusion
! 	  ----------------------------------------------------------------
! 
! 	  in fd(ii) is explicit contribution
! 	  the sign is for the term on the left side, therefore
! 	  fd(ii) must be subtracted from the right side
! 
! 	  maybe we should use real layer thickness, or even the
! 	  time dependent layer thickness

	  rvptop = difv(l-1,k) + difmol
	  rvpbot = difv(l,k) + difmol
	  hmotop =2.*rvptop*presentl(l-1)/(hcur(l-1,ii)+hcur(l,ii))
	  hmobot =2.*rvpbot*presentl(l+1)/(hcur(l,ii)+hcur(l+1,ii))
	  hmntop =2.*rvptop*presentl(l-1)/(hnew(l-1,ii)+hnew(l,ii))
	  hmnbot =2.*rvpbot*presentl(l+1)/(hnew(l,ii)+hnew(l+1,ii))

	  fd(ii) = (cl(l,ii)-cl(l+1,ii))*hmobot - &
     &		   (cl(l-1,ii)-cl(l,ii))*hmotop

	  clc(l,ii) = clc(l,ii) + as_ll * ( hmntop + hmnbot )
	  clm(l,ii) = clm(l,ii) - as_ll * ( hmntop )
	  clp(l,ii) = clp(l,ii) - as_ll * ( hmnbot )

! 	  ----------------------------------------------------------------
! 	  contributions from vertical advection
! 	  ----------------------------------------------------------------
! 
! 	  in fw(ii) is explicit contribution
! 	  the sign is for the term on the left side, therefore
! 	  fw(ii) must be subtracted from the right side
! 
! 	  if we are in last layer, w(l,ii) is zero
! 	  if we are in first layer, w(l-1,ii) is zero (see above)

	  w = wl(l-1,ii)		!top of layer
	  if( l .eq. jlevel ) w = 0.	!surface -> no transport (WZERO)
	  if( w .ge. 0. ) then
	    clc(l,ii) = clc(l,ii) + ac_ll*w
	  else
	    clm(l,ii) = clm(l,ii) + ac_ll*w
	  end if

	  w = wl(l,ii)			!bottom of layer
	  if( l .eq. ilevel ) w = 0.	!bottom -> handle flux elsewhere (WZERO)
	  if( w .gt. 0. ) then
	    clp(l,ii) = clp(l,ii) - ac_ll*w
	  else
	    clc(l,ii) = clc(l,ii) - ac_ll*w
	  end if

	  flux_tot = vflux(l-1,ii) - vflux(l,ii)

	  fw(ii) = flux_tot
	end do

! 	----------------------------------------------------------------
! 	contributions from horizontal advection (only explicit)
! 	----------------------------------------------------------------
! 
! 	f(ii) > 0 ==> flux into node ii
! 	itot=1 -> flux out of one node
! 		compute flux with concentration of this node
! 	itot=2 -> flux into one node
! 		for flux use conz. of the other two nodes and
! 		minus the sum of these nodes for the flux of this node

	if(itot.eq.1) then	!$$flux
	  fl(1)=f(1)*cl(l,isum)
	  fl(2)=f(2)*cl(l,isum)
	  fl(3)=f(3)*cl(l,isum)
	else if(itot.eq.2) then
	  isum=6-isum
	  fl(1)=f(1)*cl(l,1)
	  fl(2)=f(2)*cl(l,2)
	  fl(3)=f(3)*cl(l,3)
	  fl(isum) = 0.
	  fl(isum) = -(fl(1)+fl(2)+fl(3))
	  isum=6-isum		!reset to original value
	else			!exception	$$itot0
	  fl(1)=0.
	  fl(2)=0.
	  fl(3)=0.
	end if

! 	----------------------------------------------------------------
! 	horizontal TVD scheme start - compute fluxes fl, otherwise leave as is
! 	----------------------------------------------------------------

        if( btvd ) then
	  iext = 0
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( is_external_boundary(k) ) iext = iext + 1
	  end do

          if( iext .eq. 0 ) then
	    call tvd_fluxes(ie,l,itot,isum,dt,cl,cc,gradxv,gradyv,f,fl)
	  end if
	end if

! 	----------------------------------------------------------------
! 	horizontal TVD scheme finish
! 	----------------------------------------------------------------

! 	----------------------------------------------------------------
! 	contributions from nudging
! 	----------------------------------------------------------------

	do ii=1,3
	  fnudge_o(ii) = robs * rtau(l,ii) * cob(l,ii)  !explicit contributions
	  fnudge_c(ii) = -robs * rtau(l,ii) * cl(l,ii)  ! 
	  finu(l,ii) = robs * rtau(l,ii)		!implicit contribution:
	end do						!later set to zero

!	------------------------------------------------------
!	Set up current stage right hand side for conc.
!       F^c_l = R_l+W_l+D_l
!	rrc is explicit contribution R_l
!	fw  is vertical advection contribution W_l
!	fd  is stiffly implicit contribution D_l
!	Set up concentration matrix A^c
!	------------------------------------------------------

	do ii=1,3
	  k=kn(ii)
          hmed = haver(l)                    !new ggu   !HACK

          rrc = 3.*fl(ii) - rk3*hmed*wdiff(ii) + fnudge_c(ii)

	  cexpl = aj4 * ( &
	            hold(l,ii)*co(l,k) &
     &	              + dt *  (   c_l*hold(l,ii)*fnudge_o(ii) &
     &		                + coeff_erk(curr_stage)*rrc &
     &                          - coeff_crk(curr_stage)*fw(ii) &
     &		                - coeff_srk(curr_stage)*fd(ii) ) &
     &                  )
	  
	  !clm(1,ii) = 0.		!ERIC
	  !clp(ilevel,ii) = 0.
	  ! next check to be deleted
	  if( clm(1,ii) /= 0. .or. clp(ilevel,ii) /= 0. ) then
	    write(6,*) ie,ii,ilevel
	    write(6,*) clm(1,ii),clp(ilevel,ii)
	    stop 'error stop: assumption violated'
	  end if
	  
	  alow  = aj4 * dt * clm(l,ii)
	  ahigh = aj4 * dt * clp(l,ii)
	  adiag = aj4 * dt * clc(l,ii)  &
     &			+ aj4 * (1.+ 0.*dt*finu(l,ii)) * hnew(l,ii)
	  cn(l,k)    = cn(l,k)    + cexpl
	  clow(l,k)  = clow(l,k)  + alow
	  chigh(l,k) = chigh(l,k) + ahigh   
          cdiag(l,k) = cdiag(l,k) + adiag

!	------------------------------------------------------
!	save current stage right hand side R_l, W_l, D_l
!	------------------------------------------------------

          if (curr_stage .ne. n_rkstages) then !if (not last stage)
	    erk_reg(l,k,curr_stage) = erk_reg(l,k,curr_stage) &
     &        + aj4 * rrc
	    crk_reg(l,k,curr_stage) = crk_reg(l,k,curr_stage) &
     &        - aj4 * fw(ii)
	    srk_reg(l,k,curr_stage) = srk_reg(l,k,curr_stage) &
     &        - aj4 * fd(ii)
	  end if

	end do

	end do		! loop over l

! ----------------------------------------------------------------
!  end of loop over l
! ----------------------------------------------------------------

! 	deallocate(fw,fd,fl,fnudge)
! 	deallocate(b,c,f,wdiff)
! 	deallocate(haver,presentl)
! 	deallocate(hnew,rtau,cob)
! 	deallocate(hold,vflux,wl,cl)
! 	deallocate(clc,clm,clp,cle)
! 	
! ----------------------------------------------------------------
!  end of routine
! ----------------------------------------------------------------

      end subroutine conz3d_element

! *****************************************************************
      
       subroutine conz3d_nodes(		   &
     &			       curr_stage, &
     &			       coeff_erk,coeff_srk,coeff_crk, &
     &			       k, &
     &			       cn,cdiag,clow,chigh,cc,cbound, &
     &                         load,rload,is_explicit,dt,nlvddi, &
     &			       erk_reg,crk_reg,srk_reg)

	use mod_rungekutta, only : n_rkstages
      	use mod_bound_geom
	use mod_geom
	use mod_depth
	use mod_diff_aux
	use mod_bound_dynamic
	use mod_area
	use mod_ts
	use mod_hydro_vel
	use mod_hydro
	use evgeom
	use levels
	use basin
	use shympi
	
	implicit none
	
	integer,intent(in) :: curr_stage,k,nlvddi
        real,dimension(n_rkstages),intent(in) :: coeff_erk
        real,dimension(n_rkstages+1),intent(in) :: coeff_srk
        real,dimension(n_rkstages+1),intent(in) :: coeff_crk

	real,intent(in) :: rload
	real,dimension(nlvddi,nkn),intent(in) :: cc,cbound
	real,dimension(nlvddi,nkn),intent(inout) :: load 		!LLL
	double precision, intent(in) :: dt

	double precision,dimension(nlvddi,nkn),intent(inout) :: cn
	double precision,dimension(nlvddi),intent(inout) :: cdiag
	double precision,dimension(nlvddi),intent(inout) :: clow
	double precision,dimension(nlvddi),intent(inout) :: chigh

        real,dimension(nlvdi,nkn,n_rkstages-1),intent(in) :: erk_reg
        real,dimension(nlvdi,nkn,n_rkstages-1),intent(in) :: crk_reg
        real,dimension(nlvdi,nkn,n_rkstages-1),intent(in) :: srk_reg

	logical, intent(in) :: is_explicit

	logical :: bdry
	integer :: l,ilevel,jlevel,jstage,lstart,i,ii,ie,n,ibase
	double precision :: mflux,qflux,cconz
	double precision :: loading,aux,cload

	double precision, parameter :: d_tiny = tiny(1.d+0)
	double precision, parameter :: r_tiny = tiny(1.)
      
	logical, save :: bdebug = .false.
	character*80 aline

	logical :: is_dry_node

! ----------------------------------------------------------------
!  debug code
! ----------------------------------------------------------------

! ----------------------------------------------------------------
!  handle boundary (flux) conditions and
!  set up previous stage right hand side F^c_l = R_l+W_l+D_l
! ----------------------------------------------------------------

	  bdry = is_dry_node(k)

      	  ilevel = ilhkv(k)
          jlevel = jlhkv(k)

	  do l=jlevel,ilevel

            !mflux = cbound(l,k)		!mass flux has been passed
	    cconz = cbound(l,k)			!concentration has been passed
	    qflux = mfluxv(l,k)
	    if( bdry ) qflux = 0.
	    if( qflux .lt. 0. .and. is_boundary(k) ) cconz = cc(l,k)
	    mflux = qflux * cconz

            cn(l,k) = cn(l,k) + dt * mflux	!explicit treatment

	    loading = rload*load(l,k)
            if( loading == 0 ) then			!no loading
              !nothing
            else if ( loading < 0.d0 ) then		!excess deposition
	      cload = 0.
	      if( cn(l,k) > 0. ) then
                cload = - dt * loading
                cload = cn(l,k) * ( 1. - exp(-cload/cn(l,k)) )
	      else if( cn(l,k) < 0. ) then
	        cn(l,k) = 0.
	      end if
              if( cload > cn(l,k) ) goto 98
              loading = -cload / dt
              if( rload > 0. ) load(l,k) = loading / rload
              cn(l,k) = cn(l,k) + dt*loading
            else					!erosion
              cn(l,k) = cn(l,k) + dt*loading
            end if

	    do jstage=1,curr_stage-1 !lrp:dbg-imex
              cn(l,k) = cn(l,k) + dt*( &
     &          coeff_erk(jstage) * erk_reg(l,k,jstage) + &
     &          coeff_crk(jstage) * crk_reg(l,k,jstage) + &
     &          coeff_srk(jstage) * srk_reg(l,k,jstage) )
            end do

	    if( bdebug ) then
	     if( mflux /= 0. .or. loading /= 0. ) then
	      call get_act_timeline(aline)
	      write(666,*) trim(aline),l,k
	      write(666,*) cconz,qflux,mflux
	      write(666,*) rload,loading
	     end if
	    end if

	  end do

! ----------------------------------------------------------------
!  compute concentration for each node (solve system) A^c C = F^c
! ----------------------------------------------------------------

	if(is_explicit .or. (nlv .eq. 1)) then

	  if( nlv .gt. 1 ) then
	    write(6,*) 'conz: computing explicitly ',nlv
	  end if

	  ilevel = ilhkv(k)
	  do l=jlevel,ilevel
	    if(cdiag(l).ne.0.) then
	      cn(l,k)=cn(l,k)/cdiag(l)
	    end if
	  end do

	else

	  ilevel = ilhkv(k)
	  aux=1./cdiag(jlevel)
	  chigh(jlevel)=chigh(jlevel)*aux
	  cn(jlevel,k)=cn(jlevel,k)*aux
	  do l=jlevel+1,ilevel
	    if( cdiag(l) == 0. ) goto 99
	    aux=1./(cdiag(l)-clow(l)*chigh(l-1))
	    chigh(l)=chigh(l)*aux
	    cn(l,k)=(cn(l,k)-clow(l)*cn(l-1,k))*aux
	  end do
	  lstart = ilevel-1
	  do l=lstart,jlevel,-1	!$$LEV0 bug 14.08.1998 -> ran to 0
	    cn(l,k)=cn(l,k)-cn(l+1,k)*chigh(l)
	  end do
	end if
	
! ----------------------------------------------------------------
!  end of routine
! ----------------------------------------------------------------

	return
   98	continue
	write(6,*) 'error computing loading: ',l,k
	write(6,*) 'loading,cload,cn(l,k): ',loading,cload,cn(l,k)
	stop 'error stop conz3d_nodes: internal error (1)'
   99	continue
	write(6,*) k,l,ilevel
	write(6,*) nkn_inner,nkn_local,nkn
	write(6,*) cdiag(l),clow(l),chigh(l-1)
	stop 'error stop conz3d_nodes (omp): diag == 0'
      end subroutine conz3d_nodes

!*****************************************************************

        subroutine tsdebug(robs,cobs,rtauv)

        use basin
        use levels

        implicit none

        real robs
        real cobs(nlvdi,nkn)
        real rtauv(nlvdi,nkn)

        integer k,lmax,i
        integer ipext,ipint

        integer, parameter :: nmax=5
        integer :: nodes(nmax) = (/37806,52065,80988,65554,53926/)

        do i=1,nmax
          k = nodes(i)
          k = ipint(k)
          lmax= ilhkv(k)
          call tsdebugk(k,lmax,robs,cobs(1:lmax,k),rtauv(1:lmax,k))
        end do

        end

!*****************************************************************

        subroutine tsdebugk(k,lmax,robs,cobs,rtau)

        integer k,lmax
        real robs
        real cobs(lmax)
        real rtau(lmax)

        integer iunit,ke
        real aux(lmax)
        integer ipext,ipint
        
        iunit = 789
        aux = 0
        where ( rtau > 0 ) aux = 1./rtau
        ke = ipext(k)

        write(iunit,*) '------------------------------------'
        write(iunit,*) k,ke,lmax,robs
        write(iunit,*) minval(cobs),maxval(cobs)
        write(iunit,*) minval(rtau),maxval(rtau)
        write(iunit,*) minval(aux),maxval(aux)
        write(iunit,*) '------------------------------------'

        end

!*****************************************************************

