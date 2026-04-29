
!--------------------------------------------------------------------------
!
!    Copyright (C) 2015,2017-2019  Christian Ferrarin
!    Copyright (C) 2019  Georg Umgiesser
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

! subroutines for computing tidal potential and tidal analysis
!
!********************************************************************
!
! revision log :
!
! 30.09.2015	ccf	converted previous routines to module
! 01.10.2015	ccf	bug fix in hour, now handle negative time
! 21.10.2015	ccf	documentation included
! 29.03.2017	ccf	bug fix in the astronomical arguments chi
! 29.04.2018	ccf	new version of tidal potential and analysis
! 18.02.2019	ccf	inclued more constituents
! 13.03.2019	ggu	changed VERS_7_5_61
! 21.05.2019	ggu	changed VERS_7_5_62
! 28.04.2023	ggu	avoid using floating for **2
! 13.10.2024	ggu	bug fix for INTEL_BUG
! 15.10.2024	ggu	bug fix for INTEL_BUG (IFX)
! 13.11.2024    ggu     marked old code with INTEL_BUG_OLD
!
!********************************************************************
! DOCS  START   S_tide

! SHYFEM includes an as\-tro\-no\-mi\-cal tidal model which can be
! activated by setting the parameter |rtide| equal 1 in the |para| 
! section.
!
! The model calculates equilibrium tidal potential ($\eta$) and load 
! tides ($\beta$) and uses these to force the free surface. 
! The term $\eta$ in the momentum equations is calculated as a sum
! of the tidal potential of each tidal constituents multiplied by the
! frequency-dependent elasticity factor. The factor $\beta$ accounts 
! for the effect of the load tides, assuming that loading tides are
! in-phase with the oceanic tide. $\beta$ is function of the water 
! depth as $\beta=ltidec*H$ with |ltidec| a calibration factor to be 
! set in the str |para| section.
!
! The model cosiders the following tidal costituents:
! \begin{itemize}
! \item Semidiurnal species:
!    \begin{itemize}
!    \item M2  semi-diurnal principal lunar
!    \item S2  semi-diurnal principal solar
!    \item N2  large elliptical tide of first-order to M2
!    \item K2  semi-diurnal declination to M2
!    \item NU2 large evection tide to M2
!    \item MU2 large variation tide to M2
!    \item L2  small elliptical tide of first-order to M2
!    \item T2  large elliptical tide of first-order to S2
!    \end{itemize}
! \item Diurnal species:
!    \begin{itemize}
!    \item K1  declination luni-solar
!    \item O1  principal lunar
!    \item P1  principal solar
!    \item Q1  elliptical lunar
!    \item J1  elliptical tide of first-order to K1
!    \item OO1 evection tide to O1
!    \item S1  radiational tide 
!    \end{itemize}
! \item Long-Period species:
!    \begin{itemize}
!    \item MF  fortnightly lunar
!    \item MM  monthly lunar
!    \item MSM S0-semiannual solar
!    \item MSF Evection tide to M0
!    \item SSA semiannual solar
!    \item SA  elliptical tide of first-order to S0
!    \end{itemize}
! \end{itemize}
!
! SHYFEM also allows to perform the tidal analysis of water levels
! during the model runtime. The tidal analysis is actived by setting
! |itmtid| and |idttid|. |idttid| should be long enough to perform 
! a reliable analysis. The parameter |itmtid| can be used to start 
! the analysis after the simulation spin-up. The tidal analysis module
! write an output file .tide.shy containing amplitudes and phases
! of all tidal constituents over the computational domain (on the 
! nodes).
!
! DOCS  END
!********************************************************************
! Initialize tidal variables and parameters

        subroutine tidefini

        use tide
        use basin
        use coordinates, only : iproj
	use mod_info_output

        implicit none

        integer  	 :: isphe      !if = 1  coordinates are spherical
        double precision :: dgetpar

        zeqv = 0.
	vvk  = 0.d0
	uvk  = 0.d0
	fvk  = 0.d0
	astr = 0.d0
	ader = 0.d0

        rtide = dgetpar('rtide')

        if( rtide <= 0. ) return

        ltidec = dgetpar('ltidec')

        call get_coords_ev(isphe)

        if( isphe <= 0 .AND. iproj <= 0 ) then
            write(6,*) 'isphe,iproj: ',isphe,iproj
            write(6,*) 'for tidal potential either '
            write(6,*) 'isphe or iproj must be set'
            stop 'error stop tidefini: no lat/lon coordinates'
        end if

	if( print_not_quiet_once() ) then
          write(6,*) 'Tidal potential activated'
	end if

        end subroutine tidefini

!*********************************************************************
! This subroutine computes tidal potential on nodes

        subroutine tideforce

        use tide
        use basin, only : nkn
        use mod_depth
!$      use omp_lib

        implicit none

        integer 		:: k,iks,ikend,nchunk,nthreads

!--------------------------------------------------------
!------ compute tidal potential? ------------------------
!--------------------------------------------------------

        if( rtide <= 0 ) return

!--------------------------------------------------------
!------ computes eq. tide and load tide -----------------
!--------------------------------------------------------

!$OMP PARALLEL 
!$OMP SINGLE

        call omp_compute_chunk(nkn,nchunk)
        nthreads = 1
!$      nthreads = omp_get_num_threads()

        do k = 1,nkn,nchunk

!$OMP TASK FIRSTPRIVATE(k) PRIVATE(iks,ikend) &
!$OMP&     SHARED(nkn,nchunk) DEFAULT(NONE)

          call omp_compute_minmax(nchunk,nkn,k,ikend)

          do iks=k,ikend
            call equilibrium_tide(iks)
          end do

!$OMP END TASK

        end do

!$OMP END SINGLE
!$OMP TASKWAIT  
!$OMP END PARALLEL      

!--------------------------------------------------------
!------ end of routine ----------------------------------
!--------------------------------------------------------

        end subroutine tideforce

!*********************************************************************
! This subroutine computes the equilibrium tide as a sum of the single
! costituents using nodal modulation. Account also for the loading tide

        subroutine equilibrium_tide(k)

        use tide
        use mod_depth
        use coordinates
        use mod_hydro
        use basin, only : nkn
        use iso8601
	use shympi

        implicit none

	integer, intent(in)	 :: k

        real, allocatable        :: latf(:)
        real  			 :: llat,llon
	double precision	 :: fact,arg 			!INTEL_BUG
	double precision	 :: hdepth			!INTEL_BUG
	!real			 :: fact,arg 			!INTEL_BUG_OLD
	!real			 :: hdepth			!INTEL_BUG_OLD
        real, parameter          :: d2r=4.0*atan(1.0)/180.0
        real, parameter          :: twopi=8.d0*atan(1.d0)
	integer			 :: i,ld1
        double precision         :: eqtide  !Equilibrium tide [m] !INTEL_BUG
        !real		         :: eqtide  !Equilibrium tide [m] !INTEL_BUG_OLD
        real			 :: loadb   !loading tide factor 
					    ! [0.054,0.047,= ltidec*depth]

	integer iudb,ks
	logical bdebug
	double precision dtime
	integer ipext
	!iudb = 890 + my_id
	!ks = 90521
	!call get_act_dtime(dtime)
	!bdebug = ipext(k) == ks .and. dtime == 400.
	bdebug = .false.

!-----------------------------------------------
!       Set up constants and parameters
!-----------------------------------------------

	allocate(latf(0:2))

        llon  = xgeov(k)
        if (llon < 0.0) llon = 360. + llon
        llon = llon * d2r
        llat = ygeov(k)  * d2r

        latf(2) = cos(llat)**2         !semi-diurnal
        latf(1) = sin(2.*llat)         !diurnal
        latf(0) = 1.5*latf(2) - 1.     !long-period

!-----------------------------------------------
!       Compute equilibrium tide
!-----------------------------------------------

        eqtide = 0.
        do i = 1,ntd
          ld1  = const_ar(i)%dood(1)
          fact = const_ar(i)%mask * const_ar(i)%loven * &
     &           fvk(i,k) * const_ar(i)%amp * latf(ld1)
          arg = (vvk(i,k)+uvk(i,k))*twopi + ld1*llon
          eqtide = eqtide + fact*cos(arg)
        end do

!-----------------------------------------------
!       Compute loading tide
!-----------------------------------------------

	hdepth = dble(hkv(k))+dble(zov(k))		!INTEL_BUG
        loadb = ltidec * hdepth				!INTEL_BUG
        !loadb = ltidec * ( hkv(k) + zov(k) )		!INTEL_BUG_OLD

!-----------------------------------------------
!	Compute earth + loading
!-----------------------------------------------

        zeqv(k) = eqtide + loadb*zov(k)

	if( bdebug ) then
	  write(iudb,*) dtime,k,ipext(k),ks
	  write(iudb,*) ltidec,hkv(k),zov(k)
	  write(iudb,*) hdepth
	  write(iudb,*) ltidec,loadb,eqtide,zeqv(k)
	end if

!-----------------------------------------------
!       End of routine
!-----------------------------------------------

        end subroutine equilibrium_tide

!*********************************************************************
! This subroutine calculates the V,U,F values (stored in the tide module) 
! at time d for all constituents ntd and for all nodes. This routine is 
! called by shyfem (do_befor) if rtide > 0 .and. itidana > 0

        subroutine tide_vuf

        use tide
        use basin
        use coordinates
!$      use omp_lib

        implicit none

        double precision              :: atime
        double precision              :: d    !Universal time (days since 1989-12-31::12)
        double precision              :: lat
        integer                       :: k
        integer         	      :: iks,ikend,nchunk,nthreads

        if( rtide <= 0 .and. itidana <= -1 ) return

        call get_absolute_act_time(atime)
        call tide_time(atime,d)
        call astro(d)

!$OMP PARALLEL 
!$OMP SINGLE

        call omp_compute_chunk(nkn,nchunk)
        nthreads = 1
!$      nthreads = omp_get_num_threads()

        do k = 1,nkn,nchunk

!$OMP TASK FIRSTPRIVATE(k) PRIVATE(iks,ikend,lat) &
!$OMP&     SHARED(nkn,nchunk,ygeov,vvk,uvk,fvk) &
!$OMP&     DEFAULT(NONE)

          call omp_compute_minmax(nchunk,nkn,k,ikend)

          do iks=k,ikend
            lat = ygeov(iks)
            call setvuf(lat,vvk(:,iks),uvk(:,iks),fvk(:,iks))
          end do

!$OMP END TASK

        end do

!$OMP END SINGLE
!$OMP TASKWAIT  
!$OMP END PARALLEL      

        end subroutine tide_vuf

!*****************************************************************
! Tidal analysis of water levels during simulation runtime. 
! It accumulates values during simulation runtine and
! performs tidal analysis at selected times.

        subroutine runtime_tidana

        use tide
        use basin
        use coordinates, only : iproj
        use mod_hydro, only : znv
!$      use omp_lib

        implicit none

        double precision, save  	    :: da_out(4) = 0.d0
        double precision, save, allocatable :: acov(:,:,:)
        double precision, save, allocatable :: bcov(:,:)
        double precision, save, allocatable :: tvar(:)
        !double precision, allocatable      :: x(:)		!ggu
        !double precision, allocatable 	    :: acovl(:,:)	!ggu
        double precision                    :: x(2*ntd+1)
        double precision 	            :: acovl(2*ntd+1,2*ntd+1)
        real, allocatable                   :: tideh(:,:)
        real, allocatable                   :: tideg(:,:)
        double precision, save  	    :: dtmtid
        double precision 		    :: dtime
	real				    :: dat
        integer, save			    :: nequ
	integer, save			    :: ndat
        integer                 	    :: isphe
        integer 			    :: nvar,id,ivar
	integer 			    :: i,i1,i2,k
        logical 			    :: has_output_d,next_output_d
        character*20 	   	   	    :: aline
        integer                             :: iks,ikend,nchunk,nthreads

!-----------------------------------------------------------------
! Start of code
!-----------------------------------------------------------------

        if( itidana <= -1 ) return

! ----------------------------------------------------------
! Initialization
! ----------------------------------------------------------

        if( itidana == 0 ) then

!         --------------------------------------------------------------
!         Initialize tidal parameters
!         --------------------------------------------------------------
          call init_output_d('itmtid','idttid',da_out)
          call increase_output_d(da_out)
          if( .not. has_output_d(da_out) ) itidana = -1
          if( itidana == -1 ) return

          call get_coords_ev(isphe)
          if( isphe <= 0 .AND. iproj <= 0 ) then
            stop 'error stop runtime_tidana: proj not defined'
	  end if

          nequ = ntd*2 + 1

          allocate(acov(nequ,nequ,nkn))
          allocate(bcov(nequ,nkn))
          allocate(tvar(nkn))

          acov = 0.d0
          tvar = 0.d0
          bcov = 0.d0
	  ndat = 0

!         -------------------------------------------------------------
!         Initialize output
!         -------------------------------------------------------------

          call convert_date_d('itmtid',dtmtid)

          nvar = count(const_ar%mask /= 0) * 2

          if( has_output_d(da_out) ) then
            call shyfem_init_scalar_file('tide',nvar,.true.,id)
            da_out(4) = id
          end if

          write(6,*) 'Runtime tidal analysis of water levels activated'

          itidana = 1

        endif

! -------------------------------------------------------------------
! Normal call, time loop, accumulate values if dtime > dtmtid
! -------------------------------------------------------------------

        call get_act_dtime(dtime)

        if( dtime < dtmtid ) return

        !nequ = ntd*2 + 1		!ggu
        !allocate(acovl(nequ,nequ))	!ggu
        !allocate(x(nequ))		!ggu

	ndat = ndat + 1

!$OMP PARALLEL 
!$OMP SINGLE

        call omp_compute_chunk(nkn,nchunk)
        nthreads = 1
!$      nthreads = omp_get_num_threads()

        do k = 1,nkn,nchunk

!$OMP TASK FIRSTPRIVATE(k) PRIVATE(iks,ikend,dat,x,i1,i2,acovl) &
!$OMP&     SHARED(nequ,nkn,nchunk,znv,vvk,uvk,fvk,tvar,bcov,acov) &
!$OMP&     DEFAULT(NONE)

          call omp_compute_minmax(nchunk,nkn,k,ikend)

          do iks=k,ikend
            dat = znv(iks)
            call get_overarray(nequ,vvk(:,iks),uvk(:,iks),fvk(:,iks),x)
            tvar(iks) = tvar(iks) + dat*dat
            bcov(:,iks) = bcov(:,iks) + x*dat
            forall (i2=1:nequ,i1=1:nequ) acovl(i1,i2) = x(i1)*x(i2) 
            acov(:,:,iks) = acov(:,:,iks) + acovl
          end do

!$OMP END TASK

        end do

!$OMP END SINGLE
!$OMP TASKWAIT  
!$OMP END PARALLEL      

! -------------------------------------------------------------------
! Compute tidal analysis and write results 
! -------------------------------------------------------------------

	if( .not. next_output_d(da_out) ) return

	call get_act_timeline(aline)

        write(6,*) 'Performing tidal analysis at time: ', aline

!       -----------------------------------------------
!       Compute tidal analysis over nkn
!       -----------------------------------------------
        allocate(tideh(ntd,nkn))
        allocate(tideg(ntd,nkn))

!$OMP PARALLEL 
!$OMP SINGLE

        call omp_compute_chunk(nkn,nchunk)
        nthreads = 1
!$      nthreads = omp_get_num_threads()

        do k = 1,nkn,nchunk

!$OMP TASK FIRSTPRIVATE(k) PRIVATE(iks,ikend) &
!$OMP&     SHARED(nequ,ndat,nkn,nchunk,tvar,acov,bcov,tideh,tideg) &
!$OMP&     DEFAULT(NONE)

          call omp_compute_minmax(nchunk,nkn,k,ikend)

          do iks=k,ikend
            call tda_solve(nequ,ndat,tvar(iks),acov(:,:,iks), &
     &          bcov(:,iks),tideh(:,iks),tideg(:,iks))
          end do

!$OMP END TASK

        end do

!$OMP END SINGLE
!$OMP TASKWAIT  
!$OMP END PARALLEL      

!       -----------------------------------------------
!       Writes tidal constrants to the .tide.shy
!       -----------------------------------------------
        id = nint(da_out(4))
        ivar = 900
        do i = 1,ntd
	  if ( const_ar(i)%mask == 1 ) then
            ivar = ivar + 1
            call shy_write_scalar_record2d(id,dtime,ivar,tideh(i,:))
            ivar = ivar + 1
            call shy_write_scalar_record2d(id,dtime,ivar,tideg(i,:))
          end if
        end do

!       -----------------------------------------------
!       Reset variables
!       -----------------------------------------------

        deallocate(tideh)
        deallocate(tideg)
        acov = 0.d0
        tvar = 0.d0
        bcov = 0.d0
	ndat = 0

        end subroutine runtime_tidana

!*********************************************************************

