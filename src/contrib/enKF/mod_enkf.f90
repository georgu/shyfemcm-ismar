!===============================================================
!  mod_enkf.f90 
!---------------------------------------------------------------
!  Purpose:
!    Build the key matrices/vectors for EnKF analysis using the
!    available observations and ensemble fields:
!      - D   : perturbed observations (nobs_ok x nrens)
!      - D1  : innovation vectors per member = D - HA
!      - E   : obs perturbations (white/red noise) (nobs_ok x nrens)
!      - R   : observation error covariance (nobs_ok x nobs_ok)
!      - S   : model perturbations at obs locations S = HA - mean(HA)
!      - HA  : model (ensemble) mapped to observation space
!      - innov : innovation using ensemble mean = d - H(mean(Abk))
!
!  Notes:
!    * Double precision everywhere (dp = real64).
!    * 0D variables accept observations with stat < 2 (0 or 1).
!    * 2D currents accept only stat == 0; we enforce the same check
!      in counting and filling to keep dimensions consistent.
!
!  Copyright:
!    (C) 2017, Marco Bajo, CNR-ISMAR Venice. All rights reserved.
!    Updated comments and corrections (2026-02-13).
!===============================================================
module mod_enkf
   use iso_fortran_env, only : dp => real64
   use mod_para
   use mod_manage_obs
   use mod_ens_state
   implicit none

   !-----------------------------------------------------------------
   ! Module arrays used during the analysis step (allocated on demand)
   !-----------------------------------------------------------------
   real(dp), save, allocatable :: D(:,:)    ! perturbed measurements
   real(dp), save, allocatable :: D1(:,:)   ! perturbed innovation vectors = D - HA
   real(dp), save, allocatable :: E(:,:)    ! perturbations applied to obs
   real(dp), save, allocatable :: R(:,:)    ! obs error covariance
   real(dp), save, allocatable :: S(:,:)    ! model perturbations at obs (HA - mean(HA))
   real(dp), save, allocatable :: innov(:)  ! innovation using ensemble mean
   real(dp), save, allocatable :: HA(:,:)   ! ensemble model in obs space

contains
!===================================================================
! Build/allocate all matrices/vectors needed for the EnKF analysis.
! This routine:
!   1) counts valid observations (nobs_ok) for the active obs type
!   2) allocates D, D1, E, R, S, HA, innov accordingly
!   3) fills them by calling the proper helpers per obs type
!===================================================================
subroutine make_matrices
   implicit none
   integer :: n, ix, iy

   write(*,*) 'Building the assimilation arrays...'
   nobs_ok = 0

   !--------------------------------------------------------------
   ! Count valid observations for the ACTIVE observation type only
   ! (exactly one of the following flags must be non-zero)
   !--------------------------------------------------------------
   if (islev /= 0) then
      do n = 1, n_0dlev
         if (o0dlev(n)%stat < 2) nobs_ok = nobs_ok + 1
      end do

   else if (issalt /= 0) then
      do n = 1, n_0dsalt
         if (o0dsalt(n)%stat < 2) nobs_ok = nobs_ok + 1
      end do

   else if (istemp /= 0) then
      do n = 1, n_0dtemp
         if (o0dtemp(n)%stat < 2) nobs_ok = nobs_ok + 1
      end do

   else if (isvel /= 0) then
      ! Each valid grid-point contributes TWO observations (u and v).
      do n = 1, n_2dvel
         do iy = 1, o2dvel(n)%ny
            do ix = 1, o2dvel(n)%nx
               if (o2dvel(n)%stat(ix,iy) == 0) nobs_ok = nobs_ok + 2
            end do
         end do
      end do

   else
      error stop 'Observation type not valid.'
   end if

   write(*,*) 'Valid observations: ', nobs_ok

   !-------------------------------------
   ! Allocate analysis arrays/matrices
   !-------------------------------------
   allocate(D(nobs_ok,  nrens),  &
            E(nobs_ok,  nrens),  &
            R(nobs_ok,  nobs_ok))
   allocate(S(nobs_ok,  nrens),  innov(nobs_ok))
   allocate(HA(nobs_ok, nrens),  D1(nobs_ok, nrens))

   R = 0.0_dp

   !-------------------------------------
   ! Fill by observation type
   !-------------------------------------
   ! 0D sea level
   if (n_0dlev  > 0) then
      write(*,*) 'Assimilation of sea level'
      call fill_scalar_0d('0DLEV', n_0dlev,  o0dlev)
   end if

   ! 0D temperature
   if (n_0dtemp > 0) then
      write(*,*) 'Assimilation of temperature'
      call fill_scalar_0d('0DTEM', n_0dtemp, o0dtemp)
   end if

   ! 0D salinity
   if (n_0dsalt > 0) then
      write(*,*) 'Assimilation of salinity'
      call fill_scalar_0d('0DSAL', n_0dsalt, o0dsalt)
   end if

   ! 2D surface currents (requires 3D sim to get first-layer thickness)
   if (n_2dvel > 0) then
      write(*,*) 'Assimilation of velocities'
      call fill_scurrents(n_2dvel)
   end if
end subroutine make_matrices

!===================================================================
! Fill matrices for 0D scalar observations (sea level, temperature,
! salinity). For each file entry, if stat < 2 the observation is used.
! olabel ∈ {'0DLEV','0DTEM','0DSAL'}
!===================================================================
subroutine fill_scalar_0d(olabel, nfile, ostate)
   implicit none
   character(len=5), intent(in)    :: olabel
   integer,         intent(in)     :: nfile
   type(scalar_0d), intent(inout)  :: ostate(nfile)

   integer :: nf, ne, nook
   integer :: iemin, kmin
   real(dp) :: x, y
   real(dp) :: oval, stdv, stdm
   real(dp) :: inn1
   real(dp) :: mvalm
   real(dp) :: mval(nrens)     ! model values at obs location for each ensemble member
   real(dp) :: pvec(nrens)     ! N(0,1) (or red noise) perturbation vector

   nook = 0

   do nf = 1, nfile

      ! Create Gaussian/red-noise perturbations with mean 0 and std 1
      call make_0Dpert(olabel, nrens, nanal, ostate(nf)%id, pvec, atime_an, TTAU_0D)

      ! Skip if observation is bad (we accept stat = 0 or 1)
      if (ostate(nf)%stat > 1) cycle

      nook = nook + 1

      ! Nearest model node/element to this observation
      x = ostate(nf)%x
      y = ostate(nf)%y
      call find_el_node(x, y, iemin, kmin)
      ! if (verbose) write(*,*) 'Internal node nearest to obs: ', kmin

      ! Observation error variance on the diagonal of R
      R(nook, nook) = (ostate(nf)%std)**2

      !-------------------------------------------------------------
      ! Map the ensemble state to observation space (HA)
      ! and compute anomalies S = HA - mean(HA)
      !-------------------------------------------------------------
      select case (olabel)
      case ('0DLEV')
         do ne = 1, nrens
            mval(ne) = Abk(ne)%z(kmin)
         end do
         mvalm = Abk_m%z(kmin)

      case ('0DTEM')
         do ne = 1, nrens
            mval(ne) = Abk(ne)%t(1, kmin)
         end do
         mvalm = Abk_m%t(1, kmin)

      case ('0DSAL')
         do ne = 1, nrens
            mval(ne) = Abk(ne)%s(1, kmin)
         end do
         mvalm = Abk_m%s(1, kmin)

      end select

      ! Ensemble spread at the obs location (use dp-safe arithmetic)
      stdm = sqrt( sum(mval**2) / real(nrens, dp) - (sum(mval) / real(nrens, dp))**2 )

      S(nook, :)  = mval(:) - mvalm
      HA(nook, :) = mval(:)

      ! Innovation using ensemble mean
      oval = ostate(nf)%val
      stdv = ostate(nf)%std
      inn1 = oval - mvalm
      innov(nook) = inn1

      ! Optional consistency checks on spread/innovation
      call check_spread(inn1, stdv, mval, mvalm)

      if (verbose) then
         write(*,'(a26,1x,i4,2f8.4,1x,3f8.3)') 'ID obs,x,y,vobs,vmod,inn: ', &
              nf, x, y, oval, mvalm, inn1
      end if

      !-------------------------------------------------------------
      ! Build perturbations and innovation per member
      !   E  = stdv * pvec
      !   D  = oval + E
      !   D1 = D - HA
      !-------------------------------------------------------------
      E(nook, :)  = stdv * pvec
      D(nook, :)  = oval + E(nook, :)
      D1(nook, :) = D(nook, :) - HA(nook, :)

   end do

   ! For a single active obs type, the total good observations must match
   if (nook /= nobs_ok) error stop 'The number of good observations is wrong.'
end subroutine fill_scalar_0d

!===================================================================
! Fill matrices for 2D surface current observations.
! Each valid grid cell contributes two observations: u and v.
! We require a 3D simulation to compute first-layer thickness (hlv(1))
! to convert observed surface velocities into transports consistent
! with model layer-integrated u/v (if that is your convention).
!===================================================================
subroutine fill_scurrents(nfile)
   use levels                ! provides hlv (layer thicknesses), among others
   implicit none

   integer, intent(in) :: nfile

   integer :: nf, ix, iy, ne, nook
   integer :: iemin, kmin
   real(dp) :: x, y
   real(dp) :: uu, vv, stdv
   real(dp) :: inn1, inn2
   real(dp) :: mvalum, mvalvm
   real(dp) :: h_1st_layer
   real(dp) :: mvalu(nrens), mvalv(nrens)
   real(dp) :: pvec1(nrens), pvec2(nrens)

   ! Require 3D (at least two layers incl. free-surface contribution)
   if (size(hlv) <= 1) then
      error stop 'fill_scurrents: a 3D simulation is necessary to assimilate surface currents'
   end if

   nook = 0

   do nf = 1, nfile

      ! Independent perturbations for u and v components
      call make_0Dpert('2DVEL', nrens, nanal, o2dvel(nf)%id, pvec1, atime_an, TTAU_2D)
      call make_0Dpert('2DVEL', nrens, nanal, o2dvel(nf)%id, pvec2, atime_an, TTAU_2D)

      do iy = 1, o2dvel(nf)%ny
         do ix = 1, o2dvel(nf)%nx

            ! Accept only strictly good points: stat == 0.
            ! This matches the counting in make_matrices and prevents mismatches.
            if (o2dvel(nf)%stat(ix,iy) /= 0) cycle

            nook = nook + 2   ! two observations: u and v

            x = o2dvel(nf)%x(ix,iy)
            y = o2dvel(nf)%y(ix,iy)

            call find_el_node(x, y, iemin, kmin)
            if (verbose) write(*,*) 'Internal element nearest to obs: ', iemin

            ! Effective first-layer thickness (static + free surface)
            h_1st_layer = hlv(1) + Abk_m%z(kmin)

            ! Observation error variances (diagonal R)
            stdv = o2dvel(nf)%std(ix,iy)
            R(nook-1, nook-1) = stdv**2
            R(nook,   nook   ) = stdv**2

            ! Model values (ensemble and mean) at obs location
            do ne = 1, nrens
               mvalu(ne) = Abk(ne)%u(1, iemin)
               mvalv(ne) = Abk(ne)%v(1, iemin)
            end do
            mvalum = Abk_m%u(1, iemin)
            mvalvm = Abk_m%v(1, iemin)

            ! Ensemble anomalies and HA
            S(nook-1, :)  = mvalu(:) - mvalum
            S(nook,   :)  = mvalv(:) - mvalvm
            HA(nook-1, :) = mvalu(:)
            HA(nook,   :) = mvalv(:)

            ! Innovation using ensemble mean.
            ! If the model stores layer-transports, convert obs to transport
            ! using the first-layer thickness. Adjust this section if your
            ! model/obs are in strictly velocity units.
            uu = o2dvel(nf)%u(ix,iy) * h_1st_layer
            vv = o2dvel(nf)%v(ix,iy) * h_1st_layer
            inn1 = uu - mvalum
            inn2 = vv - mvalvm
            innov(nook-1) = inn1
            innov(nook)   = inn2

            call check_spread_speed(inn1, inn2, stdv, mvalu, mvalv, mvalum, mvalvm)

            ! Build perturbations and member-wise innovations
            ! -- u component
            E(nook-1, :)  = stdv * pvec1
            D(nook-1, :)  = uu   + E(nook-1, :)
            D1(nook-1, :) = D(nook-1, :) - HA(nook-1, :)

            ! -- v component
            E(nook, :)    = stdv * pvec2
            D(nook, :)    = vv   + E(nook, :)
            D1(nook, :)   = D(nook, :) - HA(nook, :)

         end do
      end do
   end do

   if (nook /= nobs_ok) error stop 'The number of good observations is wrong.'
end subroutine fill_scurrents

end module mod_enkf
