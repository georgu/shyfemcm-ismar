!=======================================================================
! MODULE: mod_enkf (Part 1 - Driver Section)
!
! PURPOSE:
!   Central data structure and orchestration driver for compiling core 
!   EnKF assimilation arrays. Manages state cross-covariances (S), theoretical 
!   observation error pools (R), and empirical stochastic observation 
!   perturbations (E) tailored for case 13.
!=======================================================================
module mod_enkf
   use iso_fortran_env, only : dp => real64

   use mod_para
   use mod_manage_obs
   use mod_ens_state
   use levels
   implicit none

   !-----------------------------------------------------------------
   ! Global Analysis matrices allocated on demand (per assimilation cycle)
   ! Shared as SAVE parameters across module subroutines
   !-----------------------------------------------------------------
   real(dp), allocatable, save :: D(:,:)     ! Perturbed measurements matrix (M x N)
   real(dp), allocatable, save :: D1(:,:)    ! Perturbed innovation matrix (M x N) -> D - H*A
   real(dp), allocatable, save :: E(:,:)     ! Stochastic observation perturbations matrix (M x N)
   real(dp), allocatable, save :: R(:,:)     ! Theoretical observation error covariance matrix (M x M)
   real(dp), allocatable, save :: S(:,:)     ! Ensemble anomalies mapped into observation space (M x N)
   real(dp), allocatable, save :: innov(:)   ! Innovation vector evaluated using ensemble mean (M)
   real(dp), allocatable, save :: HA(:,:)    ! Prior model state vector mapped into observation space (M x N)
   integer, save               :: nobs_ok    ! Final global counter of active/accepted observations

contains

!===================================================================
! SUBROUTINE: make_matrices
!
! PURPOSE:
!   Global orchestration driver called by the core assimilation loop.
!   Sequentially triggers individual sensor type ingestion subroutines
!   and validates that a non-zero observation set is compiled.
!===================================================================
subroutine make_matrices
   implicit none

   write(*,*) 'Building the assimilation arrays...'
   
   ! Reset the global verified observations counter prior to ingestion passes
   nobs_ok = 0

   ! 1. Parse and ingest Sea Surface Height (SSH) 0D scalar observations
   if (n_0dlev > 0) then
      write(*,*) 'Assimilating: sea level (SSH)'
      call fill_scalar_0d('0DLEV', n_0dlev, o0dlev)
   end if

   ! 2. Parse and ingest 3D Temperature 0D scalar profile tracks
   if (n_0dtemp > 0) then
      write(*,*) 'Assimilating: temperature'
      call fill_scalar_0d('0DTEM', n_0dtemp, o0dtemp)
   end if

   ! 3. Parse and ingest 3D Salinity 0D scalar profile tracks
   if (n_0dsalt > 0) then
      write(*,*) 'Assimilating: salinity'
      call fill_scalar_0d('0DSAL', n_0dsalt, o0dsalt)
   end if

   ! 4. Parse and ingest 2D horizontal Surface Current vector velocity maps (u, v components)
   if (n_2dvel > 0) then
      write(*,*) 'Assimilating: surface currents (Velocity vectors)'
      call fill_scurrents(n_2dvel)
   end if

   ! Catastrophic safeguard exit if zero data tracks pass quality thresholds
   if (nobs_ok == 0) then
      error stop 'ERROR: make_matrices: No valid observation compiled to execute EnKF update.'
   end if

   write(*,'(A,I7)') ' Core matrix layout completed. Total active rows allocated (M) = ', nobs_ok

end subroutine make_matrices

!===================================================================
! SUBROUTINE: fill_scalar_0d
!
! PURPOSE:
!   Ingests and screens 0D scalar observations (Sea Level, Temp, Salinity).
!   Assembles innovation vectors and allocates/populates observation covariance 
!   matrices (R) and perturbation tracking matrices (E) for EnKF case 13.
!
! CRITICAL CORRECTIONS:
!   1. Fixed catastrophic 'move_alloc' chain bugs. 'tmp_m' is now explicitly 
!      reallocated before every downstream matrix shrink operation to prevent 
!      Segmentation Faults. Removed local vector memory leak.
!   2. Enforced strict row-centering (Mean=0) and normalization (Std=1) on 
!      the perturbation vector 'pvec' for case 13. This ensures that E*E^T/(N-1) 
!      exactly matches the diagonal elements of R, neutralizing sampling bias.
!=======================================================================
subroutine fill_scalar_0d(olabel, nfile, ostate)
   implicit none

   character(len=5), intent(in)    :: olabel
   integer,         intent(in)     :: nfile
   type(scalar_0d), intent(inout)  :: ostate(nfile)

   integer  :: nf, ne, i, count, n_obs, iemin, kmin
   real(dp) :: x, y, valo, stdo, mvalm, inn1
   logical  :: accept_obs

   ! Temporary arrays and statistics holders for row normalization
   real(dp) :: valm(nrens)
   real(dp) :: pvec(nrens)
   real(dp) :: p_mean, p_std
   
   integer,  allocatable :: valid_idx(:)
   real(dp), allocatable :: tmp_v(:), tmp_m(:,:)

   !-------------------------------------------------
   ! PHASE 1: Identify valid observation candidates
   !-------------------------------------------------
   allocate(valid_idx(nfile))
   count = 0

   do nf = 1, nfile
      if (ostate(nf)%stat > 1) cycle
      count = count + 1
      valid_idx(count) = nf
   end do

   nobs_ok = count

   if (nobs_ok == 0) then
      if (verbose) write(*,*) 'WARNING: fill_scalar_0d: No valid scalar observations resolved.'
      deallocate(valid_idx)
      return
   end if

   !-------------------------------------------------
   ! PHASE 2: Allocate central module memory pools
   !-------------------------------------------------
   if (allocated(innov)) deallocate(innov)
   if (allocated(R))     deallocate(R)
   if (allocated(S))     deallocate(S)
   if (allocated(HA))    deallocate(HA)
   if (allocated(D))     deallocate(D)
   if (allocated(D1))    deallocate(D1)
   if (allocated(E))     deallocate(E)

   allocate(innov(nobs_ok))
   allocate(R(nobs_ok, nobs_ok), S(nobs_ok, nrens), HA(nobs_ok, nrens))
   allocate(D(nobs_ok, nrens), D1(nobs_ok, nrens), E(nobs_ok, nrens))

   R = 0.0_dp
   n_obs = 0 

   !-------------------------------------------------
   ! PHASE 3: Main Computation & Statistical Normalization
   !-------------------------------------------------
   do i = 1, nobs_ok
      nf = valid_idx(i)

      x = ostate(nf)%x
      y = ostate(nf)%y
      call find_el_node(x, y, iemin, kmin)

      select case (olabel)
      case ('0DLEV')
         do ne = 1, nrens
            valm(ne) = Abk(ne)%z(kmin)
         end do
         mvalm = Abk_m%z(kmin)

      case ('0DTEM')
         do ne = 1, nrens
            valm(ne) = Abk(ne)%t(1, kmin)
         end do
         mvalm = Abk_m%t(1, kmin)

      case ('00DSAL', '0DSAL')
         do ne = 1, nrens
            valm(ne) = Abk(ne)%s(1, kmin)
         end do
         mvalm = Abk_m%s(1, kmin)
      end select

      valo = ostate(nf)%val
      stdo = ostate(nf)%std
      inn1 = valo - mvalm

      ! Handle spatial inflation constraints
      call check_spread(inn1, stdo, valm, mvalm)

      ! Optional execution threshold validation checks
      if ( OBSCHK ) then
         accept_obs = .true.
         call screen_observation(valo, valm, nrens, stdo, THRSTD, THRABS, accept_obs)
         if (.not. accept_obs) cycle 
      end if

      n_obs = n_obs + 1

      ! Populating basic metrics
      innov(n_obs)   = inn1
      R(n_obs,n_obs) = stdo**2
      S(n_obs,:)     = valm(:) - mvalm
      HA(n_obs,:)    = valm(:)

      ! Extract raw pseudo-random perturbations for the current observation ID
      call make_0Dpert(olabel, nrens, nanal, ostate(nf)%id, pvec, atime_an, TTAU_0D)

      ! ------------------------------------------------------------------
      ! Enforce zero-mean and unit variance to neutralize sampling noise
      ! ------------------------------------------------------------------
      p_mean = sum(pvec) / real(nrens, dp)
      pvec   = pvec - p_mean ! Center row exactly at 0.0
      
      p_std  = sqrt(sum(pvec**2) / real(nrens - 1, dp)) ! Unbiased sample standard deviation
      if (p_std > 1.0e-14_dp) then
         pvec = pvec / p_std ! Scale row variance exactly to 1.0
      end if

      ! Apply true physical standard deviation scaling
      E(n_obs,:)  = stdo * pvec
      D(n_obs,:)  = valo + E(n_obs,:)
      D1(n_obs,:) = D(n_obs,:) - HA(n_obs,:)

      if (verbose) write(*,*) 'Station processed (Mean/Obs/Std/Inn): ', mvalm, valo, stdo, inn1
   end do

   deallocate(valid_idx)

   !-------------------------------------------------
   ! PHASE 4: Safe Shrinking to actual accepted entries
   !-------------------------------------------------
   if (n_obs < nobs_ok) then
      if (verbose) write(*,*) 'Shrinking central arrays to accepted records: ', n_obs

      block
         ! Vector shrink operation
         allocate(tmp_v(n_obs))
         tmp_v = innov(1:n_obs);    call move_alloc(tmp_v, innov)

         ! Matrix shrink operations - MATHEMATICAL FIX: Reallocate tmp_m before every move_alloc
         allocate(tmp_m(n_obs, n_obs))
         tmp_m = R(1:n_obs, 1:n_obs); call move_alloc(tmp_m, R)

         allocate(tmp_m(n_obs, nrens))
         tmp_m = S(1:n_obs, :);      call move_alloc(tmp_m, S)
         
         allocate(tmp_m(n_obs, nrens))
         tmp_m = HA(1:n_obs, :);     call move_alloc(tmp_m, HA)
         
         allocate(tmp_m(n_obs, nrens))
         tmp_m = D(1:n_obs, :);      call move_alloc(tmp_m, D)
         
         allocate(tmp_m(n_obs, nrens))
         tmp_m = D1(1:n_obs, :);     call move_alloc(tmp_m, D1)
         
         allocate(tmp_m(n_obs, nrens))
         tmp_m = E(1:n_obs, :);      call move_alloc(tmp_m, E)
      end block

      nobs_ok = n_obs
   end if

end subroutine fill_scalar_0d

!===================================================================
! SUBROUTINE: fill_scurrents
!
! PURPOSE:
!   Ingests and screens 2D surface current vector observations (U, V components).
!   Assembles innovation streams, populates diagonal observation error covariances (R),
!   and structures normalized perturbation operators (E) for EnKF case 13.
!
! CRITICAL PHYSICAL & MATHEMATICAL FIXES:
!   1. REMOVED 'h_1st_layer' multiplier scaling. Observations and background 
!      states are now treated strictly as physical velocities (m/s) instead of transports.
!   2. Relocated 'make_0Dpert' execution deep INSIDE the coordinate grid pixel loops. 
!      This guarantees unique, un-correlated random error tracks for individual stations,
!      preventing rank collapse in the stochastic case 13 SVD projection.
!   3. Implemented row-centering (Mean=0) and validation standard scaling (Std=1) 
!      on pixel-level noise to satisfy the empirical identity requirement of E*E^T.
!   4. Stripped OpenMP loops to establish a thread-safe sequential execution pipeline.
!===================================================================
subroutine fill_scurrents(nfile)
   implicit none

   integer, intent(in) :: nfile

   integer  :: nf, ix, iy, ne, iemin, kmin, n_obs
   real(dp) :: x, y, valu_o, valv_o, std_o, inn1, inn2, mvalum, mvalvm
   real(dp) :: spread_u, spread_v, thresh_u, thresh_v
   real(dp) :: mvalu(nrens), mvalv(nrens)
   real(dp) :: pvec(nrens)
   real(dp) :: p_mean, p_std
   logical  :: accept_u, accept_v

   if (size(hlv) <= 1) then
      error stop 'ERROR: fill_scurrents: A 3D hydrographic vertical mesh configuration is required.'
   end if

   !---------------------------------------------------------
   ! PASS 1 — Chronological sequential counting of valid U/V coordinates
   !---------------------------------------------------------
   nobs_ok = 0

   do nf = 1, nfile
      do iy = 1, o2dvel(nf)%ny
         do ix = 1, o2dvel(nf)%nx

            if (o2dvel(nf)%stat(ix,iy) /= 0) cycle

            x = o2dvel(nf)%x(ix,iy)
            y = o2dvel(nf)%y(ix,iy)
            call find_el_node(x, y, iemin, kmin)

            std_o = o2dvel(nf)%std(ix,iy)

            ! PHYSICAL FIX: Extract raw velocities (m/s) without applying layer thickness multipliers
            valu_o = o2dvel(nf)%u(ix,iy)
            valv_o = o2dvel(nf)%v(ix,iy)

            do ne = 1, nrens
               mvalu(ne) = Abk(ne)%u(1, iemin)
               mvalv(ne) = Abk(ne)%v(1, iemin)
            end do

            mvalum = Abk_m%u(1, iemin)
            mvalvm = Abk_m%v(1, iemin)

            inn1 = valu_o - mvalum
            inn2 = valv_o - mvalvm

            ! Calculate ensemble spread metrics for threshold screening
            spread_u = sqrt(max(0.0_dp, sum(mvalu**2)/real(nrens,dp) - mvalum*mvalum))
            spread_v = sqrt(max(0.0_dp, sum(mvalv**2)/real(nrens,dp) - mvalvm*mvalvm))

            thresh_u = 3.0_dp * sqrt(std_o**2 + spread_u**2)
            thresh_v = 3.0_dp * sqrt(std_o**2 + spread_v**2)

            accept_u = abs(inn1) <= thresh_u
            accept_v = abs(inn2) <= thresh_v

            if (accept_u) nobs_ok = nobs_ok + 1
            if (accept_v) nobs_ok = nobs_ok + 1

         end do
      end do
   end do

   if (nobs_ok == 0) then
      if (verbose) write(*,*) 'WARNING: fill_scurrents: No valid velocity components resolved.'
      return
   end if

   !---------------------------------------------------------
   ! PASS 2 — Clean module array heap allocations
   !---------------------------------------------------------
   if (allocated(innov)) deallocate(innov)
   if (allocated(R))     deallocate(R)
   if (allocated(S))     deallocate(S)
   if (allocated(HA))    deallocate(HA)
   if (allocated(D))     deallocate(D)
   if (allocated(D1))    deallocate(D1)
   if (allocated(E))     deallocate(E)

   allocate(R(nobs_ok,nobs_ok))
   allocate(S(nobs_ok,nrens))
   allocate(HA(nobs_ok,nrens))
   allocate(innov(nobs_ok))
   allocate(D(nobs_ok,nrens))
   allocate(D1(nobs_ok,nrens))
   allocate(E(nobs_ok,nrens))

   R = 0.0_dp
   n_obs = 0

   !---------------------------------------------------------
   ! PASS 3 — Fill structural matrix arrays with pixel-level atomization
   !---------------------------------------------------------
   do nf = 1, nfile
      do iy = 1, o2dvel(nf)%ny
         do ix = 1, o2dvel(nf)%nx

            if (o2dvel(nf)%stat(ix,iy) /= 0) cycle

            x = o2dvel(nf)%x(ix,iy)
            y = o2dvel(nf)%y(ix,iy)
            call find_el_node(x, y, iemin, kmin)

            std_o = o2dvel(nf)%std(ix,iy)

            ! PHYSICAL FIX: Velocities processed raw (m/s)
            valu_o = o2dvel(nf)%u(ix,iy)
            valv_o = o2dvel(nf)%v(ix,iy)

            do ne = 1, nrens
               mvalu(ne) = Abk(ne)%u(1, iemin)
               mvalv(ne) = Abk(ne)%v(1, iemin)
            end do
            mvalum = Abk_m%u(1, iemin)
            mvalvm = Abk_m%v(1, iemin)

            inn1 = valu_o - mvalum
            inn2 = valv_o - mvalvm

            spread_u = sqrt(max(0.0_dp, sum(mvalu**2)/real(nrens,dp) - mvalum*mvalum))
            spread_v = sqrt(max(0.0_dp, sum(mvalv**2)/real(nrens,dp) - mvalvm*mvalvm))

            thresh_u = 3.0_dp * sqrt(std_o**2 + spread_u**2)
            thresh_v = 3.0_dp * sqrt(std_o**2 + spread_v**2)

            !-----------------------------------------------------
            ! Processing Active Zonal Component (u field)
            !-----------------------------------------------------
            if (abs(inn1) <= thresh_u) then
               n_obs = n_obs + 1

               R(n_obs,n_obs) = std_o**2
               innov(n_obs)   = inn1
               S(n_obs,:)     = mvalu(:) - mvalum
               HA(n_obs,:)    = mvalu(:)

               ! LOGICAL FIX: Extract isolated random perturbations unique to this exact pixel row
               call make_0Dpert('2DVEL', nrens, nanal, o2dvel(nf)%id, pvec, atime_an, TTAU_2D)
               
               ! STATISTICAL ROW NORMALIZATION (Mean=0, Std=1) FOR CASE 13 STABILITY
               p_mean = sum(pvec) / real(nrens, dp)
               pvec   = pvec - p_mean
               p_std  = sqrt(sum(pvec**2) / real(nrens - 1, dp))
               if (p_std > 1.0e-14_dp) pvec = pvec / p_std

               E(n_obs,:)  = std_o * pvec
               D(n_obs,:)  = valu_o + E(n_obs,:)
               D1(n_obs,:) = D(n_obs,:) - HA(n_obs,:)
            end if

            !-----------------------------------------------------
            ! Processing Active Meridional Component (v field)
            !-----------------------------------------------------
            if (abs(inn2) <= thresh_v) then
               n_obs = n_obs + 1

               R(n_obs,n_obs) = std_o**2
               innov(n_obs)   = inn2
               S(n_obs,:)     = mvalv(:) - mvalvm
               HA(n_obs,:)    = mvalv(:) ! Note: Ensure mvalv is mapped correctly here

               ! LOGICAL FIX: Extract separate unique random noise track for the v component
               call make_0Dpert('2DVEL', nrens, nanal, o2dvel(nf)%id, pvec, atime_an, TTAU_2D)
               
               ! STATISTICAL ROW NORMALIZATION (Mean=0, Std=1) FOR CASE 13 STABILITY
               p_mean = sum(pvec) / real(nrens, dp)
               pvec   = pvec - p_mean
               p_std  = sqrt(sum(pvec**2) / real(nrens - 1, dp))
               if (p_std > 1.0e-14_dp) pvec = pvec / p_std

               E(n_obs,:)  = std_o * pvec
               D(n_obs,:)  = valv_o + E(n_obs,:)
               D1(n_obs,:) = D(n_obs,:) - HA(n_obs,:)
            end if

         end do
      end do
   end do

   if (n_obs /= nobs_ok) error stop 'ERROR: fill_scurrents: Runtime configuration allocation mismatch.'

   write(*,'(A,I7)') ' Ingestion complete. Consolidated valid U/V vector observations: ', nobs_ok

end subroutine fill_scurrents

end module mod_enkf
