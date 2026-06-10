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

   use mod_para          ! must define: nrens, nanal, verbose, atime_an, TTAU_0D, TTAU_2D
   use mod_manage_obs    ! must define: o0dlev, o0dtemp, o0dsalt, o2dvel, n_0dlev, ...
   use mod_ens_state     ! must define: Abk, Abk_m
   use levels            ! must define: hlv(:)
   implicit none

   !-----------------------------------------------------------------
   ! Analysis matrices allocated on demand (per assimilation cycle)
   !-----------------------------------------------------------------
   real(dp), allocatable, save :: D(:,:)     ! perturbed measurements
   real(dp), allocatable, save :: D1(:,:)    ! D - HA (perturbed innovations)
   real(dp), allocatable, save :: E(:,:)     ! perturbations applied to obs
   real(dp), allocatable, save :: R(:,:)     ! observation covariance matrix
   real(dp), allocatable, save :: S(:,:)     ! ensemble anomalies in obs space
   real(dp), allocatable, save :: innov(:)   ! innovations using ensemble mean
   real(dp), allocatable, save :: HA(:,:)    ! ensemble model in obs space
   integer, save :: nobs_ok

contains
!===================================================================
! GLOBAL DRIVER — Build all matrices for all observation types
!===================================================================
subroutine make_matrices
   implicit none

   write(*,*) 'Building the assimilation arrays...'

   if (n_0dlev > 0) then
      write(*,*) 'Assimilating: sea level'
      call fill_scalar_0d('0DLEV', n_0dlev, o0dlev)
   end if

   if (n_0dtemp > 0) then
      write(*,*) 'Assimilating: temperature'
      call fill_scalar_0d('0DTEM', n_0dtemp, o0dtemp)
   end if

   if (n_0dsalt > 0) then
      write(*,*) 'Assimilating: salinity'
      call fill_scalar_0d('0DSAL', n_0dsalt, o0dsalt)
   end if

   if (n_2dvel > 0) then
      write(*,*) 'Assimilating: surface currents'
      call fill_scurrents(n_2dvel)
   end if

   if (nobs_ok == 0) error stop 'No valid observation to assimilate'

end subroutine make_matrices


!===================================================================
! 0D SCALARS: Sea level / Temperature / Salinity
!===================================================================
subroutine fill_scalar_0d(olabel, nfile, ostate)
   implicit none

   character(len=5), intent(in)    :: olabel
   integer,         intent(in)     :: nfile
   type(scalar_0d), intent(inout)  :: ostate(nfile)

   integer :: nf, ne
   integer :: iemin, kmin
   integer :: n_obs
   real(dp) :: x, y
   real(dp) :: valo, stdo
   real(dp) :: mvalm
   real(dp) :: valm(nrens)
   real(dp) :: pvec(nrens)
   real(dp) :: inn1
   logical  :: accept_obs

   real(dp), allocatable :: val_o(:)

   integer :: i, count
   integer, allocatable :: valid_idx(:)
   ! Temporary allocatable arrays for shrinking
   real(dp), allocatable :: tmp_v(:), tmp_m(:,:)

!  Fortran 2018 code, use a temporary index for allocation.
   !-------------------------------------------------
   ! PHASE 1: Identify valid observations
   !-------------------------------------------------
   ! Allocate a temporary index array to store nf of valid observations
   allocate(valid_idx(nfile))
   count = 0

   do nf = 1, nfile
      ! Quick check on state
      if (ostate(nf)%stat > 1) cycle

      ! Increment count and store the file index
      count = count + 1
      valid_idx(count) = nf
   end do

   nobs_ok = count

   if (nobs_ok == 0) then
      write(*,*) 'WARNING: no valid scalar observations'
      deallocate(valid_idx)
      return
   end if

   !-------------------------------------------------
   ! PHASE 2: Allocation
   !-------------------------------------------------
   ! Arrays are now sized exactly to nobs_ok
   allocate(val_o(nobs_ok), innov(nobs_ok))
   allocate(R(nobs_ok, nobs_ok), S(nobs_ok, nrens), HA(nobs_ok, nrens))
   allocate(D(nobs_ok, nrens), D1(nobs_ok, nrens), E(nobs_ok, nrens))

   R = 0.0_dp
   n_obs = 0 ! Counter for the actual filled observations (after OBSCHK)

   !-------------------------------------------------
   ! PHASE 3: Main Computation
   !-------------------------------------------------
   ! Now we iterate only over the previously identified candidates
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

      case ('0DSAL')
         do ne = 1, nrens
            valm(ne) = Abk(ne)%s(1, kmin)
         end do
         mvalm = Abk_m%s(1, kmin)
      end select

      ! Get obs values
      valo = ostate(nf)%val
      stdo = ostate(nf)%std
      inn1 = valo - mvalm

      ! Adjust spread/error inflation
      call check_spread(inn1, stdo, valm, mvalm)

      ! Optional: Advanced screening
      if ( OBSCHK ) then
         accept_obs = .true.
         if (nanal > 3) call screen_observation(valo, valm, nrens, stdo, THRSTD, THRABS, accept_obs)
         if (.not. accept_obs) cycle ! If rejected here, this slot in allocated matrices will stay empty or needs packing
      end if

      ! Update global counter for valid output
      n_obs = n_obs + 1

      ! Fill structures using n_obs as the row index
      val_o(n_obs)   = valo
      innov(n_obs)   = inn1
      R(n_obs,n_obs) = stdo**2
      S(n_obs,:)     = valm(:) - mvalm
      HA(n_obs,:)    = valm(:)

      ! Perturbations and Innovations
      call make_0Dpert(olabel, nrens, nanal, ostate(nf)%id, pvec, atime_an, TTAU_0D)

      E(n_obs,:)  = stdo * pvec
      D(n_obs,:)  = valo + E(n_obs,:)
      D1(n_obs,:) = D(n_obs,:) - HA(n_obs,:)

      if (verbose) write(*,*) 'val_m, val_o, std_o, inn: ', mvalm, valo, stdo, inn1
   end do

   ! Cleanup temporary index
   deallocate(valid_idx)

   !-------------------------------------------------
   ! PHASE 4: Shrink matrices to actual accepted obs
   !-------------------------------------------------
   if (n_obs < nobs_ok) then
      if (verbose) write(*,*) 'Shrinking matrices from ', nobs_ok, ' to ', n_obs

      block

         ! Shrink vectors
         allocate(tmp_v(n_obs))

         tmp_v = val_o(1:n_obs);    call move_alloc(tmp_v, val_o)
         allocate(tmp_v(n_obs)) ! Re-allocate for next use
         tmp_v = innov(1:n_obs);    call move_alloc(tmp_v, innov)

         ! Shrink matrices (2D)
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

      ! Update nobs_ok so the rest of the code knows the actual count
      nobs_ok = n_obs
   end if

end subroutine fill_scalar_0d


!===================================================================
! SURFACE CURRENTS (U,V) with OpenMP only in pass 1
!===================================================================
subroutine fill_scurrents(nfile)
   implicit none

   integer, intent(in) :: nfile

   integer :: nf, ix, iy, ne
   integer :: iemin, kmin
   integer :: n_obs
   real(dp) :: x, y
   real(dp) :: valu_o, valv_o, std_o
   real(dp) :: inn1, inn2
   real(dp) :: mvalum, mvalvm
   real(dp) :: h_1st_layer
   real(dp) :: spread_u, spread_v
   real(dp) :: thresh_u, thresh_v
   real(dp) :: mvalu(nrens), mvalv(nrens)
   real(dp) :: pvec1(nrens), pvec2(nrens)
   logical :: accept_u, accept_v

   if (size(hlv) <= 1) then
      error stop 'fill_scurrents: a 3D simulation is required'
   end if

   !---------------------------------------------------------
   ! PASS 1 — count valid u,v observations (OpenMP safe)
   !---------------------------------------------------------
   nobs_ok = 0

!$omp parallel do default(shared) private(nf,ix,iy,x,y,iemin,kmin,ne,          &
!$omp      mvalu,mvalv,mvalum,mvalvm,h_1st_layer,std_o,valu_o,valv_o,          &
!$omp      inn1,inn2,spread_u,spread_v,thresh_u,thresh_v,accept_u,accept_v)    &
!$omp      reduction(+:nobs_ok)
   do nf = 1, nfile

      do iy = 1, o2dvel(nf)%ny
         do ix = 1, o2dvel(nf)%nx

            if (o2dvel(nf)%stat(ix,iy) /= 0) cycle

            x = o2dvel(nf)%x(ix,iy)
            y = o2dvel(nf)%y(ix,iy)
            call find_el_node(x, y, iemin, kmin)

            h_1st_layer = hlv(1) + Abk_m%z(kmin)
            std_o       = o2dvel(nf)%std(ix,iy)

            valu_o = o2dvel(nf)%u(ix,iy) * h_1st_layer
            valv_o = o2dvel(nf)%v(ix,iy) * h_1st_layer

            do ne = 1, nrens
               mvalu(ne) = Abk(ne)%u(1, iemin)
               mvalv(ne) = Abk(ne)%v(1, iemin)
            end do

            mvalum = Abk_m%u(1, iemin)
            mvalvm = Abk_m%v(1, iemin)

            inn1 = valu_o - mvalum
            inn2 = valv_o - mvalvm

            spread_u = sqrt(sum(mvalu**2)/real(nrens,dp) - mvalum*mvalum)
            spread_v = sqrt(sum(mvalv**2)/real(nrens,dp) - mvalvm*mvalvm)

            thresh_u = 3._dp * sqrt(std_o**2 + spread_u**2)
            thresh_v = 3._dp * sqrt(std_o**2 + spread_v**2)

            accept_u = abs(inn1) <= thresh_u
            accept_v = abs(inn2) <= thresh_v

            if (accept_u) nobs_ok = nobs_ok + 1
            if (accept_v) nobs_ok = nobs_ok + 1

         end do
      end do

   end do
!$omp end parallel do


   if (nobs_ok == 0) then
      write(*,*) 'WARNING: no valid current observations'
      return
   end if

   !---------------------------------------------------------
   ! PASS 2 — allocate arrays
   !---------------------------------------------------------
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
   ! PASS 3 — fill arrays
   !---------------------------------------------------------
   do nf = 1, nfile

      call make_0Dpert('2DVEL', nrens, nanal, o2dvel(nf)%id, pvec1, atime_an, TTAU_2D)
      call make_0Dpert('2DVEL', nrens, nanal, o2dvel(nf)%id, pvec2, atime_an, TTAU_2D)

      do iy = 1, o2dvel(nf)%ny
         do ix = 1, o2dvel(nf)%nx

            if (o2dvel(nf)%stat(ix,iy) /= 0) cycle

            x = o2dvel(nf)%x(ix,iy)
            y = o2dvel(nf)%y(ix,iy)
            call find_el_node(x, y, iemin, kmin)

            h_1st_layer = hlv(1) + Abk_m%z(kmin)
            std_o       = o2dvel(nf)%std(ix,iy)

            valu_o = o2dvel(nf)%u(ix,iy) * h_1st_layer
            valv_o = o2dvel(nf)%v(ix,iy) * h_1st_layer

            do ne = 1, nrens
               mvalu(ne) = Abk(ne)%u(1, iemin)
               mvalv(ne) = Abk(ne)%v(1, iemin)
            end do
            mvalum = Abk_m%u(1, iemin)
            mvalvm = Abk_m%v(1, iemin)

            inn1 = valu_o - mvalum
            inn2 = valv_o - mvalvm

            spread_u = sqrt(sum(mvalu**2)/real(nrens,dp) - mvalum*mvalum)
            spread_v = sqrt(sum(mvalv**2)/real(nrens,dp) - mvalvm*mvalvm)

            thresh_u = 3._dp * sqrt(std_o**2 + spread_u**2)
            thresh_v = 3._dp * sqrt(std_o**2 + spread_v**2)

            !-----------------------
            ! u component accepted?
            !-----------------------
            if (abs(inn1) <= thresh_u) then
               n_obs = n_obs + 1

               R(n_obs,n_obs) = std_o**2
               innov(n_obs)   = inn1

               S(n_obs,:)  = mvalu(:) - mvalum
               HA(n_obs,:) = mvalu(:)

               E(n_obs,:)  = std_o * pvec1
               D(n_obs,:)  = valu_o + E(n_obs,:)
               D1(n_obs,:) = D(n_obs,:) - HA(n_obs,:)
            end if

            !-----------------------
            ! v component accepted?
            !-----------------------
            if (abs(inn2) <= thresh_v) then
               n_obs = n_obs + 1

               R(n_obs,n_obs) = std_o**2
               innov(n_obs)   = inn2

               S(n_obs,:)  = mvalv(:) - mvalvm
               HA(n_obs,:) = mvalv(:)

               E(n_obs,:)  = std_o * pvec2
               D(n_obs,:)  = valv_o + E(n_obs,:)
               D1(n_obs,:) = D(n_obs,:) - HA(n_obs,:)
            end if

         end do
      end do
   end do

   if (n_obs /= nobs_ok) error stop 'Mismatch in surface current count'

   write(*,*) 'Number of valid 2D observations: ', nobs_ok

end subroutine fill_scurrents

end module mod_enkf
