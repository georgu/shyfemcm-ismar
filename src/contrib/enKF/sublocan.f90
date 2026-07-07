!=======================================================================
! SUBROUTINE: local_analysis
!
! PURPOSE:
!   Orchestrates spatial localization grid-point loops. Gathers multivariate 
!   sensor positions and routes states safely to per-node (locan_k) and 
!   per-element (locan_e) localized filters under OpenMP thread isolation.
!=======================================================================
subroutine local_analysis()
  use iso_fortran_env,     only : dp => real64
  use mod_ens_state,       only : nnkn, nnel, nnlv, Abk, Aan, nrens
  use mod_manage_obs,      only : o0dlev, o0dtemp, o0dsalt, o2dvel, &
                                  n_0dlev, n_0dtemp, n_0dsalt, n_2dvel
  use mod_restart,         only : ibarcl_rst
  use basin,               only : xgv, ygv, nen3v
  use mod_enkf,            only : nobs_ok, S, R, E, D, D1 
  use mod_para,            only : verbose
  implicit none

  ! ----------------------------- DECLARATIONS -----------------------------
  integer :: no, nf, ix, iy, i_obs
  real(dp), allocatable :: xobs(:), yobs(:), rho_loc(:)

  integer :: lkdim, lnedim
  real(dp), allocatable :: ex(:), ey(:)
  integer :: ne, k

  integer :: nk_l_local, ne_l_local
  logical :: local_upd_rand

  ! Per-thread master workspace pointers (Thread-safe private handles)
  real(dp), allocatable :: Ak_bk(:,:), Ak_loc(:,:)
  logical :: did_update
  real(dp), allocatable :: Ae_bk(:,:), Ae_loc(:,:)
  logical :: did_update_e
  
  integer :: k_start
  ! ------------------------------------------------------------------------

  ! 1. Establish state vector block limits depending on baroclinic physics setup
  if (ibarcl_rst == 0) then
     lkdim  = 1                 ! Barotropic: SSH (Z) only at grid nodes
     lnedim = 2*nnlv            ! 3D velocities (u,v) at grid elements
  else
     lkdim  = 1 + 2*nnlv        ! Baroclinic: Z + Temp profile + Salt profile at nodes
     lnedim = 2*nnlv            ! 3D velocity profiles (u,v) at elements
  end if

  ! ------------------------------------------------------------------------
  ! Step 2: Precompute finite element centroids (u,v element anchors)
  ! ------------------------------------------------------------------------
  allocate(ex(nnel), ey(nnel))
  do ne = 1, nnel
     ex(ne) = ( xgv(nen3v(1,ne)) + xgv(nen3v(2,ne)) + xgv(nen3v(3,ne)) ) / 3.0_dp
     ey(ne) = ( ygv(nen3v(1,ne)) + ygv(nen3v(2,ne)) + ygv(nen3v(3,ne)) ) / 3.0_dp
  end do

  ! ------------------------------------------------------------------------
  ! Step 3: Secure coordinate mapping for active multivariate observations
  ! ------------------------------------------------------------------------
  ! MATHEMATICAL FIX: Bypassed the broken PASS 1 manual count. 'nobs_ok' 
  ! is already verified and consolidated by the upstream module compilers.
  if (nobs_ok <= 0) then
     if (verbose) write(*,*) 'local_analysis: Zero active rows detected. Localization bypassed.'
     deallocate(ex, ey)
     return
  end if

  ! Safe allocation matching the active rows exactly
  allocate(xobs(nobs_ok), yobs(nobs_ok), rho_loc(nobs_ok))
  
  xobs    = 0.0_dp
  yobs    = 0.0_dp
  rho_loc = 0.0_dp
  i_obs   = 0

  ! --- Sequential Extraction: Safely mapping active dimensions ---
  
  ! A. Ingest verified Sea Level positions (0DLEV)
  if (allocated(o0dlev) .and. n_0dlev > 0) then
     do no = 1, n_0dlev
        if (o0dlev(no)%stat < 2 .and. i_obs < nobs_ok) then
           i_obs = i_obs + 1
           xobs(i_obs)    = o0dlev(no)%x
           yobs(i_obs)    = o0dlev(no)%y
           rho_loc(i_obs) = o0dlev(no)%rhol
        end if
     end do
  end if

  ! B. Ingest verified Temperature tracks (0DTEM)
  if (allocated(o0dtemp) .and. n_0dtemp > 0) then
     do no = 1, n_0dtemp
        if (o0dtemp(no)%stat < 2 .and. i_obs < nobs_ok) then
           i_obs = i_obs + 1
           xobs(i_obs)    = o0dtemp(no)%x
           yobs(i_obs)    = o0dtemp(no)%y
           rho_loc(i_obs) = o0dtemp(no)%rhol
        end if
     end do
  end if

  ! C. Ingest verified Salinity tracks (0DSAL)
  if (allocated(o0dsalt) .and. n_0dsalt > 0) then
     do no = 1, n_0dsalt
        if (o0dsalt(no)%stat < 2 .and. i_obs < nobs_ok) then
           i_obs = i_obs + 1
           xobs(i_obs)    = o0dsalt(no)%x
           yobs(i_obs)    = o0dsalt(no)%y
           rho_loc(i_obs) = o0dsalt(no)%rhol
        end if
     end do
  end if

  ! D. Ingest verified Surface Velocity pixels (2DVEL maps for U and V components)
  if (allocated(o2dvel) .and. n_2dvel > 0) then
     do nf = 1, n_2dvel
        do iy = 1, o2dvel(nf)%ny
           do ix = 1, o2dvel(nf)%nx
              if (o2dvel(nf)%stat(ix,iy) == 0) then
                 ! Component U track row
                 i_obs = i_obs + 1
                 if (i_obs <= nobs_ok) then
                    xobs(i_obs)    = o2dvel(nf)%x(ix,iy)
                    yobs(i_obs)    = o2dvel(nf)%y(ix,iy)
                    rho_loc(i_obs) = o2dvel(nf)%std(ix,iy) ! Using std as localization proxy if rhol missing
                 end if
                 ! Component V track row
                 i_obs = i_obs + 1
                 if (i_obs <= nobs_ok) then
                    xobs(i_obs)    = o2dvel(nf)%x(ix,iy)
                    yobs(i_obs)    = o2dvel(nf)%y(ix,iy)
                    rho_loc(i_obs) = o2dvel(nf)%std(ix,iy)
                 end if
              end if
           end do
        end do
     end do
  end if

  ! Final cross-check mapping validation guard
  if (i_obs /= nobs_ok) then
     write(*,'(A,I8,A,I8)') 'ERROR: local_analysis: Spatial index mismatch. Compiled: ', i_obs, ' | Expected: ', nobs_ok
     error stop
  end if

  nk_l_local = 0
  ne_l_local = 0
  local_upd_rand = .true.

  ! ========================= SEQUENTIAL FIRST STEP =========================
  ! Process node k=1 sequentially to handle random rotation generation safely.
  k_start = 1
  if (nnkn >= 1) then
     allocate(Ak_bk(lkdim,nrens), Ak_loc(lkdim,nrens))
     call type_to_kmat(ibarcl_rst, Ak_bk, 1, lkdim, nrens)
     
     call locan_k(1, lkdim, nrens, nobs_ok, xobs, yobs, rho_loc, &
                  Ak_bk, local_upd_rand, Ak_loc, did_update)

     if (did_update) nk_l_local = nk_l_local + 1
     call kmat_to_type(ibarcl_rst, Ak_loc, 1, lkdim, nrens)
     deallocate(Ak_bk, Ak_loc)
     k_start = 2 
  end if

  ! =============================== NODE PHASE =============================
  ! Multi-threaded loop over remaining grid nodes. 
  ! WORKSPACE FIX: Allocations moved INSIDE the private thread scope region.
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(k, Ak_bk, Ak_loc, did_update) &
!$OMP SHARED(lkdim, nrens, nobs_ok, xobs, yobs, rho_loc, ibarcl_rst, nk_l_local, local_upd_rand, nnkn, k_start)
    allocate(Ak_bk(lkdim,nrens), Ak_loc(lkdim,nrens))
!$OMP DO SCHEDULE(static)
    do k = k_start, nnkn
       call type_to_kmat(ibarcl_rst, Ak_bk, k, lkdim, nrens)

       call locan_k(k, lkdim, nrens, nobs_ok, xobs, yobs, rho_loc, &
                    Ak_bk, local_upd_rand, Ak_loc, did_update)

       if (did_update) then
!$OMP ATOMIC
          nk_l_local = nk_l_local + 1
       end if

       call kmat_to_type(ibarcl_rst, Ak_loc, k, lkdim, nrens)
    end do
!$OMP END DO
    deallocate(Ak_bk, Ak_loc)
!$OMP END PARALLEL

  ! ============================= ELEMENT PHASE ============================
  ! Multi-threaded loop over grid elements.
  ! WORKSPACE FIX: Allocations isolated inside threads to prevent memory thrashing.
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(ne, Ae_bk, Ae_loc, did_update_e) &
!$OMP SHARED(lnedim, nrens, nobs_ok, xobs, yobs, rho_loc, ex, ey, ne_l_local, local_upd_rand, nnel)
    allocate(Ae_bk(lnedim,nrens), Ae_loc(lnedim,nrens))
!$OMP DO SCHEDULE(static)
    do ne = 1, nnel
       call type_to_emat(Ae_bk, ne, lnedim, nrens)
       call locan_e(ne, lnedim, nrens, nobs_ok, xobs, yobs, rho_loc, &
                    ex(ne), ey(ne), Ae_bk, local_upd_rand, Ae_loc, did_update_e)

       if (did_update_e) then
!$OMP ATOMIC
          ne_l_local = ne_l_local + 1
       end if

       call emat_to_type(Ae_loc, ne, lnedim, nrens)
    end do
!$OMP END DO
    deallocate(Ae_bk, Ae_loc)
!$OMP END PARALLEL

  ! -------------------------------- Cleanup --------------------------------
  deallocate(xobs, yobs, rho_loc)
  deallocate(ex, ey)
end subroutine local_analysis

!===========================================================================
! SUBROUTINE: locan_k
!
! PURPOSE:
!   Executes localized Ensemble Kalman Filter updates for an individual grid NODE (nk).
!   Identifies active stations inside the Gaspari-Cohn radius, builds
!   compacted local matrices using perturbed innovations (D1l), and invokes the core solver.
!
! MATHEMATICAL NOTES:
!   - Empirical evidence confirms the core 'analysis' routine expects the perturbed
!     innovation matrix (D1l) at the 5th argument slot. Passing absolute values causes
!     runaway state explosions.
!   - The random rotation flag (local_upd_rand) is handled via intent(inout) to allow
!     sequential generation at k=1 and deterministic recycling across subsequent OpenMP threads.
!===========================================================================
subroutine locan_k(nk, kdim, nren, no_tot, xo, yo, rhoo, &
                   Ak_bk, local_upd_rand, Ak_loc, did_update)
  use iso_fortran_env, only: dp => real64
  use m_analysis,      only : analysis
  use basin,           only : xgv, ygv
  use mod_enkf,        only : innov, D1, S, E, R, truncation, rmode, &
                              lrandrot, lsymsqrt, inflate, infmult
  implicit none

  integer,  intent(in)    :: nk, kdim, nren, no_tot
  real(dp), intent(in)    :: Ak_bk(kdim,nren)
  real(dp), intent(in)    :: xo(no_tot), yo(no_tot), rhoo(no_tot)
  logical,  intent(inout) :: local_upd_rand
  real(dp), intent(out)   :: Ak_loc(kdim,nren)
  logical,  intent(out)   :: did_update

  integer :: no, nno, no_k
  integer,  allocatable :: ido(:)
  real(dp), allocatable :: wo(:)
  real(dp) :: dist, w, xnode, ynode
  real(dp), allocatable :: innovl(:), D1l(:,:), Sl(:,:), El(:,:), Rl(:,:)
  real(dp), parameter   :: eps_la = 1.0e-4_dp   ! Gaspari-Cohn weight cutoff threshold

  Ak_loc = Ak_bk   ! Initialize with background state
  did_update = .false.

  allocate(ido(no_tot), wo(no_tot))
  ido = 0; wo = 0.0_dp; nno = 0

  ! Extract geographic anchors for the target node
  xnode = real(xgv(nk), dp)
  ynode = real(ygv(nk), dp)

  ! 1. Identify valid local observation tracks based on Euclidean grid distance
  do no = 1, no_tot
     ! Unified Cartesian Euclidean distance mapping matching SHYFEM metric projection
     dist = sqrt((xnode - xo(no))**2 + (ynode - yo(no))**2)

     call find_weight_GC(rhoo(no), dist, w)
     if (w > eps_la) then
        nno = nno + 1
        ido(nno) = no
        wo(nno)  = w
     end if
  end do

  no_k = nno

  ! 2. If valid stations reside within the active radius, compile compacted sub-matrices
  if (no_k > 0) then
     allocate(innovl(no_k), D1l(no_k,nren), Sl(no_k,nren))
     allocate(El(no_k,nren), Rl(no_k,no_k))
     Rl = 0.0_dp   ! Clean initialization for diagonal covariance

     do no = 1, no_k
        innovl(no) = innov(ido(no)) * wo(no)
        Sl(no,:)   = S(ido(no),:)   * wo(no)
        El(no,:)   = E(ido(no),:)   * wo(no) ! Taper observation perturbations for case 13
        D1l(no,:)  = D1(ido(no),:)  * wo(no) ! REALIGNED: Pass localized perturbed innovations
        Rl(no,no)  = R(ido(no), ido(no))
     end do

     ! 3. Invoke core Subspace Solver with verified sequence interfaces
     ! REALIGNED INTERFACE: Pass D1l at the 5th argument to match the verified forward setup
     call analysis(Ak_loc, Rl, El, Sl, D1l, innovl, kdim, nren, no_k, .false., &
                   truncation, rmode, lrandrot, local_upd_rand, lsymsqrt, &
                   inflate, infmult, .true.)

     ! Enforce explicit rotation toggle lockdown after the sequential first step consumes it
     if (local_upd_rand) local_upd_rand = .false.

     deallocate(innovl, D1l, Sl, El, Rl)
     did_update = .true.
  end if

  deallocate(ido, wo)
end subroutine locan_k

!===========================================================================
! SUBROUTINE: locan_e
!
! PURPOSE:
!   Executes localized state updates for a single grid ELEMENT centroid (ne).
!===========================================================================
subroutine locan_e(ne, nedim, nren, no_tot, xo, yo, rhoo, &
                   xe, ye, Ae_bk, local_upd_rand, Ae_loc, did_update)
  use iso_fortran_env, only: dp => real64
  use m_analysis,      only : analysis
  use mod_enkf,        only : innov, D1, S, E, R, truncation, rmode, &
                              lrandrot, lsymsqrt, inflate, infmult
  implicit none

  integer,  intent(in)    :: ne, nedim, nren, no_tot
  real(dp), intent(in)    :: xo(no_tot), yo(no_tot), rhoo(no_tot)
  real(dp), intent(in)    :: xe, ye
  real(dp), intent(in)    :: Ae_bk(nedim,nren)
  logical,  intent(inout) :: local_upd_rand
  real(dp), intent(out)   :: Ae_loc(nedim,nren)
  logical,  intent(out)   :: did_update

  integer :: no, nno, no_e
  integer,  allocatable :: ido(:)
  real(dp), allocatable :: wo(:)
  real(dp) :: dist, w
  real(dp), allocatable :: innovl(:), D1l(:,:), Sl(:,:), El(:,:), Rl(:,:)
  real(dp), parameter   :: eps_la = 1.0e-4_dp

  Ae_loc = Ae_bk
  did_update = .false.

  allocate(ido(no_tot), wo(no_tot))
  ido = 0; wo = 0.0_dp; nno = 0

  ! 1. Parse observations around element centroid coordinates (xe, ye)
  do no = 1, no_tot
     dist = sqrt((xe - xo(no))**2 + (ye - yo(no))**2)
     call find_weight_GC(rhoo(no), dist, w)
     if (w > eps_la) then
        nno = nno + 1
        ido(nno) = no
        wo(nno)  = w
     end if
  end do

  no_e = nno

  ! 2. Compile compacted element observation sub-matrices
  if (no_e > 0) then
     allocate(innovl(no_e), D1l(no_e,nren), Sl(no_e,nren))
     allocate(El(no_e,nren), Rl(no_e,no_e))
     Rl = 0.0_dp

     do no = 1, no_e
        innovl(no) = innov(ido(no)) * wo(no)
        Sl(no,:)   = S(ido(no),:)   * wo(no)
        El(no,:)   = E(ido(no),:)   * wo(no)
        D1l(no,:)  = D1(ido(no),:)  * wo(no) ! REALIGNED: Pass localized perturbed innovations
        Rl(no,no)  = R(ido(no), ido(no))
     end do

     ! 3. Execute localized element analysis updates
     ! REALIGNED INTERFACE: Pass D1l at the 5th argument to match the verified forward setup
     call analysis(Ae_loc, Rl, El, Sl, D1l, innovl, nedim, nren, no_e, .false., &
                   truncation, rmode, lrandrot, local_upd_rand, lsymsqrt, &
                   inflate, infmult, .true.)

     if (local_upd_rand) local_upd_rand = .false.

     deallocate(innovl, D1l, Sl, El, Rl)
     did_update = .true.
  end if

  deallocate(ido, wo)
end subroutine locan_e

!---------------------------------------------------------------------------
! SUBROUTINE: type_to_kmat
! PURPOSE: Packs localized nodal background fields (Z, T, S) for a specific 
!          node index 'k' into a flat 2D matrix layout (state-by-member).
!---------------------------------------------------------------------------
subroutine type_to_kmat(ibrcl, Ak_bk, k, kdim, nre)
  use iso_fortran_env, only: dp => real64
  use mod_ens_state,   only : Abk, nnlv
  implicit none

  integer,  intent(in)  :: ibrcl
  integer,  intent(in)  :: k, kdim, nre
  real(dp), intent(out) :: Ak_bk(kdim,nre)
  
  integer :: n

  ! Extract Sea Surface Height (SSH) at row index 1 for all members
  do n = 1, nre
     Ak_bk(1,n) = Abk(n)%z(k)
  end do
  
  ! If running under a baroclinic configuration, extract thermodynamic tracers
  if (ibrcl > 0) then
     do n = 1, nre
        ! Continuous contiguous memory block slice assignments
        Ak_bk(2:nnlv+1,           n) = Abk(n)%t(:,k)
        Ak_bk(nnlv+2:2*nnlv+1,    n) = Abk(n)%s(:,k)
     end do
  end if
end subroutine type_to_kmat

!---------------------------------------------------------------------------
! SUBROUTINE: type_to_emat
! PURPOSE: Packs localized element velocity fields (u, v vectors) for a 
!          specific element index 'ne' into a flat 2D sub-matrix layout.
!---------------------------------------------------------------------------
subroutine type_to_emat(Ae_bk, ne, nedim, nre)
  use iso_fortran_env, only: dp => real64
  use mod_ens_state,   only : Abk, nnlv
  implicit none

  integer,  intent(in)  :: ne, nedim, nre
  real(dp), intent(out) :: Ae_bk(nedim,nre)
  
  integer :: n

  do n = 1, nre
     ! Group horizontal components contiguously (Zonal u followed by Meridional v)
     Ae_bk(1:nnlv,        n) = Abk(n)%u(:,ne)
     Ae_bk(nnlv+1:2*nnlv, n) = Abk(n)%v(:,ne)
  end do
end subroutine type_to_emat

!---------------------------------------------------------------------------
! SUBROUTINE: kmat_to_type
! PURPOSE: Unpacks analyzed localized state results from 2D sub-matrices 
!          back into the global SHYFEM ensemble array structure (Aan).
!---------------------------------------------------------------------------
subroutine kmat_to_type(ibrcl, Ak_an, k, kdim, nre)
  use iso_fortran_env, only: dp => real64
  use mod_ens_state,   only : Aan, nnlv
  implicit none

  integer,  intent(in) :: ibrcl
  integer,  intent(in) :: k, kdim, nre
  real(dp), intent(in) :: Ak_an(kdim,nre)
  
  integer :: n

  ! Commit smoothed SSH trajectories to the global member records
  do n = 1, nre
     Aan(n)%z(k) = Ak_an(1,n)
  end do
  
  ! Commit smoothed Temperature and Salinity profiles if active
  if (ibrcl > 0) then
     do n = 1, nre
        Aan(n)%t(:,k) = Ak_an(2:nnlv+1,         n)
        Aan(n)%s(:,k) = Ak_an(nnlv+2:2*nnlv+1,  n)
     end do
  end if
end subroutine kmat_to_type

!---------------------------------------------------------------------------
! SUBROUTINE: emat_to_type
! PURPOSE: Unpacks analyzed localized element velocity parameters 
!          back into the global SHYFEM ensemble array structure (Aan).
!---------------------------------------------------------------------------
subroutine emat_to_type(Ae_an, ne, nedim, nre)
  use iso_fortran_env, only: dp => real64
  use mod_ens_state,   only : Aan, nnlv
  implicit none

  integer,  intent(in) :: ne, nedim, nre
  real(dp), intent(in) :: Ae_an(nedim,nre)
  
  integer :: n

  do n = 1, nre
     ! Map flat matrix chunks back to horizontal vector grid elements
     Aan(n)%u(:,ne) = Ae_an(1:nnlv,         n)
     Aan(n)%v(:,ne) = Ae_an(nnlv+1:2*nnlv,  n)
  end do
end subroutine emat_to_type
