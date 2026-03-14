!=======================================================================
! sublocan.f90 -- Local analysis with EnKF on FEM grid (nodes & elems)
! OpenMP version with Gaspari–Cohn localization (NO MPI)
! DOUBLE PRECISION (real64 via iso_fortran_env) throughout
!
! Key features and changes:
! * Safe observation selection (no reliance on short-circuit AND)
! * Precompute element centroids (used in localization distance)
! * Initialize Rl to zero (diagonal R supported cleanly)
! * OpenMP per-thread work buffers allocated inside PARALLEL regions
! * Two-phase update: nodes (z,T,S) then elements (u,v)
! * Use iso_fortran_env real64 for all reals (dp kind)
!=======================================================================

subroutine local_analysis
  use iso_fortran_env, only: dp => real64
  use mod_ens_state           ! Abk(:), Aan(:), nnkn, nnel, nnlv, etc.
  use mod_manage_obs          ! observation pools and counters
  use mod_restart , only: ibarcl_rst
  use basin                   ! xgv(:), ygv(:), nen3v(:,:), etc.
  use mod_enkf                ! innov, D1, S, E, R, analysis, params
  implicit none

  ! ----------------------------- DECLARATIONS -----------------------------
  logical :: l_verbose
  integer :: no, nook
  real(dp), allocatable :: xobs(:), yobs(:), rho_loc(:)

  integer :: lkdim, lnedim
  real(dp), allocatable :: ex(:), ey(:)
  integer :: ne, k

  integer :: nk_l_local, ne_l_local
  logical :: local_upd_rand

  ! Per-thread workspaces (declared here, allocated inside PARALLEL):
  real(dp), allocatable :: Ak_bk(:,:), Ak_loc(:,:)
  logical :: did_update
  real(dp), allocatable :: Ae_bk(:,:), Ae_loc(:,:)
  logical :: did_update_e
  ! ------------------------------------------------------------------------

  ! ------------------------------- EXECUTION ------------------------------
  l_verbose = verbose   ! assume 'verbose' is provided by one of the used modules

  ! Local block dimensions (node and element)
  if (ibarcl_rst == 0) then
     lkdim  = 1                 ! only sea level at nodes
     lnedim = 2*nnlv            ! 3D velocities at elements (u,v)
  else
     lkdim  = 1 + 2*nnlv        ! z + T(1:nnlv) + S(1:nnlv) at nodes
     lnedim = 2*nnlv            ! u(1:nnlv), v(1:nnlv) at elements
  end if

  ! -------------------------- Precompute centroids ------------------------
  allocate(ex(nnel), ey(nnel))
  do ne = 1, nnel
     ex(ne) = ( xgv(nen3v(1,ne)) + xgv(nen3v(2,ne)) + xgv(nen3v(3,ne)) ) / 3.0_dp
     ey(ne) = ( ygv(nen3v(1,ne)) + ygv(nen3v(2,ne)) + ygv(nen3v(3,ne)) ) / 3.0_dp
  end do

  ! -------------------------- Extract observations ------------------------
  ! We assume one active obs type per call. We NEST
  ! the conditionals to avoid any reliance on non-short-circuit .AND.
  allocate(xobs(nobs_ok), yobs(nobs_ok), rho_loc(nobs_ok))
  nook = 0
  do no = 1, nobs_tot
     if (islev /= 0) then
        if (no <= n_0dlev) then
           if (o0dlev(no)%stat < 2) then
              nook = nook + 1
              xobs(nook)    = real(o0dlev(no)%x, dp)
              yobs(nook)    = real(o0dlev(no)%y, dp)
              rho_loc(nook) = real(o0dlev(no)%rhol, dp)
              if (l_verbose) write(*,*) 'LEV obs: n=',nook,' x,y,r=',xobs(nook),yobs(nook),rho_loc(nook)
           end if
        end if
     else if (istemp /= 0) then
        if (no <= n_0dtemp) then
           if (o0dtemp(no)%stat < 2) then
              nook = nook + 1
              xobs(nook)    = real(o0dtemp(no)%x, dp)
              yobs(nook)    = real(o0dtemp(no)%y, dp)
              rho_loc(nook) = real(o0dtemp(no)%rhol, dp)
              if (l_verbose) write(*,*) 'TEMP obs: n=',nook,' x,y,r=',xobs(nook),yobs(nook),rho_loc(nook)
           end if
        end if
     else if (issalt /= 0) then
        if (no <= n_0dsalt) then
           if (o0dsalt(no)%stat < 2) then
              nook = nook + 1
              xobs(nook)    = real(o0dsalt(no)%x, dp)
              yobs(nook)    = real(o0dsalt(no)%y, dp)
              rho_loc(nook) = real(o0dsalt(no)%rhol, dp)
              if (l_verbose) write(*,*) 'SALT obs: n=',nook,' x,y,r=',xobs(nook),yobs(nook),rho_loc(nook)
           end if
        end if
     else if (isvel /= 0) then
        ! Velocity-only assimilation path would go here if implemented.
     end if
  end do

  if (nook /= nobs_ok) then
     error stop 'local_analysis: mismatch in number of valid observations (nook vs nobs_ok)'
  end if

  ! --------------------------- Local counters -----------------------------
  nk_l_local = 0
  ne_l_local = 0

  ! Control flag for random-rotation update; will be flipped to .false.
  ! after the first successful local analysis.
  local_upd_rand = .true.

  ! =============================== NODE PHASE =============================
  ! Each thread processes a subset of k=1:nnkn. Per-thread private work
  ! arrays are allocated inside the PARALLEL region.
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(k, Ak_bk, Ak_loc, did_update) &
!$OMP SHARED(lkdim, nrens, nobs_ok, xobs, yobs, rho_loc, ibarcl_rst, nk_l_local, local_upd_rand, nnkn)
    allocate(Ak_bk(lkdim,nrens), Ak_loc(lkdim,nrens))
!$OMP DO SCHEDULE(static)
    do k = 1, nnkn
       ! Pack nodal variables (z, [T,S]) into Ak_bk
       call type_to_kmat(ibarcl_rst, Ak_bk, k, lkdim, nrens)

       ! Perform localized analysis for this node
       call locan_k(k, lkdim, nrens, nobs_ok, xobs, yobs, rho_loc, &
                    Ak_bk, local_upd_rand, Ak_loc, did_update)

       if (did_update) then
!$OMP ATOMIC
          nk_l_local = nk_l_local + 1
       end if

       ! Write back to analyzed state for this node
       call kmat_to_type(ibarcl_rst, Ak_loc, k, lkdim, nrens)
    end do
!$OMP END DO
    deallocate(Ak_bk, Ak_loc)
!$OMP END PARALLEL

  ! ============================= ELEMENT PHASE ============================
  ! Same pattern for velocities stored at elements. Work over ne=1:nnel.
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

  ! ------------------------------- Summary --------------------------------
  if (l_verbose) then
     write(*,*) 'Number of nodes with local analysis: ', nk_l_local
     write(*,*) 'Number of elements with local analysis: ', ne_l_local
  end if

  ! -------------------------------- Cleanup --------------------------------
  deallocate(xobs, yobs, rho_loc)
  deallocate(ex, ey)
end subroutine local_analysis

!---------------------------------------------------------------------------
!> Local analysis for a single NODE index `nk`.
!> Input Ak_bk is the background (packed); output Ak_loc is the analyzed
!> local state (packed). The routine scales innovations/anomalies by the
!> GC weight and uses a diagonal Rl (off-diagonals set to zero).
!---------------------------------------------------------------------------
subroutine locan_k(nk, kdim, nren, no_tot, xo, yo, rhoo, &
                   Ak_bk, local_upd_rand, Ak_loc, did_update)
  use iso_fortran_env, only: dp => real64
  use basin
  use mod_enkf
  implicit none

  integer, intent(in) :: nk, kdim, nren, no_tot
  real(dp), intent(in) :: Ak_bk(kdim,nren)
  real(dp), intent(in) :: xo(no_tot), yo(no_tot), rhoo(no_tot)
  logical , intent(inout) :: local_upd_rand
  real(dp), intent(out) :: Ak_loc(kdim,nren)
  logical , intent(out) :: did_update

  integer :: no, nno, no_k
  integer, allocatable :: ido(:)
  real(dp), allocatable :: wo(:)
  real(dp) :: dist, w
  real(dp), allocatable :: innovl(:), D1l(:,:), Sl(:,:), El(:,:), Rl(:,:)
  real(dp), parameter :: eps_la = 1.0e-4_dp   ! Minimum GC weight to include obs
  real(dp) :: lon_m, lat_m, rhoo_m
  integer, save :: icall = 0

  Ak_loc = Ak_bk   ! Start from background
  did_update = .false.

  allocate(ido(no_tot), wo(no_tot))
  ido = 0; wo = 0.0_dp; nno = 0

  ! Build local obs list around node position
  do no = 1, no_tot
     call deg2meters(xo(no), yo(no), rhoo(no), .false., rhoo_m)
     call deg2meters(real(xgv(nk), dp), real(ygv(nk), dp), real(xgv(nk), dp) - xo(no), .true., lon_m)
     call deg2meters(real(xgv(nk), dp), real(ygv(nk), dp), real(ygv(nk), dp) - yo(no), .false., lat_m)
     dist = sqrt( lon_m**2 + lat_m**2 )
     call find_weight_GC(rhoo_m, dist, w)
     if (w > eps_la) then
        nno = nno + 1
        ido(nno) = no
        wo(nno)  = w
     end if
  end do

  no_k = nno
  if (no_k > 0) then
     allocate(innovl(no_k), D1l(no_k,nren), Sl(no_k,nren))
     allocate(El(no_k,nren), Rl(no_k,no_k))
     Rl = 0.0_dp   ! Ensure a clean diagonal covariance

     do no = 1, no_k
        innovl(no)   = innov(ido(no)) * wo(no)
        D1l(no,:)    = D1(ido(no),:) * wo(no)
        Sl(no,:)     = S (ido(no),:) * wo(no)
        El(no,:)     = E (ido(no),:)
        Rl(no,no)    = R (ido(no), ido(no))
     end do

     call analysis(Ak_loc, Rl, El, Sl, D1l, innovl, kdim, nren, no_k, .false., &
                   truncation, rmode, lrandrot, local_upd_rand, lsymsqrt, &
                   inflate, infmult)

     deallocate(innovl, D1l, Sl, El, Rl)

!$OMP CRITICAL
     ! Flip the random-rotation update after the first successful local analysis
     if (icall == 0) then
        if (any(abs(Ak_loc - Ak_bk) > 0.0_dp)) then
           local_upd_rand = .false.
           icall = 1
        end if
     end if
!$OMP END CRITICAL

     did_update = .true.
  end if

  deallocate(ido, wo)
end subroutine locan_k

!---------------------------------------------------------------------------
!> Local analysis for a single ELEMENT index `ne`.
!> Input Ae_bk is the background (packed); output Ae_loc is the analyzed
!> local state (packed). Uses centroid (xe,ye) for localization distance.
!---------------------------------------------------------------------------
subroutine locan_e(ne, nedim, nren, no_tot, xo, yo, rhoo, &
                   xe, ye, Ae_bk, local_upd_rand, Ae_loc, did_update)
  use iso_fortran_env, only: dp => real64
  use mod_enkf
  implicit none

  integer, intent(in) :: ne, nedim, nren, no_tot
  real(dp), intent(in) :: xo(no_tot), yo(no_tot), rhoo(no_tot)
  real(dp), intent(in) :: xe, ye
  real(dp), intent(in) :: Ae_bk(nedim,nren)
  logical , intent(inout) :: local_upd_rand
  real(dp), intent(out) :: Ae_loc(nedim,nren)
  logical , intent(out) :: did_update

  integer :: no, nno, no_e
  integer, allocatable :: ido(:)
  real(dp), allocatable :: wo(:)
  real(dp) :: dist, w
  real(dp), allocatable :: innovl(:), D1l(:,:), Sl(:,:), El(:,:), Rl(:,:)
  real(dp), parameter :: eps_la = 1.0e-4_dp
  integer, save :: icall = 0

  Ae_loc = Ae_bk
  did_update = .false.

  allocate(ido(no_tot), wo(no_tot))
  ido = 0; wo = 0.0_dp; nno = 0

  do no = 1, no_tot
     dist = sqrt( (xe - xo(no))**2 + (ye - yo(no))**2 )
     call find_weight_GC(rhoo(no), dist, w)
     if (w > eps_la) then
        nno = nno + 1
        ido(nno) = no
        wo(nno)  = w
     end if
  end do

  no_e = nno
  if (no_e > 0) then
     allocate(innovl(no_e), D1l(no_e,nren), Sl(no_e,nren))
     allocate(El(no_e,nren), Rl(no_e,no_e))
     Rl = 0.0_dp

     do no = 1, no_e
        innovl(no)   = innov(ido(no)) * wo(no)
        D1l(no,:)    = D1(ido(no),:) * wo(no)
        Sl(no,:)     = S (ido(no),:) * wo(no)
        El(no,:)     = E (ido(no),:)
        Rl(no,no)    = R (ido(no), ido(no))
     end do

     call analysis(Ae_loc, Rl, El, Sl, D1l, innovl, nedim, nren, no_e, .false., &
                   truncation, rmode, lrandrot, local_upd_rand, lsymsqrt, &
                   inflate, infmult)

     deallocate(innovl, D1l, Sl, El, Rl)

!$OMP CRITICAL
     if (icall == 0) then
        if (any(abs(Ae_loc - Ae_bk) > 0.0_dp)) then
          local_upd_rand = .false.
          icall = 1
        end if
     end if
!$OMP END CRITICAL

     did_update = .true.
  end if

  deallocate(ido, wo)
end subroutine locan_e

!---------------------------------------------------------------------------
!> Pack nodal (k-th node) background into a 2D matrix (state-by-member)
!> Layout (consistent with original code):
!> Ak_bk(1,:) = z(k)
!> Ak_bk(2:nnlv+1,:) = t(:,k) [only if ibarcl>0]
!> Ak_bk(nnlv+2:2*nnlv+1,:) = s(:,k) [only if ibarcl>0]
!---------------------------------------------------------------------------
subroutine type_to_kmat(ibrcl, Ak_bk, k, kdim, nre)
  use iso_fortran_env, only: dp => real64
  use mod_ens_state
  implicit none

  integer, intent(in) :: ibrcl
  integer, intent(in) :: k, kdim, nre
  real(dp), intent(out):: Ak_bk(kdim,nre)
  integer :: n

  do n = 1, nre
     Ak_bk(1,n) = Abk(n)%z(k)
  end do
  if (ibrcl > 0) then
     do n = 1, nre
        Ak_bk(2:nnlv+1,           n) = Abk(n)%t(:,k)
        Ak_bk(nnlv+2:2*nnlv+1,    n) = Abk(n)%s(:,k)
     end do
  end if
end subroutine type_to_kmat

!---------------------------------------------------------------------------
!> Pack elemental (ne-th element) background into 2D matrix
!> Layout:
!> Ae_bk(1:nnlv,:)           = u(:,ne)
!> Ae_bk(nnlv+1:2*nnlv,:)    = v(:,ne)
!---------------------------------------------------------------------------
subroutine type_to_emat(Ae_bk, ne, nedim, nre)
  use iso_fortran_env, only: dp => real64
  use mod_ens_state
  implicit none

  integer, intent(in) :: ne, nedim, nre
  real(dp), intent(out):: Ae_bk(nedim,nre)
  integer :: n

  do n = 1, nre
     Ae_bk(1:nnlv,        n) = Abk(n)%u(:,ne)
     Ae_bk(nnlv+1:2*nnlv, n) = Abk(n)%v(:,ne)
  end do
end subroutine type_to_emat

!---------------------------------------------------------------------------
!> Unpack analyzed node back to Aan.
!---------------------------------------------------------------------------
subroutine kmat_to_type(ibrcl, Ak_an, k, kdim, nre)
  use iso_fortran_env, only: dp => real64
  use mod_ens_state
  implicit none

  integer, intent(in) :: ibrcl
  integer, intent(in) :: k, kdim, nre
  real(dp), intent(in):: Ak_an(kdim,nre)
  integer :: n

  do n = 1, nre
     Aan(n)%z(k) = Ak_an(1,n)
  end do
  if (ibrcl > 0) then
     do n = 1, nre
        Aan(n)%t(:,k) = Ak_an(2:nnlv+1,         n)
        Aan(n)%s(:,k) = Ak_an(nnlv+2:2*nnlv+1,  n)
     end do
  end if
end subroutine kmat_to_type

!---------------------------------------------------------------------------
!> Unpack analyzed element back to Aan.
!---------------------------------------------------------------------------
subroutine emat_to_type(Ae_an, ne, nedim, nre)
  use iso_fortran_env, only: dp => real64
  use mod_ens_state
  implicit none

  integer, intent(in) :: ne, nedim, nre
  real(dp), intent(in):: Ae_an(nedim,nre)
  integer :: n

  do n = 1, nre
     Aan(n)%u(:,ne) = Ae_an(1:nnlv,         n)
     Aan(n)%v(:,ne) = Ae_an(nnlv+1:2*nnlv,  n)
  end do
end subroutine emat_to_type

!=======================================================================
! End of file
!=======================================================================
