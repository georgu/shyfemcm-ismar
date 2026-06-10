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

  ! Per-thread workspaces
  real(dp), allocatable :: Ak_bk(:,:), Ak_loc(:,:)
  logical :: did_update
  real(dp), allocatable :: Ae_bk(:,:), Ae_loc(:,:)
  logical :: did_update_e
  
  integer :: k_start
  ! ------------------------------------------------------------------------

  l_verbose = .false.

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

  ! -------------------------- Dynamic Obs Counting ------------------------
  ! PASS 1: Dynamically count active observations to prevent nook vs nobs_ok mismatch
  nook = 0
  do no = 1, nobs_tot
     if (islev /= 0 .and. no <= n_0dlev) then
        if (o0dlev(no)%stat < 2) nook = nook + 1
     else if (istemp /= 0 .and. no <= n_0dtemp) then
        if (o0dtemp(no)%stat < 2) nook = nook + 1
     else if (issalt /= 0 .and. no <= n_0dsalt) then
        if (o0dsalt(no)%stat < 2) nook = nook + 1
     end if
  end do

  ! Override/Align nobs_ok with the real runtime count
  nobs_ok = nook

  ! PASS 2: Safely allocate arrays with the verified size
  allocate(xobs(nobs_ok), yobs(nobs_ok), rho_loc(nobs_ok))
  
  ! PASS 3: Populate observation arrays
  nook = 0
  do no = 1, nobs_tot
     if (islev /= 0 .and. no <= n_0dlev) then
        if (o0dlev(no)%stat < 2) then
           nook = nook + 1
           xobs(nook)    = real(o0dlev(no)%x, dp)
           yobs(nook)    = real(o0dlev(no)%y, dp)
           rho_loc(nook) = real(o0dlev(no)%rhol, dp)
        end if
     else if (istemp /= 0 .and. no <= n_0dtemp) then
        if (o0dtemp(no)%stat < 2) then
           nook = nook + 1
           xobs(nook)    = real(o0dtemp(no)%x, dp)
           yobs(nook)    = real(o0dtemp(no)%y, dp)
           rho_loc(nook) = real(o0dtemp(no)%rhol, dp)
        end if
     else if (issalt /= 0 .and. no <= n_0dsalt) then
        if (o0dsalt(no)%stat < 2) then
           nook = nook + 1
           xobs(nook)    = real(o0dsalt(no)%x, dp)
           yobs(nook)    = real(o0dsalt(no)%y, dp)
           rho_loc(nook) = real(o0dsalt(no)%rhol, dp)
        end if
     end if
  end do

  ! --------------------------- Local counters -----------------------------
  nk_l_local = 0
  ne_l_local = 0
  local_upd_rand = .true.

  ! ========================= SEQUENTIAL FIRST STEP =========================
  ! Process node k=1 sequentially to handle random rotation generation safely.
  ! This eliminates race conditions on 'local_upd_rand' inside the OpenMP region.
  k_start = 1
  if (nnkn >= 1) then
     allocate(Ak_bk(lkdim,nrens), Ak_loc(lkdim,nrens))
     call type_to_kmat(ibarcl_rst, Ak_bk, 1, lkdim, nrens)
     
     ! 'local_upd_rand' is passed with intent(inout) and will turn .false. inside
     call locan_k(1, lkdim, nrens, nobs_ok, xobs, yobs, rho_loc, &
                  Ak_bk, local_upd_rand, Ak_loc, did_update)

     if (did_update) nk_l_local = nk_l_local + 1
     call kmat_to_type(ibarcl_rst, Ak_loc, 1, lkdim, nrens)
     deallocate(Ak_bk, Ak_loc)
     k_start = 2 ! The remaining nodes will start from index 2
  end if

  ! =============================== NODE PHASE =============================
  ! Multi-threaded loop over the remaining nodes. 'local_upd_rand' is now safe as INTENT(IN)
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
  ! Multi-threaded loop over elements. 'local_upd_rand' remains read-only (.false.)
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
!> Updates the random rotation state once if local_upd_rand is .true.
!---------------------------------------------------------------------------
subroutine locan_k(nk, kdim, nren, no_tot, xo, yo, rhoo, &
                   Ak_bk, local_upd_rand, Ak_loc, did_update)
  use iso_fortran_env, only: dp => real64
  use m_analysis
  use basin
  use mod_enkf
  implicit none

  integer, intent(in) :: nk, kdim, nren, no_tot
  real(dp), intent(in) :: Ak_bk(kdim,nren)
  real(dp), intent(in) :: xo(no_tot), yo(no_tot), rhoo(no_tot)
  logical , intent(inout) :: local_upd_rand ! Updated to .false. once rotation is generated
  real(dp), intent(out) :: Ak_loc(kdim,nren)
  logical , intent(out) :: did_update

  integer :: no, nno, no_k
  integer, allocatable :: ido(:)
  real(dp), allocatable :: wo(:)
  real(dp) :: dist, w
  real(dp), allocatable :: innovl(:), D1l(:,:), Sl(:,:), El(:,:), Rl(:,:)
  real(dp), parameter :: eps_la = 1.0e-4_dp   ! Minimum GC weight to include obs
  real(dp) :: lon_m, lat_m, rhoo_m

  Ak_loc = Ak_bk   ! Initialize with background state
  did_update = .false.

  allocate(ido(no_tot), wo(no_tot))
  ido = 0; wo = 0.0_dp; nno = 0

  ! Build local observation list based on localization distance
  do no = 1, no_tot
     call deg2meters(xo(no), yo(no), rhoo(no), .false., rhoo_m)
     call deg2meters(real(xgv(nk), dp), real(ygv(nk), dp), real(xgv(nk), dp) - xo(no), .true., lon_m)
     call deg2meters(real(xgv(nk), dp), real(ygv(nk), dp), real(xgv(nk), dp) - yo(no), .false., lat_m)
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
     Rl = 0.0_dp   ! Clean initialization for diagonal covariance

     do no = 1, no_k
        innovl(no)   = innov(ido(no)) * wo(no)
        D1l(no,:)    = D1(ido(no),:) * wo(no)
        Sl(no,:)     = S (ido(no),:) * wo(no)
        El(no,:)     = E (ido(no),:)
        Rl(no,no)    = R (ido(no), ido(no))
     end do

     ! Call core EnKF analysis solver
     call analysis(Ak_loc, Rl, El, Sl, D1l, innovl, kdim, nren, no_k, .false., &
                   truncation, rmode, lrandrot, local_upd_rand, lsymsqrt, &
                   inflate, infmult, .true.)

     ! If rotation update was requested (.true.), the analysis routine has consumed it.
     ! We explicitly switch it off now to prevent downstream loops from regenerating it.
     if (local_upd_rand) local_upd_rand = .false.

     deallocate(innovl, D1l, Sl, El, Rl)
     did_update = .true.
  end if

  deallocate(ido, wo)
end subroutine locan_k

!---------------------------------------------------------------------------
!> Local analysis for a single ELEMENT index `ne`.
!---------------------------------------------------------------------------
subroutine locan_e(ne, nedim, nren, no_tot, xo, yo, rhoo, &
                   xe, ye, Ae_bk, local_upd_rand, Ae_loc, did_update)
  use iso_fortran_env, only: dp => real64
  use m_analysis
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
                   inflate, infmult, .true.)

     if (local_upd_rand) local_upd_rand = .false.

     deallocate(innovl, D1l, Sl, El, Rl)
     did_update = .true.
  end if

  deallocate(ido, wo)
end subroutine locan_e

!---------------------------------------------------------------------------
!> Pack nodal background into a 2D matrix (state-by-member)
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
!> Pack elemental background into 2D matrix
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
