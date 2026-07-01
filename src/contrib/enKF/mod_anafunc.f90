!===============================================================
! Module: mod_anafunc
! Purpose: Analysis utilities for EnKF in ocean applications
!   - Low-rank decompositions (SVD-based)
!   - Eigen-decompositions
!   - Mean-preserving square-root updates
!   - Exact diagonal inversion
!   - Inflation utilities
!   - Debug dumps
! Precision: Double precision throughout (dp)
!
! IMPROVEMENTS (v2.0):
!   - Adaptive damping in eigenvalue inversions
!   - Better numerical stability for ill-conditioned problems
!   - Enhanced diagnostics for debugging
!   - Generic support for sea level, T, S, velocity
!===============================================================
module mod_anafunc
  use iso_fortran_env, only : dp => real64
  implicit none
  private
  public :: lowrankE, eigC, eigsign, eigsign_safe, genX2, genX3, meanX5, X5sqrt
  public :: dumpX3, dumpX5, lowrankCinv, lowrankCee, svdS
  public :: exact_diag_inversion, inflationfactor, inflateA
  
  ! Parameters for numerical stability
  real(dp), parameter :: eps_eig = 1.0e-14_dp    ! Eigenvalue floor
  real(dp), parameter :: eps_damp = 1.0e-8_dp    ! Damping parameter
  
contains

!=====================================================================
! eigC: symmetric eigen-decomposition R = Z * diag(eig) * Z^T
! eigenvalues in ascending order (LAPACK DSYEVX).
!=====================================================================
subroutine eigC(R, nrobs, Z, eig)
  implicit none
  integer, intent(in) :: nrobs
  real(dp), intent(in) :: R(nrobs, nrobs)
  real(dp), intent(out) :: Z(nrobs, nrobs), eig(nrobs)
  real(dp) :: RR(nrobs, nrobs)
  real(dp) :: fwork(8*nrobs)
  integer :: iwork(5*nrobs)
  integer :: ifail(nrobs)
  real(dp) :: abstol, ddum
  integer :: idum, neig, ierr
  real(dp), external :: dlamch

  idum = 1
  abstol = 2.0_dp * dlamch('S')
  RR = R
  call dsyevx('V','A','U', nrobs, RR, nrobs, ddum, ddum, idum, idum, &
              abstol, neig, eig, Z, nrobs, fwork, 8*nrobs, iwork, ifail, ierr)
  if (ierr /= 0) then
    print *, 'eigC: dsyevx error = ', ierr
    stop
  end if
end subroutine eigC

!=====================================================================
! eigsign: invert dominant eigenvalues up to a target energy share
! and zero the tail; also (optionally) dumps spectrum.
!=====================================================================
subroutine eigsign(eig, nrobs, truncation)
  implicit none
  integer, intent(in) :: nrobs
  real(dp), intent(inout) :: eig(nrobs)
  real(dp), intent(in) :: truncation
  integer :: i, nrsigma
  real(dp) :: sigsum, sigsum1
  logical :: ex

  inquire(file='eigenvalues.dat', exist=ex)
  if (ex) then
    open(10, file='eigenvalues.dat', position='append')
    write(10,'(a,i5,a)') ' ZONE F=POINT, I=', nrobs, ' J=1 K=1'
    do i = 1, nrobs
      write(10,'(i3,g13.5)') i, eig(nrobs - i + 1)
    end do
    close(10)
  else
    open(10, file='eigenvalues.dat')
    write(10,*) 'TITLE = "Eigenvalues of C"'
    write(10,*) 'VARIABLES = "obs" "eigenvalues"'
    write(10,'(a,i5,a)') ' ZONE F=POINT, I=', nrobs, ' J=1 K=1'
    do i = 1, nrobs
      write(10,'(i3,g13.5)') i, eig(nrobs - i + 1)
    end do
    close(10)
  end if

  sigsum  = sum(eig)
  sigsum1 = 0.0_dp
  nrsigma = 0
  do i = nrobs, 1, -1
    if (sigsum1 / sigsum < truncation) then
      nrsigma = nrsigma + 1
      sigsum1 = sigsum1 + eig(i)
      eig(i)  = 1.0_dp / eig(i)
    else
      eig(1:i) = 0.0_dp
      exit
    end if
  end do
end subroutine eigsign

!=====================================================================
! eigsign_safe: IMPROVED version with adaptive damping and diagnostics
! 
! PURPOSE:
!   Apply eigenvalue truncation with Tikhonov damping for ill-conditioned
!   cases. Returns the number of eigenvalues retained for diagnostics.
!
! IMPROVEMENTS over eigsign:
!   - Adaptive damping: eig_inv = 1/(eig + damping*trace/n)
!   - Prevents excessive amplification when eig << trace
!   - Reports how many modes are retained
!   - Generic for sea level, T, S, velocity
!=====================================================================
subroutine eigsign_safe(eig, nobs, truncation, verbose)
  implicit none
  integer, intent(in) :: nobs
  real(dp), intent(inout) :: eig(nobs)
  real(dp), intent(in) :: truncation
  logical, intent(in) :: verbose
  
  integer :: i, n_retained
  real(dp) :: trace, damping_coeff, eig_threshold
  
  ! Compute trace for adaptive damping
  trace = sum(eig)
  if (nobs > 0) trace = trace / real(nobs, dp)
  
  ! Adaptive damping coefficient prevents ill-conditioning
  damping_coeff = max(eps_damp, 0.01_dp * trace)
  
  ! Truncate and invert with damping
  n_retained = 0
  eig_threshold = max(truncation, eps_eig)
  
  do i = 1, nobs
     if (eig(i) > eig_threshold) then
        ! Apply Tikhonov damping: avoid 1/x singularity for small x
        eig(i) = 1.0_dp / (eig(i) + damping_coeff)
        n_retained = n_retained + 1
     else
        eig(i) = 0.0_dp
     end if
  end do
  
  if (verbose) then
     write(*,'(a,i5,a,i5,a,e12.4,a,e12.4)') &
         'eigsign_safe: retained ', n_retained, ' of ', nobs, &
         ' eigenvalues; threshold=', eig_threshold, ' damping=', damping_coeff
  end if
  
end subroutine eigsign_safe

!=====================================================================
! lowrankE (SAFE: uses workspace-query SVD)
! Compute a low-rank basis W and diagonal spectrum eig based on
! SVD(S) and SVD(X0) with X0 = diag(sig0) * U0^T * E
! - Input:
!     S(nrobs,nrens), E(nrobs,nrens), truncation in (0,1]
!     nrmin = requested truncated rank (<= min(nrobs,nrens))
! - Output:
!     W(nrobs,nrmin), eig(nrmin) with eig(i) = 1/(1 + sigma_i^2)
!   (sigma_i are singular values of X0)
!=====================================================================
subroutine lowrankE(S, E, nrobs, nrens, nrmin, W, eig, truncation)
  implicit none
  integer, intent(in) :: nrobs, nrens, nrmin
  real(dp), intent(in) :: S(nrobs,nrens), E(nrobs,nrens)
  real(dp), intent(out) :: W(nrobs,nrmin), eig(nrmin)
  real(dp), intent(in) :: truncation
  ! Locals
  real(dp) :: U0(nrobs,nrmin), sig0(nrmin)
  real(dp), allocatable :: X0(:,:), Utmp(:,:), VT0(:,:), sval(:), work(:)
  integer :: i, j, ierr, lwork, minmn
  real(dp) :: wkopt

  ! 1) Truncated left singular vectors of S (safe)
  call svdS(S, nrobs, nrens, nrmin, U0, sig0, truncation)

  ! 2) X0 = diag(sig0) * U0^T * E (size: nrmin x nrens)
  allocate(X0(nrmin,nrens))
  call dgemm('T','N', nrmin, nrens, nrobs, 1.0_dp, U0, nrobs, E, nrobs, 0.0_dp, X0, nrmin)
  do j = 1, nrens
    do i = 1, nrmin
      X0(i,j) = sig0(i) * X0(i,j)
    end do
  end do

  ! 3) SVD(X0) with workspace query (Utmp: nrmin x min(nrmin,nrens))
  minmn = min(nrmin, nrens)
  allocate(Utmp(nrmin,minmn), sval(minmn), VT0(1,1))  ! VT not referenced for JOBVT='N'

  ! -- Workspace query
  lwork = -1
  allocate(work(1))
  call dgesvd('S','N', nrmin, nrens, X0, nrmin, sval, Utmp, nrmin, VT0, 1, work, lwork, ierr)
  wkopt = work(1)
  deallocate(work)
  if (ierr /= 0) then
    print *, 'mod_anafunc (lowrankE): DGESVD workspace query error = ', ierr
    stop
  end if

  ! -- Actual SVD
  lwork = max(1, int(wkopt))
  allocate(work(lwork))
  call dgesvd('S','N', nrmin, nrens, X0, nrmin, sval, Utmp, nrmin, VT0, 1, work, lwork, ierr)
  deallocate(work)
  if (ierr /= 0) then
    print *, 'mod_anafunc (lowrankE): dgesvd error = ', ierr
    stop
  end if

  ! 4) eig(i) = 1/(1 + sigma_i^2) with damping
  eig = 0.0_dp
  do i = 1, minmn
    ! IMPROVED: Add damping to prevent excessive inversion
    eig(i) = 1.0_dp / (1.0_dp + sval(i)**2 + eps_damp)
  end do

  ! 5) W = U0 * diag(sig0^{-1}) * Utmp
  !    Note: svdS returns sig0 = 1/sigma(S). To get diag(sig0^{-1}) we can
  !    multiply Utmp by (1/sig0)^{-1} = sigma(S). The original implementation
  !    multiplies by sig0 here, which corresponds to another stable form used
  !    in the legacy code path. We keep the legacy behavior (drop-in).
  do j = 1, nrmin
    do i = 1, nrmin
      Utmp(i,j) = sig0(i) * Utmp(i,j)
    end do
  end do
  call dgemm('N','N', nrobs, nrmin, nrmin, 1.0_dp, U0, nrobs, Utmp, nrmin, 0.0_dp, W, nrobs)

  deallocate(X0, Utmp, sval, VT0)
end subroutine lowrankE

!=====================================================================
! genX2: X2 = (I + Λ)^{-1/2} * W^T * S
!=====================================================================
subroutine genX2(nrens, nrobs, nrmin, S, W, eig, X2)
  implicit none
  integer,  intent(in)  :: nrens, nrobs, nrmin
  real(dp), intent(in)  :: W(nrobs, nrmin)
  real(dp), intent(in)  :: S(nrobs, nrens)
  real(dp), intent(in)  :: eig(nrmin)
  real(dp), intent(out) :: X2(nrmin, nrens)
  integer :: j

  ! 1. X2 = W^T * S
  ! This step computes the projection of the ensemble into the reduced space
  call dgemm('T','N', nrmin, nrens, nrobs, 1.0_dp, W, nrobs, S, nrobs, 0.0_dp, X2, nrmin)

  ! 2. Scale rows by sqrt of eigenvalues
  ! Parallelize across the ensemble members (columns) for better cache reuse
  !$omp parallel do private(j) shared(X2, eig, nrens, nrmin)
  do j = 1, nrens
    ! Multiply each column by the sqrt of the eigenvalues vector
    X2(:,j) = sqrt(max(0.0_dp, eig(:))) * X2(:,j)
  end do
  !$omp end parallel do
end subroutine genX2

!=====================================================================
! genX3: X3 = W * (diag(eig) * W^T * D)
!=====================================================================
subroutine genX3(nrens, nrobs, nrmin, eig, W, D, X3)
  implicit none
  integer,  intent(in)  :: nrens, nrobs, nrmin
  real(dp), intent(in)  :: eig(nrmin)
  real(dp), intent(in)  :: W(nrobs, nrmin)
  real(dp), intent(in)  :: D(nrobs, nrens)
  real(dp), intent(out) :: X3(nrobs, nrens)
  
  ! Temporary workspace matrices
  real(dp) :: X1(nrmin, nrobs)
  real(dp) :: X2(nrmin, nrens)
  integer  :: i

  ! 1. X1 = diag(eig) * W^T
  ! We scale each row of W^T (which is each column of W) by the eigenvalues
  !$omp parallel do private(i) shared(X1, eig, W, nrmin, nrobs)
  do i = 1, nrmin
    X1(i,:) = eig(i) * W(:,i)
  end do
  !$omp end parallel do

  ! 2. X2 = X1 * D  -> (nrmin x nrens)
  call dgemm('N','N', nrmin, nrens, nrobs, 1.0_dp, X1, nrmin, D, nrobs, 0.0_dp, X2, nrmin)

  ! 3. X3 = W * X2  -> (nrobs x nrens)
  ! Final weight matrix calculation
  call dgemm('N','N', nrobs, nrens, nrmin, 1.0_dp, W, nrobs, X2, nrmin, 0.0_dp, X3, nrobs)
end subroutine genX3

!=====================================================================
! meanX5: constructs mean-update matrix; adds 1/N term
!=====================================================================
subroutine meanX5(nrens, nrobs, nrmin, S, W, eig, innov, X5)
  implicit none
  integer,  intent(in)  :: nrens, nrobs, nrmin
  real(dp), intent(in)  :: W(nrobs, nrmin)
  real(dp), intent(in)  :: S(nrobs, nrens)
  real(dp), intent(in)  :: eig(nrmin)
  real(dp), intent(in)  :: innov(nrobs)
  real(dp), intent(out) :: X5(nrens, nrens)
  
  ! Local work arrays
  real(dp) :: y1(nrmin), y2(nrmin), y3(nrobs), y4(nrens)
  real(dp) :: inv_nrens
  integer  :: j

  ! --- Step 1: Project innovation into the reduced space ---
  if (nrobs == 1) then
    y1(1) = W(1,1) * innov(1)
    y2(1) = eig(1) * y1(1)
    y3(1) = W(1,1) * y2(1)
    y4(:) = y3(1) * S(1,:)
  else
    ! y1 = W' * innov
    call dgemv('T', nrobs, nrmin, 1.0_dp, W, nrobs, innov, 1, 0.0_dp, y1, 1)
    ! y2 = diag(eig) * y1
    y2 = eig * y1
    ! y3 = W * y2
    call dgemv('N', nrobs, nrmin, 1.0_dp, W, nrobs, y2, 1, 0.0_dp, y3, 1)
    ! y4 = S' * y3 (this represents the weights for each ensemble member)
    call dgemv('T', nrobs, nrens, 1.0_dp, S, nrobs, y3, 1, 0.0_dp, y4, 1)
  end if

  ! --- Step 2: Build the full transformation matrix X5 ---
  inv_nrens = 1.0_dp / real(nrens, dp)

  ! Parallelize the column assignment and the mean addition
  !$omp parallel do private(j) shared(X5, y4, inv_nrens, nrens)
  do j = 1, nrens
    ! Each column of X5 gets the weights (y4) plus the uniform mean weight
    X5(:, j) = inv_nrens + y4(:)
  end do
  !$omp end parallel do

end subroutine meanX5

!=====================================================================
! X5sqrt — SAFE DGESVD, preserves original mean‑preserving form
!=====================================================================
subroutine X5sqrt(X2, nrobs, nrens, nrmin, X5, lrandrot, lupdate_randrot, mode, lsymsqrt)
  use m_randrot
  use m_mean_preserving_rotation
  use iso_fortran_env, only : dp => real64
  implicit none

  ! --- Arguments ---
  integer,  intent(in)     :: nrobs, nrens
  integer,  intent(inout)  :: nrmin
  real(dp), intent(in)     :: X2(nrmin, nrens)
  real(dp), intent(inout)  :: X5(nrens, nrens)
  logical,  intent(in)     :: lrandrot, lupdate_randrot
  integer,  intent(in)     :: mode
  logical,  intent(in)     :: lsymsqrt

  ! --- Local variables ---
  real(dp), allocatable :: Utmp(:,:), VT(:,:), work(:), sig(:)
  real(dp), allocatable :: X3(:,:), X33(:,:), X4(:,:), X2loc(:,:)
  real(dp)              :: IenN(nrens, nrens), inv_n
  real(dp), save, allocatable :: rot(:,:)
  integer               :: i, j, ierr, lwork, minmn
  real(dp)              :: wkopt

  ! 1. Random Rotation management
  if (lrandrot .and. lupdate_randrot) then
    if (allocated(rot)) then
       if (size(rot,1) /= nrens) deallocate(rot)
    end if
    if (.not. allocated(rot)) allocate(rot(nrens, nrens))
    call mean_preserving_rotation(rot, nrens)
  end if

  if (mode == 21) nrmin = min(nrens, nrobs)
  minmn = min(nrmin, nrens)

  ! 2. Prepare SVD (dgesvd modifies the input matrix, so we use a copy)
  allocate(X2loc(nrmin, nrens)); X2loc = X2
  allocate(Utmp(nrmin, minmn), VT(nrens, nrens), sig(minmn))

  ! Workspace query for SVD
  lwork = -1
  allocate(work(1))
  call dgesvd('S', 'A', nrmin, nrens, X2loc, nrmin, sig, Utmp, nrmin, VT, nrens, work, lwork, ierr)
  lwork = max(1, int(work(1)))
  deallocate(work)
  allocate(work(lwork))

  ! Perform SVD: X2 = U * Sig * VT
  call dgesvd('S', 'A', nrmin, nrens, X2loc, nrmin, sig, Utmp, nrmin, VT, nrens, work, lwork, ierr)
  if (ierr /= 0) stop 'X5sqrt: dgesvd failed'
  deallocate(work, X2loc)

  ! 3. Compute the square root of the weights
  allocate(X3(nrens, nrens))
  !$omp parallel do private(j) shared(X3, VT, sig, minmn, nrens)
  do j = 1, nrens
    ! VT contains rows of V, so VT(j,i) is V(i,j). We need columns of V.
    ! Optimization: Transpose VT while applying isigma
    X3(:, j) = VT(j, :)
    if (j <= minmn) then
       X3(:, j) = X3(:, j) * sqrt(max(0.0_dp, 1.0_dp - sig(j)**2))
    end if
  end do
  !$omp end parallel do

  ! 4. Handle Symmetric and Random Rotations
  allocate(X33(nrens, nrens))
  if (lsymsqrt) then
    ! X33 = V * sqrt(I - Sig^2) * V^T
    call dgemm('N', 'N', nrens, nrens, nrens, 1.0_dp, X3, nrens, VT, nrens, 0.0_dp, X33, nrens)
  else
    X33 = X3
  end if

  allocate(X4(nrens, nrens))
  if (lrandrot .and. allocated(rot)) then
    call dgemm('N', 'N', nrens, nrens, nrens, 1.0_dp, X33, nrens, rot, nrens, 0.0_dp, X4, nrens)
  else
    X4 = X33
  end if

  ! 5. Mean-Preserving Projection: X5 = (I - 1/N * 11^T) * X4 + X5
  ! Note: X5 already contains the mean update from meanX5
  inv_n = -1.0_dp / real(nrens, dp)
  !$omp parallel do private(i, j) shared(IenN, inv_n, nrens)
  do j = 1, nrens
    do i = 1, nrens
      IenN(i, j) = inv_n
      if (i == j) IenN(i, j) = IenN(i, j) + 1.0_dp
    end do
  end do
  !$omp end parallel do

  ! Final update: A_new = Mean_Update + (Anomalies * X4)
  call dgemm('N', 'N', nrens, nrens, nrens, 1.0_dp, IenN, nrens, X4, nrens, 1.0_dp, X5, nrens)

  deallocate(X3, X33, X4, Utmp, VT, sig)
end subroutine X5sqrt

!=====================================================================
! dumpX3 / dumpX5 (debug I/O)
!=====================================================================
subroutine dumpX3(X3, S, nrobs, nrens)
  implicit none
  integer,  intent(in) :: nrens, nrobs
  real(dp), intent(in) :: X3(nrobs, nrens)
  real(dp), intent(in) :: S(nrobs, nrens)
  character(len=2)     :: tag2
  integer              :: u

  tag2 = 'X3'
  ! Use NEWUNIT to avoid conflicts with other open files
  open(newunit=u, file='X3.uf', form='unformatted', status='replace')
  write(u) tag2, nrens, nrobs
  write(u) X3
  write(u) S
  close(u)
end subroutine dumpX3

subroutine dumpX5(X5, nrens)
  implicit none
  integer,  intent(in) :: nrens
  real(dp), intent(in) :: X5(nrens, nrens)
  integer              :: j, u
  character(len=2)     :: tag2
  real(dp)             :: col_sum(nrens), row_sum(nrens)

  tag2 = 'X5'

  ! 1. Binary dump for the Smoother
  open(newunit=u, file='X5.uf', form='unformatted', status='replace')
  write(u) tag2, nrens
  write(u) X5
  close(u)

  ! 2. Diagnostics: Pre-calculate sums efficiently
  ! Column sums (Fast: memory contiguous)
  do j = 1, nrens
     col_sum(j) = sum(X5(:,j))
  end do

  ! Row sums (Calculated via sum intrinsic on the whole matrix for speed)
  do j = 1, nrens
     row_sum(j) = sum(X5(j,:)) / real(nrens, dp)
  end do

  ! 3. Write ASCII diagnostics
  open(newunit=u, file='X5col.dat', status='replace')
  do j = 1, nrens
    write(u,'(i5, f14.6)') j, col_sum(j)
  end do
  close(u)

  open(newunit=u, file='X5row.dat', status='replace')
  do j = 1, nrens
    write(u,'(i5, f14.6)') j, row_sum(j)
  end do
  close(u)
end subroutine dumpX5

!=====================================================================
! lowrankCinv — uses svdS (safe) and eigC with adaptive damping
!=====================================================================
subroutine lowrankCinv(S, R, nrobs, nrens, nrmin, W, eig, truncation)
  implicit none
  integer, intent(in) :: nrobs, nrens, nrmin
  real(dp), intent(in) :: S(nrobs, nrens), R(nrobs, nrobs)
  real(dp), intent(out) :: W(nrobs, nrmin), eig(nrmin)
  real(dp), intent(in) :: truncation
  real(dp) :: U0(nrobs, nrmin), sig0(nrmin)
  real(dp) :: B(nrmin, nrmin), Z(nrmin, nrmin)
  integer :: i, j

  call svdS(S, nrobs, nrens, nrmin, U0, sig0, truncation)
  call lowrankCee(B, nrmin, nrobs, nrens, R, U0, sig0)
  call eigC(B, nrmin, Z, eig)

  ! IMPROVED: Add damping in eigenvalue inversion
  ! This prevents division by very small numbers
  do i = 1, nrmin
    ! eig(i) = 1/(1 + eig(i)) with damping to regularize
    eig(i) = 1.0_dp / (1.0_dp + eig(i) + eps_damp)
  end do
  
  do j = 1, nrmin
    do i = 1, nrmin
      Z(i,j) = sig0(i) * Z(i,j)
    end do
  end do
  call dgemm('N','N', nrobs, nrmin, nrmin, 1.0_dp, U0, nrobs, Z, nrmin, 0.0_dp, W, nrobs)
end subroutine lowrankCinv

!=====================================================================
! lowrankCee (no DGESVD): B = sig0^{-1} * U0^T * R * U0 * sig0^{-1} * (nrens-1)
!=====================================================================
subroutine lowrankCee(B, nrmin, nrobs, nrens, R, U0, sig0)
  implicit none
  integer, intent(in) :: nrmin, nrobs, nrens
  real(dp), intent(inout) :: B(nrmin, nrmin)
  real(dp), intent(in) :: R(nrobs, nrobs), U0(nrobs, nrmin), sig0(nrmin)
  real(dp) :: X0(nrmin, nrobs)
  integer :: i, j

  call dgemm('T','N', nrmin, nrobs, nrobs, 1.0_dp, U0, nrobs, R, nrobs, 0.0_dp, X0, nrmin)
  call dgemm('N','N', nrmin, nrmin, nrobs, 1.0_dp, X0, nrmin, U0, nrobs, 0.0_dp, B, nrmin)

  do j = 1, nrmin
    do i = 1, nrmin
      B(i,j) = sig0(i) * B(i,j)
    end do
  end do
  do j = 1, nrmin
    do i = 1, nrmin
      B(i,j) = sig0(j) * B(i,j)
    end do
  end do
  B = real(nrens - 1, dp) * B
end subroutine lowrankCee

!=====================================================================
! svdS — SAFE SVD with workspace query (U0 and 1/sigma returned)
!   Inputs:
!     S(nrobs,nrens), nrmin
!   Outputs:
!     U0(nrobs,nrmin)  : left singular vectors (first nrmin columns, with truncation)
!     sig0(nrmin)      : 1/sigma for retained modes (zeros for truncated tail)
!   Notes:
!     - This keeps the drop-in interface (no nr_eff output).
!     - Truncation is based on cumulative energy within the first nrmin singular values.
!=====================================================================
subroutine svdS(S, nrobs, nrens, nrmin, U0, sig0, truncation)
  implicit none
  integer, intent(in) :: nrobs, nrens, nrmin
  real(dp), intent(in) :: S(nrobs,nrens)
  real(dp), intent(out) :: U0(nrobs,nrmin) ! left singular vectors (truncated)
  real(dp), intent(out) :: sig0(nrmin)     ! 1/sigma for retained modes, else 0
  real(dp), intent(in) :: truncation
  integer :: i, ierr, lwork, minmn, nrsigma
  real(dp) :: sigsum, sigsum1, wkopt
  real(dp), allocatable :: work(:), S0(:,:), Utmp(:,:), sval(:)
  real(dp) :: VT0(1,1) ! not referenced for JOBVT='N'

  if (nrobs <= 0 .or. nrens <= 0 .or. nrmin <= 0) then
    if (nrmin > 0) then
      U0   = 0.0_dp
      sig0 = 0.0_dp
    end if
    return
  end if

  minmn = min(nrobs, nrens)
  allocate(S0(nrobs,nrens), Utmp(nrobs,minmn), sval(minmn))
  S0   = S
  Utmp = 0.0_dp
  sval = 0.0_dp
  sig0 = 0.0_dp
  U0   = 0.0_dp

  ! -- Workspace query
  lwork = -1
  allocate(work(1))
  call dgesvd('S','N', nrobs, nrens, S0, nrobs, sval, Utmp, nrobs, VT0, 1, work, lwork, ierr)
  wkopt = work(1)
  deallocate(work)
  if (ierr /= 0) then
    print *, 'svdS: dgesvd (workspace query) error = ', ierr
    stop
  end if

  ! -- Actual SVD
  lwork = max(1, int(wkopt))
  allocate(work(lwork))
  call dgesvd('S','N', nrobs, nrens, S0, nrobs, sval, Utmp, nrobs, VT0, 1, work, lwork, ierr)
  deallocate(work)
  if (ierr /= 0) then
    print *, 'svdS: dgesvd error = ', ierr
    stop
  end if

  ! Energy within first nrmin singular values (sval is non-increasing)
  sigsum = 0.0_dp
  do i = 1, min(nrmin, minmn)
    sigsum = sigsum + sval(i)**2
  end do

  sigsum1 = 0.0_dp
  nrsigma = 0
  do i = 1, min(nrmin, minmn)
    if (sigsum > 0.0_dp .and. sigsum1/sigsum < truncation) then
      nrsigma = nrsigma + 1
      sigsum1 = sigsum1 + sval(i)**2
    else
      exit
    end if
  end do

  ! Copy first nrmin columns of Utmp into U0
  U0(:,1:min(nrmin,minmn)) = Utmp(:,1:min(nrmin,minmn))

  ! IMPROVED: Invert only retained singular values with damping
  ! This prevents ill-conditioning when singular values are small
  sig0 = 0.0_dp
  do i = 1, nrsigma
    if (sval(i) > eps_eig) then
      ! Add damping: sig_inv = 1/(sigma + damping*mean_sigma)
      sig0(i) = 1.0_dp / (sval(i) + eps_damp * (sigsum / real(max(1, nrsigma), dp)))
    else
      sig0(i) = 0.0_dp
    end if
  end do

  deallocate(S0, Utmp, sval)
end subroutine svdS

!=====================================================================
! exact_diag_inversion (no SVD)
!   Solves (I + 1/(N-1) * S^T S)^{-1} * (1/(N-1) * S^T D) in ensemble space
!   and maps back, producing X5 (transform for perturbations).
!=====================================================================
subroutine exact_diag_inversion(S, D, X5, nrens, nrobs)
  implicit none
  integer,  intent(in)  :: nrens, nrobs
  real(dp), intent(in)  :: S(nrobs, nrens), D(nrobs, nrens)
  real(dp), intent(out) :: X5(nrens, nrens)
  
  ! Temporary workspace arrays
  real(dp), allocatable :: SS(:,:), SD(:,:), ZSD(:,:), eig(:), Z(:,:)
  real(dp) :: n1, eig_threshold
  integer  :: i, j

  allocate(SS(nrens, nrens), SD(nrens, nrens), Z(nrens, nrens), &
           eig(nrens), ZSD(nrens, nrens))

  n1 = 1.0_dp / real(nrens - 1, dp)

  ! 1. Compute SS = (1/(N-1)) * S^T * S
  call dgemm('T','N', nrens, nrens, nrobs, n1, S, nrobs, S, nrobs, 0.0_dp, SS, nrens)
  
  ! Add identity matrix to SS (diagonal update)
  do i = 1, nrens
    SS(i,i) = SS(i,i) + 1.0_dp
  end do

  ! 2. Compute SD = (1/(N-1)) * S^T * D
  call dgemm('T','N', nrens, nrens, nrobs, n1, S, nrobs, D, nrobs, 0.0_dp, SD, nrens)

  ! 3. Eigen-decomposition of SS: SS = Z * diag(eig) * Z^T
  call eigC(SS, nrens, Z, eig)

  ! 4. Project SD into eigen-space: ZSD = Z^T * SD
  call dgemm('T','N', nrens, nrens, nrens, 1.0_dp, Z, nrens, SD, nrens, 0.0_dp, ZSD, nrens)

  ! 5. Scale by inverse eigenvalues with damping
  ! IMPROVED: Add damping to prevent ill-conditioning
  eig_threshold = max(eps_eig, eps_damp * sum(eig) / real(nrens, dp))
  
  !$omp parallel do private(j, i) shared(ZSD, eig, nrens, eig_threshold)
  do j = 1, nrens
    do i = 1, nrens
      if (eig(i) > eig_threshold) then
        ZSD(i,j) = (1.0_dp / (eig(i) + eps_damp)) * ZSD(i,j)
      else
        ZSD(i,j) = 0.0_dp
      end if
    end do
  end do
  !$omp end parallel do

  ! 6. Transform back: X5 = Z * ZSD
  call dgemm('N','N', nrens, nrens, nrens, 1.0_dp, Z, nrens, ZSD, nrens, 0.0_dp, X5, nrens)

  ! 7. Add Identity to the final transformation matrix
  ! This ensures A_new = A_old + A_old * X5 = A_old * (I + X5)
  do i = 1, nrens
    X5(i,i) = X5(i,i) + 1.0_dp
  end do

  deallocate(eig, Z, SS, SD, ZSD)
end subroutine exact_diag_inversion

!=====================================================================
! inflationfactor (no SVD)
!   Estimates a multiplicative inflation factor from random probes
!   after applying X5 to whitened random ensembles.
!=====================================================================
subroutine inflationfactor(X5, nrens, inffac)
  use iso_fortran_env, only : dp => real64
  use m_multa
  use m_random
  implicit none

  ! --- Arguments ---
  integer,  intent(in)  :: nrens          ! Number of ensemble members
  real(dp), intent(in)  :: X5(nrens, nrens) ! Transformation matrix
  real(dp), intent(out) :: inffac         ! Calculated inflation factor

  ! --- Local variables ---
  integer, parameter    :: ndim = 300     ! State dimension for synthetic ensemble
  real(dp)              :: verens(ndim, nrens)
  real(dp)              :: std_vec(ndim)
  real(dp)              :: ave, s_dev, total_std
  integer               :: i

  ! 1. Generate synthetic ensemble with random noise
  ! Note: Ensure your random routine is thread-safe or keep it outside parallel blocks
  call random(verens, ndim*nrens)

  ! 2. Center and Normalize the synthetic ensemble (Mean 0, Std 1)
  ! We parallelize across the state dimension (rows)
  !$omp parallel do private(i, ave, s_dev) shared(verens)
  do i = 1, ndim
    ! Calculate row mean
    ave = sum(verens(i,:)) / real(nrens, dp)
    verens(i,:) = verens(i,:) - ave
    
    ! Calculate standard deviation (N-biased)
    s_dev = sqrt( sum(verens(i,:)**2) / real(nrens, dp) )
    
    ! Scale to unit variance if not singular
    if (s_dev > 1.0e-14_dp) then
       verens(i,:) = verens(i,:) / s_dev
    end if
  end do
  !$omp end parallel do

  ! 3. Apply the analysis transformation X5 to the synthetic ensemble
  ! This mimics how the actual filter updates the ensemble perturbations
  call multa(verens, X5, ndim, nrens, ndim)

  ! 4. Measure the spread after transformation
  ! We re-center each row to isolate the perturbation variance
  !$omp parallel do private(i, ave) shared(verens, std_vec)
  do i = 1, ndim
    ave = sum(verens(i,:)) / real(nrens, dp)
    ! Calculate post-analysis standard deviation
    std_vec(i) = sqrt( sum((verens(i,:) - ave)**2) / real(nrens, dp) )
  end do
  !$omp end parallel do

  ! 5. Compute the final Inflation Factor
  ! Inflation is the inverse of the average spread contraction
  total_std = sum(std_vec) / real(ndim, dp)
  
  if (total_std > 1.0e-14_dp) then
     inffac = 1.0_dp / total_std
  else
     inffac = 1.0_dp ! Safety fallback
  end if

end subroutine inflationfactor


!=====================================================================
! inflateA (no SVD) — multiplicative inflation around ensemble mean
!=====================================================================
subroutine inflateA(ndim, nrens, A, inflation)
  use m_ensmean
  implicit none
  integer, intent(in) :: ndim, nrens
  real(dp), intent(inout) :: A(ndim, nrens)
  real(dp), intent(in) :: inflation(ndim)
  real(dp) :: ave(ndim)
  integer :: i, j

  call ensmean(A, ave, ndim, nrens)
  do j = 1, nrens
    do i = 1, ndim
      A(i,j) = ave(i) + (A(i,j) - ave(i)) * inflation(i)
    end do
  end do
end subroutine inflateA

end module mod_anafunc
