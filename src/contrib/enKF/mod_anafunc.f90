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
!===============================================================
module mod_anafunc
  use iso_fortran_env, only : dp => real64
  implicit none
  private
  public :: lowrankE, eigC, eigsign, genX2, genX3, meanX5, X5sqrt
  public :: dumpX3, dumpX5, lowrankCinv, lowrankCee, svdS
  public :: exact_diag_inversion, inflationfactor, inflateA
contains

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

  ! 4) eig(i) = 1/(1 + sigma_i^2)
  eig = 0.0_dp
  do i = 1, minmn
    eig(i) = 1.0_dp / (1.0_dp + sval(i)**2)
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
! genX2: X2 = (I + Λ)^{-1/2} * W^T * S
!=====================================================================
subroutine genX2(nrens, nrobs, nrmin, S, W, eig, X2)
  implicit none
  integer, intent(in) :: nrens, nrobs, nrmin
  real(dp), intent(in) :: W(nrobs, nrmin)
  real(dp), intent(in) :: S(nrobs, nrens)
  real(dp), intent(in) :: eig(nrmin)
  real(dp), intent(out) :: X2(nrmin, nrens)
  integer :: i, j
  call dgemm('T','N', nrmin, nrens, nrobs, 1.0_dp, W, nrobs, S, nrobs, 0.0_dp, X2, nrmin)
  do j = 1, nrens
    do i = 1, nrmin
      X2(i,j) = sqrt(eig(i)) * X2(i,j)
    end do
  end do
end subroutine genX2

!=====================================================================
! genX3: X3 = W * (diag(eig) * W^T * D)
!=====================================================================
subroutine genX3(nrens, nrobs, nrmin, eig, W, D, X3)
  implicit none
  integer, intent(in) :: nrens, nrobs, nrmin
  real(dp), intent(in) :: eig(nrmin)
  real(dp), intent(in) :: W(nrobs, nrmin)
  real(dp), intent(in) :: D(nrobs, nrens)
  real(dp), intent(out) :: X3(nrobs, nrens)
  real(dp) :: X1(nrmin, nrobs)
  real(dp) :: X2(nrmin, nrens)
  integer :: i, j

  do i = 1, nrmin
    do j = 1, nrobs
      X1(i,j) = eig(i) * W(j,i)
    end do
  end do
  call dgemm('N','N', nrmin, nrens, nrobs, 1.0_dp, X1, nrmin, D, nrobs, 0.0_dp, X2, nrmin)
  call dgemm('N','N', nrobs, nrens, nrmin, 1.0_dp, W, nrobs, X2, nrmin, 0.0_dp, X3, nrobs)
end subroutine genX3

!=====================================================================
! meanX5: constructs mean-update matrix; adds 1/N term
!=====================================================================
subroutine meanX5(nrens, nrobs, nrmin, S, W, eig, innov, X5)
  implicit none
  integer, intent(in) :: nrens, nrobs, nrmin
  real(dp), intent(in) :: W(nrobs, nrmin)
  real(dp), intent(in) :: S(nrobs, nrens)
  real(dp), intent(in) :: eig(nrmin)
  real(dp), intent(in) :: innov(nrobs)
  real(dp), intent(out) :: X5(nrens, nrens)
  real(dp) :: y1(nrmin), y2(nrmin), y3(nrobs), y4(nrens)
  integer :: i

  if (nrobs == 1) then
    y1(1) = W(1,1) * innov(1)
    y2(1) = eig(1) * y1(1)
    y3(1) = W(1,1) * y2(1)
    y4(:) = y3(1) * S(1,:)
  else
    call dgemv('T', nrobs, nrmin, 1.0_dp, W, nrobs, innov, 1, 0.0_dp, y1, 1)
    y2 = eig * y1
    call dgemv('N', nrobs, nrmin, 1.0_dp, W, nrobs, y2, 1, 0.0_dp, y3, 1)
    call dgemv('T', nrobs, nrens, 1.0_dp, S, nrobs, y3, 1, 0.0_dp, y4, 1)
  end if

  do i = 1, nrens
    X5(:,i) = y4(:)
  end do
  X5 = 1.0_dp / real(nrens, dp) + X5
end subroutine meanX5

!=====================================================================
! X5sqrt — SAFE DGESVD, preserves original mean‑preserving form
!=====================================================================
subroutine X5sqrt(X2, nrobs, nrens, nrmin, X5, lrandrot, lupdate_randrot, mode, lsymsqrt)
  use m_randrot
  use m_mean_preserving_rotation
  implicit none
  integer, intent(in) :: nrobs, nrens
  integer, intent(inout) :: nrmin
  real(dp), intent(in) :: X2(nrmin, nrens)
  real(dp), intent(inout):: X5(nrens, nrens)
  logical, intent(in) :: lrandrot, lupdate_randrot
  integer, intent(in) :: mode
  logical, intent(in) :: lsymsqrt

  real(dp), allocatable :: Utmp(:,:), VT(:,:), work(:), isigma(:)
  real(dp), allocatable :: X3(:,:), X33(:,:), X4(:,:), X2loc(:,:)
  real(dp), allocatable :: sig(:)      ! <-- allocatable (fixed)
  real(dp) :: IenN(nrens, nrens)
  real(dp), save, allocatable :: rot(:,:)
  integer :: i, j, ierr, lwork, minmn
  real(dp) :: wkopt

  if (lrandrot .and. lupdate_randrot) then
    if (allocated(rot)) deallocate(rot)
    allocate(rot(nrens, nrens))
    call mean_preserving_rotation(rot, nrens)
  end if

  if (mode == 21) nrmin = min(nrens, nrobs)

  ! SVD of X2 (nrmin x nrens)
  minmn = min(nrmin, nrens)
  allocate(X2loc(nrmin,nrens)); X2loc = X2
  allocate(Utmp(nrmin,minmn))
  allocate(VT(nrens,nrens))
  allocate(sig(minmn))

  ! -- Workspace query
  lwork = -1
  allocate(work(1))
  call dgesvd('S','A', nrmin, nrens, X2loc, nrmin, sig, Utmp, nrmin, VT, nrens, work, lwork, ierr)
  wkopt = work(1)
  deallocate(work)
  if (ierr /= 0) then
    print *, 'X5sqrt: dgesvd (workspace query) error = ', ierr
    stop
  end if

  ! -- Actual SVD
  lwork = max(1, int(wkopt))
  allocate(work(lwork))
  call dgesvd('S','A', nrmin, nrens, X2loc, nrmin, sig, Utmp, nrmin, VT, nrens, work, lwork, ierr)
  deallocate(work, X2loc)
  if (ierr /= 0) then
    print *, 'X5sqrt: dgesvd error = ', ierr
    stop
  end if

  allocate(isigma(minmn), X3(nrens,nrens))
  do j = 1, nrens
    X3(:,j) = VT(j,:)  ! columns of V
  end do
  do j = 1, minmn
    isigma(j) = sqrt( max(0.0_dp, 1.0_dp - sig(j)**2) )
    X3(:,j) = X3(:,j) * isigma(j)
  end do

  allocate(X33(nrens,nrens))
  if (lsymsqrt) then
    call dgemm('N','N', nrens, nrens, nrens, 1.0_dp, X3, nrens, VT, nrens, 0.0_dp, X33, nrens)
  else
    X33 = X3
  end if

  allocate(X4(nrens,nrens))
  if (lrandrot) then
    call dgemm('N','N', nrens, nrens, nrens, 1.0_dp, X33, nrens, rot, nrens, 0.0_dp, X4, nrens)
  else
    X4 = X33
  end if

  ! Project to zero-mean subspace: (I - 1/N 11^T)
  IenN = -1.0_dp / real(nrens, dp)
  do i = 1, nrens
    IenN(i,i) = IenN(i,i) + 1.0_dp
  end do
  call dgemm('N','N', nrens, nrens, nrens, 1.0_dp, IenN, nrens, X4, nrens, 1.0_dp, X5, nrens)

  deallocate(isigma, X3, X33, X4, Utmp, VT, sig)
end subroutine X5sqrt

!=====================================================================
! dumpX3 / dumpX5 (debug I/O)
!=====================================================================
subroutine dumpX3(X3, S, nrobs, nrens)
  implicit none
  integer, intent(in) :: nrens, nrobs
  real(dp), intent(in) :: X3(nrens, nrens)
  real(dp), intent(in) :: S(nrobs, nrens)
  character(len=2) :: tag2
  tag2 = 'X3'
  open(10, file='X5.uf', form='unformatted')
  write(10) tag2, nrens, nrobs, X3, S
  close(10)
end subroutine dumpX3

subroutine dumpX5(X5, nrens)
  implicit none
  integer, intent(in) :: nrens
  real(dp), intent(in) :: X5(nrens, nrens)
  integer :: j
  character(len=2) :: tag2
  tag2 = 'X5'
  open(10, file='X5.uf', form='unformatted')
  write(10) tag2, nrens, X5
  close(10)

  open(10, file='X5col.dat')
  do j = 1, nrens
    write(10,'(i5,f10.4)') j, sum(X5(:,j))
  end do
  close(10)

  open(10, file='X5row.dat')
  do j = 1, nrens
    write(10,'(i5,f10.4)') j, sum(X5(j,:)) / real(nrens, dp)
  end do
  close(10)
end subroutine dumpX5

!=====================================================================
! lowrankCinv — uses svdS (safe) and eigC
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

  do i = 1, nrmin
    eig(i) = 1.0_dp / (1.0_dp + eig(i))
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

  ! Invert only retained singular values
  sig0 = 0.0_dp
  do i = 1, nrsigma
    if (sval(i) > 0.0_dp) then
      sig0(i) = 1.0_dp / sval(i)
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
  integer, intent(in) :: nrens, nrobs
  real(dp), intent(in) :: S(nrobs, nrens), D(nrobs, nrens)
  real(dp), intent(out) :: X5(nrens, nrens)
  real(dp), allocatable :: SS(:,:), SD(:,:), ZSD(:,:), eig(:), Z(:,:)
  real(dp) :: n1
  integer :: i, j

  allocate(SS(nrens, nrens), SD(nrens, nrens), Z(nrens, nrens), eig(nrens), ZSD(nrens, nrens))
  n1 = 1.0_dp / real(nrens - 1, dp)
  call dgemm('T','N', nrens, nrens, nrobs, n1, S, nrobs, S, nrobs, 0.0_dp, SS, nrens)
  do i = 1, nrens
    SS(i,i) = SS(i,i) + 1.0_dp
  end do
  call dgemm('T','N', nrens, nrens, nrobs, n1, S, nrobs, D, nrobs, 0.0_dp, SD, nrens)
  call eigC(SS, nrens, Z, eig)
  call dgemm('T','N', nrens, nrens, nrens, 1.0_dp, Z, nrens, SD, nrens, 0.0_dp, ZSD, nrens)
  do j = 1, nrens
    do i = 1, nrens
      ZSD(i,j) = (1.0_dp / eig(i)) * ZSD(i,j)
    end do
  end do
  call dgemm('N','N', nrens, nrens, nrens, 1.0_dp, Z, nrens, ZSD, nrens, 0.0_dp, X5, nrens)
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
  use m_multa
  use m_random
  implicit none
  integer, intent(in) :: nrens
  real(dp), intent(in) :: X5(nrens, nrens)
  real(dp), intent(out) :: inffac
  integer, parameter :: ndim = 300
  real(dp) :: aveverens, stdverens
  real(dp) :: verens(ndim, nrens), std(ndim)
  integer :: i, j

  call random(verens, ndim*nrens) ! assumes double-precision compatible RNG

  do i = 1, ndim
    aveverens = sum(verens(i,1:nrens)) / real(nrens, dp)
    verens(i,:) = verens(i,:) - aveverens
  end do
  do i = 1, ndim
    stdverens = sqrt( sum(verens(i,:)**2) / real(nrens, dp) )
    verens(i,:) = verens(i,:) / stdverens
  end do

  call multa(verens, X5, ndim, nrens, ndim)

  do i = 1, ndim
    aveverens = sum(verens(i,1:nrens)) / real(nrens, dp)
    verens(i,:) = verens(i,:) - aveverens
  end do
  do i = 1, ndim
    std(i) = sqrt( sum(verens(i,:)**2) / real(nrens, dp) )
  end do

  stdverens = sum(std) / real(ndim, dp)
  inffac = 1.0_dp / stdverens
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
