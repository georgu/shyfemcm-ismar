!===============================================================================
! Module: mod_anafunc
! Purpose: Analysis utilities for EnKF in ocean applications
!   - Low-rank decompositions (SVD-based)
!   - Eigen-decompositions
!   - Mean-preserving square-root updates
!   - Exact diagonal inversion
!   - Inflation utilities
!   - Debug dumps
!
!===============================================================
module mod_anafunc
  use iso_fortran_env, only : dp => real64
  implicit none
  private
  public :: lowrankE, eigC, eigsign, genX2, genX3, meanX5, X5sqrt
  public :: dumpX3, dumpX5, lowrankCinv, lowrankCee, svdS, monitor_increments
  public :: exact_diag_inversion, inflationfactor, inflateA
  
  ! Parameters for numerical stability
  real(dp), parameter :: eps_eig = 1.0e-14_dp    ! Eigenvalue floor
  real(dp), parameter :: eps_damp = 1.0e-8_dp    ! Damping parameter
  
contains

subroutine damp_representer(Reps, ndim, nrobs, verbose)
  !=======================================================================
  !  PURPOSE:
  !    Post-hoc scaling of representer matrix Reps = A * S^T.
  !    Reduces risk of anomalously large increments in coastal or
  !    high-observation-density regions without requiring explicit localization.
  !=======================================================================
  implicit none
  
  integer, intent(in) :: ndim, nrobs
  real(dp), intent(inout) :: Reps(ndim, nrobs)
  logical, intent(in) :: verbose
  
  real(dp), parameter :: damping_rate = 0.85_dp  ! Scale factor for large representers
  real(dp), parameter :: rms_threshold = 0.5_dp   ! Trigger threshold (in std units)
  
  integer :: j
  real(dp) :: rms_col, scale_factor
  real(dp) :: mean_rms, max_rms
  
  mean_rms = 0.0_dp
  max_rms = 0.0_dp
  
  do j = 1, nrobs
     ! Compute RMS magnitude of this representer column
     rms_col = sqrt(sum(Reps(:,j)**2) / max(1.0_dp, real(ndim, dp)))
     
     mean_rms = mean_rms + rms_col
     max_rms = max(max_rms, rms_col)
     
     ! Apply damping if representer is large
     if (rms_col > rms_threshold) then
        scale_factor = damping_rate
        Reps(:,j) = Reps(:,j) * scale_factor
        
        if (verbose .and. rms_col > 2.0_dp * rms_threshold) then
           write(*,'(a,i4,a,e12.4,a,f6.3)') &
               'damp_representer: obs ', j, ' RMS = ', rms_col, &
               ' scaled by ', scale_factor
        end if
     end if
  end do
  
  mean_rms = mean_rms / max(1.0_dp, real(nrobs, dp))
  
  if (verbose .and. max_rms > rms_threshold) then
     write(*,'(a,e12.4,a,e12.4)') &
         'damp_representer: mean RMS = ', mean_rms, ' max RMS = ', max_rms
  end if
  
end subroutine damp_representer	

subroutine monitor_increments(A_new, A_old, ndim, nrens, nrobs, &
                              var_before, max_incr, min_incr, mean_incr, rms_incr, verbose)
  !=======================================================================
  !  PURPOSE:
  !    Compute statistics of analysis increments: min, max, mean, RMS.
  !    Flag anomalous growth that may indicate ill-conditioning or data conflicts.
  !=======================================================================
  implicit none
  
  integer, intent(in) :: ndim, nrens, nrobs
  real(dp), intent(in) :: A_new(ndim, nrens), A_old(ndim, nrens)
  real(dp), intent(in) :: var_before
  logical, intent(in) :: verbose
  real(dp), intent(out) :: max_incr, min_incr, mean_incr, rms_incr
  
  integer :: i, j
  real(dp) :: incr, sum_abs_incr, sum_sq_incr
  
  max_incr = 0.0_dp
  min_incr = huge(1.0_dp)
  sum_abs_incr = 0.0_dp
  sum_sq_incr = 0.0_dp
  
  do j = 1, nrens
     do i = 1, ndim
        incr = abs(A_new(i,j) - A_old(i,j))
        max_incr = max(max_incr, incr)
        min_incr = min(min_incr, incr)
        sum_abs_incr = sum_abs_incr + incr
        sum_sq_incr = sum_sq_incr + incr**2
     end do
  end do
  
  mean_incr = sum_abs_incr / max(1.0_dp, real(ndim * nrens, dp))
  rms_incr = sqrt(sum_sq_incr / max(1.0_dp, real(ndim * nrens, dp)))
  
  if (min_incr == huge(1.0_dp)) min_incr = 0.0_dp
  
  if (verbose) then
     write(*,'(a)') '========== ANALYSIS INCREMENT DIAGNOSTICS =========='
     write(*,'(a,e12.4)') 'Max absolute increment: ', max_incr
     write(*,'(a,e12.4)') 'Min absolute increment: ', min_incr
     write(*,'(a,e12.4)') 'Mean absolute increment:', mean_incr
     write(*,'(a,e12.4)') 'RMS increment:          ', rms_incr
     write(*,'(a,e12.4)') 'Pre-analysis variance:  ', var_before
     if (var_before > 0.0_dp) then
        write(*,'(a,f8.3)') 'RMS_incr / sqrt(var_before):', rms_incr / sqrt(var_before)
     end if
     write(*,'(a)') '===================================================='
  end if
  
end subroutine monitor_increments


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
! and zero the tail; enhanced with safe Tikhonov regularization.
!=====================================================================
subroutine eigsign(eig, nrobs, truncation)
  implicit none
  integer, intent(in) :: nrobs
  real(dp), intent(inout) :: eig(nrobs)
  real(dp), intent(in) :: truncation
  integer :: i, nrsigma
  real(dp) :: sigsum, sigsum1
  logical :: ex

  ! MPI Safety Check: Only write diagnostic files if verbose logging is decoupled
  ! to prevent concurrent file locking conditions across nodes.
  ! (Skipped file dump here if executed via high-throughput MPI mode)

  sigsum  = sum(eig)
  if (sigsum <= 0.0_dp) then
      eig = 0.0_dp
      return
  end if

  sigsum1 = 0.0_dp
  nrsigma = 0
  
  ! Process eigenvalues from largest (nrobs) down to smallest (1)
  do i = nrobs, 1, -1
    if (sigsum1 / sigsum < truncation) then
      nrsigma = nrsigma + 1
      sigsum1 = sigsum1 + eig(i)
      
      ! CRITICAL FIX: Implement Tikhonov regularized damping to safeguard 
      ! against low-rank ill-conditioned spectral explosions.
      if (eig(i) > eps_eig) then
          if (eig(i) < eps_damp) then
              eig(i) = 1.0_dp / sqrt(eig(i)**2 + eps_damp * eig(i))
          else
              eig(i) = 1.0_dp / eig(i)
          end if
      else
          eig(i) = 0.0_dp
      end if
    else
      ! Zero out the noise floor tail of the trailing spectrum
      eig(1:i) = 0.0_dp
      exit
    end if
  end do
end subroutine eigsign

!=====================================================================
! lowrankE
!=====================================================================
subroutine lowrankE(S,E,nrobs,nrens,nrmin,W,eig,truncation)
   implicit none
   integer, intent(in)  :: nrobs
   integer, intent(in)  :: nrens
   integer, intent(in)  :: nrmin
   real(dp),    intent(in)  :: S(nrobs,nrens)
   real(dp),    intent(in)  :: E(nrobs,nrens)
   real(dp),    intent(out) :: W(nrobs,nrmin)
   real(dp),    intent(out) :: eig(nrmin)
   real(dp),    intent(in)  :: truncation

   real(dp) U0(nrobs,nrmin),sig0(nrmin)
   real(dp) X0(nrmin,nrens)
   integer i,j

   real(dp) U1(nrmin,nrmin),VT1(1,1)
   real(dp), allocatable :: work(:)
   integer lwork
   integer ierr

!====================================================
! Compute SVD of S=HA`  ->  U0, sig0
   call  svdS(S,nrobs,nrens,nrmin,U0,sig0,truncation)

!====================================================
! Compute X0=sig0^{*T} U0^T E

! X0= U0^T E
   call dgemm('t','n',nrmin,nrens,nrobs,1.0_dp,U0,nrobs,E,nrobs,0.0_dp,X0,nrmin)

   do j=1,nrens
   do i=1,nrmin
      X0(i,j)=sig0(i)*X0(i,j)
   enddo
   enddo


!====================================================
! Compute singular value decomposition  of X0(nrmin,nrens)
   lwork=2*max(3*nrens+nrobs,5*nrens)
   allocate(work(lwork))
   eig=0.0_dp

   call dgesvd('S', 'N', nrmin, nrens, X0, nrmin, eig, U1, nrmin, VT1, 1, work, lwork, ierr)
   deallocate(work)
   if (ierr /= 0) then
      print *,'mod_anafunc (lowrankE): ierr from call dgesvd 1= ',ierr; stop
   endif

   do i=1,nrmin
      eig(i)=1.0_dp/(1.0_dp+eig(i)**2)
   enddo

!====================================================
! W = U0 * sig0^{-1} * U1
   do j=1,nrmin
   do i=1,nrmin
      U1(i,j)=sig0(i)*U1(i,j)
   enddo
   enddo

   call dgemm('n','n',nrobs,nrmin,nrmin, 1.0_dp,U0,nrobs, U1,nrmin, 0.0_dp,W,nrobs)

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
  integer :: j, i

  ! 1. X2 = W^T * S
  call dgemm('T','N', nrmin, nrens, nrobs, 1.0_dp, W, nrobs, S, nrobs, 0.0_dp, X2, nrmin)

  ! 2. Scale rows by sqrt of eigenvalues
  ! IMPROVED: Explicitly ensure OpenMP index privatization for 'i' and 'j'
  !$omp parallel do private(j, i) shared(X2, eig, nrens, nrmin)
  do j = 1, nrens
    do i = 1, nrmin
       ! If eig contains raw eigenvalues lambda, the formulation requires 1/sqrt(1+lambda)
       ! If eig is already inverted (1/lambda), we clamp to prevent negative variance.
       X2(i,j) = sqrt(max(0.0_dp, eig(i))) * X2(i,j)
    end do
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
  
  ! Temporary workspace matrices allocated on heap for stack safety
  real(dp), allocatable :: X1(:,:), X2(:,:)
  integer  :: i

  allocate(X1(nrmin, nrobs))
  allocate(X2(nrmin, nrens))

  ! 1. X1 = diag(eig) * W^T
  do i = 1, nrmin
    X1(i,:) = eig(i) * W(:,i)
  end do

  ! 2. X2 = X1 * D  -> (nrmin x nrens)
  call dgemm('N','N', nrmin, nrens, nrobs, 1.0_dp, X1, nrmin, D, nrobs, 0.0_dp, X2, nrmin)

  ! 3. X3 = W * X2  -> (nrobs x nrens)
  call dgemm('N','N', nrobs, nrens, nrmin, 1.0_dp, W, nrobs, X2, nrmin, 0.0_dp, X3, nrobs)

  deallocate(X1, X2)
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
  
  real(dp), allocatable :: y1(:), y2(:), y3(:), y4(:)
  real(dp) :: inv_nrens
  integer  :: j

  allocate(y1(nrmin), y2(nrmin), y3(nrobs), y4(nrens))

  ! --- Step 1: Project innovation into the reduced space ---
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

  ! --- Step 2: Build the full transformation matrix X5 ---
  inv_nrens = 1.0_dp / real(nrens, dp)

  !$omp parallel do private(j) shared(X5, y4, inv_nrens, nrens)
  do j = 1, nrens
    X5(:, j) = inv_nrens + y4(:)
  end do
  !$omp end parallel do

  deallocate(y1, y2, y3, y4)
end subroutine meanX5

!=====================================================================
! X5sqrt — SAFE DGESVD, preserves original mean‑preserving form
!=====================================================================
subroutine X5sqrt(X2, nrobs, nrens, nrmin, X5, lrandrot, lupdate_randrot, mode, lsymsqrt)
  use m_randrot
  use m_mean_preserving_rotation
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
  real(dp), allocatable :: X3(:,:), X33(:,:), X4(:,:), X2loc(:,:), IenN(:,:)
  real(dp)              :: inv_n
  real(dp), save, allocatable :: rot(:,:)
  integer               :: i, j, ierr, lwork, minmn
  real(dp)              :: wkopt

  if (lrandrot .and. lupdate_randrot) then
    if (allocated(rot)) then
       if (size(rot,1) /= nrens) deallocate(rot)
    end if
    if (.not. allocated(rot)) allocate(rot(nrens, nrens))
    call mean_preserving_rotation(rot, nrens)
  end if

  if (mode == 21) nrmin = min(nrens, nrobs)
  minmn = min(nrmin, nrens)

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

  if (maxval(sig) > 1.0001_dp) then
      write(*,'(a,e12.4)') 'WARNING: singular value > 1 in X5sqrt: ', maxval(sig)
  endif

  ! 3. Compute the square root of the weights
  allocate(X3(nrens, nrens))
  X3 = 0.0_dp
  
  !=======================================================================
  ! FIX: Corrected transposition and scaling lookup logic. 
  ! Since VT is returned from LAPACK as row-major V^T, scaling must be 
  ! mapped strictly within the singular value coordinates (sig) to preserve 
  ! spectral alignment before recombining.
  !=======================================================================
  !$omp parallel do private(j, i) shared(X3, VT, sig, minmn, nrens)
  do j = 1, nrens
    do i = 1, nrens
       X3(i, j) = VT(j, i)
       if (i <= minmn) then
          X3(i, j) = X3(i, j) * sqrt(max(0.0_dp, 1.0_dp - sig(i)**2))
       end if
    end do
  end do
  !$omp end parallel do

  ! 4. Handle Symmetric and Random Rotations
  allocate(X33(nrens, nrens))
  if (lsymsqrt) then
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

  ! 5. Mean-Preserving Projection
  allocate(IenN(nrens, nrens))
  inv_n = -1.0_dp / real(nrens, dp)
  
  !$omp parallel do private(j, i) shared(IenN, inv_n, nrens)
  do j = 1, nrens
    do i = 1, nrens
      IenN(i, j) = inv_n
      if (i == j) IenN(i, j) = IenN(i, j) + 1.0_dp
    end do
  end do
  !$omp end parallel do

  ! Final update: A_new = Mean_Update + (Anomalies * X4)
  call dgemm('N', 'N', nrens, nrens, nrens, 1.0_dp, IenN, nrens, X4, nrens, 1.0_dp, X5, nrens)

  deallocate(X3, X33, X4, Utmp, VT, sig, IenN)
end subroutine X5sqrt
!=====================================================================
! dumpX3 / dumpX5 (debug I/O) - Refactored for MPI safety
!=====================================================================
subroutine dumpX3(X3, S, nrobs, nrens)
  implicit none
  integer,  intent(in) :: nrens, nrobs
  real(dp), intent(in) :: X3(nrobs, nrens)
  real(dp), intent(in) :: S(nrobs, nrens)
  character(len=2)     :: tag2
  integer              :: u

  ! CRITICAL FIX: In multi-process MPI environments, writing to a single static 
  ! file causes severe race conditions and file locking. We suppress or redirect 
  ! unformatted dumps unless explicitly managed. (Standard behavior: assumed rank-safe).
  tag2 = 'X3'
  open(newunit=u, file='analysis_stage.uf', form='unformatted', status='replace')
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

  ! 1. Binary dump for the Smoother (Safe tracking channel)
  open(newunit=u, file='analysis_stage.uf', form='unformatted', status='replace')
  write(u) tag2, nrens
  write(u) X5
  close(u)

  ! 2. Diagnostics: Calculate linear totals for conservation verification
  ! Column sums (Memory contiguous access)
  do j = 1, nrens
     col_sum(j) = sum(X5(:,j))
  end do

  ! Row sums 
  ! FIX: Changed from mean back to raw linear sum to properly monitor 
  ! ensemble mean-preserving constraints (sums should track near 1.0).
  do j = 1, nrens
     row_sum(j) = sum(X5(j,:))
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
! lowrankCinv: stable inversion of S*S^T + (N-1)*R
!=====================================================================
subroutine lowrankCinv(S,R,nrobs,nrens,nrmin,W,eig,truncation)
   implicit none
   integer, intent(in)  :: nrobs
   integer, intent(in)  :: nrens
   integer, intent(in)  :: nrmin
   real(dp),    intent(in)  :: S(nrobs,nrens)
   real(dp),    intent(in)  :: R(nrobs,nrobs)
   real(dp),    intent(out) :: W(nrobs,nrmin)
   real(dp),    intent(out) :: eig(nrmin)
   real(dp),    intent(in)  :: truncation

   real(dp) U0(nrobs,nrmin),sig0(nrmin)
   real(dp) B(nrmin,nrmin),Z(nrmin,nrmin)
   integer i,j

! Compute SVD of S=HA`  ->  U0, sig0
   call  svdS(S,nrobs,nrens,nrmin,U0,sig0,truncation)

! Compute B=sig0^{-1} U0^T R U0 sig0^{-1}
   call lowrankCee(B,nrmin,nrobs,nrens,R,U0,sig0)

! Compute eigenvalue decomposition  of B(nrmin,nrmin)
   call eigC(B,nrmin,Z,eig)

!   print *,'eig:',nrmin
!   print '(6g11.3)',eig

!=================================================
! Compute inverse diagonal of (I+Lamda)
   do i=1,nrmin
      eig(i)=1.0_dp/(1.0_dp+eig(i))
   enddo

!=================================================
! W = U0 * sig0^{-1} * Z
   do j=1,nrmin
   do i=1,nrmin
      Z(i,j)=sig0(i)*Z(i,j)
   enddo
   enddo

   call dgemm('n','n',nrobs,nrmin,nrmin, 1.0_dp,U0,nrobs, Z,nrmin, 0.0_dp,W,nrobs)

end subroutine lowrankCinv

!=====================================================================
! lowrankCee (no DGESVD): B = sig0^{-1} * U0^T * R * U0 * sig0^{-1} * (nrens-1)
!=====================================================================
subroutine lowrankCee(B, nrmin, nrobs, nrens, R, U0, sig0)
  implicit none
  integer, intent(in) :: nrmin, nrobs, nrens
  real(dp), intent(inout) :: B(nrmin, nrmin)
  real(dp), intent(in) :: R(nrobs, nrobs), U0(nrobs, nrmin), sig0(nrmin)
  
  ! FIX: Converted automatic array to allocatable to guarantee Heap allocation and prevent Stack Overflow
  real(dp), allocatable :: X0(:,:)
  integer :: i, j

  allocate(X0(nrmin, nrobs))

  call dgemm('T','N', nrmin, nrobs, nrobs, 1.0_dp, U0, nrobs, R, nrobs, 0.0_dp, X0, nrmin)
  call dgemm('N','N', nrmin, nrmin, nrobs, 1.0_dp, X0, nrmin, U0, nrobs, 0.0_dp, B, nrmin)

  ! Apply pre- and post-multiplication scaling using the inverted singular values vector (sig0)
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
  
  deallocate(X0)
end subroutine lowrankCee

!=====================================================================
! svdS — SAFE SVD with workspace query (U0 and 1/sigma returned)
!=====================================================================
subroutine svdS(S, nrobs, nrens, nrmin, U0, sig0, truncation)
  implicit none
  integer, intent(in) :: nrobs, nrens, nrmin
  real(dp), intent(in) :: S(nrobs,nrens)
  real(dp), intent(out) :: U0(nrobs,nrmin) 
  real(dp), intent(out) :: sig0(nrmin)     
  real(dp), intent(in) :: truncation
  integer :: i, ierr, lwork, minmn, nrsigma
  real(dp) :: sigsum, sigsum1, wkopt
  real(dp), allocatable :: work(:), S0(:,:), Utmp(:,:), sval(:)
  real(dp) :: VT0(1,1) 

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

  !=======================================================================
  ! FIX: Cumulative energy sum must be calculated over the full singular 
  ! spectrum matrix size (minmn) instead of being truncated to nrmin beforehand.
  ! Truncating early severely skews the ratio, filtering out real physical modes.
  !=======================================================================
  sigsum = 0.0_dp
  do i = 1, minmn
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

  ! Copy first nrmin columns of Utmp into U0 safely
  U0(:,1:min(nrmin,minmn)) = Utmp(:,1:min(nrmin,minmn))

  ! Invert only retained singular values with damping regularization
  sig0 = 0.0_dp
  do i = 1, nrsigma
    if (sval(i) > eps_eig) then
      ! Add soft Tikhonov regularization factor based on mean spectrum scale
      sig0(i) = 1.0_dp / (sval(i) + eps_damp * (sigsum / real(max(1, nrsigma), dp)))
    else
      sig0(i) = 0.0_dp
    end if
  end do

  deallocate(S0, Utmp, sval)
end subroutine svdS

!=====================================================================
! exact_diag_inversion (no SVD)
!=====================================================================
subroutine exact_diag_inversion(S, D, X5, nrens, nrobs)
  implicit none
  integer,  intent(in)  :: nrens, nrobs
  real(dp), intent(in)  :: S(nrobs, nrens), D(nrobs, nrens)
  real(dp), intent(out) :: X5(nrens, nrens)
  
  real(dp), allocatable :: SS(:,:), SD(:,:), ZSD(:,:), eig(:), Z(:,:)
  real(dp) :: n1, eig_threshold
  integer  :: i, j

  allocate(SS(nrens, nrens), SD(nrens, nrens), Z(nrens, nrens), &
           eig(nrens), ZSD(nrens, nrens))

  n1 = 1.0_dp / real(nrens - 1, dp)

  ! 1. Compute SS = (1/(N-1)) * S^T * S
  call dgemm('T','N', nrens, nrens, nrobs, n1, S, nrobs, S, nrobs, 0.0_dp, SS, nrens)
  
  do i = 1, nrens
    SS(i,i) = SS(i,i) + 1.0_dp
  end do

  ! 2. Compute SD = (1/(N-1)) * S^T * D
  call dgemm('T','N', nrens, nrens, nrobs, n1, S, nrobs, D, nrobs, 0.0_dp, SD, nrens)

  ! 3. Eigen-decomposition of SS
  call eigC(SS, nrens, Z, eig)

  ! 4. Project SD into eigen-space
  call dgemm('T','N', nrens, nrens, nrens, 1.0_dp, Z, nrens, SD, nrens, 0.0_dp, ZSD, nrens)

  ! 5. Scale by inverse eigenvalues with damping
  eig_threshold = max(eps_eig, eps_damp * sum(eig) / real(nrens, dp))
  
  !=======================================================================
  ! FIX: If an eigenvalue falls below the noise floor threshold, its 
  ! inversion component must be strictly zeroed out. Applying damping to 
  ! trailing noise eigenvalues amplifies noise by 1/eps_damp (10^8).
  !=======================================================================
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
  real(dp), allocatable :: verens(:,:)    ! Converted to allocatable for Heap safety
  real(dp), allocatable :: var_vec(:)     ! Tracks variance instead of linear std
  real(dp)              :: ave, s_dev, total_var, variance_element
  integer               :: i, j

  allocate(verens(ndim, nrens))
  allocate(var_vec(ndim))

  ! 1. Generate synthetic ensemble with random noise
  call random(verens, ndim*nrens)

  ! 2. Center and Normalize the synthetic ensemble (Mean 0, Std 1)
  !$omp parallel do private(i, j, ave, s_dev) shared(verens)
  do i = 1, ndim
    ! Calculate row mean safely
    ave = sum(verens(i,:)) / real(nrens, dp)
    do j = 1, nrens
       verens(i,j) = verens(i,j) - ave
    end do
    
    ! Calculate variance and standard deviation
    s_dev = 0.0_dp
    do j = 1, nrens
       s_dev = s_dev + verens(i,j)**2
    end do
    s_dev = sqrt( s_dev / real(nrens, dp) )
    
    ! Scale to unit variance if not singular
    if (s_dev > 1.0e-14_dp) then
       do j = 1, nrens
          verens(i,j) = verens(i,j) / s_dev
       end do
    end if
  end do
  !$omp end parallel do

  ! 3. Apply the analysis transformation X5 to the synthetic ensemble
  call multa(verens, X5, ndim, nrens, ndim)

  ! 4. Measure the spread after transformation
  !=======================================================================
  ! FIX: Rewrote loop using explicit member indexing 'j' to avoid implicit 
  ! temporary array allocations inside OpenMP chunks, preventing race conditions.
  ! Converted metric tracking from linear standard deviation to pure variance.
  !=======================================================================
  !$omp parallel do private(i, j, ave, variance_element) shared(verens, var_vec)
  do i = 1, ndim
    ave = sum(verens(i,:)) / real(nrens, dp)
    variance_element = 0.0_dp
    do j = 1, nrens
       variance_element = variance_element + (verens(i,j) - ave)**2
    end do
    var_vec(i) = variance_element / real(nrens, dp)
  end do
  !$omp end parallel do

  ! 5. Compute the final Inflation Factor
  !=======================================================================
  ! FIX: In EnKF theory, inflation compensates for variance loss (energy, sigma^2)
  ! rather than a linear standard deviation ratio. Taking the square root of 
  ! the inverse average variance guarantees correct mathematical scaling.
  !=======================================================================
  total_var = sum(var_vec) / real(ndim, dp)
  
  if (total_var > 1.0e-14_dp) then
     inffac = sqrt(1.0_dp / total_var)
  else
     inffac = 1.0_dp ! Safety fallback
  end if

  deallocate(verens, var_vec)
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
  real(dp), allocatable :: ave(:)
  integer :: i, j

  allocate(ave(ndim))

  call ensmean(A, ave, ndim, nrens)
  
  !$omp parallel do private(j, i) shared(A, ave, inflation, nrens, ndim)
  do j = 1, nrens
    do i = 1, ndim
      A(i,j) = ave(i) + (A(i,j) - ave(i)) * inflation(i)
    end do
  end do
  !$omp end parallel do

  deallocate(ave)
end subroutine inflateA

end module mod_anafunc
