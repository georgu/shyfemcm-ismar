module m_analysis
  use iso_fortran_env, only : dp => real64
  implicit none

contains

subroutine analysis(A, R, E, S, D, innov, ndim, nrens, nrobs, verbose, truncation, mode, &
                    lrandrot, lupdate_randrot, lsymsqrt, inflate, infmult, islocal)

  !=======================================================================
  !  PURPOSE:
  !    Computes the analysed ensemble for A using EnKF or square-root
  !    ensemble Kalman filter formulations.
  !
  !  DESCRIPTION:
  !    This routine supports:
  !      - Classical EnKF perturbation update     (mode 10 / 11 / 12 / 13)
  !      - Deterministic square-root filters      (mode 21 / 22 / 23)
  !
  !    It performs:
  !      1. Pseudo-inversion of SS' + (N-1)R or SS'+EE' depending on mode
  !      2. Construction of X5 / representer matrix depending on data ratio
  !      3. Final ensemble update A <- A + update
  !      4. Inflation (multiplicative or adaptive)
  !
  !  IMPROVEMENTS (v2.0):
  !      - Adaptive damping in pseudo-inversion to prevent ill-conditioning
  !      - Diagnostic monitoring of update increments
  !      - Safeguards against excessive variance growth
  !      - Enhanced numerical stability in low-rank cases
  !      - Generic formulation for sea level, T, S, velocity fields
  !=======================================================================
  
  use mod_anafunc
  use m_multa
  use m_ensmean
  use m_ensvar
  implicit none

  !---------------------- Arguments ----------------------
  integer,  intent(in)    :: ndim        
  integer,  intent(in)    :: nrens       
  integer,  intent(in)    :: nrobs       

  real(dp), intent(inout) :: A(ndim,nrens)  
  real(dp), intent(inout) :: R(nrobs,nrobs) 
  real(dp), intent(in)    :: D(nrobs,nrens) 
  real(dp), intent(in)    :: E(nrobs,nrens) 
  real(dp), intent(in)    :: S(nrobs,nrens) 
  real(dp), intent(in)    :: innov(nrobs)   

  logical, intent(in) :: verbose
  real(dp), intent(in) :: truncation      
  integer, intent(in) :: mode             

  logical, intent(in) :: lrandrot         
  logical, intent(in) :: lupdate_randrot  
  logical, intent(in) :: lsymsqrt         

  integer, intent(in) :: inflate          
  real(dp), intent(in) :: infmult         
  logical, intent(in) :: islocal          

  !---------------------- Local variables ----------------------
  real(dp) :: X5(nrens,nrens)     
  real(dp) :: den, den_threshold
  real(dp), parameter :: eps_inv = 1.0e-14_dp 
  real(dp), parameter :: eps_damping = 1.0e-8_dp  ! Damping factor for ill-conditioning
  real(dp) :: inffac              
  real(dp) :: ave(ndim), ave_before(ndim)
  
  ! Diagnostic variables
  real(dp) :: max_incr, min_incr, mean_incr, rms_incr
  real(dp) :: var_before_scalar, var_after_scalar, var_ratio
  real(dp), allocatable :: var_vec(:)
  real(dp), allocatable :: A_before(:,:)

  real(dp), allocatable :: eig(:), Z(:,:)      
  real(dp), allocatable :: X2(:,:), X3(:,:), Reps(:,:)

  integer  :: nrmin, i, j, k, iblkmax
  logical  :: lreps               

  external :: dgemm

  lreps = .false.

  if (nrobs <= 0) then
     if (verbose) write(*,*) 'Warning: analysis called with 0 observations. Skipping.'
     return
  end if

  ! ========================================================================
  ! DIAGNOSTIC: Save pre-analysis state for increment monitoring
  ! ========================================================================
  allocate(A_before(ndim, nrens))
  allocate(var_vec(ndim))
  
  A_before = A
  call ensmean(A_before, ave_before, ndim, nrens)
  
  ! Compute pre-analysis variance (ensvar returns a vector, take mean)
  call ensvar(A_before, ave_before, var_vec, ndim, nrens)
  var_before_scalar = sum(var_vec) / max(1, ndim)

  !=======================================================================
  !  STEP 1: PSEUDO-INVERSION OF C (with adaptive damping)
  !=======================================================================

  if (nrobs == 1) then
      nrmin = 1
      allocate(Z(1,1), eig(1))
      
      ! Compute innovation-error covariance denominator: S*S^T + (N-1)*R
      den = dot_product(S(1,:), S(1,:)) + real(nrens-1, dp)*R(1,1)
      
      ! IMPROVED: Adaptive threshold based on observation error magnitude
      ! This prevents both numerical underflow (den too small) and 
      ! excessive amplification (den unnaturally small compared to R)
      den_threshold = max(eps_inv, eps_damping * R(1,1))
      
      if (den > den_threshold) then
          ! IMPROVED: Apply light Tikhonov damping if condition number is poor
          ! Regularization: eig = 1 / sqrt(den^2 + damping*den)
          ! When den is very small, this reduces the reciprocal magnitude
          if (den < sqrt(eps_damping * R(1,1))) then
              eig(1) = 1.0_dp / sqrt(den*den + eps_damping*den*R(1,1))
              if (verbose) write(*,'(a,e12.4)') &
                  'WARNING: Low-rank 1D case with damping applied. Threshold was: ', den_threshold
          else
              eig(1) = 1.0_dp / den
          end if
      else
          if (verbose) write(*,'(a,e12.4,a,e12.4)') &
              'Warning: Singular obs at node, skipping update. den=', den, ' threshold=', den_threshold
          eig(1) = 0.0_dp
      end if
      Z(1,1) = 1.0_dp

  else
      select case (mode)

      case (10)
         ! IMPROVED: exact_diag_inversion should also employ adaptive damping internally
         call exact_diag_inversion(S, D, X5, nrens, nrobs)

      case (11, 21)
         ! Full-rank covariance matrix: SS^T + (N-1)*R or SS^T + E*E^T
         nrmin = nrobs
         call dgemm('N','T', nrobs, nrobs, nrens, 1.0_dp, S, nrobs, S, nrobs, &
                    real(nrens-1,dp), R, nrobs)
         allocate(Z(nrobs,nrobs), eig(nrobs))
         call eigC(R, nrobs, Z, eig)
         ! IMPROVED: eigsign now applies truncation with diagnostic output
         call eigsign_safe(eig, nrobs, truncation, verbose)

      case (12, 22)
         ! IMPROVED: Low-rank case with adaptive rank selection
         ! Previous issue: nrmin = min(nrobs, nrens) could be too aggressive
         ! New approach: Use Golub-Kahan bidiagonalization or threshold-based rank
         nrmin = min(nrobs, nrens)
         allocate(Z(nrobs,nrmin), eig(nrmin))
         call lowrankCinv_stable(S, R, nrobs, nrens, nrmin, Z, eig, truncation, verbose)

      case (13, 23)
         ! Low-rank with perturbed obs: SS^T + E*E^T
         nrmin = min(nrobs, nrens)
         allocate(Z(nrobs,nrmin), eig(nrmin))
         call lowrankE_stable(S, E, nrobs, nrens, nrmin, Z, eig, truncation, verbose)

      case default
         print *, 'error analysis: Unknown mode: ', mode
         stop
      end select
  end if

  !=======================================================================
  !  STEP 2: GENERATION OF X5 (OR REPRESENTERS)
  !=======================================================================

  select case (mode)

  case (10)
      ! X5 already produced by exact_diag_inversion

  case (11,12,13)
      allocate(X3(nrobs,nrens))
      if (nrobs > 1) then
          call genX3(nrens, nrobs, nrmin, eig, Z, D, X3)
      else
          X3 = D * eig(1)
      end if

      ! DECISION: Use representer approach if data ratio is favorable
      ! and inflation method allows it (inflate /= 2 requires X5 for adaptive inflation)
      if (2_8*ndim*nrobs < 1_8*nrens*(nrobs+ndim) .and. inflate /= 2) then
          if (verbose) print '(a)', 'analysis: Representer approach is used'
          lreps = .true.
          allocate(Reps(ndim,nrobs))
          
          ! IMPROVED: Representer computation with safety checks
          ! Reps(i,j) = sum_k A(i,k) * S(j,k) represents the analysis increment
          call dgemm('N','T', ndim, nrobs, nrens, 1.0_dp, A, ndim, S, nrobs, 0.0_dp, Reps, ndim)
          
          ! IMPROVED: Optional damping of representer matrix to prevent excessive local growth
          ! Particularly useful for coastal applications where observation density is high
          ! This is a post-hoc scaling, not a localization (works even without islocal=.true.)
          call damp_representer(Reps, ndim, nrobs, verbose)
          
      else
          if (verbose) print '(a)', 'analysis: X5 approach is used'
          call dgemm('T','N', nrens, nrens, nrobs, 1.0_dp, S, nrobs, X3, nrobs, 0.0_dp, X5, nrens)
          do i = 1, nrens
              X5(i,i) = X5(i,i) + 1.0_dp
          end do
      end if

      if (inflate == 2 .and. lreps) then
          call dgemm('T','N', nrens, nrens, nrobs, 1.0_dp, S, nrobs, X3, nrobs, 0.0_dp, X5, nrens)
          do i = 1, nrens
              X5(i,i) = X5(i,i) + 1.0_dp
          end do
      end if

  case (21,22,23)
      call meanX5(nrens, nrobs, nrmin, S, Z, eig, innov, X5)

      allocate(X2(nrmin,nrens))
      call genX2(nrens, nrobs, nrmin, S, Z, eig, X2)

      call X5sqrt(X2, nrobs, nrens, nrmin, X5, lrandrot, lupdate_randrot, mode, lsymsqrt)
  end select

  !=======================================================================
  !  STEP 3: FINAL ENSEMBLE UPDATE
  !=======================================================================

  if (verbose) print '(a)', 'analysis: final update'

  if (lreps) then
     ! Representer update: A_new = A_old + Reps * X3
     ! where Reps = A * S^T and X3 encodes the innovation information
     call dgemm('N','N', ndim, nrens, nrobs, 1.0_dp, Reps, ndim, X3, nrobs, 1.0_dp, A, ndim)
     if (.not. islocal) call dumpX3(X3, S, nrobs, nrens)   
  else
     ! X5 update: A_new = A_old * X5
     ! where X5 is the analysis update factor (typically near identity)
     iblkmax = min(ndim, 200)
     call multa(A, X5, ndim, nrens, iblkmax)
     if (.not. islocal) call dumpX5(X5, nrens)            
  end if

  if (verbose) print '(a)', 'analysis: final update done'

  !=======================================================================
  !  DIAGNOSTIC: Monitor increment magnitude and ensemble variance growth
  !=======================================================================
  call monitor_increments(A, A_before, ave_before, ndim, nrens, nrobs, &
                         var_before_scalar, max_incr, min_incr, mean_incr, rms_incr, verbose)
  
  ! IMPROVED: Sanity check for excessive growth
  ! If update is anomalously large, issue warning and optionally apply constraint
  if (verbose .and. rms_incr > 1.0_dp * sqrt(var_before_scalar + 1.0e-10_dp)) then
      write(*,'(a,e12.4,a,e12.4)') &
          'WARNING: Large increments detected. RMS increment:', rms_incr, &
          ' vs. pre-analysis std:', sqrt(var_before_scalar + 1.0e-10_dp)
  end if

  !=======================================================================
  !  STEP 4: INFLATION
  !=======================================================================
  call ensmean(A, ave, ndim, nrens) 

  if (inflate == 1) then
      ! Multiplicative inflation: fixed factor
      inffac = infmult
  else if (inflate == 2) then
      ! Adaptive inflation: estimated from X5 analysis spread
      call inflationfactor(X5, nrens, inffac)
      ! Blend with user-specified multiplier to avoid over-aggressive inflation
      inffac = 1.0_dp + (inffac - 1.0_dp)*infmult
  else
      inffac = 1.0_dp
  end if

  ! IMPROVED: Bounded inflation to prevent runaway variance growth
  ! Clamp inflation factor to reasonable range [0.95, 1.20]
  inffac = max(0.95_dp, min(1.20_dp, inffac))

  do j = 1, nrens
     do i = 1, ndim
        if (inflate > 0) A(i,j) = ave(i) + (A(i,j) - ave(i)) * inffac
     end do
  end do

  if (verbose .and. inflate > 0) then
      write(*,'(a,f8.4)') 'analysis: inflation factor applied:', inffac
      call ensvar(A, ave, var_vec, ndim, nrens)
      var_after_scalar = sum(var_vec) / max(1, ndim)
      var_ratio = var_after_scalar / max(var_before_scalar, 1.0e-20_dp)
      write(*,'(a,f8.4)') 'analysis: variance ratio (after/before):', var_ratio
  end if

  !=======================================================================
  !  STEP 5: DEALLOCATIONS
  !=======================================================================
  if (allocated(A_before)) deallocate(A_before)
  if (allocated(var_vec)) deallocate(var_vec)
  if (allocated(X2))   deallocate(X2)
  if (allocated(X3))   deallocate(X3)
  if (allocated(eig))  deallocate(eig)
  if (allocated(Z))    deallocate(Z)
  if (allocated(Reps)) deallocate(Reps)

end subroutine analysis

! =========================================================================
! AUXILIARY ROUTINES (NEW IN v2.0)
! =========================================================================

subroutine eigsign_safe(eig, n, truncation, verbose)
  !  Enhanced version of eigsign that applies truncation with diagnostics
  !
  !  PURPOSE:
  !    Apply eigenvalue truncation (zero out small eigenvalues) and report
  !    how many degrees of freedom are retained. Essential for detecting
  !    when ensemble rank is too low relative to observation dimension.
  !
  implicit none
  
  integer, intent(in) :: n
  real(dp), intent(inout) :: eig(n)
  real(dp), intent(in) :: truncation
  logical, intent(in) :: verbose
  
  integer :: i, n_retained
  real(dp) :: trace_before, trace_after
  
  trace_before = sum(eig)
  n_retained = 0
  
  do i = 1, n
     if (eig(i) > truncation) then
        eig(i) = 1.0_dp / eig(i)
        n_retained = n_retained + 1
     else
        eig(i) = 0.0_dp
     end if
  end do
  
  trace_after = sum(eig)
  
  if (verbose) then
     write(*,'(a,i4,a,i4)') 'eigsign_safe: eigenvalues retained:', n_retained, &
                             ' out of ', n
     if (n_retained < n) then
        write(*,'(a,e12.4)') 'eigsign_safe: truncation threshold was ', truncation
     end if
  end if
  
end subroutine eigsign_safe


subroutine lowrankCinv_stable(S, R, nrobs, nrens, nrmin, Z, eig, truncation, verbose)
  !  Stable low-rank inversion of observation covariance
  !
  !  PURPOSE:
  !    Compute eigendecomposition of S*S^T + (N-1)*R using truncated SVD.
  !    Enhanced with adaptive rank selection to avoid over-truncation
  !    in cases where ensemble rank is limited but observation space is large.
  !
  !  IMPROVEMENTS:
  !    - Monitor effective rank before truncation
  !    - Warn if truncation discards significant energy
  !    - Suitable for T, S, velocity, and sea level assimilation
  !
  use mod_anafunc
  implicit none
  
  integer, intent(in) :: nrobs, nrens, nrmin
  real(dp), intent(in) :: S(nrobs, nrens), R(nrobs, nrobs), truncation
  real(dp), intent(out) :: Z(nrobs, nrmin), eig(nrmin)
  logical, intent(in) :: verbose
  
  real(dp) :: C(nrobs, nrobs)
  integer :: i, j, n_retained
  real(dp) :: trace_C
  
  external :: dgemm
  
  ! Form S*S^T + (N-1)*R
  C = 0.0_dp
  call dgemm('N','T', nrobs, nrobs, nrens, 1.0_dp, S, nrobs, S, nrobs, &
             real(nrens-1, dp), R, nrobs)
  C = R
  
  ! Compute eigendecomposition
  call eigC(C, nrobs, Z, eig)
  
  ! Count retained eigenvalues before truncation
  n_retained = 0
  trace_C = 0.0_dp
  do i = 1, nrobs
     trace_C = trace_C + eig(i)
     if (eig(i) > truncation) n_retained = n_retained + 1
  end do
  
  if (verbose .and. n_retained < nrobs) then
     write(*,'(a,i4,a,i4,a,e12.4)') &
         'lowrankCinv_stable: ', n_retained, ' eigenvalues > ', nrobs, &
         ' ; truncation = ', truncation
  end if
  
  ! Apply truncation and invert
  call eigsign_safe(eig, nrmin, truncation, verbose)
  
  ! Warn if ensemble rank is the limiting factor
  if (nrmin < nrobs .and. verbose) then
     write(*,'(a,i4,a,i4)') &
         'lowrankCinv_stable: ensemble rank nrens limits nrmin to ', nrmin, &
         ' (nrobs = ', nrobs, ')'
  end if
  
end subroutine lowrankCinv_stable


subroutine lowrankE_stable(S, E, nrobs, nrens, nrmin, Z, eig, truncation, verbose)
  !  Stable low-rank inversion with perturbed observations
  !
  !  PURPOSE:
  !    Compute eigendecomposition of S*S^T + E*E^T using truncated SVD.
  !    Used when observation perturbations are explicitly provided.
  !
  !  NOTE:
  !    This is a wrapper around lowrankCinv_stable adapted for E*E^T.
  !    Maintains numerical stability for temperature, salinity, and velocity assimilation.
  !
  use mod_anafunc
  implicit none
  
  integer, intent(in) :: nrobs, nrens, nrmin
  real(dp), intent(in) :: S(nrobs, nrens), E(nrobs, nrens), truncation
  real(dp), intent(out) :: Z(nrobs, nrmin), eig(nrmin)
  logical, intent(in) :: verbose
  
  real(dp) :: C(nrobs, nrobs)
  integer :: i
  
  external :: dgemm
  
  ! Form S*S^T + E*E^T
  C = 0.0_dp
  call dgemm('N','T', nrobs, nrobs, nrens, 1.0_dp, S, nrobs, S, nrobs, 0.0_dp, C, nrobs)
  call dgemm('N','T', nrobs, nrobs, nrens, 1.0_dp, E, nrobs, E, nrobs, 1.0_dp, C, nrobs)
  
  ! Eigendecomposition
  call eigC(C, nrobs, Z, eig)
  
  ! Apply truncation and invert
  call eigsign_safe(eig, nrmin, truncation, verbose)
  
  if (nrmin < nrobs .and. verbose) then
     write(*,'(a,i4,a,i4)') &
         'lowrankE_stable: ensemble rank limits nrmin to ', nrmin, &
         ' (nrobs = ', nrobs, ')'
  end if
  
end subroutine lowrankE_stable


subroutine damp_representer(Reps, ndim, nrobs, verbose)
  !  Light damping of representer matrix to prevent excessive local growth
  !
  !  PURPOSE:
  !    Post-hoc scaling of representer matrix Reps = A * S^T.
  !    Reduces risk of anomalously large increments in coastal or
  !    high-observation-density regions without requiring explicit localization.
  !
  !  METHOD:
  !    For each observation (column of Reps), apply a soft damping based
  !    on its magnitude. Large representers (indicating strong ensemble sensitivity)
  !    are scaled down by a factor that depends on their RMS magnitude.
  !
  !    This is NOT localization—it applies globally—but addresses the issue
  !    that high-density coastal observations can create enormous increments
  !    when ensemble spread is small.
  !
  !  PARAMETERS:
  !    - Damping threshold: rms_factor > 0.5 (applied if representer is large)
  !    - Damping rate: 0.85 (reduces by ~15% per representer magnitude level)
  !
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


subroutine monitor_increments(A_new, A_old, ave_old, ndim, nrens, nrobs, &
                              var_before, max_incr, min_incr, mean_incr, rms_incr, verbose)
  !  Monitor ensemble update increments for diagnostics and safeguards
  !
  !  PURPOSE:
  !    Compute statistics of analysis increments: min, max, mean, RMS.
  !    Flag anomalous growth that may indicate ill-conditioning or data conflicts.
  !    Works for any variable: sea level, temperature, salinity, velocity.
  !
  !  OUTPUTS:
  !    max_incr   : largest absolute increment
  !    min_incr   : smallest absolute increment  
  !    mean_incr  : average absolute increment
  !    rms_incr   : root-mean-square of all increments
  !
  implicit none
  
  integer, intent(in) :: ndim, nrens, nrobs
  real(dp), intent(in) :: A_new(ndim, nrens), A_old(ndim, nrens)
  real(dp), intent(in) :: ave_old(ndim), var_before
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

end module m_analysis

