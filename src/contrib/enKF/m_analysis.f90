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
      den_threshold = max(eps_inv, eps_damping * R(1,1))
      
      if (den > den_threshold) then
          ! IMPROVED: Apply light Tikhonov damping if condition number is poor
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
         !=======================================================================
         ! FIX: Allocate a dedicated temporary covariance matrix 'C' to prevent
         ! permanent destruction/overwriting of the original 'R' intent(inout) argument.
         ! This preserves observational error variances for adjacent local nodes or loops.
         !=======================================================================
         nrmin = nrobs
         allocate(Z(nrobs,nrobs), eig(nrobs))
         
         block
            real(dp), allocatable :: C(:,:)
            allocate(C(nrobs,nrobs))
            
            ! Initialize the innovation covariance matrix with the observation error term
            if (mode == 11) then
               ! Classical Perturbed EnKF: (N-1)*R
               C = real(nrens-1, dp) * R
            else
               ! Deterministic / Square-Root filter modes (e.g., Mode 21)
               C = real(nrens-1, dp) * R
            end if
            
            ! Compute C = S*S^T + C via BLAS dgemm matrix multiplication
            call dgemm('N','T', nrobs, nrobs, nrens, 1.0_dp, S, nrobs, S, nrobs, &
                       1.0_dp, C, nrobs)
            
            ! Perform spectral decomposition on the innovation covariance matrix C (instead of R)
            call eigC(C, nrobs, Z, eig)
         end block

         ! Apply the eigenvalue truncation threshold filter to handle low-rank approximations
         call eigsign(eig, nrobs, truncation)

         !=======================================================================
         ! FIX: Explicitly compute the reciprocal of valid eigenvalues (Pseudo-Inverse).
         ! Omitting this step causes the EnKF update step to multiply by the covariance
         ! instead of dividing, leading to unphysical explosive increment magnitudes.
         !=======================================================================
         do i = 1, nrmin
            if (eig(i) > eps_inv) then
                ! If the eigenvalue approaches the noise floor, apply adaptive Tikhonov regularization
                if (eig(i) < eps_damping) then
                    eig(i) = 1.0_dp / sqrt(eig(i)**2 + eps_damping*eig(i))
                else
                    eig(i) = 1.0_dp / eig(i)
                end if
            else
                eig(i) = 0.0_dp
            end if
         end do

      case (12, 22)
         ! IMPROVED: Low-rank case with adaptive rank selection
         nrmin = min(nrobs, nrens)
         allocate(Z(nrobs,nrmin), eig(nrmin))
         call lowrankCinv(S, R, nrobs, nrens, nrmin, Z, eig, truncation, verbose)

      case (13, 23)
         ! Low-rank with perturbed obs: SS^T + E*E^T
         nrmin = min(nrobs, nrens)
         allocate(Z(nrobs,nrmin), eig(nrmin))
         call lowrankE(S, E, nrobs, nrens, nrmin, Z, eig, truncation, verbose)

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

      !=======================================================================
      ! FIX/CRITICAL CRITERIA: Stochastic perturbed EnKF modes (11, 12, 13) 
      ! mathematically require an additive update formulation (A = A + Reps * X3).
      ! The multiplicative X5 formulation (A = A * X5) is strictly designed 
      ! for deterministic square-root filters (21, 22, 23) to avoid adding 
      ! explicit noise. Forcing X5 on stochastic modes destroys the ensemble 
      ! perturbations, collapsing the eigenvalues. We force lreps = .true. here.
      !=======================================================================
      lreps = .true.
      if (verbose) print '(a)', 'analysis: Representer approach forced for stochastic mode'
      
      allocate(Reps(ndim,nrobs))
      
      ! Compute representer matrix via BLAS dgemm: Reps = A * S^T
      ! This maps the ensemble anomalies into the model state space
      call dgemm('N','T', ndim, nrobs, nrens, 1.0_dp, A, ndim, S, nrobs, 0.0_dp, Reps, ndim)
      
      ! Apply damping to representer matrix to avoid unphysical local growth 
      ! near boundaries or dense observation clusters
      call damp_representer(Reps, ndim, nrobs, verbose)

      ! If adaptive inflation requires X5 diagnostics later in the code, 
      ! we compute X5 as a secondary tracking variable without using it for the update
      if (inflate == 2) then
          call dgemm('T','N', nrens, nrens, nrobs, 1.0_dp, S, nrobs, X3, nrobs, 0.0_dp, X5, nrens)
          do i = 1, nrens
              X5(i,i) = X5(i,i) + 1.0_dp
          end do
      end if

  case (21,22,23)
      ! Deterministic square-root modes correctly use the multiplicative X5 weighting matrix
      lreps = .false.
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
     !=======================================================================
     ! Correct additive update for stochastic EnKF (Modes 11, 12, 13)
     ! A_new = A_old + 1.0 * Reps * X3
     !=======================================================================
     call dgemm('N','N', ndim, nrens, nrobs, 1.0_dp, Reps, ndim, X3, nrobs, 1.0_dp, A, ndim)
     if (.not. islocal) call dumpX3(X3, S, nrobs, nrens)   
  else
     !=======================================================================
     ! Correct multiplicative update for deterministic Square-Root EnKF (Modes 21, 22, 23)
     ! A_new = A_old * X5
     !=======================================================================
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
  
  ! SANITY CHECK: Warn if analysis updates are larger than background spread
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

  ! Bounded inflation to prevent runaway variance growth
  ! Clamp inflation factor to reasonable stable range [0.95, 1.20]
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

end module m_analysis

