subroutine analysis(A, R, E, S, D, innov, ndim, nrens, nrobs, verbose, truncation, mode, &
                    lrandrot, lupdate_randrot, lsymsqrt, inflate, infmult)

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
  !      3. Final ensemble update A ← A + update
  !      4. Optional inflation (multiplicative or adaptive)
  !
  !    All operations are performed in double precision.
  !=======================================================================

  use iso_fortran_env, only : dp => real64
  use mod_anafunc
  use m_multa
  use m_ensmean
  use m_ensvar
  implicit none

  !---------------------- Arguments ----------------------
  integer,  intent(in)    :: ndim        ! Dimension of model state vector
  integer,  intent(in)    :: nrens       ! Number of ensemble members
  integer,  intent(in)    :: nrobs       ! Number of observations

  real(dp), intent(inout) :: A(ndim,nrens)  ! Ensemble matrix (modified)
  real(dp), intent(inout) :: R(nrobs,nrobs) ! R is updated in some modes → must be INOUT
  real(dp), intent(in)    :: D(nrobs,nrens) ! Perturbed measurement innovations (d + E - H*A)
  real(dp), intent(in)    :: E(nrobs,nrens) ! Observation perturbations (used in mode 13 / 23)
  real(dp), intent(in)    :: S(nrobs,nrens) ! HA' (observation operator applied to ensemble)
  real(dp), intent(in)    :: innov(nrobs)   ! Innovation vector: d - H*mean(A)

  logical, intent(in) :: verbose
  real(dp), intent(in) :: truncation      ! Fraction of variance retained in pseudo-inversion
  integer, intent(in) :: mode             ! 10/11/12/13 (EnKF) or 21/22/23 (sqrt)

  logical, intent(in) :: lrandrot         ! Random rotation for square-root updates
  logical, intent(in) :: lupdate_randrot  ! True only when rotation is updated at this point
  logical, intent(in) :: lsymsqrt         ! Use symmetric square-root (Sakov)

  integer, intent(in) :: inflate          ! 0=Off, 1=Multiplicative, 2=Adaptive (Evensen 2009)
  real(dp), intent(in) :: infmult         ! Inflation multiplier / adjustment

  !---------------------- Local variables ----------------------
  real(dp) :: X5(nrens,nrens)     ! Transformation matrix for ensemble perturbations
  real(dp) :: inffac              ! Inflation factor
  real(dp) :: ave(ndim)           ! Ensemble mean

  ! Deferred-shape arrays must be allocatable to be allocated at runtime.
  real(dp), allocatable :: eig(:), Z(:,:)      ! Eigenvalues/eigenvectors or SVD components
  real(dp), allocatable :: X2(:,:), X3(:,:), Reps(:,:)

  integer  :: nrmin, i, j, iblkmax
  logical  :: lreps               ! True → representer approach is used for few observations

  external :: dgemm

  lreps = .false.

  !=======================================================================
  !  STEP 1: PSEUDO-INVERSION OF C = S*S' + (N-1)*R  OR  S*S' + E*E'
  !
  !  Modes:
  !    10  = exact diagonal case
  !    11  = eigenvalue pseudo-inversion   (C = SS' + (N-1)R)
  !    12  = low-rank SVD pseudo-inversion (C = SS' + (N-1)R)
  !    13  = low-rank SVD pseudo-inversion (C = SS' + EE')
  !    21/22/23 = same as 11/12/13, but for deterministic square-root updates
  !=======================================================================

  if (nrobs == 1) then
      ! C is scalar → invert directly
      nrmin = 1
      allocate(Z(1,1), eig(1))
      eig(1) = dot_product(S(1,:), S(1,:)) + real(nrens-1,dp)*R(1,1)
      eig(1) = 1.0_dp / eig(1)
      Z(1,1) = 1.0_dp

  else
      select case (mode)

      case (10)
         ! Exact diagonal inversion (special case filter)
         call exact_diag_inversion(S, D, X5, nrens, nrobs)

      case (11, 21)
         ! Full eigen-decomposition of C = SS' + (N-1)R
         nrmin = nrobs
         call dgemm('N','T', nrobs, nrobs, nrens, 1.0_dp, S, nrobs, S, nrobs, &
                    real(nrens-1,dp), R, nrobs)
         allocate(Z(nrobs,nrobs), eig(nrobs))
         call eigC(R, nrobs, Z, eig)
         call eigsign(eig, nrobs, truncation)

      case (12, 22)
         ! Low-rank pseudo-inversion using SVD of C = SS' + (N-1)R
         nrmin = min(nrobs, nrens)
         allocate(Z(nrobs,nrmin), eig(nrmin))
         call lowrankCinv(S, R, nrobs, nrens, nrmin, Z, eig, truncation)

      case (13, 23)
         ! SVD pseudo-inversion using EE' instead of R
         nrmin = min(nrobs, nrens)
         allocate(Z(nrobs,nrmin), eig(nrmin))
         call lowrankE(S, E, nrobs, nrens, nrmin, Z, eig, truncation)

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
      ! Build X3 from eig/Z/D
      allocate(X3(nrobs,nrens))
      if (nrobs > 1) then
          call genX3(nrens, nrobs, nrmin, eig, Z, D, X3)
      else
          X3 = D * eig(1)
      end if

      ! Choose representers vs X5 based on work estimate (avoid overflow with 64-bit literals)
      if (2_8*ndim*nrobs < 1_8*nrens*(nrobs+ndim) .and. inflate /= 2) then
          if (verbose) print '(a)', 'analysis: Representer approach is used'
          lreps = .true.

          allocate(Reps(ndim,nrobs))
          call dgemm('N','T', ndim, nrobs, nrens, 1.0_dp, A, ndim, S, nrobs, 0.0_dp, Reps, ndim)

      else
          if (verbose) print '(a)', 'analysis: X5 approach is used'
          call dgemm('T','N', nrens, nrens, nrobs, 1.0_dp, S, nrobs, X3, nrobs, 0.0_dp, X5, nrens)
          do i = 1, nrens
              X5(i,i) = X5(i,i) + 1.0_dp
          end do
      end if

  case (21,22,23)
      ! Mean update
      call meanX5(nrens, nrobs, nrmin, S, Z, eig, innov, X5)

      ! Generate perturbation matrix X2
      allocate(X2(nrmin,nrens))
      call genX2(nrens, nrobs, nrmin, S, Z, eig, X2)

      ! Construct full square-root transformation
      call X5sqrt(X2, nrobs, nrens, nrmin, X5, lrandrot, lupdate_randrot, mode, lsymsqrt)
  end select

  !=======================================================================
  !  STEP 3: FINAL ENSEMBLE UPDATE
  !=======================================================================

  if (verbose) print '(a)', 'analysis: final update'

  if (lreps) then
      ! Representer update: A += Reps * X3
      call dgemm('N','N', ndim, nrens, nrobs, 1.0_dp, Reps, ndim, X3, nrobs, 1.0_dp, A, ndim)
      if (ndim > 1000) call dumpX3(X3, S, nrobs, nrens)   ! no write for local analysis (omp) - mbj
  else
      ! Full X5 transform: A = A * X5
      iblkmax = min(ndim, 200)
      call multa(A, X5, ndim, nrens, iblkmax)
      if (ndim > 1000) call dumpX5(X5, nrens)            ! no write for local analysis (omp) - mbj
  end if

  if (verbose) print '(a)', 'analysis: final update done'

  !=======================================================================
  !  STEP 4: INFLATION (optional)
  !=======================================================================
  if (inflate == 1) then
      inffac = infmult
  elseif (inflate == 2) then
      call inflationfactor(X5, nrens, inffac)              ! adaptive inflation factor
      inffac = 1.0_dp + (inffac - 1.0_dp)*infmult          ! adjust adaptive factor
  end if

  if (inflate > 0) then
      if (verbose) print '(a,f10.4)', 'analysis: inflation update with inflation factor= ', inffac
      call ensmean(A, ave, ndim, nrens)
      do j = 1, nrens
          do i = 1, ndim
              A(i,j) = ave(i) + (A(i,j) - ave(i))*inffac
          end do
      end do
  end if

  !=======================================================================
  !  STEP 5: DEALLOCATIONS
  !=======================================================================
  if (allocated(X2))   deallocate(X2)
  if (allocated(X3))   deallocate(X3)
  if (allocated(eig))  deallocate(eig)
  if (allocated(Z))    deallocate(Z)
  if (allocated(Reps)) deallocate(Reps)

end subroutine analysis
