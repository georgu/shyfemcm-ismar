! ======================================================================
!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights
! reserved.
!
! ======================================================================
!  Module: mod_obs_states
!
!  PURPOSE:
!    Define lightweight data structures for observations used by
!    the EnKF data assimilation system:
!      - files      : entry in the observation file list
!      - scalar_0d  : 0-D point observations (e.g., tide gauges)
!      - vector_2d  : 2-D gridded vector fields (e.g., currents)
!
!  STATUS CODES (kept as in your system):
!    0 = normal obs (assimilated)
!    1 = super-observation (assimilated)
!    2 = merged into a super-observation (not assimilated)
!    3 = out of range (not assimilated)
!    4 = flagged value (not assimilated)
!
!  PRECISION POLICY:
!    - Time fields (t) are DOUBLE PRECISION, matching downstream
!      time handling in analysis modules.
!    - Spatial coordinates (x, y, z) are REAL (single) to preserve
!      binary I/O compatibility with FEM readers.
!    - Observation VALUES (val, u, v) are REAL (single) for
!      consistency with model fields.
!    - Observation STD and rhol are REAL (single) for practical
!      observation error specification.
!
!  IMPROVEMENTS (v2.0):
!    - Explicit documentation of covariance matrix generation for mode 13
!    - Safety notes on observation error handling
!    - Clear distinction between model and observation precision
!    - Enhanced comments on matrix construction in analysis
!
!  IMPORTANT NOTES ON COVARIANCE MATRICES:
!    =====================================================
!    Mode 12 (SS^T + (N-1)*R):
!      - R should be a DIAGONAL matrix of observation errors
!      - Pass diag(std(:)^2) as the observation error covariance
!      - Standard Kalman formulation
!
!    Mode 13 (SS^T + E*E^T):
!      - E is the PERTURBED OBSERVATION matrix
!      - Each row of E represents one perturbed observation
!      - E has shape (nrobs, nrens), one perturbation per member
!      - Construction in mod_anafunc:
!        * If you have obs errors std(:), create:
!          E(i,j) = sqrt(0.5 * std(i)^2) * random_normal(0,1)
!        * Then E*E^T ≈ 0.5 * diag(std(:)^2) (averaged over ensemble)
!        * This adds STOCHASTIC perturbations, not fixed covariance
!      - CAVEAT: E*E^T is not deterministic like R!
!        * Different ensemble runs give different E*E^T
!        * For reproducibility, fix random seed BEFORE obs perturbations
!      - Better suited for large observation sets where full R is expensive
!
!    Mode 10 or 11 (Full-rank, explicit R):
!      - Use if you need full observation error covariance matrix
!      - Allows off-diagonal correlations
!      - More expensive computationally
!    =====================================================
!
!=======================================================================
module mod_obs_states
  use iso_fortran_env, only : dp => real64
  implicit none
  private

  ! Make the types available to users of this module
  public :: files, scalar_0d, vector_2d

  !-------------------------------------------------------------
  ! Observation file entry
  !   ty   : type tag (e.g., '0DLEV', '0DTEM', '0DSAL', '2DVEL')
  !   name : path or basename of the file
  !-------------------------------------------------------------
  type files
     character(len=5)  :: ty       ! Type of observations (0DLEV, 0DTEM, 0DSAL, 2DVEL, etc.)
     character(len=80) :: name     ! Name of the file (path or basename)
     real(dp)              :: x        ! X coordinate (for single-point files, unused for 2D)
     real(dp)              :: y        ! Y coordinate (for single-point files, unused for 2D)
     real(dp)              :: z        ! Z coordinate (depth/height, unused for surface observations)
     real(dp)              :: std      ! Observation standard deviation (assumed uniform)
     real(dp)              :: rhol     ! Localization radius for local analysis (in model units)

  end type files

  !-------------------------------------------------------------
  ! 0-D SCALAR OBSERVATION AT A SINGLE LOCATION
  !
  ! USAGE:
  !   Represents point measurements: tide gauges, thermometers,
  !   conductivity probes, etc.
  !
  ! PRECISION NOTES:
  !   - t      : DOUBLE because time arithmetic downstream
  !   - x, y, z: REAL (single) for FEM mesh compatibility
  !   - val    : REAL (single) to match model variables
  !   - std    : REAL (single) observation error estimate
  !   - rhol   : REAL (single) localization radius
  !
  ! COVARIANCE IMPLICATIONS:
  !   For a single observation, the obs error covariance R is a
  !   1x1 scalar: R(1,1) = std^2
  !
  !   In analysis (mod_anafunc):
  !   - Mode 12: R is read as a 1x1 matrix with R(1,1) = std^2
  !   - Mode 13: E is a 1 x nrens matrix of perturbed obs:
  !              E(1,j) = sqrt(std^2) * random_normal_j
  !              Then E*E^T ≈ std^2 (on average)
  !-------------------------------------------------------------
  type scalar_0d
     real(dp)          :: t       ! Time (absolute, for temporal collocation)
     real(dp)              :: x       ! X coordinate (longitude or projected)
     real(dp)              :: y       ! Y coordinate (latitude or projected)
     real(dp)              :: z       ! Z coordinate (depth: negative or positive, or NaN)
     real(dp)              :: val     ! Observed value
     real(dp)              :: std     ! Observation standard deviation
     integer           :: stat    ! Status (0=normal, 1=super-obs, 2=merged, 3=out-of-range, 4=flagged)
     integer           :: id      ! ID number of source file
     real(dp)              :: rhol    ! Localization radius (0 = no localization)

  end type scalar_0d

  !-------------------------------------------------------------
  ! 2-D VECTOR FIELD OBSERVATION ON A REGULAR GRID
  !
  ! USAGE:
  !   Represents gridded vector fields: satellite SST/SSH maps,
  !   gridded current measurements, etc.
  !
  ! PRECISION NOTES:
  !   - t        : DOUBLE because time arithmetic
  !   - x(:,:), y(:,:): REAL (single) grid coordinates
  !   - u(:,:), v(:,:): REAL (single) vector components
  !   - std(:,:) : REAL (single) per-point observation errors
  !   - z        : REAL (single) constant depth
  !
  ! COVARIANCE IMPLICATIONS:
  !   For 2D gridded observations, the full obs error covariance R
  !   would be nrobs x nrobs, where nrobs = nx * ny * ncomps.
  !
  !   DIAGONAL assumption (commonly used):
  !      R(i,i) = std(i)^2  for all i
  !      R(i,j) = 0         for i ≠ j
  !   This is efficient but ignores spatial correlations.
  !
  !   In analysis (mod_anafunc):
  !   - Mode 12: R constructed as diag(reshape(std(:,:)^2, [nrobs]))
  !              CAVEAT: This assumes INDEPENDENT spatial errors
  !              which may be unrealistic for satellite data.
  !
  !   - Mode 13: E is nrobs x nrens matrix:
  !              E(i,j) = std_reshaped(i) * random_normal(0,1)
  !              Each perturbed observation is independent
  !              CAVEAT: Correlations are NOT preserved across
  !              ensemble members (stochastic EnKF limitation).
  !
  !   BETTER PRACTICE for spatially correlated obs:
  !   - Use mode 11 with FULL R matrix
  !   - Construct R with off-diagonal correlations
  !   - E.g., R(i,j) = std(i)*std(j)*exp(-d_ij^2/L_c^2)
  !   - This is expensive but realistic
  !
  !-------------------------------------------------------------
  type vector_2d
     real(dp)          :: t              ! Time of field (absolute)
     integer           :: nx, ny         ! Grid dimensions (nx columns, ny rows)
     real(dp), allocatable :: x(:,:)         ! X coordinates (nx x ny)
     real(dp), allocatable :: y(:,:)         ! Y coordinates (nx x ny)
     real(dp)              :: z              ! Depth/height (constant for entire field)
     real(dp), allocatable :: u(:,:)         ! U component (nx x ny)
     real(dp), allocatable :: v(:,:)         ! V component (nx x ny)
     real(dp), allocatable :: std(:,:)       ! Observation error std (nx x ny)
     integer, allocatable :: stat(:,:)   ! Status flags per grid point (0..4)
     integer           :: id             ! ID number of source file

  end type vector_2d

contains

!=======================================================================
! UTILITY: Check consistency of scalar_0d observation
! (Optional: for debugging/validation)
!=======================================================================
subroutine check_scalar_obs(obs)
   implicit none
   type(scalar_0d), intent(in) :: obs

   if (obs%std <= 0.0) then
      write(*,'(a)') 'WARNING: check_scalar_obs: std <= 0'
      write(*,'(a,e12.4)') '  std = ', obs%std
   end if

   if (isnan(obs%val)) then
      write(*,'(a)') 'WARNING: check_scalar_obs: value is NaN'
   end if

   if (obs%rhol < 0.0) then
      write(*,'(a,e12.4)') 'WARNING: check_scalar_obs: negative rhol =', obs%rhol
   end if

end subroutine check_scalar_obs

!=======================================================================
! UTILITY: Check consistency of vector_2d observation
!=======================================================================
subroutine check_vector_obs(obs)
   implicit none
   type(vector_2d), intent(in) :: obs
   integer :: i, j, nnan, nbad

   if (.not. allocated(obs%u) .or. .not. allocated(obs%v)) then
      write(*,'(a)') 'WARNING: check_vector_obs: u or v not allocated'
      return
   end if

   nnan = 0
   nbad = 0

   do j = 1, obs%ny
      do i = 1, obs%nx
         if (isnan(obs%u(i,j)) .or. isnan(obs%v(i,j))) nnan = nnan + 1
         if (obs%std(i,j) <= 0.0) nbad = nbad + 1
      end do
   end do

   if (nnan > 0) write(*,'(a,i8,a,i8,a,i8)') &
       'WARNING: check_vector_obs: Found ', nnan, ' NaN values in ', &
       obs%nx*obs%ny, ' total grid points'

   if (nbad > 0) write(*,'(a,i8,a,i8,a,i8)') &
       'WARNING: check_vector_obs: Found ', nbad, ' bad std values in ', &
       obs%nx*obs%ny, ' total grid points'

end subroutine check_vector_obs

end module mod_obs_states

