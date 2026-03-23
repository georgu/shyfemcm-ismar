!======================================================================
!  Module: mod_para
!
!  Purpose:
!    Central configuration parameters for:
!      - EnKF analysis modes
!      - Ensemble inflation
!      - Localisation
!      - Initial ensemble generation
!      - Observation handling
!
!  Notes:
!    * Time-related tolerances and decay scales use DOUBLE PRECISION.
!    * Observation QC bounds and flags remain REAL to match obs arrays.
!======================================================================
module mod_para
  use iso_fortran_env, only : dp => real64
  implicit none

  !--------------------------------------------------------------------
  ! Analysis method selector
  !--------------------------------------------------------------------
  integer, save :: rmode = 13
  !   13 : EnKF with SVD pseudo-inversion of SS' + EE'
  !   22 : Square-root method with SVD pseudo-inversion of SS' + (N-1)R
  !   23 : Square-root method with SVD pseudo-inversion of SS' + EE'
  !   10 : Exact update for diagonal observation-error covariance

  ! Diagnostics
  logical, parameter :: verbose = .true.

  !--------------------------------------------------------------------
  ! Ensemble inflation
  !--------------------------------------------------------------------
  integer, parameter :: inflate = 2    ! 0=off, 1=multiplicative, 2=adaptive
  real(dp),    parameter :: infmult = 1.0  ! Inflation multiplier when inflate=1

  !--------------------------------------------------------------------
  ! Local analysis (localisation)
  !--------------------------------------------------------------------
  integer, save :: is_local = 0        ! 0=off, 1=enabled

  !--------------------------------------------------------------------
  ! Innovation limiter
  !--------------------------------------------------------------------
  integer, parameter :: mode_an   = 0
  real(dp),    parameter :: inn_alpha = 5.0

  !--------------------------------------------------------------------
  ! Initial ensemble generation
  !--------------------------------------------------------------------
  integer, parameter :: fmult_init   = 10   ! supersampling factor
  real(dp),    parameter :: theta_init   = 0.0  ! rotation angle (degrees)
  real(dp),    parameter :: sigma_init_z = 0.03 ! std of free surface
  real(dp),    parameter :: sigma_init_t = 1.0  ! std of temperature
  real(dp),    parameter :: sigma_init_s = 1.0  ! std of salinity
  logical, parameter :: sample_fix_init = .true.  ! fixed sampling seed

  !--------------------------------------------------------------------
  ! Observation perturbations (temporal correlation)
  !--------------------------------------------------------------------
  ! Negative value means white noise (no temporal correlation)
  real(dp), parameter :: TTAU_0D = -1.0d0
  real(dp), parameter :: TTAU_2D = 3.0d0 * 3600.0d0

  !--------------------------------------------------------------------
  ! Observation quality-control and processing
  !--------------------------------------------------------------------
  real(dp), parameter :: KSTD = 2.0       ! std limiter factor (<=0 disables)

  ! Makes superobservations, depending on the local analysis radius
  logical, parameter  :: SUPEROBS = .true.

  ! Maximum biases to keep
  logical, parameter  :: NOBIAS = .false.
  real(dp), parameter :: ZBIAS_MAX = 0.2
  real(dp), parameter :: TBIAS_MAX = 2.
  real(dp), parameter :: SBIAS_MAX = 4.

  ! Check on observations
  logical, parameter  :: OBSCHK = .true.

  ! SVD truncation and perturbation options
  real(dp),    parameter :: truncation     = 0.99
  logical, parameter :: lrandrot       = .false.
  logical, parameter :: lupdate_randrot = .true.
  logical, parameter :: lsymsqrt       = .true.

  ! Time window (seconds) to accept obs near analysis time
  real(dp), parameter :: TEPS = 300.0d0

  ! Missing-observation flag
  real(dp), parameter :: OFLAG = -999.0

  ! Min/max QC bounds for observations
  real(dp), parameter :: TEM_MIN = -20.0, TEM_MAX = 60.0
  real(dp), parameter :: SAL_MIN =  -0.5, SAL_MAX = 60.0
  real(dp), parameter :: SSH_MIN =  -8.0, SSH_MAX =  8.0
  real(dp), parameter :: VEL_MIN = -90000.0, VEL_MAX = 90000.0

end module mod_para
