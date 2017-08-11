module mod_para

!------------
! Settings for the analysis step (enkf_analysis)
!
  integer :: rmode = 13 ! Ensemble Kalman Filter
  !integer :: rmode = 23 ! Square root algorithm

  logical :: verbose = .true. ! Prints diagnostic output

  ! do not touch these
  !
  real :: truncation = 0.995 ! truncation of the SVD eigenvalues
  logical :: update_randrot = .true. ! False for local analysis

!------------
! Settings for the observation perturbations (mod_enkf)
!
  ! decay time for the red noise of the observations (sec). 
  ! Set lower than 0 to disable (white noise)
  !
  double precision, parameter :: TTAU_0DLEV = -1
  double precision, parameter :: TTAU_2DVEL = 3*3600.

!------------
! Settings to manage the observations (mod_manage_obs)
!
  ! multiplication factor for the ens std, to resize
  ! the obs std (default = 2). Set <= 0 to disable it.
  ! See Sakov et al. 2012.
  !
  real, parameter :: KSTD = 2

  ! time interval (sec) to select an observation
  ! with respect to the analysis time
  !
  double precision, parameter :: TEPS = 300.

  ! standard flag for a missing observation
  !
  real, parameter :: OFLAG = -999.

  ! min-max values for the observation and model check
  !
  real, parameter :: TEM_MIN = 0.0d0
  real, parameter :: TEM_MAX = 50.0d0
  real, parameter :: SAL_MIN = 0.0d0
  real, parameter :: SAL_MAX = 40.0d0
  real, parameter :: SSH_MIN = -3.0d0
  real, parameter :: SSH_MAX = 3.0d0
  real, parameter :: VEL_MIN = 0.0d0
  real, parameter :: VEL_MAX = 4.0d0	!used also for u and v, min with -

end module mod_para
