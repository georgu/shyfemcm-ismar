!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
!
module mod_para

  integer :: rmode = 13 ! Ensemble Kalman Filter with SVD pseudo inversion of SS'+ EE'
  !integer :: rmode = 22 ! Square root algorithm with SVD pseudo inversion of SS'+(N-1)R
  !integer :: rmode = 23 ! Square root algorithm with SVD pseudo inversion of SS'+ EE'
  !integer :: rmode = 21

  logical :: verbose = .true. ! Prints diagnostic output

  integer, parameter :: type_infl = 1  ! 1: RTPS inflation (see WHITAKER, 2012)
                                       ! 2: Multiplication inflation. Seems better with
				       !    uniform observation, but can explode if there 
                                       !    are grid areas without observations (spread 
                                       !    is not reduced)
  real, parameter :: alpha_infl = 0.02 ! type_infl = 1 -> ~ 0.01 
                                       ! type_infl = 2 -> ~ 0.01 (lower with SQRT method)

  ! Set these parameters for local analysis. Important.				       
  !
  integer, parameter :: is_local = 1 !Local analysis. 0 disable, 1 local analysis
  real, parameter :: rho_loc = 6.   !Radius for local analysis (use the same coords of the grid)

  ! set this to 1 to include the model errors in the analysis (to test) or equal 2 to
  ! include model parameters (todo)
  !
  integer, parameter :: mode_an = 0

!------------
! Settings for the observation perturbations (mod_enkf)
!
  ! decay time for the red noise of the observations (sec). 
  ! Set lower than 0 to disable (white noise)
  !
  double precision, parameter :: TTAU_0D = -1
  double precision, parameter :: TTAU_2D = 3*3600.

!------------
! Settings to manage the observations (mod_manage_obs)
!
  ! multiplication factor for the ens std, to resize
  ! the obs std (default = 2). Set <= 0 to disable it.
  ! See Sakov et al. 2012.
  !
  real, parameter :: KSTD = 2

  ! parameters for the analysis
  real :: truncation = 0.995 ! truncation of the SVD eigenvalues
  logical :: update_randrot = .true. 
  !logical :: update_randrot = .false. ! False to avoid randomness

  ! time interval (sec) to select an observation
  ! with respect to the analysis time
  !
  double precision, parameter :: TEPS = 300.

  ! standard flag for a missing observation
  !
  real, parameter :: OFLAG = -999.

  ! min-max values for the observation and model check
  !
  real, parameter :: TEM_MIN = -5.0d0
  real, parameter :: TEM_MAX = 50.0d0
  real, parameter :: SAL_MIN = 0.0d0
  real, parameter :: SAL_MAX = 40.0d0
  real, parameter :: SSH_MIN = -3.0d0
  real, parameter :: SSH_MAX = 3.0d0
  real, parameter :: VEL_MIN = 0.0d0
  real, parameter :: VEL_MAX = 4.0d0	!used also for u and v, min with -

end module mod_para
