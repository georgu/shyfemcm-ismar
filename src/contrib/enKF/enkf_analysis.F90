!------------------------------------------------------------------------------
! Runs the Ensemble Kalman filter or the Square Root analysis by using the 
! Evenseen's routines
!------------------------------------------------------------------------------
! Parameters for the routine analysis.F90
! Depending on the type of analysis, different input
! matrices are necessary for the observations. However 
! the program computes all the matrices for both the
! cases
!------------------------------------------------------------------------------
  program enKF_analysis

  use mod_states
  use mod_enkf
  implicit none

  integer :: rmode = 13 ! Ensemble Kalman Filter
  !integer :: rmode = 23 ! Square root algorithm
  real :: truncation = 0.995 ! truncation of the SVD eigenvalues
  logical :: update_randrot = .true. ! False for local analysis

  logical :: verbose = .true. ! Prints diagnostic output

!----------------------------------------------------
! Opens info file
!----------------------------------------------------
  call read_info
  
!----------------------------------------------------
! Read basin
!----------------------------------------------------
  call read_basin

!----------------------------------------------------
! Init restart variables
!----------------------------------------------------
  call init_rst_vars

!----------------------------------------------------
! Read the ensemble
!----------------------------------------------------
  call read_ensemble

!----------------------------------------------------
! Makes the mean of A ans saves a restart file with it
! before the analysis
!----------------------------------------------------
  call average_mat(-1)

!----------------------------------------------------
! Read observations and makes D1, E and R, S and d-H*mean(A)
!----------------------------------------------------
  call make_matrices

!--------------------------------
! Call the analysis routine
!--------------------------------
! Note that the ENKF scheme needs D1, the observation innovations, not the
! perturbed measurements. See analysis2.F90, not analysis.F90, for a correct
! description.
  ! Decide if use an augmented state with the model error
  if( is_mod_err.eq.0 ) then

    call analysis(A,R,E,S,D1,innov,global_ndim,nrens,nobs_tot,verbose,truncation,rmode,update_randrot)
    !call analysis6c(A,E,S,innov,global_ndim,nrens,nobs_tot,verbose)	!SQRT alg
    !call analysis2(A,D1,R,S,global_ndim,nrens,nobs_tot,verbose)	!ENKF alg

  else if( is_mod_err.eq.1 ) then

    call push_aug
    call analysis(Aaug,R,E,S,D1,innov,2*global_ndim,nrens,nobs_tot,verbose,truncation,rmode,update_randrot)
    !call analysis6c(Aaug,E,S,innov,2*global_ndim,nrens,nobs_tot,verbose)	!SQRT alg
    !call analysis2(Aaug,D1,R,S,2*global_ndim,nrens,nobs_tot,verbose)		!ENKF alg
    call pull_aug

  else

    error stop 'enkf_analysis: invalid option for is_mod_err'

  end if

  write(*,*) 'Analysis done'

!--------------------------------
! Save the output in different restart files
!--------------------------------
  write(*,*) 'Writing the restart files with the analysis states...'
  call write_ensemble

!--------------------------------
! Save the average state
!--------------------------------
  call average_mat(-2)

  end program enKF_analysis
