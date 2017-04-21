!------------------------------------------------------------------------------
! Run the ensemble Kalman filter analysis using the Evenseen's routines
!------------------------------------------------------------------------------
  program enKF_analysis

  use mod_states
  use mod_enkf
  implicit none

  ! Analysis parameters for the routine analysis.F90
  integer :: rmode = 13 ! mode to run the program
  !integer :: rmode = 23 ! mode to run the program
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
! Read the ensemble
!----------------------------------------------------
  call read_ensemble

!----------------------------------------------------
! Makes the mean of A ans saves a restart file with it
! before the analysis
!----------------------------------------------------
  call average_mat(-1)

!----------------------------------------------------
! Read observations and makes D, E and R, S and d-H*mean(A)
!----------------------------------------------------
  call make_matrices

!--------------------------------
! Call the analysis routine
!--------------------------------
  ! Decide if use an augmented state with the model error
  if( is_mod_err.eq.0 ) then

    call analysis(A,R,E,S,D,innov,global_ndim,nrens,nobs_tot,verbose,truncation,rmode,update_randrot)
    !call analysis6c(A,E,S,innov,global_ndim,nrens,nobs_tot,verbose)
    !call analysis2(A,D,R,S,global_ndim,nrens,nobs_tot,verbose)

  else if( is_mod_err.eq.1 ) then

    call push_aug
    call analysis(Aaug,R,E,S,D,innov,2*global_ndim,nrens,nobs_tot,verbose,truncation,rmode,update_randrot)
    !call analysis6c(Aaug,E,S,innov,2*global_ndim,nrens,nobs_tot,verbose)
    !call analysis2(Aaug,D,R,S,2*global_ndim,nrens,nobs_tot,verbose)
    call pull_aug

  else

    stop 'enkf_analysis: invalid option for is_mod_err'

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
