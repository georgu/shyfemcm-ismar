!------------------------------------------------------------------------------
! Run the ensemble Kalman filter analysis using the Evenseen's routines
!------------------------------------------------------------------------------
  program enKF_analysis

  use mod_states
  use mod_enkf
  use m_read_ensemble
  use m_write_ensemble
  implicit none

  ! Analysis parameters
  integer :: rmode = 23 ! mode to run the program (see analysis.F90)
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
! Read observations and makes D, E and R
!----------------------------------------------------
! Better to read D directly and to make it outside.
  call read_obs
  call make_D_E_R

!--------------------------------
! Make S(nobs,nrens), matrix holding HA`, and innov(nobs), innovation vector holding 
! d-H*mean(A) 
  call make_S_innov
!--------------------------------

!--------------------------------
! Call the analysis routine
!--------------------------------
  call analysis(A,R,E,S,D,innov,global_ndim,nrens,nobs,verbose,truncation,rmode,update_randrot)
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
