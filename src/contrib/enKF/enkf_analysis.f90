!------------------------------------------------------------------------------
! Run the ensemble Kalman filter analysis using the Evenseen's routines
!------------------------------------------------------------------------------
  program enKF_analysis

  use mod_enKF

  implicit none

  ! Analysis parameters
  integer :: rmode = 23 ! mode to run the program (see analysis.F90)
  real :: truncation = 0.995 ! truncation of the SVD eigenvalues
  logical :: update_randrot = .true. ! False for local analysis
  logical :: verbose = .true. ! Prints diagnostic output

  integer ne

!----------------------------------------------------
! Opens info file
!----------------------------------------------------
  call read_info
  
!----------------------------------------------------
! Read basin
!----------------------------------------------------
  call read_basin(basfile)

!----------------------------------------------------
! Load the ensemble state A
!----------------------------------------------------
  write(*,*) 'Reading the restart files with the background states...'
  do ne = 0,nens-1
     call rst_read(ne,na,tobs)
     call push_state(ne+1)	!matrix columns start from 1
  end do
  write(*,*) 'Dimension of the model state: ',sdim
  write(*,*) 'Number of ensemble members: ',nens

!----------------------------------------------------
! Makes the mean of A ans saves a restart file with it
! before the analysis
!----------------------------------------------------
  call average_mat(A,Am,sdim,nens)
  call pull_av_state
  call rst_write(-1,na,tobs)	! -1 to write with average backgr label

!----------------------------------------------------
! Read observations and makes D, E and R
!----------------------------------------------------
! Better to read D directly and to make it outside.
  call read_obs
  call make_D_E_R

!--------------------------------
! Make S(nobs,nens), matrix holding HA`, and innov(nobs), innovation vector holding 
! d-H*mean(A) 
  call make_S_innov
!--------------------------------

!--------------------------------
! Call the analysis routine
!--------------------------------
  call analysis(A,R,E,S,D,innov,sdim,nens,nobs,verbose,truncation,rmode,update_randrot)
  write(*,*) 'Analysis done'

!--------------------------------
! Save the output in different restart files
!--------------------------------
  write(*,*) 'Writing the restart files with the analysis states...'
  do ne = 0,nens-1
     call rst_read(ne,na,tobs)
     call pull_state(ne+1)
     call rst_write(ne,na,tobs)
  end do

!--------------------------------
! Save the average state
!--------------------------------
! Warning! The variables not used in the analysis are the last stored,
! i.e., of the last ens state. This must be corrected.
  call average_mat(A,Am,sdim,nens)
  call pull_av_state
  call rst_write(-2,na,tobs)	! -2 to write with average analysis label

  end program enKF_analysis
