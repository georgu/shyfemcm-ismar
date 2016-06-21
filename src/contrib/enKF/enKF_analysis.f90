!------------------------------------------------------------------------------
! Run the ensemble Kalman filter analysis using the Evenseen's routines
!------------------------------------------------------------------------------
  program enKF_analysis

  use basin

! to store the restart vars
  use mod_conz
  use mod_geom_dynamic
  use mod_ts
  use mod_hydro_vel
  use mod_hydro
  use levels, only : nlvdi,nlv
  use mod_restart

  use mod_enKF

  implicit none

  integer ne

!----------------------------------------------------
! Opens info file
!----------------------------------------------------
  call read_info
  
!----------------------------------------------------
! Read basin
!----------------------------------------------------
  call read_basin

!----------------------------------------------------
! Load the ensemble state A
!----------------------------------------------------
  do ne = 1,nens
     call rst_read(ne)
     call push_state(ne)
  end do

!----------------------------------------------------
! Makes the mean of A ans saves a restart file with it
! before the analysis
!----------------------------------------------------
  call average_mat(A,Am,sdim,nens)
  call pull_av_state
  call rst_write(-1)	! -1 to write with average backgr label

!----------------------------------------------------
! Read observations and store in D
!----------------------------------------------------
! Better to read D directly and to make it outside.
  call read_obs
  call make_D_E

!--------------------------------
! Make S(nobs,nens), matrix holding HA`, and innov(nobs), innovation vector holding 
! d-H*mean(A) 
  call make_S_innov
!--------------------------------

!--------------------------------
! Make R(nobs,nobs) matrix holding R (only used if mode=?1 or ?2) (no for low-rank sq root)
!--------------------------------
  allocate(R(nobs,nobs))
  R = 0.02	! R seems used also in the case nobs==1

!--------------------------------
! Call the analysis routine
!--------------------------------
! See the file analysis.F90 for an explanation
  call analysis(A,R,E,S,D,innov,sdim,nens,nobs,.true.,0.99,23,.true.)
  write(*,*) 'Analysis done'

!--------------------------------
! Save the output in different restart files
!--------------------------------
  do ne = 1,nens
     call rst_read(ne)
     call pull_state(ne)
     call rst_write(ne)
  end do

!--------------------------------
! Save the average state
!--------------------------------
! Warning! The variables not used in the analysis are the last stored,
! i.e., of the last ens state. This must be corrected.
  call average_mat(A,Am,sdim,nens)
  call pull_av_state
  call rst_write(-2)	! -2 to write with average analysis label

  end program enKF_analysis
