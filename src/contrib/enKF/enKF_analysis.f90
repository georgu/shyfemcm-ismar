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
! TODO: the init state should contain all the variables of the restart file
  do ne = 1,nens
     call rst_read(ne)
     call push_state(ne)
  end do

!----------------------------------------------------
! Makes the mean of A
!----------------------------------------------------
  call average_mat(A,Am,sdim,nens)

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
  call analysis(A,R,E,S,D,innov,ndim,nens,nobs,verbose,truncation,mode,update_randrot)
  write(*,*) 'Analysis done'

!--------------------------------
! Save the output in different restart files
!--------------------------------
  do ne = 1,nens
     call rst_read(ne)
     call pull_state(ne)
     call rst_write(ne)
  end do

  end program enKF_analysis
