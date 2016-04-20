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

  !tmp to subtitute with modules
  !include 'param.h'
  !include 'nlevel.h'
  !include 'hydro.h'
  !include 'geom_dynamic.h'
  !include 'basin.h'
  !include 'ts.h'
  !include 'hydro_vel.h'
  !include 'conz.h'

!----------------------------------------------------
! Opens info file
!----------------------------------------------------
  call read_info
  
!----------------------------------------------------
! Read basin
!----------------------------------------------------
  call read_basin

!----------------------------------------------------
! Load the initial ensemble state A
!----------------------------------------------------
  do ne = 1,nrens
     call rst_read(23,trim(rstfile(ne)))
     call store_state(ne)
  end do

!----------------------------------------------------
! Make the mean state of A
!----------------------------------------------------
  call average_mat(A,Am,xdim,nrens)

!----------------------------------------------------
! Read observations and store in D
!----------------------------------------------------
  call read_obs(tobs,nanl,baseobs)



  end program enKF_analysis
