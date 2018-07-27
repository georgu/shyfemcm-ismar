!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
!
!------------------------------------------------------------------------------
! Runs the Ensemble Kalman filter or the Square Root analysis by using the 
! Evenseen's routines
!------------------------------------------------------------------------------
! Parameters for the routine analysis.F90
! Different types of analysis need different input
! arguments. However the program computes all the 
! matrices for both enkf and sqrt cases
!
! Remember that all the model variables and the shyfem
! routines LINKED to this code work in real*4, while this
! code is compiled in double precision. This is better for
! the analysis (Topaz does the same).
!------------------------------------------------------------------------------
program enKF_analysis

  use mod_para
  use mod_mod_err
  use mod_enkf
  use mod_ens_state

  implicit none

  write(*,*) "***"
  write(*,*) "*** Geir's analysis routine, method: ",rmode
  write(*,*) "***"

!----------------------------------------------------
! Opens info file
!----------------------------------------------------
  call read_info

!----------------------------------------------------
! Set shyfem variables and init modules
!----------------------------------------------------
  call set_model_params

!----------------------------------------------------
! Read the ensemble of model states
!----------------------------------------------------
  call read_ensemble

!----------------------------------------------------
! compute the prior mean and std
!----------------------------------------------------
  call make_mean_std('b',nanal)

!----------------------------------------------------
! Read observations and pre-process them
!----------------------------------------------------
  call read_obs

!----------------------------------------------------
! makes D1, E and R, S and d-H*mean(A)
!----------------------------------------------------
  call make_matrices

!--------------------------------
! Call the analysis routine
!--------------------------------
  select case(bmod_err)

   case(0) !no model error
   
    call analysis(A,R,E,S,D1,innov,global_ndim,nrens,nobs_tot,verbose,truncation,rmode,update_randrot)
    !call analysis6c(A,E,S,innov,global_ndim,nrens,nobs_tot,verbose)	!SQRT alg
    !call analysis2(A,D1,R,S,global_ndim,nrens,nobs_tot,verbose)	!ENKF alg

   case(1) !add a model error to the state (augmented state. Not working)

    call info_moderr
    call push_aug
    deallocate(A)
    call analysis(Aaug,R,E,S,D1,innov,2*global_ndim,nrens,nobs_tot,&
                  verbose,truncation,rmode,update_randrot)
    call pull_aug

  end select

!--------------------------------
! Local analysis
!--------------------------------
  select case (is_local)

   case(1)	! localisation
 
    continue	!TODO
    !call analysis_nkn
    !call analysis_nel
     
    ! save a total X5 for enKS
    !call save_X5('local',atime) !TODO

   case default  ! no local analysis
       
    ! save a total X5 for enKS
    call save_X5('global',atime)

  end select
   

!--------------------------------
!  compute the posterior mean and std
!--------------------------------
  call make_mean_std('a',nanal)

!--------------------------------
!  state inflation
!--------------------------------
  call inflate_state

!--------------------------------
! Save the output in different restart files
!--------------------------------
  call write_ensemble

end program enKF_analysis
