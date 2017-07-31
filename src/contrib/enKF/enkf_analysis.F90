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

  use mod_para
  use mod_mod_err
  use mod_enkf
  use mod_ens_state

  implicit none

  integer :: date,time !date and time in the rst files

!----------------------------------------------------
! Opens info file
!----------------------------------------------------
  call read_info
  
!----------------------------------------------------
! Set shyfem variables and init modules
!----------------------------------------------------
  call set_model_params

!----------------------------------------------------
! Read the ensemble
!----------------------------------------------------
  call read_ensemble(date,time)

!----------------------------------------------------
! Makes the mean of A and saves a restart file with it
! before the analysis
!----------------------------------------------------
  call average_mat(date,time,-1)

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
! Note that the ENKF scheme needs D1, the observation innovations, not the
! perturbed measurements. See analysis2.F90, not analysis.F90, for a correct
! description.
  ! Decide if use an augmented state with the model error
  !
  select case(bmod_err)

   case(0)
   
    call analysis(A,R,E,S,D1,innov,global_ndim,nrens,nobs_tot,verbose,truncation,rmode,update_randrot)
    !call analysis6c(A,E,S,innov,global_ndim,nrens,nobs_tot,verbose)	!SQRT alg
    !call analysis2(A,D1,R,S,global_ndim,nrens,nobs_tot,verbose)	!ENKF alg

    call check_values

   case(1)

    call init_moderr
    call push_aug(A)
    deallocate(A)
    call analysis(Aaug,R,E,S,D1,innov,2*global_ndim,nrens,nobs_tot,&
                  verbose,truncation,rmode,update_randrot)
    allocate(A(nrens))
    call pull_aug(A)

  end select

!--------------------------------
! Save the output in different restart files
!--------------------------------
  call write_ensemble(date,time)

!--------------------------------
! Save the average state
!--------------------------------
  call average_mat(date,time,-2)

  end program enKF_analysis
