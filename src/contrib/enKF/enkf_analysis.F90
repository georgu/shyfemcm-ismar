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

  character(len=3) :: nal
  character(len=80) :: filin

!----------------------------------------------------
! Opens info file
!----------------------------------------------------
  call read_info
  call num2str(nanal,nal)
  
!----------------------------------------------------
! Set shyfem variables and init modules
!----------------------------------------------------
  call set_model_params

!----------------------------------------------------
! Read the ensemble
!----------------------------------------------------
  call read_ensemble

!----------------------------------------------------
! Makes the mean of A and saves a restart file with it
! before the analysis
!----------------------------------------------------
  call mean_state
  call std_state('old')

  filin = 'an'//nal//'_mean_state_b.rst'
  call write_state(Am,filin)
  filin = 'an'//nal//'_std_states_b.rst'
  call write_state(Astd_old,filin)

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

   case(0) !no model error
   
    call analysis(A,R,E,S,D1,innov,global_ndim,nrens,nobs_tot,verbose,truncation,rmode,update_randrot)
    !call analysis6c(A,E,S,innov,global_ndim,nrens,nobs_tot,verbose)	!SQRT alg
    !call analysis2(A,D1,R,S,global_ndim,nrens,nobs_tot,verbose)	!ENKF alg

   case(1) !add a model error to the state

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
!  the average state
!--------------------------------
  call mean_state
  call std_state('new')

  filin = 'an'//nal//'_mean_state_a.rst'
  call write_state(Am,filin)
  filin = 'an'//nal//'_std_states_a.rst'
  call write_state(Astd_new,filin)

!--------------------------------
!  state inflation
!--------------------------------
  call inflate_state

!--------------------------------
! Save the output in different restart files
!--------------------------------
  call write_ensemble


  end program enKF_analysis
