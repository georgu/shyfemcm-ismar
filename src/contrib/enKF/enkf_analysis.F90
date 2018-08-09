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

  use mod_mod_states
  use mod_para
  use mod_mod_err
  use mod_enkf
  use mod_ens_state

  implicit none
  integer :: dim_tot
  real, allocatable :: Amat(:,:)

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
! makes D1, E and R, S and d-H*mean(Ashy)
!----------------------------------------------------
  call make_matrices

!--------------------------------
! Make the state for the analysis (model error, parameters)
!--------------------------------
  select case(mode_an)

   case(0) !normal call
 
     allocate(Amat(global_ndim,nrens))
     call tystate_to_matrix(nrens,Ashy,Amat)
     deallocate(Ashy)
     dim_tot = global_ndim

   case(1) !state with model error

     call info_moderr
     call push_aug
     allocate(Amat(2*global_ndim,nrens))
     call tyqstate_to_matrix(nrens,Ashy_aug,Amat)
     deallocate(Ashy_aug)
     dim_tot = 2*global_ndim

    case(2) !state with model parameter

     error stop 'TODO'
     ! call load_model_params(npar,nrens_p,pars)
     ! if (nrens_p /= nrens) stop error 'bad parameter ensemble'
     ! allocate(Amat(global_ndim+npar,nrens))
     ! call typstate_to_matrix(nrens,npar,pars,Ashy,Amat)
     ! deallocate(Ashy)
     ! dim_tot = global_ndim + npar

  end select

!--------------------------------
! Call the analysis routine
!--------------------------------
  call analysis(Amat,R,E,S,D1,innov,dim_tot,nrens,nobs_tot,verbose,&
                  truncation,rmode,update_randrot)
  !call analysis6c(Amat,E,S,innov,dim_tot,nrens,nobs_tot,verbose)	!SQRT alg
  !call analysis2(Amat,D1,R,S,dim_tot,nrens,nobs_tot,verbose)	!ENKF alg
    
!--------------------------------
! Do after... 
!--------------------------------
  select case(mode_an)

   case(0) !normal call
   
     allocate(Ashy(nrens))
     call matrix_to_tystate(nrens,Amat,Ashy)
     deallocate(Amat)

   case(1) !state with model error

     allocate(Ashy_aug(nrens))
     call matrix_to_tyqstate(nrens,Amat,Ashy_aug)
     deallocate(Amat)
     call pull_aug

   case(2) !state with model parameter
     
     !allocate(Ashy(nrens))
     !call matrix_to_typstate(nrens,npar,Amat,pars,Ashy)
     !call save_model_params(npar,nrens,pars)
     !deallocate(Amat)

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
