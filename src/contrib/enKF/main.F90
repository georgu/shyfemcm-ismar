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
program main

  use mod_mod_states
  use mod_para
  use mod_enkf
  use mod_ens_state
  use mod_mod_err
  use mod_restart , only : ibarcl_rst

  implicit none
  integer :: ndim
  real, allocatable :: Amat(:,:)

  write(*,*) "***"
  write(*,*) "*** Geir's analysis routine, method: ",rmode
  write(*,*) "***"

!----------------------------------------------------
! Opens info file
!----------------------------------------------------
  call read_info

!----------------------------------------------------
! Read the ensemble of model states
!----------------------------------------------------
  call read_ensemble

!----------------------------------------------------
! compute the prior mean and std
!----------------------------------------------------
  call make_mean_std('b')

!----------------------------------------------------
! Read observations and pre-process them
!----------------------------------------------------
  call read_obs

!----------------------------------------------------
! makes D1, E and R, S and d-H*mean(Abk)
!----------------------------------------------------
  call make_matrices

!--------------------------------
! Make the state for the analysis (model error, parameters)
!--------------------------------
  ! define the dimension of the state
  if (ibarcl_rst == 0) then
     ndim = nnkn + 2*nnel*nnlv
  else
     ndim = nnkn + 2*nnel*nnlv + 2*nnkn*nnlv
  end if

  select case(mode_an)

   case(0) !normal call

     allocate(Amat(ndim,nrens))
     call tystate_to_matrix(ibarcl_rst,nrens,ndim,Abk,Amat)

   case(1) !state with model error

     call info_moderr
     call push_aug
     allocate(Amat(2*ndim,nrens))
     call tyqstate_to_matrix(ibarcl_rst,nrens,2*ndim,Abk_aug,Amat)

    case(2) !state with model parameter

     error stop 'TODO'
     ! call load_model_params(npar,nrens_p,pars)
     ! if (nrens_p /= nrens) stop error 'bad parameter ensemble'
     ! allocate(Amat(ndim+npar,nrens))
     ! call typstate_to_matrix(ibarcl_rst,nrens,ndim+npar,pars,Abk,Amat)

  end select

!--------------------------------
! Call the analysis routine
!--------------------------------
  
  select case(mode_an)

   case(0) !normal call

     if (is_local == 0) then

        call analysis(Amat,R,E,S,D1,innov,ndim,nrens,nobs_ok,verbose,&
                  truncation,rmode,lrandrot,lupdate_randrot,lsymsqrt,&
		  inflate,infmult)
	  
        allocate(Aan(nrens))
        call matrix_to_tystate(ibarcl_rst,nrens,ndim,Amat,Aan)
        deallocate(Amat)

        ! save a total X5 for enKS
        call save_X5('global',atime_an)

     elseif (is_local == 1) then

	write(*,*) 'Running local analysis...'
        allocate(Aan(nrens))
        call matrix_to_tystate(ibarcl_rst,nrens,ndim,Amat,Aan)
        deallocate(Amat)

	call local_analysis
	
     end if
   
   case(1) !state with model error

     call matrix_to_tyqstate(ibarcl_rst,nrens,2*ndim,Amat,Abk_aug)
     deallocate(Amat)
     call pull_aug

   case(2) !state with model parameter
     
     !call matrix_to_typstate(ibarcl_rst,nrens,ndim+npar,Amat,pars,Aan)
     !call save_model_params(npar,nrens,pars)
     !deallocate(Amat)

  end select

!--------------------------------
!  check and correct all the out-of-range and nans in analysis
!--------------------------------
  call check_and_correct(Abk,Aan)

!--------------------------------
!  compute the posterior mean and std. This is needed for inflation
!--------------------------------
  call make_mean_std('a')

!--------------------------------
! Save the output in different restart files
!--------------------------------
  call write_ensemble

end program main
