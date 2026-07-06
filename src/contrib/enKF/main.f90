!------------------------------------------------------------------------------
! Copyright (C) 2017,
! Marco Bajo, CNR-ISMAR Venice. All rights reserved.
!
! PURPOSE:
!   Central execution orchestrator for the Ensemble Kalman Filter (EnKF forward pass).
!   Ingests ensemble ocean backgrounds, streams concurrent multivariate data channels,
!   solves the state-space update (Global matrix or OpenMP Localized blocks), and 
!   manages state-augmentation pipelines for parameter estimation.
! ==========================================================================
program main
   use iso_fortran_env, only : dp => real64
   use m_set_random_seed2
   use mod_enkf
   use mod_mod_err
   use mod_restart
   use m_analysis
   use mod_pestimation
   use mod_ens_state
   implicit none

   real(dp), allocatable :: Amat(:,:)
   integer :: istat, ndim
   integer :: ndim_pe = 0

   ! 1. Initialize persistent random generation seed sequence
   call init_random_seed_persistent('random_seed.dat', .true.)

   ! 2. Ingest centralized framework parameters and execution switches
   call read_info

   write(*,*) '========================================'
   write(*,'(A,I4)') ' EnKF Forward Analysis Solver | Mode: ', rmode
   write(*,*) '========================================'

   ! Ingest active parameter estimation configuration maps if active
   call read_pe_files(mode_an, ndim_pe)

   ! Ingest standard background restart structures for all active members
   call read_ensemble

   ! Compute and cache baseline prior ensemble spread statistics
   call make_mean_std('b')

   ! Ingest and screen concurrent multivariate sensor arrays
   call read_observations

   ! Compile centralized observation matrices (R, E, S, D, D1 layers)
   call make_matrices

   ! 3. Map physical state layout constraints to compute dimension 'ndim'
   if (ibarcl_rst == 0) then
      ndim = nnkn + 2*nnel*nnlv ! Barotropic layout size
   else
      ndim = nnkn + 2*nnel*nnlv + 2*nnkn*nnlv ! Baroclinic layout size
   end if

   if (ndim <= 0)  error stop 'ERROR: Main: Physical state dimension ndim <= 0.'
   if (nrens <= 0) error stop 'ERROR: Main: Active ensemble parameter size nrens <= 0.'
   if (nobs_ok < 0) error stop 'ERROR: Main: Verified observation count rows nobs_ok < 0.'

   ! Append active augmentation parameters dimension bounds
   ndim = ndim + ndim_pe

   ! 4. Structure and materialize state-space tracking matrix 'Amat'
   select case (mode_an)
   case (0)   
      ! Standard Filter State tracking matrix allocation
      allocate(Amat(ndim, nrens), stat=istat)
      if (istat /= 0) error stop 'ERROR: Main: Memory allocation failure on Amat (Mode 0).'
      call tystate_to_matrix(ibarcl_rst, nrens, ndim, Abk, ndim_pe, Amat=Amat)

   case (1)   
      ! Augmented State tracking matrix allocation (Includes model error biases)
      call info_moderr
      call push_aug
      allocate(Amat(2*ndim, nrens), stat=istat)
      if (istat /= 0) error stop 'ERROR: Main: Memory allocation failure on Amat (Mode 1).'
      call tyqstate_to_matrix(ibarcl_rst, nrens, 2*ndim, Abk_aug, Amat)

   case (2)  
      ! Parameter Estimation State tracking matrix allocation
      allocate(Amat(ndim, nrens), stat=istat)
      if (istat /= 0) error stop 'ERROR: Main: Memory allocation failure on Amat (Mode 2).'
      call tystate_to_matrix(ibarcl_rst, nrens, ndim, Abk, ndim_pe, pe_mat, Amat)
   end select

   ! 5. Launch the targeted Subspace Analysis Solver Path
   select case (mode_an)
   case (0)   
      if (is_local == 0) then
         write(*,*) ' Executing Global Subspace Analysis Update...'
         call analysis(Amat, R, E, S, D1, innov, ndim, nrens, nobs_ok, verbose, &
                       truncation, rmode, lrandrot, lupdate_randrot, lsymsqrt, &
                       inflate, infmult, .false.)

         call matrix_to_tystate(ibarcl_rst, nrens, ndim, Amat, ndim_pe, A=Aan)
         deallocate(Amat)
         
         ! Archive transition weights for downstream Smoother execution blocks
         call save_X5(atime_an)
      else
         write(*,*) ' Executing Localized Subspace Analysis Update loops...'
         call matrix_to_tystate(ibarcl_rst, nrens, ndim, Amat, ndim_pe, A=Aan)
         deallocate(Amat)
         call local_analysis
      end if

   case (1)   
      write(*,*) ' Executing Augmented Model Error Subspace Analysis Update...'
      call analysis(Amat, R, E, S, D1, innov, 2*ndim, nrens, nobs_ok, verbose, &
                    truncation, rmode, lrandrot, lupdate_randrot, lsymsqrt, &
                    inflate, infmult, .false.)

      call matrix_to_tyqstate(ibarcl_rst, nrens, 2*ndim, Amat, Abk_aug)
      deallocate(Amat)
      call pull_aug

   case (2)  
      if (is_local == 0) then
         write(*,*) ' Executing Global Parameter Estimation Update...'
         call analysis(Amat, R, E, S, D1, innov, ndim, nrens, nobs_ok, verbose, &
                       truncation, rmode, lrandrot, lupdate_randrot, lsymsqrt, &
                       inflate, infmult, .false.)

         call matrix_to_tystate(ibarcl_rst, nrens, ndim, Amat, ndim_pe, pe_mat, A=Aan)
         deallocate(Amat)
         call save_X5(atime_an)
      else
         write(*,*) ' WARNING: Local Parameter Estimation requested. Validating convergence boundaries.'
         call matrix_to_tystate(ibarcl_rst, nrens, ndim, Amat, ndim_pe, pe_mat, A=Aan)
         deallocate(Amat)
         call local_analysis
      end if
      call check_pe_bounds(ndim_pe)
   end select

   ! 6. Run physical Quality Control validation screening checks
   call val_check_correct

   ! 7. Compute and cache posterior analytical ensemble spread statistics
   call make_mean_std('a')

   ! 8. Serialize updated member tracking structures to output storage restarts
   call write_ensemble

   ! If active, export revised parameter estimation state sets to external ledger
   call write_pe_files(ndim_pe)

   write(*,*) '========================================'
   write(*,*) ' EnKF Forward Cycle Executed Successfully.'
   write(*,*) '========================================'

end program main
