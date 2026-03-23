!------------------------------------------------------------------------------
! Copyright (C) 2017,
! Marco Bajo, CNR-ISMAR Venice. All rights reserved.
!
! Main driver for Ensemble Kalman Filter (EnKF) and Square-Root analysis.
!
! Workflow:
!  - Read global settings and the input ensemble
!  - Compute background (prior) mean and standard deviation
!  - Read and prepare observations
!  - Build analysis matrices (D1, E, R, S, innovation)
!  - Optionally augment state with model-error components
!  - Execute global or local analysis
!  - Apply boundary/physical-value checks
!  - Compute analysis (posterior) statistics
!  - Write the updated ensemble
!
! All floating-point computations use double precision (dp = real64).
!------------------------------------------------------------------------------

program main
   use iso_fortran_env, only : dp => real64
   use m_set_random_seed2
   use mod_enkf
   use mod_mod_err
   use mod_restart, only : ibarcl_rst
   implicit none

   real(dp), allocatable :: Amat(:,:)
   integer :: istat, ndim


   !--------------------------------------------------------------------------
   ! Init a random seed and save in random_seed.dat, or read it from this file
   !--------------------------------------------------------------------------
   !call init_random_seed_persistent('random_seed.dat',.true.)
   call init_random_seed_constant

   !--------------------------------------------------------------------------
   ! Read configuration, parameters, runtime flags (sets rmode, mode_an, etc.)
   !--------------------------------------------------------------------------
   call read_info

   write(*,*) '***'
   write(*,*) '*** Analysis method code: ', rmode
   write(*,*) '***'

   ! Basic sanity checks on control flags (defensive programming).
   if (is_local /= 0 .and. is_local /= 1) then
      error stop 'main: is_local must be 0 (global) or 1 (local)'
   end if

   !--------------------------------------------------------------------------
   ! Load prior ensemble into module arrays (e.g., Abk)
   !--------------------------------------------------------------------------
   call read_ensemble

   !--------------------------------------------------------------------------
   ! Compute background (prior) statistics
   !--------------------------------------------------------------------------
   call make_mean_std('b')

   !--------------------------------------------------------------------------
   ! Read and preprocess observations
   !--------------------------------------------------------------------------
   call read_obs

   !--------------------------------------------------------------------------
   ! Build analysis matrices: D1, E, R, S, and innovation
   !--------------------------------------------------------------------------
   call make_matrices

   !--------------------------------------------------------------------------
   ! Determine state dimension depending on restart layout.
   !  - ibarcl_rst = 0 → standard state size
   !  - ibarcl_rst ≠ 0 → larger state (baroclinic content included)
   ! Note that the dimension of the model type is always with T/S
   !--------------------------------------------------------------------------
   if (ibarcl_rst == 0) then
      ndim = nnkn + 2*nnel*nnlv
   else
      ndim = nnkn + 2*nnel*nnlv + 2*nnkn*nnlv
   end if

   if (ndim <= 0) error stop 'main: computed ndim <= 0'
   if (nrens <= 0) error stop 'main: nrens must be > 0'
   if (nobs_ok < 0) error stop 'main: nobs_ok must be >= 0'

   !--------------------------------------------------------------------------
   ! Prepare analysis matrix Amat according to analysis mode.
   !
   ! mode_an:
   !   0 → standard state EnKF
   !   1 → augmented state (state + model error)
   !   2 → parameters (not implemented)
   !--------------------------------------------------------------------------
   select case (mode_an)

   case (0)   ! Standard EnKF on state only
      allocate(Amat(ndim, nrens), stat=istat)
      if (istat /= 0) error stop 'main: allocation failed (Amat, mode 0)'

      call tystate_to_matrix(ibarcl_rst, nrens, ndim, Abk, Amat)

   case (1)   ! Augmented state with model error
      call info_moderr
      call push_aug

      allocate(Amat(2*ndim, nrens), stat=istat)
      if (istat /= 0) error stop 'main: allocation failed (Amat, mode 1)'

      call tyqstate_to_matrix(ibarcl_rst, nrens, 2*ndim, Abk_aug, Amat)

   case (2)
      error stop 'main: analysis mode 2 (parameters) is TODO'

   end select

   !--------------------------------------------------------------------------
   ! Run analysis (global or local) and convert matrix back to state arrays.
   !--------------------------------------------------------------------------
   select case (mode_an)

   case (0)   ! Standard EnKF
      if (is_local == 0) then
         ! Global analysis on Amat
         call analysis(Amat, R, E, S, D1, innov, ndim, nrens, nobs_ok, verbose, &
                       truncation, rmode, lrandrot, lupdate_randrot, lsymsqrt, &
                       inflate, infmult)

         call matrix_to_tystate(ibarcl_rst, nrens, ndim, Amat, Aan)
         deallocate(Amat)

         call save_X5('global', atime_an)
      else
         ! Local analysis: convert first, then run location-wise update
         write(*,*) 'Running local analysis...'

         call matrix_to_tystate(ibarcl_rst, nrens, ndim, Amat, Aan)
         deallocate(Amat)

         call local_analysis
      end if

   case (1)   ! Augmented state path
      call matrix_to_tyqstate(ibarcl_rst, nrens, 2*ndim, Amat, Abk_aug)
      deallocate(Amat)
      call pull_aug

   case (2)
      ! Not implemented
      continue

   end select

   !--------------------------------------------------------------------------
   ! Apply boundary-condition corrections and physical-range checks
   !--------------------------------------------------------------------------
   call bc_val_check_correct

   !--------------------------------------------------------------------------
   ! Compute analysis (posterior) statistics
   !--------------------------------------------------------------------------
   call make_mean_std('a')

   !--------------------------------------------------------------------------
   ! Write the updated ensemble to output/restart
   !--------------------------------------------------------------------------
   call write_ensemble

end program main
