module m_analysis
  implicit none

contains

subroutine analysis(A, R, E, S, D, innov, ndim, nrens, nrobs, verbose, truncation, mode, &
                    lrandrot, lupdate_randrot, lsymsqrt, inflate, infmult, islocal)

  !=======================================================================
  !  PURPOSE:
  !    Computes the analysed ensemble for A using EnKF or square-root
  !    ensemble Kalman filter formulations.
  !
  !  DESCRIPTION:
  !    This routine supports:
  !      - Classical EnKF perturbation update     (mode 10 / 11 / 12 / 13)
  !      - Deterministic square-root filters      (mode 21 / 22 / 23)
  !
  !    It performs:
  !      1. Pseudo-inversion of SS' + (N-1)R or SS'+EE' depending on mode
  !      2. Construction of X5 / representer matrix depending on data ratio
  !      3. Final ensemble update A <- A + update
  !      4. Inflation (multiplicative or adaptive) and safety clipping
  !
  !  THREAD SAFETY NOTE:
  !    This routine is executed within an OpenMP parallel region during 
  !    local analysis. All local variables and deferred-shape arrays (Z, eig)
  !    are allocated on a per-thread basis, ensuring thread isolation.
  !    Nested OpenMP parallel directives inside Step 4 have been removed 
  !    to avoid thread collisions and severe performance degradation.
  !=======================================================================
  
  use iso_fortran_env, only : dp => real64
  use mod_anafunc
  use m_multa
  use m_ensmean
  use m_ensvar
  implicit none

  !---------------------- Arguments ----------------------
  integer,  intent(in)    :: ndim        
  integer,  intent(in)    :: nrens       
  integer,  intent(in)    :: nrobs       

  real(dp), intent(inout) :: A(ndim,nrens)  
  real(dp), intent(inout) :: R(nrobs,nrobs) 
  real(dp), intent(in)    :: D(nrobs,nrens) 
  real(dp), intent(in)    :: E(nrobs,nrens) 
  real(dp), intent(in)    :: S(nrobs,nrens) 
  real(dp), intent(in)    :: innov(nrobs)   

  logical, intent(in) :: verbose
  real(dp), intent(in) :: truncation      
  integer, intent(in) :: mode             

  logical, intent(in) :: lrandrot         
  logical, intent(in) :: lupdate_randrot  
  logical, intent(in) :: lsymsqrt         

  integer, intent(in) :: inflate          
  real(dp), intent(in) :: infmult         
  logical, intent(in) :: islocal          

  !---------------------- Local variables ----------------------
  real(dp) :: X5(nrens,nrens)     
  real(dp) :: den
  real(dp), parameter :: eps_inv = 1.0e-14_dp 
  real(dp) :: inffac              
  real(dp) :: ave(ndim)           

  real(dp), allocatable :: eig(:), Z(:,:)      
  real(dp), allocatable :: X2(:,:), X3(:,:), Reps(:,:)

  integer  :: nrmin, i, j, iblkmax
  logical  :: lreps               

  external :: dgemm

  lreps = .false.

  if (nrobs <= 0) then
     if (verbose) write(*,*) 'Warning: analysis called with 0 observations. Skipping.'
     return
  end if

  !=======================================================================
  !  STEP 1: PSEUDO-INVERSION OF C
  !=======================================================================

  if (nrobs == 1) then
      nrmin = 1
      allocate(Z(1,1), eig(1))
      
      den = dot_product(S(1,:), S(1,:)) + real(nrens-1, dp)*R(1,1)
      if (den > eps_inv) then
          eig(1) = 1.0_dp / den
      else
          if (verbose) write(*,*) 'Warning: Singular obs at node, skipping update'
          eig(1) = 0.0_dp
      end if
      Z(1,1) = 1.0_dp

  else
      select case (mode)

      case (10)
         call exact_diag_inversion(S, D, X5, nrens, nrobs)

      case (11, 21)
         nrmin = nrobs
         call dgemm('N','T', nrobs, nrobs, nrens, 1.0_dp, S, nrobs, S, nrobs, &
                    real(nrens-1,dp), R, nrobs)
         allocate(Z(nrobs,nrobs), eig(nrobs))
         call eigC(R, nrobs, Z, eig)
         call eigsign(eig, nrobs, truncation)

      case (12, 22)
         nrmin = min(nrobs, nrens)
         allocate(Z(nrobs,nrmin), eig(nrmin))
         call lowrankCinv(S, R, nrobs, nrens, nrmin, Z, eig, truncation)

      case (13, 23)
         nrmin = min(nrobs, nrens)
         allocate(Z(nrobs,nrmin), eig(nrmin))
         call lowrankE(S, E, nrobs, nrens, nrmin, Z, eig, truncation)

      case default
         print *, 'error analysis: Unknown mode: ', mode
         stop
      end select
  end if

  !=======================================================================
  !  STEP 2: GENERATION OF X5 (OR REPRESENTERS)
  !=======================================================================

  select case (mode)

  case (10)
      ! X5 already produced by exact_diag_inversion

  case (11,12,13)
      allocate(X3(nrobs,nrens))
      if (nrobs > 1) then
          call genX3(nrens, nrobs, nrmin, eig, Z, D, X3)
      else
          X3 = D * eig(1)
      end if

      if (2_8*ndim*nrobs < 1_8*nrens*(nrobs+ndim) .and. inflate /= 2) then
          if (verbose) print '(a)', 'analysis: Representer approach is used'
          lreps = .true.
          allocate(Reps(ndim,nrobs))
          call dgemm('N','T', ndim, nrobs, nrens, 1.0_dp, A, ndim, S, nrobs, 0.0_dp, Reps, ndim)
      else
          if (verbose) print '(a)', 'analysis: X5 approach is used'
          call dgemm('T','N', nrens, nrens, nrobs, 1.0_dp, S, nrobs, X3, nrobs, 0.0_dp, X5, nrens)
          do i = 1, nrens
              X5(i,i) = X5(i,i) + 1.0_dp
          end do
      end if

      if (inflate == 2 .and. lreps) then
          call dgemm('T','N', nrens, nrens, nrobs, 1.0_dp, S, nrobs, X3, nrobs, 0.0_dp, X5, nrens)
          do i = 1, nrens
              X5(i,i) = X5(i,i) + 1.0_dp
          end do
      end if

  case (21,22,23)
      call meanX5(nrens, nrobs, nrmin, S, Z, eig, innov, X5)

      allocate(X2(nrmin,nrens))
      call genX2(nrens, nrobs, nrmin, S, Z, eig, X2)

      call X5sqrt(X2, nrobs, nrens, nrmin, X5, lrandrot, lupdate_randrot, mode, lsymsqrt)
  end select

  !=======================================================================
  !  STEP 3: FINAL ENSEMBLE UPDATE
  !=======================================================================

  if (verbose) print '(a)', 'analysis: final update'

  if (lreps) then
     call dgemm('N','N', ndim, nrens, nrobs, 1.0_dp, Reps, ndim, X3, nrobs, 1.0_dp, A, ndim)
     if (.not. islocal) call dumpX3(X3, S, nrobs, nrens)   
  else
     iblkmax = min(ndim, 200)
     call multa(A, X5, ndim, nrens, iblkmax)
     if (.not. islocal) call dumpX5(X5, nrens)            
  end if

  if (verbose) print '(a)', 'analysis: final update done'

  !=======================================================================
  !  STEP 4: INFLATION AND SAFETY CLIPPING
  !=======================================================================
  call ensmean(A, ave, ndim, nrens) 

  if (inflate == 1) inffac = infmult
  if (inflate == 2) then
     call inflationfactor(X5, nrens, inffac)
     inffac = 1.0_dp + (inffac - 1.0_dp)*infmult
  end if

  do j = 1, nrens
     do i = 1, ndim
        if (inflate > 0) A(i,j) = ave(i) + (A(i,j) - ave(i)) * inffac

        ! --- SELECTIVE SAFETY CLIPPING ---
        ! 1. Remove NaNs or extreme overflows
        if (A(i,j) /= A(i,j) .or. abs(A(i,j)) > 1.0e15_dp) then
           A(i,j) = ave(i)
           cycle
        end if

        ! 2. Variable-dependent Clipping (Logic for T and S)
        if (abs(A(i,j)) > 1.0e5_dp) A(i,j) = ave(i) 
     end do
  end do

  !=======================================================================
  !  STEP 5: DEALLOCATIONS
  !=======================================================================
  if (allocated(X2))   deallocate(X2)
  if (allocated(X3))   deallocate(X3)
  if (allocated(eig))  deallocate(eig)
  if (allocated(Z))    deallocate(Z)
  if (allocated(Reps)) deallocate(Reps)

end subroutine analysis

endmodule m_analysis
