module m_pseudo1D
   !!
   !! 1D pseudo-random field generator using FFTW3 and the Newton-based
   !! correlation-length calibration. Updated version:
   !!   • fully double precision
   !!   • FFTW-safe
   !!   • robust for small grids (1° resolution, etc.)
   !!   • no longer uses r1 = 3/rx (unsafe for small scales)
   !!
   use iso_fortran_env, only: dp => real64
   use m_random
   use m_newton1D
   use mod_fftw3
   implicit none
contains

   subroutine pseudo1D(A, nx, nrfields, rx, dx, n1)
      implicit none

      !==================== Arguments ======================
      integer,  intent(in)  :: nx
      integer,  intent(in)  :: nrfields
      real(dp), intent(out) :: A(nx,nrfields)
      real(dp), intent(in)  :: rx, dx
      integer,  intent(in)  :: n1

      !==================== Local variables =================
      real(dp) :: r1, r12
      logical  :: cnv
      integer  :: l, i
      real(dp) :: c
      real(dp) :: kappa2, kappa
      real(dp) :: pi2, deltak, summ
      real(dp), dimension(0:n1/2,2) :: fampl
      real(dp), dimension(0:n1/2)   :: phi
      real(dp) :: tt
      logical, save :: diag = .false.

      real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp

      ! FFTW plan (double precision)
      integer(C_INT) :: plan
      complex(dp), dimension(n1/2+1) :: arrayC
      real(dp),    dimension(n1)     :: y

      !==================== Precompute geometry =================
      pi2    = 2.0_dp * pi
      deltak = pi2 / ( real(n1,dp) * dx )
      kappa  = pi2 / ( real(n1,dp) * dx )
      kappa2 = kappa * kappa

      !==================== FFTW plan ===========================
      call dfftw_plan_dft_c2r_1d(plan, n1, arrayC, y, FFTW_ESTIMATE)

      !==================== NEW robust initial guess ==============
      !!
      !! Old: r1 = 3.0_dp / rx   (unsafe for small grids)
      !! New: scale-independent guess inside Newton’s basin
      !! Newton will explore 10 scaled values anyway.
      !!
      r1 = 0.2_dp        ! universal safe guess, independent of rx

      if (diag) print *, 'Initial r1 guess = ', r1
      call newton1D(r1, n1, dx, rx, cnv)
      if (.not. cnv) stop 'pseudo1D: Newton did not converge.'

      if (diag) print *, 'Newton result r1 = ', r1

      r12 = r1*r1

      !==================== Normalization constant =================
      summ = 0.0_dp
      do l = -n1/2+1, n1/2
         summ = summ + exp( -2.0_dp * kappa2 * real(l*l,dp) / r12 )
      end do
      summ = summ - 1.0_dp
      c    = sqrt( 1.0_dp / ( deltak * summ ) )

      !==================== Generate random fields =================
      do i = 1, nrfields

         call random_number(phi)
         phi = pi2 * phi

         do l = 1, n1/2
!$omp simd
            tt = kappa2 * real(l*l,dp) / r12
            fampl(l,1) = exp(-tt) * cos(phi(l)) * sqrt(deltak) * c
            fampl(l,2) = exp(-tt) * sin(phi(l)) * sqrt(deltak) * c
         end do

         fampl(0,1) = 0.0_dp
         fampl(0,2) = 0.0_dp

         arrayC(:) = cmplx(fampl(:,1), fampl(:,2), kind=dp)

         call dfftw_execute(plan)
         A(1:nx,i) = y(1:nx)

      end do

      !==================== Cleanup ==========================
      call dfftw_destroy_plan(plan)

   end subroutine pseudo1D
end module m_pseudo1D
