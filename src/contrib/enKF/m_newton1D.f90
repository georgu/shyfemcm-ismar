module m_newton1D
   !!
   !! Double-precision 1D Newton solver with multi-start strategy.
   !! Tries j = 1..max_iter_outer scaled initial guesses and performs
   !! damped Newton iterations for each. Optional OpenMP parallelism
   !! explores the different starts concurrently.
   !!
   use iso_fortran_env, only: dp => real64
   use m_newtonfunc1D                 ! Must be double precision compatible
   implicit none
contains

   subroutine newton1D(r1, n1, dx, rx, lconv)
      !!
      !! INPUTS:
      !!   n1      - problem/context for newtonfunc1D
      !!   dx, rx  - parameters forwarded to newtonfunc1D
      !! INOUT:
      !!   r1      - initial guess on entry, solution on exit (if converged)
      !! OUTPUT:
      !!   lconv   - .true. if converged for at least one start
      !!
      implicit none

      !----------------------------- arguments
      integer,  intent(in)     :: n1
      real(dp), intent(inout)  :: r1
      real(dp), intent(in)     :: dx, rx
      logical,  intent(out)    :: lconv

      !----------------------------- parameters (tunable)
      integer, parameter :: max_iter_outer = 10     ! number of initial scales
      integer, parameter :: max_iter_inner = 100    ! Newton steps per start
      real(dp), parameter :: tol_rel = 1.0e-5_dp    ! relative tolerance
      real(dp), parameter :: rmin = 1.0e-6_dp       ! lower bound on r
      real(dp), parameter :: rmax = 1.0_dp          ! upper bound on r
      real(dp), parameter :: eps  = 1.0e-12_dp      ! safeguard for divisions
      logical, parameter :: verbose = .false.       ! set .true. for prints

      !----------------------------- locals
      real(dp) :: r1ini, gamma, f, f1, inc1, err1, rtest
      integer  :: i, j

      ! Parallel exploration bookkeeping
      logical              :: conv_arr(max_iter_outer)
      integer              :: iters_arr(max_iter_outer)
      real(dp)             :: r_arr(max_iter_outer)
      real(dp)             :: fres_arr(max_iter_outer)

      ! Best-candidate selection variables
      integer :: jbest, itmin
      real(dp) :: fmin

      !----------------------------- init
      lconv  = .false.
      r1ini  = max(rmin, min(rmax, r1))   ! sanitize initial guess

      conv_arr  = .false.
      iters_arr = huge(0)                 ! largest representable integer
      r_arr     = r1ini
      fres_arr  = huge(1.0_dp)            ! largest representable real(dp)

      !===========================================================
      ! Explore j=1..max_iter_outer starts (optionally in parallel)
      !===========================================================
!$omp parallel default(shared) private(j, i, gamma, rtest, f, f1, inc1, err1) if (max_iter_outer>1)
!$omp do schedule(static)
      do j = 1, max_iter_outer
         rtest = real(j, dp) * r1ini                  ! scaled initial guess
         rtest = max(rmin, min(rmax, rtest))
         gamma = 1.25_dp - 0.25_dp / real(j, dp)     ! damping as in original

         do i = 1, max_iter_inner
            call newtonfunc1D(f, f1, rtest, n1, dx, rx)

            ! Guard against zero/ill-conditioned derivative
            if (abs(f1) <= eps) then
               exit  ! cannot take a reliable Newton step
            else
               inc1 = f / f1
            end if

            ! Damped Newton update with bounds
            rtest = rtest - gamma * inc1
            rtest = max(rmin, min(rmax, rtest))

            ! Relative correction measure
            err1  = inc1 / (abs(rtest) + eps)

            if (abs(err1) < tol_rel) then
               conv_arr(j)   = .true.
               iters_arr(j)  = i
               r_arr(j)      = rtest
               fres_arr(j)   = abs(f)
               exit
            end if
         end do
      end do
!$omp end do
!$omp end parallel

      !===========================================================
      ! Pick the "best" converged candidate:
      !   1) fewest iterations, 2) smallest |f|, else keep original r1
      !===========================================================
      if (any(conv_arr)) then
         itmin = huge(0)
         fmin  = huge(1.0_dp)
         jbest = 1
         do j = 1, max_iter_outer
            if (conv_arr(j)) then
               if ( (iters_arr(j) < itmin) .or.  &
                    ((iters_arr(j) == itmin) .and. (fres_arr(j) < fmin)) ) then
                  itmin = iters_arr(j)
                  fmin  = fres_arr(j)
                  jbest = j
               end if
            end if
         end do

         r1    = r_arr(jbest)
         lconv = .true.

         if (verbose) then
            write(*,'(a,i3,a,i4,a,es12.5)') ' Newton converged (best j= ', jbest,  &
                 ') in ', itmin, ' iters, |f|=', fmin
         end if
      else
         if (verbose) write(*,*) 'Newton did not converge'
         lconv = .false.
         ! r1 remains the sanitized initial value
      end if

   end subroutine newton1D

end module m_newton1D
