module m_newton2D
   !!
   !! Double-precision 2D Newton solver with multi-start strategy.
   !! Tries j = 1..max_iter_outer scaled initial guesses and performs
   !! damped Newton iterations for each. Optional OpenMP parallelism
   !! explores the different starts concurrently when verbose=.false.
   !!
   use iso_fortran_env, only: dp => real64
   use m_newtonfunc2D                 ! Must be double precision compatible
   implicit none
contains

   subroutine newton2D(r1, r2, n1, n2, dx, dy, rx, ry, lconv, verbose)
      !!
      !! Solves F(r1,r2) = (0,0) by Newton's method.
      !!
      !! INPUTS:
      !!   n1, n2  - problem/context for newtonfunc2D
      !!   dx, dy  - parameters forwarded to newtonfunc2D
      !!   rx, ry  - additional parameters forwarded to newtonfunc2D
      !!   verbose - if .true., write detailed iteration log to 'newton.cnv'
      !! INOUT:
      !!   r1, r2  - initial guesses on entry, solutions on exit (if converged)
      !! OUTPUT:
      !!   lconv   - .true. if converged for at least one start
      !!
      implicit none

      !----------------------------- arguments
      integer,  intent(in)     :: n1, n2
      real(dp), intent(inout)  :: r1, r2
      real(dp), intent(in)     :: dx, dy, rx, ry
      logical,  intent(out)    :: lconv
      logical,  intent(in)     :: verbose

      !----------------------------- parameters (tunable)
      integer, parameter :: max_iter_outer = 10     ! number of initial scales
      integer, parameter :: max_iter_inner = 100    ! Newton steps per start
      real(dp), parameter :: tol_rel = 1.0e-5_dp    ! relative tol on step
      real(dp), parameter :: rmin = 1.0e-6_dp       ! lower bounds (clipping)
      real(dp), parameter :: rmax = 1.0_dp          ! upper bounds (clipping)
      real(dp), parameter :: eps  = 1.0e-12_dp      ! safeguard for divisions
      logical, parameter :: do_log_header = .true.  ! header in log

      !----------------------------- locals
      real(dp) :: r1ini, r2ini, gamma
      real(dp) :: f, g, f1, g1, f2, g2, det
      real(dp) :: inc1, inc2, err1, err2
      real(dp) :: r1test, r2test
      integer  :: i, j, ulog

      ! Parallel exploration bookkeeping (used when verbose=.false.)
      logical              :: conv_arr(max_iter_outer)
      integer              :: iters_arr(max_iter_outer)
      real(dp)             :: r1_arr(max_iter_outer), r2_arr(max_iter_outer)
      real(dp)             :: res_arr(max_iter_outer)   ! abs(err1)+abs(err2)

      ! Best-candidate selection variables
      integer :: jbest, itmin
      real(dp) :: resmin

      !----------------------------- init
      lconv  = .false.
      r1ini  = max(rmin, min(rmax, r1))   ! sanitize initial guesses
      r2ini  = max(rmin, min(rmax, r2))

      if (verbose) then
         !======================== SERIAL (logging) path =======================
         ulog = 10
         open(ulog, file='newton.cnv', status='replace', action='write')
         if (do_log_header) then
            write(ulog,'(a)') '# i   r1          r2          f           g' // &
                              '           f1          g1          f2          g2' // &
                              '          err1        err2'
         end if

         gamma = 1.25_dp
         do j = 1, max_iter_outer
            r1 = real(j, dp) * r1ini
            r2 = real(j, dp) * r2ini
            r1 = max(rmin, min(rmax, r1))
            r2 = max(rmin, min(rmax, r2))
            gamma = 1.25_dp - 0.25_dp / real(j, dp)   ! as in original

            do i = 1, max_iter_inner
               call newtonfunc2D(f, g, f1, g1, f2, g2, r1, r2, n1, n2, dx, dy, rx, ry)

               det = f1*g2 - f2*g1
               if (abs(det) <= eps) then
                  ! Ill-conditioned Jacobian: cannot take stable step
                  exit
               end if

               ! Newton step via 2x2 inverse (Cramer's rule)
               inc1 = ( f*g2 - f2*g ) / det
               inc2 = ( f1*g - f*g1 ) / det

               ! Damped update + bounding
               r1 = max(rmin, min(rmax, r1 - gamma*inc1))
               r2 = max(rmin, min(rmax, r2 - gamma*inc2))

               ! Relative corrections
               err1 = inc1 / (abs(r1) + eps)
               err2 = inc2 / (abs(r2) + eps)

               write(ulog,'(i5,10(1x,es13.5))') i, r1, r2, f, g, f1, g1, f2, g2, err1, err2

               if (abs(err1) + abs(err2) < tol_rel) then
                  lconv = .true.
                  exit
               end if
            end do

            if (lconv) exit
         end do

         close(ulog)

      else
         !======================== PARALLEL (no logging) path ==================
         conv_arr  = .false.
         iters_arr = huge(0)
         r1_arr    = r1ini
         r2_arr    = r2ini
         res_arr   = huge(1.0_dp)

!$omp parallel default(shared) &
!$omp& private(j, i, gamma, r1test, r2test, f, g, f1, g1, f2, g2, det, inc1, inc2, err1, err2) &
!$omp& if (max_iter_outer > 1)
!$omp do schedule(static)
         do j = 1, max_iter_outer
            r1test = max(rmin, min(rmax, real(j, dp) * r1ini))
            r2test = max(rmin, min(rmax, real(j, dp) * r2ini))
            gamma  = 1.25_dp - 0.25_dp / real(j, dp)

            do i = 1, max_iter_inner
               call newtonfunc2D(f, g, f1, g1, f2, g2, r1test, r2test, n1, n2, dx, dy, rx, ry)

               det = f1*g2 - f2*g1
               if (abs(det) <= eps) exit

               inc1 = ( f*g2 - f2*g ) / det
               inc2 = ( f1*g - f*g1 ) / det

               r1test = max(rmin, min(rmax, r1test - gamma*inc1))
               r2test = max(rmin, min(rmax, r2test - gamma*inc2))

               err1 = inc1 / (abs(r1test) + eps)
               err2 = inc2 / (abs(r2test) + eps)

               if (abs(err1) + abs(err2) < tol_rel) then
                  conv_arr(j)  = .true.
                  iters_arr(j) = i
                  r1_arr(j)    = r1test
                  r2_arr(j)    = r2test
                  res_arr(j)   = abs(err1) + abs(err2)
                  exit
               end if
            end do
         end do
!$omp end do
!$omp end parallel

         ! Select the best converged candidate: fewest iterations, then smallest residual
         if (any(conv_arr)) then
            itmin = huge(0)
            resmin = huge(1.0_dp)
            jbest = 1
            do j = 1, max_iter_outer
               if (conv_arr(j)) then
                  if ( (iters_arr(j) < itmin) .or.  &
                       ((iters_arr(j) == itmin) .and. (res_arr(j) < resmin)) ) then
                     itmin = iters_arr(j)
                     resmin = res_arr(j)
                     jbest = j
                  end if
               end if
            end do
            r1 = r1_arr(jbest)
            r2 = r2_arr(jbest)
            lconv = .true.
         else
            lconv = .false.
            r1 = r1ini
            r2 = r2ini
         end if

      end if

      if (.not. lconv .and. verbose) then
         write(*,*) 'Newton did not converge.'
         write(*,*) 'Probably an error in input parameters.'
         write(*,*) 'Is dx (and dy) resolving the rx (and ry) scales?'
      end if

   end subroutine newton2D

end module m_newton2D
