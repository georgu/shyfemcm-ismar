module m_newton2D

   use iso_fortran_env, only : dp => real64
   use m_newtonfunc2D

   implicit none

contains

   subroutine newton2D(r1, r2, n1, n2, dx, dy, rx, ry, lconv, verbose)

      implicit none

      integer,  intent(in)    :: n1, n2
      logical,  intent(in)    :: verbose

      real(dp), intent(inout) :: r1, r2
      real(dp), intent(in)    :: dx, dy
      real(dp), intent(in)    :: rx, ry

      logical,  intent(out)   :: lconv

      integer, parameter :: max_iter_outer = 10
      integer, parameter :: max_iter_inner = 100

      real(dp), parameter :: tol_rel = 1.0e-5_dp
      real(dp), parameter :: eps     = 1.0e-12_dp
      real(dp), parameter :: rmin    = 1.0e-6_dp
      real(dp), parameter :: rmax    = 1.0_dp

      integer  :: i, j, ulog

      real(dp) :: f, g
      real(dp) :: f1, g1, f2, g2
      real(dp) :: det

      real(dp) :: inc1, inc2
      real(dp) :: err1, err2

      real(dp) :: r1ini, r2ini
      real(dp) :: gamma

      lconv = .false.

      r1ini = max(rmin, min(rmax, r1))
      r2ini = max(rmin, min(rmax, r2))

      ulog = -1

      if (verbose) then
         ulog = 10
         open(ulog, file='newton.cnv', status='replace', action='write')
      endif

      gamma = 1.25_dp

      do j = 1, max_iter_outer

         r1 = real(j, dp) * r1ini
         r2 = real(j, dp) * r2ini

         r1 = max(rmin, min(rmax, r1))
         r2 = max(rmin, min(rmax, r2))

         gamma = gamma - 0.25_dp / real(j, dp)

         do i = 1, max_iter_inner

            call newtonfunc2D( &
                 f, g, &
                 f1, g1, f2, g2, &
                 r1, r2, &
                 n1, n2, dx, dy, rx, ry )

            det = f1*g2 - f2*g1

            ! Protection against singular Jacobian
            if (abs(det) <= eps) exit

            inc1 = (f*g2 - f2*g) / det
            inc2 = (f1*g  - f*g1) / det

            r1 = r1 - gamma*inc1
            r2 = r2 - gamma*inc2

            r1 = max(rmin, min(rmax, r1))
            r2 = max(rmin, min(rmax, r2))

            err1 = inc1 / (abs(r1) + tol_rel)
            err2 = inc2 / (abs(r2) + tol_rel)

            if (verbose) then
               write(ulog,'(i5,10g13.5)') &
                    i, r1, r2, &
                    f, g, &
                    f1, g1, f2, g2, &
                    err1, err2
            endif

            if (abs(err1) + abs(err2) < tol_rel) then

               lconv = .true.

               if (verbose) then
                  write(*,*) 'Newton converged in iteration ', j, i
               endif

               exit

            endif

         enddo

         if (lconv) exit

      enddo

      if (verbose) close(ulog)

      if (.not. lconv) then

         r1 = r1ini
         r2 = r2ini

         write(*,*) 'Newton did not converge.'
         write(*,*) 'Probably an error in input parameters.'
         write(*,*) 'Is dx (and dy) resolving the rx (and ry) scales?'

      endif

   end subroutine newton2D

end module m_newton2D
