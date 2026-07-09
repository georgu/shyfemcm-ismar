module m_newton1D

   use iso_fortran_env, only : dp => real64
   use m_newtonfunc1D

   implicit none

contains

   subroutine newton1D(r1, n1, dx, rx, lconv)

      implicit none

      integer,  intent(in)    :: n1
      real(dp), intent(inout) :: r1
      real(dp), intent(in)    :: dx
      real(dp), intent(in)    :: rx
      logical,  intent(out)   :: lconv

      integer, parameter :: max_iter_outer = 10
      integer, parameter :: max_iter_inner = 100

      real(dp), parameter :: tol_rel = 1.0e-5_dp
      real(dp), parameter :: rmin    = 1.0e-6_dp
      real(dp), parameter :: rmax    = 1.0_dp
      real(dp), parameter :: eps     = 1.0e-12_dp

      logical, parameter :: verbose = .false.

      real(dp) :: f
      real(dp) :: f1
      real(dp) :: inc1
      real(dp) :: err1
      real(dp) :: gamma
      real(dp) :: r1ini

      integer :: i
      integer :: j

      lconv = .false.

      r1ini = max(rmin, min(rmax, r1))

      do j = 1, max_iter_outer

         r1 = real(j, dp) * r1ini
         r1 = max(rmin, min(rmax, r1))

         gamma = 1.25_dp - 0.25_dp / real(j, dp)

         do i = 1, max_iter_inner

            call newtonfunc1D(f, f1, r1, n1, dx, rx)

            ! Protection against singular derivative
            if (abs(f1) <= eps) exit

            inc1 = f / f1

            r1 = r1 - gamma * inc1

            r1 = max(rmin, min(rmax, r1))

            err1 = inc1 / (abs(r1) + eps)

            if (abs(err1) < tol_rel) then

               lconv = .true.

               if (verbose) then
                  write(*,'(a,i3,a,i4)') &
                     'Newton converged at start ', j, &
                     ' after iterations ', i
               endif

               exit

            endif

         enddo

         if (lconv) exit

      enddo

      if (.not. lconv) then

         if (verbose) then
            write(*,*) 'Newton did not converge'
         endif

         r1 = r1ini

      endif

   end subroutine newton1D

end module m_newton1D
