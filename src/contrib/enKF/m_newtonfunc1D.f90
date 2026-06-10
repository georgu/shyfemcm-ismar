module m_newtonfunc1D
   !!
   !! Provides the function and its derivative for the 1D Newton solver:
   !!   Inputs:  r1, n1, dx, rx
   !!   Outputs: f(r1), f1(r1) = d f / d r1
   !!
   use iso_fortran_env, only: dp => real64
   implicit none
contains

   subroutine newtonfunc1D(f, f1, r1, n1, dx, rx)
      !!
      !! Computes:
      !!   f(r1)  = sum_{l ≠ 0} exp(-2 * kappa^2 * l^2 / r1^2) *
      !!             [ cos(kappa * l * rx) - exp(-1) ]
      !!   f1(r1) = d f / d r1 (analytical)
      !!
      !! where kappa = 2*pi / (n1*dx), and l runs from -n1/2+1 to n1/2, l≠0.
      !!
      implicit none

      !----------------------------- arguments
      real(dp), intent(out) :: f, f1     ! function value and derivative
      real(dp), intent(in)  :: r1        ! parameter (assumed > 0)
      integer,  intent(in)  :: n1        ! grid count (assumed > 0)
      real(dp), intent(in)  :: dx, rx    ! spacing and coordinate

      !----------------------------- local constants
      real(dp), parameter :: pi     = 3.14159265358979323846264338327950288_dp
      real(dp), parameter :: expm1  = 2.7182818284590452353602874713526625_dp**(-1) ! exp(-1)
      real(dp), parameter :: tiny_r = 1.0e-12_dp  ! safeguard for r1

      !----------------------------- locals
      integer  :: l
      real(dp) :: two_pi, kappa, kappa2, r1s, e, l2

      !----------------------------- precompute geometry
      two_pi = 2.0_dp * pi
      ! Use real(n1,dp) to avoid integer division; n1>0 is expected.
      kappa  = two_pi / ( real(n1, dp) * dx )
      kappa2 = kappa * kappa

      ! Safeguard r1 against zero/negative values
      r1s = max(r1, tiny_r)

      f  = 0.0_dp
      f1 = 0.0_dp

      ! Sum over integer modes, skipping l=0
      ! The loop bounds match your original: l = -n1/2+1 .. n1/2.
      ! Optional SIMD hint helps the compiler vectorize.
!$omp simd
      do l = -n1/2 + 1, n1/2
         if (l == 0) cycle
         l2 = real(l*l, dp)
         ! e = exp( -2 * kappa^2 * l^2 / r1^2 )
         e  = exp( -2.0_dp * kappa2 * l2 / (r1s*r1s) )

         ! f contribution
         f  = f  + e * ( cos( kappa * real(l, dp) * rx ) - expm1 )

         ! df/dr1 contribution:
         !   d/dr1 exp(-2 k^2 l^2 / r1^2) = exp(...) * (4 k^2 l^2 / r1^3)
         !   multiplied by the same bracket [cos(...) - exp(-1)]
         f1 = f1 + e * ( 4.0_dp * kappa2 * l2 / (r1s*r1s*r1s) ) *  &
                   ( cos( kappa * real(l, dp) * rx ) - expm1 )
      end do

   end subroutine newtonfunc1D

end module m_newtonfunc1D
