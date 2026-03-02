module m_newtonfunc2D
   !!
   !! Provides the 2D functions and their partial derivatives for the Newton solver.
   !! Inputs:  r1, r2, n1, n2, dx, dy, rx, ry
   !! Outputs: f, g, f1=∂f/∂r1, g1=∂g/∂r1, f2=∂f/∂r2, g2=∂g/∂r2
   !!
   use iso_fortran_env, only: dp => real64
   implicit none
contains

   subroutine newtonfunc2D(f, g, f1, g1, f2, g2, r1, r2, n1, n2, dx, dy, rx, ry)
      implicit none

      !----------------------------- arguments
      real(dp), intent(out) :: f, g, f1, g1, f2, g2
      real(dp), intent(in)  :: r1, r2
      integer,  intent(in)  :: n1, n2
      real(dp), intent(in)  :: dx, dy, rx, ry

      !----------------------------- local constants
      real(dp), parameter :: pi     = 3.14159265358979323846264338327950288_dp
      real(dp), parameter :: expm1  = 2.7182818284590452353602874713526625_dp**(-1) ! exp(-1)
      real(dp), parameter :: tiny_r = 1.0e-12_dp  ! safeguard for r1, r2

      !----------------------------- locals
      integer  :: l, p
      real(dp) :: two_pi, kappa, kappa2, lambda, lambda2
      real(dp) :: r1s, r2s
      real(dp) :: e, l2, p2

      !----------------------------- precompute geometry
      two_pi = 2.0_dp * pi
      kappa  = two_pi / ( real(n1, dp) * dx )
      kappa2 = kappa * kappa
      lambda = two_pi / ( real(n2, dp) * dy )
      lambda2= lambda * lambda

      ! Safeguard radii against non-positive values
      r1s = max(r1, tiny_r)
      r2s = max(r2, tiny_r)

      f  = 0.0_dp
      g  = 0.0_dp
      f1 = 0.0_dp
      g1 = 0.0_dp
      f2 = 0.0_dp
      g2 = 0.0_dp

      ! The original loops span p = -n2/2+1..n2/2 and l = -n1/2+1..n1/2,
      ! skipping only the joint (p=0,l=0) mode. Keep that behavior.
      !
      ! Optional SIMD hint; harmless without OpenMP and useful to vectorize.
!$omp simd collapse(2)
      do p = -n2/2 + 1, n2/2
         do l = -n1/2 + 1, n1/2
            if (p == 0 .and. l == 0) cycle

            l2 = real(l*l, dp)
            p2 = real(p*p, dp)

            ! Exponential envelope:
            ! e = exp( -2 * ( kappa^2 l^2 / r1^2 + lambda^2 p^2 / r2^2 ) )
            e  = exp( -2.0_dp * ( kappa2*l2/(r1s*r1s) + lambda2*p2/(r2s*r2s) ) )

            ! Function values
            f  = f  + e * ( cos( kappa * real(l, dp) * rx )   - expm1 )
            g  = g  + e * ( cos( lambda* real(p, dp) * ry )   - expm1 )

            ! Partial derivatives wrt r1 and r2
            ! d/dr1 exp(...) = exp(...) * ( 4 kappa^2 l^2 / r1^3 )
            f1 = f1 + e * ( 4.0_dp * kappa2 * l2 / (r1s*r1s*r1s) ) *  &
                      ( cos( kappa * real(l, dp) * rx )   - expm1 )

            g1 = g1 + e * ( 4.0_dp * kappa2 * l2 / (r1s*r1s*r1s) ) *  &
                      ( cos( lambda* real(p, dp) * ry )   - expm1 )

            ! d/dr2 exp(...) = exp(...) * ( 4 lambda^2 p^2 / r2^3 )
            f2 = f2 + e * ( 4.0_dp * lambda2 * p2 / (r2s*r2s*r2s) ) *  &
                      ( cos( kappa * real(l, dp) * rx )   - expm1 )

            g2 = g2 + e * ( 4.0_dp * lambda2 * p2 / (r2s*r2s*r2s) ) *  &
                      ( cos( lambda* real(p, dp) * ry )   - expm1 )
         end do
      end do

   end subroutine newtonfunc2D

end module m_newtonfunc2D
