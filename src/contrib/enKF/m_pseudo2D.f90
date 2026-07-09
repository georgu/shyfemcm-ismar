module m_pseudo2D
   !
   !
   use iso_fortran_env, only : dp => real64
   use mod_fftw3
   use m_newton2D
   implicit none
contains

!======================================================================
subroutine pseudo2D(Amat, nx, ny, lde, rx, ry, dx, dy, n1, n2, theta, verbose)
!======================================================================
   !
   ! Generates lde two‑dimensional pseudo‑random fields with an
   ! anisotropic Gaussian spectrum following the Evensen (1994) method.
   !
   ! This is an optimized rewrite of the legacy version, keeping
   ! identical output statistics and behaviour.
   !
   logical, intent(in) :: verbose
   integer, intent(in) :: nx, ny          ! physical domain size
   integer, intent(in) :: lde             ! number of fields to generate
   real(dp), intent(out) :: Amat(nx,ny,lde)
   real(dp), intent(in)  :: rx, ry        ! decorrelation lengths
   real(dp), intent(in)  :: dx, dy        ! grid spacing
   real(dp), intent(in)  :: theta         ! rotation angle (deg, CCW)
   integer, intent(in)   :: n1, n2        ! FFT grid sizes (>= nx, >= ny)

   integer(kind=8) :: plan                ! FFTW plan handle
   real(dp), allocatable :: y(:,:)        ! real output field
   complex(dp), allocatable :: x(:,:)     ! complex spectrum
   real(dp) :: r1, r2, c
   real(dp) :: kappa, kappa2, lambda, lambda2
   real(dp) :: pi2, deltak, summ
   real(dp) :: a11tmp, a22tmp, a11, a22, a12, torad
   logical :: cnv
   integer :: i, m, j, l, p
   real(dp), parameter :: pi = 3.141592653589793238_dp

   !------------------ Basic checks
   if (lde < 1) stop 'pseudo2D: lde < 1'
   if (rx <= 0.0_dp) stop 'pseudo2D: rx <= 0'
   if (ry <= 0.0_dp) stop 'pseudo2D: ry <= 0'
   if (n1 < nx) stop 'pseudo2D: n1 < nx'
   if (n2 < ny) stop 'pseudo2D: n2 < ny'

   !------------------ Allocate FFTW‑compatible arrays
   !
   !   Real array : y(0:n1-1, 0:n2-1)
   !   Complex array (half spectrum): x(0:n1/2, 0:n2-1)
   !
   allocate(y(0:n1-1, 0:n2-1))
   allocate(x(0:n1/2, 0:n2-1))

   !------------------ Spectral geometry
   pi2 = 2.0_dp*pi
   deltak = (pi2*pi2) / ( real(n1*n2,dp)*dx*dy )

   kappa  = pi2 / ( real(n1,dp)*dx )
   kappa2 = kappa*kappa
   lambda = pi2 / ( real(n2,dp)*dy )
   lambda2 = lambda*lambda

   !------------------ Create FFTW plan once (safe and fast)
   call dfftw_plan_dft_c2r_2d(plan, n1, n2, x, y, FFTW_ESTIMATE)

   !------------------ Initial guess as in old version
   r1 = 3.0_dp/rx
   r2 = 3.0_dp/ry
   if (verbose) write(*,'(a,2f10.4)') 'Calling Newton with initial r1,r2=', r1, r2

   ! Find calibrated correlation radii
   call newton2D(r1,r2,n1,n2,dx,dy,rx,ry,cnv,verbose)
   if (.not.cnv) stop 'pseudo2D: Newton did not converge'

   !------------------ Normalisation constant (kept identical)
   summ = 0.0_dp
   do p = -n2/2+1, n2/2
      do l = -n1/2+1, n1/2
         summ = summ + exp( -2.0_dp * ( kappa2*real(l*l,dp)/(r1*r1) +  &
                                       lambda2*real(p*p,dp)/(r2*r2) ) )
      end do
   end do
   summ = summ - 1.0_dp
   c = sqrt( 1.0_dp / (deltak*summ) )

   !------------------ Rotation matrix in frequency space
   a11tmp = 1.0_dp/(r1*r1)
   a22tmp = 1.0_dp/(r2*r2)
   torad  = -pi/180.0_dp    ! negative = anticlockwise as in old code

   a11 = a11tmp*cos(theta*torad)**2 + a22tmp*sin(theta*torad)**2
   a22 = a11tmp*sin(theta*torad)**2 + a22tmp*cos(theta*torad)**2
   a12 = (a22tmp - a11tmp)*cos(theta*torad)*sin(theta*torad)

!======================================================================
   do j=1, lde

      ! Generate spectrum (preserves exactly old behaviour)
      call wave_amp(n1,n2,pi2,a11,a12,a22,kappa,kappa2,lambda,lambda2, &
                     deltak,c,x)

      ! Inverse FFT: complex spectrum -> real field
      call dfftw_execute_dft_c2r(plan, x, y)

      ! Extract the physical subdomain (nx,ny)
      do m=1,ny
         do i=1,nx
            Amat(i,m,j) = y(i-1,m-1)
         end do
      end do

   end do

   !------------------ Cleanup
   call dfftw_destroy_plan(plan)
   deallocate(x,y)

end subroutine pseudo2D
!======================================================================



!======================================================================
subroutine wave_amp(n1,n2,pi2,a11,a12,a22,kappa,kappa2,lambda,lambda2, &
                     deltak,c,x)
!======================================================================
   !
   ! Constructs the anisotropic Gaussian spectrum with random phases.
   !
   integer, intent(in) :: n1, n2
   real(dp), intent(in) :: pi2, a11, a12, a22
   real(dp), intent(in) :: kappa, kappa2, lambda, lambda2
   real(dp), intent(in) :: deltak, c
   complex(dp), intent(inout) :: x(0:n1/2,0:n2-1)

   real(dp), allocatable :: phi(:,:)
   real(dp), allocatable :: fampl(:,:,:)
   real(dp) :: e
   integer :: p, l

   !------------------ Allocate exact n2 frequencies:
   ! p = -n2/2 ... n2/2 - 1  (exactly n2 points)
   allocate(phi(0:n1/2, -n2/2 : n2/2-1))
   allocate(fampl(0:n1/2, -n2/2 : n2/2-1, 2))

   !------------------ Draw random phases in [0, 2π)
   call random_number(phi)
   phi = pi2 * phi

   !------------------ Build anisotropic Gaussian envelope
   do p = -n2/2, n2/2-1
      do l = 0, n1/2
         e = exp( - ( a11*kappa2*l*l                                &
                 + 2.0_dp*a12*kappa*lambda*real(l*p,dp)             &
                 + a22*lambda2*p*p ) )
         fampl(l,p,1) = e * cos(phi(l,p)) * sqrt(deltak) * c
         fampl(l,p,2) = e * sin(phi(l,p)) * sqrt(deltak) * c
      end do
   end do

   ! Zero the DC component (required for zero‑mean fields)
   fampl(0,0,1) = 0.0_dp
   fampl(0,0,2) = 0.0_dp

   !------------------ Pack to FFTW "half‑spectrum" layout
   !
   ! First half-columns (positive frequencies)
   do p = 0, n2/2-1
      x(:,p) = cmplx( fampl(:,p,1), fampl(:,p,2), dp )
   end do

   ! Second half-columns (negative frequencies wrapped by FFTW)
   do p = n2/2, n2-1
      x(:,p) = cmplx( fampl(:, p-n2, 1), fampl(:, p-n2, 2), dp )
   end do

   deallocate(phi, fampl)
end subroutine wave_amp
!======================================================================

end module m_pseudo2D
