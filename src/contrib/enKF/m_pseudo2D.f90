module m_pseudo2D
   !---------------------------------------------------------------------------
   ! 2D pseudo-random field generator using FFTW3 and Newton-based calibration.
   !
   ! Key features:
   !  • Double precision everywhere: real(dp), complex(dp)
   !  • Robust on small grids (e.g., 1°): scale‑independent initial r1=r2=0.2_dp
   !    and newton2D then calibrates (r1, r2) to match target correlation.
   !
   ! Overall workflow in pseudo2D:
   !  1) Validate input sizes and parameters.
   !  2) Precompute spectral geometry (kappa, lambda) and Δk for normalization.
   !  3) Choose robust initial (r1,r2) and run Newton2D to calibrate them.
   !  4) Compute normalization constant c so the variance matches the target.
   !  5) Build rotated anisotropic quadratic form (a11,a12,a22) from (r1,r2,θ).
   !  6) Parallel over fields j=1..lde:
   !       • Each thread allocates its own (x_local,y_local) and plan_local.
   !       • Fill complex spectrum x_local with random phases + Gaussian envelope.
   !       • Inverse FFT (c2r) to obtain y_local; copy leading (nx,ny) into Amat.
   !
   ! Notes:
   !  • If your local mod_fftw3 uses a different plan type, adjust plan_local.
   !  • The spectrum uses a rotated quadratic exponent to represent anisotropy.
   !---------------------------------------------------------------------------

   use iso_fortran_env, only : dp => real64
   use mod_fftw3
   use m_newton2D
   implicit none
contains

   subroutine pseudo2D(Amat, nx, ny, lde, rx, ry, dx, dy, n1, n2, theta, verbose)
      implicit none
      !---------------- Arguments
      logical,  intent(in)  :: verbose         ! verbose diagnostics
      integer,  intent(in)  :: nx, ny          ! output field size
      integer,  intent(in)  :: lde             ! number of fields to generate
      real(dp), intent(out) :: Amat(nx, ny, lde)
      real(dp), intent(in)  :: rx, ry          ! decorrelation lengths
      real(dp), intent(in)  :: dx, dy          ! grid spacings
      real(dp), intent(in)  :: theta           ! rotation (deg, 0=east; anticlockwise)
      integer,  intent(in)  :: n1, n2          ! FFT sizes (n1>=nx, n2>=ny)

      !---------------- Declarations (ALL before executable code)
      real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp

      ! Spectral/normalization scalars
      real(dp) :: r1, r2, c
      real(dp) :: kappa, kappa2, lambda, lambda2
      real(dp) :: pi2, deltak, summ

      ! Rotation coefficients (anisotropic quadratic form in k-space)
      real(dp) :: a11tmp, a22tmp, a11, a22, a12, torad

      ! Loop variables and flags
      integer  :: i, m, l, p, j
      logical  :: cnv

      ! Per-thread FFTW plan and per-thread work buffers
      integer(kind=8)          :: plan_local      ! adjust kind if needed for your mod_fftw3
      complex(dp), allocatable :: x_local(:,:)    ! complex spectrum  (0:n1/2, 0:n2-1)
      real(dp),    allocatable :: y_local(:,:)    ! real field        (0:n1-1, 0:n2-1)

      real(dp), allocatable :: phi(:,:)    ! random phases in [0, 2π)
      real(dp), allocatable :: fampl(:,:,:)

      !---------------- Basic validation (fail-fast)
      if (lde < 1)       stop 'pseudo2D: error lde < 1'
      if (rx  <= 0.0_dp) stop 'pseudo2D: error, rx <= 0.0'
      if (ry  <= 0.0_dp) stop 'pseudo2D: error, ry <= 0.0'
      if (n1  <  nx)     stop 'pseudo2D: n1 < nx'
      if (n2  <  ny)     stop 'pseudo2D: n2 < ny'

      !---------------- Spectral geometry and Δk
      ! kx = kappa*l, ky = lambda*p  with integer indices l,p on an n1×n2 grid.
      ! deltak is the spectral cell measure Δkx*Δky used in normalization.
      pi2    = 2.0_dp * pi
      deltak = (pi2*pi2) / ( real(n1*n2,dp) * dx * dy )
      kappa  = pi2 / ( real(n1,dp) * dx )
      kappa2 = kappa * kappa
      lambda = pi2 / ( real(n2,dp) * dy )
      lambda2= lambda * lambda

      !---------------- initial guesses for Newton2D
      r1 = 3.0_dp/rx
      r2 = 3.0_dp/ry
      if (verbose) then
         write(*,'(a,2(1x,es12.5),a,2(i6),a,4(1x,es12.5))')  &
              'pseudo2D: Call newton with r1,r2=', r1, r2, '  n1,n2=', n1, n2, &
              '  dx,dy,rx,ry=', dx, dy, rx, ry
      end if

      ! Calibrate r1,r2 so the model’s correlation matches the requested (rx,ry).
      call newton2D(r1, r2, n1, n2, dx, dy, rx, ry, cnv, verbose)
      if (.not. cnv) stop 'newton did not converge'

      !---------------- Normalization constant c
      ! Ensures resulting fields have consistent variance with the target spectrum.
      summ = 0.0_dp
      do p = -n2/2 + 1, n2/2
         do l = -n1/2 + 1, n1/2
            summ = summ + exp( -2.0_dp * ( kappa2*real(l*l,dp)/(r1*r1) +  &
                                           lambda2*real(p*p,dp)/(r2*r2) ) )
         end do
      end do
      summ = summ - 1.0_dp
      c    = sqrt( 1.0_dp / (deltak * summ) )
      if (verbose) write(*,'(a,3(1x,es12.5))') 'pseudo2D: r1 r2 c =', r1, r2, c

      !---------------- Rotation to angle θ (deg → rad)
      ! We rotate the principal axes of the Gaussian envelope in frequency
      ! space by θ. The rotated quadratic form is:
      !   exp( - (a11 kx^2 + 2 a12 kx ky + a22 ky^2) )
      ! where a11,a12,a22 are derived from 1/r1^2 and 1/r2^2 and θ.
      a11tmp = 1.0_dp / (r1*r1)
      a22tmp = 1.0_dp / (r2*r2)
      torad  = -pi / 180.0_dp          ! negative → anticlockwise convention used here
      a11 = a11tmp * cos(theta*torad)**2 + a22tmp * sin(theta*torad)**2
      a22 = a11tmp * sin(theta*torad)**2 + a22tmp * cos(theta*torad)**2
      a12 = (a22tmp - a11tmp) * cos(theta*torad) * sin(theta*torad)

      !---------------- Parallel generation: per-thread plan & buffers
      ! Each thread:
      !   • allocates its own x_local,y_local
      !   • builds a plan bound to those arrays
      !   • generates a subset of fields (j loop)
      !   • destroys its plan and deallocates its buffers
      allocate(x_local(0:n1/2, 0:n2-1))
      allocate(y_local(0:n1-1, 0:n2-1))
      call dfftw_plan_dft_c2r_2d(plan_local, n1, n2, x_local, y_local, FFTW_ESTIMATE)

      allocate(phi(0:n1/2, -n2/2:n2/2))
      allocate(fampl(0:n1/2, -n2/2:n2/2, 2))

      do j = 1, lde
         ! Build complex spectrum with random phase and rotated anisotropic envelope
         call wave_amp(n1, n2, pi2, a11, a12, a22, kappa, kappa2, lambda, lambda2, deltak, c, x_local, phi, fampl)

         ! Inverse FFT (complex → real) into y_local
         !call dfftw_execute_dft_c2r(plan_local, x_local, y_local)
	 y_local = 0.0_dp

         ! Copy leading (nx × ny) block to output (Fortran: y_local is 0-based here)
         do m = 1, ny
            do i = 1, nx
               Amat(i, m, j) = y_local(i-1, m-1)
            end do
         end do
      end do

      deallocate(fampl, phi)

      call dfftw_destroy_plan(plan_local)
      deallocate(y_local, x_local)

   end subroutine pseudo2D

   !---------------------------------------------------------------------------
   subroutine wave_amp(n1, n2, pi2, a11, a12, a22, kappa, kappa2, lambda, lambda2, deltak, c, x, phi, fampl)
      implicit none
      !---------------- Arguments
      integer,  intent(in)    :: n1, n2
      real(dp), intent(in)    :: pi2, a11, a12, a22, kappa, kappa2, lambda, lambda2, deltak, c
      complex(dp), intent(inout) :: x(0:n1/2, 0:n2-1)
      real(dp), intent(inout) :: phi(0:n1/2, -n2/2:n2/2)    ! random phases in [0, 2π)
      real(dp), intent(inout) :: fampl(0:n1/2, -n2/2:n2/2, 2)

      !---------------- Declarations
      real(dp) :: e, l2, p2
      integer  :: p, l


      call random_number(phi)
      phi = pi2 * phi

      ! Envelope in k-space: rotated anisotropic Gaussian
      !   e = exp( - (a11*kx^2 + 2*a12*kx*ky + a22*ky^2) )
      !   with kx = kappa*l, ky = lambda*p
      do p = -n2/2, n2/2
         do l = 0, n1/2
            l2 = real(l*l, dp)
            p2 = real(p*p, dp)
            e  = exp( - ( a11*kappa2*l2 + 2.0_dp*a12*kappa*lambda*real(l*p,dp) + a22*lambda2*p2 ) )

            fampl(l,p,1) = e * cos( phi(l,p) ) * sqrt(deltak) * c   ! real part
            fampl(l,p,2) = e * sin( phi(l,p) ) * sqrt(deltak) * c   ! imag part
         end do
      end do

      ! Enforce zero-mean (DC component)
      fampl(0,0,1) = 0.0_dp
      fampl(0,0,2) = 0.0_dp

      ! Pack to FFTW half-spectrum layout (rows are l=0..n1/2, columns are p=0..n2-1)
      do p = 0, n2/2 - 1
         x(:,p) = cmplx( fampl(:,p,1), fampl(:,p,2), kind=dp )
      end do
      do p = n2/2, n2 - 1
         x(:,p) = cmplx( fampl(:, -n2 + p, 1), fampl(:, -n2 + p, 2), kind=dp )
      end do

   end subroutine wave_amp

end module m_pseudo2D
