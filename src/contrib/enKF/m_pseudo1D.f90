module m_pseudo1D
   !
   !  • Arrays are clearly dimensioned according to FFTW requirements.
   !  • The FFTW execution call is explicit (dft_c2r), not generic.
   !  • Added full English documentation for each computational step.
   !  • Added safety checks and clarified assumptions (e.g., n1 must be even).
   !
   use iso_fortran_env, only : dp => real64
   use m_random
   use m_newton1D
   use mod_fftw3
   implicit none
contains
!======================================================================
subroutine pseudo1D(A, nx, nrfields, rx, dx, n1)
!======================================================================
   !
   ! Generates nrfields 1D pseudo-random fields with a specified
   ! correlation length rx using a Gaussian spectrum and FFTW.
   !
   integer, intent(in) :: nx        ! physical vector length
   integer, intent(in) :: nrfields  ! number of fields to produce
   real(dp), intent(out) :: A(nx,nrfields)
   real(dp), intent(in) :: rx       ! correlation length
   real(dp), intent(in) :: dx       ! grid spacing (>0)
   integer, intent(in) :: n1        ! FFT length (MUST be even, n1 >= nx)

   ! Internal variables
   real(dp) :: r1, r12, c
   real(dp) :: kappa, kappa2
   real(dp) :: pi2, deltak, summ
   real(dp) :: tt
   integer :: l, i
   logical :: cnv
   logical, save :: diag = .false.
   real(dp), parameter :: pi = 3.141592653589793238_dp

   ! FFTW plan and buffers
   integer(kind=8) :: plan
   complex(dp), allocatable :: arrayC(:)  ! complex half-spectrum
   real(dp),    allocatable :: y(:)       ! real output field

   ! Spectrum arrays
   real(dp), allocatable :: fampl(:,:)    ! (0:n1/2, 2)
   real(dp), allocatable :: phi(:)        ! random phases

   !--------------------------- Safety checks
   if (n1 < nx) stop 'pseudo1D: n1 must be >= nx'
   if (mod(n1,2) /= 0) stop 'pseudo1D: n1 must be even (FFT requirement)'
   if (rx <= 0.0_dp) stop 'pseudo1D: rx must be > 0'
   if (dx <= 0.0_dp) stop 'pseudo1D: dx must be > 0'

   !--------------------------- FFTW-compatible dimensions
   allocate(arrayC(0:n1/2))     ! complex, n1/2+1 elements
   allocate(y(0:n1-1))          ! full-size real field

   !--------------------------- Spectral constants
   pi2    = 2.0_dp*pi
   deltak = pi2 / ( real(n1,dp)*dx )
   kappa  = pi2 / ( real(n1,dp)*dx )
   kappa2 = kappa*kappa

   ! Plan creation (complex-to-real)
   call dfftw_plan_dft_c2r_1d(plan, n1, arrayC, y, FFTW_ESTIMATE)

   !--------------------------- Initial guess for decorrelation radius
   r1  = 3.0_dp/rx
   if (diag) write(*,*) 'Newton initial guess r1=', r1

   ! Newton solver for matching target correlation length
   call newton1D(r1, n1, dx, rx, cnv)
   if (.not.cnv) stop 'pseudo1D: Newton did not converge'

   r12 = r1*r1

   !--------------------------- Normalisation constant (kept identical)
   summ = 0.0_dp
   do l = -n1/2+1, n1/2
      summ = summ + exp( -2.0_dp*(kappa2*l*l)/r12 )
   end do
   summ = summ - 1.0_dp
   c = sqrt( 1.0_dp / (deltak*summ) )

   !--------------------------- Allocate per-field spectral arrays
   allocate(phi(0:n1/2))
   allocate(fampl(0:n1/2, 2))

   !==================================================================
   ! Main generation loop (no OpenMP here because FFTW plan is shared
   ! and we keep exact behaviour of the old version)
   !==================================================================
   do i = 1, nrfields

      ! Draw random phases uniformly in [0, 2π)
      call random_number(phi)
      phi = pi2 * phi

      ! Build spectral Gaussian envelope
      do l = 1, n1/2
	 tt = kappa2 * real(l*l,dp) / r12
         fampl(l,1) = exp(-tt) * cos(phi(l)) * sqrt(deltak) * c
         fampl(l,2) = exp(-tt) * sin(phi(l)) * sqrt(deltak) * c
      end do

      ! Zero DC component: enforces zero-mean fields
      fampl(0,1) = 0.0_dp
      fampl(0,2) = 0.0_dp

      ! Pack into FFTW half-spectrum
      arrayC = cmplx(fampl(:,1), fampl(:,2), kind=dp)

      ! Inverse FFT
      call dfftw_execute_dft_c2r(plan, arrayC, y)

      ! Extract physical subvector
      A(:,i) = y(0:nx-1)

   end do

   !--------------------------- Cleanup
   call dfftw_destroy_plan(plan)
   deallocate(arrayC, y, phi, fampl)

end subroutine pseudo1D
!======================================================================

end module m_pseudo1D
