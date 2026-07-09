module m_sample1D

   use iso_fortran_env, only : dp => real64
   use m_pseudo1D
   use m_fixsample1D

   implicit none

contains

subroutine sample1D(A2, n, nrens, nre, dx, rh, samp_fix, periodic)

   implicit none

   !==================== Arguments ====================
   integer,  intent(in) :: n
   integer,  intent(in) :: nrens
   integer,  intent(in) :: nre

   real(dp), intent(in) :: dx
   real(dp), intent(in) :: rh

   logical, intent(in) :: samp_fix
   logical, intent(in) :: periodic

   real(dp), intent(out) :: A2(n, nrens)

   !==================== Locals =======================
   integer :: ns
   integer :: msx
   integer :: nsx
   integer :: n1

   integer :: i, j
   integer :: ierr
   integer :: lwork

   real(dp) :: summ

   logical :: debug

   real(dp), allocatable :: A(:,:)
   real(dp), allocatable :: A0(:,:)

   real(dp), allocatable :: U(:,:)
   real(dp), allocatable :: VT1(:,:)

   real(dp), allocatable :: sig(:)
   real(dp), allocatable :: work(:)

   real(dp), allocatable :: mean(:)
   real(dp), allocatable :: var(:)

   real(dp) :: VTdummy(1,1)

   debug = .false.

   !===================================================
   ! Determine FFT size
   !===================================================
   if (periodic) then
      n1 = n
   else
      n1 = int(real(n,dp)*1.2_dp)
   endif

   do i = 1, 100
      if (2**i >= n1) then
         n1 = 2**i
         exit
      endif
   enddo

   if (periodic .and. n1 /= n) then
      write(*,'(a,i6,a,i6)') &
           'For periodic grid: n=', n, ' but n1=', n1
      stop 'm_sample1D: adjust model grid for periodic sampling'
   endif

   !===================================================
   ! Ensemble dimensions
   !===================================================
   ns  = nre * nrens
   msx = min(ns, n)
   nsx = min(nrens, n)

   !===================================================
   ! Standard Monte-Carlo sampling
   !===================================================
   if (nre == 1) then

      call pseudo1D(A2, n, nrens, rh, dx, n1)

      if (samp_fix) then
         call fixsample1D(A2, n, nrens)
      endif

      return

   endif

   !===================================================
   ! Oversampled ensemble
   !===================================================
   lwork = 2 * max(3*ns + max(n,ns), 5*ns)

   allocate(work(lwork))

   allocate(A(n,ns))

   call pseudo1D(A, n, ns, rh, dx, n1)

   !===================================================
   ! Generate orthogonal mixing matrix VT1
   !===================================================
   allocate(A0(nsx,nsx))
   allocate(U(nsx,nsx))
   allocate(sig(nsx))
   allocate(VT1(nsx,nsx))

   call pseudo1D(A0, nsx, nsx, rh, dx, nsx)

   call dgesvd('N','S', &
               nsx, nsx, &
               A0, nsx, &
               sig, &
               U, nsx, &
               VT1, nsx, &
               work, lwork, ierr)

   if (ierr /= 0) then
      write(*,*) 'sample1D: dgesvd(A0) ierr=', ierr
      stop
   endif

   deallocate(A0)
   deallocate(sig)
   deallocate(U)

   !===================================================
   ! SVD of oversized ensemble
   !===================================================
   allocate(U(n,msx))
   allocate(sig(msx))

   call dgesvd('S','N', &
               n, ns, &
               A, n, &
               sig, &
               U, n, &
               VTdummy, 1, &
               work, lwork, ierr)

   if (ierr /= 0) then
      write(*,*) 'sample1D: dgesvd(A) ierr=', ierr
      stop
   endif

   !===================================================
   ! Optional diagnostics
   !===================================================
   if (debug) then

      open(10,file='sigma_.dat')

      summ = 0.0_dp

      do i = 1, msx
         summ = summ + sig(i)**2

         write(10,'(i4,3e12.4)') &
              i, &
              sig(i)/sig(1), &
              sig(i)**2/sig(1)**2, &
              summ/real(n*ns,dp)
      enddo

      close(10)

   endif

   !===================================================
   ! Construct final ensemble
   !===================================================
   A2 = 0.0_dp

   do j = 1, nsx
      do i = 1, nsx

         A2(:,j) = A2(:,j) + &
                   U(:,i) * &
                   (sig(i)/sqrt(real(nre,dp))) * &
                   VT1(i,j)

      enddo
   enddo

   deallocate(U)
   deallocate(sig)
   deallocate(VT1)
   deallocate(A)

   !===================================================
   ! Mean removal / variance correction
   !===================================================
   if (samp_fix) then
      call fixsample1D(A2, n, nrens)
   endif

   !===================================================
   ! Debug check
   !===================================================
   if (debug) then

      allocate(mean(n))
      allocate(var(n))

      mean = 0.0_dp

      do j = 1, nrens
         mean(:) = mean(:) + A2(:,j)
      enddo

      mean = mean / real(nrens,dp)

      var = 0.0_dp

      do j = 1, nrens
         do i = 1, n
            var(i) = var(i) + A2(i,j)**2
         enddo
      enddo

      open(10,file='check.dat')

      do i = 1, n
         write(10,'(i5,2g13.5)') i, mean(i), var(i)
      enddo

      close(10)

      deallocate(mean)
      deallocate(var)

      stop 'DEBUG mode'

   endif

   if (allocated(work)) deallocate(work)

end subroutine sample1D

end module m_sample1D
