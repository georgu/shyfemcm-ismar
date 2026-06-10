module m_sample1D
   !! ------------------------------------------------------------------------
   !! 1D ensemble sampler with “improved conditioning” via oversampling + SVD.
   !!
   !! Purpose:
   !!   Generate nrens 1D random fields of length n with spatial correlation
   !!   defined by rh, dx — using m_pseudo1D (double‑precision version).
   !!
   !!   If nre > 1:
   !!      1) Generate ns = nre*nrens independent pseudo random fields.
   !!      2) Compute their SVD:  A = U * diag(sig) * VT
   !!      3) Keep only the nsx = min(n,nrens) leading singular vectors.
   !!      4) Generate an orthogonal mixing matrix VT1 via a second random sample,
   !!         project the dominant SVD modes into nsx final samples.
   !!
   !!   This produces an ensemble with much better conditioning, reduced
   !!   spurious correlations, and more isotropic sampling of the dominant modes.
   !!
   !! Inputs:
   !!   n       – physical grid size
   !!   nrens   – final ensemble size
   !!   nre     – oversampling factor (nre=1 → plain Monte Carlo)
   !!   dx,rh   – correlation parameters forwarded to pseudo1D
   !!   samp_fix – whether to impose zero mean and unit variance (fixsample1D)
   !!   periodic – periodic grid option
   !!
   !! Output:
   !!   A2(n,nrens) – final ensemble
   !!
   !! Dependencies:
   !!   use m_pseudo1D   – DP pseudo random field generator
   !!   use m_fixsample1D – optional diagnostic correction
   !!
   !! ------------------------------------------------------------------------
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
   real(dp), intent(in) :: dx, rh
   logical, intent(in) :: samp_fix
   logical, intent(in) :: periodic
   real(dp), intent(out) :: A2(n, nrens)

   !==================== Locals ========================
   integer :: ns, msx, nsx, n1
   integer :: i, j, ierr, lwork
   real(dp) :: summ
   logical :: debug
   real(dp), allocatable :: A(:,:), A0(:,:), U(:,:), VT(:,:), VT1(:,:)
   real(dp), allocatable :: sig(:), work(:)
   real(dp), allocatable :: mean(:), var(:)

   debug = .false.

   !--------------------------------------------------------------------
   ! Determine FFT size n1
   ! If periodic → n1 must equal n
   ! Otherwise choose smallest power of 2 ≥ 1.2*n
   !--------------------------------------------------------------------
   if (periodic) then
      n1 = n
   else
      n1 = int(real(n,dp)*1.2_dp)
   end if

   do i = 1, 100
      if (2**i >= n1) then
         n1 = 2**i
         exit
      end if
   end do

   if (periodic .and. n1 /= n) then
      write(*,'(a,i6,a,i6)') 'For periodic grid: n=',n,' but n1=',n1
      stop 'm_sample1D: adjust model grid for periodic sampling'
   end if

   !--------------------------------------------------------------------
   ! Ensemble dimensions
   ! ns = full oversampled ensemble size
   ! nsx = min(n, nrens) — rank truncation size
   ! msx = min(ns, n) — number of usable singular values
   !--------------------------------------------------------------------
   ns  = nre * nrens
   nsx = min(n, nrens)
   msx = min(ns, n)

   !--------------------------------------------------------------------
   ! Case 1: plain Monte Carlo (no SVD conditioning)
   !--------------------------------------------------------------------
   if (nre == 1) then

      call pseudo1D(A2, n, nrens, rh, dx, n1)
      if (samp_fix) call fixsample1D(A2, n, nrens)
      return

   end if

   !--------------------------------------------------------------------
   ! Case 2: oversampling (nre > 1)
   !
   ! Step 1 – Generate oversized ensemble A(n × ns)
   !--------------------------------------------------------------------
   lwork = 2 * max( 3*ns + max(n,ns), 5*ns )
   allocate(work(lwork))

   allocate(A(n, ns))
   call pseudo1D(A, n, ns, rh, dx, n1)

   !--------------------------------------------------------------------
   ! Step 2 – Generate orthogonal mixing matrix VT1(nsx × nsx)
   !
   ! We generate a random square pseudo1D sample A0, then take its SVD.
   ! The left singular vectors (U) or right singular vectors (VT1)
   ! produce an orthonormal basis used to mix dominant modes.
   !--------------------------------------------------------------------
   allocate( A0(nsx, nsx), U(nsx, nsx), sig(nsx), VT1(nsx, nsx) )

   call pseudo1D(A0, nsx, nsx, rh, dx, nsx)

   ! SVD of A0: keep singular vectors in VT1
   call dgesvd('N', 'S', nsx, nsx, A0, nsx, sig, U, nsx, VT1, nsx, work, lwork, ierr)
   if (ierr /= 0) write(*,*) 'sample1D: dgesvd(A0) ierr=', ierr

   deallocate(A0, sig, U)

   !--------------------------------------------------------------------
   ! Step 3 – SVD of oversized ensemble A(n × ns)
   !          Keep only top msx = min(n,ns) modes for conditioning
   !--------------------------------------------------------------------
   allocate(U(n, msx), sig(msx), VT(msx, msx))

   call dgesvd('S', 'N', n, ns, A, n, sig, U, n, VT, msx, work, lwork, ierr)
   if (ierr /= 0) write(*,*) 'sample1D: dgesvd(A) ierr=', ierr

   !--------------------------------------------------------------------
   ! Optional diagnostic output of spectra
   !--------------------------------------------------------------------
   open(10,file='sigma_.dat')
   summ = 0.0_dp
   do i = 1, msx
      summ = summ + sig(i)**2
      write(10,'(i4,3e12.4)') i, sig(i)/sig(1), sig(i)**2/sig(1)**2, summ/real(n*ns,dp)
   end do
   close(10)

   !--------------------------------------------------------------------
   ! Step 4 – Construct final ensemble
   !
   ! A2(:,j) = Σ_i  U(:,i) * sig(i)/√nre * VT1(i,j)
   !
   ! This compresses the ensemble to the spanning subspace defined
   ! by the leading singular vectors, then mixes via VT1.
   !--------------------------------------------------------------------
   A2 = 0.0_dp

   do j = 1, nsx
      do i = 1, nsx
         A2(:,j) = A2(:,j) + U(:,i) * ( sig(i) / sqrt(real(nre,dp)) ) * VT1(i,j)
      end do
   end do

   deallocate(U, VT, sig, VT1)
   deallocate(A)

   !--------------------------------------------------------------------
   ! Step 5 – Fix mean and variance (optional)
   !--------------------------------------------------------------------
   if (samp_fix) call fixsample1D(A2, n, nrens)

   !--------------------------------------------------------------------
   ! Optional debugging: check mean & variance
   !--------------------------------------------------------------------
   if (debug) then
      allocate(mean(n), var(n))
      mean = 0.0_dp
      do j = 1, nrens
         mean(:) = mean(:) + A2(:,j)
      end do
      mean = mean / real(nrens,dp)

      var = 0.0_dp
      do j = 1, nrens
         do i = 1, n
            var(i) = var(i) + A2(i,j)**2
         end do
      end do

      open(10,file='check.dat')
      do i = 1, n
         write(10,'(i5,2g13.5)') i, mean(i), var(i)
      end do
      close(10)

      stop 'DEBUG mode'
   end if

end subroutine sample1D

end module m_sample1D
