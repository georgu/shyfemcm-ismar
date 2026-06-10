module m_sample2D
!------------------------------------------------------------------------------
!  Purpose
!  --------
!  Sample 2-D random fields. If nre>1, draw an oversized ensemble (ns=nre*nrens),
!  compute its SVD, and build a reduced ensemble with improved independence/
!  orthogonality by mixing leading left singular vectors with a random
!  orthogonal matrix.
!
!  Precision: double everywhere (wp = real64).
!  Notes:
!    * We keep LAPACK DGESVD with an implicit interface (external dgesvd).
!    * OpenMP is used for embarrassingly-parallel loops (flatten/reconstruction).
!------------------------------------------------------------------------------

   use iso_fortran_env, only : wp => real64
   use m_pseudo2D
   use m_randrot
   use m_fixsample2D
   implicit none
   private
   public :: sample2D

contains

   subroutine sample2D(A2, nx, ny, nrens, nre, dx, dy, rx, ry, theta, samp_fix, verbose)
      implicit none

      ! Arguments
      integer,  intent(in)  :: nx, ny
      integer,  intent(in)  :: nrens
      integer,  intent(in)  :: nre
      logical,  intent(in)  :: samp_fix, verbose
      real(wp), intent(out) :: A2(nx,ny,nrens)
      real(wp), intent(in)  :: dx, dy, rx, ry, theta

      ! Sizes and counters
      integer :: n1, n2, n, ns, msx, nsx
      integer :: i, j, ierr, pow2, iens, lwork

      ! Misc
      real(wp) :: summ
      logical, parameter :: debug = .false.

      ! Workspace / temporaries
      real(wp), allocatable :: A(:,:,:)
      real(wp), allocatable :: Aflat(:,:), Eflat(:,:)
      real(wp), allocatable :: U(:,:), VT(:,:), VT1(:,:), mean(:,:), var(:,:)
      real(wp), allocatable :: sig(:), work(:)
      real(wp)              :: work_query(1)     ! LAPACK workspace query wants an array(1)

      ! Thread-private buffers (allocated inside parallel regions)
      ! They are declared allocatable here and made PRIVATE in the OMP clauses.
      real(wp), allocatable :: accum_vec(:)
      real(wp), allocatable :: local_mean(:,:), local_var(:,:)

      external :: dgesvd     ! implicit interface to LAPACK DGESVD

      !========================
      ! Power-of-two pads n1,n2 >= 1.2*nx,ny
      !========================
      n1 = int( real(nx,wp)*1.2_wp + 0.5_wp )
      n2 = int( real(ny,wp)*1.2_wp + 0.5_wp )

      do pow2 = 1, 100
         if (2**pow2 >= n1) then
            n1 = 2**pow2
            exit
         end if
      end do
      do pow2 = 1, 100
         if (2**pow2 >= n2) then
            n2 = 2**pow2
            exit
         end if
      end do

      if (verbose) print *, 'nx=',nx, ' ny=',ny, ' n1=',n1, ' n2=',n2

      n   = nx*ny
      ns  = nre*nrens
      msx = min(ns, n)
      nsx = min(nrens, n)

      if (nre == 1) then
         !---------------------------------------------------------------
         ! Standard Monte Carlo sampling
         !---------------------------------------------------------------
         if (verbose) print *,'sample2D: calling pseudo2D'
         call pseudo2D(A2, nx, ny, nrens, rx, ry, dx, dy, n1, n2, theta, verbose)
         if (verbose) print *,'sample2D: pseudo2D done'

      else if (nre > 1) then
         !---------------------------------------------------------------
         ! Improved sampling via SVD of an oversized ensemble
         !---------------------------------------------------------------
         if (verbose) print *, 'sample2D with nre=', nre

         allocate(A(nx,ny,ns))
         if (verbose) print *,'sample2D: calling pseudo2D (oversized)'
         call pseudo2D(A, nx, ny, ns, rx, ry, dx, dy, n1, n2, theta, verbose)
         if (verbose) print *,'sample2D: pseudo2D done'

         ! Random orthogonal nsx×nsx mixing matrix
         allocate(VT1(nsx,nsx))
         if (verbose) print *,'sample2D: calling randrot'
         call randrot(VT1, nsx)
         if (verbose) print *,'sample2D: randrot done'

         ! Flatten A(x,y,k) -> Aflat(n,k). Parallelize safely over k.
         allocate(Aflat(n,ns))
!$omp parallel do default(none) private(j) shared(Aflat,A,n,nx,ny,ns)
         do j = 1, ns
            Aflat(:,j) = reshape( A(:,:,j), (/ n /) )
         end do
!$omp end parallel do

         ! Thin SVD: Aflat = U * diag(sig) * VT (VT not needed with JOBVT='N')
         allocate(U(n,msx), sig(msx))
         allocate(VT(1,1))  ! dummy because JOBVT='N' ignores VT

         ! Workspace query (lwork = -1) with WORK as length-1 array
         lwork = -1
         call dgesvd('S','N', n, ns, Aflat, n, sig, U, n, VT, 1, work_query, lwork, ierr)
         if (ierr /= 0) error stop 'DGESVD(work query) failed in sample2D'
         lwork = max(1, int(work_query(1)))
         allocate(work(lwork))

         call dgesvd('S','N', n, ns, Aflat, n, sig, U, n, VT, 1, work, lwork, ierr)
         if (ierr /= 0) error stop 'DGESVD failed in sample2D'
         if (verbose) print *,'sample2D: SVD done'

         ! Build improved ensemble:
         ! A2(:,:,j) = sum_{i=1..nsx} reshape(U(:,i),nx,ny) * (sig(i)/sqrt(nre) * VT1(i,j))
         A2 = 0.0_wp

!$omp parallel default(none) &
!$omp shared(U,sig,VT1,A2,n,nx,ny,nsx,nre) private(j,i,accum_vec)
         allocate(accum_vec(n))
!$omp do schedule(static)
         do j = 1, nsx
            accum_vec = 0.0_wp
            do i = 1, nsx
               accum_vec = accum_vec + U(:,i) * ( sig(i) / sqrt(real(nre,wp)) * VT1(i,j) )
            end do
            A2(:,:,j) = reshape(accum_vec, (/ nx, ny /) )
         end do
!$omp end do
         deallocate(accum_vec)
!$omp end parallel

         ! Cleanup SVD temporaries and oversized sample
         deallocate(U, sig, VT1, Aflat, VT, work, A)

         ! Optional debug: singular spectrum of the improved ensemble
         if (debug) then
            allocate(Eflat(n,nsx))
!$omp parallel do default(none) private(j) shared(Eflat,A2,n,nx,ny,nsx)
            do j = 1, nsx
               Eflat(:,j) = reshape( A2(:,:,j), (/ n /) )
            end do
!$omp end parallel do

            allocate(U(n,nsx), sig(nsx), VT(nsx,nsx))
            lwork = -1
            call dgesvd('S','S', n, nsx, Eflat, n, sig, U, n, VT, nsx, work_query, lwork, ierr)
            if (ierr /= 0) error stop 'DGESVD(work query) failed (debug)'
            lwork = max(1, int(work_query(1)))
            allocate(work(lwork))
            call dgesvd('S','S', n, nsx, Eflat, n, sig, U, n, VT, nsx, work, lwork, ierr)
            if (ierr /= 0) error stop 'DGESVD failed (debug)'

            open(unit=10, file='sigma2.dat', status='replace', action='write')
            summ = 0.0_wp
            do i = 1, nsx
               summ = summ + sig(i)**2
               write(10,'(i6,3e16.8)') i, sig(i)/sig(1), (sig(i)/sig(1))**2, &
                                      summ/real(n*nsx,wp)
            end do
            close(10)
            deallocate(U, sig, VT, work, Eflat)
            error stop 'Debug stop'
         end if

      else
         error stop 'sample2D: invalid value for nre'
      end if

      ! Optional ensemble mean/variance fix
      if (samp_fix) call fixsample2D(A2, nx, ny, nrens)

      ! Optional diagnostics (OpenMP-safe reductions)
      if (debug) then
         allocate(mean(nx,ny), var(nx,ny))
         mean = 0.0_wp
         var  = 0.0_wp

!$omp parallel default(none) shared(A2,mean,var,nx,ny,nrens) private(local_mean,local_var,iens)
         allocate(local_mean(nx,ny), local_var(nx,ny))
         local_mean = 0.0_wp
         local_var  = 0.0_wp
!$omp do
         do iens = 1, nrens
            local_mean = local_mean + A2(:,:,iens)
            local_var  = local_var  + A2(:,:,iens)**2
         end do
!$omp end do
!$omp critical
         mean = mean + local_mean
         var  = var  + local_var
!$omp end critical
         deallocate(local_mean, local_var)
!$omp end parallel

         mean = mean / real(nrens,wp)

         open(unit=10, file='check.dat', status='replace', action='write')
         do j = 1, ny
            do i = 1, nx
               write(10,'(2i6,2e16.8)') i, j, mean(i,j), var(i,j)
            end do
         end do
         close(10)
         deallocate(mean, var)
         error stop 'Debug stop'
      end if

   end subroutine sample2D

end module m_sample2D
