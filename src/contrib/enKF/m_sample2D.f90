module m_sample2D

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

   !--------------------------------------------------
   ! Arguments
   !--------------------------------------------------
   integer,  intent(in)  :: nx
   integer,  intent(in)  :: ny
   integer,  intent(in)  :: nrens
   integer,  intent(in)  :: nre

   logical,  intent(in)  :: samp_fix
   logical,  intent(in)  :: verbose

   real(wp), intent(out) :: A2(nx,ny,nrens)

   real(wp), intent(in)  :: dx
   real(wp), intent(in)  :: dy
   real(wp), intent(in)  :: rx
   real(wp), intent(in)  :: ry
   real(wp), intent(in)  :: theta

   !--------------------------------------------------
   ! Local variables
   !--------------------------------------------------
   integer :: n1
   integer :: n2
   integer :: n

   integer :: ns
   integer :: msx
   integer :: nsx

   integer :: i
   integer :: j
   integer :: ierr
   integer :: iens
   integer :: lwork
   integer :: pow2

   real(wp) :: summ

   logical, parameter :: debug = .false.

   !--------------------------------------------------
   ! Arrays
   !--------------------------------------------------
   real(wp), allocatable :: A(:,:,:)
   real(wp), allocatable :: UU(:,:,:)

   real(wp), allocatable :: Aflat(:,:)
   real(wp), allocatable :: Eflat(:,:)

   real(wp), allocatable :: U(:,:)
   real(wp), allocatable :: VT(:,:)
   real(wp), allocatable :: VT1(:,:)

   real(wp), allocatable :: mean(:,:)
   real(wp), allocatable :: var(:,:)

   real(wp), allocatable :: sig(:)
   real(wp), allocatable :: work(:)

   real(wp) :: work_query(1)

   ! Dummy matrix because JOBVT='N'
   real(wp) :: VTdummy(1,1)

   external :: dgesvd

   !--------------------------------------------------
   ! Compute FFT dimensions
   !--------------------------------------------------
   n1 = int(real(nx,wp)*1.2_wp + 0.5_wp)
   n2 = int(real(ny,wp)*1.2_wp + 0.5_wp)

   do pow2 = 1,100
      if (2**pow2 >= n1) then
         n1 = 2**pow2
         exit
      endif
   enddo

   do pow2 = 1,100
      if (2**pow2 >= n2) then
         n2 = 2**pow2
         exit
      endif
   enddo

   if (verbose) then
      print *,'nx=',nx
      print *,'ny=',ny
      print *,'n1=',n1
      print *,'n2=',n2
   endif

   n   = nx*ny
   ns  = nre*nrens
   msx = min(ns,n)
   nsx = min(nrens,n)

   !---------------------------------------------------------------
   ! Standard Monte Carlo sampling
   !---------------------------------------------------------------
   if (nre == 1) then

      if (verbose) print *,'sample2D: calling pseudo2D'

      call pseudo2D(A2, nx, ny, nrens, &
                    rx, ry, dx, dy, &
                    n1, n2, theta, verbose)

      if (verbose) print *,'sample2D: pseudo2D done'

   !---------------------------------------------------------------
   ! Improved sampling via SVD of oversized ensemble
   !---------------------------------------------------------------
   else if (nre > 1) then

      if (verbose) print *,'sample2D with nre=',nre

      allocate(A(nx,ny,ns))

      if (verbose) print *,'sample2D: calling pseudo2D'

      call pseudo2D(A, nx, ny, ns, &
                    rx, ry, dx, dy, &
                    n1, n2, theta, verbose)

      if (verbose) print *,'sample2D: pseudo2D done'

      allocate(VT1(nsx,nsx))

      if (verbose) print *,'sample2D: calling randrot'

      call randrot(VT1, nsx)

      if (verbose) print *,'sample2D: randrot done'

      !------------------------------------------------------------
      ! Flatten A(nx,ny,ns) --> Aflat(n,ns)
      !------------------------------------------------------------
      allocate(Aflat(n,ns))

      do j = 1, ns
         Aflat(:,j) = reshape(A(:,:,j), (/n/) )
      end do

      !------------------------------------------------------------
      ! SVD of oversized ensemble
      !------------------------------------------------------------
      allocate(U(n,msx))
      allocate(sig(msx))

      VTdummy = 0.0_wp

      lwork = -1

      call dgesvd('S','N', &
                  n, ns, &
                  Aflat, n, &
                  sig, &
                  U, n, &
                  VTdummy, 1, &
                  work_query, lwork, ierr)

      if (ierr /= 0) stop 'sample2D: DGESVD workspace query failed'

      lwork = max(1, int(work_query(1)))

      allocate(work(lwork))

      call dgesvd('S','N', &
                  n, ns, &
                  Aflat, n, &
                  sig, &
                  U, n, &
                  VTdummy, 1, &
                  work, lwork, ierr)

      if (ierr /= 0) stop 'sample2D: DGESVD failed'

      if (verbose) print *,'sample2D: svd done'

      !------------------------------------------------------------
      ! Generate improved ensemble (same algorithm as original)
      !------------------------------------------------------------
      allocate(UU(nx,ny,nsx))

      UU = reshape(U(:,1:nsx), (/nx,ny,nsx/) )

      A2 = 0.0_wp

      do j = 1, nsx
         do i = 1, nsx

            A2(:,:,j) = A2(:,:,j) + &
                         UU(:,:,i) * &
                         sig(i) / sqrt(real(nre,wp)) * &
                         VT1(i,j)

         end do
      end do

      if (verbose) print *,'sample2D: improved ensemble done'

      deallocate(UU)
      deallocate(U)
      deallocate(sig)
      deallocate(VT1)

      deallocate(Aflat)

      if (allocated(work)) deallocate(work)

      deallocate(A)

   else

      stop 'sample2D: invalid value of nre'

   endif

      ! Optional ensemble mean/variance fix
      if (samp_fix) call fixsample2D(A2, nx, ny, nrens)

      ! Optional diagnostics
      if (debug) then

         allocate(mean(nx,ny))
         allocate(var(nx,ny))

         mean = 0.0_wp
         var  = 0.0_wp

         do iens = 1, nrens
            mean = mean + A2(:,:,iens)
         end do

         mean = mean / real(nrens,wp)

         do iens = 1, nrens
            var = var + A2(:,:,iens)**2
         end do

         open(unit=10, file='check.dat', status='replace', action='write')

         do j = 1, ny
            do i = 1, nx
               write(10,'(2i6,2e16.8)') &
                    i, j, mean(i,j), var(i,j)
            end do
         end do

         close(10)

         deallocate(mean)
         deallocate(var)

         stop 'DEBUG mode'

      end if

   end subroutine sample2D

end module m_sample2D
