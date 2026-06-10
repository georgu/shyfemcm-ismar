module m_randrot
   ! Generate a random real orthogonal matrix Q ∈ R^{nrens×nrens} (Haar measure)
   ! Method (Mezzadri, 2007):
   !   1) Draw A with i.i.d. N(0,1) entries (Box–Muller).
   !   2) QR factorization A = Q R via DGEQRF.
   !   3) Form explicit Q via DORGQR.
   !   4) Flip signs of columns of Q so diag(R) ≥ 0.
   use iso_fortran_env, only : dp => real64
   implicit none
contains

   subroutine randrot(Q, nrens)
      implicit none
      ! Arguments
      integer,  intent(in)  :: nrens
      real(dp), intent(out) :: Q(nrens, nrens)

      ! Locals
      real(dp), allocatable :: A(:,:), U1(:,:), U2(:,:)   ! Box–Muller uniforms
      real(dp), allocatable :: TAU(:)                     ! Householder scalars
      real(dp), allocatable :: WORK(:)                    ! LAPACK workspace
      real(dp), allocatable :: diagR(:)                   ! signs of diag(R)
      real(dp) :: pi, two_pi
      integer  :: lwork, info, i

      ! Constants
      pi     = 3.14159265358979323846264338327950288_dp
      two_pi = 2.0_dp * pi

      if (nrens <= 0) return

      ! Allocate
      allocate(A(nrens, nrens), U1(nrens, nrens), U2(nrens, nrens))
      allocate(TAU(nrens))
      allocate(diagR(nrens))

      ! Draw N(0,1) via vectorized Box–Muller
      call random_number(U1)
      call random_number(U2)
      U1 = max(U1, tiny(1.0_dp))                         ! avoid log(0)
      A  = sqrt(-2.0_dp * log(U1)) * cos(two_pi * U2)
      deallocate(U1, U2)

      ! -------------------- QR: A = Q*R (Q in A, R in upper triangle) --------------------
      ! Workspace query
      lwork = -1
      allocate(WORK(1))
      call dgeqrf(nrens, nrens, A, nrens, TAU, WORK, lwork, info)
      if (info /= 0) write(*,*) 'randrot: dgeqrf(query) info=', info
      lwork = max(1, int(WORK(1)))     ! optimal size returned in WORK(1)
      deallocate(WORK)
      allocate(WORK(lwork))

      ! Actual factorization
      call dgeqrf(nrens, nrens, A, nrens, TAU, WORK, lwork, info)
      if (info /= 0) write(*,*) 'randrot: dgeqrf info=', info
      deallocate(WORK)

      ! Capture sign of diag(R) from A(i,i)
      do i = 1, nrens
         if (A(i,i) > 0.0_dp) then
            diagR(i) =  1.0_dp
         else
            diagR(i) = -1.0_dp    ! also covers zero
         end if
      end do

      ! -------------------- Form explicit Q via DORGQR (overwrites A with Q) -------------
      lwork = -1
      allocate(WORK(1))
      call dorgqr(nrens, nrens, nrens, A, nrens, TAU, WORK, lwork, info)
      if (info /= 0) write(*,*) 'randrot: dorgqr(query) info=', info
      lwork = max(1, int(WORK(1)))
      deallocate(WORK)
      allocate(WORK(lwork))

      call dorgqr(nrens, nrens, nrens, A, nrens, TAU, WORK, lwork, info)
      if (info /= 0) write(*,*) 'randrot: dorgqr info=', info
      deallocate(WORK, TAU)

      ! Apply column sign corrections so diag(R) is nonnegative
      Q = A
      do i = 1, nrens
         if (diagR(i) < 0.0_dp) Q(:, i) = -Q(:, i)
      end do
      deallocate(A, diagR)
   end subroutine randrot

end module m_randrot
