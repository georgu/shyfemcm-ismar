module m_random
   !---------------------------------------------------------------------------
   ! Double-precision random normal generator (standard normal N(0,1)).
   !
   ! Subroutine:
   !   random(work1, n)
   !     • work1(n): output vector of length n with N(0,1) samples
   !     • n       : length of the vector (no output if n <= 0)
   !
   ! Method:
   !   Vectorized Box–Muller transform:
   !     U1, U2 ~ i.i.d. Uniform(0,1)
   !     Z = sqrt(-2 ln U1) * cos(2π U2)   ~ N(0,1)
   !
   ! Notes:
   !   • Fully in double precision via iso_fortran_env (dp = real64).
   !   • Adds a guard: U1 = max(U1, tiny(1.0_dp)) to avoid log(0).
   !   • Keeps the original interface so existing callers need no changes.
   !---------------------------------------------------------------------------

   use iso_fortran_env, only : dp => real64
   implicit none
contains

   subroutine random(work1, n)
      implicit none
      !---------------- Arguments
      integer,  intent(in)  :: n
      real(dp), intent(out) :: work1(n)

      !---------------- Locals
      real(dp), allocatable :: work2(:)     ! second uniform for Box–Muller
      real(dp), parameter   :: pi     = 3.14159265358979323846264338327950288_dp
      real(dp), parameter   :: two_pi = 2.0_dp * pi

      !---------------- Quick return for empty size
      if (n <= 0) return

      ! Draw two independent uniform(0,1) vectors
      allocate(work2(n))
      call random_number(work1)
      call random_number(work2)

      ! Safeguard to prevent log(0): clamp U1 away from 0
      work1 = max(work1, tiny(1.0_dp))

      ! Vectorized Box–Muller: Z = sqrt(-2 ln U1) * cos(2π U2)
      work1 = sqrt(-2.0_dp * log(work1)) * cos(two_pi * work2)

      ! Cleanup
      deallocate(work2)
   end subroutine random

end module m_random
