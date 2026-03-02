module m_multa
   use iso_fortran_env, only: dp => real64
   implicit none
contains

   subroutine multa(A, X, ndim, nrens, iblkmax)
      use omp_lib
      implicit none

      integer, intent(in)    :: ndim, nrens, iblkmax
      real(dp), intent(in)   :: X(nrens, nrens)
      real(dp), intent(inout):: A(ndim, nrens)

      ! Work array becomes THREAD-PRIVATE
      real(dp) :: v(iblkmax, nrens)

      integer :: ia, ib, nrow

!$omp parallel default(shared) private(ia, ib, nrow, v)
!$omp do schedule(static)
      do ia = 1, ndim, iblkmax

         ib   = min(ia + iblkmax - 1, ndim)
         nrow = ib - ia + 1

         ! Copy local block
         v(1:nrow,1:nrens) = A(ia:ib,1:nrens)

         ! Multiply block: A(ia:ib,:) = v * X
         call dgemm('N','N',                &
                    nrow, nrens, nrens,     &
                    1.0_dp,                 &
                    v,     iblkmax,         &
                    X,     nrens,           &
                    0.0_dp,                 &
                    A(ia,1), ndim)

      end do
!$omp end do
!$omp end parallel

   end subroutine multa

end module m_multa
