module m_multa
   use iso_fortran_env, only : dp => real64
   implicit none
contains

subroutine multa(A, X, ndim, nrens, iblkmax)

   implicit none

   integer, intent(in) :: ndim
   integer, intent(in) :: nrens
   integer, intent(in) :: iblkmax

   real(dp), intent(in)    :: X(nrens,nrens)
   real(dp), intent(inout) :: A(ndim,nrens)

   real(dp), allocatable :: v(:,:)

   integer :: ia, ib, nrow

   allocate(v(iblkmax,nrens))

   do ia = 1, ndim, iblkmax

      ib   = min(ia+iblkmax-1,ndim)
      nrow = ib-ia+1

      v(1:nrow,1:nrens) = A(ia:ib,1:nrens)

      call dgemm('N','N', &
                 nrow, nrens, nrens, &
                 1.0_dp, &
                 v, iblkmax, &
                 X, nrens, &
                 0.0_dp, &
                 A(ia,1), ndim)

   enddo

   deallocate(v)

end subroutine multa

end module m_multa
