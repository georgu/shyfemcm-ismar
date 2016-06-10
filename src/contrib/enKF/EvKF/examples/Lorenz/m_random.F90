module m_random
contains
subroutine random(work1)
!  Returns a vector of random values N(variance=1,mean=0)
   use mod_dimensions
   implicit none
   real, intent(out) :: work1(3*ndim)
   real, allocatable :: work2(:)

   allocate (work2(3*ndim))

   call random_number(work1)
   call random_number(work2)
   work1= sqrt(-2.0*log(work1))*cos(2.0*pi*work2)

   deallocate(work2)
end subroutine random
end module m_random
