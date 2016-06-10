module m_random3
contains
subroutine random3(work1,n)
!  Returns a vector of random values N(variance=1,mean=0)
   use mod_dimensions
   implicit none
   integer, intent(in) :: n
   real,   intent(out) :: work1(n)
   integer i,j
   real t1,t2

   do i=1,n
      do j=1,100
         call random_number(t1)
         call random_number(t2)
         work1(i)= sqrt(-2.0*log(t1))*cos(2.0*pi*t2)
         if (abs(work1(i)) < 2.0) exit
      enddo
   enddo
end subroutine random3
end module m_random3
