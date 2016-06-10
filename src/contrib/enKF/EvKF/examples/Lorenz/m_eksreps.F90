module m_eksreps
contains
subroutine eksreps(S,A,RR,Reps)
   use mod_dimensions
   use m_ranmean2
   implicit none
   real, intent(inout):: A(3,nrsamp)
   real, intent(in)   :: S(nrtobs,nrsamp)
   real, intent(in)   :: RR(nrtobs,nrtobs)
   real, intent(out)  :: Reps(nrtobs,3)

   integer i,j,m,ii

   real, allocatable :: w(:)

   allocate (w(3))

   call ranmean2(A,w)

   do i=1,nrsamp
   A(:,i)=A(:,i)-w(:)
   enddo

   Reps=matmul(S,transpose(A))

   Reps=Reps/float(nrsamp)

   do i=1,nrsamp
   A(:,i)=A(:,i)+w(:)
   enddo

   Reps=matmul(transpose(RR),Reps)

   deallocate (w)

end subroutine eksreps
end module m_eksreps
