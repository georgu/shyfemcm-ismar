module m_eksscalc_lag
contains
subroutine eksscalc_lag(p,A,S,bb,zdat,idat,vars)
   use mod_dimensions
   use m_random2
   use m_ranmean2
   implicit none
   real, intent(in)    :: p(3)
   real, intent(inout) :: A(3,nrsamp)
   real, intent(inout) :: S(nrobs,nrsamp)
   real, intent(inout) :: bb(nrobs,0:nrsamp)
   real, intent(inout) :: zdat(nrobs,nrmes)
   integer, intent(in) :: idat
   type(variances) vars

   integer i,j,m

   real, allocatable :: H(:,:)
   real, allocatable :: w(:)
   allocate (H(nrobs,3))
   allocate (w(3))

! Subtract the emsemble mean
   call ranmean2(A,w)
   do i=1,nrsamp
      A(:,i)=A(:,i)-w(:)
   enddo

   H=0.0
   do m=1,3
      H(m,m)=1.0
   enddo

   S(1:nrobs,1:nrsamp)=matmul(H,A)
   

! Add back the emsemble mean 
   do i=1,nrsamp
      A(:,i)=A(:,i)+w(:)
   enddo

! Store data misfit 
   bb(1:nrobs,0)=zdat(1:nrobs,idat)-matmul(H,p)
   do i=1,nrsamp
      call random2(w,3)
      bb(1:nrobs,i)=zdat(1:nrobs,idat)+sqrt(vars%mes)*w-matmul(H,A(1:3,i))
   enddo

   deallocate (H, w)

end subroutine eksscalc_lag
end module
