module m_eksscalc
contains
subroutine eksscalc(p,A,S,bb,zdat,idat,vars)
   use mod_dimensions
   use m_ranmean2
   use m_random2
   implicit none
   real, intent(in)    :: p(3)
   real, intent(inout) :: A(3,nrsamp)
   real, intent(inout) :: S(nrtobs,nrsamp)
   real, intent(inout) :: bb(nrobs,nrmes,0:nrsamp)
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
   print *,'idat',idat
   print *,'Rows in S',(idat-1)*nrobs+1,idat*nrobs
   S((idat-1)*nrobs+1:idat*nrobs,1:nrsamp)=matmul(H,A)
   

! Add back the emsemble mean 
   do i=1,nrsamp
      A(:,i)=A(:,i)+w(:)
   enddo

! Store data misfit 
   bb(1:nrobs,idat,0)=zdat(1:nrobs,idat)-matmul(H,p)
   do i=1,nrsamp
      call random2(w,3)
      bb(1:nrobs,idat,i)=zdat(1:nrobs,idat)+sqrt(vars%mes)*w-matmul(H,A(1:3,i))
   enddo

   deallocate (H, w)

end subroutine eksscalc
end module
