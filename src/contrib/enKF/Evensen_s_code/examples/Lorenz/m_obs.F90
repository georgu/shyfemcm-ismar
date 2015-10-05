module m_obs
contains
subroutine obs(sol,d,ih,v,dt,vars,deltaobs,inobs)
   use mod_dimensions
   use m_random2
   implicit none
   type(variances), intent(in) :: vars
   integer, intent(out) :: ih(nrmes)     ! Measurement matrix
   real, intent(out) :: v(nrmes,nrmes)   ! Measurement weigth matrix
   real, intent(out) :: d(3,nrmes)       ! data
   real, intent(in)  :: sol(3,ndim)      ! state variable
   real, intent(in)  :: dt
   integer, intent(in)  :: inobs         ! read observations if 1
   real, intent(out) :: deltaobs

   integer, parameter :: naux=33*nrmes
   real, allocatable :: aux(:)

   real det(2),rcond,t1,t2
   real h(1,1)
   integer i,ii

   allocate(aux(naux))

! Generating the measurement matrix
   h = 0.0

   deltaobs=float(ndim-1)*dt/float(nrmes)
   print *,'delta_obs=',deltaobs

   do i=1,nrmes
      ii=nint(deltaobs*float(i)/dt+1.0)
!H      h(i,ii)=1.0
      ih(i)=ii
   enddo
   if (ih(nrmes) > ndim) then
      print *,'OBS: ih(nrmes) >ndim',ih(nrmes),ndim
   endif

   call random2(d,3*nrmes)
   do i=1,nrmes
      d(:,i)=sqrt(vars%mes)*d(:,i)+sol(:,ih(i))
   enddo
!H   d=matmul(sol,transpose(H))

   if (inobs == 1) then
      open(10,file='mesA.dat')
         do i=1,nrmes
            read(10,*)t1,d(1:3,i),t2
         enddo
      close(10)
   endif

! Generating the measurement error covariance matrix
   v = 0.0
   do i = 1,nrmes
      v(i,i)=vars%mes
   enddo

   deallocate(aux)

end subroutine obs
end module m_obs
