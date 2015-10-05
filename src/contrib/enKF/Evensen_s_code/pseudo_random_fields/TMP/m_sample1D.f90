












module m_sample1D
! This routine samples pseudo random fields with improved coonditioning.
! This is done by first drawing a large sample and then construct the
! final sample using the dominant singular vectors of the large sample.

contains
subroutine sample1D(A2,n,nrens,nre,dx,rh,samp_fix,periodic)
   use m_pseudo1D
   use m_fixsample1D
   implicit none
   integer, intent(in)     ::  n
   integer, intent(in)     ::  nrens
   integer, intent(in)     ::  nre
   real,    intent(in)     ::  rh
   real,    intent(in)     ::  dx
   logical, intent(in)     ::  samp_fix
   logical, intent(in)     ::  periodic
   real,    intent(out)    ::  A2(n,nrens)
   logical, parameter :: new=.true.

   integer ns,msx,i,j,nsx,n1
   integer lwork,ierr
   real summ
   real, allocatable, dimension(:,:) :: A,U,VT,A0,VT1
   real, allocatable, dimension(:)   :: sig,work,mean,var

   logical :: debug=.false.

   if (periodic) then
      n1=n
   else
      n1=nint(real(n)*1.2)
   endif


   do i=1,100
      if (2**i >= n1) then
         n1=2**i
         exit
      endif
   enddo
   if (periodic .and. n /= n1) then
      print '(a,i5,a,i5)','Modified:  n=',n,' n1=',n1
      stop 'm_sample1D: You have to change model grid size for periodic samples'
   endif


   ns=nre*nrens
   msx=min(ns,n)
   nsx=min(nrens,n)

   if (nre == 1) then
! Standard Monte Carlo sampling
      call pseudo1D(A2,n,nrens,rh,dx,n1)

   elseif (nre > 1) then
! Start with oversized ensemble of ns members
      lwork=2*max(3*ns+max(n,ns),5*ns)
      allocate(work(lwork))

      allocate(A(n,ns))
      call pseudo1D(A,n,ns,rh,dx,n1)

! make an orthogonal VT1 used as linear combination for final ensemble
      allocate (A0(nsx,nsx), U(nsx,nsx), sig(nsx), VT1(nsx,nsx) )
      call pseudo1D(A0,nsx,nsx,rh,dx,nsx)
      call dgesvd('N', 'S', nsx, nsx, A0, nsx, sig, U, nsx, VT1, nsx, work, lwork, ierr)
      if (ierr /= 0) print *, 'ierr',ierr
      deallocate(A0, sig, U)


! Compute SVD of oversized ensemble
      allocate( U(n,msx), sig(msx), VT(msx,msx) )
      call dgesvd('S', 'N', n, ns, A, n, sig, U, n, VT, ns, work, lwork, ierr)
      if (ierr /= 0) print *, 'ierr',ierr

      open(10,file='sigma_.dat')
         summ=0.0
         do i=1,msx
            summ=summ+sig(i)**2
            write(10,'(i4,3e12.4)')i,sig(i)/sig(1),sig(i)**2/sig(1)**2,summ/real(n*ns)
         enddo
      close(10)

! Generate first nrens or nsx members
      A2=0.0
      do j=1,nsx
         do i=1,nsx
            A2(:,j)=A2(:,j)+U(:,i)*sig(i)/sqrt(real(nre))*VT1(i,j)
         enddo
      enddo
      deallocate(U, VT, sig, VT1)


      if (debug) then
! SVD of new ensemble
         allocate (U(n,nsx))
         allocate (sig(nsx))
         allocate (VT(nsx,nsx))
         sig=0.0
         call dgesvd('S', 'S', n, nsx, A2, n, sig, U, n, VT, nsx, work, lwork, ierr)
         if (ierr /= 0) print *, 'ierr',ierr

         open(10,file='sigma2.dat')
            summ=0.0
            do i=1,nsx
               summ=summ+sig(i)**2
               write(10,'(i4,3e12.4)')i,sig(i)/sig(1),sig(i)**2/sig(1)**2,summ/real(n*ns)
            enddo
         close(10)
         deallocate(U, VT, sig)
         stop
      endif
   endif



! subtract mean and correct variance
   
   if (samp_fix) call fixsample1D(A2,n,nrens)

   if (debug) then
      allocate(mean(n))
      allocate(var(n))
      mean=0.0
      do j=1,nrens
         mean(:)=mean(:)+A2(:,j)
      enddo
      mean=(1.0/real(nrens))*mean


      var=0.0
      do j=1,nrens
         do i=1,n
            var(i)=var(i)+ A2(i,j)**2
         enddo
      enddo

      open(10,file='check.dat')
         do i=1,n
            write(10,'(i5,2g13.5)')i,mean(i),var(i)
         enddo
      close(10)

      deallocate(mean)
      deallocate(var)

      stop
   endif

end subroutine sample1D
end module m_sample1D
