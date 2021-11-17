module m_sample2D
! This routine samples pseudo random fields with improved independence
! or orthogonality.  This is done by first drawing a large sample
! and construct the final sample using the dominant singular vectors
! of the large sample.

contains
subroutine sample2D(A2,nx,ny,nrens,nre,dx,dy,rx,ry,theta,samp_fix,verbose)
   use m_pseudo2D
   use m_randrot
   use m_fixsample2D
   implicit none
   integer, intent(in)     :: nx,ny
   integer, intent(in)     :: nrens
   integer, intent(in)     :: nre
   logical, intent(in)     :: samp_fix , verbose
   real,    intent(out)    :: A2(nx,ny,nrens)
   real,    intent(in)     :: dx,dy,rx,ry,theta

   integer n1,n2,n
   integer ns,msx,i,j,nsx
   integer lwork,ierr,iens
   integer ishape3(3)
   character(len=3) tag3
   real summ
   real, allocatable, dimension(:,:,:) :: A,UU
   real, allocatable, dimension(:,:)   :: U,VT,VT1,mean,var
   real, allocatable, dimension(:)     :: sig,work

   logical :: debug=.false.



   n1=nint(real(nx)*1.2)
   n2=nint(real(ny)*1.2)
   if (verbose) then
      print *,'nx=',nx
      print *,'ny=',ny
      print *,'n1=',n1
      print *,'n2=',n2
   end if

   do i=1,100
      if (2**i >= n1) then
         n1=2**i
         exit
      endif
   enddo
   do i=1,100
      if (2**i >= n2) then
         n2=2**i
         exit
      endif
   enddo

   if (verbose) then
      print *,'n1=',n1
      print *,'n2=',n2
   end if
   
   n=nx*ny
   ns=nre*nrens
   msx=min(ns,n)
   nsx=min(nrens,n)


   if (nre == 1) then
! Standard Monte Carlo sampling
      if (verbose) print *,'sample2d: calling pseudo2d'
      call pseudo2D(A2,nx,ny,nrens,rx,ry,dx,dy,n1,n2,theta,verbose)
      if (verbose) print *,'sample2d: pseudo2d done'

   elseif (nre > 1) then
! Start with oversized ensemble of ns members
      print *,'sample2d with nre=',nre
      lwork=2*max(3*ns+max(n,ns),5*ns)
      allocate(work(lwork))

      allocate(A(nx,ny,ns))
      print *,'sample2d: calling pseudo2d'
      call pseudo2D(A,nx,ny,ns,rx,ry,dx,dy,n1,n2,theta , verbose)
      print *,'sample2d: pseudo2d done'

      print *,'sample2d: calling randrot'
      allocate(VT1(nsx,nsx))
      call randrot(VT1,nsx)
      print *,'sample2d: randrot done'

! Compute SVD of oversized ensemble
      print *,'sample2d: svd of oversized ensemble'
      allocate( U(n,msx), sig(msx), VT(msx,msx) )
      call dgesvd('S', 'N', n, ns, A, n, sig, U, n, VT, msx, work, lwork, ierr)
      if (ierr /= 0) print *, 'ierr',ierr
      print *,'sample2d: svd done'

      write(tag3,'(i3.3)')ns


!      inquire(unit=11,opened=lopen)
!      open(11,file='sigma_'//tag3//'.dat')
!         summ=0.0
!         do i=1,ns
!            summ=summ+sig(i)**2
!            write(10,'(i4,3e12.4)')i,sig(i)/sig(1),sig(i)**2/sig(1)**2,summ/real(n*ns)
!         enddo
!      close(11)

! Generate first min(nrens,n) members
      print *,'sample2d: Generation of improved ensemble'
      ishape3=(/nx,ny,nsx/)
      allocate(UU(nx,ny,nsx))
      UU=reshape(U(:,1:nsx),ishape3)
      A2=0.0
      do j=1,nsx
         do i=1,nsx
            A2(:,:,j)=A2(:,:,j)+UU(:,:,i)*sig(i)/sqrt(real(nre))*VT1(i,j)
         enddo
      enddo
      deallocate(U, UU, VT, sig, VT1)
      print *,'sample2d: improved ensemble done'


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
   else
      print *,'invalid value for nre=',nre
      stop 'm_sample2D'
   endif



! subtract mean and correct variance
   
   if (samp_fix) call fixsample2D(A2,nx,ny,nrens)

   if (debug) then
      allocate(mean(nx,ny))
      allocate(var(nx,ny))
      mean=0.0
      do j=1,nrens
         mean(:,:)=mean(:,:)+A2(:,:,j)
      enddo
      mean=(1.0/real(nrens))*mean


      var=0.0
      do iens=1,nrens
         do j=1,ny
         do i=1,nx
            var(i,j)=var(i,j)+A2(i,j,iens)**2
         enddo
         enddo
      enddo

      open(10,file='check.dat')
         do j=1,ny
         do i=1,nx
            write(10,'(2i5,2g13.5)')i,j,mean(i,j),var(i,j)
         enddo
         enddo
      close(10)

      deallocate(mean)
      deallocate(var)

      stop
   endif

end subroutine sample2D
end module m_sample2D
