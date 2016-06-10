subroutine analysis6c(A, E, S, d, ndim, nrens, nrobs, verbose)
! Computes the analysed ensemble for A using the fast low rank square
! root formulation from Evensen 2004.  I.e. the sophisitcated inversion 
! algorithm is used.  This algorith does not use the full R but
! assumes (nrens-1) R = E E^T and E is specified.

   use m_multa
   use m_randrot
   implicit none
   integer, intent(in) :: ndim             ! dimension of model state
   integer, intent(in) :: nrens            ! number of ensemble members
   integer, intent(in) :: nrobs            ! number of observations
   
   real, intent(inout) :: A(ndim,nrens)    ! ensemble matrix
   real, intent(in)    :: S(nrobs,nrens)   ! matrix holding HA' 
   real, intent(in)    :: d(nrobs,1)         ! vector holding d-HA
   real, intent(in)    :: E(nrobs,nrens)   ! matrix holding observation perturbations
   logical, intent(in) :: verbose


   real ave(ndim)       ! ensemble mean

   real S0(nrobs,nrens)

   real VT(nrens,nrens)
   real, allocatable :: U0(:,:),sig0(:),VT0(:,:)
   real, allocatable :: U1(:,:),sig1(:),VT1(:,:)
   real, allocatable :: U2(:,:),sig2(:),VT2(:,:)
   integer iwork(nrobs)
   real    fwork(8*nrobs)
   real, allocatable :: X0(:,:),X1(:,:),X2(:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real, allocatable, dimension(:)   :: work,isigma

   real X3(nrens,nrens)
   real, allocatable :: y1(:)
   real, allocatable :: y2(:)
   real y3(nrobs)
   real y4(nrens) 

   real sigsum,sigsum1
   integer ierr,nrsigma,i,j,lwork, nrmin,m
   integer iblkmax
   character(len=2) tag2

   real X4(nrens,nrens)
   real IenN(nrens,nrens)
   real X5(nrens,nrens)

   nrmin=min(nrobs,nrens)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute SVD of S=HA'  ->  U0, sig0
   S0=S
   allocate (U0(nrobs,nrmin)  )
   allocate (sig0(nrmin))
   allocate (VT0(1,1))
   lwork=2*max(3*nrens+nrobs,5*nrens)
   allocate(work(lwork))
   sig0=0.0
   call dgesvd('S', 'N', nrobs, nrens, S0, nrobs, sig0, U0, nrobs, VT0, nrens, work, lwork, ierr)
   deallocate(work,VT0)
   if (ierr /= 0) then
      print *,'analysis6: ierr from call dgesvd 0= ',ierr; stop
   endif

   sigsum=sum( sig0(1:nrmin) )
   sigsum1=0.0
! Significant eigenvalues.
   nrsigma=0
   do i=1,nrmin                       
      if (sigsum1/sigsum < 0.9999) then
         nrsigma=nrsigma+1
         sigsum1=sigsum1+sig0(i)
      else
         sig0(i:nrmin)=0.0
         exit
      endif
   enddo

   if (verbose) then
      write(*,'(a,i5,g13.5)') ' dominant sing. values and share ',nrsigma,sigsum1/sigsum
      write(*,'(5g11.3)')sig0
   endif

   do i=1,nrsigma
       sig0(i) = 1.0/sig0(i)
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute X0=sig0^{*T} U0^T E 

! X0= U0^T R
   allocate(X0(nrmin,nrens))
   call dgemm('t','n',nrmin,nrens,nrobs, 1.0,U0,nrobs, E,nrobs, 0.0,X0,nrmin)


   do j=1,nrens
   do i=1,nrmin
      X0(i,j)=sig0(i)*X0(i,j)
   enddo
   enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute singular value decomposition  of X0(nrmin,nrens)
   allocate (U1(nrmin,nrmin)  )
   allocate (sig1(nrmin))
   allocate (VT1(1,1))
   lwork=2*max(3*nrens+nrobs,5*nrens)
   allocate(work(lwork))
   sig1=0.0

   call dgesvd('S', 'N', nrmin, nrens, X0, nrmin, sig1, U1, nrmin, VT1, 1, work, lwork, ierr)
   deallocate(work,VT1)
   if (ierr /= 0) then
      print *,'analysis6: ierr from call dgesvd 1= ',ierr; stop
   endif

   do i=1,nrmin
      sig1(i)=1.0/(1.0+sig1(i)**2)
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! X1 = U0 * sig0^{-1} * U1
   do j=1,nrmin
   do i=1,nrmin
      U1(i,j)=sig0(i)*U1(i,j)
   enddo
   enddo

   allocate(X1(nrobs,nrmin))
   call dgemm('n','n',nrobs,nrmin,nrmin, 1.0,U0,nrobs, U1,nrmin, 0.0,X1,nrobs)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update mean
   allocate(y1(1:nrmin))
   call dgemv('t',nrobs,nrmin,1.0,X1,nrobs,d,1,0.0,y1 ,1)
   allocate(y2(1:nrmin))
   y2=sig1*y1  
   call dgemv('n',nrobs,nrmin,1.0,X1,nrobs,y2,1,0.0,y3 ,1)
   call dgemv('t',nrobs,nrens,1.0,S ,nrobs,y3,1,0.0,y4 ,1)
   deallocate(y1,y2)

   do i=1,nrens
      X5(:,i)=y4(:)
   enddo

! X5=enN + (I - enN) X5  = enN + X5
   X5=1.0/float(nrens) + X5



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! X2 = (I+sig1^2)^{-0.5} * X1^T * S  
   allocate(X2(nrmin,nrens))
   call dgemm('t','n',nrmin,nrens,nrobs,1.0,X1,nrobs, S,nrobs, 0.0,X2,nrmin)


   do j=1,nrens
   do i=1,nrmin
      X2(i,j)=sqrt(sig1(i))*X2(i,j)
   enddo
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SVD of X2
   lwork=2*max(3*nrens+nrens,5*nrens)
   allocate (U2(nrmin,nrmin), sig2(nrmin), VT2(nrens,nrens), work(lwork))
   sig2=0.0
   call dgesvd('N', 'A', nrmin, nrens, X2, nrmin, sig2, U2, nrmin, VT2, nrens, work, lwork, ierr)
   deallocate(work)
   if (ierr /= 0) then
      print *,'ierr from call dgesvd 2 = ',ierr
      stop
   endif

   print *,'svd2'

   allocate(isigma(nrmin))
   isigma=1.0
   do i=1,nrmin
      if ( sig2(i) > 1.0 ) print *,'WARNING (analysis 5): sig2 > 1',i,sig2(i)
      isigma(i)=sqrt( max(1.0-sig2(i)**2,0.0) )
   enddo

   do j=1,nrens
      X3(:,j)=VT2(j,:)
   enddo


   do j=1,nrmin
      X3(:,j)=X3(:,j)*isigma(j)
   enddo

   print '(a)','A6: sig2: '
   print '(5g11.3)',sig2(1:nrmin)

   call randrot(VT,nrens)
   call dgemm('n','n',nrens,nrens,nrens,1.0,X3,nrens,VT,nrens,0.0,X4,nrens)
   
   IenN=-1.0/float(nrens)
   do i=1,nrens
      IenN(i,i)=  IenN(i,i) + 1.0
   enddo
   call dgemm('n','n',nrens,nrens,nrens,1.0,IenN,nrens,X4,nrens,1.0,X5,nrens)

! Final ensemble update
   iblkmax=min(ndim,200)
   call  multa(A, X5, ndim, nrens, iblkmax )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Analysis done
   tag2(1:2)='X5.uf'
   open(10,file='X5.uf',form='unformatted')
      write(10)tag2,nrens,X5
   close(10)

   open(10,file='X5col.dat')
      do j=1,nrens
         write(10,'(i5,f10.4)')j,sum(X5(:,j))
      enddo
   close(10)

   open(10,file='X5row.dat')
      do j=1,nrens
         write(10,'(i5,f10.4)')j,sum(X5(j,:))/float(nrens)
       enddo
   close(10)



   
end subroutine analysis6c
