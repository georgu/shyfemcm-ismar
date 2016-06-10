subroutine analysis4c(A, R, S, d, ndim, nrens, nrobs, verbose)
! Computes the analysed ensemble for A using the square root formulation
! Algorithm from Section 3 in Evensen 2004, but combining the 
! update of the mean and the peterurbations into one equation.

   use m_multa
   use m_randrot
   implicit none
   integer, intent(in) :: ndim             ! dimension of model state
   integer, intent(in) :: nrens            ! number of ensemble members
   integer, intent(in) :: nrobs            ! number of observations
   
   real, intent(inout) :: A(ndim,nrens)    ! ensemble matrix
   real, intent(in)    :: R(nrobs,nrobs)   ! matrix holding R
   real, intent(in)    :: S(nrobs,nrens)   ! matrix holding HA' 
   real, intent(in)    :: d(nrobs)         ! vector holding d-HA
   logical, intent(in) :: verbose


   real ave(ndim)                          ! ensemble mean

! eigenvalue decomposition using RSM
   real eig(nrobs),Z(nrobs,nrobs),fwork(8*nrobs)  
   integer iwork(nrobs)

   real ZZ(nrobs,nrobs),y22(nrobs)
   

   real X1(nrobs,nrobs),X2(nrobs,nrens),X3(nrens,nrens),X4(nrens,nrens)
   real y1(nrobs),y2(nrobs),y3(nrobs),y4(nrens) 

   real, allocatable :: U(:,:),V(:,:),VT(:,:),sig(:),work(:)
   real, allocatable, dimension(:)   :: isigma

   real, allocatable :: ap(:)

   real sigsum,sigsum1,oneobs(1,1)
   integer ierr,nrsigma,i,j,lwork, nrmin,m,k
   integer iblkmax
   character(len=2) tag2

   real IenN(nrens,nrens)
   real X5(nrens,nrens)


! Evaluate R= S*S' + (nrens)*R
   call dgemm('n','t',nrobs,nrobs,nrens, &
                 1.0, S, nrobs, &
                      S, nrobs, &
        float(nrens), R, nrobs)



! Compute eigenvalue decomposition of R -> Z*eig*Z'
#ifdef IBM
! Upper packed storage as in ESSL manual
   allocate (ap(nrobs*(nrobs+1)/2) )
   k=0
   do j=1,nrobs
   do i=1,j
       k=k+1
       ap(k)=R(i,j)
   enddo
   enddo
   call dspev(21,ap,eig,Z,nrobs,nrobs,fwork,2*nrobs)
   deallocate(ap)
#else
   call rsm(nrobs, nrobs, R, eig, nrobs, Z, fwork, iwork, ierr)
   if (ierr /= 0) then
      print *,'analysis5: ierr from call rsm 1= ',ierr; stop
   endif
#endif



! Significant eigenvalues
   sigsum=sum( eig(1:nrobs) )
   sigsum1=0.0
   nrsigma=0
   do i=nrobs,1,-1
      print *,'Eigen values: ',i,eig(i)
      if (sigsum1/sigsum < 0.999) then
         nrsigma=nrsigma+1
         sigsum1=sigsum1+eig(i)
         eig(i) = 1.0/eig(i)
      else
         eig(1:i)=0.0
         exit
      endif
   enddo
   if (verbose) then
      write(*,'(a,i5,g13.5)') ' dominant eigenvalues and share ',nrsigma,sigsum1/sigsum
   endif

   open(10,file='ana4Ceig.dat')
      do i=1,nrobs
         write(10,'(i3,f13.5)')i,eig(i)
      enddo
   close(10)
   open(10,file='ana4Cvec.dat')
      write(10,'(6f13.5)')Z
   close(10)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update mean
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call dgemv('t',nrobs,nrobs,1.0,Z,nrobs,d,1,0.0,y1 ,1)
   y2=eig*y1
   call dgemv('n',nrobs,nrobs,1.0,Z ,nrobs,y2,1,0.0,y3 ,1)
   call dgemv('t',nrobs,nrens,1.0,S ,nrobs,y3,1,0.0,y4 ,1)

   do i=1,nrens
      X5(:,i)=y4(:)
   enddo

! X5=enN + (I - enN) X5  = enN + X5
   X5=1.0/float(nrens) + X5


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update perturbations
! Creates  X2= eig^-0.5 U^T S
   do i=1,nrobs
      eig(i)=sqrt(eig(i))
   enddo
   print '(a)','A4: eig: '
   print '(5g11.3)',eig(1:nrobs)
   call dgemm('t','n',nrobs,nrens,nrobs,1.0,Z,nrobs, S,nrobs, 0.0,X2,nrobs)
   do j=1,nrens
   do i=1,nrobs
      X2(i,j)=eig(i)*X2(i,j)
   enddo
   enddo


! SVD of X2
   nrmin=min(nrens,nrobs)
   lwork=2*max(3*nrens+nrens,5*nrens)
   allocate (U(nrobs,1), sig(nrmin), VT(nrens,nrens), work(lwork))
   sig=0.0
   call dgesvd('N', 'A', nrobs, nrens, X2, nrobs, sig, U, nrobs, VT, nrens, work, lwork, ierr)
   deallocate(work)
   if (ierr /= 0) then
      print *,'ierr from call dgesvd 2 = ',ierr
      stop
   endif

   allocate(isigma(nrmin))
   isigma=1.0
   do i=1,nrmin
      if ( sig(i) > 1.0 ) print *,'WARNING (analysis 4): sigma > 1',i,sig(i)
      isigma(i)=sqrt( max(1.0-sig(i)**2,0.0) )
   enddo

   do j=1,nrens
      X3(:,j)=VT(j,:)
   enddo

   do j=1,nrmin
      X3(:,j)=X3(:,j)*isigma(j)
   enddo

   print '(a)','A4: sig: '
   print '(5g11.3)',sig(1:nrmin)

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

   deallocate(U,sig,VT,isigma)


end subroutine analysis4c
