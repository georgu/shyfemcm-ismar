subroutine analysis4(A, R, S, d, ndim, nrens, nrobs, verbose)
! Computes the analysed ensemble for A using the square root formulation
! Algorithm from Appendix A1 in Evensen 2004

   use m_multa
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
   

   real X1(nrobs,nrobs),X2(nrobs,nrens),X3(nrens,nrens)
   real y1(nrobs),y2(nrobs),y3(nrobs),y4(nrens) 

   real, allocatable :: U(:,:),V(:,:),VT(:,:),sig(:),work(:)
   real, allocatable, dimension(:)   :: isigma

   real, allocatable :: ap(:)

   real sigsum,sigsum1,oneobs(1,1)
   integer ierr,nrsigma,i,j,lwork, nrmin,m,k
   integer iblkmax
   character(len=2) tag2


! Subtract mean from ensemble
   ave(:)=A(:,1)
   do i=2,nrens
      ave(:)=ave(:)+A(:,i)
   enddo
   ave=(1.0/real(nrens))*ave
   do i=1,nrens
      A(:,i)=A(:,i)-ave(:)
   enddo



! Evaluate R= S*S' + (nrens)*R
   call dgemm('n','t',nrobs,nrobs,nrens, &
                 1.0, S, nrobs, &
                      S, nrobs, &
        float(nrens-1), R, nrobs)



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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update mean
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call dgemv('t',nrobs,nrobs,1.0,Z,nrobs,d,1,0.0,y1 ,1)
   y2=eig*y1
   call dgemv('n',nrobs,nrobs,1.0,Z ,nrobs,y2,1,0.0,y3 ,1)
   call dgemv('t',nrobs,nrens,1.0,S ,nrobs,y3,1,0.0,y4 ,1)
   call dgemv('n',ndim ,nrens,1.0,A ,ndim ,y4,1,1.0,ave,1)


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

! Final ensemble perturbations
   iblkmax=min(ndim,200)
   call  multa(A, X3, ndim, nrens, iblkmax )

   
   do i=1,nrens
      A(:,i)=A(:,i)+ave(:)
   enddo

   deallocate(U,sig,VT,isigma)

end subroutine analysis4
