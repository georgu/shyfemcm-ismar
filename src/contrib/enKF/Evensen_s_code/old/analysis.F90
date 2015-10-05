subroutine analysis(A, D, E, S, ndim, nrens, nrobs)
! This routine is obsolete!
! It should be replaced by the analysis6.F90 routine
! which is equally efficient but do not use the approximate
! factorization from Evensen 2003.

! Computes the analysed ensemble in the EnKF
! Written by G. Evensen (Geir.Evensen@nersc.no)
! This routine uses subroutines from BLAS and EISPACK
! and calls the additional multiplication routine multa.
   use m_multa
   implicit none

! Dimension of model state
   integer, intent(in) :: ndim            

! Number of ensemble members
   integer, intent(in) :: nrens            

! Number of observations
   integer, intent(in) :: nrobs            
   
! Ensemble matrix
   real, intent(inout) :: A(ndim,nrens)    

! Matrix holding innovations
   real, intent(in)    :: D(nrobs,nrens)   

! Matrix holding HA' 
   real, intent(in)    :: S(nrobs,nrens)   

! Matrix holding observation perturbations
   real, intent(in)    :: E(nrobs,nrens)   


! Local variables
   real,    allocatable, dimension(:,:) :: &
       X1, X2, U, X4, Reps

   real,    allocatable, dimension(:)   :: &
       sig, work

   real ES(nrobs,nrens), X3(nrobs,nrens), V(nrens,nrens)
   real sigsum,sigsum1
   integer ierr, nrsigma, i, j, lwork, nrmin, iblkmax



   print *,'analysis.F90: please replace me with the new fast analysis6.F90 routine.'
   stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do nothing if only one measurement
   if (nrobs == 1) then
      print *,'analysis: no support for nrobs=1'
      return
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Minimum of nrobs and nrens
   nrmin=min(nrobs,nrens)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute HA'+E
   ES=S+E

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute SVD of HA'+E  ->  U and sig, using Eispack
   allocate (U(nrobs,nrmin))
   allocate (sig(nrmin))
   lwork=2*max(3*nrens+nrobs,5*nrens)
   allocate(work(lwork))

   sig=0.0
   call dgesvd('S', 'N', nrobs, nrens, ES, nrobs, sig,  &
               U, nrobs, V, nrens, work, lwork, ierr)
   deallocate(work)
   if (ierr /= 0) then
      print *,'ierr from call dgesvd= ',ierr
      stop
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Convert to eigenvalues
   do i=1,nrmin
      sig(i)=sig(i)**2
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute number of significant eigenvalues
   sigsum=sum(sig(1:nrmin))
   sigsum1=0.0
   nrsigma=0
   do i=1,nrmin                
      if (sigsum1/sigsum < 0.999) then
         nrsigma=nrsigma+1
         sigsum1=sigsum1+sig(i)
         sig(i) = 1.0/sig(i)
      else
         sig(i:nrmin)=0.0
         exit
      endif
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute X1
   allocate (X1(nrmin,nrobs))
   do j=1,nrobs
   do i=1,nrmin
      X1(i,j)=sig(i)*U(j,i)
   enddo
   enddo
   deallocate(sig)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute X2=X1*D
   allocate (X2(nrmin,nrens))
   call dgemm('n', 'n', nrmin, nrens, nrobs, 1.0, X1, &
              nrmin, D, nrobs, 0.0, X2, nrmin)
   deallocate(X1) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute X3=U*X2
   call dgemm('n', 'n', nrobs, nrens, nrmin, 1.0, U,  &
              nrobs, X2, nrmin, 0.0, X3, nrobs)
   deallocate(U)
   deallocate(X2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute final analysis 
   if (2*ndim*nrobs > nrens*(nrobs+ndim)) then
!  Case with nrobs 'large'

!     Compute X4=(HA')^T * X3
      allocate(X4(nrens,nrens))
      call dgemm('t', 'n', nrens, nrens, nrobs, 1.0,  &
                 S, nrobs, X3, nrobs, 0.0, X4, nrens)
 
!     Compute X5=X4+I (stored in X4)
      do i=1,nrens
         X4(i,i)=X4(i,i)+1.0  
      enddo

!     Compute A=A*X5
      iblkmax=min(ndim,200)
      call  multa(A, X4, ndim, nrens, iblkmax)
      deallocate(X4)

   else
!  Case with nrobs 'small'

!     Compute representers Reps=A'*S^T
      allocate (Reps(ndim,nrobs)) 
      call dgemm('n', 't', ndim, nrobs, nrens, 1.0, A, &
                 ndim, S, nrobs, 0.0, Reps, ndim)

!     Compute A=A+Reps*X3
      call dgemm('n', 'n', ndim, nrens, nrobs, 1.0,    &
                 Reps, ndim, X3, nrobs, 1.0, A, ndim)
      deallocate(Reps)
   endif 


end subroutine analysis
