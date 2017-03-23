subroutine analysis2(A, D, R, S, ndim, nrens, nrobs, verbose)
! Computes the analysed ensemble for A 
! Uses the exact algorithm from Evensen 2003 where the
! error covariance matrix is specified.
! If the error covariance matrix is represented by a low rank
! ensemble it will lead to a loss of rank when
! m > N.  Includes the perturbations of measurements.
   use m_multa
   implicit none
   integer, intent(in) :: ndim             ! dimension of model state
   integer, intent(in) :: nrens            ! number of ensemble members
   integer, intent(in) :: nrobs            ! number of observations
   
   real, intent(inout) :: A(ndim,nrens)   ! ensemble matrix
   real, intent(in)    :: D(nrobs,nrens)  ! vector holding observation innovations
   real, intent(in)    :: S(nrobs,nrens)  ! matrix holding HA' 
   real, intent(inout) :: R(nrobs,nrobs)  ! Error covariance matrix for observations
   logical, intent(in) :: verbose


   real, allocatable, dimension(:,:) :: X1,X2,U,X4,Reps,V
   real X3(nrobs,nrens)
   real, allocatable, dimension(:)   :: sig,work

   real sigsum,sigsum1,oneobs(1,1)
   integer ierr,nrsigma,i,j,lwork, m
   integer iblkmax
   character(len=2) tag2



   if (nrobs > 1) then
!      R=float(nrens)*R+matmul(S,transpose(S))
      call dgemm('n','t', nrobs, nrobs, nrens, &
                 1.0, S, nrobs, &
                      S, nrobs, &
        float(nrens-1), R, nrobs)



      allocate (U(nrobs,nrobs)  )
      allocate (V(nrobs,nrobs)  )
      allocate (sig(nrobs)  )
      lwork=2*max(3*nrobs+nrobs,5*nrobs)
      allocate(work(lwork))
      sig=0.0
      call dgesvd('A', 'A', nrobs, nrobs, R, nrobs, sig, U, nrobs, V, nrobs, work, lwork, ierr)
      deallocate(work)
      if (ierr /= 0) then
         print *,'ierr from call dgesvd= ',ierr
         stop
      endif

      open(10,file='sigma.dat')
         do i=1,nrobs
            write(10,'(i5,g12.3)')i,sig(i)
         enddo
      close(10)

      sigsum=sum( sig(1:nrobs) )
      sigsum1=0.0
   ! Significant eigenvalues.
      nrsigma=0
      do i=1,nrobs                 ! singular values are in descending order
         if (sigsum1/sigsum < 0.999) then
            nrsigma=nrsigma+1
            sigsum1=sigsum1+sig(i)
            sig(i) = 1.0/sig(i)
         else
            sig(i:nrobs)=0.0
            exit
         endif
      enddo

      if (verbose) then
         write(*,'(a,i5,g13.5)') ' dominant sing. values and share ',nrsigma,sigsum1/sigsum
         write(*,'(8g12.3)')1./sig(1:nrsigma), sig(nrsigma+1:nrobs)
      endif

      allocate (X1(nrobs,nrobs))
      do i=1,nrobs
      do j=1,nrobs
         X1(i,j)=sig(i)*U(j,i)
      enddo
      enddo
      deallocate(sig)

      allocate (X2(nrobs,nrens))
!     X2=matmul(X1,D)
      call dgemm('n','n',nrobs,nrens,nrobs,1.0,X1,nrobs,D ,nrobs,0.0,X2,nrobs)
      deallocate(X1) 

!     X3=matmul(V,X2)
      call dgemm('t','n',nrobs,nrens,nrobs,1.0,V ,nrobs,X2,nrobs,0.0,X3,nrobs)
      deallocate(V)
      deallocate(X2)

   else
      oneobs=matmul(S,transpose(S))+float(nrens-1)*R
      print *,'oneobs: ',oneobs(1,1)
      X3=D/oneobs(1,1)
   endif

   if (2_8*ndim*nrobs < 1_8*nrens*(nrobs+ndim)) then
!    Code for few observations ( m<nN/(2n-N) )
      if (verbose) print * ,'analysis: Representer approach is used'
      allocate (Reps(ndim,nrobs))

!    Reps=matmul(A,transpose(S))
      call dgemm('n','t',ndim,nrobs,nrens,1.0,A,ndim,S,nrobs,0.0,Reps,ndim)
!      call printreps(Reps,nrobs)

!    A=A+matmul(Reps,X3)
      call dgemm('n','n',ndim,nrens,nrobs,1.0,Reps,ndim,X3,nrobs,1.0,A,ndim)
      deallocate(Reps)

      tag2(1:2)='X3'
      open(10,file='X.uf',form='unformatted')
         write(10)tag2,nrens,nrobs,X3,S
      close(10)

   else
      if (verbose) print * ,'analysis: X5 appraoch is used'
      allocate(X4(nrens,nrens))
!      X4=matmul(transpose(S),X3)
      call dgemm('t','n',nrens,nrens,nrobs,1.0,S,nrobs,X3,nrobs,0.0,X4,nrens)
      do i=1,nrens
         X4(i,i)=X4(i,i)+1.0  ! Addition with Af
      enddo

      iblkmax=min(ndim,200)
      call  multa(A, X4, ndim, nrens, iblkmax )

      tag2(1:2)='X5'
      open(10,file='X.uf',form='unformatted')
         write(10)tag2,nrens,X4
      close(10)

      open(10,file='X5.dat')
         write(10,*)'TITLE = "X5"'
         write(10,*)'VARIABLES = "iens_f" "iens_a" "X5"'
         write(10,'(2(a,i5),a)')' ZONE  F=BLOCK, I=',nrens,', J=',nrens,', K=1'
         write(10,'(30I5)')((i,i=1,nrens),j=1,nrens)
         write(10,'(30I5)')((j,i=1,nrens),j=1,nrens)
         write(10,'(10(1x,e12.5))')((X4(i,j),i=1,nrens),j=1,nrens)
      close(10)

      open(10,file='X5col.dat')
         write(10,'(a)')'These should all be ones'
         do j=1,nrens
            write(10,'(i5,f10.4)')j,sum(X4(:,j))  
         enddo
      close(10)

      open(10,file='X5row.dat')
         do j=1,nrens
            write(10,'(i5,f10.4)')j,sum(X4(j,:))/float(nrens)
          enddo
      close(10)

      deallocate(X4)
   endif 


end subroutine analysis2
