module m_pseudo2D

contains

subroutine pseudo2D(Amat,nx,ny,lde,rx,ry,dx,dy,n1,n2,theta,verbose)
! This routine calculates the pseudo random filds using
! the procedure outlined in Evensen (1994) \cite{eve94a}.


   use mod_fftw3
   use m_newton2D
   implicit none
   logical, Intent(IN) :: verbose
   integer, intent(in) :: nx,ny           ! horizontal dimensions
   integer, intent(in) :: lde             ! number of fields stored at the time
   real, intent(out)   :: Amat(nx,ny,lde) ! generated random fields
   real, intent(in)    :: rx,ry           ! Horizontal decorrelation lengths (rx in dir theta)
   real, intent(in)    :: dx,dy           ! grid spacing
   real, intent(in)    :: theta           ! rotation angle in deg (theta=0 is east, rotation anticlocwise)
   integer, intent(in) :: n1,n2           ! horizontal dimensions in fft grid

   real r1,r2,c

   integer(kind=8) plan

   integer l,p,j,m,i
   real kappa2,lambda2,kappa,lambda
   real pi2,deltak,summ,scale
   real a11tmp,a22tmp,a11,a22,a12,torad

   real, allocatable    :: fampl(:,:,:)
   real, allocatable    :: phi(:,:)
   real, allocatable    :: y(:,:)   ! Physical field
   complex, allocatable :: x(:,:)   ! Fourier amplitudes

   real, parameter :: pi=3.141592653589
   logical cnv
   real e

   if (lde < 1)    stop 'pseudo2D: error lde < 1'
   if (rx <= 0.0)  stop 'pseudo2D: error, rx <= 0.0'
   if (ry <= 0.0)  stop 'pseudo2D: error, ry <= 0.0'
   if (n1 < nx)    stop 'pseudo2D: n1 < nx'
   if (n2 < ny)    stop 'pseudo2D: n2 < ny'

   allocate(fampl(0:n1/2,-n2/2:n2/2,2))
   allocate(phi(0:n1/2,-n2/2:n2/2))
   allocate(y(0:n1+1,0:n2-1))
   allocate(x(0:n1/2,0:n2-1))

   pi2=2.0*pi
   deltak=pi2**2/(real(n1*n2)*dx*dy)
   kappa=pi2/(real(n1)*dx)
   kappa2=kappa**2
   lambda=pi2/(real(n2)*dy)
   lambda2=lambda**2
   scale=1.0

   if (allocated(y)) deallocate(y)
   allocate(y(0:n1-1,0:n2-1))
   call dfftw_plan_dft_c2r_2d(plan,n1,n2,x,y,FFTW_ESTIMATE)


! computing the coefficients r1, r2, and c
   r1=3.0/rx
   r2=3.0/ry

   if (verbose) print '(a,2f12.5,2i6,4f10.2)','pseudo2D: Call newton with ',r1,r2,n1,n2,dx,dy,rx,ry
   call newton2D(r1,r2,n1,n2,dx,dy,rx,ry,cnv,verbose)
   if (.not.cnv) then
      stop 'newton did not converge'
   endif
   
   summ=0.0
   do p=-n2/2+1,n2/2
   do l=-n1/2+1,n1/2
      summ=summ+exp(-2.0*(kappa2*real(l*l)/r1**2 + lambda2*float(p*p)/r2**2))
   enddo
   enddo

   summ=summ-1.0

   c=sqrt(1.0/(deltak*summ))

   if (verbose) then
      print *,'pseudo2D: r1=  ',r1
      print *,'pseudo2D: r2=  ',r2
      print *,'pseudo2D:  c=  ',c
   end if

! Rotation to angle theta
   a11tmp=1.0/r1**2
   a22tmp=1.0/r2**2
   torad=-pi/180.0
   a11=a11tmp*cos(theta*torad)**2 + a22tmp*sin(theta*torad)**2
   a22=a11tmp*sin(theta*torad)**2 + a22tmp*cos(theta*torad)**2
   a12=(a22tmp-a11tmp)*cos(theta*torad)*sin(theta*torad)

   do j=1,lde
      ! Calculating the random wave phases
      call random_number(phi)
      phi=pi2*phi



      ! Calculating the wave amplitues
      do p=-n2/2,n2/2
      do l=0,n1/2 
!         e=exp(-(kappa2*real(l*l)/r1**2+lambda2*float(p*p)/r2**2))
         e=exp(-( a11*kappa2*real(l*l) + 2.0*a12*kappa*lambda*float(l*p) + a22*lambda2*float(p*p) ))
         fampl(l,p,1)=e*cos(phi(l,p))*sqrt(deltak)*c
         fampl(l,p,2)=e*sin(phi(l,p))*sqrt(deltak)*c
      enddo
      enddo
      fampl(0,0,1)=0.0
      fampl(0,0,2)=0.0

      do p=0,n2/2-1
         x(:,p)=cmplx(fampl(:,p,1),fampl(:,p,2))
      enddo

      do p=n2/2,n2-1
         x(:,p)=cmplx(fampl(:,-n2+p,1),fampl(:,-n2+p,2))
      enddo

      call dfftw_execute(plan)

      do m=1,ny
      do i=1,nx
         Amat(i,m,j)=y(i-1,m-1)
      enddo
      enddo

   enddo

   call dfftw_destroy_plan(plan)

   deallocate(fampl, phi, y, x)

end subroutine pseudo2D
end module m_pseudo2D
