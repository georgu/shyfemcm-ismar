module m_pseudo2D

#ifdef IBM
   integer :: naux1,naux2,naux3,nn,s1
   double precision, dimension(:), allocatable,save :: aux1,aux2,aux3
#endif


contains

subroutine pseudo2D(Amat,nx,ny,lde,rx,ry,dx,dy,n1,n2,theta,verbose)
! This routine calculates the pseudo random filds using
! the procedure outlined in Evensen (1994) \cite{eve94a}.

#ifdef DEC
!   use mydxml
#endif

#ifdef LINUX
   use mod_fftw3
#endif
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

#ifdef CRAY
   real, allocatable, dimension(:), save :: work 
   real, allocatable, dimension(:), save :: table 
#endif

#ifdef DEC
   integer status
   record /dxml_d_fft_structure_2d/ fft_struct
#endif

#ifdef SGI
   real, allocatable, dimension(:), save :: coeff
#endif

   real r1,r2,c

#ifdef LINUX
   integer*8 plan
#endif

   integer l,p,j,m,i
   real kappa2,lambda2,kappa,lambda
   real pi2,deltak,summ,scale
   real a11tmp,a22tmp,a11,a22,a12,torad
#ifdef IBM
   logical ln1,ln2
#endif

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

#ifndef CRAY
#ifndef IBM
#ifndef DEC
#ifndef SGI
#ifndef LINUX
   print *,'ranfield is only running on the following machines:'
   print *,'   CRAY'
   print *,'   IBM having essl'
   print *,'   DEC having dxml'
   print *,'   SGI'
   print *,'   LINUX having FFTW3'
   stop
#endif
#endif
#endif
#endif
#endif

   pi2=2.0*pi
   deltak=pi2**2/(real(n1*n2)*dx*dy)
   kappa=pi2/(real(n1)*dx)
   kappa2=kappa**2
   lambda=pi2/(real(n2)*dy)
   lambda2=lambda**2
   scale=1.0


#ifdef SGI
   if (mod(n1,2) /= 0) print *,'pseudo2D: mod(n1,2) must be zero, n1=',n1
   if (mod(n2,2) /= 0) print *,'pseudo2D: mod(n2,2) must be zero, n2=',n2
   if(allocated(coeff)) deallocate(coeff)
   allocate( coeff((n1+15) + 2*(n2+15)) )
   call dzfft2dui(n1,n2,coeff)
#endif

#ifdef CRAY
   if(allocated(work)) deallocate(work)
   allocate(work(512*n1))

   if(allocated(table)) deallocate(table)
   allocate(table(100+2*(n1+n2)))

   call scfft2d(0,n1,n2,scale,x,n1/2+1,y,n1+2,table,work,0)
#endif

#ifdef IBM
   ln1=.false.
   ln2=.false.
   do i=1,20
      if (n1 == 2**i) ln1=.true.
      if (n2 == 2**i) ln2=.true.
   enddo

   if (.not.ln1 .or. .not.ln2) then
      print *,'pseudo2D: For IBM ESSL routines n1 and n2 must be equal to 2^i'
      stop
   endif

   nn   = max(n1/2,n2)

   if (nn<=2048) then
      naux1= 42000 
   else
      naux1= ceiling(40000+1.64*n1+2.28*n2)
   end if

   if (n1 <= 4096 ) then
      naux2 = 20000
   else if (n1 > 4096 ) then 
      naux2 = ceiling(20000+1.14*n1)
   end if

   if ( n2 > 252) then 
      s1 = min(64, 1+n1/2)
      naux2 = naux2 + ceiling((2*n2+256)*(2.28+s1))
   end if

   naux3=1

   if (allocated(aux1)) deallocate(aux1); allocate(aux1(naux1))
   if (allocated(aux2)) deallocate(aux2); allocate(aux2(naux2))
   if (allocated(aux3)) deallocate(aux3); allocate(aux3(naux3))

   call dcrft2(1,x,n1/2+1,y,n1+2,n1,n2,-1,scale,aux1,naux1,aux2,naux2,aux3,naux3)
#endif

#ifdef DEC
   if (mod(n1,2) /= 0) then
      print *,'ranfield: n1 is not even. n1=',n1
   endif
   status=dfft_init_2d(n1,n2,fft_struct,.true.)
   if (status /= 0 ) print *,'status: dfft_init_2d',status
#endif

#ifdef LINUX
   if (allocated(y)) deallocate(y)
   allocate(y(0:n1-1,0:n2-1))
   call dfftw_plan_dft_c2r_2d(plan,n1,n2,x,y,FFTW_ESTIMATE)
#endif


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

#ifdef CRAY
      call csfft2d(-1,n1,n2,scale,x,n1/2+1,y,n1+2,table,work,0)
#endif
#ifdef SGI
      call zdfft2du(-1,n1,n2,x,n1+2,coeff)
      y=reshape(transfer(x,(/0.0 /) ),(/n1+2,n2/))
#endif

#ifdef IBM
      call dcrft2(0,x,n1/2+1,y,n1+2,n1,n2,-1,scale,aux1,naux1,aux2,naux2,aux3,naux3)
#endif

#ifdef DEC
      status=dfft_apply_2d('C','R','B',x,y,n1+2,fft_struct,1,1)
      if (status /= 0 ) print *,'status: dfft_apply_2d',status
      y=y*real(n1*n2)
#endif

#ifdef LINUX
      call dfftw_execute(plan)
#endif

      do m=1,ny
      do i=1,nx
         Amat(i,m,j)=y(i-1,m-1)
      enddo
      enddo

   enddo

#ifdef LINUX
   call dfftw_destroy_plan(plan)
#endif

   deallocate(fampl, phi, y, x)
!#ifdef IBM
!   if(allocated(aux1)) deallocate(aux1)
!   if(allocated(aux2)) deallocate(aux2)
!   if(allocated(aux3)) deallocate(aux3)
!#endif

#ifdef DEC
!   status=dfft_exit_2d(fft_struct)
!   print *,'status: dfft_exit_2d',status
#endif

end subroutine pseudo2D
end module m_pseudo2D
