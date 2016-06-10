program sample
! Code for testing sampling of pseudo random fields on different architectures.
   use m_sample1D
   use m_sample2D
   use m_tecfld
   use m_set_random_seed2
   implicit none


   real :: xlength=10000.0
   real :: ylength=10000.0
   integer, parameter :: nrens=250
   integer, parameter :: nx=100
   integer, parameter :: ny=100
   integer, parameter :: nre=1

   real :: cor1=1500.0
   real :: cor2=600.0
   real dir

   integer i
   real, parameter :: pi=3.1415927

   real dx,dy
   real, allocatable :: Amat(:,:,:),B(:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Set a variable random seed
   call set_random_seed2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   cor1=cor1/sqrt(3.0)
   cor2=cor2/sqrt(3.0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   dx=xlength/float(nx-1)
   dy=ylength/float(ny-1)

   print *,'dx=',dx
   print *,'dy=',dy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Sample 2D pseudo random fields
   allocate(Amat(nx,ny,nrens))

   do i=1,1
   dir=20.0
   print '(a,f10.2)','calling sample2D with dir=',dir
   call sample2D(Amat,nx,ny,nrens,nre,dx,dy,cor1,cor2,dir,.true.,.true.)
   call tecfld('ens2DA',nx,ny,min(10,nrens),Amat)

   dir=60.0
   print '(a,f10.2)','calling sample2D with dir=',dir
   call sample2D(Amat,nx,ny,nrens,nre,dx,dy,cor1,cor2,dir,.true.,.true.)
   call tecfld('ens2DB',nx,ny,min(10,nrens),Amat)
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Sample all random fields and interpolate to model grid
   allocate(B(nx,nrens))

   print '(a,f10.2)','calling sample1D'
   call sample1D(B,nx,nrens,nre,dx,cor1,.true.,.false.)
   open(10,file='tec_ens1D.dat')
   do i=1,nx
      write(10,'(i3,20f12.2)')i,(i-1)*dx,B(i,1:min(nrens,10))
   enddo
   close(10)

end program
