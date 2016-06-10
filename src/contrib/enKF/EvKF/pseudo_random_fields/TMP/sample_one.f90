












program sample
! Code for testing sampling of pseudo random fields on different architectures.
   use m_sample1D
   use m_tecfld
   use m_set_random_seed2
   implicit none


!   real :: xlength=3600.*8760
   integer, parameter :: nrens=4
   integer, parameter :: nx=8760
   integer, parameter :: nre=1

   real :: cor1=1500.0

   integer i
   real, parameter :: pi=3.1415927

   real dx
   real, allocatable :: B(:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Set a variable random seed
   !call set_random_seed2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   cor1=cor1/sqrt(3.0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   dx=3600.

   print *,'dx=',dx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Sample all random fields and interpolate to model grid
   allocate(B(nx,nrens))

   print '(a,f10.2)','calling sample1D'
   call sample1D(B,nx,nrens,nre,dx,cor1,.true.,.false.)
   open(10,file='ens1D_nad.dat')
   do i=1,nx
      write(10,'(i3,20f12.2)')i,(i-1)*dx,B(i,1:min(nrens,10))
   enddo
   close(10)

end program
