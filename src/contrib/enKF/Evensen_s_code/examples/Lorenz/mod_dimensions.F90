module mod_dimensions
   integer, parameter :: nrsamp=1000     ! Ensemble size
   integer, parameter :: ndim=801        ! Dimension of state
   integer, parameter :: nfft=1260       ! (tab p 647)
!   integer, parameter :: ndim=2401       !801        ! Dimension of state
!   integer, parameter :: nfft=3840       ! 1260 (tab p 647)
   integer, parameter :: nrmes=80        !80 160    ! number of measurements
   integer, parameter :: nrobs=3         ! number of mes at each time
   integer, parameter :: nrtobs=nrmes*nrobs   ! total number of mes 
   integer, parameter :: n0=ndim
   integer, parameter :: n1=ndim-1
   integer, parameter :: n2=ndim-2
   integer, parameter :: n3=ndim-3
   integer, parameter :: n4=ndim-4
   integer, parameter :: n5=ndim-5
   integer, parameter :: n6=ndim-6
   real, parameter   ::  pi=3.141592653589
   real, parameter   ::  sigma=10.0          ! Prantl number
   real, parameter   ::  r=28.0              ! Normalized Rayleigh number
   real, parameter   ::  b=8.0/3.0           ! Nondimensional wave number

   type variances
      real dyn
      real mes
      real ini
   end type variances

end module mod_dimensions
