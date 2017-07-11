module mod_observations
! observation files
!
   type files
       integer :: id	! id number
       character(len=5) :: ty	!type of file
       character(len=80) :: name	!name of the file
   end type files
! sea level obs (not correlated)
!
   type levels
      double precision :: t          ! time
      real :: x                      ! x coord
      real :: y                      ! y coord
      real :: val                    ! value
      real :: std                    ! std value
      integer :: status		     ! = 0 can be assimilated
                                     ! = 1 bad value not to be ass
                                     ! = 2 flag value not to be ass
      integer  :: id                 ! id number of the file
   end type levels

! velocity obs
!
   type currentf
      double precision  :: t                         ! time of the field
      integer           :: nx,ny                     ! dimensions
      real, allocatable :: x(:)                      ! x coords
      real, allocatable :: y(:)                      ! y coords
      real              :: z                         ! z coord
      real, allocatable :: u(:,:)                    ! u value
      real, allocatable :: v(:,:)                    ! v value
      real, allocatable :: std(:,:)                  ! std value
      integer, allocatable :: status(:,:)	     ! = 0 can be assimilated
                                                     ! = 1 bad value not to be ass
                                                     ! = 2 flag value not to be ass
      integer           :: id			     ! id number of the file
   end type currentf

end module mod_observations
