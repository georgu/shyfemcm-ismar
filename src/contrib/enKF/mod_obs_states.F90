!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
!
module mod_obs_states

! stat of an observation:
! 0 = normal obs (assimilated)
! 1 = super-observation (assimilated)
! 2 = observation merged into a super-observation (not-assimilated)
! 3 = observation out of range (not-assimilated)
! 4 = observation with a flag value (not-assimilated)

! observation files
!
   type files
       character(len=5) :: ty	!type of file
       character(len=80) :: name	!name of the file
   end type files
   
! scalar variable in one location (0D)
!
   type scalar_0d
      double precision :: t          ! time
      real :: x                      ! x coord
      real :: y                      ! y coord
      real :: z                      ! z coord
      real :: val                    ! value
      real :: std                    ! std value
      integer :: stat		     ! = 0,1,2,3,4
      integer  :: id                 ! id number of the file
   end type scalar_0d

! vector variable in a field (2D)
!
   type vector_2d
      double precision  :: t                         ! time of the field
      integer           :: nx,ny                     ! dimensions
      real, allocatable :: x(:,:)                      ! x coords
      real, allocatable :: y(:,:)                      ! y coords
      real              :: z                         ! z coord
      real, allocatable :: u(:,:)                    ! u value
      real, allocatable :: v(:,:)                    ! v value
      real, allocatable :: std(:,:)                  ! std value
      integer, allocatable :: stat(:,:)	     ! = 0,1,2,3,4
      integer           :: id			     ! id number of the file
   end type vector_2d

end module mod_obs_states
