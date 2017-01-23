module mod_observations
! sea level obs (not correlated)
   type levels
      double precision, allocatable :: t(:)          ! time
      real, allocatable :: x(:)                      ! x coord
      real, allocatable :: y(:)                      ! y coord
      real, allocatable :: val(:)                    ! value
      real, allocatable :: std(:)                    ! std value
   end type levels
end module mod_observations
