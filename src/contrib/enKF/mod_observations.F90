module mod_observations
! sea level obs (not correlated)
   type obs_level
      double precision t          ! Time
      real x                      ! x coord
      real y                      ! y coord
      real val                    ! value
      real std                ! std value
   end type obs_level

end module mod_observations
