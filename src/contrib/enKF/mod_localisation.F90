!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
!
module mod_localisation

  use mod_dimension
  use mod_para
  use mod_ens_state

  implicit none

contains

!*************************************************************

  subroutine analysis_nkn

  implicit none
  
  integer ik_list(nnkn) ! list of nodes inside the radius
  real rho(nnkn)	! value of the distance function
  integer k

  do k = 1,nnkn
     !call find_nodes(r_local,ik_list,rho)
     !call make_local
     !call analysis(??)
     !call add_local
  end do

  end subroutine analysis_nkn


!*************************************************************

  subroutine analysis_nel

  end subroutine analysis_nel

end module mod_localisation

