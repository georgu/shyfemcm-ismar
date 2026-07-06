! ======================================================================
! DUMBER STUBS / DUMMY SUBROUTINES
! Required to resolve missing diagnostic symbols in SHYFEM's elems_dealing
! ======================================================================
subroutine check_set_unit(u)
   implicit none
   integer, intent(in) :: u
   ! Dummy routine - Do nothing
end subroutine check_set_unit

subroutine check_elem(e)
   implicit none
   integer, intent(in) :: e
   ! Dummy routine - Do nothing
end subroutine check_elem

subroutine check_nodes_in_elem(e)
   implicit none
   integer, intent(in) :: e
   ! Dummy routine - Do nothing
end subroutine check_nodes_in_elem

subroutine check_node(n)
   implicit none
   integer, intent(in) :: n
   ! Dummy routine - Do nothing
end subroutine check_node

subroutine check_elems_around_node(n)
   implicit none
   integer, intent(in) :: n
   ! Dummy routine - Do nothing
end subroutine check_elems_around_node

