module mod_states
! Model state definition for HYCOM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module defines a "structure" or type definition which is used
! to hold, and operate on, a full model state.  It uses the grid dimensions
! as defined in "mod_dimensions".   It also "overloads" the operators, +, -,
! *, and = (which in Fortran90 by default only work between reals and integers),
! to also work between variables of the type states.   This makes it possible
! to define a variable and write simple code as follows.
!    type(states) member(2), average, product
!    average=member(1)+member(2)
!    product=member(1)*member(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use mod_dimensions
   type states
      real u(nx,ny,nz)     ! 3-D u-velocity
      real v(nx,ny,nz)     ! 3-D v-velocity
      real d(nx,ny,nz)     ! 3-D Layer thickness
      real t(nx,ny,nz)     ! 3-D Temperature
      real s(nx,ny,nz)     ! 3-D Temperature 
      real ub(nx,ny)       ! 2-D barotropic u-velocity
      real vb(nx,ny)       ! 2-D barotropic v-velocity
      real pb(nx,ny)       ! 2-D barotropic pressure
   end type states
   integer, parameter ::  global_ndim = 5*nx*ny*nz+3*nx*ny  ! Dimension of states

! Overloaded and generic operators
   interface operator(+)
      module procedure add_states
   end interface

   interface operator(-)
      module procedure subtract_states
   end interface

   interface operator(*)
      module procedure states_real_mult,&
                       real_states_mult,&
                       states_states_mult
   end interface

   interface assignment(=)
      module procedure assign_states
   end interface


contains
   function add_states(A,B)
      type(states) add_states
      type(states), intent(in) :: A
      type(states), intent(in) :: B
       add_states%u = A%u + B%u
       add_states%v = A%v + B%v
       add_states%d = A%d + B%d
       add_states%t = A%t + B%t
       add_states%s = A%s + B%s
       add_states%ub = A%ub + B%ub
       add_states%vb = A%vb + B%vb
       add_states%pb = A%pb + B%pb
   end function add_states

   function subtract_states(A,B)
      type(states) subtract_states
      type(states), intent(in) :: A
      type(states), intent(in) :: B
       subtract_states%u = A%u - B%u
       subtract_states%v = A%v - B%v
       subtract_states%d = A%d - B%d
       subtract_states%t = A%t - B%t
       subtract_states%s = A%s - B%s
       subtract_states%ub = A%ub - B%ub
       subtract_states%vb = A%vb - B%vb
       subtract_states%pb = A%pb - B%pb
   end function subtract_states

   function states_real_mult(A,B)
      type(states) states_real_mult
      type(states), intent(in) :: A
      real, intent(in) :: B
       states_real_mult%u = B*A%u
       states_real_mult%v = B*A%v
       states_real_mult%d = B*A%d
       states_real_mult%t = B*A%t
       states_real_mult%s = B*A%s
       states_real_mult%ub = B*A%ub
       states_real_mult%vb = B*A%vb
       states_real_mult%pb = B*A%pb
   end function states_real_mult

   function real_states_mult(B,A)
      type(states) real_states_mult
      type(states), intent(in) :: A
      real, intent(in) :: B
       real_states_mult%u = B*A%u
       real_states_mult%v = B*A%v
       real_states_mult%d = B*A%d
       real_states_mult%t = B*A%t
       real_states_mult%s = B*A%s
       real_states_mult%ub = B*A%ub
       real_states_mult%vb = B*A%vb
       real_states_mult%pb = B*A%pb
   end function real_states_mult

   function states_states_mult(A,B)
      type(states) states_states_mult
      type(states), intent(in) :: A
      type(states), intent(in) :: B
       states_states_mult%u = A%u * B%u
       states_states_mult%v = A%v * B%v
       states_states_mult%d = A%d * B%d
       states_states_mult%t = A%t * B%t
       states_states_mult%s = A%s * B%s
       states_states_mult%ub = A%ub * B%ub
       states_states_mult%vb = A%vb * B%vb
       states_states_mult%pb = A%pb * B%pb
   end function states_states_mult

   subroutine assign_states(A,r)
      type(states), intent(out) :: A
      real, intent(in) :: r
       A%u = r
       A%v = r
       A%d = r
       A%t = r
       A%s = r
       A%ub = r
       A%vb = r
       A%pb = r
   end subroutine assign_states

end module mod_states

