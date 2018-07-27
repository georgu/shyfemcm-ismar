!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
!
module mod_mod_states

   use mod_dimensions
   
! standard model state
   type states
      real u(nnlv,nnel)                      ! 3-D u-velocity
      real v(nnlv,nnel)                      ! 3-D v-velocity
      real z(nnkn) 	                     ! 2-D water level 
      real t(nnlv,nnkn)                      ! 3-D Temperature
      real s(nnlv,nnkn)                      ! 3-D Salinity 
   end type states
   integer, save ::  global_ndim = 2*nnlv*nnel + nnkn + 2*nnlv*nnkn
                      

 ! Dimension of states

! single precision model state (used for read and write to files)
   type states4
      real*4 u(nnlv,nnel)                      ! 3-D u-velocity
      real*4 v(nnlv,nnel)                      ! 3-D v-velocity
      real*4 z(nnkn) 	                     ! 2-D water level 
      real*4 t(nnlv,nnkn)                      ! 3-D Temperature
      real*4 s(nnlv,nnkn)                      ! 3-D Salinity 
   end type states4

! Since we have a staggered grid we must use two local states.
! It should be still a good assumption...
! model state at one node (used in cov localisation)
   type sub_states_nkn
      real z 	                     	       ! 0-D water level 
      real t(nnlv)                             ! 1-D Temperature
      real s(nnlv)                             ! 1-D Salinity
   end type sub_states_nkn
   integer, parameter ::  local_ndim_nkn = 2*nnlv + 1
! model state at one element (used in cov localisation)
   type sub_states_nel
      real u(nnlv)                             ! 1-D u-velocity
      real v(nnlv)                             ! 1-D v-velocity
   end type sub_states_nel
   integer, parameter ::  local_ndim_nel = 2*nnlv 

! double state, with the model errors
     type qstate
      real qu(nnlv,nnel)                       ! 3-D u-velocity error
      real qv(nnlv,nnel)                       ! 3-D v-velocity error
      real qz(nnkn) 	                     ! 2-D water level 
      real qt(nnlv,nnkn)                       ! 3-D Temperature error
      real qs(nnlv,nnkn)                       ! 3-D Salinity error

      real u(nnlv,nnel)                        ! 3-D u-velocity
      real v(nnlv,nnel)                        ! 3-D v-velocity
      real z(nnkn) 	                     ! 2-D water level 
      real t(nnlv,nnkn)                        ! 3-D Temperature
      real s(nnlv,nnkn)                        ! 3-D Salinity 
   end type qstate
 
!-----------
! Overloaded and generic operators
   interface operator(+)
      module procedure add_states,&
		       states_real_add,&
		       real_states_add
   end interface

   interface operator(-)
      module procedure subtract_states
   end interface

   interface operator(*)
      module procedure states_real_mult,&
                       real_states_mult,&
                       states_states_mult
   end interface

   interface operator(/)
      module procedure states_states_div
   end interface

   interface assignment(=)
      module procedure assign_states
      module procedure states4to8
      module procedure states8to4
   end interface


contains


   function add_states(A,B)
      type(states) add_states
      type(states), intent(in) :: A
      type(states), intent(in) :: B
       add_states%u = A%u + B%u
       add_states%v = A%v + B%v
       add_states%z = A%z + B%z
       add_states%t = A%t + B%t
       add_states%s = A%s + B%s
   end function add_states

   function subtract_states(A,B)
      type(states) subtract_states
      type(states), intent(in) :: A
      type(states), intent(in) :: B
       subtract_states%u = A%u - B%u
       subtract_states%v = A%v - B%v
       subtract_states%z = A%z - B%z
       subtract_states%t = A%t - B%t
       subtract_states%s = A%s - B%s
   end function subtract_states

   function states_real_mult(A,B)
      type(states) states_real_mult
      type(states), intent(in) :: A
      real, intent(in) :: B
       states_real_mult%u = B*A%u
       states_real_mult%v = B*A%v
       states_real_mult%z = B*A%z
       states_real_mult%t = B*A%t
       states_real_mult%s = B*A%s
   end function states_real_mult

   function real_states_mult(B,A)
      type(states) real_states_mult
      type(states), intent(in) :: A
      real, intent(in) :: B
       real_states_mult%u = B*A%u
       real_states_mult%v = B*A%v
       real_states_mult%z = B*A%z
       real_states_mult%t = B*A%t
       real_states_mult%s = B*A%s
   end function real_states_mult

   function states_real_add(A,B)
      type(states) states_real_add
      type(states), intent(in) :: A
      real, intent(in) :: B
       states_real_add%u = B + A%u
       states_real_add%v = B + A%v
       states_real_add%z = B + A%z
       states_real_add%t = B + A%t
       states_real_add%s = B + A%s
   end function states_real_add


   function real_states_add(B,A)
      type(states) real_states_add
      type(states), intent(in) :: A
      real, intent(in) :: B
       real_states_add%u = B + A%u
       real_states_add%v = B + A%v
       real_states_add%z = B + A%z
       real_states_add%t = B + A%t
       real_states_add%s = B + A%s
   end function real_states_add


   function states_states_mult(A,B)
      type(states) states_states_mult
      type(states), intent(in) :: A
      type(states), intent(in) :: B
       states_states_mult%u = A%u * B%u
       states_states_mult%v = A%v * B%v
       states_states_mult%z = A%z * B%z
       states_states_mult%t = A%t * B%t
       states_states_mult%s = A%s * B%s
   end function states_states_mult

   function states_states_div(A,B)
      type(states) states_states_div
      type(states), intent(in) :: A
      type(states), intent(in) :: B
       states_states_div%u = A%u / B%u
       states_states_div%v = A%v / B%v
       states_states_div%z = A%z / B%z
       states_states_div%t = A%t / B%t
       states_states_div%s = A%s / B%s
   end function states_states_div

   function root_state(A)
      type(states) root_state
      type(states), intent(in) :: A
       root_state%u = sqrt(A%u)
       root_state%v = sqrt(A%v)
       root_state%z = sqrt(A%z)
       root_state%t = sqrt(A%t)
       root_state%s = sqrt(A%s)
   end function root_state

   subroutine assign_states(A,r)
      type(states), intent(out) :: A
      real, intent(in) :: r
       A%u = r
       A%v = r
       A%z = r
       A%t = r
       A%s = r
   end subroutine assign_states

   subroutine states4to8(A,B)
      type(states), intent(out) :: A
      type(states4), intent(in)  :: B
      A%u=DBLE(B%u)
      A%v=DBLE(B%v)
      A%z=DBLE(B%z)
      A%t=DBLE(B%t)
      A%s=DBLE(B%s)
   end subroutine states4to8

   subroutine states8to4(A,B)
      type(states), intent(in)  :: B
      type(states4),  intent(out) :: A
      A%u=real(B%u)
      A%v=real(B%v)
      A%z=real(B%z)
      A%t=real(B%t)
      A%s=real(B%s)
   end subroutine states8to4

   subroutine push_qstate(A,B,C)
      implicit none
      type(states), intent(in)  :: A,B
      type(qstate), intent(out) :: C
      C%qu=B%u
      C%qv=B%v
      C%qz=B%z
      C%qt=B%t
      C%qs=B%s
      C%u=A%u
      C%v=A%v
      C%z=A%z
      C%t=A%t
      C%s=A%s
   end subroutine push_qstate

   subroutine pull_qstate(A,B,C)
      implicit none
      type(qstate), intent(in) :: C
      type(states), intent(out)  :: A,B
      B%u=C%qu
      B%v=C%qv
      B%z=C%qz
      B%t=C%qt
      B%s=C%qs
      A%u=C%u
      A%v=C%v
      A%z=C%z
      A%t=C%t
      A%s=C%s
   end subroutine pull_qstate

end module mod_mod_states
