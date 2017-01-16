module mod_states
   use mod_dimensions
! standard model state
   type states
      real u(nnlv,nnel)                      ! 3-D u-velocity
      real v(nnlv,nnel)                      ! 3-D v-velocity
      !real z(nnkn)                           ! 2-D water level
      real ze(3,nnel)                        ! 2-D water level at vertices
      real t(nnlv,nnkn)                      ! 3-D Temperature
      real s(nnlv,nnkn)                      ! 3-D Salinity 
   end type states
   integer, save ::  global_ndim = 2*nnlv*nnel + 3*nnel + 2*nnlv*nnkn  ! Dimension of states
   !integer, save ::  global_ndim = 2*nnlv*nnel + nnkn + 3*nnel + 2*nnlv*nnkn  ! Dimension of states

! single precision model state (used for read and write to files)
   type states4
      real*4 u(nnlv,nnel)                      ! 3-D u-velocity
      real*4 v(nnlv,nnel)                      ! 3-D v-velocity
      !real*4 z(nnkn)                          ! 2-D water level
      real*4 ze(3,nnel)                        ! 2-D water level at vertices
      real*4 t(nnlv,nnkn)                      ! 3-D Temperature
      real*4 s(nnlv,nnkn)                      ! 3-D Salinity 
   end type states4

! model state at one grid point (used in local analysis)
   type sub_states
      real u(nnlv)                               ! 1-D u-velocity
      real v(nnlv)                               ! 1-D v-velocity
      !real z                                   ! 0-D water level
      real ze(3)                                   ! 0-D water level
      real t(nnlv)                               ! 1-D Temperature
      real s(nnlv)                               ! 1-D Salinity
   end type sub_states
   integer, parameter ::  local_ndim=4*nnlv + 1    ! Dimension of sub_states


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

!   interface operator(/)
!      module procedure divide_states
!   end interface

   interface assignment(=)
      module procedure assign_states
      module procedure states4to8
      module procedure states8to4
   end interface


contains
   type (sub_states) function getA(A,k,ie)
      implicit none
      type(states), intent(in)     :: A
      integer, intent(in) :: k,ie
      getA%u(:)=A%u(:,ie)
      getA%v(:)=A%v(:,ie)
      !getA%z=A%z(k)
      getA%ze=A%ze(:,ie)
      getA%t(:)=A%t(:,k)
      getA%s(:)=A%s(:,k)
   end function getA

   subroutine putA(subA,A,k,ie)
      implicit none
      type(sub_states), intent(in) :: subA
      type(states), intent(inout)  :: A
      integer, intent(in) :: k,ie
      A%u(:,ie)=subA%u(:)
      A%v(:,ie)=subA%v(:)
      !A%z(k)=subA%z
      A%ze(:,ie)=subA%ze
      A%t(:,k)=subA%t(:)
      A%s(:,k)=subA%s(:)
   end subroutine putA



   function add_states(A,B)
      type(states) add_states
      type(states), intent(in) :: A
      type(states), intent(in) :: B
       add_states%u = A%u + B%u
       add_states%v = A%v + B%v
       !add_states%z = A%z + B%z
       add_states%ze = A%ze + B%ze
       add_states%t = A%t + B%t
       add_states%s = A%s + B%s
   end function add_states

   function subtract_states(A,B)
      type(states) subtract_states
      type(states), intent(in) :: A
      type(states), intent(in) :: B
       subtract_states%u = A%u - B%u
       subtract_states%v = A%v - B%v
       !subtract_states%z = A%z - B%z
       subtract_states%ze = A%ze - B%ze
       subtract_states%t = A%t - B%t
       subtract_states%s = A%s - B%s
   end function subtract_states

   function states_real_mult(A,B)
      type(states) states_real_mult
      type(states), intent(in) :: A
      real, intent(in) :: B
       states_real_mult%u = B*A%u
       states_real_mult%v = B*A%v
       !states_real_mult%z = B*A%z
       states_real_mult%ze = B*A%ze
       states_real_mult%t = B*A%t
       states_real_mult%s = B*A%s
   end function states_real_mult

   function real_states_mult(B,A)
      type(states) real_states_mult
      type(states), intent(in) :: A
      real, intent(in) :: B
       real_states_mult%u = B*A%u
       real_states_mult%v = B*A%v
       !real_states_mult%z = B*A%z
       real_states_mult%ze = B*A%ze
       real_states_mult%t = B*A%t
       real_states_mult%s = B*A%s
   end function real_states_mult

   function states_states_mult(A,B)
      type(states) states_states_mult
      type(states), intent(in) :: A
      type(states), intent(in) :: B
       states_states_mult%u = A%u * B%u
       states_states_mult%v = A%v * B%v
       !states_states_mult%z = A%z * B%z
       states_states_mult%ze = A%ze * B%ze
       states_states_mult%t = A%t * B%t
       states_states_mult%s = A%s * B%s
   end function states_states_mult


   subroutine assign_states(A,r)
      type(states), intent(out) :: A
      real, intent(in) :: r
       A%u = r
       A%v = r
       !A%z = r
       A%ze = r
       A%t = r
       A%s = r
   end subroutine assign_states

   subroutine states4to8(A,B)
      type(states), intent(out) :: A
      type(states4), intent(in)  :: B
      A%u=DBLE(B%u)
      A%v=DBLE(B%v)
      !A%z=DBLE(B%z)
      A%ze=DBLE(B%ze)
      A%t=DBLE(B%t)
      A%s=DBLE(B%s)
   end subroutine states4to8

   subroutine states8to4(A,B)
      type(states), intent(in)  :: B
      type(states4),  intent(out) :: A
      A%u=real(B%u)
      A%v=real(B%v)
      !A%z=real(B%z)
      A%ze=real(B%ze)
      A%t=real(B%t)
      A%s=real(B%s)
   end subroutine states8to4

end module mod_states

