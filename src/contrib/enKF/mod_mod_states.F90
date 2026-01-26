!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
!
! This routine determines the state of the model and some operations on it
!
module mod_mod_states

   integer,save :: nnkn = 0 !number of nodes
   integer,save :: nnel = 0 !number of elements

! standard model state
   type states
      real,allocatable :: u(:,:)                      ! 3-D u-velocity
      real,allocatable :: v(:,:)                      ! 3-D v-velocity
      real,allocatable :: z(:) 	                     ! 2-D water level 
      real,allocatable :: t(:,:)                      ! 3-D Temperature
      real,allocatable :: s(:,:)                      ! 3-D Salinity 
   end type states

! single precision model state (used for read and write to files)
   type states4
      real*4,allocatable :: u(:,:)                      ! 3-D u-velocity
      real*4,allocatable :: v(:,:)                      ! 3-D v-velocity
      real*4,allocatable :: z(:) 	                     ! 2-D water level 
      real*4,allocatable :: t(:,:)                      ! 3-D Temperature
      real*4,allocatable :: s(:,:)                      ! 3-D Salinity 
   end type states4

! double state, with the model errors
     type qstates
      real,allocatable :: qu(:,:)                       ! 3-D u-velocity error
      real,allocatable :: qv(:,:)                       ! 3-D v-velocity error
      real,allocatable :: qz(:) 	                     ! 2-D water level 
      real,allocatable :: qt(:,:)                       ! 3-D Temperature error
      real,allocatable :: qs(:,:)                       ! 3-D Salinity error

      real,allocatable :: u(:,:)                        ! 3-D u-velocity
      real,allocatable :: v(:,:)                        ! 3-D v-velocity
      real,allocatable :: z(:) 	                     ! 2-D water level 
      real,allocatable :: t(:,:)                        ! 3-D Temperature
      real,allocatable :: s(:,:)                        ! 3-D Salinity 
   end type qstates
 
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

   subroutine mod_mod_init(nk,ne)
    implicit none
    integer,intent(in) :: nk,ne

    nnkn = nk
    nnel = ne

   end subroutine mod_mod_init

   subroutine allocate_states(AA,nkn_dim,nel_dim,nlv_dim)
    implicit none
    type(states) :: AA
    integer :: nkn_dim,nel_dim,nlv_dim

       allocate(AA%u(nlv_dim,nel_dim))
       allocate(AA%v(nlv_dim,nel_dim))
       allocate(AA%z(nkn_dim))
       allocate(AA%t(nlv_dim,nkn_dim))
       allocate(AA%s(nlv_dim,nkn_dim))
       AA%u(:,:) = 0.
       AA%v(:,:) = 0.
       AA%z(:) = 0.
       AA%t(:,:) = 0.
       AA%s(:,:) = 0.

   end subroutine allocate_states

   subroutine allocate_states4(AA,nkn_dim,nel_dim,nlv_dim)
    implicit none
    type(states4) :: AA
    integer :: nkn_dim,nel_dim,nlv_dim

       allocate(AA%u(nlv_dim,nel_dim))
       allocate(AA%v(nlv_dim,nel_dim))
       allocate(AA%z(nkn_dim))
       allocate(AA%t(nlv_dim,nkn_dim))
       allocate(AA%s(nlv_dim,nkn_dim))
       AA%u(:,:) = 0.
       AA%v(:,:) = 0.
       AA%z(:) = 0.
       AA%t(:,:) = 0.
       AA%s(:,:) = 0.

   end subroutine allocate_states4

   subroutine allocate_qstates(AA,nkn_dim,nel_dim,nlv_dim)
    implicit none
    type(qstates) :: AA
    integer :: nkn_dim,nel_dim,nlv_dim

       allocate(AA%qu(nlv_dim,nel_dim))
       allocate(AA%qv(nlv_dim,nel_dim))
       allocate(AA%qz(nkn_dim))
       allocate(AA%qt(nlv_dim,nkn_dim))
       allocate(AA%qs(nlv_dim,nkn_dim))

       allocate(AA%u(nlv_dim,nel_dim))
       allocate(AA%v(nlv_dim,nel_dim))
       allocate(AA%z(nkn_dim))
       allocate(AA%t(nlv_dim,nkn_dim))
       allocate(AA%s(nlv_dim,nkn_dim))

       AA%qu(:,:) = 0.
       AA%qv(:,:) = 0.
       AA%qz(:) = 0.
       AA%qt(:,:) = 0.
       AA%qs(:,:) = 0.

       AA%u(:,:) = 0.
       AA%v(:,:) = 0.
       AA%z(:) = 0.
       AA%t(:,:) = 0.
       AA%s(:,:) = 0.

   end subroutine allocate_qstates

   subroutine deallocate_states(AA)
    implicit none
    type(states) :: AA

       deallocate(AA%u)
       deallocate(AA%v)
       deallocate(AA%z)
       deallocate(AA%t)
       deallocate(AA%s)

   end subroutine deallocate_states

   subroutine deallocate_states4(AA,nkn_dim,nel_dim,nlv_dim)
    implicit none
    type(states4) :: AA
    integer :: nkn_dim,nel_dim,nlv_dim

       deallocate(AA%u)
       deallocate(AA%v)
       deallocate(AA%z)
       deallocate(AA%t)
       deallocate(AA%s)

   end subroutine deallocate_states4

   subroutine deallocate_qstates(AA,nkn_dim,nel_dim,nlv_dim)
    implicit none
    type(qstates) :: AA
    integer :: nkn_dim,nel_dim,nlv_dim

       deallocate(AA%qu)
       deallocate(AA%qv)
       deallocate(AA%qz)
       deallocate(AA%qt)
       deallocate(AA%qs)

       deallocate(AA%u)
       deallocate(AA%v)
       deallocate(AA%z)
       deallocate(AA%t)
       deallocate(AA%s)

   end subroutine deallocate_qstates


!-------------------------------------------------------------------
! Functions
!-------------------------------------------------------------------

   function add_states(A,B)
      implicit none
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
      implicit none
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
      implicit none
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
      implicit none
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
      implicit none
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
      implicit none
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
      implicit none
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
      implicit none
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
      implicit none
      type(states) root_state
      type(states), intent(in) :: A
       root_state%u = sqrt(A%u)
       root_state%v = sqrt(A%v)
       root_state%z = sqrt(A%z)
       root_state%t = sqrt(A%t)
       root_state%s = sqrt(A%s)
   end function root_state

!-------------------------------------------------------------------
! subroutines
!-------------------------------------------------------------------

   subroutine assign_states(A,r)
      implicit none
      type(states), intent(out) :: A
      real, intent(in) :: r
       A%u = r
       A%v = r
       A%z = r
       A%t = r
       A%s = r
   end subroutine assign_states

   subroutine states4to8(A,B)
      implicit none
      type(states), intent(out) :: A
      type(states4), intent(in)  :: B
      A%u=DBLE(B%u)
      A%v=DBLE(B%v)
      A%z=DBLE(B%z)
      A%t=DBLE(B%t)
      A%s=DBLE(B%s)
   end subroutine states4to8

   subroutine states8to4(A,B)
      implicit none
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
      type(qstates), intent(out) :: C
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
      type(qstates), intent(in) :: C
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

   subroutine mean_state(n,A,B)
      implicit none
      integer, intent(in) :: n
      type(states), intent(in) :: A(n)
      type(states), intent(inout) :: B
      integer i

      if (n <= 0) error stop 'mean_state: invalid n'

      B%u(:,:) = 0.
      B%v(:,:) = 0.
      B%z(:) = 0.
      B%t(:,:) = 0.
      B%s(:,:) = 0.
      do i = 1,n
       B%u = B%u + A(i)%u
       B%v = B%v + A(i)%v
       B%z = B%z + A(i)%z
       B%t = B%t + A(i)%t
       B%s = B%s + A(i)%s
      end do
      B%u = B%u / float(n)
      B%v = B%v / float(n)
      B%z = B%z / float(n)
      B%t = B%t / float(n)
      B%s = B%s / float(n)
   end subroutine mean_state

   subroutine std_state(n,nk,ne,nl,A,B)
      implicit none
      integer, intent(in) :: n,nk,ne,nl
      type(states), intent(in) :: A(n)
      type(states), intent(inout) :: B
      type(states), allocatable :: Aaux
      integer i

      if (n <= 0) error stop 'mean_state: invalid n'

      allocate(Aaux)
      call allocate_states(Aaux,nk,ne,nl)

      ! This is not zero in case some vars are always zero
      Aaux%u(:,:) = 0.
      Aaux%v(:,:) = 0.
      Aaux%z(:) = 0.
      Aaux%t(:,:) = 0.
      Aaux%s(:,:) = 0.
      do i = 1,n
       Aaux%u = Aaux%u + A(i)%u
       Aaux%v = Aaux%v + A(i)%v
       Aaux%z = Aaux%z + A(i)%z
       Aaux%t = Aaux%t + A(i)%t
       Aaux%s = Aaux%s + A(i)%s
      end do
      Aaux%u = Aaux%u / float(n)
      Aaux%v = Aaux%v / float(n)
      Aaux%z = Aaux%z / float(n)
      Aaux%t = Aaux%t / float(n)
      Aaux%s = Aaux%s / float(n)

      B%u(:,:) = 0.
      B%v(:,:) = 0.
      B%z(:) = 0.
      B%t(:,:) = 0.
      B%s(:,:) = 0.
      do i = 1,n
       B%u = B%u + (A(i)%u - Aaux%u)**2
       B%v = B%v + (A(i)%v - Aaux%v)**2
       B%z = B%z + (A(i)%z - Aaux%z)**2
       B%t = B%t + (A(i)%t - Aaux%t)**2
       B%s = B%s + (A(i)%s - Aaux%s)**2
      end do
      B%u = sqrt(B%u / float(n-1))
      B%v = sqrt(B%v / float(n-1))
      B%z = sqrt(B%z / float(n-1))
      B%t = sqrt(B%t / float(n-1))
      B%s = sqrt(B%s / float(n-1))

      deallocate(Aaux)

   end subroutine std_state

   subroutine rtps_inflation(alpha,n,A,Amean,Astdo,Astdn)
      implicit none
   ! spatially variable rtps mult inflation

      real, intent(in) :: alpha
      integer, intent(in) :: n
      type(states), intent(inout) :: A(n)
      type(states), intent(in) :: Amean,Astdo,Astdn
      type(states), allocatable :: Aaux
      integer i

      allocate(Aaux)
      Aaux%u = (alpha * (Astdo%u - Astdn%u) / Astdn%u) + 1.
      Aaux%v = (alpha * (Astdo%v - Astdn%v) / Astdn%v) + 1.
      Aaux%z = (alpha * (Astdo%z - Astdn%z) / Astdn%z) + 1.
      Aaux%t = (alpha * (Astdo%t - Astdn%t) / Astdn%t) + 1.
      Aaux%s = (alpha * (Astdo%s - Astdn%s) / Astdn%s) + 1.

      do i = 1,n
         A(i)%u = Amean%u + (A(i)%u - Amean%u) * Aaux%u
         A(i)%v = Amean%v + (A(i)%v - Amean%v) * Aaux%v
         A(i)%z = Amean%z + (A(i)%z - Amean%z) * Aaux%z
         A(i)%t = Amean%t + (A(i)%t - Amean%t) * Aaux%t
         A(i)%s = Amean%s + (A(i)%s - Amean%s) * Aaux%s
      end do
      deallocate(Aaux)

   end subroutine rtps_inflation

   subroutine mult_inflation(alpha,n,A,Amean)
      implicit none
   ! spatially constant mult inflation

      real, intent(in) :: alpha
      integer, intent(in) :: n
      type(states), intent(inout) :: A(n)
      type(states), intent(in) :: Amean
      integer i

      do i = 1,n
         A(i)%u = Amean%u + (A(i)%u - Amean%u) * (alpha + 1.)
         A(i)%v = Amean%v + (A(i)%v - Amean%v) * (alpha + 1.)
         A(i)%z = Amean%z + (A(i)%z - Amean%z) * (alpha + 1.)
         A(i)%t = Amean%t + (A(i)%t - Amean%t) * (alpha + 1.)
         A(i)%s = Amean%s + (A(i)%s - Amean%s) * (alpha + 1.)
      end do

   end subroutine mult_inflation


end module mod_mod_states
