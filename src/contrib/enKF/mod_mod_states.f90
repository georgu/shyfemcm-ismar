module mod_mod_states

  use iso_fortran_env, only: dp => real64
  use mod_init_enkf,  only: nnlv, nnkn, nnel
  implicit none

  !==============================================================
  ! Type definitions
  !==============================================================
  type states
     ! NOTE: all components are allocatable to allow shape-conform operations
     real(dp), allocatable :: u(:,:)  ! (nnlv, nnel)
     real(dp), allocatable :: v(:,:)  ! (nnlv, nnel)
     real(dp), allocatable :: z(:)    ! (nnkn)
     real(dp), allocatable :: t(:,:)  ! (nnlv, nnkn)
     real(dp), allocatable :: s(:,:)  ! (nnlv, nnkn)
  end type states

  type states4
     ! Single-precision container (I/O or memory-saving scenarios)
     real,    allocatable :: u(:,:)
     real,    allocatable :: v(:,:)
     real,    allocatable :: z(:)
     real,    allocatable :: t(:,:)
     real,    allocatable :: s(:,:)
  end type states4

  type qstates
     ! Augmented container: physical fields + their error fields (q*)
     real(dp), allocatable :: u(:,:)  ! (nnlv, nnel)
     real(dp), allocatable :: v(:,:)  ! (nnlv, nnel)
     real(dp), allocatable :: z(:)    ! (nnkn)
     real(dp), allocatable :: t(:,:)  ! (nnlv, nnkn)
     real(dp), allocatable :: s(:,:)  ! (nnlv, nnkn)
     real(dp), allocatable :: qu(:,:) ! (nnlv, nnel)
     real(dp), allocatable :: qv(:,:) ! (nnlv, nnel)
     real(dp), allocatable :: qz(:)   ! (nnkn)
     real(dp), allocatable :: qt(:,:) ! (nnlv, nnkn)
     real(dp), allocatable :: qs(:,:) ! (nnlv, nnkn)
  end type qstates

  !==============================================================
  ! Public operator overloading for type(states)
  !==============================================================
  interface operator(+)
     module procedure add_states          ! states + states
     module procedure states_real_add     ! states + real(dp)
     module procedure real_states_add     ! real(dp) + states
  end interface

  interface operator(-)
     module procedure subtract_states     ! states - states
  end interface

  interface operator(*)
     module procedure states_real_mult    ! states * real(dp)
     module procedure real_states_mult    ! real(dp) * states
     module procedure states_states_mult  ! states * states (Hadamard/element-wise)
  end interface

  ! Defined assignment: states = real(dp)
  interface assignment(=)
     module procedure assign_states_real
  end interface

contains

!==============================================================
subroutine allocate_states(A, nk, ne, nl)
  ! Allocate and zero-initialize all components of a states object.
  type(states), intent(inout) :: A
  integer, intent(in) :: nk, ne, nl
  allocate(A%u(nl,ne))
  allocate(A%v(nl,ne))
  allocate(A%z(nk))
  allocate(A%t(nl,nk))
  allocate(A%s(nl,nk))
  A%u = 0.0_dp; A%v = 0.0_dp; A%z = 0.0_dp
  A%t = 0.0_dp; A%s = 0.0_dp
end subroutine allocate_states

!==============================================================
subroutine allocate_qstates(A, nk, ne, nl)
  ! Allocate and zero-initialize an augmented (q*) container.
  type(qstates), intent(inout) :: A
  integer, intent(in) :: nk, ne, nl
  allocate(A%u(nl,ne)); allocate(A%v(nl,ne)); allocate(A%z(nk))
  allocate(A%t(nl,nk)); allocate(A%s(nl,nk))
  allocate(A%qu(nl,ne)); allocate(A%qv(nl,ne)); allocate(A%qz(nk))
  allocate(A%qt(nl,nk)); allocate(A%qs(nl,nk))
  A%u = 0.0_dp; A%v = 0.0_dp; A%z = 0.0_dp
  A%t = 0.0_dp; A%s = 0.0_dp
  A%qu = 0.0_dp; A%qv = 0.0_dp; A%qz = 0.0_dp
  A%qt = 0.0_dp; A%qs = 0.0_dp
end subroutine allocate_qstates

!==============================================================
subroutine deallocate_qstates(A)
  ! Explicit deallocation of all components in qstates.
  type(qstates), intent(inout) :: A
  deallocate(A%u, A%v, A%z, A%t, A%s, A%qu, A%qv, A%qz, A%qt, A%qs)
end subroutine deallocate_qstates

!==============================================================
subroutine push_qstate(A,B,C)
  ! Pack (A,B) into an augmented container C = (A, q=B).
  implicit none
  type(states), intent(in)  :: A,B
  type(qstates), intent(out):: C
  C%qu=B%u; C%qv=B%v; C%qz=B%z; C%qt=B%t; C%qs=B%s
  C%u =A%u; C%v =A%v; C%z =A%z; C%t =A%t; C%s =A%s
end subroutine push_qstate

!==============================================================
subroutine pull_qstate(A,B,C)
  ! Unpack an augmented container: C -> (A,B), where B holds q* fields.
  implicit none
  type(qstates), intent(in)  :: C
  type(states),  intent(out) :: A,B
  B%u=C%qu; B%v=C%qv; B%z=C%qz; B%t=C%qt; B%s=C%qs
  A%u=C%u;  A%v=C%v;  A%z=C%z;  A%t=C%t;  A%s=C%s
end subroutine pull_qstate

!==============================================================
subroutine allocate_states4(A, nk, ne, nl)
  ! Single-precision allocation helper.
  type(states4), intent(inout) :: A
  integer, intent(in) :: nk, ne, nl
  allocate(A%u(nl,ne)); allocate(A%v(nl,ne)); allocate(A%z(nk))
  allocate(A%t(nl,nk)); allocate(A%s(nl,nk))
  A%u = 0.0; A%v = 0.0; A%z = 0.0
  A%t = 0.0; A%s = 0.0
end subroutine allocate_states4

!==============================================================
subroutine states8to4(A4, A)
  ! Convert double-precision states -> single-precision states4.
  type(states4), intent(inout) :: A4
  type(states),  intent(in)    :: A
  A4%u = real(A%u); A4%v = real(A%v); A4%z = real(A%z)
  A4%t = real(A%t); A4%s = real(A%s)
end subroutine states8to4

!==============================================================
subroutine states4to8(A, A4)
  ! Convert single-precision states4 -> double-precision states.
  type(states),  intent(inout) :: A
  type(states4), intent(in)    :: A4
  A%u = real(A4%u,dp); A%v = real(A4%v,dp); A%z = real(A4%z,dp)
  A%t = real(A4%t,dp); A%s = real(A4%s,dp)
end subroutine states4to8

!==============================================================
! states + states
function add_states(A,B) result(C)
  type(states)             :: C
  type(states), intent(in) :: A,B
  call allocate_states(C, nnkn, nnel, nnlv)
  C%u = A%u + B%u; C%v = A%v + B%v; C%z = A%z + B%z
  C%t = A%t + B%t; C%s = A%s + B%s
end function add_states

!==============================================================
! states - states
function subtract_states(A,B) result(C)
  type(states)             :: C
  type(states), intent(in) :: A,B
  call allocate_states(C, nnkn, nnel, nnlv)
  C%u = A%u - B%u; C%v = A%v - B%v; C%z = A%z - B%z
  C%t = A%t - B%t; C%s = A%s - B%s
end function subtract_states

!==============================================================
! real(dp) + states
function real_states_add(B, A) result(C)
  real(dp),     intent(in) :: B
  type(states), intent(in) :: A
  type(states)             :: C
  call allocate_states(C, nnkn, nnel, nnlv)
  C%u = B + A%u; C%v = B + A%v; C%z = B + A%z
  C%t = B + A%t; C%s = B + A%s
end function real_states_add

!==============================================================
! states + real(dp)
function states_real_add(A, B) result(C)
  type(states), intent(in) :: A
  real(dp),     intent(in) :: B
  type(states)             :: C
  call allocate_states(C, nnkn, nnel, nnlv)
  C%u = A%u + B; C%v = A%v + B; C%z = A%z + B
  C%t = A%t + B; C%s = A%s + B
end function states_real_add

!==============================================================
! real(dp) * states
function real_states_mult(B, A) result(C)
  real(dp),     intent(in) :: B
  type(states), intent(in) :: A
  type(states)             :: C
  call allocate_states(C, nnkn, nnel, nnlv)
  C%u = B * A%u; C%v = B * A%v; C%z = B * A%z
  C%t = B * A%t; C%s = B * A%s
end function real_states_mult

!==============================================================
! states * real(dp)
function states_real_mult(A, B) result(C)
  type(states), intent(in) :: A
  real(dp),     intent(in) :: B
  type(states)             :: C
  call allocate_states(C, nnkn, nnel, nnlv)
  C%u = A%u * B; C%v = A%v * B; C%z = A%z * B
  C%t = A%t * B; C%s = A%s * B
end function states_real_mult

!==============================================================
! states * states (Hadamard / element-wise product)
function states_states_mult(A,B) result(C)
  type(states)             :: C
  type(states), intent(in) :: A,B
  call allocate_states(C, nnkn, nnel, nnlv)
  C%u = A%u * B%u; C%v = A%v * B%v; C%z = A%z * B%z
  C%t = A%t * B%t; C%s = A%s * B%s
end function states_states_mult

!==============================================================
! Defined assignment: states = real(dp)
subroutine assign_states_real(A, B)
  ! Sets all components of A to scalar B (requires components to be allocatable).
  type(states), intent(inout) :: A
  real(dp),     intent(in)    :: B
  if (allocated(A%u)) A%u = B
  if (allocated(A%v)) A%v = B
  if (allocated(A%z)) A%z = B
  if (allocated(A%t)) A%t = B
  if (allocated(A%s)) A%s = B
end subroutine assign_states_real

!==============================================================
! Ensemble mean of type(states)
!   Amean = (1/nrens) * sum_{k=1}^{nrens} Aens(k)
! Requirements:
!   - Aens(:) length is nrens
!   - Operator overloading for states (+, *, assignment) available
!==============================================================
subroutine mean_state(Amean, Aens, nrens)
  implicit none
  type(states),               intent(out) :: Amean
  type(states), dimension(:), intent(in)  :: Aens
  integer,                    intent(in)  :: nrens
  integer :: k

  if (nrens <= 0) then
     ! Allocate and return zeros for safety
     call allocate_states(Amean, nnkn, nnel, nnlv)
     Amean = 0.0_dp
     return
  end if

  ! Allocate output and initialize to zero
  call allocate_states(Amean, nnkn, nnel, nnlv)
  Amean = 0.0_dp

  ! Accumulate sum
  do k = 1, nrens
     Amean = Amean + Aens(k)
  end do

  ! Scale by 1/n
  Amean = Amean * (1.0_dp / real(nrens, dp))
end subroutine mean_state

!==============================================================
! Ensemble standard deviation of type(states) (unbiased)
!   Astd = sqrt( sum_{k=1}^{nrens} (Aens(k) - Amean)^2 / (nrens - 1) )
! Inputs:
!   - Aens(:)  : ensemble members
!   - Amean    : ensemble mean (compute with mean_state)
! Notes:
!   - If nrens <= 1, returns zeros.
!   - Uses element-wise (Hadamard) operations via overloaded operators.
!==============================================================
subroutine std_state(Astd, Aens, Amean, nrens)
  implicit none
  type(states),               intent(out) :: Astd
  type(states), dimension(:), intent(in)  :: Aens
  type(states),               intent(in)  :: Amean
  integer,                    intent(in)  :: nrens
  integer :: k
  type(states) :: Atmp

  ! Guard for small ensembles
  call allocate_states(Astd, nnkn, nnel, nnlv)
  if (nrens <= 1) then
     Astd = 0.0_dp
     return
  end if

  ! Working buffer
  call allocate_states(Atmp, nnkn, nnel, nnlv)
  Astd = 0.0_dp

  ! Sum squared deviations
  do k = 1, nrens
     Atmp = Aens(k) - Amean     ! deviation
     Atmp = Atmp * Atmp         ! square (element-wise)
     Astd = Astd + Atmp
  end do

  ! Divide by (n-1)
  Astd = Astd * (1.0_dp / real(nrens - 1, dp))

  ! Element-wise square root
  Astd%u = sqrt(Astd%u)
  Astd%v = sqrt(Astd%v)
  Astd%z = sqrt(Astd%z)
  Astd%t = sqrt(Astd%t)
  Astd%s = sqrt(Astd%s)
end subroutine std_state


!==============================================================
end module mod_mod_states
