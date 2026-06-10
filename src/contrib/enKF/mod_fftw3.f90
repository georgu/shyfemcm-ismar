!===============================================================
!  Module: mod_fftw3
!  Purpose:
!    Provide named integer parameters (flags and directions)
!    compatible with FFTW3 usage from Fortran code.
!
!  Notes:
!    - This module only defines constants; no procedures are
!      required for double precision. Double precision control
!      happens in your arrays and FFTW plan/execution calls
!      (e.g., using double-precision real/complex data and the
!      proper FFTW interfaces).
!    - Constants are kept identical to your original file to
!      preserve compatibility with existing code.
!
!  Double precision guidance:
!    Use selected_real_kind(15, 307) or real64 kind from
!    ISO_FORTRAN_ENV for your data and plans elsewhere, e.g.:
!       use, intrinsic :: iso_fortran_env, only: real64
!       real(real64), allocatable :: x(:)
!    This module is agnostic to the floating-point kind.
!
!  Maintainer notes:
!    - Added IMPLICIT NONE and an explicit PUBLIC list to avoid
!      namespace leaks and accidental symbol shadowing.
!    - Added clear English comments for each flag.
!===============================================================
module mod_fftw3
  implicit none
  private

  !-----------------------------
  ! Transform type identifiers
  !-----------------------------
  integer, parameter, public :: FFTW_R2HC    = 0   ! Real-to-halfcomplex (legacy/utility transform)
  integer, parameter, public :: FFTW_HC2R    = 1   ! Halfcomplex-to-real
  integer, parameter, public :: FFTW_DHT     = 2   ! Discrete Hartley Transform

  !-----------------------------
  ! Real even/odd DCT/DST kinds
  ! (REDFT = DCT, RODFT = DST)
  !-----------------------------
  integer, parameter, public :: FFTW_REDFT00 = 3   ! DCT-I
  integer, parameter, public :: FFTW_REDFT01 = 4   ! DCT-II
  integer, parameter, public :: FFTW_REDFT10 = 5   ! DCT-III
  integer, parameter, public :: FFTW_REDFT11 = 6   ! DCT-IV

  integer, parameter, public :: FFTW_RODFT00 = 7   ! DST-I
  integer, parameter, public :: FFTW_RODFT01 = 8   ! DST-II
  integer, parameter, public :: FFTW_RODFT10 = 9   ! DST-III
  integer, parameter, public :: FFTW_RODFT11 = 10  ! DST-IV

  !-----------------------------
  ! Transform directions
  !-----------------------------
  integer, parameter, public :: FFTW_FORWARD  = -1 ! Forward transform sign
  integer, parameter, public :: FFTW_BACKWARD = +1 ! Backward (inverse) transform sign

  !-----------------------------
  ! Planner flags
  ! (combine with bitwise OR in C; in Fortran you typically pass
  !  a single integer containing the OR-combined flags as needed)
  !-----------------------------
  integer, parameter, public :: FFTW_MEASURE          = 0     ! Measure timings to find a good plan
  integer, parameter, public :: FFTW_DESTROY_INPUT    = 1     ! Allow planner to overwrite input
  integer, parameter, public :: FFTW_UNALIGNED        = 2     ! Permit unaligned arrays (less strict)
  integer, parameter, public :: FFTW_CONSERVE_MEMORY  = 4     ! Use less memory at possible speed cost
  integer, parameter, public :: FFTW_EXHAUSTIVE       = 8     ! Try many more plans (slow planning)
  integer, parameter, public :: FFTW_PRESERVE_INPUT   = 16    ! Do not overwrite input during planning
  integer, parameter, public :: FFTW_PATIENT          = 32    ! More thorough than MEASURE
  integer, parameter, public :: FFTW_ESTIMATE         = 64    ! No runtime measurements (fast planning)
  integer, parameter, public :: FFTW_ESTIMATE_PATIENT = 128   ! Estimate but more plan choices
  integer, parameter, public :: FFTW_BELIEVE_PCOST    = 256   ! Trust planner cost heuristic
  integer, parameter, public :: FFTW_DFT_R2HC_ICKY    = 512   ! Internal/legacy flag (rarely used)
  integer, parameter, public :: FFTW_NONTHREADED_ICKY = 1024  ! Internal/legacy flag (rarely used)
  integer, parameter, public :: FFTW_NO_BUFFERING     = 2048  ! Avoid extra buffering
  integer, parameter, public :: FFTW_NO_INDIRECT_OP   = 4096  ! Avoid certain algorithmic choices
  integer, parameter, public :: FFTW_ALLOW_LARGE_GENERIC = 8192 ! Permit generic codelets for large sizes
  integer, parameter, public :: FFTW_NO_RANK_SPLITS   = 16384 ! Disable rank split strategies
  integer, parameter, public :: FFTW_NO_VRANK_SPLITS  = 32768 ! Disable vector-rank split strategies
  integer, parameter, public :: FFTW_NO_VRECURSE      = 65536 ! Disable vector recursion strategies
  integer, parameter, public :: FFTW_NO_SIMD          = 131072! Disable SIMD usage in planner/kernels

end module mod_fftw3
