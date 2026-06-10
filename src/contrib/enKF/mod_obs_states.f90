!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights
! reserved.
!
!===============================================================
!  Module: mod_obs_states
!  Purpose:
!    Define lightweight data structures for observations used by
!    the EnKF data assimilation system:
!      - files      : entry in the observation file list
!      - scalar_0d  : 0-D point observations (e.g., tide gauges)
!      - vector_2d  : 2-D gridded vector fields (e.g., currents)
!
!  Status codes (kept as in your system):
!    0 = normal obs (assimilated)
!    1 = super-observation (assimilated)
!    2 = merged into a super-observation (not assimilated)
!    3 = out of range (not assimilated)
!    4 = flagged value (not assimilated)
!
!  Precision policy:
!    - Time fields (t) are double precision, matching your original
!      design and downstream time handling.
!    - Spatial coordinates and values remain REAL as in the original
!      to preserve binary I/O compatibility with existing FEM/data
!      readers and other modules.
!
!  Notes:
!    - No procedures are defined here; this module only contains
!      type declarations and is intentionally minimal.
!===============================================================
module mod_obs_states
  use iso_fortran_env, only : dp => real64
  implicit none
  private

  ! Make the types available to users of this module
  public :: files, scalar_0d, vector_2d

  !-------------------------------------------------------------
  ! Observation file entry
  !   ty   : type tag (e.g., '0DLEV', '0DTEM', '0DSAL', '2DVEL')
  !   name : path or basename of the file
  !-------------------------------------------------------------
  type files
     character(len=5)  :: ty    ! type of file
     character(len=80) :: name  ! name of the file (path or basename)
     real(dp)             :: x     ! x coordinate
     real(dp)             :: y     ! y coordinate
     real(dp)             :: z     ! z coordinate
     real(dp)             :: std   ! observation std
     real(dp)             :: rhol  ! radius for local analysis

  end type files

  !-------------------------------------------------------------
  ! 0-D scalar observation at a single location
  !-------------------------------------------------------------
  type scalar_0d
     real(dp)             :: t     ! time (absolute)
     real(dp)             :: x     ! x coordinate
     real(dp)             :: y     ! y coordinate
     real(dp)             :: z     ! z coordinate
     real(dp)             :: val   ! observed value
     real(dp)             :: std   ! observation std
     integer          :: stat  ! status = 0,1,2,3,4
     integer          :: id    ! id number of the source file
     real(dp)             :: rhol  ! radius for local analysis
  end type scalar_0d

  !-------------------------------------------------------------
  ! 2-D vector field observation on a regular grid
  !-------------------------------------------------------------
  type vector_2d
     real(dp)               :: t           ! time of the field (absolute)
     integer                :: nx, ny      ! dimensions
     real(dp), allocatable  :: x(:,:)      ! x coordinates
     real(dp), allocatable  :: y(:,:)      ! y coordinates
     real(dp)               :: z           ! depth (if applicable)
     real(dp), allocatable  :: u(:,:)      ! u component
     real(dp), allocatable  :: v(:,:)      ! v component
     real(dp), allocatable  :: std(:,:)    ! observation std
     integer, allocatable :: stat(:,:) ! status flags (0..4)
     integer            :: id          ! id number of the source file
  end type vector_2d

end module mod_obs_states
