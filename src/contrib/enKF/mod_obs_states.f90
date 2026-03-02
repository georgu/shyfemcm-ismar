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
  end type files

  !-------------------------------------------------------------
  ! 0-D scalar observation at a single location
  !   t     : observation time (absolute), double precision
  !   x,y,z : coordinates (model grid / projected units)
  !   val   : observed value
  !   std   : observation standard deviation
  !   stat  : status flag (0..4, see header)
  !   id    : originating file index (as read from list)
  !   rhol  : localisation radius (if localisation is used)
  !-------------------------------------------------------------
  type scalar_0d
     double precision :: t     ! time (absolute)
     real             :: x     ! x coordinate
     real             :: y     ! y coordinate
     real             :: z     ! z coordinate
     real             :: val   ! observed value
     real             :: std   ! observation std
     integer          :: stat  ! status = 0,1,2,3,4
     integer          :: id    ! id number of the source file
     real             :: rhol  ! radius for local analysis
  end type scalar_0d

  !-------------------------------------------------------------
  ! 2-D vector field observation on a regular grid
  !   t        : field time (absolute), double precision
  !   nx, ny   : grid dimensions
  !   x, y     : 2-D arrays of grid coordinates (nx, ny)
  !   z        : representative depth/level (scalar)
  !   u, v     : vector components on the grid (nx, ny)
  !   std      : per-point std on the grid (nx, ny)
  !   stat     : per-point status flag (nx, ny)
  !   id       : originating file index (as read from list)
  !-------------------------------------------------------------
  type vector_2d
     double precision   :: t           ! time of the field (absolute)
     integer            :: nx, ny      ! dimensions
     real, allocatable  :: x(:,:)      ! x coordinates
     real, allocatable  :: y(:,:)      ! y coordinates
     real               :: z           ! depth (if applicable)
     real, allocatable  :: u(:,:)      ! u component
     real, allocatable  :: v(:,:)      ! v component
     real, allocatable  :: std(:,:)    ! observation std
     integer, allocatable :: stat(:,:) ! status flags (0..4)
     integer            :: id          ! id number of the source file
  end type vector_2d

end module mod_obs_states
