!======================================================================
!  Module: mod_init_enkf
!
!  Purpose:
!   Read all initial settings needed to initialize an EnKF cycle.
!   The information is taken from the file 'analysis.info'.
!
!
!  Copyright:
!    (C) 2017, Marco Bajo, CNR-ISMAR Venice. All rights reserved.
!    Updated comments and corrections (2026-02-13).
!======================================================================
module mod_init_enkf
  implicit none

  !-----------------------------
  ! Control parameters read from analysis.info
  !-----------------------------
  integer            :: nnkn       ! number of nodes
  integer            :: nnel       ! number of elements
  integer            :: nnlv       ! number of vertical levels
  integer            :: nrens      ! number of ensemble members
  integer            :: nanal      ! analysis step index
  character(len=80)  :: basfile    ! name of the basin file (no extension)
  character(len=80)  :: obsfile    ! name of the observation file list
  character(len=80)  :: ostring    ! date/time string for observation time
  double precision   :: atime_an   ! analysis time in absolute units

contains
!======================================================================
!  Subroutine: read_info
!
!  Purpose:
!    Read all configuration parameters from 'analysis.info',
!    check their validity, and convert time strings to absolute time.
!
!  Behavior:
!    - Reads parameters in strict order.
!    - Performs validity checks for rmode, is_local, and nrens.
!    - Converts the observation time string into absolute time.
!======================================================================
  subroutine read_info
    use iso8601                ! for string2date and absolute time conversion
    use mod_para               ! provides rmode and is_local
    implicit none

    integer :: ierr
    integer :: date, time
    integer :: ios

    !------------------------------------------------------------
    ! Open and read configuration file
    !------------------------------------------------------------
    open(20, file='analysis.info', status='old', action='read', iostat=ios)
    if (ios /= 0) then
       error stop 'read_info: cannot open analysis.info'
    end if

    read(20, *, iostat=ios) nnkn
    if (ios /= 0) error stop 'read_info: error reading nnkn'

    read(20, *, iostat=ios) nnel
    if (ios /= 0) error stop 'read_info: error reading nnel'

    read(20, *, iostat=ios) nnlv
    if (ios /= 0) error stop 'read_info: error reading nnlv'

    read(20, *, iostat=ios) nrens
    if (ios /= 0) error stop 'read_info: error reading nrens'

    read(20, *, iostat=ios) nanal
    if (ios /= 0) error stop 'read_info: error reading nanal'

    read(20, *, iostat=ios) basfile
    if (ios /= 0) error stop 'read_info: error reading basfile'

    read(20, *, iostat=ios) ostring
    if (ios /= 0) error stop 'read_info: error reading ostring'

    read(20, *, iostat=ios) obsfile
    if (ios /= 0) error stop 'read_info: error reading obsfile'

    ! rmode and is_local come from mod_para
    read(20, *, iostat=ios) rmode
    if (ios /= 0) error stop 'read_info: error reading rmode'

    read(20, *, iostat=ios) is_local
    if (ios /= 0) error stop 'read_info: error reading is_local'

    close(20)

    !------------------------------------------------------------
    ! Validity checks
    !------------------------------------------------------------

    ! Check allowed analysis modes
    if ( (rmode /= 11) .and. (rmode /= 12) .and. (rmode /= 13) &
         .and. (rmode /= 21) .and. (rmode /= 22) .and. (rmode /= 23) ) then
       error stop 'Wrong analysis method.'
    end if

    ! Localisation option must be 0 or 1
    if ( (is_local /= 0) .and. (is_local /= 1) ) then
       error stop 'Wrong localisation option.'
    end if

    ! EnKF constraint: nrens must be odd (control member first)
    if (mod(nrens, 2) == 0) then
       error stop 'read_info: number of ensemble members must be odd.'
    end if

    !------------------------------------------------------------
    ! Convert date string to absolute time
    !------------------------------------------------------------
    call string2date(trim(ostring), date, time, ierr)
    if (ierr /= 0) then
       error stop 'read_info: invalid date string'
    end if

    call dts_to_abs_time(date, time, atime_an)

    !------------------------------------------------------------
    ! User feedback
    !------------------------------------------------------------
    write(*,*) 'Time of the analysis step: ', trim(ostring)
    write(*,*) 'Number of ensemble members: ', nrens

  end subroutine read_info

end module mod_init_enkf
