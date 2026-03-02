!======================================================================
!  Program:   make_antime_list
!  Purpose:   Build a list of analysis times that have at least one
!             valid observation across a set of input observation files.
!
!  Usage:
!     ./make_antime_list [date1] [date2] [dt] [filelist]
!
!     date1    : initial date  (yyyy-mm-dd::HH:MM:SS)
!     date2    : final date    (yyyy-mm-dd::HH:MM:SS)
!     dt       : minimum time step between observations (seconds)
!     filelist : text file; each non-empty, non-comment line contains:
!                  <flag> <filename>
!                The <flag> token is read and ignored here, <filename>
!                is the path of an observation file to scan.
!
!  Output:
!     antime_list.txt  (formatted)
!       For every time with at least one valid observation:
!         line 1: yyyy-mm-dd::HH:MM:SS <nfiles>
!         line 2: a row of 0/1 flags (length = nfiles), where 1 marks
!                 that file contains a valid observation at that time.
!
!  Notes:
!    - All floating-point quantities use real64 (double precision).
!    - Time conversions are provided by module iso8601:
!        string2date, dts_to_abs_time, dts_from_abs_time, date2string.
!    - We compare times using integer rounding:
!        nint(t1 - t2) == 0
!      assuming absolute times are in seconds with integer resolution.
!    - Missing values use FLAG = -999.0_dp.
!
!  Updated:   2026-02-13
!======================================================================
program make_antime_list
  use iso_fortran_env, only : dp => real64
  use iso8601
  implicit none

  !-------------------------
  ! CLI arguments
  !-------------------------
  character(len=20) :: date1, date2
  character(len=80) :: filelist
  character(len=8)  :: dtstr

  !-------------------------
  ! Time control
  !-------------------------
  real(dp) :: dt                    ! time step in seconds
  integer      :: d, t, ierr
  real(dp) :: atime1, atime2, atime
  integer      :: n, nsteps, nn

  !-------------------------
  ! File list & per-file data
  !-------------------------
  integer                 :: k, ktot
  character(len=256)      :: line
  character(len=5)        :: oflag
  character(len=256)      :: ofile
  character(len=20)       :: dstr
  integer                 :: nr

  ! Work arrays for a single file
  real(dp), allocatable :: atval(:), val(:)

  ! Per-file/per-time storage
  integer,      allocatable :: nrec(:), isfile(:)
  real(dp), allocatable :: atval_tot(:,:), val_tot(:,:)

  ! Counters per time
  integer :: iok, iflag

  ! Constants
  real(dp), parameter :: FLAG = -999.0_dp

  !--------------------------------------------------------------------
  ! Read command-line arguments (no interactivity; fail fast on errors)
  !--------------------------------------------------------------------
  call get_command_argument(1, date1)
  call get_command_argument(2, date2)
  call get_command_argument(3, dtstr)
  call get_command_argument(4, filelist)

  if (trim(filelist) == '') then
    write(*,*) ''
    write(*,*) 'Make a list of times with at least one observation'
    write(*,*) ''
    write(*,*) 'Usage:'
    write(*,*) ''
    write(*,*) './make_antime_list [date1] [date2] [dt] [filelist]'
    write(*,*) ''
    write(*,*) 'date1: initial date (yyyy-mm-dd::HH:MM:SS)'
    write(*,*) 'date2: final date   (yyyy-mm-dd::HH:MM:SS)'
    write(*,*) 'dt:    minimum timestep between observations (sec)'
    write(*,*) 'filelist: list of observation files'
    write(*,*) ''
    stop
  end if

  ! date1
  call string2date(trim(date1), d, t, ierr)
  if (ierr /= 0) error stop 'Invalid date1'
  call dts_to_abs_time(d, t, atime1)

  ! date2
  call string2date(trim(date2), d, t, ierr)
  if (ierr /= 0) error stop 'Invalid date2'
  call dts_to_abs_time(d, t, atime2)

  ! time step
  read(dtstr, *, iostat=ierr) dt
  if (ierr /= 0) error stop 'Invalid timestep'
  if (dt <= 0.0_dp .or. dt > 864000.0_dp) error stop 'Bad timestep'
  if (atime2 - atime1 < 2.0_dp*dt) error stop 'Bad times: window too short'

  ! Number of steps (include both endpoints if aligned)
  nsteps = nint( (atime2 - atime1)/dt ) + 1

  !--------------------------------------------------------------------
  ! Read the file list robustly: skip blank and comment lines.
  ! Count valid entries first, then allocate and read again.
  !--------------------------------------------------------------------
  ktot = 0
  open(unit=20, file=trim(filelist), iostat=ierr, status='old', action='read')
  if (ierr /= 0) error stop 'Error opening the file list'

  do
    read(20, '(A)', iostat=ierr) line
    if (ierr < 0) exit                 ! EOF
    if (ierr > 0) cycle                ! read error: skip line
    if (len_trim(line) == 0) cycle     ! blank line
    if (line(1:1) == '#') cycle        ! comment

    ! Probe if the line has the two tokens we expect
    read(line, *, iostat=ierr) oflag, ofile
    if (ierr == 0 .and. len_trim(ofile) > 0) ktot = ktot + 1
  end do

  if (ktot <= 0) error stop 'Empty file list (no valid entries)'

  rewind(20)

  ! Allocate storage
  allocate(nrec(ktot), isfile(ktot))
  allocate(atval_tot(ktot, nsteps), val_tot(ktot, nsteps))
  allocate(atval(nsteps), val(nsteps))

  atval_tot = FLAG
  val_tot   = FLAG

  ! Second pass: ingest files
  k = 0
  do
    read(20, '(A)', iostat=ierr) line
    if (ierr < 0) exit
    if (ierr > 0) cycle
    if (len_trim(line) == 0) cycle
    if (line(1:1) == '#') cycle

    read(line, *, iostat=ierr) oflag, ofile
    if (ierr /= 0 .or. len_trim(ofile) == 0) cycle  ! skip malformed

    k = k + 1
    write(*,*) 'Reading: ', trim(ofile)

    ! Initialize work arrays to FLAG for safety
    atval = FLAG
    val   = FLAG

    call read_file_obs(trim(ofile), atime1, atime2, nsteps, nr, atval, val, dt)

    nrec(k)            = nr
    atval_tot(k, :)    = atval
    val_tot(k,  :)     = val

    if (k == ktot) exit
  end do
  close(20)

  if (k /= ktot) then
    write(*,*) 'Warning: counted files=', ktot, ' but read only k=', k
  end if

  !--------------------------------------------------------------------
  ! Build the output list: for each time, mark which files have data.
  !--------------------------------------------------------------------
  open(unit=30, file='antime_list.txt', form='formatted', status='replace', iostat=ierr)
  if (ierr /= 0) error stop 'Cannot open output file antime_list.txt'

  do n = 1, nsteps
    atime = atime1 + real(n - 1, dp) * dt

    isfile = 0
    iflag  = 0
    iok    = 0

    ! Loop on files
    do k = 1, ktot
      ! Loop on records stored for that file
      do nn = 1, nrec(k)
        ! Time match using integer rounding (sec resolution)
        if (nint(atime - atval_tot(k, nn)) == 0) then
          if (nint(val_tot(k, nn) - FLAG) == 0) then
            iflag = iflag + 1         ! encountered a flagged value
          else
            iok        = iok + 1
            isfile(k)  = 1            ! file k contributes at this time
          end if
        end if
      end do
    end do

    ! Output this time only if at least one valid obs is present
    if (iok > 0) then
      call dts_from_abs_time(d, t, atime)
      call date2string(d, t, dstr)
      write(30, '(a20,1x,i5)') dstr, ktot
      write(30, '(9999i2)')    isfile(1:ktot)
    end if
  end do

  close(30)

end program make_antime_list

!======================================================================
!  Subroutine: read_file_obs
!  Purpose   : Read an observation file and keep only those records whose
!              time falls exactly on the regular grid defined by:
!                atime1 + j*dt,  j = 0..nsteps-1
!              Matching is done with nint(atime - ctime) == 0 (sec grid).
!
!  Input:
!    ofile    : observation file path
!    atime1   : initial absolute time (sec)
!    atime2   : final   absolute time (sec)
!    nsteps   : number of time steps on the regular grid
!    dt       : time step (sec)
!
!  Output:
!    nrec     : number of kept records (<= nsteps)
!    atval(:) : absolute times of kept records; FLAG elsewhere
!    val(:)   : observed values aligned to the grid; FLAG elsewhere
!
!  File format (list-directed):
!    Each data line:  <datestr> <value>
!      where <datestr> is 'yyyy-mm-dd::HH:MM:SS'
!
!  Errors:
!    - Invalid or unreadable file -> error stop
!    - Invalid date in a line     -> error stop
!    - nrec > nsteps              -> error stop (should not happen)
!======================================================================
subroutine read_file_obs(ofile, atime1, atime2, nsteps, nrec, atval, val, dt)
  use iso_fortran_env, only : dp => real64
  use iso8601
  implicit none

  character(len=*), intent(in)    :: ofile
  integer,          intent(in)    :: nsteps
  real(dp),     intent(in)    :: atime1, atime2, dt
  integer,          intent(out)   :: nrec
  real(dp),     intent(inout) :: atval(nsteps), val(nsteps)

  integer          :: ierr
  character(len=20):: dstring
  integer          :: d, t
  real(dp)     :: atime, v, ctime
  integer          :: k, j

  open(unit=21, file=trim(ofile), status='old', action='read', iostat=ierr)
  if (ierr /= 0) error stop 'Error opening obs file: '//trim(ofile)

  k = 0

  do
    ! Read a line; tolerate blank lines by checking iostat and trimming
    read(21, *, iostat=ierr) dstring, v
    if (ierr < 0) exit               ! EOF
    if (ierr > 0) then
      ! Skip malformed/blank lines quietly
      cycle
    end if

    call string2date(trim(dstring), d, t, ierr)
    if (ierr /= 0) error stop 'Invalid date in file: '//trim(ofile)//' line='//trim(dstring)

    call dts_to_abs_time(d, t, atime)

    ! Keep only values that fall on the regular grid within 0.5 s (nint test)
    do j = 0, nsteps - 1
      ctime = atime1 + real(j, dp) * dt
      if (nint(atime - ctime) == 0) then
        k         = k + 1
        if (k > nsteps) error stop 'wrong nrec: exceeds nsteps for file '//trim(ofile)
        atval(k)  = atime
        val(k)    = v
        exit
      end if
    end do
  end do

  close(21)

  nrec = k

end subroutine read_file_obs
