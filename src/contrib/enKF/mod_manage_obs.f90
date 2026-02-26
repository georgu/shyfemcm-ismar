!======================================================================
!  Module: mod_manage_obs
!
!  Purpose:
!    Prepare and manage observations to be assimilated by the EnKF DA.
!    - Reads observation file list (obsfile) and detects available types.
!    - Loads 0D scalar time series (sea level, temperature, salinity).
!    - Loads 2D velocity fields from FEM files.
!    - Performs quality control and creates super-observations.
!
!  Precision policy:
!    - Time variables and time differences handled in double precision.
!    - Distances for grouping (super-obs) computed in double precision.
!    - Observation containers keep their original kinds (compatibility).
!
!  Status codes:
!    0 = normal obs (assimilated)
!    1 = super-observation (assimilated)
!    2 = merged into a super-observation (not-assimilated)
!    3 = out of range (not-assimilated)
!    4 = flagged value (not-assimilated)
!
!  Copyright:
!    (C) 2017, Marco Bajo, CNR-ISMAR Venice. All rights reserved.
!    Updated comments and corrections (2026-02-13).
!======================================================================
module mod_manage_obs
  use mod_para
  use mod_init_enkf
  use mod_obs_states
  implicit none

  type(files),      allocatable, private :: ofile(:)

  integer :: islev, isvel, istemp, issalt

  integer :: n_0dlev, n_1dlev, n_2dlev
  integer :: n_0dvel,           n_2dvel
  integer :: n_0dtemp, n_1dtemp, n_2dtemp
  integer :: n_0dsalt, n_1dsalt, n_2dsalt

  integer :: nobs_tot
  integer :: nobs_ok

  type(scalar_0d), allocatable :: o0dlev(:),  o0dtemp(:), o0dsalt(:)
  type(vector_2d), allocatable :: o2dvel(:)

contains

!======================================================================
!  read_obs
!  -------
!  Two-pass logic:
!    (1) scan to detect types & counts;
!    (2) allocate & store; build super-obs.
!======================================================================
  subroutine read_obs
    implicit none

    integer            :: ios
    character(len=80)  :: line
    integer            :: n, nfile
    integer            :: kinit, kend
    logical            :: linit
    integer            :: nobs
    double precision   :: tobs
    real               :: xobs, yobs, zobs
    real               :: vobs, stdobs, rho
    integer            :: statobs

    ! Reset type flags
    islev  = 0
    isvel  = 0
    istemp = 0
    issalt = 0

    ! Reset counters
    n_0dlev  = 0; n_1dlev  = 0; n_2dlev  = 0
    n_0dtemp = 0; n_1dtemp = 0; n_2dtemp = 0
    n_0dsalt = 0; n_1dsalt = 0; n_2dsalt = 0
    n_0dvel  = 0; n_2dvel  = 0

    write(*,*) 'Reading observations...'

    !-------------------------------
    ! Count lines in obsfile (no GOTO)
    !-------------------------------
    open(25, file=obsfile, status='old', form='formatted', iostat=ios)
    if (ios /= 0) error stop 'read_obs: error opening file list'

    nfile = 0
    do
      read(25, '(A)', iostat=ios) line
      if (ios < 0) exit                ! EOF
      if (ios > 0) error stop 'read_obs: error reading file list'
      nfile = nfile + 1
    end do

    rewind(25)
    if (nfile <= 0) error stop 'read_obs: empty file list'

    allocate(ofile(nfile))
    do n = 1, nfile
      read(25, *, iostat=ios) ofile(n)%ty, ofile(n)%name
      if (ios /= 0) error stop 'read_obs: malformed list entry'
    end do
    close(25)

    !-------------------------------
    ! First pass: detect types & counts
    !-------------------------------
    do n = 1, nfile
      select case (trim(ofile(n)%ty))

      case ('0DLEV')
        linit = .true.;  kinit = n_0dlev
        call read_scalar_0d('0DLEV', linit, trim(ofile(n)%name), TEPS, &
                            kinit, kend, tobs, xobs, yobs, zobs, vobs, stdobs, statobs, rho)
        if (kend > kinit) then
          n_0dlev = n_0dlev + 1
          islev   = 1
        end if

      case ('0DTEM')
        linit = .true.;  kinit = n_0dtemp
        call read_scalar_0d('0DTEM', linit, trim(ofile(n)%name), TEPS, &
                            kinit, kend, tobs, xobs, yobs, zobs, vobs, stdobs, statobs, rho)
        if (kend > kinit) then
          n_0dtemp = n_0dtemp + 1
          istemp   = 1
        end if

      case ('0DSAL')
        linit = .true.;  kinit = n_0dsalt
        call read_scalar_0d('0DSAL', linit, trim(ofile(n)%name), TEPS, &
                            kinit, kend, tobs, xobs, yobs, zobs, vobs, stdobs, statobs, rho)
        if (kend > kinit) then
          n_0dsalt = n_0dsalt + 1
          issalt   = 1
        end if

      case ('2DVEL')
        n_2dvel = n_2dvel + 1
        isvel   = 1

      case default
        error stop 'read_obs: Unknown file type'
      end select
    end do

    ! Enforce single-type assimilation per cycle (as per original note)
    if ((islev + isvel + istemp + issalt) > 1) then
      write(*,*) 'islev ',  islev
      write(*,*) 'isvel ',  isvel
      write(*,*) 'istemp ', istemp
      write(*,*) 'issalt ', issalt
      error stop 'Different observation types present. Assimilate them at different times.'
    end if

    !-------------------------------
    ! Allocate containers
    !-------------------------------
    if (n_0dlev  > 0) allocate(o0dlev(n_0dlev))
    if (n_0dtemp > 0) allocate(o0dtemp(n_0dtemp))
    if (n_0dsalt > 0) allocate(o0dsalt(n_0dsalt))
    if (n_2dvel  > 0) allocate(o2dvel(n_2dvel))

    !-------------------------------
    ! Second pass: read & store
    !-------------------------------
    n_2dvel  = 0
    kinit    = 0
    nobs_tot = 0

    do n = 1, nfile
      select case (trim(ofile(n)%ty))

      case ('0DLEV')
        linit = .false.
        call read_scalar_0d('0DLEV', linit, trim(ofile(n)%name), TEPS, &
                            kinit, kend, tobs, xobs, yobs, zobs, vobs, stdobs, statobs, rho)
        if (kend > kinit) then
          if (verbose) write(*,*) 'Station n. ', n
          o0dlev(kend)%t    = tobs
          o0dlev(kend)%x    = xobs
          o0dlev(kend)%y    = yobs
          o0dlev(kend)%z    = zobs
          o0dlev(kend)%val  = vobs
          o0dlev(kend)%std  = stdobs
          o0dlev(kend)%stat = statobs
          o0dlev(kend)%id   = n
          o0dlev(kend)%rhol = rho
          nobs_tot = nobs_tot + 1
        end if
        kinit = kend

      case ('0DTEM')
        linit = .false.
        call read_scalar_0d('0DTEM', linit, trim(ofile(n)%name), TEPS, &
                            kinit, kend, tobs, xobs, yobs, zobs, vobs, stdobs, statobs, rho)
        if (kend > kinit) then
          if (verbose) write(*,*) 'Station n. ', n
          o0dtemp(kend)%t    = tobs
          o0dtemp(kend)%x    = xobs
          o0dtemp(kend)%y    = yobs
          o0dtemp(kend)%z    = zobs
          o0dtemp(kend)%val  = vobs
          o0dtemp(kend)%std  = stdobs
          o0dtemp(kend)%stat = statobs
          o0dtemp(kend)%id   = n
          o0dtemp(kend)%rhol = rho
          nobs_tot = nobs_tot + 1
        end if
        kinit = kend

      case ('0DSAL')
        linit = .false.
        call read_scalar_0d('0DSAL', linit, trim(ofile(n)%name), TEPS, &
                            kinit, kend, tobs, xobs, yobs, zobs, vobs, stdobs, statobs, rho)
        if (kend > kinit) then
          if (verbose) write(*,*) 'Station n. ', n
          o0dsalt(kend)%t    = tobs
          o0dsalt(kend)%x    = xobs
          o0dsalt(kend)%y    = yobs
          o0dsalt(kend)%z    = zobs
          o0dsalt(kend)%val  = vobs
          o0dsalt(kend)%std  = stdobs
          o0dsalt(kend)%stat = statobs
          o0dsalt(kend)%id   = n
          o0dsalt(kend)%rhol = rho
          nobs_tot = nobs_tot + 1
        end if
        kinit = kend

      case ('2DVEL')
        n_2dvel = n_2dvel + 1
        call read_2dvel(trim(ofile(n)%name), n, n_2dvel, TEPS, nobs)
        nobs_tot = nobs_tot + 2 * nobs  ! u and v

      case default
        error stop 'read_obs: Unknown file type'
      end select
    end do

    ! Super-observations
    call make_super_2dvel
    call make_super_1dlev

    if (nobs_tot < 1) error stop 'No valid observations, stopping.'
  end subroutine read_obs

!======================================================================
!  read_scalar_0d
!  --------------
!======================================================================
  subroutine read_scalar_0d(olabel, linit, filin, eps, kinit, kend, atime_obs, xv, yv, zv, vv, stdvv, &
                            ostatusv, rho)
    use iso8601
    implicit none
    character(len=*), intent(in)  :: olabel
    logical,          intent(in)  :: linit
    character(len=*), intent(in)  :: filin
    double precision, intent(in)  :: eps
    integer,          intent(in)  :: kinit
    integer,          intent(out) :: kend
    double precision, intent(out) :: atime_obs
    real,             intent(out) :: xv, yv, zv, vv, stdvv, rho
    integer,          intent(out) :: ostatusv

    integer :: ios
    real    :: x, y, z, v, stdv
    integer :: ostatus
    character(len=80) :: dstring
    integer :: ierr
    integer :: date, time
    integer :: k
    logical :: stored

    ! Defaults
    xv = -999.0; yv = -999.0; zv = -999.0
    vv = -999.0; stdvv = -999.0
    ostatusv = -999

    ! Meta info (coordinates, std, rho)
    open(27, file=trim(filin)//'.info', status='old', form='formatted', iostat=ios)
    if (ios /= 0) error stop 'read_scalar_0d: error opening info file'
    read(27, *, iostat=ios) x, y, z, stdv, rho
    if (ios /= 0) error stop 'read_scalar_0d: malformed info file'
    close(27)

    ! Time series
    open(26, file=trim(filin), status='old', form='formatted', iostat=ios)
    if (ios /= 0) error stop 'read_scalar_0d: error opening series file'

    k = kinit
    stored = .false.

    do
      read(26, *, iostat=ios) dstring, v
      if (ios < 0) exit
      if (ios > 0) error stop 'read_scalar_0d: read error'

      if (isnan(v)) v = OFLAG
      call string2date(trim(dstring), date, time, ierr)
      if (ierr /= 0) error stop 'read_scalar_0d: bad date string'
      call dts_to_abs_time(date, time, atime_obs)

      if (abs(atime_obs - atime_an) < eps) then
        ostatus = 0
        call check_obs(olabel, v, v, OFLAG, ostatus)
        k = k + 1
        if (.not. linit .and. .not. stored) then
          xv = x; yv = y; zv = z
          vv = v; stdvv = stdv
          ostatusv = ostatus
          stored = .true.          ! store just the first valid obs
        end if
      end if
    end do

    close(26)
    kend = k
  end subroutine read_scalar_0d

!======================================================================
!  read_2dvel
!  ----------
!  Read a 2-D velocity field (u, v) from a FEM file at analysis time.
!  Memory-safe (temporary arrays deallocated in all paths).
!======================================================================
  subroutine read_2dvel(filin, fid, nrec, eps, nobs)
    implicit none
    character(len=*), intent(in)  :: filin
    integer,          intent(in)  :: fid
    integer,          intent(in)  :: nrec
    double precision, intent(in)  :: eps
    integer,          intent(out) :: nobs

    integer :: ios
    integer :: np, iformat, iunit
    integer :: irec, i, ii, jj
    integer :: nvers, lmax, nvar, ntype
    double precision :: tt
    integer :: datetime(2), ierr, nlvddi
    real*4, allocatable :: hhlv(:)
    real*4              :: regpar(7)
    integer :: nx, ny
    real    :: flag, dx, dy, x0, y0
    integer, allocatable :: ilhkv(:)
    real*4,  allocatable :: hd(:)
    real*4,  allocatable :: idata(:,:,:)
    character(len=50) :: string
    integer :: ostatus
    logical :: bdata
    real    :: uu, vv, ostd, rho
    integer :: ix, iy
    double precision :: atime_obs

    nobs  = 0
    bdata = .false.

    ! Standard deviation and rho from companion .info
    open(27, file=trim(filin)//'.info', status='old', form='formatted', iostat=ios)
    if (ios /= 0) error stop 'read_2dvel: error opening info file'
    read(27, *, iostat=ios) ostd, rho
    if (ios /= 0) error stop 'read_2dvel: malformed info file'
    close(27)

    ! Open FEM file
    np = 0
    write(*,*) 'Opening velocity FEM file: ', trim(filin)
    call fem_file_read_open(trim(filin), np, iformat, iunit)

    irec = 0
    do
      irec = irec + 1

      ! Headers
      call fem_file_read_params(iformat, iunit, tt, nvers, np, lmax, nvar, ntype, datetime, ierr)
      if (ierr < 0) exit

      call dts_convert_to_atime(datetime, tt, atime_obs)

      allocate(hhlv(lmax))
      nlvddi = lmax
      call fem_file_read_2header(iformat, iunit, ntype, lmax, hhlv, regpar, ierr)

      nx = nint(regpar(1)); ny = nint(regpar(2))
      x0 = regpar(3); y0 = regpar(4)
      dx = regpar(5); dy = regpar(6)
      flag = regpar(7)

      if (flag /= OFLAG) then
        deallocate(hhlv, stat=ios)
        error stop 'read_2dvel: bad flag value in header'
      end if

      ! If not the target time, skip data safely and continue
      if (abs(atime_obs - atime_an) > eps) then
        do i = 1, nvar
          call fem_file_skip_data(iformat, iunit, nvers, np, lmax, string, ierr)
          if (ierr /= 0) then
            deallocate(hhlv, stat=ios)
            error stop 'read_2dvel: error skipping data'
          end if
        end do
        deallocate(hhlv, stat=ios)
        cycle
      end if

      ! Matching time: proceed to read
      deallocate(hhlv, stat=ios)
      allocate(ilhkv(np), hd(np), idata(1, nx, ny))
      allocate(o2dvel(nrec)%x(nx,ny), o2dvel(nrec)%y(nx,ny), &
               o2dvel(nrec)%u(nx,ny), o2dvel(nrec)%v(nx,ny), &
               o2dvel(nrec)%std(nx,ny), o2dvel(nrec)%stat(nx,ny))

      do i = 1, nvar
        idata = flag
        call fem_file_read_data(iformat, iunit, nvers, np, lmax, string, ilhkv, hd, nlvddi, idata, ierr)
        if (ierr /= 0) then
          deallocate(ilhkv, hd, idata, stat=ios)
          error stop 'read_2dvel: error reading data'
        end if
        select case (i)
        case (1)
          o2dvel(nrec)%u = idata(1,:,:)
        case (2)
          o2dvel(nrec)%v = idata(1,:,:)
        case default
          deallocate(ilhkv, hd, idata, stat=ios)
          error stop 'read_2dvel: too many variables'
        end select
      end do
      deallocate(ilhkv, hd, idata, stat=ios)

      ! Coordinates
      do ii = 1, nx
        do jj = 1, ny
          o2dvel(nrec)%x(ii,jj) = x0 + dx * real(ii-1)
          o2dvel(nrec)%y(ii,jj) = y0 + dy * real(jj-1)
        end do
      end do
      o2dvel(nrec)%z   = 0.0
      o2dvel(nrec)%nx  = nx
      o2dvel(nrec)%ny  = ny
      o2dvel(nrec)%std = ostd   ! broadcast scalar std
      o2dvel(nrec)%id  = fid

      ! QC per point
      do ix = 1, nx
        do iy = 1, ny
          uu = o2dvel(nrec)%u(ix,iy)
          vv = o2dvel(nrec)%v(ix,iy)
          call check_obs('2DVEL', uu, vv, flag, ostatus)
          o2dvel(nrec)%stat(ix,iy) = ostatus
          if (ostatus == 0) then
            bdata = .true.
            nobs  = nobs + 1
          end if
        end do
      end do

      exit   ! first matching field is enough
    end do

    close(iunit)
    if (.not. bdata) error stop 'read_2dvel: 2DVEL file without valid data'
  end subroutine read_2dvel

!======================================================================
!  check_obs
!======================================================================
  subroutine check_obs(ty, v1, v2, flag, stat)
    implicit none
    character(len=*), intent(in)  :: ty
    real,            intent(in)   :: v1, v2, flag
    integer,         intent(out)  :: stat
    real :: vmin, vmax

    stat = 0

    ! Flag value check (keep original semantics)
    if (nint(v1 - flag) == 0 .or. nint(v2 - flag) == 0) then
      stat = 4
      return
    end if

    select case (trim(ty))
    case ('0DLEV'); vmin = SSH_MIN; vmax = SSH_MAX
    case ('0DTEM'); vmin = TEM_MIN; vmax = TEM_MAX
    case ('0DSAL'); vmin = SAL_MIN; vmax = SAL_MAX
    case ('2DVEL'); vmin = VEL_MIN; vmax = VEL_MAX
    case default
      write(*,*) 'Observation not implemented yet'
      error stop
    end select

    if (v1 <= vmin .or. v1 >= vmax) stat = 3
    if (v2 <= vmin .or. v2 >= vmax) stat = 3
  end subroutine check_obs

!======================================================================
!  make_super_1dlev
!  ----------------
!  Distance in double precision; unchanged interface.
!======================================================================
  subroutine make_super_1dlev
    implicit none
    integer, allocatable :: near_sts(:), near_sts_mat(:,:)
    integer, allocatable :: id_sorted(:), near_sts_sort(:)
    double precision     :: dist
    real                 :: x, y, vstd, val, vrhol
    real, parameter      :: mult_coeff = 0.25
    integer :: i, j, kk, nid
    integer :: knorm, ksup, kbad

    if (n_0dlev < 2) return

    allocate(near_sts(n_0dlev), near_sts_mat(n_0dlev,n_0dlev))
    allocate(near_sts_sort(n_0dlev), id_sorted(n_0dlev))

    write(*,*) 'Making superobs for the sea level...'
    near_sts     = 0
    near_sts_mat = 0

    do i = 1, n_0dlev
      do j = 1, n_0dlev
        dist = sqrt( dble(o0dlev(i)%x - o0dlev(j)%x)**2 + dble(o0dlev(i)%y - o0dlev(j)%y)**2 )
        if ( (dist < dble(o0dlev(i)%rhol)*dble(mult_coeff)) .or. &
             (dist < dble(o0dlev(j)%rhol)*dble(mult_coeff)) ) then
          near_sts(i)       = near_sts(i) + 1
          near_sts_mat(i,j) = 1
        end if
      end do
    end do

    call dsort(n_0dlev, near_sts, near_sts_sort, id_sorted)

    do i = 1, n_0dlev
      nid  = id_sorted(i)
      kk   = 0;  x = 0.0; y = 0.0; vstd = 0.0; vrhol = 0.0; val = 0.0

      do j = 1, n_0dlev
        if ((near_sts_mat(nid,j) == 1) .and. (o0dlev(j)%stat == 0)) then
          x     = x     + o0dlev(j)%x
          y     = y     + o0dlev(j)%y
          vstd  = vstd  + o0dlev(j)%std
          vrhol = vrhol + o0dlev(j)%rhol
          val   = val   + o0dlev(j)%val
          if (nid /= j) o0dlev(j)%stat = 2
          kk = kk + 1
        end if
      end do

      if (kk > 1) then
        write(*,*) 'Making sea-level superobservation...'
        o0dlev(nid)%x    = x    / kk
        o0dlev(nid)%y    = y    / kk
        o0dlev(nid)%std  = vstd / kk
        o0dlev(nid)%val  = val  / kk
        o0dlev(nid)%rhol = vrhol/ kk
        o0dlev(nid)%stat = 1
      end if
    end do

    knorm = 0; ksup = 0; kbad = 0
    write(*,*) 'Final sea-level observations:'
    do i = 1, n_0dlev
      write(*,'(a24,i1,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3)')  &
         'status,x,y,val,std,rho: ', o0dlev(i)%stat, o0dlev(i)%x, o0dlev(i)%y, &
         o0dlev(i)%val, o0dlev(i)%std, o0dlev(i)%rhol
      if (o0dlev(i)%stat == 0) knorm = knorm + 1
      if (o0dlev(i)%stat == 1) ksup  = ksup  + 1
      if (o0dlev(i)%stat >  1) kbad  = kbad  + 1
    end do
    write(*,*) 'Normal sea-level observations: ', knorm
    write(*,*) 'Super sea-level observations:  ', ksup
    write(*,*) 'Bad or merged sea-level obs:   ', kbad
  end subroutine make_super_1dlev

!======================================================================
!  dsort
!  -----
!  Sort integer array A (descending). 
!======================================================================
  subroutine dsort(ndim, A, B, IA)
    implicit none
    integer, intent(in)  :: ndim
    integer, intent(in)  :: A(ndim)
    integer, intent(out) :: B(ndim), IA(ndim)

    integer :: i, j, pos, best_val

    ! Initialize index array
    do i = 1, ndim
      IA(i) = i
    end do

    ! Simple selection sort on indices (descending by A)
    do i = 1, ndim
      pos = i
      best_val = A(IA(i))
      do j = i+1, ndim
        if (A(IA(j)) > best_val) then
          pos = j
          best_val = A(IA(j))
        end if
      end do
      if (pos /= i) then
        call swap_int(IA(i), IA(pos))
      end if
      B(i) = A(IA(i))
    end do
  contains
    subroutine swap_int(a, b)
      integer, intent(inout) :: a, b
      integer :: t
      t = a; a = b; b = t
    end subroutine swap_int
  end subroutine dsort

!======================================================================
!  make_super_2dvel
!  ----------------
!  Correct linear indexing (idx = (jj-1)*nx + ii).
!======================================================================
  subroutine make_super_2dvel
    implicit none
    integer :: nobs, nx, ny
    real,    allocatable :: x(:), y(:), v1(:), v2(:)
    integer, allocatable :: stat(:)
    integer :: n, ii, jj, idx

    if (n_2dvel < 1) return

    do n = 1, n_2dvel
      nx = o2dvel(n)%nx
      ny = o2dvel(n)%ny
      nobs = nx * ny

      allocate(x(nobs), y(nobs), v1(nobs), v2(nobs), stat(nobs))

      do jj = 1, ny
        do ii = 1, nx
          idx      = (jj-1)*nx + ii
          x(idx)   = o2dvel(n)%x(ii,jj)
          y(idx)   = o2dvel(n)%y(ii,jj)
          stat(idx)= o2dvel(n)%stat(ii,jj)
          v1(idx)  = o2dvel(n)%u(ii,jj)
          v2(idx)  = o2dvel(n)%v(ii,jj)
        end do
      end do

      call superobs_horiz_el(nobs, x, y, stat, v1, v2)

      do jj = 1, ny
        do ii = 1, nx
          idx = (jj-1)*nx + ii
          o2dvel(n)%stat(ii,jj) = stat(idx)
          o2dvel(n)%u(ii,jj)    = v1(idx)
          o2dvel(n)%v(ii,jj)    = v2(idx)
        end do
      end do

      deallocate(x, y, v1, v2, stat)
    end do
  end subroutine make_super_2dvel

!======================================================================
!  check_spread
!======================================================================
  subroutine check_spread(d, stdv, mval, mvalm)
    use levels
    use mod_ens_state
    implicit none
    real, intent(in)    :: d
    real, intent(in)    :: mvalm, mval(nrens)
    real, intent(inout) :: stdv
    integer :: ne
    integer, save :: icall = 0
    real :: ens_std, stdv_new

    if (icall == 0) write(*,*) 'Obs error variance limited to KSTD times the ensemble spread (Sakov 2012).'
    if (KSTD <= 0.0) return

    ens_std = 0.0
    do ne = 1, nrens
      ens_std = ens_std + (mval(ne) - mvalm)**2
    end do
    ens_std = sqrt( ens_std / real(nrens - 1) )

    if (abs(d) > KSTD * ens_std) then
      stdv_new = sqrt( sqrt( (ens_std**2 + stdv**2)**2 + ( (1.0/KSTD) * ens_std * d )**2 ) - ens_std**2 )
      write(*,'(a20,1x,3f8.3)') 'STDe,STDo,STDo_new: ', ens_std, stdv, stdv_new
      stdv = stdv_new
    end if
    icall = 1
  end subroutine check_spread

!======================================================================
!  check_spread_speed
!======================================================================
  subroutine check_spread_speed(d1, d2, stdv, um, vm, umm, vmm)
    use mod_ens_state
    implicit none
    real, intent(in)    :: d1, d2
    real, intent(inout) :: stdv
    real, intent(in)    :: um(nrens), vm(nrens), umm, vmm
    real :: cs(nrens), csm, inn
    integer :: ne
    real :: ens_std, stdv_new

    if (KSTD <= 0.0) return

    cs  = sqrt( um**2 + vm**2 )
    csm = sqrt( umm**2 + vmm**2 )

    ens_std = 0.0
    do ne = 1, nrens
      ens_std = ens_std + (cs(ne) - csm)**2
    end do
    ens_std = sqrt( ens_std / real(nrens - 1) )

    inn = sqrt( d1**2 + d2**2 )

    if (abs(inn) > KSTD * ens_std) then
      stdv_new = sqrt( sqrt( (ens_std**2 + stdv**2)**2 + ( (1.0/KSTD) * ens_std * inn )**2 ) - ens_std**2 )
      if (verbose) write(*,'(a18,2f8.4)') ' changing obs std ', stdv, stdv_new
      stdv = stdv_new
    end if
  end subroutine check_spread_speed

end module mod_manage_obs
