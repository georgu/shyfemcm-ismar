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
  use iso_fortran_env, only : dp => real64
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
    real(dp)   :: tobs
    real(dp)               :: xobs, yobs, zobs
    real(dp)               :: vobs, stdobs, rho, vmean
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
    ! Count lines in obsfile
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
                            kinit, kend, tobs, xobs, yobs, zobs, vobs, stdobs, vmean, statobs, rho)
        if (kend > kinit) then
          n_0dlev = n_0dlev + 1
          islev   = 1
        end if

      case ('0DTEM')
        linit = .true.;  kinit = n_0dtemp
        call read_scalar_0d('0DTEM', linit, trim(ofile(n)%name), TEPS, &
                            kinit, kend, tobs, xobs, yobs, zobs, vobs, stdobs, vmean, statobs, rho)
        if (kend > kinit) then
          n_0dtemp = n_0dtemp + 1
          istemp   = 1
        end if

      case ('0DSAL')
        linit = .true.;  kinit = n_0dsalt
        call read_scalar_0d('0DSAL', linit, trim(ofile(n)%name), TEPS, &
                            kinit, kend, tobs, xobs, yobs, zobs, vobs, stdobs, vmean, statobs, rho)
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

    ! Enforce single-type assimilation per cycle
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
                            kinit, kend, tobs, xobs, yobs, zobs, vobs, stdobs, vmean, statobs, rho)
        if (kend > kinit) then
          call check_mean('0DLEV',vmean,zobs,vobs)
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
                            kinit, kend, tobs, xobs, yobs, zobs, vobs, stdobs, vmean, statobs, rho)
        if (kend > kinit) then
	  ! not very correct for temp
          !call check_mean('0DTEM',vmean,zobs,vobs)
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
                            kinit, kend, tobs, xobs, yobs, zobs, vobs, stdobs, vmean, statobs, rho)
        if (kend > kinit) then
	  ! not very correct for salt
          !call check_mean('0DSAL',vmean,zobs,vobs)
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
    if ( SUPEROBS ) call make_super_1dlev
    if ( SUPEROBS ) call make_super_2dvel

    if (nobs_tot < 1) error stop 'No valid observations, stopping.'
  end subroutine read_obs

!======================================================================
!  read_scalar_0d
!  --------------
!======================================================================
  subroutine read_scalar_0d(olabel, linit, filin, eps, kinit, kend, atime_obs, xv, yv, zv, vv, stdvv, &
                            vmean, ostatusv, rho)
    use iso8601
    implicit none
    character(len=*), intent(in)  :: olabel
    logical,          intent(in)  :: linit
    character(len=*), intent(in)  :: filin
    real(dp),         intent(in)  :: eps
    integer,          intent(in)  :: kinit
    integer,          intent(out) :: kend
    real(dp),         intent(out) :: atime_obs
    real(dp),         intent(out) :: xv, yv, zv, vv, stdvv, rho
    integer,          intent(out) :: ostatusv
    real(dp),         intent(out) :: vmean

    integer :: ios
    real(dp)    :: x, y, z, v, stdv
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

    vmean = 0.
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
        vmean = vmean + v
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
    vmean = vmean / kend
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
    real(dp), intent(in)  :: eps
    integer,          intent(out) :: nobs

    integer :: ios
    integer :: np, iformat, iunit
    integer :: irec, i, ii, jj
    integer :: nvers, lmax, nvar, ntype
    real(dp) :: tt
    integer :: datetime(2), ierr, nlvddi
    real*4, allocatable :: hhlv(:)
    real*4              :: regpar(7)
    integer :: nx, ny
    real(dp)    :: flag, dx, dy, x0, y0
    integer, allocatable :: ilhkv(:)
    real*4,  allocatable :: hd(:)
    real*4,  allocatable :: idata(:,:,:)
    character(len=50) :: string
    integer :: ostatus
    logical :: bdata
    real(dp)    :: uu, vv, ostd, rho
    integer :: ix, iy
    real(dp) :: atime_obs

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
          o2dvel(nrec)%x(ii,jj) = x0 + dx * real(ii-1, dp)
          o2dvel(nrec)%y(ii,jj) = y0 + dy * real(jj-1, dp)
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
    real(dp),            intent(in)   :: v1, v2, flag
    integer,         intent(out)  :: stat
    real(dp) :: vmin, vmax

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
!  Create superobservations for sea-level scalar observations.
!  Distances are computed in meters using deg2meters.
!======================================================================
subroutine make_super_1dlev
  implicit none
  integer, allocatable :: near_sts(:), near_sts_mat(:,:)
  integer, allocatable :: id_sorted(:), near_sts_sort(:)
  real(dp) :: dist_mx, dist_my, rho_m1, rho_m2, dist_m
  real(dp) :: x, y, vstd, val, vrhol
  real(dp), parameter :: mult_coeff = 1.0_dp
  integer :: i, j, nid, kk
  integer :: knorm, ksup, kbad

  if (n_0dlev < 2) return

  ! 1. Allocation
  allocate(near_sts(n_0dlev), near_sts_mat(n_0dlev,n_0dlev))
  allocate(near_sts_sort(n_0dlev), id_sorted(n_0dlev))

  near_sts     = 0
  near_sts_mat = 0

  ! 2. Identify neighbors based on rhol (horizontal correlation length)
  do i = 1, n_0dlev
    do j = 1, n_0dlev
      if (i == j) cycle

      ! Calculate distances in meters
      call deg2meters(o0dlev(i)%x, o0dlev(i)%y, o0dlev(i)%x - o0dlev(j)%x, .true.,  dist_mx)
      call deg2meters(o0dlev(i)%x, o0dlev(i)%y, o0dlev(i)%y - o0dlev(j)%y, .false., dist_my)
      
      ! Convert correlation lengths from degrees to meters
      call deg2meters(o0dlev(i)%x, o0dlev(i)%y, o0dlev(i)%rhol, .false., rho_m1)
      call deg2meters(o0dlev(j)%x, o0dlev(j)%y, o0dlev(j)%rhol, .false., rho_m2)

      dist_m = sqrt(dist_mx**2 + dist_my**2)
      
      if ( (dist_m < rho_m1 * mult_coeff) .or. (dist_m < rho_m2 * mult_coeff) ) then
        near_sts(i)       = near_sts(i) + 1
        near_sts_mat(i,j) = 1
      end if
    end do
  end do

  ! 3. Sort observations by number of neighbors to process dense areas first
  call dsort(n_0dlev, near_sts, near_sts_sort, id_sorted)

  ! 4. Merge observations
  do i = 1, n_0dlev
    nid = id_sorted(i)
    
    ! Skip if this observation has already been merged into another (stat=2)
    if (o0dlev(nid)%stat /= 0) cycle

    ! Initialize super-observation sums with the current "pivot" (nid)
    kk = 1
    x     = o0dlev(nid)%x
    y     = o0dlev(nid)%y
    vstd  = o0dlev(nid)%std
    val   = o0dlev(nid)%val
    vrhol = o0dlev(nid)%rhol

    do j = 1, n_0dlev
      if (j == nid) cycle
      
      ! If j is a neighbor and still available (stat=0)
      if ((near_sts_mat(nid,j) == 1) .and. (o0dlev(j)%stat == 0)) then
        x     = x     + o0dlev(j)%x
        y     = y     + o0dlev(j)%y
        vstd  = vstd  + o0dlev(j)%std
        val   = val   + o0dlev(j)%val
        vrhol = vrhol + o0dlev(j)%rhol
        
        o0dlev(j)%stat = 2 ! Mark as merged/inactive
        kk = kk + 1
      end if
    end do

    ! If neighbors were found, update the pivot and set it as a super-observation (stat=1)
    if (kk > 1) then
      o0dlev(nid)%x    = x / real(kk, dp)
      o0dlev(nid)%y    = y / real(kk, dp)
      o0dlev(nid)%std  = vstd / real(kk, dp)
      o0dlev(nid)%val  = val / real(kk, dp)
      o0dlev(nid)%rhol = vrhol / real(kk, dp)
      o0dlev(nid)%stat = 1
    end if
  end do

  ! 5. Final count and cleanup
  knorm = 0; ksup = 0; kbad = 0
  do i = 1, n_0dlev
    if (o0dlev(i)%stat == 0) knorm = knorm + 1
    if (o0dlev(i)%stat == 1) ksup  = ksup  + 1
    if (o0dlev(i)%stat >  1) kbad  = kbad  + 1
  end do
  
  write(*,*) 'Normal sea-level observations: ', knorm
  write(*,*) 'Super sea-level observations:  ', ksup
  write(*,*) 'Bad or merged sea-level obs:   ', kbad

  deallocate(near_sts, near_sts_mat, near_sts_sort, id_sorted)

end subroutine make_super_1dlev

!======================================================================
!  make_super_2dvel
!  ----------------
!  Flatten 2D velocity fields and send the observations to
!  superobs_horiz_el, which performs FEM-element horizontal clustering.
!======================================================================
subroutine make_super_2dvel
  implicit none
  integer :: nobs, nx, ny
  real(dp), allocatable :: x(:), y(:), u(:), v(:)
  integer,  allocatable :: stat(:)
  integer :: n

  if (n_2dvel < 1) return

  write(*,*) 'Processing 2D velocity super-observations...'

  do n = 1, n_2dvel
    nx = o2dvel(n)%nx
    ny = o2dvel(n)%ny
    nobs = nx * ny

    ! Allocate temporary 1D buffers
    allocate(x(nobs), y(nobs), u(nobs), v(nobs), stat(nobs))

    ! Flatten 2D fields into 1D vectors using Fortran intrinsic RESHAPE/PACK 
    ! or simple array assignment (column-major order)
    x    = reshape(o2dvel(n)%x,    [nobs])
    y    = reshape(o2dvel(n)%y,    [nobs])
    u    = reshape(o2dvel(n)%u,    [nobs])
    v    = reshape(o2dvel(n)%v,    [nobs])
    stat = reshape(o2dvel(n)%stat, [nobs])

    ! Perform super-observation logic
    call superobs_horiz_el(nobs, x, y, stat, u, v)

    ! Copy back 1D results into the original 2D structure
    o2dvel(n)%u    = reshape(u,    [nx, ny])
    o2dvel(n)%v    = reshape(v,    [nx, ny])
    o2dvel(n)%stat = reshape(stat, [nx, ny])

    ! Explicit deallocation inside the loop to prevent memory accumulation
    deallocate(x, y, u, v, stat)
  end do

end subroutine make_super_2dvel

!----------------------------------------------------------------------
!  superobs_horiz_el  (versione compatibile con find_element in real*4)
!----------------------------------------------------------------------
subroutine superobs_horiz_el(no, x, y, ostatus, val1, val2)
  use basin
  implicit none

  integer,  intent(in)    :: no
  real(dp), intent(in)    :: x(no), y(no)
  integer,  intent(inout) :: ostatus(no)
  real(dp), intent(inout) :: val1(no), val2(no)

  integer :: n, ie, nn, omax
  integer, allocatable :: ieobs(:)
  integer, allocatable :: oindex(:,:)
  integer, allocatable :: nobs(:)     ! Changed to allocatable to handle memory better
  real(dp) :: av1, av2
  real(4)  :: x4, y4                  ! Consistent with find_element requirements

  if (no <= 0 .or. nel <= 0) return

  ! Use allocate for nobs if nel is large to avoid stack overflow
  allocate(nobs(nel), ieobs(no))
  ieobs = -999
  nobs  = 0

  ! 1. Identify containing FEM element for each observation
  do n = 1, no
    if (ostatus(n) /= 0) cycle ! Process only 'normal' observations

    x4 = real(x(n), 4)
    y4 = real(y(n), 4)

    call find_element(x4, y4, ie)

    if (ie >= 1 .and. ie <= nel) then
       ieobs(n) = ie
       nobs(ie) = nobs(ie) + 1
    end if
  end do

  ! 2. Determine maximum occupancy and handle empty elements
  omax = maxval(nobs)
  if (omax <= 0) then
    deallocate(ieobs, nobs)
    return
  end if

  ! 3. Build observation index per element
  allocate(oindex(0:omax, nel))
  oindex = 0

  do n = 1, no
    ie = ieobs(n)
    if (ie < 1 .or. ie > nel) cycle

    oindex(0, ie) = oindex(0, ie) + 1
    nn = oindex(0, ie)
    oindex(nn, ie) = n
  end do

  ! 4. Build super-observations (average per element)
  do ie = 1, nel
    nn = oindex(0, ie)
    if (nn <= 0) cycle

    ! Compute average using vector subscripting
    av1 = sum(val1(oindex(1:nn, ie))) / real(nn, dp)
    av2 = sum(val2(oindex(1:nn, ie))) / real(nn, dp)

    ! Update first observation in element as the 'Super-Observation'
    val1(oindex(1, ie)) = av1
    val2(oindex(1, ie)) = av2
    ostatus(oindex(1, ie)) = 1

    ! Mark remaining observations as merged (stat=2)
    if (nn > 1) ostatus(oindex(2:nn, ie)) = 2
  end do

  deallocate(ieobs, oindex, nobs)
end subroutine superobs_horiz_el

!======================================================================
!  dsort
!  -----
!  Sort integer array A in descending order.
!  Output:
!     IA(i) = index of i-th sorted element of A
!     B(i)  = sorted values (B = A(IA))
!======================================================================
subroutine dsort(ndim, A, B, IA)
  implicit none
  integer, intent(in)  :: ndim
  integer, intent(in)  :: A(ndim)
  integer, intent(out) :: B(ndim), IA(ndim)

  integer :: i, j, pos, best_val

  if (ndim <= 0) return

  do i = 1, ndim
    IA(i) = i
  end do

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
! Check the mean of the observations with respect to the model
!======================================================================
subroutine check_mean(obs_type, vmean, z, v)
  use mod_ens_state
  implicit none
  
  ! Argument declarations
  character(len=*), intent(in) :: obs_type  ! Type of observation string
  real(dp), intent(in)         :: vmean     ! Observed value
  real(dp), intent(in)         :: z         ! Depth/vertical coordinate
  real(dp), intent(inout)      :: v         ! Value to be corrected
  
  ! Local variables
  real(dp) :: mmean
  integer  :: idx_near, i, k

  if ( NOBIAS ) return

  if (.not. allocated(Abk)) then
     write(*,*) 'ERROR: Abk is not allocated in check_mean!'
     return
  end if

  select case (trim(obs_type))

  case ('0DLEV')
     ! Calculate ensemble mean for sea level (z-component)
     mmean = 0._dp
     do i = 1, nnkn
        mmean = mmean + Abk(1)%z(i)
     end do
     mmean = mmean / real(nnkn, dp)

     if ( abs(vmean - mmean) > ZBIAS_MAX ) then
        write(*,*) 'Warning: bias in observations, correcting (mod/obs): ', mmean, vmean
        v = v - vmean + mmean
     end if

  case ('0DTEM')
     ! Find nearest vertical index for temperature (t-component)
     idx_near = nearest_index(Abk(1)%t(:,1), z)
     
     ! Calculate ensemble mean at that specific depth using a loop
     mmean = 0._dp
     k = 0
     do i = 1, nnkn
       if (Abk(1)%t(idx_near,i) > -90.) then
          mmean = mmean + Abk(1)%t(idx_near,i)
          k = k + 1
       end if
     end do
     mmean = mmean / real(k, dp)

     if ( abs(vmean - mmean) > TBIAS_MAX ) then
        write(*,*) 'Warning: bias in observations, correcting (mod/obs): ', mmean, vmean
        v = v - vmean + mmean
     end if

  case ('0DSAL')
     ! Find nearest vertical index for salinity (s-component)
     idx_near = nearest_index(Abk(1)%s(:,1), z)
     
     ! Calculate ensemble mean at that specific depth using a loop
     mmean = 0._dp
     k = 0
     do i = 1, nnkn
       if (Abk(1)%s(idx_near,i) > -90.) then
          mmean = mmean + Abk(1)%s(idx_near,i)
          k = k + 1
       end if
     end do
     mmean = mmean / real(k, dp)

     if ( abs(vmean - mmean) > SBIAS_MAX ) then
        write(*,*) 'Warning: bias in observations, correcting (mod/obs): ', mmean, vmean
        v = v - vmean + mmean
     end if

  case default
     return

  end select

contains

    function nearest_index(vec, z0_local) result(idx)
        ! Internal function to find the index of the closest value in a vector
        real(dp), intent(in) :: vec(:)
        real(dp), intent(in) :: z0_local
        integer              :: idx
        idx = minloc(abs(vec - z0_local), 1)
    end function nearest_index

end subroutine check_mean

!======================================================================
!  screen_observation
!======================================================================
subroutine screen_observation(obs, x_ens, nmem, obs_std, k_std, k_rel, accept_obs)
  implicit none

  !------------------------------------------------------------------
  ! Inputs
  !------------------------------------------------------------------
  real(dp), intent(in) :: obs              ! observed scalar value
  real(dp), intent(in) :: x_ens(nmem)      ! ensemble model values
  integer, intent(in)  :: nmem             ! number of ensemble members
  real(dp), intent(in) :: obs_std          ! observation std deviation
  real(dp), intent(in) :: k_std            ! threshold multiplier for std (e.g., 3.0)
  real(dp), intent(in) :: k_rel            ! relative threshold (0.1 - 0.5 typical)

  !------------------------------------------------------------------
  ! Output
  !------------------------------------------------------------------
  logical, intent(inout) :: accept_obs       ! .true. = keep observation

  !------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------
  real(dp) :: innovation
  real(dp) :: mean_model
  real(dp) :: scale_value
  real(dp) :: spread
  real(dp) :: thresh

  !==================================================================
  ! Compute ensemble mean
  !==================================================================
  mean_model = sum(x_ens) / real(nmem, dp)

  !==================================================================
  ! Innovation = obs - H(x)
  !==================================================================
  innovation = obs - mean_model

  !==================================================================
  ! Physical scale for relative screening
  ! (prevents rejecting deep-water Z or large T/S values)
  !==================================================================
  scale_value = max(abs(mean_model), abs(obs))

!  !==================================================================
!  ! 1) Std-based check
!  !==================================================================
!  if (abs(innovation) > k_std * obs_std) then
!     write(*,*) 'Std-based check not passed: ', abs(innovation), k_std*obs_std
!     accept_obs = .false.
!     return
!  end if

  !==================================================================
  ! 2) Relative-scale check
  !==================================================================
  if (scale_value > 0.0_dp) then
     if (abs(innovation) > k_rel * scale_value) then
	write(*,*) 'Relative-scale check not passed: ', abs(innovation), k_rel*scale_value
        accept_obs = .false.
        return
     end if
  end if

  !==================================================================
  ! 3) Spread check
  !==================================================================
  spread = sqrt(sum(x_ens**2)/real(nmem, dp) - mean_model * mean_model)
  thresh = 3._dp * sqrt(obs_std**2 + spread**2)
  if (abs(innovation) > thresh) then
	write(*,*) 'Spread check not passed: ', abs(innovation), thresh
	accept_obs = .false.
	return
  end if

  ! Passed both checks -> accept
  accept_obs = .true.

end subroutine screen_observation

!======================================================================
!  check_spread
!======================================================================
  subroutine check_spread(d, stdv, mval, mvalm)
    use levels
    use mod_ens_state
    implicit none
    real(dp), intent(in)    :: d
    real(dp), intent(in)    :: mvalm, mval(nrens)
    real(dp), intent(inout) :: stdv
    integer :: ne
    integer, save :: icall = 0
    real(dp) :: ens_std, stdv_new

    if (icall == 0) write(*,*) 'Obs error variance limited to KSTD times the ensemble spread (Sakov 2012).'
    if (KSTD <= 0.0) return

    ens_std = 0.0
    do ne = 1, nrens
      ens_std = ens_std + (mval(ne) - mvalm)**2
    end do
    ens_std = sqrt( ens_std / real(nrens - 1, dp) )

    if (abs(d) > KSTD * ens_std) then
      stdv_new = sqrt( sqrt( (ens_std**2 + stdv**2)**2 + ( (1.0/KSTD) * ens_std * d )**2 ) - ens_std**2 )
      write(*,'(a22,1x,3f8.3)') 'STDe, STDo, STDo_new: ', ens_std, stdv, stdv_new
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
    real(dp), intent(in)    :: d1, d2
    real(dp), intent(inout) :: stdv
    real(dp), intent(in)    :: um(nrens), vm(nrens), umm, vmm
    real(dp) :: cs(nrens), csm, inn
    integer :: ne
    real(dp) :: ens_std, stdv_new

    if (KSTD <= 0.0) return

    cs  = sqrt( um**2 + vm**2 )
    csm = sqrt( umm**2 + vmm**2 )

    ens_std = 0.0
    do ne = 1, nrens
      ens_std = ens_std + (cs(ne) - csm)**2
    end do
    ens_std = sqrt( ens_std / real(nrens - 1, dp) )

    inn = sqrt( d1**2 + d2**2 )

    if (abs(inn) > KSTD * ens_std) then
      stdv_new = sqrt( sqrt( (ens_std**2 + stdv**2)**2 + ( (1.0/KSTD) * ens_std * inn )**2 ) - ens_std**2 )
      if (verbose) write(*,'(a18,2f8.4)') ' changing obs std ', stdv, stdv_new
      stdv = stdv_new
    end if
  end subroutine check_spread_speed

end module mod_manage_obs
