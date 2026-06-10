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
subroutine read_observations
    use iso_fortran_env, only: error_unit, dp => real64
    use mod_obs_states 
    implicit none

    integer            :: ios, u, n, nfile
    integer            :: kend, kinit_loc, nobs_file
    character(len=80)  :: line
    real(dp)           :: tobs, vobs, vmean
    integer            :: statobs

    ! Reset counters
    islev = 0; isvel = 0; istemp = 0; issalt = 0
    n_0dlev = 0; n_0dtemp = 0; n_0dsalt = 0; n_2dvel = 0
    nobs_tot = 0

    ! 1. Open and count entries in the main list file
    open(newunit=u, file=obsfile, status='old', form='formatted', iostat=ios)
    if (ios /= 0) error stop 'read_observations: error opening '//trim(obsfile)

    nfile = 0
    do
        read(u, '(A)', iostat=ios) line
        if (ios < 0) exit
        nfile = nfile + 1
    end do
    rewind(u)

    allocate(ofile(nfile))

    ! 2. Load metadata and pre-count types (No file opening here yet)
    do n = 1, nfile
        read(u, *) ofile(n)%ty, ofile(n)%name, ofile(n)%x, ofile(n)%y, &
                   ofile(n)%z, ofile(n)%std, ofile(n)%rhol
        
        select case (trim(ofile(n)%ty))
        case ('0DLEV'); n_0dlev = n_0dlev + 1; islev = 1
        case ('0DTEM'); n_0dtemp = n_0dtemp + 1; istemp = 1
        case ('0DSAL'); n_0dsalt = n_0dsalt + 1; issalt = 1
        case ('2DVEL'); n_2dvel = n_2dvel + 1; isvel = 1
        end select
    end do
    close(u)

    ! 3. Consistency check
    if ((islev + isvel + istemp + issalt) > 1) error stop 'Multiple obs types'

    ! 4. Allocation (Exact size from metadata pre-count)
    if (n_0dlev  > 0) allocate(o0dlev(n_0dlev))
    if (n_0dtemp > 0) allocate(o0dtemp(n_0dtemp))
    if (n_0dsalt > 0) allocate(o0dsalt(n_0dsalt))
    if (n_2dvel  > 0) allocate(o2dvel(n_2dvel))

    ! 5. SINGLE PASS DATA READING
    ! We reset counters to use them as current indices for the arrays
    n_0dlev = 0; n_0dtemp = 0; n_0dsalt = 0; n_2dvel = 0

    do n = 1, nfile
        select case (trim(ofile(n)%ty))
        case ('0DLEV', '0DTEM', '0DSAL')
            ! Call the reading routine ONLY ONCE per file
            ! kinit=0, linit=.false. means: read and return the data immediately
            call read_scalar_0d(ofile(n)%ty, .false., trim(ofile(n)%name), TEPS, &
                                0, kend, tobs, vobs, vmean, statobs)
            
            if (kend > 0) then
                nobs_tot = nobs_tot + 1
                select case (trim(ofile(n)%ty))
                case ('0DLEV')
                    n_0dlev = n_0dlev + 1
                    call store_data(o0dlev(n_0dlev), n, tobs, vobs, statobs)
                    call check_mean('0DLEV', vmean, ofile(n)%z, vobs)
                case ('0DTEM')
                    n_0dtemp = n_0dtemp + 1
                    call store_data(o0dtemp(n_0dtemp), n, tobs, vobs, statobs)
                case ('0DSAL')
                    n_0dsalt = n_0dsalt + 1
                    call store_data(o0dsalt(n_0dsalt), n, tobs, vobs, statobs)
                end select
            end if

        case ('2DVEL')
            n_2dvel = n_2dvel + 1
            call read_2dvel(trim(ofile(n)%name), n, n_2dvel, TEPS, nobs_file, ofile(n)%std)
            nobs_tot = nobs_tot + 2 * nobs_file
        end select
    end do

contains

    subroutine store_data(target, idx, t, v, s)
        type(scalar_0d), intent(out) :: target
        integer, intent(in) :: idx, s
        real(dp), intent(in) :: t, v
        target%t = t; target%x = ofile(idx)%x; target%y = ofile(idx)%y
        target%z = ofile(idx)%z; target%val = v; target%std = ofile(idx)%std
        target%stat = s; target%id = idx; target%rhol = ofile(idx)%rhol
    end subroutine
    
end subroutine read_observations

!======================================================================
    subroutine read_scalar_0d(olabel, linit, filin, eps, kinit, kend, atime_obs, vv, &
                            vmean, ostatusv)
    use iso_fortran_env, only: dp => real64
    use iso8601
    implicit none

    ! Arguments
    character(len=*), intent(in)  :: olabel    ! Observation type label
    logical,          intent(in)  :: linit     ! Keep for compatibility (not strictly needed now)
    character(len=*), intent(in)  :: filin     ! Filename to read
    real(dp),         intent(in)  :: eps       ! Time tolerance (assimilation window)
    integer,          intent(in)  :: kinit     ! Initial offset (usually 0 in single-pass)
    integer,          intent(out) :: kend      ! Number of valid observations found
    real(dp),         intent(out) :: atime_obs ! Absolute time of the stored observation
    real(dp),         intent(out) :: vv        ! Value of the stored observation
    real(dp),         intent(out) :: vmean     ! Mean value of all valid obs in file
    integer,          intent(out) :: ostatusv  ! Status of the stored observation

    ! Local variables
    integer           :: ios, u, ierr, date, time
    integer           :: k_count
    real(dp)          :: v, t_tmp
    character(len=80) :: dstring
    logical           :: stored

    ! Initialization
    ostatusv = -999
    vmean    = 0.0_dp
    vv       = 0.0_dp
    atime_obs = 0.0_dp
    k_count  = 0
    stored   = .false.

    ! Open the data file using a safe unit number
    open(newunit=u, file=trim(filin), status='old', form='formatted', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Warning: read_scalar_0d could not open file: ', trim(filin)
        kend = 0
        return
    end if

    do
      ! Read date string and value
      read(u, *, iostat=ios) dstring, v
      if (ios < 0) exit ! End of file
      if (ios > 0) then
          write(*,*) 'Error: reading data in ', trim(filin)
          exit
      end if

      ! Check for NaN (v /= v is true only if v is NaN)
      if (v /= v) v = OFLAG

      ! Convert ISO8601 string to absolute time
      call string2date(trim(dstring), date, time, ierr)
      if (ierr /= 0) then
          write(*,*) 'Warning: bad date string in ', trim(filin)
          cycle
      end if
      call dts_to_abs_time(date, time, t_tmp)

      ! Check if the observation falls within the assimilation window [atime_an - eps, atime_an + eps]
      if (abs(t_tmp - atime_an) <= eps) then
        
        k_count = k_count + 1
        vmean = vmean + v
        
        ! Store only the first valid observation encountered
        if (.not. stored) then
          vv = v
          atime_obs = t_tmp
          ostatusv  = 0 ! Initialize status
          ! Check quality/bounds
          call check_obs(olabel, vv, vv, OFLAG, ostatusv)
          stored = .true.
        end if
      end if
    end do
    close(u)

    ! Finalize outputs
    kend = k_count
    
    if (k_count > 0) then
        vmean = vmean / real(k_count, dp)
    else
        vmean = OFLAG
        ostatusv = -999
    end if

  end subroutine read_scalar_0d

!======================================================================
  subroutine read_2dvel(filin, fid, nrec, eps, nobs, ostd)
    use iso_fortran_env, only: dp => real64, sp => real32
    implicit none
    character(len=*), intent(in)  :: filin
    integer,          intent(in)  :: fid
    integer,          intent(in)  :: nrec
    real(dp),         intent(in)  :: eps, ostd
    integer,          intent(out) :: nobs

    integer :: ios, u, i, ii, jj, ix, iy
    integer :: np, iformat, iunit, irec
    integer :: nvers, lmax, nvar, ntype
    integer :: nlvddi, nx, ny
    integer :: datetime(2), ierr
    
    real(dp) :: tt, flag, dx, dy, x0, y0, uu, vv, atime_obs
    character(len=50)    :: string
    integer              :: ostatus
    logical              :: bdata
    
    ! Explicitly typed temporary arrays
    real(sp),  allocatable :: hhlv(:), hd(:)
    real(sp)               :: regpar(7)
    integer,   allocatable :: ilhkv(:)
    real(sp),  allocatable :: idata(:,:,:)

    nobs  = 0
    bdata = .false.
    np    = 0

    write(*,*) 'Opening velocity FEM file: ', trim(filin)
    call fem_file_read_open(trim(filin), np, iformat, iunit)

    irec = 0
    do
      irec = irec + 1

      ! Read record headers
      call fem_file_read_params(iformat, iunit, tt, nvers, np, lmax, nvar, ntype, datetime, ierr)
      if (ierr < 0) exit ! EOF

      call dts_convert_to_atime(datetime, tt, atime_obs)

      allocate(hhlv(lmax))
      nlvddi = lmax
      call fem_file_read_2header(iformat, iunit, ntype, lmax, hhlv, regpar, ierr)

      nx = nint(regpar(1)); ny = nint(regpar(2))
      x0 = regpar(3); y0 = regpar(4)
      dx = regpar(5); dy = regpar(6)
      flag = real(regpar(7), dp)

      ! Time window check
      if (abs(atime_obs - atime_an) > eps) then
        do i = 1, nvar
          call fem_file_skip_data(iformat, iunit, nvers, np, lmax, string, ierr)
        end do
        if (allocated(hhlv)) deallocate(hhlv)
        cycle
      end if

      ! Found matching time: allocate and read
      if (allocated(hhlv)) deallocate(hhlv)
      allocate(ilhkv(np), hd(np), idata(1, nx, ny))
      
      ! Allocation of the global structure for this record
      allocate(o2dvel(nrec)%x(nx,ny), o2dvel(nrec)%y(nx,ny), &
               o2dvel(nrec)%u(nx,ny), o2dvel(nrec)%v(nx,ny), &
               o2dvel(nrec)%std(nx,ny), o2dvel(nrec)%stat(nx,ny))

      do i = 1, nvar
        idata = real(flag, sp)
        call fem_file_read_data(iformat, iunit, nvers, np, lmax, string, ilhkv, hd, nlvddi, idata, ierr)
        if (ierr /= 0) then
          if (allocated(ilhkv)) deallocate(ilhkv, hd, idata)
          error stop 'read_2dvel: error reading data from '//trim(filin)
        end if
        
        select case (i)
        case (1); o2dvel(nrec)%u = real(idata(1,:,:), dp)
        case (2); o2dvel(nrec)%v = real(idata(1,:,:), dp)
        end select
      end do

      ! Vectorized coordinate calculation (faster than nested loops)
      do jj = 1, ny
         o2dvel(nrec)%y(:,jj) = y0 + dy * real(jj-1, dp)
      end do
      do ii = 1, nx
         o2dvel(nrec)%x(ii,:) = x0 + dx * real(ii-1, dp)
      end do

      o2dvel(nrec)%z   = 0.0_dp
      o2dvel(nrec)%nx  = nx
      o2dvel(nrec)%ny  = ny
      o2dvel(nrec)%std = ostd
      o2dvel(nrec)%id  = fid

      ! Quality Control and counting
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

      ! Cleanup and exit after first valid time match
      if (allocated(ilhkv)) deallocate(ilhkv, hd, idata)
      exit 
    end do

    close(iunit)
    if (.not. bdata) error stop 'read_2dvel: No valid data found in '//trim(filin)

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
!  superobs_horiz_el  
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
  integer, allocatable :: nobs(:)     
  real(dp) :: av1, av2
  real(4)  :: x4, y4                  

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
  real(dp) :: ens_spread
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
  ! 1) Std-based check
  !==================================================================
  if (abs(innovation) > k_std * obs_std) then
     write(*,*) 'Std-based check not passed (inn > k*stdo): ', abs(innovation), k_std*obs_std
     accept_obs = .false.
     return
  end if

  !==================================================================
  ! 2) Relative-scale check
  !==================================================================
  scale_value = max(abs(mean_model), abs(obs))
  if (scale_value > 0.0_dp) then
     if (abs(innovation) > k_rel * scale_value) then
	write(*,*) 'Relative-scale check not passed (inn > k*scale): ', abs(innovation), k_rel*scale_value
        accept_obs = .false.
        return
     end if
  end if

  !==================================================================
  ! 3) Spread check
  !==================================================================
  ens_spread = sqrt(sum(x_ens**2)/real(nmem, dp) - mean_model * mean_model)
  thresh = 3._dp * sqrt(obs_std**2 + ens_spread**2)
  if (abs(innovation) > thresh) then
	write(*,*) 'Spread check not passed (inn > thresh): ', abs(innovation), thresh
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
    
    ! Input/Output
    real(dp), intent(in)    :: d          ! Innovation (Observation - Ensemble Mean)
    real(dp), intent(in)    :: mvalm      ! Ensemble mean at the specific node
    real(dp), intent(in)    :: mval(nrens)! Values for all ensemble members
    real(dp), intent(inout) :: stdv       ! Observation error (Standard Deviation)
    
    ! Local variables
    integer :: ne
    integer, save :: icall = 0
    real(dp) :: var_e, var_o, var_o_new, d2
    real(dp), parameter :: tiny_spread = 1.0e-12_dp ! Guard against zero ens_spread

    if (icall == 0) then
        write(*,*) 'EnKF: R-adaptive inflation active (Sakov 2012), KSTD = ', KSTD
        icall = 1
    end if

    ! Exit if KSTD is not set or invalid
    if (KSTD <= 0.0_dp) return

    ! 1. Calculate Ensemble Variance (var_e)
    var_e = 0.0_dp
    do ne = 1, nrens
        var_e = var_e + (mval(ne) - mvalm)**2
    end do
    var_e = var_e / real(nrens - 1, dp)

    ! 2. Square the innovation and get current observation variance (var_o)
    d2 = d**2
    var_o = stdv**2

    ! 3. Sakov (2012) logic: If innovation exceeds expected ens_spread
    ! Condition: d^2 > KSTD^2 * var_e
    if (d2 > (KSTD**2 * var_e)) then
        
        ! Robustness: Ensure var_e is not zero to avoid numerical instability
        if (var_e < tiny_spread) var_e = tiny_spread
        
        ! Correct formula for inflated observation variance (var_o_new)
        ! This limits the influence of outliers by increasing R
        var_o_new = sqrt( (var_e + var_o)**2 + (d2 * var_e / KSTD**2) ) - var_e
        
        ! Safety cap: Prevent the error from growing indefinitely at unstable nodes
        if (var_o_new > 100.0_dp * var_o) var_o_new = 100.0_dp * var_o
        
        ! Update the standard deviation for the analysis step
        stdv = sqrt(max(var_o_new, var_o))
    end if

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
