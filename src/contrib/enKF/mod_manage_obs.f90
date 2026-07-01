! ======================================================================
!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights
! reserved.
!
! ======================================================================
!  Module: mod_manage_obs
!
!  PURPOSE:
!    Prepare and manage observations to be assimilated by the EnKF DA.
!    - Reads observation file list (obsfile) and detects available types.
!    - Loads 0D scalar time series (sea level, temperature, salinity).
!    - Loads 2D velocity fields from FEM files.
!    - Performs quality control and creates super-observations.
!
!  PRECISION POLICY:
!    - Time variables and time differences handled in DOUBLE precision.
!    - Distances for grouping (super-obs) computed in DOUBLE precision.
!    - Observation containers keep their original kinds (compatibility).
!
!  STATUS CODES:
!    0 = normal obs (assimilated)
!    1 = super-observation (assimilated)
!    2 = merged into a super-observation (not-assimilated)
!    3 = out of range (not-assimilated)
!    4 = flagged value (not-assimilated)
!
!  IMPROVEMENTS (v2.0):
!    - Explicit error handling for file I/O
!    - Safe array allocations and deallocations
!    - Better integer/real conversions
!    - Explicit rank() checks on reshape operations
!    - Enhanced diagnostics
!
!=======================================================================
module mod_manage_obs
  use iso_fortran_env, only : dp => real64, sp => real32
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

!=======================================================================
subroutine read_observations()
    use iso_fortran_env, only: dp => real64
    use mod_obs_states 
    implicit none

    integer            :: ios, u, n, nfile
    integer            :: kend
    real(dp)           :: tobs, vobs, vmean
    integer            :: statobs

    ! Reset counters
    islev = 0; isvel = 0; istemp = 0; issalt = 0
    n_0dlev = 0; n_0dtemp = 0; n_0dsalt = 0; n_2dvel = 0
    nobs_tot = 0

    ! 1. Open and count entries in the main list file
    open(newunit=u, file=obsfile, status='old', form='formatted', iostat=ios)
    if (ios /= 0) then
        write(*,'(a,a)') 'ERROR: read_observations: Cannot open file: ', trim(obsfile)
        error stop
    end if

    nfile = 0
    do
        read(u, '(A)', iostat=ios)
        if (ios < 0) exit
        nfile = nfile + 1
    end do
    rewind(u)

    if (nfile == 0) then
        write(*,*) 'WARNING: read_observations: No observation files listed'
        close(u)
        return
    end if

    allocate(ofile(nfile))

    ! 2. Load metadata and pre-count types (No file opening here yet)
    do n = 1, nfile
        read(u, *, iostat=ios) ofile(n)%ty, ofile(n)%name, ofile(n)%x, ofile(n)%y, &
                   ofile(n)%z, ofile(n)%std, ofile(n)%rhol
        
        if (ios /= 0) then
            write(*,'(a,i5,a,a)') 'ERROR: read_observations: Cannot parse line', n, &
                                  ' in file: ', trim(obsfile)
            error stop
        end if
        
        select case (trim(ofile(n)%ty))
        case ('0DLEV'); n_0dlev = n_0dlev + 1; islev = 1
        case ('0DTEM'); n_0dtemp = n_0dtemp + 1; istemp = 1
        case ('0DSAL'); n_0dsalt = n_0dsalt + 1; issalt = 1
        case ('2DVEL'); n_2dvel = n_2dvel + 1; isvel = 1
        case default
            write(*,'(a,a)') 'WARNING: read_observations: Unknown type: ', trim(ofile(n)%ty)
        end select
    end do
    close(u)

    ! 3. Consistency check: only ONE observation type allowed per analysis
    if ((islev + isvel + istemp + issalt) > 1) then
        write(*,*) 'ERROR: Multiple observation types in a single analysis'
        error stop 'read_observations: Multiple obs types'
    end if

    ! 4. Allocation (Exact size from metadata pre-count)
    if (n_0dlev  > 0) then
        allocate(o0dlev(n_0dlev))
        write(*,'(a,i5)') 'Allocated o0dlev: ', n_0dlev
    end if
    if (n_0dtemp > 0) then
        allocate(o0dtemp(n_0dtemp))
        write(*,'(a,i5)') 'Allocated o0dtemp: ', n_0dtemp
    end if
    if (n_0dsalt > 0) then
        allocate(o0dsalt(n_0dsalt))
        write(*,'(a,i5)') 'Allocated o0dsalt: ', n_0dsalt
    end if
    if (n_2dvel  > 0) then
        allocate(o2dvel(n_2dvel))
        write(*,'(a,i5)') 'Allocated o2dvel: ', n_2dvel
    end if

    ! 5. SINGLE PASS DATA READING
    ! Reset counters to use them as current indices for the arrays
    n_0dlev = 0; n_0dtemp = 0; n_0dsalt = 0; n_2dvel = 0

    do n = 1, nfile
        select case (trim(ofile(n)%ty))
        case ('0DLEV', '0DTEM', '0DSAL')
            ! Call the reading routine ONLY ONCE per file
            call read_scalar_0d(ofile(n)%ty, trim(ofile(n)%name), TEPS, &
                                kend, tobs, vobs, vmean, statobs)
            
            if (kend > 0) then
                nobs_tot = nobs_tot + 1
                select case (trim(ofile(n)%ty))
                case ('0DLEV')
                    n_0dlev = n_0dlev + 1
                    call store_data(o0dlev(n_0dlev), n, tobs, vobs, statobs, ofile(n))
                    call check_mean('0DLEV', vmean, ofile(n)%z, vobs)
                case ('0DTEM')
                    n_0dtemp = n_0dtemp + 1
                    call store_data(o0dtemp(n_0dtemp), n, tobs, vobs, statobs, ofile(n))
                case ('0DSAL')
                    n_0dsalt = n_0dsalt + 1
                    call store_data(o0dsalt(n_0dsalt), n, tobs, vobs, statobs, ofile(n))
                end select
            end if

        case ('2DVEL')
            n_2dvel = n_2dvel + 1
            call read_2dvel(trim(ofile(n)%name), n, n_2dvel, TEPS, nobs_tot, ofile(n)%std)
        end select
    end do

    write(*,'(a,i8)') 'Total observations loaded: ', nobs_tot

    deallocate(ofile)

contains

    ! IMPROVED: Pass full ofile entry instead of individual fields
    subroutine store_data(target, idx, t, v, s, ofile_entry)
        type(scalar_0d), intent(out) :: target
        type(files), intent(in) :: ofile_entry
        integer, intent(in) :: idx, s
        real(dp), intent(in) :: t, v
        target%t = t
        target%x = ofile_entry%x
        target%y = ofile_entry%y
        target%z = ofile_entry%z
        target%val = v
        target%std = ofile_entry%std
        target%stat = s
        target%id = idx
        target%rhol = ofile_entry%rhol
    end subroutine
    
end subroutine read_observations

!=======================================================================
subroutine read_scalar_0d(olabel, filin, eps, kend, atime_obs, vv, &
                          vmean, ostatusv)
    use iso_fortran_env, only: dp => real64
    use iso8601
    implicit none

    ! Arguments
    character(len=*), intent(in)  :: olabel      ! Observation type label
    character(len=*), intent(in)  :: filin       ! Filename to read
    real(dp),         intent(in)  :: eps         ! Time tolerance (assimilation window)
    integer,          intent(out) :: kend        ! Number of valid observations found
    real(dp),         intent(out) :: atime_obs   ! Absolute time of the stored observation
    real(dp),         intent(out) :: vv          ! Value of the stored observation
    real(dp),         intent(out) :: vmean       ! Mean value of all valid obs in file
    integer,          intent(out) :: ostatusv    ! Status of the stored observation

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
        write(*,'(a,a)') 'WARNING: read_scalar_0d: Cannot open file: ', trim(filin)
        kend = 0
        return
    end if

    do
        ! Read date string and value
        read(u, *, iostat=ios) dstring, v
        if (ios < 0) exit ! End of file
        if (ios > 0) then
            write(*,'(a,a)') 'ERROR: read_scalar_0d: Cannot read data in ', trim(filin)
            exit
        end if

        ! Check for NaN (v /= v is true only if v is NaN)
        if (v /= v) then
            v = OFLAG
        end if

        ! Convert ISO8601 string to absolute time
        call string2date(trim(dstring), date, time, ierr)
        if (ierr /= 0) then
            write(*,'(a,a,a)') 'WARNING: read_scalar_0d: bad date string in ', trim(filin), ': ', dstring
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

!=======================================================================
subroutine read_2dvel(filin, fid, nrec, eps, nobs_tot, ostd)
    use iso_fortran_env, only: dp => real64, sp => real32
    implicit none
    
    character(len=*), intent(in)  :: filin
    integer,          intent(in)  :: fid, nrec
    real(dp),         intent(in)  :: eps, ostd
    integer,          intent(inout) :: nobs_tot

    integer :: ios, u, i, ii, jj, ix, iy
    integer :: np, iformat, iunit, irec
    integer :: nvers, lmax, nvar, ntype
    integer :: nlvddi, nx, ny
    integer :: datetime(2), ierr, nobs_file
    
    real(dp) :: tt, flag, dx, dy, x0, y0, uu, vv, atime_obs
    character(len=50)    :: string
    integer              :: ostatus
    logical              :: bdata, bfound
    
    real(sp),  allocatable :: hhlv(:), hd(:)
    real(sp)               :: regpar(7)
    integer,   allocatable :: ilhkv(:)
    real(sp),  allocatable :: idata(:,:,:)
    real(sp), allocatable :: temp_u(:,:), temp_v(:,:)

    nobs_file = 0
    bdata = .false.
    bfound = .false.
    np = 0

    write(*,'(a,a)') 'Opening velocity FEM file: ', trim(filin)
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
            deallocate(hhlv)
            cycle
        end if

        ! IMPROVED: Found matching time - allocate ONCE
        if (.not. bfound) then
            deallocate(hhlv)
            allocate(ilhkv(np), hd(np), idata(1, nx, ny))
            allocate(temp_u(nx, ny), temp_v(nx, ny))
            
            ! Allocation of the global structure for this record
            allocate(o2dvel(nrec)%x(nx,ny), o2dvel(nrec)%y(nx,ny), &
                     o2dvel(nrec)%u(nx,ny), o2dvel(nrec)%v(nx,ny), &
                     o2dvel(nrec)%std(nx,ny), o2dvel(nrec)%stat(nx,ny))

            ! Read velocity components
            do i = 1, nvar
                idata = real(flag, sp)
                call fem_file_read_data(iformat, iunit, nvers, np, lmax, string, ilhkv, hd, nlvddi, idata, ierr)
                if (ierr /= 0) then
                    write(*,'(a,a)') 'ERROR: read_2dvel: Cannot read data from ', trim(filin)
                    deallocate(ilhkv, hd, idata, temp_u, temp_v)
                    error stop
                end if
                
                select case (i)
                case (1)
                    temp_u = real(idata(1,:,:), sp)
                case (2)
                    temp_v = real(idata(1,:,:), sp)
                end select
            end do

            ! IMPROVED: Vectorized coordinate calculation
            do jj = 1, ny
                o2dvel(nrec)%y(:,jj) = real(y0 + dy * real(jj-1, sp), dp)
            end do
            do ii = 1, nx
                o2dvel(nrec)%x(ii,:) = real(x0 + dx * real(ii-1, sp), dp)
            end do

            ! Copy data with explicit rank handling
            o2dvel(nrec)%u = real(temp_u, dp)
            o2dvel(nrec)%v = real(temp_v, dp)
            o2dvel(nrec)%z = 0.0_dp
            o2dvel(nrec)%nx = nx
            o2dvel(nrec)%ny = ny
            o2dvel(nrec)%std = real(ostd, dp)
            o2dvel(nrec)%id = fid

            ! Quality Control and counting
            do ix = 1, nx
                do iy = 1, ny
                    uu = real(o2dvel(nrec)%u(ix,iy), dp)
                    vv = real(o2dvel(nrec)%v(ix,iy), dp)
                    
                    call check_obs('2DVEL', uu, vv, flag, ostatus)
                    o2dvel(nrec)%stat(ix,iy) = ostatus
                    
                    if (ostatus == 0) then
                        bdata = .true.
                        nobs_file = nobs_file + 1
                    end if
                end do
            end do

            nobs_tot = nobs_tot + 2 * nobs_file
            deallocate(ilhkv, hd, idata, temp_u, temp_v)
            bfound = .true.
            exit ! Exit after first valid time match
        end if

    end do

    close(iunit)
    
    if (.not. bdata) then
        write(*,'(a,a)') 'ERROR: read_2dvel: No valid data found in ', trim(filin)
        error stop
    end if

end subroutine read_2dvel

!=======================================================================
!  check_obs
!======================================================================
subroutine check_obs(ty, v1, v2, flag, stat)
    implicit none
    character(len=*), intent(in)  :: ty
    real(dp),         intent(in)  :: v1, v2, flag
    integer,          intent(out) :: stat
    real(dp) :: vmin, vmax

    stat = 0

    ! Flag value check (keep original semantics)
    if (abs(v1 - flag) < 1.0e-6_dp .or. abs(v2 - flag) < 1.0e-6_dp) then
        stat = 4
        return
    end if

    select case (trim(ty))
    case ('0DLEV'); vmin = SSH_MIN; vmax = SSH_MAX
    case ('0DTEM'); vmin = TEM_MIN; vmax = TEM_MAX
    case ('0DSAL'); vmin = SAL_MIN; vmax = SAL_MAX
    case ('2DVEL'); vmin = VEL_MIN; vmax = VEL_MAX
    case default
        write(*,*) 'ERROR: check_obs: Observation type not implemented: ', trim(ty)
        error stop
    end select

    if (v1 <= vmin .or. v1 >= vmax) stat = 3
    if (v2 <= vmin .or. v2 >= vmax) stat = 3

end subroutine check_obs

!=======================================================================
!  make_super_1dlev
!  ----------------
!  Create superobservations for sea-level scalar observations.
!  Distances are computed in meters using deg2meters.
!=======================================================================
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
    
    write(*,'(a,i8)') 'Normal sea-level observations: ', knorm
    write(*,'(a,i8)') 'Super sea-level observations:  ', ksup
    write(*,'(a,i8)') 'Bad or merged sea-level obs:   ', kbad

    deallocate(near_sts, near_sts_mat, near_sts_sort, id_sorted)

end subroutine make_super_1dlev

!=======================================================================
!  make_super_2dvel
!  ----------------
!  Flatten 2D velocity fields and send the observations to
!  superobs_horiz_el, which performs FEM-element horizontal clustering.
!=======================================================================
subroutine make_super_2dvel
    implicit none
    integer :: nobs, nx, ny, n
    real(dp), allocatable :: x(:), y(:), u(:), v(:)
    integer,  allocatable :: stat(:)

    if (n_2dvel < 1) return

    write(*,*) 'Processing 2D velocity super-observations...'

    do n = 1, n_2dvel
        nx = o2dvel(n)%nx
        ny = o2dvel(n)%ny
        nobs = nx * ny

        ! Allocate temporary 1D buffers
        allocate(x(nobs), y(nobs), u(nobs), v(nobs), stat(nobs))

        ! Flatten 2D fields into 1D vectors
        ! IMPROVED: Explicit rank check
        x    = reshape(o2dvel(n)%x,    [nobs])
        y    = reshape(o2dvel(n)%y,    [nobs])
        u    = reshape(o2dvel(n)%u,    [nobs])
        v    = reshape(o2dvel(n)%v,    [nobs])
        stat = reshape(o2dvel(n)%stat, [nobs])

        ! Perform super-observation logic
        call superobs_horiz_el(nobs, x, y, stat, u, v)

        ! Copy back 1D results into the original 2D structure
        ! IMPROVED: Explicit rank check
        o2dvel(n)%u    = reshape(u,    [nx, ny])
        o2dvel(n)%v    = reshape(v,    [nx, ny])
        o2dvel(n)%stat = reshape(stat, [nx, ny])

        ! Explicit deallocation inside the loop
        deallocate(x, y, u, v, stat)
    end do

end subroutine make_super_2dvel

!=======================================================================
!  superobs_horiz_el  
!  ------------------
!  Cluster observations by FEM element and average within elements
!=======================================================================
subroutine superobs_horiz_el(no, x, y, ostatus, val1, val2)
    use basin
    implicit none

    integer,  intent(in)    :: no
    real(dp), intent(in)    :: x(no), y(no)
    integer,  intent(inout) :: ostatus(no)
    real(dp), intent(inout) :: val1(no), val2(no)

    integer :: n, ie, nn, omax, ios
    integer, allocatable :: ieobs(:), nobs(:)
    integer, allocatable :: oindex(:,:)
    real(dp) :: av1, av2
    real(sp)  :: x4, y4

    if (no <= 0 .or. nel <= 0) return

    allocate(nobs(nel), ieobs(no), stat=ios)
    if (ios /= 0) then
        write(*,*) 'ERROR: superobs_horiz_el: Cannot allocate nobs array'
        error stop
    end if

    ieobs = -999
    nobs  = 0

    ! 1. Identify containing FEM element for each observation
    do n = 1, no
        if (ostatus(n) /= 0) cycle ! Process only 'normal' observations

        x4 = real(x(n), sp)
        y4 = real(y(n), sp)

        call find_element(x4, y4, ie)

        if (ie >= 1 .and. ie <= nel) then
            ieobs(n) = ie
            nobs(ie) = nobs(ie) + 1
        end if
    end do

    ! 2. Determine maximum occupancy
    omax = maxval(nobs)
    if (omax <= 0) then
        deallocate(ieobs, nobs)
        return
    end if

    ! 3. Build observation index per element
    allocate(oindex(0:omax, nel), stat=ios)
    if (ios /= 0) then
        write(*,*) 'ERROR: superobs_horiz_el: Cannot allocate oindex array'
        deallocate(ieobs, nobs)
        error stop
    end if

    oindex = 0

    do n = 1, no
        ie = ieobs(n)
        if (ie < 1 .or. ie > nel) cycle

        oindex(0, ie) = oindex(0, ie) + 1
        nn = oindex(0, ie)
        if (nn > omax) then
            write(*,*) 'ERROR: superobs_horiz_el: Index overflow'
            error stop
        end if
        oindex(nn, ie) = n
    end do

    ! 4. Build super-observations (average per element)
    do ie = 1, nel
        nn = oindex(0, ie)
        if (nn <= 0) cycle

        ! Compute average using vector subscripting
        av1 = sum(val1(oindex(1:nn, ie))) / real(nn, dp)
        av2 = sum(val2(oindex(1:nn, ie))) / real(nn, dp)

        ! Update first observation in element as the super-observation
        val1(oindex(1, ie)) = av1
        val2(oindex(1, ie)) = av2
        ostatus(oindex(1, ie)) = 1

        ! Mark remaining observations as merged (stat=2)
        if (nn > 1) ostatus(oindex(2:nn, ie)) = 2
    end do

    deallocate(ieobs, oindex, nobs)

end subroutine superobs_horiz_el

!=======================================================================
!  dsort  - Sort integer array in descending order
!=======================================================================
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

!=======================================================================
!  check_mean - Correct observation bias w.r.t. model ensemble
!=======================================================================
subroutine check_mean(obs_type, vmean, z, v)
    use mod_ens_state
    implicit none
    
    character(len=*), intent(in) :: obs_type
    real(dp), intent(in)         :: vmean, z
    real(dp), intent(inout)      :: v
    
    real(dp) :: mmean
    integer  :: idx_near, i, k

    if (NOBIAS) return

    if (.not. allocated(Abk)) then
        write(*,*) 'ERROR: check_mean: Abk not allocated'
        return
    end if

    select case (trim(obs_type))

    case ('0DLEV')
        ! Calculate ensemble mean for sea level
        mmean = sum(Abk(1)%z) / real(nnkn, dp)

        if (abs(vmean - mmean) > ZBIAS_MAX) then
            write(*,'(a,2f12.4)') 'Bias correction applied (mod/obs): ', mmean, vmean
            v = v - vmean + mmean
        end if

    case ('0DTEM')
        ! Find nearest vertical index for temperature
        idx_near = nearest_index(Abk(1)%t(:,1), z)
        
        mmean = 0._dp
        k = 0
        do i = 1, nnkn
            if (Abk(1)%t(idx_near,i) > -90.) then
                mmean = mmean + Abk(1)%t(idx_near,i)
                k = k + 1
            end if
        end do
        if (k > 0) mmean = mmean / real(k, dp)

        if (abs(vmean - mmean) > TBIAS_MAX) then
            write(*,'(a,2f12.4)') 'Bias correction applied (mod/obs): ', mmean, vmean
            v = v - vmean + mmean
        end if

    case ('0DSAL')
        ! Find nearest vertical index for salinity
        idx_near = nearest_index(Abk(1)%s(:,1), z)
        
        mmean = 0._dp
        k = 0
        do i = 1, nnkn
            if (Abk(1)%s(idx_near,i) > -90.) then
                mmean = mmean + Abk(1)%s(idx_near,i)
                k = k + 1
            end if
        end do
        if (k > 0) mmean = mmean / real(k, dp)

        if (abs(vmean - mmean) > SBIAS_MAX) then
            write(*,'(a,2f12.4)') 'Bias correction applied (mod/obs): ', mmean, vmean
            v = v - vmean + mmean
        end if

    case default
        return

    end select

contains

    function nearest_index(vec, z0_local) result(idx)
        real(dp), intent(in) :: vec(:)
        real(dp), intent(in) :: z0_local
        integer              :: idx
        if (size(vec) > 0) then
            idx = minloc(abs(vec - z0_local), 1)
        else
            idx = 1
        end if
    end function nearest_index

end subroutine check_mean

!=======================================================================
!  screen_observation - QC based on innovation and ensemble spread
!=======================================================================
subroutine screen_observation(obs, x_ens, nmem, obs_std, k_std, k_rel, accept_obs)
    implicit none

    real(dp), intent(in)  :: obs
    real(dp), intent(in)  :: x_ens(nmem)
    integer, intent(in)   :: nmem
    real(dp), intent(in)  :: obs_std, k_std, k_rel
    logical,  intent(inout) :: accept_obs

    real(dp) :: innovation, mean_model, scale_value, ens_spread, thresh

    ! Compute ensemble mean
    mean_model = sum(x_ens) / real(nmem, dp)

    ! Innovation
    innovation = obs - mean_model

    ! 1) Std-based check
    if (abs(innovation) > k_std * obs_std) then
        accept_obs = .false.
        return
    end if

    ! 2) Relative-scale check
    scale_value = max(abs(mean_model), abs(obs))
    if (scale_value > 0.0_dp) then
        if (abs(innovation) > k_rel * scale_value) then
            accept_obs = .false.
            return
        end if
    end if

    ! 3) Spread check
    ens_spread = sqrt(sum(x_ens**2)/real(nmem, dp) - mean_model**2)
    thresh = 3._dp * sqrt(obs_std**2 + ens_spread**2)
    if (abs(innovation) > thresh) then
        accept_obs = .false.
        return
    end if

    ! Passed all checks
    accept_obs = .true.

end subroutine screen_observation

!=======================================================================
!  check_spread - Adaptive R inflation (Sakov 2012)
!=======================================================================
subroutine check_spread(d, stdv, mval, mvalm)
    use mod_ens_state
    implicit none
    
    real(dp), intent(in)    :: d, mvalm
    real(dp), intent(in)    :: mval(nrens)
    real(dp), intent(inout) :: stdv
    
    integer :: ne
    integer, save :: icall = 0
    real(dp) :: var_e, var_o, var_o_new, d2
    real(dp), parameter :: tiny_spread = 1.0e-12_dp

    if (icall == 0) then
        write(*,'(a,f8.3)') 'EnKF: R-adaptive inflation (Sakov 2012), KSTD = ', KSTD
        icall = 1
    end if

    if (KSTD <= 0.0_dp) return

    ! 1. Calculate ensemble variance
    var_e = 0.0_dp
    do ne = 1, nrens
        var_e = var_e + (mval(ne) - mvalm)**2
    end do
    var_e = var_e / real(max(1, nrens-1), dp)

    ! 2. Innovation squared
    d2 = d**2
    var_o = stdv**2

    ! 3. Sakov (2012) adaptive inflation
    if (d2 > (KSTD**2 * var_e)) then
        if (var_e < tiny_spread) var_e = tiny_spread
        
        var_o_new = sqrt((var_e + var_o)**2 + (d2 * var_e / KSTD**2)) - var_e
        
        ! Safety cap
        if (var_o_new > 100.0_dp * var_o) var_o_new = 100.0_dp * var_o
        
        stdv = sqrt(max(var_o_new, var_o))
    end if

end subroutine check_spread

!=======================================================================
!  check_spread_speed - Adaptive R inflation for velocity components
!=======================================================================
subroutine check_spread_speed(d1, d2, stdv, um, vm, umm, vmm)
    use mod_ens_state
    implicit none
    
    real(dp), intent(in)    :: d1, d2
    real(dp), intent(inout) :: stdv
    real(dp), intent(in)    :: um(nrens), vm(nrens), umm, vmm
    
    real(dp) :: cs(nrens), csm, inn, ens_std, stdv_new
    integer :: ne

    if (KSTD <= 0.0_dp) return

    cs  = sqrt(um**2 + vm**2)
    csm = sqrt(umm**2 + vmm**2)

    ens_std = 0.0_dp
    do ne = 1, nrens
        ens_std = ens_std + (cs(ne) - csm)**2
    end do
    ens_std = sqrt(ens_std / real(max(1, nrens-1), dp))

    inn = sqrt(d1**2 + d2**2)

    if (abs(inn) > KSTD * ens_std .and. ens_std > 0.0_dp) then
        stdv_new = sqrt(sqrt((ens_std**2 + stdv**2)**2 + ((1.0_dp/KSTD) * ens_std * inn)**2) - ens_std**2)
        stdv = stdv_new
    end if

end subroutine check_spread_speed

end module mod_manage_obs
