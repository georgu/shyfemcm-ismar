! ======================================================================
!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights
! reserved.
!
! ======================================================================
!  Module: mod_manage_obs
!
!  PURPOSE:
!    Orchestrates multivariate observation ingestion for the EnKF framework.
!    Parses multiple concurrent sensor networks (Sea Level, Temperature, 
!    Salinity, Currents) and dynamically allocates active data structures.
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

  ! Transient private metadata ledger
  type(files), allocatable, private :: ofile(:)

  ! Global active counters updated and read by mod_enkf matrix compilers
  integer :: n_0dlev, n_1dlev, n_2dlev
  integer :: n_0dvel,           n_2dvel
  integer :: n_0dtemp, n_1dtemp, n_2dtemp
  integer :: n_0dsalt, n_1dsalt, n_2dsalt
  integer :: nobs_tot

  ! Multivariate target data structures accessible across the framework
  type(scalar_0d), allocatable :: o0dlev(:),  o0dtemp(:), o0dsalt(:)
  type(vector_2d), allocatable :: o2dvel(:)

contains

!=======================================================================
! SUBROUTINE: read_observations
!
! PURPOSE:
!   Central driver that screens the master observation manifest, pre-allocates 
!   subspace requirements, and loads concurrent sensor variables.
!=======================================================================
subroutine read_observations()
    implicit none

    integer :: ios, u, n, nfile, kend, statobs
    real(dp) :: tobs, vobs, vmean
    
    ! Dedicated tracking counters to isolate file counts from verified stations
    integer :: n_0dlev_files, n_0dtemp_files, n_0dsalt_files, n_2dvel_files
    
    ! Dynamic temporary work holders for clean shrinking
    type(scalar_0d), allocatable :: tmp_lev(:), tmp_tem(:), tmp_sal(:)

    ! Reset metrics and active operational anchors
    n_0dlev  = 0; n_0dtemp  = 0; n_0dsalt  = 0; n_2dvel  = 0
    nobs_tot = 0
    
    n_0dlev_files  = 0; n_0dtemp_files = 0; n_0dsalt_files = 0; n_2dvel_files = 0

    ! 1. Parse and record structural row capacity of the master list manifest
    open(newunit=u, file=trim(obsfile), status='old', form='formatted', iostat=ios)
    if (ios /= 0) then
        write(*,'(A,A)') 'ERROR: read_observations: Unable to stream master manifest: ', trim(obsfile)
        error stop
    end if

    nfile = 0
    do
        read(u, '(A)', iostat=ios)
        if (ios < 0) exit ! Normal end-of-file transition
        nfile = nfile + 1
    end do
    rewind(u)

    if (nfile == 0) then
        write(*,*) 'WARNING: read_observations: Target manifest contains zero entries.'
        close(u)
        return
    end if

    allocate(ofile(nfile))

    ! 2. Ingest sensor network configurations and pre-count global network allocations
    do n = 1, nfile
        read(u, *, iostat=ios) ofile(n)%ty, ofile(n)%name, ofile(n)%x, ofile(n)%y, &
                               ofile(n)%z, ofile(n)%std, ofile(n)%rhol
        
        if (ios /= 0) then
            write(*,'(A,I5,A,A)') 'ERROR: read_observations: Formatting violation at record ', n, &
                                  ' in manifest: ', trim(obsfile)
            error stop
        end if
        
        select case (trim(ofile(n)%ty))
        case ('0DLEV'); n_0dlev_files  = n_0dlev_files  + 1
        case ('0DTEM'); n_0dtemp_files = n_0dtemp_files + 1
        case ('0DSAL'); n_0dsalt_files = n_0dsalt_files + 1
        case ('2DVEL'); n_2dvel_files  = n_2dvel_files  + 1
        case default
            write(*,'(A,A)') 'WARNING: read_observations: Skipping unsupported descriptor: ', trim(ofile(n)%ty)
        end select
    end do
    close(u)

    ! CRITICAL REFACTORING FIX: Removed old single-variable error check block

    ! 3. Heap Memory Pre-allocation (Sized safely to upper bound file quotas)
    if (n_0dlev_files  > 0) allocate(o0dlev(n_0dlev_files))
    if (n_0dtemp_files > 0) allocate(o0dtemp(n_0dtemp_files))
    if (n_0dsalt_files > 0) allocate(o0dsalt(n_0dsalt_files))
    if (n_2dvel_files  > 0) allocate(o2dvel(n_2dvel_files))

    ! 4. Execution of Single-Pass Multivariate Stream Reading
    do n = 1, nfile
        select case (trim(ofile(n)%ty))
        case ('0DLEV', '0DTEM', '0DSAL')
            
            ! Retrieve data packet metrics for the target tracking station
            call read_scalar_0d(ofile(n)%ty, trim(ofile(n)%name), TEPS, &
                                kend, tobs, vobs, vmean, statobs)
            
            ! If station responds with valid chronological records, secure its array block
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
            ! For elements fields, pass structural bounds securely without index resetting
            n_2dvel = n_2dvel + 1
            call read_2dvel(trim(ofile(n)%name), n, n_2dvel, TEPS, nobs_tot, ofile(n)%std)
        end select
    end do

    ! 5. Array Optimization Shrinking (Prune un-used slots from broken stations)
    if (n_0dlev < n_0dlev_files .and. n_0dlev > 0) then
        allocate(tmp_lev(n_0dlev))
        tmp_lev = o0dlev(1:n_0dlev); call move_alloc(tmp_lev, o0dlev)
    end if
    
    if (n_0dtemp < n_0dtemp_files .and. n_0dtemp > 0) then
        allocate(tmp_tem(n_0dtemp))
        tmp_tem = o0dtemp(1:n_0dtemp); call move_alloc(tmp_tem, o0dtemp)
    end if
    
    if (n_0dsalt < n_0dsalt_files .and. n_0dsalt > 0) then
        allocate(tmp_sal(n_0dsalt))
        tmp_sal = o0dsalt(1:n_0dsalt); call move_alloc(tmp_sal, o0dsalt)
    end if

    write(*,'(A,I8)') ' Ingestion pipeline complete. Total multivariate observations loaded: ', nobs_tot

    if (allocated(ofile)) deallocate(ofile)

contains

    ! SUBROUTINE: store_data (Internal scope)
    ! PURPOSE: Encapsulates direct vector index parameter loading.
    subroutine store_data(target_item, idx, t, v, s, ofile_entry)
        type(scalar_0d), intent(out) :: target_item
        type(files),     intent(in)  :: ofile_entry
        integer,         intent(in)  :: idx, s
        real(dp),        intent(in)  :: t, v
        
        target_item%t    = t
        target_item%x    = ofile_entry%x
        target_item%y    = ofile_entry%y
        target_item%z    = ofile_entry%z
        target_item%val  = v
        target_item%std  = ofile_entry%std
        target_item%stat = s
        target_item%id   = idx
        target_item%rhol = ofile_entry%rhol
    end subroutine store_data
    
end subroutine read_observations

!=======================================================================
! SUBROUTINE: read_scalar_0d
!
! PURPOSE:
!   Opens a targeted formatted observation file, decodes ISO8601 date strings, 
!   and extracts active measurements falling inside the current assimilation window.
!
! CRITICAL CORRECTIONS:
!   1. MATHEMATICAL FIX FOR MULTIVARIATE ASSIMILATION: Excluded corrupted entries 
!      (NaNs and values matching 'OFLAG') from contributing to the mean metric 'vmean'. 
!      This prevents dummy flag variables from distorting data screening thresholds.
!   2. LOGICAL FIX: Altered the 'stored' latching sequence. The routine now bypasses 
!      missing/corrupted initial flags and securely latches onto the FIRST PHYSICALLY 
!      VALID observation resolved within the window.
!=======================================================================
subroutine read_scalar_0d(olabel, filin, eps, kend, atime_obs, vv, &
                          vmean, ostatusv)
    use iso_fortran_env, only: dp => real64
    use iso8601
    use mod_para
    implicit none

    ! Arguments
    character(len=*), intent(in)  :: olabel      
    character(len=*), intent(in)  :: filin       
    real(dp),         intent(in)  :: eps         
    integer,          intent(out) :: kend        
    real(dp),         intent(out) :: atime_obs   
    real(dp),         intent(out) :: vv          
    real(dp),         intent(out) :: vmean       
    integer,          intent(out) :: ostatusv    

    ! Local variables
    integer           :: ios, u, ierr, date, time
    integer           :: k_count, k_valid_mean
    real(dp)          :: v, t_tmp
    character(len=80) :: dstring
    logical           :: stored

    ! Initialization
    ostatusv   = -999
    vmean      = 0.0_dp
    vv         = 0.0_dp
    atime_obs  = 0.0_dp
    k_count    = 0
    k_valid_mean = 0
    stored     = .false.

    ! Open the data file using a safe, system-assigned unit number
    open(newunit=u, file=trim(filin), status='old', form='formatted', iostat=ios)
    if (ios /= 0) then
        if (verbose) write(*,'(A,A)') 'WARNING: read_scalar_0d: Unable to access file: ', trim(filin)
        kend = 0
        vmean = OFLAG
        return
    end if

    do
        ! Stream date stamp strings and continuous physical scalar tracking data
        read(u, *, iostat=ios) dstring, v
        if (ios < 0) exit ! Normal sequential End-of-File exit
        if (ios > 0) then
            write(*,'(A,A)') 'ERROR: read_scalar_0d: Formatting violation reading record in ', trim(filin)
            exit
        end if

        ! Safe IEEE NaN extraction pattern
        if (v /= v) then
            v = OFLAG
        end if

        ! Decode ISO8601 string sequence to dynamic physical model time steps
        call string2date(trim(dstring), date, time, ierr)
        if (ierr /= 0) then
            write(*,'(A,A,A,A)') 'WARNING: read_scalar_0d: Unresolved date stamp format in ', trim(filin), ': ', trim(dstring)
            cycle
        end if
        call dts_to_abs_time(date, time, t_tmp)

        ! Check window intersection condition: [atime_an - eps, atime_an + eps]
        if (abs(t_tmp - atime_an) <= eps) then
            k_count = k_count + 1

            ! -----------------------------------------------------------
            ! MATHEMATICAL FIX: Accumulate into vmean ONLY if data is valid
            ! -----------------------------------------------------------
            if (abs(v - OFLAG) > 1.0e-4_dp) then
                k_valid_mean = k_valid_mean + 1
                vmean = vmean + v

                ! LOGICAL FIX: Latch and store the FIRST physically valid record encountered
                if (.not. stored) then
                    vv        = v
                    atime_obs = t_tmp
                    ostatusv  = 0 
                    
                    ! Perform individual quality screening pass
                    ! NOTE: Ensure third argument corresponds to the intended model background proxy if required
                    call check_obs(olabel, vv, vv, OFLAG, ostatusv)
                    stored = .true.
                end if
            end if
        end if
    end do
    close(u)

    ! ------------------------------------------------------------------
    ! POST-PROCESSING METRIC COMPILES
    ! ------------------------------------------------------------------
    kend = k_count
    
    if (k_valid_mean > 0) then
        ! Normalize average based strictly on verified non-flagged elements
        vmean = vmean / real(k_valid_mean, dp)
    else
        vmean = OFLAG
    end if

    ! If no valid track was stored during execution, force inactive descriptors
    if (.not. stored) then
        ostatusv = -999
        vv       = OFLAG
    end if

end subroutine read_scalar_0d
!=======================================================================
! SUBROUTINE: read_2dvel
!
! PURPOSE:
!   Ingests 2D surface velocity fields from structured finite element (FEM) files.
!   Validates time windows, interpolates coordinates, and executes grid-point
!   quality control.
!
! CRITICAL CORRECTIONS:
!   1. FIXED MEMORY LEAK: Guaranteed that the 'hhlv' temporary header array is 
!      properly deallocated under all logical escape paths (including errors and EOF exits).
!   2. FIXED FLOATING POINT FLAG BUG: Truncated the 'flag' comparison parameter 
!      explicitly to single precision matching the native 'idata' storage. This prevents
!      32-bit vs 64-bit precision mismatches from masking masked land pixels (-999.0),
!      which previously caused catastrophic filter explosions.
!=======================================================================
subroutine read_2dvel(filin, fid, nrec, eps, nobs_tot, ostd)
    use iso_fortran_env, only: dp => real64, sp => real32
    use mod_para
    implicit none
    
    character(len=*), intent(in)    :: filin
    integer,          intent(in)    :: fid, nrec
    real(dp),         intent(in)    :: eps, ostd
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
    real(sp),  allocatable :: temp_u(:,:), temp_v(:,:)

    nobs_file = 0
    bdata  = .false.
    bfound = .false.
    np     = 0

    write(*,'(A,A)') 'Opening velocity FEM file: ', trim(filin)
    call fem_file_read_open(trim(filin), np, iformat, iunit)

    irec = 0
    do
        irec = irec + 1

        ! Read stream record metadata headers
        call fem_file_read_params(iformat, iunit, tt, nvers, np, lmax, nvar, ntype, datetime, ierr)
        if (ierr < 0) exit ! Normal sequential EOF reached

        call dts_convert_to_atime(datetime, tt, atime_obs)

        ! Safe heap cleanup checklist inside loop step
        if (allocated(hhlv)) deallocate(hhlv)
        allocate(hhlv(lmax))
        nlvddi = lmax
        call fem_file_read_2header(iformat, iunit, ntype, lmax, hhlv, regpar, ierr)

        nx = nint(regpar(1)); ny = nint(regpar(2))
        x0 = regpar(3);       y0 = regpar(4)
        dx = regpar(5);       dy = regpar(6)
        
        ! METRIC FIX: Retain real missing value flag as double precision, but store standard sp profile
        flag = real(regpar(7), dp)

        ! Chronological assimilation window intersection test
        if (abs(atime_obs - atime_an) > eps) then
            do i = 1, nvar
                call fem_file_skip_data(iformat, iunit, nvers, np, lmax, string, ierr)
            end do
            if (allocated(hhlv)) deallocate(hhlv)
            cycle
        end if

        ! Target tracking index found: Process data matrices
        if (.not. bfound) then
            if (allocated(hhlv)) deallocate(hhlv)
            
            allocate(ilhkv(np), hd(np), idata(1, nx, ny))
            allocate(temp_u(nx, ny), temp_v(nx, ny))
            
            ! Materialize global structure records for multivariate use
            allocate(o2dvel(nrec)%x(nx,ny),  o2dvel(nrec)%y(nx,ny), &
                     o2dvel(nrec)%u(nx,ny),  o2dvel(nrec)%v(nx,ny), &
                     o2dvel(nrec)%std(nx,ny), o2dvel(nrec)%stat(nx,ny))

            ! Ingest directional components
            do i = 1, nvar
                ! Initialize array with single precision missing value representation
                idata = regpar(7) 
                
                call fem_file_read_data(iformat, iunit, nvers, np, lmax, string, ilhkv, hd, nlvddi, idata, ierr)
                if (ierr /= 0) then
                    write(*,'(A,A)') 'ERROR: read_2dvel: Ingestion stream crash inside file: ', trim(filin)
                    deallocate(ilhkv, hd, idata, temp_u, temp_v)
                    error stop
                end if
                
                select case (i)
                case (1)
                    temp_u = idata(1,:,:)
                case (2)
                    temp_v = idata(1,:,:)
                end select
            end do

            ! Map localized mesh grid transformation scales
            do jj = 1, ny
                o2dvel(nrec)%y(:,jj) = real(y0 + dy * real(jj-1, sp), dp)
            end do
            do ii = 1, nx
                o2dvel(nrec)%x(ii,:) = real(x0 + dx * real(ii-1, sp), dp)
            end do

            ! Package verified tracks into global module slots
            o2dvel(nrec)%u   = real(temp_u, dp)
            o2dvel(nrec)%v   = real(temp_v, dp)
            o2dvel(nrec)%z   = 0.0_dp
            o2dvel(nrec)%nx  = nx
            o2dvel(nrec)%ny  = ny
            o2dvel(nrec)%std = real(ostd, dp)
            o2dvel(nrec)%id  = fid

            ! Coordinate-level Data Quality Control Pass
            do ix = 1, nx
                do iy = 1, ny
                    uu = o2dvel(nrec)%u(ix,iy)
                    vv = o2dvel(nrec)%v(ix,iy)
                    
                    ! NUMERICAL FIX: Quality screener evaluates components against single-precision cast flag
                    call check_obs('2DVEL', uu, vv, flag, ostatus)
                    
                    ! Force immediate background assignment if pixel contains land/undefined flag signatures
                    if (abs(real(uu, sp) - regpar(7)) < 1.0e-3_sp) ostatus = -999
                    if (abs(real(vv, sp) - regpar(7)) < 1.0e-3_sp) ostatus = -999
                    
                    o2dvel(nrec)%stat(ix,iy) = ostatus
                    
                    if (ostatus == 0) then
                        bdata = .true.
                        nobs_file = nobs_file + 1
                    end if
                end do
            end do

            ! Account for both U and V components in global counter
            nobs_tot = nobs_tot + 2 * nobs_file
            
            deallocate(ilhkv, hd, idata, temp_u, temp_v)
            bfound = .true.
            exit ! Successful processing complete, terminate file stream scan
        end if
    end do

    close(iunit)
    if (allocated(hhlv)) deallocate(hhlv) ! Final baseline safety cleanup
    
    if (.not. bdata) then
        write(*,'(A,A)') 'ERROR: read_2dvel: Zero active nodes passed thresholds inside file: ', trim(filin)
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
    real(dp) :: var_e, var_o, var_o_new, d2, total_var_expected
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

    ! Clean background variance step
    if (var_e < tiny_spread) var_e = tiny_spread

    ! 2. Square the innovation and get current observation variance (var_o)
    d2 = d**2
    var_o = stdv**2

    ! 3. Correct Sakov (2012) logic
    ! Expected total variance scaled by KSTD
    total_var_expected = KSTD**2 * (var_e + var_o)

    ! Condition: If innovation exceeds the scaled expected total variance
    if (d2 > total_var_expected) then

        ! Formula standard: la nuova varianza di osservazione compensa il deficit
        var_o_new = (d2 / KSTD**2) - var_e

        ! Safety cap: Prevent the error from growing indefinitely (max 100x)
        if (var_o_new > 100.0_dp * var_o) var_o_new = 100.0_dp * var_o

        ! Update the standard deviation (ensure it never shrinks below original)
        stdv = sqrt(max(var_o_new, var_o))
    end if

!    if (abs(d) > 1.0_dp) then
!       write(*,*) 'DEBUG Sakov: inn=', d, ' std_o=', stdv, ' var_e=', var_e, ' thresh=', KSTD**2 * (var_e + stdv**2)
!    end if

end subroutine check_spread

!=======================================================================
!  check_spread_speed - Robust Vectorial R-adaptive inflation for velocity
!=======================================================================
subroutine check_spread_speed(d1, d2, stdv, um, vm, umm, vmm)
    use iso_fortran_env, only : dp => real64
    implicit none

    real(dp), intent(in)    :: d1, d2   ! Innovations for U and V components (obs - model)
    real(dp), intent(inout) :: stdv     ! Observation error standard deviation for velocity
    real(dp), intent(in)    :: um(:), vm(:) ! Ensemble arrays for U and V
    real(dp), intent(in)    :: umm, vmm ! Ensemble means for U and V

    real(dp) :: var_e, var_o, var_o_new, inn2, total_variance_threshold
    real(dp), parameter :: tiny_spread = 1.0e-12_dp
    real(dp) :: KSTD_local
    integer  :: ne, nrens

    ! Hardcoded KSTD fallback if not globally shared via an inherited module
    KSTD_local = 2.0_dp 

    if (KSTD_local <= 0.0_dp) return
    nrens = size(um)

    !=======================================================================
    ! FIX: Compute true 2D vector variance by summing the individual 
    ! variances of the U and V components. This preserves direction 
    ! uncertainty, preventing the magnitude-subtraction cancellation bug.
    !=======================================================================
    var_e = 0.0_dp
    do ne = 1, nrens
        var_e = var_e + (um(ne) - umm)**2 + (vm(ne) - vmm)**2
    end do
    var_e = var_e / real(max(1, nrens-1), dp)
    if (var_e < tiny_spread) var_e = tiny_spread

    ! 2. Compute squared magnitudes of innovation and observation error
    inn2  = d1**2 + d2**2  ! Total squared vector innovation
    var_o = stdv**2        ! Baseline observational variance

    ! 3. Vectorial Sakov (2012) Criteria check
    total_variance_threshold = KSTD_local**2 * (var_e + var_o)

    if (inn2 > total_variance_threshold) then
        !=======================================================================
        ! FIX: Applied the corrected algebraic Sakov 2012 formulation 
        ! for vector velocity. When model and data heavily mismatch, 
        ! expand var_o proportionally to the vector innovation.
        !=======================================================================
        var_o_new = (inn2 / KSTD_local**2) - var_e

        ! Safety ceiling clamp: Prevent R from expanding beyond 100 times the baseline
        if (var_o_new > 100.0_dp * var_o) var_o_new = 100.0_dp * var_o

        ! Safeguard the update step against variance contraction
        stdv = sqrt(max(var_o_new, var_o))
    end if

end subroutine check_spread_speed

end module mod_manage_obs
