module mod_stoch
    use iso_fortran_env, only : dp => real64, sp => real32, error_unit
    implicit none

contains

    subroutine generate_gaussian_noise(vec)
        real(dp), intent(out) :: vec(:)
        real(dp) :: u1, u2, r, theta
        real(dp), parameter :: PI = 3.141592653589793_dp
        integer  :: i, n
        
        n = size(vec)
        do i = 1, n, 2
            call random_number(u1)
            call random_number(u2)
            u1 = max(u1, tiny(u1))
            
            r     = sqrt(-2.0_dp * log(u1))
            theta = 2.0_dp * PI * u2
            
            vec(i) = r * cos(theta)
            if (i + 1 <= n) vec(i+1) = r * sin(theta)
        end do
    end subroutine generate_gaussian_noise

subroutine perturbe_0d(n_size, pvec, dt, tau)
    implicit none

    ! --- Arguments ---
    integer,  intent(in)    :: n_size
    real(dp), intent(out)   :: pvec(:)
    real(dp), intent(in)    :: dt      ! Time step between records (seconds)
    real(dp), intent(in)    :: tau     ! Decorrelation time scale (seconds)

    ! --- Local Variables ---
    real(dp), allocatable   :: raw_noise(:)
    real(dp)                :: alpha, ave
    integer                 :: n

    ! 1. Calculate Alpha based on physical time scales (seconds)
    ! alpha = exp(-dt / tau). If tau is very large, alpha -> 1 (very red)
    if (tau > 0.0_dp) then
        alpha = exp(-dt / tau)
    else
        alpha = 0.0_dp ! White noise
    end if

    allocate(raw_noise(n_size))
    call generate_gaussian_noise(raw_noise)
    
    ! 2. Red Noise AR(1) Process
    ! We use sqrt(1 - alpha**2) to keep the theoretical variance = 1
    pvec(1) = raw_noise(1)
    do n = 2, n_size
        pvec(n) = alpha * pvec(n-1) + sqrt(1.0_dp - alpha**2) * raw_noise(n)
    end do
    
    ! 3. Final Bias Correction (Mean = 0)
    ! IMPORTANT: We ONLY remove the mean. We DO NOT force std_dev = 1 
    ! because that would overwrite the AR(1) statistical structure.
    ave = sum(pvec) / real(n_size, dp)
    pvec = pvec - ave
    
    deallocate(raw_noise)
end subroutine perturbe_0d
	
subroutine read_ts(linit, filein, kdim, v, vtime, dstring)
    use iso8601          ! Assuming this provides string2date and dts_to_abs_time
    implicit none

    ! --- Arguments ---
    logical,          intent(in)    :: linit    ! .true. to initialize/count, .false. to read next
    character(len=*), intent(in)    :: filein   ! Input filename
    integer,          intent(inout) :: kdim     ! Total number of records found
    real(dp),         intent(out)   :: v        ! Current value read
    real(dp),         intent(out)   :: vtime    ! Absolute time (double precision)
    character(len=80),intent(out)   :: dstring  ! Date string from file

    ! --- Local Variables ---
    integer :: ios, ierr, date, time
    integer, save :: unit_ts = 26  ! Save the unit number across calls

    ! Default null values
    v = -999.0_sp
    vtime = -999.0_dp

    if (linit) then
        ! --- CASE: INITIALIZATION (Count records) ---
        open(unit=unit_ts, file=trim(filein), status='old', form='formatted', iostat=ios)
        if (ios /= 0) then
            write(error_unit, '(A,A)') "Error: read_ts could not open file: ", trim(filein)
            error stop
        end if

        kdim = 0
        do
            read(unit_ts, *, iostat=ios) dstring, v
            if (ios < 0) exit ! End of file reached
            if (ios > 0) error stop "read_ts: format error during counting"

            ! Check for NaNs (using intrinsic ieee_is_nan if available, or simple check)
            if (v /= v) error stop "read_ts: input file contains NaNs"

            call string2date(trim(dstring), date, time, ierr)
            if (ierr /= 0) error stop "read_ts: invalid date string format"
            
            kdim = kdim + 1
        end do

        rewind(unit_ts)

    else
        ! --- CASE: SEQUENTIAL READING ---
        read(unit_ts, *, iostat=ios) dstring, v
        
        if (ios < 0) then
            ! End of file: close and return
            close(unit_ts)
            return
        else if (ios > 0) then
            error stop "read_ts: error reading data line"
        end if

        ! Process the valid line
        call string2date(trim(dstring), date, time, ierr)
        if (ierr /= 0) error stop "read_ts: error converting string to date"
        
        call dts_to_abs_time(date, time, vtime)
    end if

end subroutine read_ts

! ======================================================================
! Subroutine to apply perturbations to a scalar time-series value
! ======================================================================
subroutine apply_ts_perturbation(v_nom, std_obs, pvec, nrens, v_pert)
    implicit none

    ! --- Arguments ---
    real(dp), intent(in)  :: v_nom        ! Nominal value from time-series (v)
    real(dp), intent(in)  :: std_obs      ! Standard deviation (amplitude) of error
    real(dp), intent(in)  :: pvec(:)      ! Stochastic vector from perturbe_0d (size nrens)
    integer,  intent(in)  :: nrens        ! Number of ensemble members
    real(dp), intent(out) :: v_pert(:)    ! Perturbed values for each member

    ! --- Local Variables ---
    integer :: ne

    ! 1. Initialize ensemble perturbed vector
    if (size(v_pert) < nrens) error stop "apply_ts_perturbation: output vector too small"

    ! 2. Apply perturbation logic
    ! Member 1 is the control (unperturbed)
    v_pert(1) = v_nom

    ! Members 2 to nrens are perturbed
    do ne = 2, nrens
        ! The formula: perturbed = nominal + (normalized_noise * amplitude)
        ! pvec is assumed to have mean=0 and std_dev ~1
        v_pert(ne) = v_nom + (pvec(ne-1) * std_obs)
    end do

end subroutine apply_ts_perturbation


end module mod_stoch

!=================================================================================================
! Main
!=================================================================================================

program perturbe_ts
    use iso_fortran_env, only : dp => real64, sp => real32, error_unit
    use mod_stoch
    implicit none

    ! --- Variables ---
    character(len=255)      :: filename, basename, out_file
    integer                 :: dot_pos
    integer                 :: nrens, kdim, i, ne
    real(dp)                :: std_obs, tau, vmin, vmax
    real(dp)                :: v_nom
    real(dp)                :: vtime
    character(len=80)       :: dstring
    character(len=3)        :: nel

    ! Arrays for full time-series storage
    real(dp), allocatable   :: ts_values(:), pvec(:), v_pert(:)
    real(dp), allocatable   :: ts_times(:)
    character(len=80), allocatable :: ts_dates(:)
    real(dp)                :: dt

    ! --- 1. CLI Arguments Parsing ---
    if (command_argument_count() /= 6) then
        write(error_unit, '(A)') "Usage: ./perturbe_ts <filename> <nrens> <STD> <Tau> <MIN> <MAX>"
        write(error_unit, '(A)') "  Tau: Correlation scale (0.0: white noise, >0.0: red noise)"
        stop 1
    end if

    call get_command_argument(1, filename)
    dot_pos = index(filename, '.', back=.true.)
    if (dot_pos > 0) then
	    basename = filename(1:dot_pos-1)
    else
	    basename = trim(filename)
    end if

    block
        character(len=62) :: arg
        call get_command_argument(2, arg); read(arg, *) nrens
        call get_command_argument(3, arg); read(arg, *) std_obs
        call get_command_argument(4, arg); read(arg, *) tau
        call get_command_argument(5, arg); read(arg, *) vmin
        call get_command_argument(6, arg); read(arg, *) vmax
    end block

    ! --- 2. Read Time-Series (Initialization & Loading) ---
    ! Phase A: Count records
    call read_ts(.true., filename, kdim, v_nom, vtime, dstring)

    allocate(ts_values(kdim), ts_times(kdim), ts_dates(kdim))

    ! Phase B: Load data into memory
    do i = 1, kdim
        call read_ts(.false., filename, kdim, ts_values(i), ts_times(i), ts_dates(i))
    end do
    dt = (ts_times(kdim) - ts_times(1)) / (kdim - 1)
    write(*, '(A,I0,A)') "Loaded ", kdim, " records from " // trim(filename)

    ! --- 3. Generate Noise and Apply Perturbations ---
    allocate(pvec(kdim))    ! We need one noise value per time step
    allocate(v_pert(nrens)) ! Buffer for ensemble at each step

    ! Loop over each ensemble member to create its own file
    do ne = 1, nrens
        ! Format filename: inputname_001.ts
        write(nel, '(I3.3)') ne - 1
        out_file = trim(basename) // "_" // nel // ".dat"

        open(unit=30, file=trim(out_file), status='replace', form='formatted')

        ! Generate a new stochastic path for this specific member (except for member 1)
        if (ne == 1) then
            pvec = 0.0_sp  ! Control run: no noise
        else
            ! Generate red noise vector for the entire time-series length
            call perturbe_0d(kdim, pvec, dt, tau)
        end if

        ! Write the perturbed time series
        do i = 1, kdim
            ! Apply noise (pvec(i)) scaled by std_obs to the nominal value ts_values(i)
            ! We use a simplified inline version of apply_ts_perturbation for efficiency
            v_nom = ts_values(i) + (pvec(i) * std_obs)

	    v_nom = max(vmin, min(vmax, v_nom))

            write(30, '(A20,2X,F12.5)') ts_dates(i), v_nom
        end do

        close(30)
        write(*, '(A,A)') "Created: ", trim(out_file)
    end do

    write(*, '(A)') "Ensemble generation complete."

end program perturbe_ts

