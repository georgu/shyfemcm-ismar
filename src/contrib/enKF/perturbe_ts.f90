module mod_stoch
    use iso_fortran_env, only : dp => real64, error_unit
    implicit none

contains

! ======================================================================
! Subroutine to generate Gaussian noise using Box-Muller transform
! ======================================================================
subroutine generate_gaussian_noise(vec)
    use iso_fortran_env, only : dp => real64
    implicit none

    ! --- Arguments ---
    real(dp), intent(out) :: vec(:)

    ! --- Local Variables ---
    real(dp) :: u1, u2, r, theta
    real(dp), parameter :: PI = 3.141592653589793_dp
    real(dp), parameter :: LIMIT = 4.0_dp  ! Clipping threshold (4 standard deviations)
    integer  :: i, n

    n = size(vec)
    do i = 1, n, 2
        call random_number(u1)
        call random_number(u2)

        ! Ensure u1 is not zero to avoid log(0)
        u1 = max(u1, tiny(u1))

        ! Box-Muller transformation: generates two independent Gaussian values
        r     = sqrt(-2.0_dp * log(u1))
        theta = 2.0_dp * PI * u2

        vec(i) = r * cos(theta)
        if (i + 1 <= n) vec(i+1) = r * sin(theta)
    end do

    ! --- Extreme values control (Clipping) ---
    ! Prevents extreme outliers that could destabilize the model
    vec = max(-LIMIT, min(LIMIT, vec))

end subroutine generate_gaussian_noise

! ======================================================================
! Subroutine to create a Red Noise (AR1) stochastic vector
! ======================================================================
subroutine perturbe_0d(n_size, pvec, dt, tau)
    use iso_fortran_env, only : dp => real64
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

    ! 1. Calculate Alpha based on physical time scales
    ! If tau is very large, alpha -> 1 (strong temporal correlation)
    if (tau > 0.0_dp) then
        alpha = exp(-dt / tau)
    else
        alpha = 0.0_dp ! Becomes white noise if tau is zero
    end if

    allocate(raw_noise(n_size))
    call generate_gaussian_noise(raw_noise)

    ! 2. Red Noise AR(1) Process
    ! Scaling with sqrt(1 - alpha**2) maintains theoretical variance = 1
    pvec(1) = raw_noise(1)
    do n = 2, n_size
        pvec(n) = alpha * pvec(n-1) + sqrt(1.0_dp - alpha**2) * raw_noise(n)
    end do

    ! 3. Final Bias Correction (Ensure mean = 0)
    ave = sum(pvec) / real(n_size, dp)
    pvec = pvec - ave

    deallocate(raw_noise)
end subroutine perturbe_0d

! ======================================================================
! Subroutine to apply perturbations to a scalar value for EnKF
! ======================================================================
subroutine apply_ts_perturbation(v_nom, std_obs, pvec, nrens, v_pert)
    use iso_fortran_env, only : dp => real64
    implicit none

    ! --- Arguments ---
    real(dp), intent(in)  :: v_nom        ! Nominal (background) value
    real(dp), intent(in)  :: std_obs      ! Standard deviation of the error (amplitude)
    real(dp), intent(in)  :: pvec(:)      ! Stochastic noise vector (size should be >= nrens-1)
    integer,  intent(in)  :: nrens        ! Number of ensemble members
    real(dp), intent(out) :: v_pert(:)    ! Output perturbed ensemble members

    ! --- Local Variables ---
    integer :: ne

    ! 1. Safety checks
    if (size(v_pert) < nrens) error stop "apply_ts_perturbation: output vector too small"
    if (size(pvec) < nrens - 1) error stop "apply_ts_perturbation: pvec size insufficient"

    ! 2. Member 1 is the control (unperturbed)
    v_pert(1) = v_nom

    ! 3. Apply Gaussian perturbation to members 2 to nrens
    do ne = 2, nrens
        ! Formula: Background + (Normalized_Noise * Amplitude)
        v_pert(ne) = v_nom + (pvec(ne-1) * std_obs)
    end do

end subroutine apply_ts_perturbation

end module mod_stoch

	
subroutine read_ts(linit, filein, kdim, v, vtime, dstring)
    use iso_fortran_env, only : dp => real64, error_unit
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
    v = -999.0_dp
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


!=================================================================================================
! Main
!=================================================================================================
program perturbe_ts
    use iso_fortran_env, only : dp => real64, error_unit
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
    real(dp), allocatable   :: ts_values(:), pvec(:)
    real(dp), allocatable   :: ts_times(:)
    character(len=80), allocatable :: ts_dates(:)
    real(dp)                :: dt

    ! --- 1. Initialization ---
    ! Initialize the random number generator to ensure unique noise for each run
    call random_seed()

    ! --- 2. CLI Arguments Parsing ---
    if (command_argument_count() /= 6) then
        write(error_unit, '(A)') "Usage: ./perturbe_ts <filename> <nrens> <STD> <Tau> <MIN> <MAX>"
        write(error_unit, '(A)') "  Tau: Correlation scale (0.0: white noise, >0.0: red noise)"
        stop 1
    end if

    call get_command_argument(1, filename)

    ! Extract basename by removing the file extension
    dot_pos = index(filename, '.', back=.true.)
    if (dot_pos > 0) then
        basename = filename(1:dot_pos-1)
    else
        basename = trim(filename)
    end if

    ! Parse numeric arguments from command line
    block
        character(len=62) :: arg
        call get_command_argument(2, arg); read(arg, *) nrens
        call get_command_argument(3, arg); read(arg, *) std_obs
        call get_command_argument(4, arg); read(arg, *) tau
        call get_command_argument(5, arg); read(arg, *) vmin
        call get_command_argument(6, arg); read(arg, *) vmax
    end block

    ! --- 3. Read Time-Series (Initialization & Loading) ---
    ! Phase A: Count records to determine array size
    call read_ts(.true., filename, kdim, v_nom, vtime, dstring)

    allocate(ts_values(kdim), ts_times(kdim), ts_dates(kdim))

    ! Phase B: Load data into memory
    do i = 1, kdim
        call read_ts(.false., filename, kdim, ts_values(i), ts_times(i), ts_dates(i))
    end do

    ! Calculate the average time step (dt) for the AR(1) process
    dt = (ts_times(kdim) - ts_times(1)) / real(kdim - 1, dp)
    write(*, '(A,I0,A)') "Loaded ", kdim, " records from " // trim(filename)

    ! --- 4. Generate Noise and Apply Perturbations ---
    allocate(pvec(kdim))    ! Stochastic noise vector for the entire time-series

    ! Loop over each ensemble member to create separate output files
    do ne = 1, nrens
        ! Format member index: 000 (Control), 001...N (Perturbed)
        write(nel, '(I3.3)') ne - 1
        out_file = trim(basename) // "_" // nel // ".dat"

        ! Generate a new stochastic path for this member (Member 1 is control)
        if (ne == 1) then
            pvec = 0.0_dp  ! Control run remains unperturbed
        else
            ! Generate a Red Noise path (AR1) unique to this member
            call perturbe_0d(kdim, pvec, dt, tau)
        end if

        open(unit=30, file=trim(out_file), status='replace', form='formatted')

        ! Apply noise and write the perturbed time-series
        do i = 1, kdim
            ! Apply noise (mean 0, std 1) scaled by the observed standard deviation
            v_nom = ts_values(i) + (pvec(i) * std_obs)

            ! Physical clipping to keep values within the [vmin, vmax] range
            v_nom = max(vmin, min(vmax, v_nom))

            ! Write date string and the perturbed value
            write(30, '(A20,2X,F12.5)') trim(ts_dates(i)), v_nom
        end do

        close(30)
        write(*, '(A,A)') "Created: ", trim(out_file)
    end do

    write(*, '(A)') "Ensemble generation complete."

end program perturbe_ts
