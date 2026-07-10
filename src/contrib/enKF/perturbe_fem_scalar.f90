!============================================================================================
subroutine generate_ensemble_perturbations(nx, ny, nrens, &
    dx_deg, dy_deg, y0_lat, rx_km, ry_km, alpha, pmat)

    use iso_fortran_env, only : dp => real64
    use m_sample2D, only : sample2D
    implicit none

    ! Arguments
    integer, intent(in)            :: nx, ny, nrens
    real(dp), intent(in)           :: dx_deg, dy_deg, y0_lat
    real(dp), intent(in)           :: rx_km, ry_km, alpha
    real(dp), intent(inout)        :: pmat(nx, ny, nrens)   ! Persistent AR1 state

    ! Local variables
    real(dp) :: dx_m, dy_m, rx_m, ry_m, lat_rad
    real(dp), allocatable :: amat(:,:,:)
    real(dp), parameter   :: PI = 3.141592653589793_dp
    real(dp), parameter   :: R_EARTH = 6371000.0_dp
    integer  :: nre_factor

    ! --- 1. Convert Degrees to Meters ---
    lat_rad = y0_lat * PI / 180.0_dp
    dx_m = R_EARTH * cos(lat_rad) * (dx_deg * PI / 180.0_dp)
    dy_m = R_EARTH * (dy_deg * PI / 180.0_dp)

    ! Convert correlation scales from km to meters
    rx_m = rx_km * 1000.0_dp
    ry_m = ry_km * 1000.0_dp

    ! --- 2. Generate Spatial White Noise (Innovation) ---
    ! nre=5 is a typical value for the oversized ensemble in sample2D
    nre_factor = 5
    allocate(amat(nx, ny, nrens))

    ! Calling your module's routine
    call sample2D(amat, nx, ny, nrens, nre_factor, dx_m, dy_m, rx_m, ry_m, &
                  0.0_dp, .true., .false.)

    ! --- 3. Apply AR1 Temporal Evolution (Red Noise) ---
    ! pmat is the noise field
    ! pmat(t) = alpha * pmat(t-1) + sqrt(1 - alpha^2) * NewNoise
    pmat = alpha * pmat + sqrt(1.0_dp - alpha**2) * amat

    deallocate(amat)

end subroutine generate_ensemble_perturbations

!============================================================================================
subroutine generate_random_pert(nrens, nx, ny, pmat, alpha)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    
    ! Arguments
    integer, intent(in)           :: nrens, nx, ny
    real(dp), intent(inout)       :: pmat(nx, ny, nrens)
    real(dp), intent(in)          :: alpha
    
    ! Local variables
    real(dp) :: temp_vec(nrens)
    real(dp) :: amat(nx, ny, nrens)
    real(dp) :: u1, u2, r, factor, v_mean
    integer :: i, j, k
    real(dp), parameter :: limit = 4.0_dp ! Clipping threshold

    ! 1. Generate Gaussian random vector (Mean 0, Std 1) using Box-Muller
    k = 1
    do while (k <= nrens)
        do
            call random_number(u1)
            call random_number(u2)
            u1 = 2.0_dp * u1 - 1.0_dp
            u2 = 2.0_dp * u2 - 1.0_dp
            r = u1**2 + u2**2
            if (r > 0.0_dp .and. r < 1.0_dp) exit
        end do
        
        factor = sqrt(-2.0_dp * log(r) / r)
        temp_vec(k) = u1 * factor
        if (k + 1 <= nrens) temp_vec(k+1) = u2 * factor
        k = k + 2
    end do

    ! 2. Optional: Force zero mean to prevent bias in EnKF
    v_mean = sum(temp_vec) / real(nrens, dp)
    temp_vec = temp_vec - v_mean

    ! 3. Control extreme values (Clipping)
    where (temp_vec > limit)
        temp_vec = limit
    elsewhere (temp_vec < -limit)
        temp_vec = -limit
    end where

    ! 4. Distribute vector to the matrix
    do j = 1, ny
        do i = 1, nx
            amat(i, j, :) = temp_vec
        end do
    end do

    ! --- 3. Apply AR1 Temporal Evolution (Red Noise) ---
    ! pmat is the noise field
    ! pmat(t) = alpha * pmat(t-1) + sqrt(1 - alpha^2) * NewNoise
    pmat = alpha * pmat + sqrt(1.0_dp - alpha**2) * amat


end subroutine generate_random_pert


!============================================================================================
subroutine make_pert_field(nvar, ne, nrens, nx, ny, lmax, hlv, var3d, var3d_ens, &
    std_r, vmin, vmax, pmat)

    use iso_fortran_env, only : dp => real64, sp => real32
    implicit none

    ! Arguments
    integer,  intent(in)    :: nvar, ne, nrens, nx, ny, lmax
    real(dp), intent(in)    :: pmat(nx, ny, nrens)
    real(sp), intent(in)    :: hlv(lmax)
    real(dp), intent(in)    :: var3d(nvar, nx, ny, lmax)
    real(dp), intent(out)   :: var3d_ens(nvar, nx, ny, lmax)
    real(dp), intent(in)    :: std_r, vmin, vmax

    ! Local variables
    real(dp) :: scaling_factor
    integer  :: l, i

    ! 1. Initialize output with background values (entire 4D volume)
    var3d_ens = var3d 

    do l = 1, lmax

       ! 2. Calculate vertical scaling (hyperbolic dampening: 1/z)
       ! Protect against division by zero if hlv(l) is not yet set
       if (hlv(l) > 0.0_sp) then
           scaling_factor = real(hlv(1), dp) / real(hlv(l), dp)
       else
           scaling_factor = 1.0_dp
       end if

      do i = 1, nvar  ! Only one STD for all vars!!!

       ! 3. Apply perturbation to the first variable (index 1) at level l
       ! The scaling_factor reduces the std_r as depth increases
       var3d_ens(i,:,:,l) = var3d(i,:,:,l) + ( std_r * pmat(:,:,ne) * scaling_factor )

       ! 4. Final clipping to physical range [vmin, vmax]
       var3d_ens(i,:,:,l) = max(vmin, min(vmax, var3d_ens(i,:,:,l)))

      end do
    end do

end subroutine make_pert_field

!============================================================================================
program perturbe_fem_scalar
!============================================================================================
    use iso_fortran_env, only : dp => real64, sp => real32, error_unit
    use m_set_random_seed2
    implicit none

    ! --- CLI and File variables ---
    character(len=255)                :: filename, arg
    character(len=255), allocatable   :: out_file(:)
    integer                           :: nrens, dot_pos, ptype
    real(dp)                          :: std_r, tau, vmin, vmax

    ! --- FEM format variables ---
    integer                 :: fem_size, iformat, fem_unit, nvers, np, lmax, nvar, ntype, nlvddi
    integer                 :: datetime(2), ierr
    real(dp)                :: dtime, atime, atime_old
    real(sp)                :: regpar(7) 
    integer                 :: nx, ny
    real(sp)                :: x0, y0, dx, dy, flag
    integer, allocatable    :: fid(:)
    character(len=255)      :: basename

    ! --- Allocatable arrays ---
    real(sp), allocatable   :: hlv(:), hd(:), femdata(:,:,:)
    character(len=50), allocatable :: vstring(:)
    integer, allocatable    :: ilhkv(:)

    ! Variables for calculation
    real(dp), allocatable   :: var3d(:,:,:,:), var3d_ens(:,:,:,:)
    real(dp), allocatable   :: pmat(:,:,:) 
    real(dp)                :: alpha, dt_sec

    integer                 :: i, j, k, l, n, ne

    ! CLI Arguments Parsing ---
    if (command_argument_count() /= 7) then
	write(error_unit, '(A)')
        write(error_unit, '(A)') "Usage: ./perturbe_fem_scalar <filename> <nrens> <ptype> <STD_r> <vmin> <vmax> <Tau>"
	write(error_unit, '(A)')
	write(error_unit, '(A)') "filename: FEM file"
	write(error_unit, '(A)') "nrens: number of ensemble members, control included"
	write(error_unit, '(A)') "ptype: type of perturbation, scalar (1), 2D pseudo-Gaussian (2)"
	write(error_unit, '(A)') "STD_r: error (standard deviation) for the scalar"
	write(error_unit, '(A)') "vmin: minimum value for the scalar"
	write(error_unit, '(A)') "vmax: maximum value for the scalar"
	write(error_unit, '(A)') "Tau: Decay time (s) for the red noise"
	write(error_unit, '(A)')
        stop 1
    end if
        
    call get_command_argument(1, filename)
    call get_command_argument(2, arg); read(arg, *) nrens
    call get_command_argument(3, arg); read(arg, *) ptype
    call get_command_argument(4, arg); read(arg, *) std_r
    call get_command_argument(5, arg); read(arg, *) vmin
    call get_command_argument(6, arg); read(arg, *) vmax
    call get_command_argument(7, arg); read(arg, *) tau ! Tau in seconds

    dot_pos = index(filename, '.', back=.true.)
    basename = merge(filename(1:dot_pos-1), trim(filename), dot_pos > 0)

    ! Check arguments
    write(*,*) 'Input arguments: ', trim(filename), nrens, std_r, vmin, vmax, tau
    if ((nrens < 2).or.(std_r < 0).or.(tau < 0)) &
	    error stop 'Bad input arguments, stopping.'
    if (mod(nrens, 2) == 0) error stop 'The ensemble members must be odd with control.'

    call init_random_seed_persistent('random_seed.dat', .true.)

    ! Open Input FEM File ---
    write(*,*) 'Opening input FEM file...'
    call fem_file_read_open(trim(filename), fem_size, iformat, fem_unit)   

    ! Allocation write units ---
    allocate(fid(nrens), out_file(nrens))
    do ne = 1, nrens
       fid(ne) = fem_unit + ne + 10
       write(arg, '(I3.3)') ne - 1
       out_file(ne) = trim(basename) // "_" // trim(arg) // ".fem"
    end do

    n = 0
    atime_old = -1.0_dp ! Initialize to flag first step
    
    read_loop: do 
        n = n + 1
        
        write(*,*)
        write(*,*) 'Time loop, record: ', n

        call fem_file_read_params(iformat, fem_unit, dtime, nvers, np, lmax, nvar, ntype, datetime, ierr)
        if (ierr < 0) exit read_loop 
        
        call dts_convert_to_atime(datetime, dtime, atime) ! atime in seconds

        ! --- Lazy Allocation ---
        if (.not. allocated(hlv))     allocate(hlv(lmax))
        if (.not. allocated(ilhkv))   allocate(ilhkv(np))
        if (.not. allocated(hd))      allocate(hd(np))
        if (.not. allocated(vstring)) allocate(vstring(nvar))

        call fem_file_read_2header(iformat, fem_unit, ntype, lmax, hlv, regpar, ierr)
        nx = nint(regpar(1)); ny = nint(regpar(2)); x0 = regpar(3)
        y0 = regpar(4); dx = regpar(5); dy = regpar(6); flag = regpar(7)

        if (n == 1) write(*,*) 'FEM parameters: ', nx, ny, x0, y0, dx, dy, flag

        if (.not. allocated(femdata))   allocate(femdata(lmax, nx, ny))
        if (.not. allocated(var3d))     allocate(var3d(nvar, nx, ny, lmax))
        if (.not. allocated(var3d_ens)) allocate(var3d_ens(nvar, nx, ny, lmax))
        if (.not. allocated(pmat)) then
            allocate(pmat(nx, ny, nrens-1))
            pmat = 0.0_dp ! Zero start for AR1
        end if

        ! Read Nominal Variables ---
        do i = 1, nvar
            femdata = flag
	    nlvddi = lmax
            call fem_file_read_data(iformat, fem_unit, nvers, np, lmax, &
                                    vstring(i), ilhkv, hd, nlvddi, femdata, ierr)
	    do concurrent (l=1:lmax, j=1:nx, k=1:ny)
               var3d(i, j, k, l) = real(femdata(l, j, k), dp)
            end do
        end do

	! Variable check
	if (n == 1) write(*,*) 'FEM variable: ', vstring(1)

        ! --- 4. AR1 Coefficient Calculation ---
        if (atime_old < 0.0_dp .or. tau <= 0.0_dp) then
            alpha = 0.0_dp
        else
            dt_sec = abs(atime - atime_old)
            alpha = exp(-dt_sec / tau)
        end if

        ! Perturbation Logic ---
        select case (ptype)
	! ------------------------------------
        case(1) ! Scalar. 
            
            if (n == 1) write(*,*) 'Scalar field'

            call generate_random_pert(nrens-1, nx, ny, pmat, alpha)

            do ne = 1, nrens

                ! Open member file only once at first time step
                if (n == 1) call fem_file_write_open(trim(out_file(ne)), iformat, fid(ne))

                call fem_file_write_header(iformat, fid(ne), dtime, nvers, np, lmax, &
                       nvar, ntype, nlvddi, hlv, datetime, regpar)

                if (ne > 1) then
                   ! Compute perturbed field for this member
		   call make_pert_field(nvar, ne-1, nrens-1, nx, ny, lmax, hlv, var3d, var3d_ens, &
                        std_r, vmin, vmax, pmat)
                else
                      var3d_ens = var3d
                end if

                do i = 1, nvar
	            do concurrent (l=1:lmax, j=1:nx, k=1:ny)
                       femdata(l, j, k) = real(var3d_ens(i, j, k, l), sp)
                    end do
                    call fem_file_write_data(iformat, fid(ne), nvers, np, lmax, &
                           vstring(i), ilhkv, hd, nlvddi, femdata)
                end do

            end do

	! ------------------------------------
        case(2) ! 2D pseudo-Gaussian field. 
            
            if (n == 1) write(*,*) '2D pseudo-Gaussian field. One for all variables'

            ! Generate spatial innovation and update AR1 state (pmat)
            ! rx_km/ry_km set to 500km as placeholder
            write(*,*) 'Making random - red noise perturbations...'
            call generate_ensemble_perturbations(nx, ny, nrens-1, &
                        real(dx,dp), real(dy,dp), real(y0,dp), &
                        500.0_dp, 500.0_dp, alpha, pmat)

            do ne = 1, nrens

                ! Open member file only once at first time step
                if (n == 1) call fem_file_write_open(trim(out_file(ne)), iformat, fid(ne))

                call fem_file_write_header(iformat, fid(ne), dtime, nvers, np, lmax, &
                       nvar, ntype, nlvddi, hlv, datetime, regpar)

                if (ne > 1) then
                   ! Compute perturbed field for this member
		   call make_pert_field(nvar, ne-1, nrens-1, nx, ny, lmax, hlv, var3d, var3d_ens, &
                        std_r, vmin, vmax, pmat)
                else
                      var3d_ens = var3d
                end if

                do i = 1, nvar
	            do concurrent (l=1:lmax, j=1:nx, k=1:ny)
                       femdata(l, j, k) = real(var3d_ens(i, j, k, l), sp)
                    end do
                    call fem_file_write_data(iformat, fid(ne), nvers, np, lmax, &
                           vstring(i), ilhkv, hd, nlvddi, femdata)
                end do

            end do

	! ------------------------------------
        case default
            error stop 'Unsupported perturbation type.'

        end select

        atime_old = atime
        write(*,'(A,F16.2,A,F6.4)') ' Processed atime: ', atime, ' | Alpha: ', alpha
    end do read_loop

    ! Close all open units
    close(fem_unit)
    do ne = 1, nrens
        close(fid(ne))
    end do

end program perturbe_fem_scalar
