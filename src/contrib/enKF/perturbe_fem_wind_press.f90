
!============================================================================================
subroutine generate_ensemble_perturbations(nx, ny, nrens, &
    dx_deg, dy_deg, y0_lat, rx_km, ry_km, alpha, pmat)

    use iso_fortran_env, only : dp => real64
    use m_sample2D, only : sample2D
    implicit none

    ! Arguments
    integer, intent(in)          :: nx, ny, nrens
    real(dp), intent(in)         :: dx_deg, dy_deg, y0_lat
    real(dp), intent(in)         :: rx_km, ry_km, alpha
    real(dp), intent(inout)      :: pmat(nx, ny, nrens)   ! Persistent AR1 state

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
subroutine make_geo_field(nvar, ne, nrens, nx, ny, lmax, dx_deg, dy_deg, &
                          pmat, var3d, var3d_ens, flag, std_p_pa, y0_lat)
    use iso_fortran_env, only : dp => real64
    implicit none

    ! Arguments
    integer,  intent(in)    :: nvar, ne, nrens, nx, ny, lmax
    real(dp), intent(in)    :: dx_deg, dy_deg, y0_lat, flag, std_p_pa
    real(dp), intent(in)    :: pmat(nx, ny, nrens)      ! Normalized noise
    real(dp), intent(in)    :: var3d(nvar, nx, ny, lmax) 
    real(dp), intent(out)   :: var3d_ens(nvar, nx, ny, lmax)
    
    ! Local variables
    integer  :: ix, iy, i_u, i_v, i_p
    real(dp) :: lat, f_coriolis, dx_m, dy_m, dp_dx, dp_dy, w_speed
    real(dp), parameter :: PI = 3.141592653589793_dp
    real(dp), parameter :: R_EARTH = 6371000.0_dp
    real(dp), parameter :: OMEGA = 7.2921e-5_dp 
    real(dp), parameter :: RHO_AIR = 1.225_dp    
    real(dp), parameter :: V_MAX = 75.0_dp 

    ! Index mapping: 1=U, 2=V, 3=P
    i_u = 1; i_v = 2; i_p = 3
    
    var3d_ens = var3d
    dy_m = R_EARTH * (dy_deg * PI / 180.0_dp)

    do iy = 1, ny
        lat = (y0_lat + (iy-1) * dy_deg) * PI / 180.0_dp
        f_coriolis = 2.0_dp * OMEGA * sin(lat)
        dx_m = R_EARTH * cos(lat) * (dx_deg * PI / 180.0_dp)

        ! Avoid issues near equator
        if (abs(f_coriolis) < 2.0e-5_dp) f_coriolis = 2.0e-5_dp 

        do ix = 1, nx
            ! Land mask check
            if (abs(var3d(i_p, ix, iy, 1) - flag) < 1.0e-3_dp) then
                var3d_ens(:, ix, iy, 1) = flag
                cycle 
            end if

            ! --- A. Pressure Perturbation ---
            ! pmat is scaled directly by std_p_pa (Pascal)
            var3d_ens(i_p, ix, iy, 1) = var3d(i_p, ix, iy, 1) + (pmat(ix, iy, ne) * std_p_pa)

            ! --- B. Gradients of the perturbation (Pa/m) ---
            if (ix == 1) then
                dp_dx = (pmat(ix+1, iy, ne) - pmat(ix, iy, ne)) * std_p_pa / dx_m
            else if (ix == nx) then
                dp_dx = (pmat(ix, iy, ne) - pmat(ix-1, iy, ne)) * std_p_pa / dx_m
            else
                dp_dx = (pmat(ix+1, iy, ne) - pmat(ix-1, iy, ne)) * std_p_pa / (2.0_dp * dx_m)
            end if

            if (iy == 1) then
                dp_dy = (pmat(ix, iy+1, ne) - pmat(ix, iy, ne)) * std_p_pa / dy_m
            else if (iy == ny) then
                dp_dy = (pmat(ix, iy, ne) - pmat(ix, iy-1, ne)) * std_p_pa / dy_m
            else
                dp_dy = (pmat(ix, iy+1, ne) - pmat(ix, iy-1, ne)) * std_p_pa / (2.0_dp * dy_m)
            end if

            ! --- C. Apply Geostrophic Balance ---
            ! u_g = -1/(f*rho) * dP/dy  |  v_g = 1/(f*rho) * dP/dx
            var3d_ens(i_u, ix, iy, 1) = var3d(i_u, ix, iy, 1) - (1.0_dp / (f_coriolis * RHO_AIR)) * dp_dy
            var3d_ens(i_v, ix, iy, 1) = var3d(i_v, ix, iy, 1) + (1.0_dp / (f_coriolis * RHO_AIR)) * dp_dx

            ! Wind speed capping
            w_speed = sqrt(var3d_ens(i_u, ix, iy, 1)**2 + var3d_ens(i_v, ix, iy, 1)**2)
            if (w_speed > V_MAX) then
                var3d_ens(i_u, ix, iy, 1) = var3d_ens(i_u, ix, iy, 1) * (V_MAX / w_speed)
                var3d_ens(i_v, ix, iy, 1) = var3d_ens(i_v, ix, iy, 1) * (V_MAX / w_speed)
            end if
        end do
    end do
end subroutine make_geo_field


!============================================================================================
program perturbe_fem_wind_press
!============================================================================================
    use iso_fortran_env, only : dp => real64, sp => real32, error_unit
    implicit none

    ! --- CLI and File variables ---
    character(len=255)      :: filename, out_file, arg
    integer                 :: pert_type, nrens, dot_pos
    real(dp)                :: std_p, tau

    ! --- FEM format variables ---
    integer                 :: fem_size, iformat, fem_unit, nvers, np, lmax, nvar, ntype
    integer                 :: datetime(2), ierr, i, n
    real(dp)                :: dtime, atime, atime_old
    real(sp)                :: regpar(7) 
    integer                 :: nx, ny
    real(sp)                :: x0, y0, dx, dy, flag
    integer                 :: ne, fid
    character(len=255)      :: basename

    ! --- Allocatable arrays ---
    real(sp), allocatable   :: hlv(:), hd(:), femdata(:,:,:)
    character(len=50), allocatable :: vstring(:)
    integer, allocatable    :: ilhkv(:)

    ! Variables for calculation
    real(dp), allocatable   :: var3d(:,:,:,:), var3d_ens(:,:,:,:)
    real(dp), allocatable   :: pmat(:,:,:) 
    real(dp)                :: alpha, dt_sec

    ! --- 1. CLI Arguments Parsing ---
    if (command_argument_count() /= 5) then
        write(error_unit, '(A)') "Usage: ./perturbe_fem_wind_press <filename> <nrens> <pert_type> <STD_p> <Tau>"
        stop 1
    end if
        
    call get_command_argument(1, filename)
    call get_command_argument(2, arg); read(arg, *) nrens
    call get_command_argument(3, arg); read(arg, *) pert_type
    call get_command_argument(4, arg); read(arg, *) std_p
    call get_command_argument(5, arg); read(arg, *) tau ! Tau in seconds

    dot_pos = index(filename, '.', back=.true.)
    basename = merge(filename(1:dot_pos-1), trim(filename), dot_pos > 0)

    ! Check arguments
    write(*,*) 'Input arguments: ', trim(filename), nrens, pert_type, std_p, tau
    if ((nrens < 2).or.(pert_type > 2).or.(std_p < 0).or.(tau < 0)) error stop 'Bad input arguments, stopping.'
    if (mod(nrens, 2) == 0) error stop 'The ensemble members must be odd with control.'

    ! --- 2. Open Input FEM File ---
    write(*,*) 'Opening input FEM file...'
    call fem_file_read_open(trim(filename), fem_size, iformat, fem_unit)   

    n = 0
    atime_old = -1.0_dp ! Initialize to flag first step
    
    read_loop: do 
        n = n + 1
        
        write(*,*)
        write(*,*) 'Time loop, record: ', n

        call fem_file_read_params(iformat, fem_unit, dtime, nvers, np, lmax, nvar, ntype, datetime, ierr)
        if (ierr < 0) exit read_loop 
        if (lmax > 1) error stop 'In meteo FEM files lmax must be 1.'
        
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
            allocate(pmat(nx, ny, nrens))
            pmat = 0.0_dp ! Zero start for AR1
        end if

        ! --- 3. Read Nominal Variables ---
        do i = 1, nvar
            femdata = flag
            call fem_file_read_data(iformat, fem_unit, nvers, np, lmax, &
                                    vstring(i), ilhkv, hd, 1, femdata, ierr)
            var3d(i, :, :, 1) = femdata(1, :, :)
        end do

        ! --- 4. AR1 Coefficient Calculation ---
        if (atime_old < 0.0_dp .or. tau <= 0.0_dp) then
            alpha = 0.0_dp
        else
            dt_sec = abs(atime - atime_old)
            alpha = exp(-dt_sec / tau)
        end if

        ! --- 5. Perturbation Logic ---
        select case (pert_type)
        case(1) ! Pressure + Geostrophic Wind
            
            if (n == 1) write(*,*) 'Pressure + Geostrophic Wind '

            ! Generate spatial innovation and update AR1 state (pmat)
            ! rx_km/ry_km set to 500km as placeholder
            write(*,*) 'Making random - red noise perturbations...'
            call generate_ensemble_perturbations(nx, ny, nrens-1, &
                        real(dx,dp), real(dy,dp), real(y0,dp), &
                        500.0_dp, 500.0_dp, alpha, pmat)

            do ne = 1, nrens

                ! Construct filename and handle file I/O
                write(arg, '(I3.3)') ne - 1
                out_file = trim(basename) // "_" // trim(arg) // ".fem"
                fid = fem_unit + 10 + ne

                ! Open member file only once at first time step
                if (n == 1) call fem_file_write_open(trim(out_file), iformat, fid)

                call fem_file_write_header(iformat, fid, dtime, nvers, np, lmax, &
                       nvar, ntype, 1, hlv, datetime, regpar)

                if (ne > 1) then
                   ! Compute perturbed field for this member
		   call make_geo_field(nvar, ne-1, nrens-1, nx, ny, lmax, &      ! 1-6
                      real(dx,dp), real(dy,dp), &                                ! 7-8
                      pmat, var3d, var3d_ens, &                                  ! 9-11
                      real(flag,dp), std_p, real(y0,dp))                         ! 12-14
                else
                      var3d_ens = var3d
                end if

                do i = 1, nvar
                    femdata(1, :, :) = real(var3d_ens(i,:,:,1), sp)
                    call fem_file_write_data(iformat, fid, nvers, np, lmax, &
                           vstring(i), ilhkv, hd, 1, femdata)
                end do

            end do

        case default
            error stop 'Unsupported perturbation type.'
        end select

        atime_old = atime
        write(*,'(A,F16.2,A,F6.4)') ' Processed atime: ', atime, ' | Alpha: ', alpha
    end do read_loop

    ! Close all open units
    close(fem_unit)
    do ne = 1, nrens
        close(fem_unit + 10 + ne)
    end do

end program perturbe_fem_wind_press
