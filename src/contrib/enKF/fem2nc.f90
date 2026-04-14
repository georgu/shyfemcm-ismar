program fem2nc
    use iso_fortran_env, only : dp => real64, sp => real32, error_unit
    use netcdf
    implicit none

    ! --- FEM Variables ---
    character(len=255) :: filename, nc_out
    integer :: fem_unit, iformat, fem_size, nvers, np, lmax, nvar, ntype, nlvddi
    integer :: datetime(2), ierr, i, t_idx
    real(dp) :: dtime, atime
    real(sp) :: regpar(7)
    integer :: nx, ny
    real(sp) :: x0, y0, dx, dy, flag

    ! --- Allocatable arrays ---
    real(sp), allocatable :: hlv(:), femdata(:,:,:), hd(:)
    character(len=50), allocatable :: vstring(:)
    integer, allocatable :: ilhkv(:)

    ! --- NetCDF IDs ---
    integer :: ncid, x_dim, y_dim, z_dim, t_dim
    integer :: x_varid, y_varid, z_varid, t_varid
    integer, allocatable :: data_varids(:)

    ! --- Dynamic Dimension Management ---
    integer :: active_dims(4), dim_count, s_ptr(4), c_ptr(4)
    real(sp), allocatable :: x_coords(:), y_coords(:)

    if (command_argument_count() /= 2) then
        write(error_unit, '(A)') "Usage: ./fem2nc <input_file.fem> <output_file.nc>"
        stop 1
    end if

    call get_command_argument(1, filename)
    call get_command_argument(2, nc_out)

    ! 1. INITIAL SCAN
    call fem_file_read_open(trim(filename), fem_size, iformat, fem_unit)
    call fem_file_read_params(iformat, fem_unit, dtime, nvers, np, lmax, nvar, ntype, datetime, ierr)

    nlvddi = lmax

    allocate(hlv(lmax))
    call fem_file_read_2header(iformat, fem_unit, ntype, lmax, hlv, regpar, ierr)

    nx = nint(regpar(1)); ny = nint(regpar(2))
    x0 = regpar(3); y0 = regpar(4); dx = regpar(5); dy = regpar(6); flag = regpar(7)

    allocate(vstring(nvar), ilhkv(np), hd(np), data_varids(nvar), femdata(lmax, nx, ny))
    allocate(x_coords(nx), y_coords(ny))
    do i = 1, nx; x_coords(i) = x0 + (i-1)*dx; end do
    do i = 1, ny; y_coords(i) = y0 + (i-1)*dy; end do

    ! 2. NETCDF DEFINE MODE
    call check( nf90_create(trim(nc_out), NF90_CLOBBER, ncid) )

    ! Define Dimensions
    dim_count = 0
    if (nx > 1) then
        call check( nf90_def_dim(ncid, "x", nx, x_dim) )
        dim_count = dim_count + 1; active_dims(dim_count) = x_dim
    end if
    if (ny > 1) then
        call check( nf90_def_dim(ncid, "y", ny, y_dim) )
        dim_count = dim_count + 1; active_dims(dim_count) = y_dim
    end if
    if (lmax > 1) then
        call check( nf90_def_dim(ncid, "z", lmax, z_dim) )
        dim_count = dim_count + 1; active_dims(dim_count) = z_dim
    end if
    call check( nf90_def_dim(ncid, "time", NF90_UNLIMITED, t_dim) )
    dim_count = dim_count + 1; active_dims(dim_count) = t_dim

    ! Define Coordinate Variables with CF attributes
    if (nx > 1) then
        call check( nf90_def_var(ncid, "x", NF90_FLOAT, (/ x_dim /), x_varid) )
        !call check( nf90_put_att(ncid, x_varid, "units", "meters") )
        call check( nf90_put_att(ncid, x_varid, "axis", "X") )
    end if
    if (ny > 1) then
        call check( nf90_def_var(ncid, "y", NF90_FLOAT, (/ y_dim /), y_varid) )
        !call check( nf90_put_att(ncid, y_varid, "units", "meters") )
        call check( nf90_put_att(ncid, y_varid, "axis", "Y") )
    end if
    if (lmax > 1) then
        call check( nf90_def_var(ncid, "z", NF90_FLOAT, (/ z_dim /), z_varid) )
        call check( nf90_put_att(ncid, z_varid, "units", "level") )
    end if
    call check( nf90_def_var(ncid, "time", NF90_DOUBLE, (/ t_dim /), t_varid) )
    call check( nf90_put_att(ncid, t_varid, "units", "seconds since 0001-01-01 00:00:00") )
    call check( nf90_put_att(ncid, t_varid, "calendar", "proleptic_gregorian") )

    ! Define Data Variables
    close(fem_unit)
    call fem_file_read_open(trim(filename), fem_size, iformat, fem_unit)
    call fem_file_read_params(iformat, fem_unit, dtime, nvers, np, lmax, nvar, ntype, datetime, ierr)
    call fem_file_read_2header(iformat, fem_unit, ntype, lmax, hlv, regpar, ierr)

    do i = 1, nvar
        call fem_file_read_data(iformat, fem_unit, nvers, np, lmax, vstring(i), ilhkv, hd, nlvddi, femdata, ierr)
	if (vstring(i) == "") then
	   write(vstring(i), '(A,I0)') 'Var_', i
	end if
        call check( nf90_def_var(ncid, trim(vstring(i)), NF90_FLOAT, active_dims(1:dim_count), data_varids(i)) )
        call check( nf90_put_att(ncid, data_varids(i), "_FillValue", flag) )
    end do

    call check( nf90_enddef(ncid) )

    ! 3. DATA WRITING
    if (nx > 1) call check( nf90_put_var(ncid, x_varid, x_coords) )
    if (ny > 1) call check( nf90_put_var(ncid, y_varid, y_coords) )
    if (lmax > 1) call check( nf90_put_var(ncid, z_varid, hlv) )

    close(fem_unit)
    call fem_file_read_open(trim(filename), fem_size, iformat, fem_unit)

    t_idx = 0
    read_loop: do
        call fem_file_read_params(iformat, fem_unit, dtime, nvers, np, lmax, nvar, ntype, datetime, ierr)
        if (ierr < 0) exit read_loop

        t_idx = t_idx + 1
        call dts_convert_to_atime(datetime, dtime, atime)
        call check( nf90_put_var(ncid, t_varid, atime, start=(/ t_idx /)) )
        call fem_file_read_2header(iformat, fem_unit, ntype, lmax, hlv, regpar, ierr)

        do i = 1, nvar
                call fem_file_read_data(iformat, fem_unit, nvers, np, lmax, vstring(i), ilhkv, hd, nlvddi, femdata, ierr)

                dim_count = 0
                if (nx > 1) then; dim_count = dim_count + 1; s_ptr(dim_count) = 1; c_ptr(dim_count) = nx; end if
                if (ny > 1) then; dim_count = dim_count + 1; s_ptr(dim_count) = 1; c_ptr(dim_count) = ny; end if
                if (lmax >= 1) then; dim_count = dim_count + 1; s_ptr(dim_count) = 1; c_ptr(dim_count) = lmax; end if

                dim_count = dim_count + 1
                s_ptr(dim_count) = t_idx; c_ptr(dim_count) = 1

                call check( nf90_put_var(ncid, data_varids(i), real(femdata, dp), &
                            start=s_ptr(1:dim_count), count=c_ptr(1:dim_count)) )
                !call check( nf90_put_var(ncid, data_varids(i), real(reshape(femdata, (/ nx, ny, lmax /), order=(/ 2, 3, 1 /)), dp), &
                !            start=s_ptr(1:dim_count), count=c_ptr(1:dim_count)) )
        end do
    end do read_loop

    call check( nf90_close(ncid) )
    close(fem_unit)
    write(*,*) "Success: NetCDF created. Steps:", t_idx

contains
    subroutine check(status)
        integer, intent(in) :: status
        if (status /= nf90_noerr) then
            write(error_unit, '(A)') trim(nf90_strerror(status))
            stop "NetCDF Error."
        end if
    end subroutine check
end program fem2nc
