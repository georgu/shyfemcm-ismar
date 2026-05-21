program fem2nc
    use iso_fortran_env, only : dp => real64, sp => real32, error_unit
    use netcdf
    implicit none

    ! --- FEM Variables ---
    character(len=255) :: filename, nc_out
    character(len=10) :: lnrec
    integer :: nrec
    integer :: fem_unit, iformat, nvers, np, lmax, nvar, ntype, nlvddi
    integer :: datetime(2), ierr, t_idx
    real(dp) :: dtime, atime, offset_1970
    real(sp) :: regpar(7)
    integer :: nx, ny
    real(sp) :: x0, y0, dx, dy, flag
    character(len=25) :: x_units, y_units, x_name, y_name

    ! --- Allocatable arrays ---
    real(sp), allocatable :: hlv(:), femdata(:,:,:), hd(:)
    character(len=50), allocatable :: vstring(:)
    character(len=50) :: vvstring
    integer, allocatable :: ilhkv(:)

    ! --- NetCDF IDs ---
    integer :: ncid, x_dim, y_dim, z_dim, t_dim
    integer :: x_varid, y_varid, z_varid, t_varid
    integer, allocatable :: data_varids(:)

    integer :: active_dims(4), dim_count, s_ptr(4), c_ptr(4)
    real(sp), allocatable :: x_coords(:), y_coords(:)
    integer :: i, ii

    if (command_argument_count() /= 3) then
        write(error_unit, '(A)') "Usage: ./fem2nc <input_file.fem> <output_file.nc> <nrec>"
        write(error_unit, '(A)') "nrec: number of records (all = -1)"
        stop 1
    end if

    call get_command_argument(1, filename)
    call get_command_argument(2, nc_out)
    call get_command_argument(3, lnrec)

    read(lnrec, *) nrec

    ! --- 1. DYNAMIC OFFSET CALCULATION ---
    ! exact seconds for 1970-01-01
    ! datetime(1) = YYYYMMDD, datetime(2) = HHMMSS
    datetime(1) = 19700101
    datetime(2) = 0
    call dts_convert_to_atime(datetime, 0.0_dp, offset_1970)

    ! --- 2. INITIAL SCAN ---
    np = 0
    call fem_file_read_open(trim(filename), np, iformat, fem_unit)
    call fem_file_read_params(iformat, fem_unit, dtime, nvers, np, lmax, nvar, ntype, datetime, ierr)
    if (ierr < 0) error stop 'Error: Cannot read FEM file.'

    allocate(hlv(lmax), ilhkv(np), hd(np), vstring(nvar), data_varids(nvar))

    call fem_file_read_2header(iformat, fem_unit, ntype, lmax, hlv, regpar, ierr)
    nx = nint(regpar(1)); ny = nint(regpar(2)); x0 = regpar(3); y0 = regpar(4)
    dx = regpar(5); dy = regpar(6); flag = regpar(7)

    if (abs(x0) <= 360.0 .and. abs(y0) <= 90.0 .and. abs(dx) < 1.0) then
        x_units = "degrees_east";  y_units = "degrees_north"
        x_name  = "longitude";     y_name  = "latitude"
    else
        x_units = "meters";        y_units = "meters"
        x_name  = "x";             y_name  = "y"
    end if

    allocate(femdata(lmax, nx, ny), x_coords(nx), y_coords(ny))
    do i = 1, nx; x_coords(i) = x0 + (i-1)*dx; end do
    do i = 1, ny; y_coords(i) = y0 + (i-1)*dy; end do

    ! --- 3. NETCDF DEFINE MODE ---
    call check( nf90_create(trim(nc_out), NF90_CLOBBER, ncid) )

    dim_count = 0
    if (lmax > 1) then
        call check( nf90_def_dim(ncid, "z", lmax, z_dim) )
        dim_count = dim_count + 1; active_dims(dim_count) = z_dim
    end if
    if (nx > 1) then
        call check( nf90_def_dim(ncid, trim(x_name), nx, x_dim) )
        dim_count = dim_count + 1; active_dims(dim_count) = x_dim
    end if
    if (ny > 1) then
        call check( nf90_def_dim(ncid, trim(y_name), ny, y_dim) )
        dim_count = dim_count + 1; active_dims(dim_count) = y_dim
    end if
    call check( nf90_def_dim(ncid, "time", NF90_UNLIMITED, t_dim) )
    dim_count = dim_count + 1; active_dims(dim_count) = t_dim

    if (nx > 1) then
        call check( nf90_def_var(ncid, trim(x_name), NF90_FLOAT, (/ x_dim /), x_varid) )
        call check( nf90_put_att(ncid, x_varid, "units", trim(x_units)) )
    end if
    if (ny > 1) then
        call check( nf90_def_var(ncid, trim(y_name), NF90_FLOAT, (/ y_dim /), y_varid) )
        call check( nf90_put_att(ncid, y_varid, "units", trim(y_units)) )
    end if
    if (lmax > 1) then
        call check( nf90_def_var(ncid, "z", NF90_FLOAT, (/ z_dim /), z_varid) )
        call check( nf90_put_att(ncid, z_varid, "units", "meters") )
    end if
    
    call check( nf90_def_var(ncid, "time", NF90_DOUBLE, (/ t_dim /), t_varid) )
    call check( nf90_put_att(ncid, t_varid, "units", "seconds since 1970-01-01 00:00:00") )
    call check( nf90_put_att(ncid, t_varid, "calendar", "standard") )

    nlvddi = lmax
    do i = 1, nvar
        call fem_file_read_data(iformat, fem_unit, nvers, np, lmax, vstring(i), ilhkv, hd, nlvddi, femdata, ierr)

	! String correction
	vvstring = vstring(i)
        if (vvstring == "") write(vvstring, '(A,I0)') 'Var_', i
	do ii = 1,len_trim(vvstring)
         if (vvstring(ii:ii) == " ") then
            vvstring(ii:ii) = "_"
         end if
        end do
	vstring(i) = vvstring

        call check( nf90_def_var(ncid, trim(vstring(i)), NF90_FLOAT, active_dims(1:dim_count), data_varids(i)) )
        call check( nf90_put_att(ncid, data_varids(i), "_FillValue", flag) )
    end do
    call check( nf90_enddef(ncid) )

    ! --- 4. DATA WRITING LOOP ---
    if (nx > 1) call check( nf90_put_var(ncid, x_varid, x_coords) )
    if (ny > 1) call check( nf90_put_var(ncid, y_varid, y_coords) )
    if (lmax > 1) call check( nf90_put_var(ncid, z_varid, hlv) )

    rewind(fem_unit)
    t_idx = 0
    read_loop: do
        call fem_file_read_params(iformat, fem_unit, dtime, nvers, np, lmax, nvar, ntype, datetime, ierr)
        if (ierr < 0) exit read_loop 

        t_idx = t_idx + 1
        call dts_convert_to_atime(datetime, dtime, atime)

	write(*,*) "Processing record: ", t_idx, atime
        
        ! Shift using the dynamically calculated offset and round
        atime = anint(atime, kind=dp) - offset_1970
        
        call check( nf90_put_var(ncid, t_varid, (/ atime /), start=(/ t_idx /)) )
        call fem_file_read_2header(iformat, fem_unit, ntype, lmax, hlv, regpar, ierr)

        do i = 1, nvar
            call fem_file_read_data(iformat, fem_unit, nvers, np, lmax, vstring(i), ilhkv, hd, nlvddi, femdata, ierr)
            dim_count = 0
            if (lmax > 1) then; dim_count = dim_count + 1; s_ptr(dim_count) = 1; c_ptr(dim_count) = lmax; end if
            if (nx > 1) then;   dim_count = dim_count + 1; s_ptr(dim_count) = 1; c_ptr(dim_count) = nx;   end if
            if (ny > 1) then;   dim_count = dim_count + 1; s_ptr(dim_count) = 1; c_ptr(dim_count) = ny;   end if
            dim_count = dim_count + 1; s_ptr(dim_count) = t_idx; c_ptr(dim_count) = 1

            call check( nf90_put_var(ncid, data_varids(i), femdata, start=s_ptr(1:dim_count), count=c_ptr(1:dim_count)) )
        end do

	if ((nrec /= -1).and.(t_idx == nrec)) exit

    end do read_loop

    call check( nf90_close(ncid) )
    close(fem_unit)
    write(*,*) "Success! NetCDF created."

contains
    subroutine check(status)
        integer, intent(in) :: status
        if (status /= nf90_noerr) then
            write(error_unit, '(A)') trim(nf90_strerror(status))
            stop "NetCDF Error."
        end if
    end subroutine check
end program fem2nc
