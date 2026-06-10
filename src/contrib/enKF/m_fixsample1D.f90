module m_fixsample1D
  !
  ! Row-wise mean removal and standardization across m samples
  ! using double precision (real64) from iso_fortran_env.
  !
  use iso_fortran_env, only: real64, int32
  implicit none

  ! Define double precision alias
  integer, parameter :: dp  = real64
  integer, parameter :: ip  = int32

contains

  subroutine fixsample1D(E, n, m, verbose)
    !!
    !! Standardizes each row of E(n, m) across the m samples:
    !!   1. Compute row-wise means
    !!   2. Remove means (anomalies)
    !!   3. Compute row-wise variance (population)
    !!   4. Scale each row to unit standard deviation
    !!
    !! Arguments:
    !!   E(n,m) [inout] : real(dp)   → data matrix
    !!   n     [in]     : integer(ip)
    !!   m     [in]     : integer(ip)
    !!   verbose[in,opt]: logical    → diagnostics
    !!

    integer(ip),               intent(in)    :: n, m
    real(dp),                  intent(inout) :: E(n, m)
    logical,          optional,intent(in)    :: verbose

    real(dp), allocatable :: avg(:)
    real(dp), allocatable :: var_row(:)
    real(dp)              :: global_var
    real(dp)              :: inv_std_i
    integer(ip)           :: i, j, istat
    logical               :: do_print

    do_print = .false.
    if (present(verbose)) do_print = verbose

    if (n <= 0 .or. m <= 0) then
       if (do_print) write(*,*) 'fixsample1D: n and m must be positive.'
       return
    end if

    ! Allocate workspace
    allocate(avg(n), var_row(n), stat=istat)
    if (istat /= 0) then
       if (do_print) write(*,*) 'Allocation error, stat=', istat
       return
    end if

    avg     = 0.0_dp
    var_row = 0.0_dp

    ! ---- 1) Row mean
    do j = 1, m
       avg(:) = avg(:) + E(:, j)
    end do
    avg = avg / real(m, dp)

    ! ---- 2) Remove mean (center anomalies)
    do j = 1, m
       E(:, j) = E(:, j) - avg(:)
    end do

    ! ---- 3) Row-wise variance (population)
    do j = 1, m
       var_row(:) = var_row(:) + E(:, j)**2
    end do

    global_var = sum(var_row) / real(n*m, dp)
    if (do_print) write(*,*) 'fixsample1D: average variance = ', global_var

    var_row = var_row / real(m, dp)

    ! ---- 4) Row-wise standardization
    do i = 1, n
       inv_std_i = 1.0_dp / sqrt( max(var_row(i), tiny(1.0_dp)) )
       E(i, :)   = inv_std_i * E(i, :)
    end do

    ! ---- Cleanup
    deallocate(avg, var_row, stat=istat)
    if (istat /= 0 .and. do_print) write(*,*) 'Deallocation error, stat=', istat

  end subroutine fixsample1D

end module m_fixsample1D
