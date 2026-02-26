module m_fixsample2D
  !!
  !! Standardize a 2-D field across an ensemble:
  !!   E(nx, ny, nrens) → remove ensemble mean at each (i,j)
  !!   and scale to unit standard deviation (population).
  !!
  !! Implementation uses double precision (real64) via iso_fortran_env.
  !!
  use iso_fortran_env, only: real64, int32
  implicit none

  integer, parameter :: dp = real64
  integer, parameter :: ip = int32

contains

  subroutine fixsample2D(E, nx, ny, nrens)
    !!
    !! Arguments (unchanged from the original):
    !!   E     [inout] : real(dp) :: E(nx, ny, nrens)
    !!   nx,ny [in]    : integer(ip)
    !!   nrens [in]    : integer(ip)
    !!
    !! Notes:
    !!   - Uses population variance (divide by nrens), like the original.
    !!   - Protects against divide-by-zero with tiny(1.0_dp).
    !!   - Standard output prints are kept.
    !!   - The 'sampvar.dat' write block from the original is left in place
    !!     but commented out, per your request.
    !!

    integer(ip),               intent(in)    :: nx, ny, nrens
    real(dp),                  intent(inout) :: E(nx, ny, nrens)

    real(dp), allocatable :: average(:,:)   ! ensemble mean per (i,j)
    real(dp), allocatable :: variance(:,:)  ! population variance per (i,j) (temp)
    real(dp)              :: global_var     ! global average variance (diagnostic)
    integer(ip)           :: i, j, k, istat
    ! integer(ip)           :: unitno       ! used for file I/O in the original

    if (nx <= 0 .or. ny <= 0 .or. nrens <= 0) then
       write(*,*) 'fixsample2D: nx, ny, and nrens must be positive.'
       return
    end if

    allocate(average(nx,ny), variance(nx,ny), stat=istat)
    if (istat /= 0) then
       write(*,*) 'fixsample2D: allocation error; stat=', istat
       return
    end if

    average  = 0.0_dp
    variance = 0.0_dp

    ! 1) Ensemble mean at each grid point (i,j)
    do k = 1, nrens
       average(:,:) = average(:,:) + E(:,:,k)
    end do
    average = average / real(nrens, dp)

    ! 2) Remove the mean (compute anomalies)
    do k = 1, nrens
      E(:,:,k) = E(:,:,k) - average(:,:)
    end do

    ! 3) Accumulate sum of squares for variance
    do k = 1, nrens
      variance(:,:) = variance(:,:) + E(:,:,k)**2
    end do

    ! Global average variance across all points and members
    global_var = sum(variance) / real(nx*ny*nrens, dp)
    write(*,*) 'variance '
    write(*,*) '2D var=', global_var
    write(*,*)

    ! open(10, file='sampvar.dat', position='append')
    ! write(10,'(f12.4)') global_var
    ! close(10)

    ! Convert sum of squares → population variance per grid point
    variance = variance / real(nrens, dp)

    ! 4) Standardize: E(i,j,:) ← E(i,j,:) / sqrt( variance(i,j) )
    !    Use a tiny floor to avoid division by zero for flat fields.
    do j = 1, ny
      do i = 1, nx
        E(i,j,:) = E(i,j,:) / sqrt( max( variance(i,j), tiny(1.0_dp) ) )
      end do
    end do

    deallocate(average, variance, stat=istat)
    if (istat /= 0) write(*,*) 'fixsample2D: deallocation error; stat=', istat

  end subroutine fixsample2D

end module m_fixsample2D
