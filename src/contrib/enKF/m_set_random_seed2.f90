module m_set_random_seed2
  !! High‑quality seed initialization suitable for double‑precision RNG use
  use, intrinsic :: iso_fortran_env, only : int32, int64, real64
  implicit none
  public :: set_random_seed2
  public :: init_random_seed_persistent
contains

  pure function mix64(x) result(y)
    !! SplitMix64 mixing function (fast, good diffusion)
    !! Pure: does not modify its argument.
    integer(int64), intent(in) :: x
    integer(int64) :: y
    integer(int64) :: z
    ! Work on a local copy to preserve purity
    z = x + int(z'9E3779B97F4A7C15', int64)
    y = z
    y = ieor(y, ishft(y, -30)); y = y * int(z'BF58476D1CE4E5B9', int64)
    y = ieor(y, ishft(y, -27)); y = y * int(z'94D049BB133111EB', int64)
    y = ieor(y, ishft(y, -31))
  end function mix64

  subroutine set_random_seed2()
    !! Initialize the Fortran RNG with a strong seed suitable for double‑precision runs.
    integer :: sze, i, istat
    integer, allocatable :: seed(:)
    integer :: vals(8)
    integer :: cnt, rate, cntmax
    integer(int64) :: state, acc

    ! Query required seed size
    call random_seed(size = sze)
    if (sze < 1) sze = 1      ! Robust fallback (rare)
    allocate(seed(sze), stat=istat)
    if (istat /= 0) then
       write(*,*) 'ERROR(set_random_seed2): allocation failed for seed size ', sze
       return
    end if

    ! --- Collect entropy ---
    call date_and_time(values = vals)   ! vals: year,month,day,tz,h,m,s,ms
    call system_clock(cnt, rate, cntmax)

    ! --- Build a 64-bit entropy accumulator ---
    acc = int(vals(1), int64)                          ! year
    acc = ieor(acc, ishft(int(vals(2),int64),  6))     ! month
    acc = ieor(acc, ishft(int(vals(3),int64), 12))     ! day
    acc = ieor(acc, ishft(int(vals(5),int64), 18))     ! hour
    acc = ieor(acc, ishft(int(vals(6),int64), 24))     ! minute
    acc = ieor(acc, ishft(int(vals(7),int64), 30))     ! second
    acc = ieor(acc, ishft(int(vals(8),int64), 36))     ! millisecond
    acc = ieor(acc, ishft(int(cnt   ,int64),  1))      ! cycle counter
    acc = ieor(acc, ishft(int(rate  ,int64), 43))      ! counter rate
    acc = ieor(acc, ishft(int(cntmax,int64), 51))      ! counter range
    if (acc == 0_int64) acc = 1_int64                  ! avoid degenerate state

    ! --- Produce as many 32‑bit seeds as required ---
    state = acc
    do i = 1, sze
       state   = mix64(state)
       seed(i) = int( iand(state, int(huge(0_int32), int64)), int32 )
       seed(i) = abs(seed(i)) + 1    ! avoid zeros in the seed vector
    end do

    call random_seed(put = seed)
    deallocate(seed)
  end subroutine set_random_seed2


  !---------------------------------------------------------------------------
  !> Initialize RNG using a persistent seed file if available; otherwise
  !> create a new strong seed (via set_random_seed2), save it, and use it.
  !>
  !> @param[in] fname   (optional) File name; default 'random_seed.dat'
  !> @param[in] verbose (optional) If .true., print informative messages
  !---------------------------------------------------------------------------
  subroutine init_random_seed_persistent(fname, verbose)
    character(len=*), intent(in), optional :: fname
    logical,          intent(in), optional :: verbose

    character(len=256) :: filename
    integer, allocatable :: seed(:)
    integer :: sze, istat, unit
    logical :: exists
    logical :: chatty

    ! ---- Settings ----
    if (present(fname)) then
       filename = adjustl(fname)
    else
       filename = 'random_seed.dat'
    end if
    chatty = .false.; if (present(verbose)) chatty = verbose

    ! ---- Determine seed vector size and allocate ----
    call random_seed(size = sze)
    if (sze < 1) sze = 1
    allocate(seed(sze), stat=istat)
    if (istat /= 0) then
       if (chatty) write(*,*) 'ERROR(init_random_seed_persistent): cannot allocate seed of size ', sze
       return
    end if

    ! ---- Check if seed file exists ----
    inquire(file=filename, exist=exists)
    unit = 77  ! fixed unit number (Fortran 90); keep unique within your app

    if (exists) then
       ! Try to read seed from existing file
       open(unit=unit, file=filename, status='old', iostat=istat)
       if (istat == 0) then
          read(unit, *, iostat=istat) seed
          close(unit)
          if (istat == 0) then
             call random_seed(put=seed)
             if (chatty) write(*,*) 'Random seed read from file: ', trim(filename)
          else
             ! Corrupted content: fall back to new seed and overwrite
             if (chatty) write(*,*) 'WARNING: corrupted seed file. Creating and saving a new seed.'
             call set_random_seed2()
             call random_seed(get=seed)
             open(unit=unit, file=filename, status='replace', iostat=istat)
             if (istat == 0) then
                write(unit, *) seed
                close(unit)
                if (chatty) write(*,*) 'New random seed written to: ', trim(filename)
             else
                if (chatty) write(*,*) 'ERROR: could not write seed file: ', trim(filename)
             end if
          end if
       else
          ! Could not open for read: create new seed and try writing
          if (chatty) write(*,*) 'WARNING: could not open seed file. Creating a new seed.'
          call set_random_seed2()
          call random_seed(get=seed)
          open(unit=unit, file=filename, status='replace', iostat=istat)
          if (istat == 0) then
             write(unit, *) seed
             close(unit)
             if (chatty) write(*,*) 'New random seed written to: ', trim(filename)
          else
             if (chatty) write(*,*) 'ERROR: could not write seed file: ', trim(filename)
          end if
       end if

    else
       ! No seed file: create a new one
       call set_random_seed2()
       call random_seed(get=seed)
       open(unit=unit, file=filename, status='replace', iostat=istat)
       if (istat == 0) then
          write(unit, *) seed
          close(unit)
          if (chatty) write(*,*) 'Created random seed file: ', trim(filename)
       else
          if (chatty) write(*,*) 'ERROR: could not write seed file: ', trim(filename)
       end if
    end if

    deallocate(seed, stat=istat)
  end subroutine init_random_seed_persistent

  subroutine init_random_seed_constant

        implicit none
        integer n,i
        integer, allocatable :: seed(:)
        integer, save :: icall = 0

        if( icall > 0 ) return
        icall = icall + 1

        call RANDOM_SEED(size=n)

        write(6,*) 'seed size = ',n
        allocate(seed(n))
        do i=1,n
          seed(i) = i*2
        end do
        call RANDOM_SEED(put=seed)

  end subroutine init_random_seed_constant

end module m_set_random_seed2
