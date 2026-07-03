module m_set_random_seed2
  !! High‑quality seed initialization suitable for double‑precision RNG use
  !! Optimized strictly for Pure OpenMP thread execution architectures.
  use, intrinsic :: iso_fortran_env, only : int32, int64, real64
  implicit none
  public :: set_random_seed2
  public :: init_random_seed_persistent
contains

  pure function mix64(x) result(y)
    integer(int64), intent(in) :: x
    integer(int64) :: y
    integer(int64) :: z
    z = x + int(z'9E3779B97F4A7C15', int64)
    y = z
    y = ieor(y, ishft(y, -30)); y = y * int(z'BF58476D1CE4E5B9', int64)
    y = ieor(y, ishft(y, -27)); y = y * int(z'94D049BB133111EB', int64)
    y = ieor(y, ishft(y, -31))
  end function mix64

  subroutine set_random_seed2()
    !! Strong serial seed generator utilizing both time and CPU cycle entropy
    integer :: sze, i, istat
    integer, allocatable :: seed(:)
    integer :: vals(8)
    integer :: cnt, rate, cntmax
    integer(int64) :: state, acc

    call random_seed(size = sze)
    if (sze < 1) sze = 1      
    allocate(seed(sze), stat=istat)
    if (istat /= 0) return

    call date_and_time(values = vals)   
    call system_clock(cnt, rate, cntmax)

    acc = int(vals(1), int64)                          
    acc = ieor(acc, ishft(int(vals(2),int64),  6))     
    acc = ieor(acc, ishft(int(vals(3),int64), 12))     
    acc = ieor(acc, ishft(int(vals(5),int64), 18))     
    acc = ieor(acc, ishft(int(vals(6),int64), 24))     
    acc = ieor(acc, ishft(int(vals(7),int64), 30))     
    acc = ieor(acc, ishft(int(vals(8),int64), 36))     
    acc = ieor(acc, ishft(int(cnt   ,int64),  1))      
    acc = ieor(acc, ishft(int(rate  ,int64), 43))      
    acc = ieor(acc, ishft(int(cntmax,int64), 51))      
    if (acc == 0_int64) acc = 1_int64                  

    state = acc
    do i = 1, sze
       state   = mix64(state)
       seed(i) = int( iand(state, int(huge(0_int32), int64)), int32 )
       seed(i) = abs(seed(i)) + 1    
    end do

    call random_seed(put = seed)
    deallocate(seed)
  end subroutine set_random_seed2

  subroutine init_random_seed_persistent(fname, verbose)
    character(len=*), intent(in), optional :: fname
    logical,          intent(in), optional :: verbose

    character(len=256) :: filename
    integer, allocatable :: seed(:)
    integer :: sze, istat, u
    logical :: exists
    logical :: chatty

    if (present(fname)) then
       filename = adjustl(fname)
    else
       filename = 'random_seed.dat'
    end if
    chatty = .false.; if (present(verbose)) chatty = verbose

    call random_seed(size = sze)
    if (sze < 1) sze = 1
    allocate(seed(sze), stat=istat)
    if (istat /= 0) return

    ! In OpenMP, this check runs inside the single master execution thread
    inquire(file=filename, exist=exists)

    if (exists) then
       open(newunit=u, file=filename, status='old', iostat=istat)
       if (istat == 0) then
          read(u, *, iostat=istat) seed
          close(u)
          if (istat == 0) then
             call random_seed(put=seed)
             if (chatty) write(*,*) 'Random seed successfully loaded from single file: ', trim(filename)
          else
             if (chatty) write(*,*) 'WARNING: corrupted seed file. Regenerating.'
             call set_random_seed2()
             call random_seed(get=seed)
             open(newunit=u, file=filename, status='replace', iostat=istat)
             if (istat == 0) then
                write(u, *) seed
                close(u)
             end if
          end if
       else
          call set_random_seed2()
          call random_seed(get=seed)
          open(newunit=u, file=filename, status='replace', iostat=istat)
          if (istat == 0) then
             write(u, *) seed
             close(u)
          end if
       end if

    else
       ! Fresh initialization
       call set_random_seed2()
       call random_seed(get=seed)
       open(newunit=u, file=filename, status='replace', iostat=istat)
       if (istat == 0) then
          write(u, *) seed
          close(u)
          if (chatty) write(*,*) 'Created standalone high-entropy random seed registry: ', trim(filename)
       end if
    end if

    deallocate(seed, stat=istat)
  end subroutine init_random_seed_persistent

end module m_set_random_seed2
