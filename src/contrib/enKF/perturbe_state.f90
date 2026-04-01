module mod_perturbations
    use iso_fortran_env, only : sp => real32, dp => real64
    implicit none

    ! Buffers for perturbed states (to be saved as real*4)
    real(sp), allocatable, public :: zz_out(:,:)
    real(sp), allocatable, public :: tt_out(:,:,:), ss_out(:,:,:)

contains

    subroutine allocate_pert_fields(dimk, nnlv, nrens, baroclinic)
        integer, intent(in) :: dimk, nnlv, nrens
        logical, intent(in) :: baroclinic

        if (allocated(zz_out)) deallocate(zz_out)
        allocate(zz_out(dimk, nrens))

        if (baroclinic) then
            if (allocated(tt_out)) deallocate(tt_out)
            if (allocated(ss_out)) deallocate(ss_out)
            allocate(tt_out(nnlv, dimk, nrens))
            allocate(ss_out(nnlv, dimk, nrens))
        end if
    end subroutine

end module mod_perturbations

! ======================================================================
! SUBROUTINES
! ======================================================================
! ======================================================================

! ======================================================================
subroutine init_shyfem(basinf, nnlv)
   use basin
   use levels
   use mod_geom_dynamic
   use mod_hydro
   use mod_hydro_vel
   use mod_ts
   use mod_restart
   use shympi
   use mod_gotm_aux
   implicit none
      
   character(len=*), intent(in) :: basinf
   integer, intent(in) :: nnlv
   integer :: ios
      
   open(21,file=basinf,status="old",form="unformatted",iostat=ios)
   if (ios/=0) error stop "Cannot open basin"
   call basin_read_by_unit(21)
   close(21)
      
   nlv = nnlv
   nlvdi = nnlv

   call mod_geom_dynamic_init(nkn,nel)
   call mod_hydro_init(nkn,nel,nlv)
   call mod_hydro_vel_init(nkn,nel,nlv)
   call mod_ts_init(nkn,nlv)
   call levels_init(nkn,nel,nlv)
   call mod_gotm_aux_init(nkn, nlv)
   call shympi_set_hlv(nlv, hlv)
   call shympi_init(.false.)

   nkn_global = nkn
   nel_global = nel
   nlv_global = nlv
end subroutine init_shyfem

! ======================================================================
subroutine read_rst(rstname, atimea)
   use iso_fortran_env, only : dp=>real64
   use mod_restart
   use levels, only : nlvdi, nlv, hlv, ilhv, ilhkv
   use shympi
   implicit none

   character(len=*), intent(in) :: rstname
   real(dp), intent(in) :: atimea

   integer :: ierr, iflag, ios
   real(dp) :: atimef
   integer, save :: icall = 0
   real*4 :: ibarcl4, iconz4, imerc4, iturb4, iwvert_rst4, ieco_rst4, zero4
   zero4 = 0.0

   open(24,file=trim(rstname),status='old',form='unformatted',action='read',iostat=ios)
   if (ios /= 0) error stop "read_rst: cannot open"

   do
      call rst_read_record(24, atimef, iflag, ierr)
      if (ierr /= 0) then
         close(24)
         error stop "read_rst: time not found"
      end if
      if (abs(atimef - atimea) < epsilon(atimea)) exit
   end do
   close(24)

   if (icall==0) then
      hlv        = hlvrst
      hlv_global = hlvrst
      ilhv       = ilhrst
      ilhkv      = ilhkrst

      ibarcl4     = ibarcl_rst
      iwvert_rst4 = iwvert_rst
      ieco_rst4   = ieco_rst
      iconz4      = iconz_rst
      imerc4      = imerc_rst
      iturb4      = iturb_rst

      call addpar("ibarcl",ibarcl4)
      call addpar("iconz", iconz4)
      call addpar("iwvert",iwvert_rst4)
      call addpar("ieco"  ,ieco_rst4)
      call addpar("ibio"  ,zero4)
      call addpar("ibfm"  ,zero4)
      call addpar("imerc" ,imerc4)
      call addpar("iturb" ,iturb4)

      call daddpar("date",0.0_dp)
      call daddpar("time",0.0_dp)

      write(*,*) 'SHYFEM flags from restart:'
      write(*,*) 'nvers = ', nvers_rst
      write(*,*) 'nvmax = ', nvmax
      write(*,*) 'ibarcl = ', ibarcl_rst
      write(*,*) 'iconz  = ', iconz_rst
      write(*,*) 'iwvert = ', iwvert_rst
      write(*,*) 'ieco   = ', ieco_rst
      write(*,*) 'imerc  = ', imerc_rst
      write(*,*) 'iturb  = ', iturb_rst
      write(*,*) 'nlv   = ', nlv
      !write(*,*) 'hlvrst = ', hlvrst(1:nlv)
      write(*,*) 'hlv    = ', hlv(1:nlv)

   end if

   icall = icall + 1
end subroutine read_rst

! ======================================================================
subroutine make_perturbations(nrens, nnlv, err_z, err_t, err_s)
    use iso_fortran_env, only : dp => real64, sp => real32
    use basin
    use levels, only : hlv
    use mod_hydro
    use mod_hydro_vel
    use mod_ts
    use mod_restart,   only : ibarcl_rst
    use mod_perturbations ! Access to zz_out, tt_out, ss_out

    implicit none

    integer,  intent(in) :: nrens, nnlv
    real(dp), intent(in) :: err_z, err_t, err_s

    integer               :: dimk, ne, nl
    real(dp)              :: weight
    real(dp), allocatable :: kpert(:,:)

    dimk = size(znv)

    ! Prepare the output module storage
    call allocate_pert_fields(dimk, nnlv, nrens, (ibarcl_rst /= 0))

    write(*,*) '===================================================================='
    write(*,*) 'Perturbing Z...'
    write(*,*) '===================================================================='
    allocate(kpert(dimk, nrens-1))
    call make_2Dpert(kpert, dimk, nrens-1)

    ! --- Z Perturbation ---
    zz_out(:, 1) = real(znv(:), sp)
    do ne = 2, nrens
       ! Calculation in DP, storage in SP
       zz_out(:, ne) = real(znv(:) + kpert(:, ne-1) * err_z, sp)
    end do

    ! --- T/S Perturbations ---
    if (ibarcl_rst /= 0) then

	write(*,*) '===================================================================='
	write(*,*) 'Baroclinic simulation, perturbing T and S, with depth smoothing...'
	write(*,*) '===================================================================='

        tt_out(:, :, 1) = real(tempv(:, :), sp)
        ss_out(:, :, 1) = real(saltv(:, :), sp)

        call make_2Dpert(kpert, dimk, nrens-1)
        do ne = 2, nrens
            do nl = 1, nnlv
                weight = real(hlv(nl), dp) / real(hlv(1), dp)
                tt_out(nl, :, ne) = real(tempv(nl, :) + kpert(:, ne-1) * err_t * weight, sp)
            end do
        end do

        call make_2Dpert(kpert, dimk, nrens-1)
        do ne = 2, nrens
            do nl = 1, nnlv
                weight = real(hlv(nl), dp) / real(hlv(1), dp)
                ss_out(nl, :, ne) = real(saltv(nl, :) + kpert(:, ne-1) * err_s * weight, sp)
            end do
        end do

    end if

    deallocate(kpert)
end subroutine make_perturbations

! ======================================================================
! Subroutine to write perturbed ensemble states to restart files
subroutine write_states(rstf, atime, nrens)
    use mod_perturbations
    use iso_fortran_env, only: dp=>real64, sp=>real32
    use mod_hydro      ! contain znv
    use mod_ts         ! contain saltv, tempv
    use mod_restart,   only: ibarcl_rst
    implicit none

    ! --- Arguments ---
    character(*), intent(in) :: rstf
    real(dp),     intent(in) :: atime
    integer,      intent(in) :: nrens

    ! --- Local Variables ---
    character(len=255) :: ens_rstf
    character(len=3)   :: nel
    integer            :: ne, ios, base_len

    base_len = index(rstf, '.rst', back=.true.)
    if (base_len == 0) base_len = len_trim(rstf)

    do ne = 1, nrens

       write(nel, '(I3.3)') ne - 1
       ens_rstf = rstf(1:base_len-1) // '_' // nel // '.rst'

       ! Update global fields (casting back to the type expected by modules)
       znv(:) = zz_out(:, ne)
       if (ibarcl_rst /= 0) then
          saltv(:, :) = ss_out(:, :, ne)
          tempv(:, :) = tt_out(:, :, ne)
       end if

       ! Open and write the restart record
       open(unit=22, file=trim(ens_rstf), form='unformatted', status='replace', iostat=ios)

       if (ios == 0) then
          write(*, '(A,A)') "Writing restart: ", trim(ens_rstf)
          call rst_write_record(atime, 22)
          close(22)
       else
          write(*, '(A,A)') "Error: Could not open file for writing: ", trim(ens_rstf)
       end if

    end do
end subroutine write_states

! ======================================================================
! PROGRAM perturbe_state
!
! Perturbes a state, saved in RST file in an nrens number of states.
! ======================================================================
program perturbe_state
    use iso_fortran_env, only: dp => real64, error_unit
    use iso8601 
    implicit none

    ! --- Variable Declarations ---
    character(len=255)      :: basinf, rstf
    character(len=20)       :: datestr
    character(len=32)       :: arg_buffer
    integer                 :: nnlv, nrens
    real(dp)                :: errz, errt, errs
    integer                 :: num_args


    ! internal variables
    integer                 :: date, time
    real(dp)                :: atime
    integer                 :: ierr

    ! --- 1. CLI Argument Validation ---
    num_args = command_argument_count()

    if (num_args /= 8) then
        write(error_unit, '(A)') "Error: Wrong number of arguments."
        write(error_unit, '(A)') "Usage: ./perturbe_state <basinf> <rstf> <date> <nlv> <nrens> <errz> <errt> <errs>"
        write(error_unit, '(A)') "Parameters:"
        write(error_unit, '(A)') "  basinf : basin file"
        write(error_unit, '(A)') "  rstf   : restart/state file"
        write(error_unit, '(A)') "  date   : Reference date (e.g., yyyy-mm-dd::HH:MM:SS)"
        write(error_unit, '(A)') "  nlv    : Number of vertical levels (integer)"
        write(error_unit, '(A)') "  nrens  : Number of ensemble members (must be odd, including control)"
        write(error_unit, '(A)') "  errz   : Perturbation amplitude for sea surface height (real)"
        write(error_unit, '(A)') "  errt   : Perturbation amplitude for temperature (real)"
        write(error_unit, '(A)') "  errs   : Perturbation amplitude for salinity (real)"
        stop 1
    end if

    ! --- 2. Parsing Arguments ---
    call get_command_argument(1, basinf)
    call get_command_argument(2, rstf)
    call get_command_argument(3, datestr)
    
    ! Parsing integer and real values from strings
    call get_command_argument(4, arg_buffer); read(arg_buffer, *) nnlv
    call get_command_argument(5, arg_buffer); read(arg_buffer, *) nrens
    call get_command_argument(6, arg_buffer); read(arg_buffer, *) errz
    call get_command_argument(7, arg_buffer); read(arg_buffer, *) errt
    call get_command_argument(8, arg_buffer); read(arg_buffer, *) errs

    ! --- 3. Consistency Checks ---
    ! Ensemble members must be odd to ensure a symmetric distribution around the control state
    if (mod(nrens, 2) == 0) then
        write(error_unit, '(A,I0,A)') "Error: nrens (", nrens, ") is even."
        error stop "The number of ensemble members must be odd to include the control state."
    end if

    ! --- 4. Initialization Summary ---
    write(*, '(/,A)') "--- perturbe_state setup ---"
    write(*, '(A,A)') "Basin file:    ", trim(basinf)
    write(*, '(A,A)') "Restart file:  ", trim(rstf)
    write(*, '(A,I0)') "Ens. members:  ", nrens
    write(*, '(A,A)') "Target date:   ", datestr
    write(*, '(A,3(F8.4,1X))') "Perturbations (z,t,s): ", errz, errt, errs
    write(*, '(A,/)') "----------------------------"

    ! Init shyfem basin and modules
    call init_shyfem(basinf, nnlv) ! Initialize SHYFEM grids

    ! Load a state from a RST file
    call string2date(datestr, date, time, ierr)
    if (ierr /= 0) error stop 'read_info: invalid date string'
    call dts_to_abs_time(date, time, atime)
    call read_rst(rstf, atime)

    ! perturbe state
    call make_perturbations(nrens, nnlv, errz, errt, errs)

    ! Write final states
    call write_states(rstf, atime, nrens)

end program perturbe_state
