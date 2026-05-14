! ======================================================================
!  MODULE: mod_enks_data  (used internally)
! ======================================================================
module mod_enks_data
   use iso_fortran_env, only: dp=>real64
   implicit none
   save

   ! Structure to store each analysis step's transformation
   type t_x5_record
      real(dp)              :: tt        ! Timestamp
      character(len=2)      :: tag       ! 'X3' or 'X5'
      integer               :: nrens, nrobs
      real(dp), allocatable :: mat(:,:)  ! Holds X5 (nrens,nrens) or X3 (nrobs,nrens)
      real(dp), allocatable :: S(:,:)    ! Only used if tag is 'X3' (nrobs,nrens)
   end type t_x5_record

   ! Global storage for all records found in X5_tot.uf
   type(t_x5_record), allocatable :: all_x5(:)
   integer  :: total_x5_records = 0
   real(dp) :: dt_x5 = -1.0_dp

contains

! ======================================================================
! LOAD ALL RECORDS FROM THE CUMULATIVE FILE
! ======================================================================
subroutine load_all_x5(nrens_expected)
   integer, intent(in) :: nrens_expected
   integer :: u, ios, count
   real(dp) :: tmp_tt
   character(len=2)  :: tmp_tag
   integer :: tmp_nrens, tmp_nrobs

   ! 1. First pass: count records to allocate the array
   open(newunit=u, file='X5_tot.uf', form='unformatted', status='old', action='read', iostat=ios)
   if (ios /= 0) error stop "load_all_x5: Cannot open X5_tot.uf. Check if analysis generated it."

   count = 0
   do
      read(u, iostat=ios) tmp_tt, tmp_tag
      if (ios /= 0) exit ! End of file reached

      if (tmp_tag == 'X3') then
         read(u) tmp_nrens, tmp_nrobs
         read(u) ! Skip X3 data
         read(u) ! Skip S data
      else
         read(u) tmp_nrens
         read(u) ! Skip X5 data
      end if
      count = count + 1
   end do

   total_x5_records = count
   if (total_x5_records == 0) error stop "load_all_x5: No records found in X5_tot.uf"

   if (allocated(all_x5)) deallocate(all_x5)
   allocate(all_x5(total_x5_records))
   rewind(u)

   ! 2. Second pass: load data into memory
   write(*,*) "load_all_x5: Loading", total_x5_records, " records..."

   do count = 1, total_x5_records
      read(u) all_x5(count)%tt, all_x5(count)%tag

      if (all_x5(count)%tag == 'X3') then
         read(u) all_x5(count)%nrens, all_x5(count)%nrobs
         allocate(all_x5(count)%mat(all_x5(count)%nrobs, all_x5(count)%nrens)) ! X3 matrix
         allocate(all_x5(count)%S(all_x5(count)%nrobs, all_x5(count)%nrens))   ! S matrix
         read(u) all_x5(count)%mat
         read(u) all_x5(count)%S
      else
         read(u) all_x5(count)%nrens
         all_x5(count)%nrobs = 0
         allocate(all_x5(count)%mat(all_x5(count)%nrens, all_x5(count)%nrens)) ! X5 matrix
         read(u) all_x5(count)%mat
      end if

      ! Consistency check
      if (all_x5(count)%nrens /= nrens_expected) then
         write(*,*) "Warning: record", count, " has nrens =", all_x5(count)%nrens, &
                    " but expected", nrens_expected
      end if
   end do

   close(u)
   write(*,*) "load_all_x5: Done."

   ! Sort records by time to be safe
   call check_sort_x5()
   ! Set time step
   dt_x5 = estimate_dt_x5()
end subroutine load_all_x5

! ======================================================================
! SORT RECORDS BY TIMESTAMP (Selection Sort)
! ======================================================================
subroutine check_sort_x5()
   implicit none
   integer :: i, j, kmin
   type(t_x5_record) :: tmp

   if (total_x5_records <= 1) return

   do i=1, total_x5_records-1
      kmin = i
      do j=i+1, total_x5_records
         if (all_x5(j)%tt < all_x5(kmin)%tt) kmin=j
      end do
      if (kmin /= i) then
         tmp = all_x5(i)
         all_x5(i) = all_x5(kmin)
         all_x5(kmin) = tmp
      end if
   end do
end subroutine check_sort_x5

! ======================================================================
! ESTIMATE TIME STEP BETWEEN ANALYSES
! ======================================================================
real(dp) function estimate_dt_x5()
   implicit none
   integer :: i
   estimate_dt_x5 = 0.0_dp
   if (total_x5_records < 2) return

   do i = 2, total_x5_records
      if (abs(all_x5(i)%tt - all_x5(i-1)%tt) > 1.0e-8_dp) then
         estimate_dt_x5 = all_x5(i)%tt - all_x5(i-1)%tt
         return
      end if
   end do
end function estimate_dt_x5

! ======================================================================
! FLOATING POINT TIME COMPARISON
! ======================================================================
logical function equal_time(t1,t2)
   real(dp), intent(in) :: t1,t2
   equal_time = abs(t1-t2) < 1.0e-6_dp
end function equal_time

! ======================================================================
! INITIALIZE STATE VECTORS BASED ON SHYFEM CONFIG
! ======================================================================
subroutine allocate_states(A, Amean, Astd, n, nrens)
   use mod_hydro
   use mod_hydro_vel
   use mod_ts
   use mod_restart,   only : ibarcl_rst
   implicit none

   real(dp), allocatable, intent(out) :: A(:,:), Amean(:), Astd(:)
   integer, intent(out)               :: n
   integer, intent(in)                :: nrens

   ! Calculate total state dimension (Elevation + Velocity U/V)
   n = size(utlnv) + size(vtlnv) + size(znv)

   ! Add Temperature and Salinity if baroclinic
   if (ibarcl_rst /= 0) n = n + size(tempv) + size(saltv)

   if (allocated(A))     deallocate(A)
   if (allocated(Amean)) deallocate(Amean)
   if (allocated(Astd))  deallocate(Astd)

   allocate(A(n, nrens))
   allocate(Amean(n), Astd(n))

   write(*,*) "allocate_states: Total state dimension n =", n
end subroutine allocate_states

end module mod_enks_data

! ======================================================================
! SUBROUTINES
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
subroutine num2str(num, str)
   implicit none
   integer, intent(in)        :: num
   character(len=5), intent(out) :: str

   if (num < 0 .or. num > 99999) then 
      error stop "num2str: argument out of range"
   end if   
            
   write(str, '(I5.5)') num
end subroutine num2str

! ======================================================================
subroutine read_rst(rstname, atimea)
   use iso_fortran_env, only : dp=>real64
   use mod_restart
   use levels, only : nlvdi, nlv, hlv, ilhv, ilhkv
   use shympi
   use mod_enks_data, only : equal_time
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
      if (equal_time(atimef, atimea)) exit
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
! Write a restart record at time atimea.
subroutine rst_write_rec(atimea, iunit)
  use iso_fortran_env, only : dp => real64
  use basin, only : nel, nkn
  use levels, only : nlv
  use mod_hydro
  use mod_hydro_vel
  use mod_layer_thickness
  implicit none
  integer, intent(in) :: iunit
  real(dp),        intent(in) :: atimea
  integer :: ios

  call mod_layer_thickness_init(nkn,nel,nlv)

  call setzev_enkf

  zov = znv
  utlov = utlnv
  vtlov = vtlnv
  zeov = zenv

  call rst_write_record(atimea, iunit)

end subroutine rst_write_rec

!----------------------------------------------------------------------
! Duplicate of setzev with some modifications
!----------------------------------------------------------------------
subroutine setzev_enkf

  use mod_geom_dynamic
  use mod_hydro
  use basin

  implicit none

! local
  integer ie,ii

  iwetv = 0

  do ie=1,nel
      do ii=1,3
        zenv(ii,ie)=znv(nen3v(ii,ie))
        if ((zenv(ii,ie) + hm3v(ii,ie)) < 0.05) then
	  iwegv(ie) = iwegv(ie) + 1
	  zenv(ii,ie) = 0.02 - hm3v(ii,ie)
          znv(nen3v(ii,ie)) = zenv(ii,ie)
	end if
      end do
  end do

end subroutine setzev_enkf


! ======================================================================
subroutine push_matrix(sdim, nrens, nre, Amat)
   use iso_fortran_env, only: dp=>real64
   use mod_hydro
   use mod_hydro_vel
   use mod_ts
   use mod_restart,   only: ibarcl_rst
   implicit none

   integer, intent(in) :: sdim, nrens, nre
   real(dp), intent(inout) :: Amat(sdim,nrens)
   integer :: d_uv, d_z, d_ts

   d_uv = size(utlnv)
   d_z  = size(znv)

   ! Correct mapping for SHYFEM state vector
   Amat(1:d_uv, nre) = reshape(real(utlnv,dp), [d_uv])
   Amat(d_uv+1:2*d_uv, nre) = reshape(real(vtlnv,dp), [d_uv])
   Amat(2*d_uv+1:2*d_uv+d_z, nre) = real(znv,dp)

   if (ibarcl_rst /= 0) then
      d_ts = size(tempv)
      Amat(2*d_uv+d_z+1 : 2*d_uv+d_z+d_ts, nre) = reshape(real(tempv,dp), [d_ts])
      Amat(2*d_uv+d_z+d_ts+1 : 2*d_uv+d_z+2*d_ts, nre) = reshape(real(saltv,dp), [d_ts])
   end if
end subroutine push_matrix

! ======================================================================
subroutine pull_matrix(sdim, nrens, nre, Amat)
   use iso_fortran_env, only: dp=>real64
   use mod_hydro
   use mod_hydro_vel
   use mod_ts
   use mod_restart,   only: ibarcl_rst
   implicit none

   integer, intent(in) :: sdim, nrens, nre
   real(dp), intent(in) :: Amat(sdim,nrens)
   integer :: d_uv, d_z, d_ts, nnkn, nnel, nnlv

   d_z   = size(znv)
   d_uv  = size(utlnv)
   nnkn  = size(znv)
   nnel  = size(utlnv, 2)
   nnlv  = size(utlnv, 1)

   utlnv = reshape(real(Amat(1:d_uv, nre), dp), [nnlv, nnel])
   vtlnv = reshape(real(Amat(d_uv+1:2*d_uv, nre), dp), [nnlv, nnel])
   znv   = real(Amat(2*d_uv+1:2*d_uv+d_z, nre), dp)

   if (ibarcl_rst /= 0) then
      d_ts = size(tempv)
      tempv = reshape(real(Amat(2*d_uv+d_z+1:2*d_uv+d_z+d_ts, nre), dp), [nnlv, nnkn])
      saltv = reshape(real(Amat(2*d_uv+d_z+d_ts+1:2*d_uv+d_z+2*d_ts, nre), dp), [nnlv, nnkn])
   end if
end subroutine pull_matrix

! ======================================================================
! ======================================================================
! Lagged smoother analysis routine
! ======================================================================
subroutine make_analysis(atime, sdim, nrens, Amat, nlag)
   use iso_fortran_env, only : dp=>real64
   use mod_enks_data
   implicit none

   real(dp), intent(in)    :: atime
   integer, intent(in)     :: sdim, nrens, nlag
   real(dp), intent(inout) :: Amat(sdim,nrens)

   real(dp) :: t_end
   integer :: i, count, j
   real(dp), allocatable :: Xacc(:,:), Xtmp(:,:), Amat_tmp(:,:), X5_equiv(:,:)

   ! 1. Determine the end of the smoothing window
   if (nlag == -1) then
      t_end = huge(1.0_dp)
   else
      t_end = atime + real(nlag, dp) * dt_x5 + 1.0e-6_dp
   end if

   allocate(Xacc(nrens,nrens), Xtmp(nrens,nrens))

   ! Initialise Xacc as Identity matrix
   Xacc = 0.0_dp
   do i=1,nrens
      Xacc(i,i) = 1.0_dp
   end do

   count = 0
   ! 2. ACCUMULATE all analysis matrices within the lag window
   do i=1, total_x5_records
      if (all_x5(i)%tt >= atime .and. all_x5(i)%tt <= t_end) then
         
         ! If the record is X3, we must convert it to X5-equivalent first
         if (all_x5(i)%tag == 'X3') then
            allocate(X5_equiv(nrens, nrens))
            X5_equiv = 0.0_dp
            do j = 1, nrens
               X5_equiv(j, j) = 1.0_dp
            end do
            call dgemm('T', 'N', nrens, nrens, all_x5(i)%nrobs, 1.0_dp, &
                       all_x5(i)%S, all_x5(i)%nrobs, all_x5(i)%mat, all_x5(i)%nrobs, &
                       1.0_dp, X5_equiv, nrens)
            
            call dgemm("N", "N", nrens, nrens, nrens, 1.0_dp, Xacc, nrens, &
                       X5_equiv, nrens, 0.0_dp, Xtmp, nrens)
            deallocate(X5_equiv)
         else
            ! Standard X5 multiplication
            call dgemm("N", "N", nrens, nrens, nrens, 1.0_dp, Xacc, nrens, &
                       all_x5(i)%mat, nrens, 0.0_dp, Xtmp, nrens)
         end if
         
         Xacc = Xtmp
         count = count + 1
      end if
   end do

   ! 3. APPLY the accumulated transformation ONCE to the state ensemble
   if (count > 0) then
      allocate(Amat_tmp(sdim, nrens))
      call dgemm("N", "N", sdim, nrens, nrens, 1.0_dp, Amat, sdim, &
                 Xacc, nrens, 0.0_dp, Amat_tmp, sdim)
      Amat = Amat_tmp
      deallocate(Amat_tmp)
   end if

   deallocate(Xacc, Xtmp)
end subroutine make_analysis

! ======================================================================
! Stats
subroutine make_mn_std(n, nens, Amat, mean_vec, std_vec)
    use iso_fortran_env, only: dp => real64
    implicit none

    ! Input/Output
    integer, intent(in) :: n, nens
    real(dp), intent(in) :: Amat(n, nens)
    real(dp), intent(out) :: mean_vec(n)
    real(dp), intent(out) :: std_vec(n)

    integer :: i

    do i = 1, n
        mean_vec(i) = sum(Amat(i, :)) / real(nens, dp)
    end do

    do i = 1, n
        std_vec(i) = sqrt(sum((Amat(i, :) - mean_vec(i))**2) / real(nens, dp))
    end do
end subroutine make_mn_std

! ======================================================================
! PROGRAM enKF2enKS
!
! Compute the ensemble mean and std of the KS from the KF and the X5 matrices.
! All the members are necessary.
! It is easy to implement the output for each member.
! ======================================================================
program enKF2enKS
    use iso_fortran_env, only: dp=>real64
    use mod_enks_data
    implicit none

    character(len=80) :: basinf
    character(len=3)  :: lnnlv
    character(len=6)  :: lnrens, lnlag
    integer           :: nnlv, nrens, nlag
    integer           :: rrec, nre, sdim
    character(len=5)  :: nrel
    real(dp)          :: atime
    real(dp), allocatable :: Astate(:,:), AmeanKS(:), AstdKS(:)

    ! 1. CLI Arguments Parsing
    call get_command_argument(1, basinf)
    call get_command_argument(2, lnnlv)
    call get_command_argument(3, lnrens)
    call get_command_argument(4, lnlag)

    if (lnlag == '') stop "Usage: enKF2enKS [basinf] [nlv] [nrens] [nlag]"

    read(lnnlv, *) nnlv
    read(lnrens, *) nrens
    read(lnlag, *) nlag

    write(*,*) "EnKS setup: Lag =", nlag, " Members =", nrens

    ! 2. Setup Data
    call load_all_x5(nrens)      ! Load X5_tot.uf into memory
    call init_shyfem(basinf, nnlv) ! Initialize SHYFEM grids

    ! 3. Opening Output Restarts
    open(18, file="analKS_mean.rst", status="replace", form="unformatted")
    open(19, file="analKS_std.rst",  status="replace", form="unformatted")

    ! ======================================================================
    ! MAIN TIME LOOP
    ! ======================================================================
    do rrec = 1, total_x5_records
       atime = all_x5(rrec)%tt

       ! -----------------------------------
       ! Read all ensemble members for current time
       ! -----------------------------------
       do nre = 1, nrens
          call num2str(nre-1, nrel)

          ! Read the specific restart for this member and time
          call read_rst("analKF_en"//nrel//".rst", atime)

          ! Allocate matrices on the first successful read
          if (.not. allocated(Astate)) then
              call allocate_states(Astate, AmeanKS, AstdKS, sdim, nrens)
          end if

          ! Map SHYFEM variables into the Astate matrix column
          call push_matrix(sdim, nrens, nre, Astate)
       end do

       ! -----------------------------------
       ! Execute Smoother Analysis
       ! Apply future weights: A = A * (X5_t * X5_t+1 * ... * X5_t+lag)
       ! -----------------------------------
       call make_analysis(atime, sdim, nrens, Astate, nlag)

       ! -----------------------------------
       ! Statistics & Output Generation
       ! -----------------------------------
       call make_mn_std(sdim, nrens, Astate, AmeanKS, AstdKS)

       ! Write Mean Restart
       call pull_matrix(sdim, nrens, 1, AmeanKS) ! Overwrites SHYFEM globals
       call rst_write_rec(atime, 18)

       ! Write Std Dev Restart
       call pull_matrix(sdim, nrens, 1, AstdKS)  ! Overwrites SHYFEM globals
       call rst_write_rec(atime, 19)

       write(*,*) "Processed Record:", rrec, " Time:", atime
    end do

    close(18)
    close(19)
    write(*,*) "Smoothing completed successfully."

end program enKF2enKS
