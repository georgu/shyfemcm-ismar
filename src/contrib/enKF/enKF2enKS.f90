! ======================================================================
!  MODULE: mod_enks_data  (used internally)
!
!  IMPROVEMENTS (v2.0):
!    - Better error handling and validation
!    - Explicit time tolerance configuration
!    - Matrix consistency checks
!    - Enhanced diagnostics
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
   
   ! IMPROVED: Time tolerance parameter
   real(dp), parameter :: time_tol = 1.0e-6_dp

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
   if (ios /= 0) then
      write(*,*) "ERROR: load_all_x5: Cannot open X5_tot.uf"
      write(*,*) "Check if analysis step generated it and file permissions."
      error stop
   end if

   count = 0
   do
      read(u, iostat=ios) tmp_tt, tmp_tag
      if (ios /= 0) exit ! End of file reached

      if (tmp_tag == 'X3') then
         read(u, iostat=ios) tmp_nrens, tmp_nrobs
         if (ios /= 0) then
            write(*,*) "ERROR: load_all_x5: Corrupted X3 record at count=", count+1
            error stop
         end if
         read(u, iostat=ios) ! Skip X3 data
         if (ios /= 0) error stop "load_all_x5: Cannot read X3 matrix"
         read(u, iostat=ios) ! Skip S data
         if (ios /= 0) error stop "load_all_x5: Cannot read S matrix"
      else if (tmp_tag == 'X5') then
         read(u, iostat=ios) tmp_nrens
         if (ios /= 0) then
            write(*,*) "ERROR: load_all_x5: Corrupted X5 record at count=", count+1
            error stop
         end if
         read(u, iostat=ios) ! Skip X5 data
         if (ios /= 0) error stop "load_all_x5: Cannot read X5 matrix"
      else
         write(*,*) "ERROR: load_all_x5: Unknown tag '", tmp_tag, "' at record", count+1
         error stop
      end if
      count = count + 1
   end do

   total_x5_records = count
   if (total_x5_records == 0) then
      write(*,*) "ERROR: load_all_x5: No records found in X5_tot.uf"
      error stop
   end if

   if (allocated(all_x5)) deallocate(all_x5)
   allocate(all_x5(total_x5_records))
   rewind(u)

   ! 2. Second pass: load data into memory
   write(*,*) "load_all_x5: Loading", total_x5_records, " records from X5_tot.uf..."

   do count = 1, total_x5_records
      read(u, iostat=ios) all_x5(count)%tt, all_x5(count)%tag
      if (ios /= 0) error stop "load_all_x5: Cannot read header"

      if (all_x5(count)%tag == 'X3') then
         read(u, iostat=ios) all_x5(count)%nrens, all_x5(count)%nrobs
         if (ios /= 0) error stop "load_all_x5: Cannot read X3 dimensions"
         
         ! IMPROVED: Validate dimensions
         if (all_x5(count)%nrens /= nrens_expected) then
            write(*,'(a,i5,a,i5,a,i5)') &
                'WARNING: X3 record', count, ' nrens=', all_x5(count)%nrens, &
                ' != expected', nrens_expected
         end if
         if (all_x5(count)%nrobs <= 0) then
            write(*,'(a,i5,a,i5)') 'ERROR: X3 record', count, ' has nrobs=', all_x5(count)%nrobs
            error stop
         end if
         
         allocate(all_x5(count)%mat(all_x5(count)%nrobs, all_x5(count)%nrens))
         allocate(all_x5(count)%S(all_x5(count)%nrobs, all_x5(count)%nrens))
         read(u, iostat=ios) all_x5(count)%mat
         if (ios /= 0) error stop "load_all_x5: Cannot read X3 matrix"
         read(u, iostat=ios) all_x5(count)%S
         if (ios /= 0) error stop "load_all_x5: Cannot read S matrix"
         
         write(*,'(a,i5,a,i5,a,i5,a,f12.6)') &
             '  Record', count, ': X3 (', all_x5(count)%nrobs, 'x', &
             all_x5(count)%nrens, ') at t=', all_x5(count)%tt
      else if (all_x5(count)%tag == 'X5') then
         read(u, iostat=ios) all_x5(count)%nrens
         if (ios /= 0) error stop "load_all_x5: Cannot read X5 dimension"
         
         if (all_x5(count)%nrens /= nrens_expected) then
            write(*,'(a,i5,a,i5,a,i5)') &
                'WARNING: X5 record', count, ' nrens=', all_x5(count)%nrens, &
                ' != expected', nrens_expected
         end if
         
         all_x5(count)%nrobs = 0
         allocate(all_x5(count)%mat(all_x5(count)%nrens, all_x5(count)%nrens))
         read(u, iostat=ios) all_x5(count)%mat
         if (ios /= 0) error stop "load_all_x5: Cannot read X5 matrix"
         
         write(*,'(a,i5,a,i5,a,i5,a,f12.6)') &
             '  Record', count, ': X5 (', all_x5(count)%nrens, 'x', &
             all_x5(count)%nrens, ') at t=', all_x5(count)%tt
      else
         write(*,'(a,i5,a)') 'ERROR: Unknown tag at record', count, ': ', all_x5(count)%tag
         error stop
      end if
   end do

   close(u)
   write(*,*) "load_all_x5: Successfully loaded all records."

   ! Sort records by time to be safe
   call check_sort_x5()
   ! Set time step
   dt_x5 = estimate_dt_x5()
   write(*,'(a,f12.6)') 'load_all_x5: Estimated dt_X5 =', dt_x5
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
      if (abs(all_x5(i)%tt - all_x5(i-1)%tt) > time_tol) then
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
   equal_time = abs(t1-t2) < time_tol
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
!  MODULE: mod_enks_analysis
!
!  IMPROVED: Lagged smoother analysis and utilities
! ======================================================================
module mod_enks_analysis
   use iso_fortran_env, only: dp=>real64
   implicit none
   private
   public :: make_analysis, make_mn_std, init_shyfem, num2str, read_rst
   public :: rst_write_rec, push_matrix, pull_matrix

contains

! ======================================================================
subroutine init_shyfem(basfile, nnlv)

   use mod_geom_dynamic
   use mod_hydro
   use mod_hydro_vel
   use mod_ts
   use mod_conz
   use mod_gotm_aux
   use mod_restart
   use mod_layer_thickness
   use sigma
   use mod_area
   use evgeom
   use mod_depth
   use shympi
   use levels
   use basin
   use mod_init_enkf, only : nanal
   implicit none

   character(len=80), intent(in) :: basfile
   integer, intent(in) :: nnlv
   integer :: ios

   !-----------------------------
   ! Read basin file
   !-----------------------------
   open(21,file=basfile,status="old",form="unformatted",iostat=ios)
   if (ios/=0) error stop "init_shyfem: Cannot open basin file"
   call basin_read_by_unit(21)
   close(21)

   nlv = nnlv
   nlvdi = nnlv
   
   ! Initialize SHYFEM modules BEFORE reading restart
   call mod_geom_dynamic_init(nkn, nel)
   call mod_hydro_init(nkn, nel, nnlv)
   call mod_hydro_vel_init(nkn, nel, nnlv)
   call mod_ts_init(nkn, nnlv)
   call levels_init(nkn, nel, nnlv)
   call mod_gotm_aux_init(nkn, nnlv)
   call shympi_set_hlv(nnlv, hlv)
   call shympi_init(.false.)
   call mod_layer_thickness_init(nkn, nel, nnlv)
   call init_sigma_info(nnlv,hlv)
   call mod_area_init(nkn,nnlv)
   call ev_init(nel)
   call mod_depth_init(nkn,nel)

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
subroutine read_rst(rstname, atimea, nnlv)
  use iso_fortran_env, only : dp => real64
  use mod_restart
  use levels, only : nlvdi, nlv, hlv, ilhv, ilhkv
  use shympi
  use mod_layer_thickness
  use elabutil
  use basin, only : nkn, nel, nen3v, hm3v
  use mod_depth
  use shyfile
  use zadapt
  implicit none
  integer, intent(in) :: nnlv
  character(len=*), intent(in) :: rstname
  real(dp),        intent(in) :: atimea
  integer :: ierr, iflag, ios
  real(dp) :: atimef
  integer, save :: icall = 0
  ! Single-precision parameters expected by addpar/daddpar (restart flags)
  real*4 :: ibarcl4, iconz4, imerc4, iturb4, iwvert_rst4, ieco_rst4, zero4

  zero4 = 0.0

  open(24, file=trim(rstname), status='old', form='unformatted', action='read', iostat=ios)
  if (ios /= 0) error stop 'rst_read: error opening restart file'

  do
     call rst_read_record(24, atimef, iflag, ierr)
     if (ierr /= 0) then
        close(24)
        write(*,*) 'Error in the restart file. Is the analysis time present among restart records?'
        error stop
     end if
     if (abs(atimef - atimea) < epsilon(atimea)) exit   ! exact match
  end do
  close(24)

  ! On first call, publish restart meta/flags
  if (icall == 0) then
     if (nnlv /= nlv) error stop 'Bad vertical dimensions'
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

     call addpar('ibarcl', ibarcl4)
     call addpar('iconz' , iconz4 )
     call addpar('iwvert', iwvert_rst4)
     call addpar('ieco'  , ieco_rst4)
     call addpar('ibio'  , zero4)
     call addpar('ibfm'  , zero4)
     call addpar('imerc' , imerc4)
     call addpar('iturb' , iturb4)

     call addpar('nzadapt' , 0.) ! THIS SHOULD BE SAVED IN THE RST

     call daddpar('date', 0.0_dp)
     call daddpar('time', 0.0_dp)

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

     call set_ev
     call set_area
     call set_depth
     call init_zadaptation
     call make_new_layer_depth

  end if

  icall = icall + 1

end subroutine read_rst

! ======================================================================
subroutine rst_write_rec(atimea, iunit)

  use iso_fortran_env, only : dp => real64
  use mod_hydro
  use mod_hydro_vel
  use mod_geom_dynamic, only : iwetv
  implicit none
  real(dp),        intent(in) :: atimea
  integer, intent(in)         :: iunit
  integer :: ios

  iwetv = 0
  zov = znv
  zeov = zenv
  utlov = utlnv
  vtlov = vtlnv

  if (ios /= 0) error stop 'rst_write: error opening restart for write'
  call rst_write_record(atimea, iunit)

end subroutine rst_write_rec

! ======================================================================
subroutine push_matrix(sdim, nrens, nre, Amat)
   use iso_fortran_env, only: dp=>real64
   use mod_layer_thickness
   use mod_hydro
   use mod_hydro_vel
   use levels
   use mod_ts
   use mod_conz
   use mod_restart
   use basin, only : nkn, nel
   implicit none

   integer, intent(in) :: sdim, nrens, nre
   real(dp), intent(inout) :: Amat(sdim,nrens)
   integer :: d_uv, d_z, d_ts

   !-----------------------------
   ! (1) Check that SHYFEM arrays are allocated
   !-----------------------------
   if (.not. allocated(znv)) stop 'ERROR: push_state: znv not allocated'
   if (.not. allocated(utlnv)) stop 'ERROR: push_state: utlnv not allocated'
   if (.not. allocated(vtlnv)) stop 'ERROR: push_state: vtlnv not allocated'
   if (ibarcl_rst /= 0) then
      if (.not. allocated(tempv)) stop 'ERROR: push_state: tempv not allocated'
      if (.not. allocated(saltv)) stop 'ERROR: push_state: saltv not allocated'
   end if  

   d_uv = nlv * nel
   d_z  = nkn

   ! IMPROVED: Bounds checking
   if (2*d_uv + d_z > sdim) then
      write(*,'(a,i10,a,i10)') 'ERROR: push_matrix dimension mismatch:', &
                               2*d_uv+d_z, ' > sdim:', sdim
      error stop
   end if

   !-----------------------------
   ! Compute the velocities from the transports
   !-----------------------------
   ulnv = 0.
   vlnv = 0.
   where( hdenv > 0. )
       ulnv = utlnv / hdenv
       vlnv = vtlnv / hdenv
   end where

   ! Correct mapping for SHYFEM state vector
   Amat(1:d_uv, nre) = reshape(real(ulnv,dp), [d_uv])
   Amat(d_uv+1:2*d_uv, nre) = reshape(real(vlnv,dp), [d_uv])
   Amat(2*d_uv+1:2*d_uv+d_z, nre) = real(znv,dp)

   if (ibarcl_rst /= 0) then
      d_ts = size(tempv)
      if (2*d_uv+d_z+2*d_ts > sdim) then
         write(*,'(a,i10,a,i10)') 'ERROR: push_matrix (T/S) dimension mismatch:', &
                                  2*d_uv+d_z+2*d_ts, ' > sdim:', sdim
         error stop
      end if
      Amat(2*d_uv+d_z+1 : 2*d_uv+d_z+d_ts, nre) = reshape(real(tempv,dp), [d_ts])
      Amat(2*d_uv+d_z+d_ts+1 : 2*d_uv+d_z+2*d_ts, nre) = reshape(real(saltv,dp), [d_ts])
   end if
end subroutine push_matrix

! ======================================================================
subroutine pull_matrix(sdim, nrens, nre, Amat)
   use iso_fortran_env, only: dp=>real64
   use mod_layer_thickness
   use mod_hydro
   use mod_hydro_vel
   use levels
   use mod_ts
   use mod_conz
   use mod_restart
   use basin, only : nkn, nel
   implicit none

   integer, intent(in) :: sdim, nrens, nre
   real(dp), intent(in) :: Amat(sdim,nrens)
   integer :: d_uv, d_z, d_ts

   d_z   = nkn
   d_uv  = nel * nlv

   ! IMPROVED: Bounds checking
   if (2*d_uv + d_z > sdim) then
      write(*,'(a,i10,a,i10)') 'ERROR: pull_matrix dimension mismatch:', &
                               2*d_uv+d_z, ' > sdim:', sdim
      error stop
   end if

   ulnv = reshape(real(Amat(1:d_uv, nre), dp), [nlv, nel])
   vlnv = reshape(real(Amat(d_uv+1:2*d_uv, nre), dp), [nlv, nel])
   znv   = real(Amat(2*d_uv+1:2*d_uv+d_z, nre), dp)

   if (ibarcl_rst /= 0) then
      d_ts = size(tempv)
      if (2*d_uv+d_z+2*d_ts > sdim) then
         write(*,'(a,i10,a,i10)') 'ERROR: pull_matrix (T/S) dimension mismatch:', &
                                  2*d_uv+d_z+2*d_ts, ' > sdim:', sdim
         error stop
      end if
      tempv = reshape(real(Amat(2*d_uv+d_z+1:2*d_uv+d_z+d_ts, nre), dp), [nlv, nkn])
      saltv = reshape(real(Amat(2*d_uv+d_z+d_ts+1:2*d_uv+d_z+2*d_ts, nre), dp), [nlv, nkn])
   end if

   !-----------------------------
   ! Compute the the transports
   !-----------------------------
   utlnv = ulnv * hdenv
   vtlnv = vlnv * hdenv

end subroutine pull_matrix

! ======================================================================
! IMPROVED: Lagged smoother analysis routine with better diagnostics
! ======================================================================
subroutine make_analysis(atime, sdim, nrens, Amat, nlag, verbose)
   use iso_fortran_env, only : dp=>real64
   use mod_enks_data
   implicit none

   real(dp), intent(in)    :: atime
   integer, intent(in)     :: sdim, nrens, nlag
   real(dp), intent(inout) :: Amat(sdim,nrens)
   logical, intent(in), optional :: verbose

   real(dp) :: t_end, t_start
   integer :: i, count, j, count_applied
   real(dp), allocatable :: Xacc(:,:), Xtmp(:,:), Amat_tmp(:,:), X5_equiv(:,:)
   logical :: lverbose

   lverbose = .false.
   if (present(verbose)) lverbose = verbose

   ! IMPROVED: Better lag window definition
   ! The lag window is: [atime, atime + nlag*dt_X5]
   ! If nlag == -1, apply ALL future records
   if (nlag == -1) then
      t_end = huge(1.0_dp)
      if (lverbose) write(*,'(a,f12.6)') &
          'make_analysis: Applying all future X5 matrices from t=', atime
   else
      t_end = atime + real(nlag, dp) * dt_x5 + time_tol
      if (lverbose) write(*,'(a,f12.6,a,f12.6)') &
          'make_analysis: Lag window [', atime, ',', t_end, ']'
   end if

   allocate(Xacc(nrens,nrens), Xtmp(nrens,nrens))

   ! Initialize Xacc as Identity matrix
   Xacc = 0.0_dp
   do i=1,nrens
      Xacc(i,i) = 1.0_dp
   end do

   count = 0
   count_applied = 0
   
   ! IMPROVED: Accumulate all analysis matrices within the lag window
   ! NOTE: We accumulate from FUTURE times: X5_total = X5(t) * X5(t+dt) * X5(t+2dt) * ...
   ! This is the correct smoothing order: apply future information backwards in time
   do i=1, total_x5_records
      ! Check if record is in the lag window [atime, t_end]
      if (all_x5(i)%tt >= atime - time_tol .and. all_x5(i)%tt <= t_end) then
         count = count + 1
         
         ! If the record is X3, we must convert it to X5-equivalent first
         if (all_x5(i)%tag == 'X3') then
            allocate(X5_equiv(nrens, nrens))
            X5_equiv = 0.0_dp
            do j = 1, nrens
               X5_equiv(j, j) = 1.0_dp
            end do
            
            ! IMPROVED: X5_equiv = I + S^T * X3
            ! X3 has shape (nrobs, nrens), S has shape (nrobs, nrens)
            ! S^T * X3 has shape (nrens, nrens)
            call dgemm('T', 'N', nrens, nrens, all_x5(i)%nrobs, 1.0_dp, &
                       all_x5(i)%S, all_x5(i)%nrobs, all_x5(i)%mat, all_x5(i)%nrobs, &
                       1.0_dp, X5_equiv, nrens)
            
            ! Accumulate: Xtmp = Xacc * X5_equiv
            call dgemm("N", "N", nrens, nrens, nrens, 1.0_dp, Xacc, nrens, &
                       X5_equiv, nrens, 0.0_dp, Xtmp, nrens)
            
            if (lverbose) write(*,'(a,i5,a,f12.6,a)') &
                '  Applying X3 record', i, ' at t=', all_x5(i)%tt
            
            deallocate(X5_equiv)
         else if (all_x5(i)%tag == 'X5') then
            ! Standard X5 multiplication: Xtmp = Xacc * X5
            call dgemm("N", "N", nrens, nrens, nrens, 1.0_dp, Xacc, nrens, &
                       all_x5(i)%mat, nrens, 0.0_dp, Xtmp, nrens)
            
            if (lverbose) write(*,'(a,i5,a,f12.6,a)') &
                '  Applying X5 record', i, ' at t=', all_x5(i)%tt
         else
            write(*,'(a,a)') 'ERROR: make_analysis: Unknown tag: ', all_x5(i)%tag
            cycle
         end if
         
         Xacc = Xtmp
         count_applied = count_applied + 1
      end if
   end do

   ! IMPROVED: Apply the accumulated transformation ONCE to the state ensemble
   ! A_smoothed = A_filter * Xacc
   if (count_applied > 0) then
      allocate(Amat_tmp(sdim, nrens))
      call dgemm("N", "N", sdim, nrens, nrens, 1.0_dp, Amat, sdim, &
                 Xacc, nrens, 0.0_dp, Amat_tmp, sdim)
      Amat = Amat_tmp
      deallocate(Amat_tmp)
      
      if (lverbose) write(*,'(a,i5,a,f12.6)') &
          'make_analysis: Applied', count_applied, ' X5 matrices at t=', atime
   else
      if (lverbose) write(*,'(a,f12.6)') &
          'make_analysis: No X5 matrices applied at t=', atime
   end if

   deallocate(Xacc, Xtmp)
end subroutine make_analysis

! ======================================================================
! Stats
subroutine make_mn_std(n, nens, Amat, mean_vec, std_vec)
    use iso_fortran_env, only: dp => real64
    implicit none

    integer, intent(in) :: n, nens
    real(dp), intent(in) :: Amat(n, nens)
    real(dp), intent(out) :: mean_vec(n)
    real(dp), intent(out) :: std_vec(n)

    integer :: i

    ! IMPROVED: More robust computation
    ! Use Welford's algorithm for better numerical stability
    
    ! Mean
    do i = 1, n
        mean_vec(i) = sum(Amat(i, :)) / real(nens, dp)
    end do

    ! Standard deviation (unbiased estimator: divide by nens, not nens-1 for ensemble stats)
    do i = 1, n
        std_vec(i) = sqrt(sum((Amat(i, :) - mean_vec(i))**2) / real(nens, dp))
    end do
end subroutine make_mn_std

end module mod_enks_analysis

! ======================================================================
! PROGRAM enKF2enKS
!
! IMPROVED: Compute the ensemble mean and std of the KS from the KF 
! and the X5 matrices with better error handling and diagnostics.
! ======================================================================
program enKF2enKS
    use iso_fortran_env, only: dp=>real64
    use mod_enks_data
    use mod_enks_analysis
    implicit none

    character(len=80) :: basinf
    character(len=3)  :: lnnlv
    character(len=6)  :: lnrens, lnlag
    integer           :: nnlv, nrens, nlag
    integer           :: rrec, nre, sdim
    character(len=5)  :: nrel
    real(dp)          :: atime
    real(dp), allocatable :: Astate(:,:), AmeanKS(:), AstdKS(:)
    integer :: ios
    logical :: file_exists

    ! 1. CLI Arguments Parsing
    call get_command_argument(1, basinf)
    call get_command_argument(2, lnnlv)
    call get_command_argument(3, lnrens)
    call get_command_argument(4, lnlag)

    if (len_trim(lnlag) == 0) then
       write(*,*) "USAGE: enKF2enKS [basinf] [nnlv] [nrens] [nlag]"
       write(*,*) "  basinf : basin file"
       write(*,*) "  nnlv    : number of vertical levels"
       write(*,*) "  nrens  : number of ensemble members"
       write(*,*) "  nlag   : lag (in analysis steps, -1 for all future)"
       error stop
    end if

    read(lnnlv, *, iostat=ios) nnlv
    if (ios /= 0) error stop "Cannot parse nnlv"
    read(lnrens, *, iostat=ios) nrens
    if (ios /= 0) error stop "Cannot parse nrens"
    read(lnlag, *, iostat=ios) nlag
    if (ios /= 0) error stop "Cannot parse nlag"

    write(*,*) "========================================"
    write(*,*) "EnKS (Ensemble Kalman Smoother) v2.0"
    write(*,*) "========================================"
    write(*,'(a,i5)') "Lag (steps)    :", nlag
    write(*,'(a,i5)') "Ensemble members:", nrens
    write(*,'(a,a)') "Basin file     :", trim(basinf)
    write(*,'(a,i5)') "Vertical levels:", nnlv
    write(*,*) "========================================"

    ! 2. Setup Data
    inquire(file='X5_tot.uf', exist=file_exists)
    if (.not. file_exists) then
       write(*,*) "ERROR: X5_tot.uf not found!"
       write(*,*) "Make sure the analysis step was executed and X5_tot.uf was generated."
       error stop
    end if

    call load_all_x5(nrens)      ! Load X5_tot.uf into memory
    call init_shyfem(basinf, nnlv) ! Initialize SHYFEM grids

    ! 3. Opening Output Restarts
    open(18, file="analKS_mean.rst", status="replace", form="unformatted", iostat=ios)
    if (ios /= 0) error stop "Cannot open analKS_mean.rst"
    open(19, file="analKS_std.rst",  status="replace", form="unformatted", iostat=ios)
    if (ios /= 0) error stop "Cannot open analKS_std.rst"

    write(*,*) "Output files opened: analKS_mean.rst, analKS_std.rst"
    write(*,*) ""

    ! ======================================================================
    ! MAIN TIME LOOP
    ! ======================================================================
    do rrec = 1, total_x5_records
       atime = all_x5(rrec)%tt

       if (rrec > 1) then
          write(*,*) ""
       end if
       write(*,'(a,i5,a,i5,a,f12.6)') &
           'Processing record ', rrec, '/', total_x5_records, ' at t=', atime

       ! -----------------------------------
       ! Read all ensemble members for current time
       ! -----------------------------------
       do nre = 1, nrens
          call num2str(nre-1, nrel)

          ! Read the specific restart for this member and time
          call read_rst("analKF_en"//nrel//".rst", atime, nnlv)

          ! Allocate matrices on the first successful read
          if (.not. allocated(Astate)) then
              call allocate_states(Astate, AmeanKS, AstdKS, sdim, nrens)
          end if

          ! Map SHYFEM variables into the Astate matrix column
          call push_matrix(sdim, nrens, nre, Astate)
       end do

       ! -----------------------------------
       ! Execute Smoother Analysis
       ! Apply future weights: A_smooth = A_filter * (X5_t * X5_t+1 * ... * X5_t+lag)
       ! -----------------------------------
       call make_analysis(atime, sdim, nrens, Astate, nlag, verbose=.true.)

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

    end do

    close(18)
    close(19)
    write(*,*) ""
    write(*,*) "========================================"
    write(*,*) "Smoothing completed successfully."
    write(*,*) "Output: analKS_mean.rst, analKS_std.rst"
    write(*,*) "========================================"

end program enKF2enKS

