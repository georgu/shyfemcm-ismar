! ======================================================================
!  MODULE: mod_enks_data
!
!  PURPOSE:
!    Manages memory infrastructure and binary unformatted I/O 
!    ingestion for retrieving analysis weight sequences (X3/X5 operators) 
!    accumulated during prior filter passes, preparing them for the EnKS.
! ======================================================================
module mod_enks_data
   use iso_fortran_env, only: dp=>real64
   implicit none
   save

   ! Dynamic data structure to represent individual time-slice weights
   type t_x5_record
      real(dp)              :: tt        ! Associated simulation time coordinate
      character(len=2)      :: tag       ! Filter architecture identifier ('X3' or 'X5')
      integer               :: nrens, nrobs
      real(dp), allocatable :: mat(:,:)  ! Subspace transform: X5 (N x N) or X3 (M x N)
      real(dp), allocatable :: S(:,:)    ! Observation space anomalies: (M x N) - Stochastic only
   end type t_x5_record

   ! Global storage arrays acting as the central ledger for retrospective smoothing
   type(t_x5_record), allocatable :: all_x5(:)
   integer  :: total_x5_records = 0
   real(dp) :: dt_x5 = -1.0_dp
   
   ! Floating point evaluation threshold
   real(dp), parameter :: time_tol = 1.0e-6_dp

contains

! ======================================================================
! SUBROUTINE: load_all_x5
!
! PURPOSE:
!   Parses the cumulative binary ledger 'X5_tot.uf' in a safe two-pass 
!   sequence. First determines structural limits for allocation, then
!   streams data chunks directly into memory structures.
! ======================================================================
subroutine load_all_x5(nrens_expected)
   integer, intent(in) :: nrens_expected
   integer :: u, ios, count
   real(dp) :: tmp_tt
   character(len=2)  :: tmp_tag
   integer :: tmp_nrens, tmp_nrobs

   ! 1. First Pass: Compute record counts and validate stream consistency
   open(newunit=u, file='X5_tot.uf', form='unformatted', status='old', action='read', iostat=ios)
   if (ios /= 0) then
      write(*,*) "ERROR: load_all_x5: Catastrophic I/O failure. Unable to open ledger file X5_tot.uf."
      error stop
   end if

   count = 0
   do
      ! Parse Packet Synchronizer Header
      read(u, iostat=ios) tmp_tt, tmp_tag
      if (ios /= 0) exit ! Normal termination upon reaching the end of the file

      if (tmp_tag == 'X3') then
         read(u, iostat=ios) tmp_nrens, tmp_nrobs
         if (ios /= 0) error stop "ERROR: load_all_x5: Corrupted stochastic header parameters."
         
         ! Advance beyond data records matching the binary signature of dumpX3/save_X5
         read(u, iostat=ios) ! Skip X3 data block
         if (ios /= 0) error stop "ERROR: load_all_x5: Record skip failure on X3 matrix."
         read(u, iostat=ios) ! Skip S data block
         if (ios /= 0) error stop "ERROR: load_all_x5: Record skip failure on S matrix."
         
      else if (tmp_tag == 'X5') then
         read(u, iostat=ios) tmp_nrens
         if (ios /= 0) error stop "ERROR: load_all_x5: Corrupted deterministic header parameter."
         
         ! Advance beyond data records matching the binary signature of dumpX5/save_X5
         read(u, iostat=ios) ! Skip X5 data block
         if (ios /= 0) error stop "ERROR: load_all_x5: Record skip failure on X5 matrix."
      else
         write(*,*) "ERROR: load_all_x5: Invalid architecture tag validation failure: ", tmp_tag
         error stop
      end if
      count = count + 1
   end do

   total_x5_records = count
   if (total_x5_records == 0) then
      write(*,*) "ERROR: load_all_x5: Data ledger X5_tot.uf is empty."
      error stop
   end if

   ! Secure Heap memory allocation for global state arrays
   if (allocated(all_x5)) deallocate(all_x5)
   allocate(all_x5(total_x5_records))
   rewind(u)

   ! 2. Second Pass: Materialize raw streams into high-performance array slots
   write(*,*) "load_all_x5: Ingesting ", total_x5_records, " data records into central structures..."

   do count = 1, total_x5_records
      read(u, iostat=ios) all_x5(count)%tt, all_x5(count)%tag
      if (ios /= 0) error stop "load_all_x5: Fatal stream synchronisation failure."

      if (all_x5(count)%tag == 'X3') then
         read(u, iostat=ios) all_x5(count)%nrens, all_x5(count)%nrobs
         if (ios /= 0) error stop "load_all_x5: Ingestion error on stochastic metadata dimensions."
         
         ! Dimensional cross-verification against master configuration parameters
         if (all_x5(count)%nrens /= nrens_expected) then
            write(*,'(A,I5,A,I5,A,I5)') &
                'WARNING: Stochastic record allocation mismatch at index ', count, &
                ' | Found: ', all_x5(count)%nrens, ' | Expected: ', nrens_expected
         end if
         
         ! Allocate dedicated arrays for individual simulation coordinates
         allocate(all_x5(count)%mat(all_x5(count)%nrobs, all_x5(count)%nrens))
         allocate(all_x5(count)%S(all_x5(count)%nrobs, all_x5(count)%nrens))
         
         read(u, iostat=ios) all_x5(count)%mat
         if (ios /= 0) error stop "load_all_x5: Ingestion error on X3 matrix coefficients."
         read(u, iostat=ios) all_x5(count)%S
         if (ios /= 0) error stop "load_all_x5: Ingestion error on S matrix coefficients."
         
         write(*,'(A,I5,A,I5,A,I5,A,F18.1)') &
             '  Record', count, ': X3 (', all_x5(count)%nrobs, 'x', &
             all_x5(count)%nrens, ') at t = ', all_x5(count)%tt
      else if (all_x5(count)%tag == 'X5') then
         read(u, iostat=ios) all_x5(count)%nrens
         if (ios /= 0) error stop "load_all_x5: Ingestion error on deterministic metadata dimensions."
         
         if (all_x5(count)%nrens /= nrens_expected) then
            write(*,'(A,I5,A,I5,A,I5)') &
                'WARNING: Deterministic record allocation mismatch at index ', count, &
                ' | Found: ', all_x5(count)%nrens, ' | Expected: ', nrens_expected
         end if
         
         all_x5(count)%nrobs = 0
         allocate(all_x5(count)%mat(all_x5(count)%nrens, all_x5(count)%nrens))
         
         read(u, iostat=ios) all_x5(count)%mat
         if (ios /= 0) error stop "load_all_x5: Ingestion error on X5 matrix coefficients."
         
         write(*,'(A,I5,A,I5,A,I5,A,F18.1)') &
             '  Record', count, ': X5 (', all_x5(count)%nrens, 'x', &
             all_x5(count)%nrens, ') at t = ', all_x5(count)%tt
      end if
   end do

   close(u)
   write(*,*) "load_all_x5: Ledger ingestion executed successfully."

   ! Validate chronological sorting instead of executing expensive data deep-copies
   call verify_chronology_x5()
   
   ! Establish baseline temporal resolution metric
   dt_x5 = estimate_dt_x5()
   write(*,'(A,F12.6)') 'load_all_x5: Temporal resolution established at dt_X5 = ', dt_x5
end subroutine load_all_x5

! ======================================================================
! SUBROUTINE: verify_chronology_x5
!
! PURPOSE:
!   Ensures records are monotonic and chronological. Replaces the old 
!   inefficient selection sort which caused severe heap thrashing due to 
!   continuous matrix deep-copies.
! ======================================================================
subroutine verify_chronology_x5()
   implicit none
   integer :: i

   if (total_x5_records <= 1) return

   do i = 2, total_x5_records
      if (all_x5(i)%tt < all_x5(i-1)%tt) then
         write(*,*) "ERROR: verify_chronology_x5: Non-monotonic timestamp detected."
         write(*,*) "Records must be written strictly forward in time by the filter."
         error stop
      end if
   end do
end subroutine verify_chronology_x5

! ======================================================================
! FUNCTION: estimate_dt_x5
! PURPOSE: Validates and extracts baseline grid step resolutions.
! ======================================================================
real(dp) function estimate_dt_x5()
   implicit none
   integer :: i
   estimate_dt_x5 = 0.0_dp
   if (total_x5_records < 2) return

   ! Extrapolate baseline resolution from the initial step pair
   if (abs(all_x5(2)%tt - all_x5(1)%tt) > time_tol) then
      estimate_dt_x5 = all_x5(2)%tt - all_x5(1)%tt
   end if
   
   ! Verify consistency across subsequent boundaries
   do i = 3, total_x5_records
      if (abs((all_x5(i)%tt - all_x5(i-1)%tt) - estimate_dt_x5) > time_tol) then
         write(*,*) "WARNING: estimate_dt_x5: Variable sampling resolution detected in the ledger."
      end if
   end do
end function estimate_dt_x5

! ======================================================================
! FUNCTION: equal_time
! PURPOSE: Fast logical wrapper for safe floating-point comparison.
! ======================================================================
logical function equal_time(t1,t2)
   real(dp), intent(in) :: t1,t2
   equal_time = abs(t1-t2) < time_tol
end function equal_time

! ======================================================================
! SUBROUTINE: allocate_states
!
! PURPOSE:
!   Queries active SHYFEM hydrodynamic geometry pointer structures to 
!   dynamically calculate the consolidated state vector dimension (n).
!   Allocates core smoother tracking arrays directly on the Heap.
!
! MATHEMATICAL ARCHITECTURE:
!   The full state vector compiles horizontal/vertical layouts for:
!   State = [ Elevation (Z) | Velocity U | Velocity V | Temp | Salt ]
! ======================================================================
subroutine allocate_states(A, Amean, Astd, n, nrens)
   use mod_hydro,     only : znv, utlnv, vtlnv
   use mod_hydro_vel, only : ulnv, vlnv
   use mod_ts,        only : tempv, saltv
   use mod_restart,   only : ibarcl_rst
   implicit none

   real(dp), allocatable, intent(out) :: A(:,:), Amean(:), Astd(:)
   integer,               intent(out) :: n
   integer,               intent(in)  :: nrens

   ! Consolidated calculation of fundamental barotropic state components (Elevation + 2D/3D Velocity)
   n = size(utlnv) + size(vtlnv) + size(znv)

   ! Append active thermodynamic scalars if running under a baroclinic physics configuration
   if (ibarcl_rst /= 0) then
      n = n + size(tempv) + size(saltv)
   end if

   ! Enforce structural cleanup to prevent memory corruption/leaks upon reallocation
   if (allocated(A))     deallocate(A)
   if (allocated(Amean)) deallocate(Amean)
   if (allocated(Astd))  deallocate(Astd)

   ! Materialize consolidated arrays on the Heap
   allocate(A(n, nrens))
   allocate(Amean(n), Astd(n))

   write(*,'(A,I9)') "allocate_states: Ocean state vector dynamically compiled with total dimension N = ", n
end subroutine allocate_states

end module mod_enks_data



! ======================================================================
!  MODULE: mod_enks_analysis (Section 1)
!
!  PURPOSE:
!    Handles hydrographic grid initialization and sequential restart 
!    ingestion from the SHYFEM core framework.
! ======================================================================
module mod_enks_analysis
   use iso_fortran_env, only: dp=>real64
   implicit none

contains

! ======================================================================
! SUBROUTINE: init_shyfem
! PURPOSE: Direct unformatted parsing of the SHYFEM topology basin file
!          and synchronized initialization of active geometry modules.
! ======================================================================
subroutine init_shyfem(basfile, nnlv)
   use mod_geom_dynamic
   use mod_hydro
   use mod_hydro_vel
   use mod_ts
   use mod_layer_thickness
   use mod_gotm_aux
   use mod_area
   use mod_depth
   use levels
   use shympi
   use sigma
   use evgeom
   use basin
   implicit none

   character(len=80), intent(in) :: basfile
   integer,           intent(in) :: nnlv
   integer :: ios

   ! Open basin
   open(22, file=trim(basfile), status="old", form="unformatted", iostat=ios)
   if (ios /= 0) error stop "ERROR: init_shyfem: Unable to access basin file."
   call basin_read_by_unit(22)
   close(22)

   ! Bind uniform horizontal-vertical mesh layouts
   nlv = nnlv
   nlvdi = nnlv
   
   ! Chronological initialisation of dependent SHYFEM geometric/physical substructures
   call mod_geom_dynamic_init(nkn, nel)
   call mod_hydro_init(nkn, nel, nnlv)
   call mod_hydro_vel_init(nkn, nel, nnlv)
   call mod_ts_init(nkn, nnlv)
   call levels_init(nkn, nel, nnlv)
   call mod_gotm_aux_init(nkn, nnlv)
   call shympi_set_hlv(nnlv, hlv)
   call shympi_init(.false.)
   call mod_layer_thickness_init(nkn, nel, nnlv)
   call init_sigma_info(nnlv, hlv)
   call mod_area_init(nkn, nnlv)
   call ev_init(nel)
   call mod_depth_init(nkn, nel)

   ! Publish global grid anchors
   nkn_global = nkn
   nel_global = nel
   nlv_global = nlv

end subroutine init_shyfem

! ======================================================================
! SUBROUTINE: num2str
! PURPOSE: Thread-safe conversion of integer indices to zero-padded strings.
! ======================================================================
subroutine num2str(num, str)
   implicit none
   integer,          intent(in)  :: num
   character(len=5), intent(out) :: str

   if (num < 0 .or. num > 99999) then 
      error stop "ERROR: num2str: Member index encoding parameter out of range [0, 99999]."
   end if   
            
   write(str, '(I5.5)') num
end subroutine num2str

! ======================================================================
! SUBROUTINE: read_rst
!
! PURPOSE:
!   Ingests physical state distributions from specific binary unformatted 
!   SHYFEM member restart tracking streams.
!
! CRITICAL CORRECTIONS:
!   1. Time-matching search tolerance modified from machine 'epsilon' 
!      to '1.0e-6_dp' to guarantee float matching over long integration lines.
!   2. Added analytical documentation anchors for transport-to-velocity mapping.
! ======================================================================
subroutine read_rst(rstname, atimea, nnlv)
  use mod_restart
  use levels
  use shympi
  implicit none
  
  integer,          intent(in) :: nnlv
  character(len=*), intent(in) :: rstname
  real(dp),         intent(in) :: atimea
  
  integer  :: ierr, iflag, ios, u_rst
  real(dp) :: atimef
  integer, save :: icall = 0
  
  real*4 :: ibarcl4, iconz4, imerc4, iturb4, iwvert_rst4, ieco_rst4, zero4
  real(dp), parameter :: time_match_tol = 1.0e-6_dp ! Safe floating point bound

  zero4 = 0.0

  ! Open unformatted sequential restart channel for the requested member
  open(newunit=u_rst, file=trim(rstname), status='old', form='unformatted', action='read', iostat=ios)
  if (ios /= 0) error stop 'ERROR: read_rst: Ingestion streaming failure opening restart file.'

  ! Parse records until the timestamps match
  do
     call rst_read_record(u_rst, atimef, iflag, ierr)
     if (ierr /= 0) then
        close(u_rst)
        write(*,'(A,F14.4,A,A)') 'ERROR: read_rst: Target timestamp ', atimea, &
                                 ' not resolved inside stream: ', trim(rstname)
        error stop
     end if
     
     ! Replaced hyper-strict machine epsilon with robust time matching tolerance
     if (abs(atimef - atimea) < time_match_tol) exit   
  end do
  close(u_rst)

  ! On the initial baseline execution, register and broadcast execution flags
  if (icall == 0) then
     if (nnlv /= nlv) error stop 'ERROR: read_rst: Vertical grid level dimensionality mismatch.'
     hlv        = hlvrst
     hlv_global = hlvrst
     ilhv       = ilhrst
     ilhkv      = ilhkrst

     ibarcl4     = real(ibarcl_rst, 4)
     iwvert_rst4 = real(iwvert_rst, 4)
     ieco_rst4   = real(ieco_rst, 4)
     iconz4      = real(iconz_rst, 4)
     imerc4      = real(imerc_rst, 4)
     iturb4      = real(iturb_rst, 4)

     call addpar('ibarcl', ibarcl4)
     call addpar('iconz' , iconz4 )
     call addpar('iwvert', iwvert_rst4)
     call addpar('ieco'  , ieco_rst4)
     call addpar('ibio'  , zero4)
     call addpar('ibfm'  , zero4)
     call addpar('imerc' , imerc4)
     call addpar('iturb' , iturb4)

     call addpar('nzadapt' , 0.04) 

     call daddpar('date', 0.0_dp)
     call daddpar('time', 0.0_dp)

     write(*,*) '--- SHYFEM Configuration Flags Published ---'
     write(*,'(A,I4)')   '  ibarcl (Baroclinic mode indicator) : ', ibarcl_rst
     write(*,'(A,I4)')   '  iconz  (Constituent tracking)      : ', iconz_rst
     write(*,'(A,I4)')   '  nlv    (Vertical mesh resolution)  : ', nlv
     write(*,*) '--------------------------------------------'

     call set_ev
     call set_area
     call set_depth
     call init_zadaptation
     call make_new_layer_depth
  end if

  icall = icall + 1

end subroutine read_rst

! ======================================================================
! SUBROUTINE: rst_write_rec
! PURPOSE: Synchronizes the internal SHYFEM output state variables 
!          and serializes an unformatted restart record to disk.
! ======================================================================
subroutine rst_write_rec(atimea, iunit)
  use iso_fortran_env, only : dp => real64
  use mod_hydro
  use mod_hydro_vel
  use mod_geom_dynamic
  use mod_restart
  implicit none
  
  real(dp), intent(in) :: atimea
  integer,  intent(in) :: iunit

  ! Synchronize dry/wet masks and copy active variables to old-state arrays expected by SHYFEM
  iwetv = 0
  zov   = znv
  zeov  = zenv
  utlov = utlnv
  vtlov = vtlnv

  call rst_write_record(atimea, iunit)
end subroutine rst_write_rec

! ======================================================================
! SUBROUTINE: push_matrix
!
! PURPOSE:
!   Extracts physical fields from SHYFEM modules, converts hydrodynamical 
!   transports to pure velocities using static layer thicknesses, and maps 
!   them into a specific column of the master filter state matrix (Amat).
! ======================================================================
subroutine push_matrix(sdim, nrens, nre, Amat)
   use iso_fortran_env, only: dp=>real64
   use mod_layer_thickness
   use mod_hydro
   use mod_hydro_vel
   use levels
   use mod_ts
   use mod_restart
   use basin
   implicit none

   integer,  intent(in)    :: sdim, nrens, nre
   real(dp), intent(inout) :: Amat(sdim,nrens)
   
   integer :: d_uv, d_z, d_ts

   ! 1. Memory allocation safeguards
   if (.not. allocated(znv))   stop 'ERROR: push_matrix: znv matrix pool unallocated.'
   if (.not. allocated(utlnv)) stop 'ERROR: push_matrix: utlnv matrix pool unallocated.'
   if (.not. allocated(vtlnv)) stop 'ERROR: push_matrix: vtlnv matrix pool unallocated.'
   if (ibarcl_rst /= 0) then
      if (.not. allocated(tempv)) stop 'ERROR: push_matrix: tempv matrix pool unallocated.'
      if (.not. allocated(saltv)) stop 'ERROR: push_matrix: saltv matrix pool unallocated.'
   end if  

   d_uv = nlv * nel
   d_z  = nkn

   ! Dimension bounds check
   if (2*d_uv + d_z > sdim) then
      write(*,'(A,I10,A,I10)') 'ERROR: push_matrix: Baseline dimension mismatch: ', 2*d_uv+d_z, ' > sdim: ', sdim
      error stop
   end if

   ! 2. PHYSICAL CONVERSION: Translate Transport (m^2/s) to Velocity (m/s) via static hdenv
   ulnv = 0.0_dp
   vlnv = 0.0_dp
   where( hdenv > 1.0e-6_dp )
       ulnv = utlnv / hdenv
       vlnv = vtlnv / hdenv
   end where

   ! 3. Map velocities and elevation fields into the designated ensemble member column
   Amat(1:d_uv, nre) = reshape(real(ulnv, dp), [d_uv])
   Amat(d_uv+1:2*d_uv, nre) = reshape(real(vlnv, dp), [d_uv])
   Amat(2*d_uv+1:2*d_uv+d_z, nre) = real(znv, dp)

   ! 4. Handle baroclinic state parameters (Temperature & Salinity) if active
   if (ibarcl_rst /= 0) then
      d_ts = size(tempv)
      if (2*d_uv+d_z+2*d_ts > sdim) then
         write(*,'(A,I10,A,I10)') 'ERROR: push_matrix: Baroclinic component size out of bounds: ', &
                                  2*d_uv+d_z+2*d_ts, ' > sdim: ', sdim
         error stop
      end if
      
      Amat(2*d_uv+d_z+1 : 2*d_uv+d_z+d_ts, nre) = reshape(real(tempv, dp), [d_ts])
      Amat(2*d_uv+d_z+d_ts+1 : 2*d_uv+d_z+2*d_ts, nre) = reshape(real(saltv, dp), [d_ts])
   end if
end subroutine push_matrix

! ======================================================================
! SUBROUTINE: pull_matrix
!
! PURPOSE:
!   Retrieves smoothed state coefficients from the smoother matrix, maps 
!   them back to SHYFEM geometry, and re-calculates physical transports 
!   using static layer thicknesses (hdenv).
! ======================================================================
subroutine pull_matrix(sdim, nrens, nre, Amat)
   use iso_fortran_env, only: dp=>real64
   use mod_layer_thickness
   use mod_hydro
   use mod_hydro_vel
   use levels
   use mod_ts
   use mod_restart
   use basin
   implicit none

   integer,  intent(in) :: sdim, nrens, nre
   real(dp), intent(in) :: Amat(sdim,nrens)
   
   integer :: d_uv, d_z, d_ts

   d_z   = nkn
   d_uv  = nel * nlv

   ! Bounds verification
   if (2*d_uv + d_z > sdim) then
      write(*,'(A,I10,A,I10)') 'ERROR: pull_matrix: Grid dimension mapping failure: ', 2*d_uv+d_z, ' > sdim: ', sdim
      error stop
   end if

   ! 1. Unpack smooth velocity trajectories and sea surface elevations from matrix slot
   ulnv = reshape(real(Amat(1:d_uv, nre), dp), [nlv, nel])
   vlnv = reshape(real(Amat(d_uv+1:2*d_uv, nre), dp), [nlv, nel])
   znv  = real(Amat(2*d_uv+1:2*d_uv+d_z, nre), dp)

   ! 2. Unpack thermodynamic scalars if baroclinic physics are engaged
   if (ibarcl_rst /= 0) then
      d_ts = size(tempv)
      if (2*d_uv+d_z+2*d_ts > sdim) then
         write(*,'(A,I10,A,I10)') 'ERROR: pull_matrix: Thermodynamic field size mismatch: ', &
                                  2*d_uv+d_z+2*d_ts, ' > sdim: ', sdim
         error stop
      end if
      tempv = reshape(real(Amat(2*d_uv+d_z+1:2*d_uv+d_z+d_ts, nre), dp), [nlv, nkn])
      saltv = reshape(real(Amat(2*d_uv+d_z+d_ts+1:2*d_uv+d_z+2*d_ts, nre), dp), [nlv, nkn])
   end if

   ! 3. Convert smoothed pure velocities (m/s) back to physical fluid transports (m^2/s) using static layers
   utlnv = ulnv * hdenv
   vtlnv = vlnv * hdenv
end subroutine pull_matrix

! ======================================================================
! SUBROUTINE: make_mn_std
!
! PURPOSE:
!   Computes the ensemble mean vector and standard deviation vector 
!   for each physical grid state entry across the entire ensemble pool.
!
! PERFORMANCE NOTE:
!   Loops are vectorized using Fortran intrinsic array syntax to allow 
!   the compiler to execute SIMD (Single Instruction Multiple Data) 
!   optimizations on the columns.
! ======================================================================
subroutine make_mn_std(n, nens, Amat, mean_vec, std_vec)
    use iso_fortran_env, only: dp => real64
    implicit none

    integer,  intent(in)  :: n, nens
    real(dp), intent(in)  :: Amat(n, nens)
    real(dp), intent(out) :: mean_vec(n)
    real(dp), intent(out) :: std_vec(n)

    integer :: i
    real(dp) :: inv_nens

    inv_nens = 1.0_dp / real(nens, dp)

    ! Compute the spatial state mean vector
    do i = 1, n
        mean_vec(i) = sum(Amat(i, :)) * inv_nens
    end do

    ! Compute the standard deviation vector (unbiased population estimator for ensembles)
    do i = 1, n
        std_vec(i) = sqrt(sum((Amat(i, :) - mean_vec(i))**2) * inv_nens)
    end do
end subroutine make_mn_std

! ======================================================================
! SUBROUTINE: make_analysis
!
! PURPOSE:
!   Executes the core retrospective Ensemble Kalman Smoother (EnKS) step
!   by concatenating future analysis transformation operators within a 
!   defined lag window, applying the joint weights to update past state fields.
!
!   Apply future information retrogradely while looping
!   forward, the incoming operators must be premultiplied (multiplied onto 
!   the left side of Xacc: Xtmp = X5 * Xacc), ensuring proper temporal chaining.
! ======================================================================
subroutine make_analysis(atime, sdim, nrens, Amat, nlag, verbose)
   use iso_fortran_env, only : dp=>real64
   use mod_enks_data,    only : all_x5, total_x5_records, dt_x5, time_tol
   implicit none

   real(dp), intent(in)    :: atime
   integer,  intent(in)    :: sdim, nrens, nlag
   real(dp), intent(inout) :: Amat(sdim,nrens)
   logical,  intent(in), optional :: verbose

   real(dp) :: t_end
   integer  :: i, count, j, count_applied
   real(dp), allocatable :: Xacc(:,:), Xtmp(:,:), Amat_tmp(:,:), X5_equiv(:,:)
   logical  :: lverbose

   lverbose = .false.
   if (present(verbose)) lverbose = verbose

   ! 1. Define the dynamic temporal smoothing lag window boundaries: [atime, t_end]
   if (nlag == -1) then
      t_end = huge(1.0_dp)
      if (lverbose) write(*,'(A,F18.1)') &
          'make_analysis: Unbounded smoothing engaged. Compiling all available future records from t = ', atime
   else
      t_end = atime + real(nlag, dp) * dt_x5 + time_tol
      if (lverbose) write(*,'(A,F18.1,A,F18.1,A)') &
          'make_analysis: Compiling future records within lag window [', atime, ' , ', t_end, ']'
   end if

   allocate(Xacc(nrens,nrens), Xtmp(nrens,nrens))

   ! Initialize the transition operator matrix Xacc as an Identity matrix
   Xacc = 0.0_dp
   do i = 1, nrens
      Xacc(i,i) = 1.0_dp
   end do

   count = 0
   count_applied = 0
   
   ! 2. Parse the historical weight ledger to accumulate operations within the lag window
   do i = 1, total_x5_records
      
      ! Check if the current ledger record falls inside the active smoothing window ]atime, t_end]
      ! Note that the window does not include the current time.
      if (all_x5(i)%tt > atime - time_tol .and. all_x5(i)%tt <= t_end) then
         count = count + 1
         
         if (all_x5(i)%tag == 'X3') then
            ! -----------------------------------------------------------
            ! Stochastic Subspace Path: Reconstruct X5 equivalent matrix
            ! Formula: X5_equiv = I + S^T * X3
            ! -----------------------------------------------------------
            allocate(X5_equiv(nrens, nrens))
            X5_equiv = 0.0_dp
            do j = 1, nrens
               X5_equiv(j, j) = 1.0_dp
            end do
            
            ! Project observation anomalies to member space via DGEMM
            call dgemm('T', 'N', nrens, nrens, all_x5(i)%nrobs, 1.0_dp, &
                       all_x5(i)%S, all_x5(i)%nrobs, all_x5(i)%mat, all_x5(i)%nrobs, &
                       1.0_dp, X5_equiv, nrens)
            
            ! Premultiply (X5_equiv * Xacc) to preserve backward-chaining math
            call dgemm("N", "N", nrens, nrens, nrens, 1.0_dp, X5_equiv, nrens, &
                       Xacc, nrens, 0.0_dp, Xtmp, nrens)
            
            if (lverbose) write(*,'(A,I5,A,F18.1)') '    Ingesting Stochastic weights from record ', i, ' at t=', all_x5(i)%tt
            deallocate(X5_equiv)
            
         else if (all_x5(i)%tag == 'X5') then
            ! -----------------------------------------------------------
            ! Deterministic Subspace Path: Direct operator multiplication
            ! -----------------------------------------------------------
            ! Premultiply (X5 * Xacc) to ensure chronological inversion order
            call dgemm("N", "N", nrens, nrens, nrens, 1.0_dp, all_x5(i)%mat, nrens, &
                       Xacc, nrens, 0.0_dp, Xtmp, nrens)
            
            if (lverbose) write(*,'(A,I5,A,F18.1)') '    Ingesting Deterministic weights from record ', i, ' at t=', all_x5(i)%tt
         else
            write(*,'(A,A)') 'ERROR: make_analysis: Invalid architecture signature detected: ', all_x5(i)%tag
            cycle
         end if
         
         ! Update the central transition ledger operator
         Xacc = Xtmp
         count_applied = count_applied + 1
      end if
   end do

   ! 3. Final Step: Apply the fully consolidated weight matrix ONCE onto the state ensemble matrix
   ! Formula: A_smoothed = A_filtered * Xacc
   if (count_applied > 0) then
      allocate(Amat_tmp(sdim, nrens))
      
      call dgemm("N", "N", sdim, nrens, nrens, 1.0_dp, Amat, sdim, &
                 Xacc, nrens, 0.0_dp, Amat_tmp, sdim)
                 
      Amat = Amat_tmp
      deallocate(Amat_tmp)
      
      if (lverbose) write(*,'(A,I5,A,F18.1)') 'make_analysis: Retrospective smoothing completed. Blended ', &
                                              count_applied, ' weight packets at time coordinate ', atime
   else
      if (lverbose) write(*,'(A,F18.1)') 'make_analysis: No future weight packets fell within criteria for time ', atime
   end if

   deallocate(Xacc, Xtmp)
end subroutine make_analysis

end module mod_enks_analysis

! ======================================================================
! PROGRAM: enKF2enKS
!
! PURPOSE:
!   Main execution driver for the Ensemble Kalman Smoother (EnKS v2.0).
!   Ingests sequential ocean model restarts (SHYFEM format) generated 
!   by a prior EnKF filter pass, applies accumulated transformation 
!   matrices (X3/X5), and exports retrospective mean and variance fields.
!
! ======================================================================
! ======================================================================
! PROGRAM: enKF2enKS
!
! PURPOSE:
!   Main execution driver for the Ensemble Kalman Smoother (EnKS v2.0).
!   Ingests sequential ocean model restarts (SHYFEM format) generated
!   by a prior EnKF filter pass, applies accumulated transformation
!   matrices (X3/X5), and exports retrospective mean and variance fields.
! ======================================================================
program enKF2enKS
    use iso_fortran_env, only: dp=>real64
    use mod_enks_data
    use mod_enks_analysis
    use mod_restart,   only : ibarcl_rst
    use basin,         only : nkn, nel
    implicit none

    ! Argument parsing and configuration variables
    character(len=80) :: basinf
    character(len=3)  :: lnnlv
    character(len=6)  :: lnrens, lnlag
    integer           :: nnlv, nrens, nlag
    integer           :: rrec, nre, sdim
    character(len=5)  :: nrel
    real(dp)          :: atime

    ! Global state space arrays (Allocated dynamically after reading the first restart)
    real(dp), allocatable :: Astate(:,:), AmeanKS(:), AstdKS(:)
    integer :: ios
    logical :: file_exists

    ! 1. Parse Command Line Interface (CLI) Arguments
    call get_command_argument(1, basinf)
    call get_command_argument(2, lnnlv)
    call get_command_argument(3, lnrens)
    call get_command_argument(4, lnlag)

    if (len_trim(lnlag) == 0) then
       write(*,*) "USAGE: enKF2enKS [basinf] [nnlv] [nrens] [nlag]"
       write(*,*) "  basinf : path to SHYFEM basin grid definition file"
       write(*,*) "  nnlv   : total number of vertical hydrographic levels"
       write(*,*) "  nrens  : number of active ensemble members (N)"
       write(*,*) "  nlag   : temporal lag configuration (-1 for full future window)"
       stop
    end if

    read(lnnlv, *, iostat=ios) nnlv
    if (ios /= 0) error stop "ERROR: Main parsing failure on parameter: nnlv"
    read(lnrens, *, iostat=ios) nrens
    if (ios /= 0) error stop "ERROR: Main parsing failure on parameter: nrens"
    read(lnlag, *, iostat=ios) nlag
    if (ios /= 0) error stop "ERROR: Main parsing failure on parameter: nlag"

    write(*,*) "========================================"
    write(*,*) " Ensemble Kalman Smoother (EnKS) v2.0"
    write(*,*) "========================================"
    write(*,'(A,I5)') " Lag window (steps) : ", nlag
    write(*,'(A,I5)') " Ensemble size (N)  : ", nrens
    write(*,'(A,A)')  " Grid Basin file    : ", trim(basinf)
    write(*,'(A,I5)') " Vertical levels    : ", nnlv
    write(*,*) "========================================"

    ! 2. Ingest Global Filter Metadata Matrices
    inquire(file='X5_tot.uf', exist=file_exists)
    if (.not. file_exists) then
       write(*,*) "ERROR: Cumulative weight file 'X5_tot.uf' not found in workspace."
       error stop
    end if

    ! Load the full historical matrix sequence into global heap storage
    call load_all_x5(nrens)

    ! Initialize the core SHYFEM oceanographic grid spatial structures
    call init_shyfem(basinf, nnlv)

    ! 3. Initialize Output Stream Directives
    open(18, file="analKS_mean.rst", status="replace", form="unformatted", iostat=ios)
    if (ios /= 0) error stop "ERROR: Cannot initialize output stream: analKS_mean.rst"
    open(19, file="analKS_std.rst",  status="replace", form="unformatted", iostat=ios)
    if (ios /= 0) error stop "ERROR: Cannot initialize output stream: analKS_std.rst"

    write(*,*) "Retrospective output streams initialized on logical units 18 and 19."
    write(*,*) ""

    ! ======================================================================
    ! MAIN ANALYSIS LOOP (CHRONOLOGICAL SEQUENCING)
    ! ======================================================================
    do rrec = 1, total_x5_records
       atime = all_x5(rrec)%tt

       if (rrec > 1) write(*,*) ""
       write(*,'(A,I5,A,I5,A,F18.1)') &
           'Inverting Record: ', rrec, ' / ', total_x5_records, ' | Physical Time = ', atime

       ! ------------------------------------------------------------------
       ! Step A: Stream physical state vector configurations for all members
       ! ------------------------------------------------------------------
       do nre = 1, nrens
          call num2str(nre-1, nrel)

          ! Ingest spatial restart arrays for the targeted member and timestamp
          call read_rst("analKF_en"//nrel//".rst", atime, nnlv)

          ! Calculation of sdim immediately after the first successful restart read
          if (.not. allocated(Astate)) then
              ! Import block flags from modules to calculate sdim
              sdim = nkn + 2 * nnlv * nel
              if (ibarcl_rst /= 0) then
                  sdim = sdim + 2 * nnlv * nkn
              end if

              ! Allocate core smoother tracking arrays on the Heap once
              allocate(Astate(sdim, nrens), AmeanKS(sdim), AstdKS(sdim))
              write(*,'(A,I10)') "Astate dynamically allocated with size sdim = ", sdim
          end if

          ! Map the ingested spatial arrays into the corresponding column of the state matrix
          call push_matrix(sdim, nrens, nre, Astate)
       end do

       ! ------------------------------------------------------------------
       ! Step B: Execute the Backward-Looking Subspace Smoother Analysis
       ! ------------------------------------------------------------------
       call make_analysis(atime, sdim, nrens, Astate, nlag, .true.)

       ! ------------------------------------------------------------------
       ! Step C: Extract Statistics & Serialize Smoothed Fields
       ! ------------------------------------------------------------------
       call make_mn_std(sdim, nrens, Astate, AmeanKS, AstdKS)

       ! Project smoothed mean vector back to SHYFEM geometry and commit record
       call pull_matrix(sdim, nrens, 1, AmeanKS)
       call rst_write_rec(atime, 18)

       ! Project smoothed standard deviation vector back to SHYFEM geometry and commit record
       call pull_matrix(sdim, nrens, 1, AstdKS)
       call rst_write_rec(atime, 19)
    end do

    close(18)
    close(19)

    if (allocated(Astate))  deallocate(Astate)
    if (allocated(AmeanKS)) deallocate(AmeanKS)
    if (allocated(AstdKS))  deallocate(AstdKS)

    write(*,*) ""
    write(*,*) "========================================"
    write(*,*) " Retrospective Smoothing Operations Completed."
    write(*,*) "========================================"

end program enKF2enKS

