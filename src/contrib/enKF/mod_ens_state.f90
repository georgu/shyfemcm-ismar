! ======================================================================
!  MODULE: mod_ens_state
!
!  PURPOSE: Manage ensemble of states (background, analysis, mean, std)
!           with physical validation, boundary corrections, and conversions
!           between state types and matrix representations.
!
!  IMPROVEMENTS (v2.0):
!    - Better numerical stability in increment bounds checking
!    - Proper array dimension validation
!    - Enhanced diagnostics for physical constraint violations
!    - Fixed rank issues in allocations
!    - Safe handling of optional parameters
!    - Thread-safe boundary corrections
! ======================================================================

module mod_ens_state

   use iso_fortran_env, only: dp => real64
   use mod_init_enkf
   use mod_mod_states
   use mod_para
   use basin
   use shympi, only: shympi_set_hlv, shympi_init
   use levels, only: nlv, hlv, levels_init

   implicit none

   ! Ensemble of states
   type(states), allocatable :: Abk(:), Aan(:)
   type(states) :: Abk_m, Aan_m
   type(states) :: Abk_std, Aan_std

contains

!=======================================================================
subroutine read_ensemble()
! Read all ensemble restart files and initialize model fields
!=======================================================================
   use mod_geom_dynamic
   use mod_hydro
   use mod_hydro_vel
   use mod_ts
   use mod_conz
   use mod_gotm_aux
   use mod_restart
   implicit none

   character(len=5) :: nrel, nal
   character(len=200) :: rstname
   integer :: ne
   logical :: bexist

   !-----------------------------
   ! Read basin file
   !-----------------------------
   open(21, file=basfile, status='old', form='unformatted')
   call basin_read_by_unit(21)
   close(21)

   ! check dimensions
   if (( nnkn /= nkn ).or.( nnel /= nel )) &
      error stop 'Horizontal dimensions of restart and basin differ.'

   ! Initialize SHYFEM modules BEFORE reading restart
   nlv  = nnlv
   call mod_geom_dynamic_init(nnkn, nnel)
   call mod_hydro_init(nnkn, nnel, nnlv)
   call mod_hydro_vel_init(nnkn, nnel, nnlv)
   call mod_ts_init(nnkn, nnlv)
   call levels_init(nnkn, nnel, nnlv)
   call mod_gotm_aux_init(nnkn, nnlv)
   call shympi_set_hlv(nnlv, hlv)
   call shympi_init(.false.)

   ! Allocate ensemble
   call allocate_all()

   call num2str(nanal, nal)

   write(*,*) 'Loading an ensemble of initial states (RSTs)...'
   do ne = 1, nrens
      call num2str(ne-1, nrel)
      rstname = 'an'//nal//'_'//'en'//nrel//'b.rst'
      inquire(file=rstname, exist=bexist)
      if( .not. bexist ) then
        write(6,*) 'restart file does not exists: ',trim(rstname)
        error stop 'error stop read_ensemble: no restart file'
      end if
      call read_state(Abk(ne), rstname)
   end do

end subroutine read_ensemble
!=======================================================================

!=======================================================================
! Write the whole analyzed ensemble to restart files
!=======================================================================
subroutine write_ensemble()
   implicit none

   character(len=5)  :: nrel, nal
   character(len=200):: rstname
   integer :: ne

   write(*,*) 'Writing analysis restart files...'
   call num2str(nanal, nal)

   do ne = 1, nrens
      call num2str(ne-1, nrel)
      rstname = 'an'//nal//'_'//'en'//nrel//'a.rst'
      call write_state(Aan(ne), rstname)
   end do

   rstname = 'an'//nal//'_mean_a.rst'
   call write_state(Aan_m, rstname)

   rstname = 'an'//nal//'_std_a.rst'
   call write_state(Aan_std, rstname)
end subroutine write_ensemble

!=======================================================================
! Check and correct ensemble fields near boundaries and ensure values
! remain physical (no NaN, no large increments, no out-of-range).
!
! IMPROVEMENTS (v2.0):
!   - Better increment bounds based on variable type
!   - Explicit diagnostics for each constraint type
!   - Safer handling of edge cases (small background values)
!   - OpenMP reduction for thread-safe aggregation
! ======================================================================
subroutine bc_val_check_correct()
   implicit none

   integer :: ne, k, nl, ie
   integer :: zout, uvout, sout, tout
   integer :: znan, uvnan, snan, tnan
   integer :: zbig, uvbig, sbig, tbig
   integer :: otot, ntot, btot
   logical :: file_exists
   integer :: nbc
   integer, allocatable :: bcid(:)
   real(dp), allocatable :: bcrho(:)
   real(dp) :: w

   ! Physical limits for model stability
   real(dp), parameter :: UV_MAX = 3.0_dp
   real(dp), parameter :: UV_MIN = -3.0_dp

   nbc = 1
   inquire(file='lbound.dat', exist=file_exists)

   if (file_exists) then
      allocate(bcid(nbc), bcrho(nbc))
      call read_bc_file(0, 'lbound.dat', nbc, bcid, bcrho)
      deallocate(bcid, bcrho)
      allocate(bcid(nbc), bcrho(nbc))
      call read_bc_file(1, 'lbound.dat', nbc, bcid, bcrho)
      write(*,*) 'Boundary value correction active.'
   else
      write(*,*) 'No boundary correction applied.'
      allocate(bcid(1), bcrho(1))
      bcid = 0
      bcrho = 0.0_dp
   end if

   otot = 0 ; ntot = 0 ; btot = 0

   !===================================================================
   ! Loop over ensemble members - OMP parallelized at highest level
   !===================================================================
!$OMP PARALLEL DO PRIVATE(ne, k, nl, ie, w) &
!$OMP PRIVATE(znan, zout, zbig, snan, sout, sbig, tnan, tout, tbig) &
!$OMP PRIVATE(uvnan, uvout, uvbig) &
!$OMP SHARED(nrens, Abk, Aan, file_exists, nbc, bcid, bcrho) &
!$OMP SHARED(nnkn, nnel, nnlv) &
!$OMP REDUCTION(+:otot, ntot, btot)
   do ne = 1, nrens

      zout = 0 ; uvout = 0 ; sout = 0 ; tout = 0
      znan = 0 ; uvnan = 0 ; snan = 0 ; tnan = 0
      zbig = 0 ; uvbig = 0 ; sbig = 0 ; tbig = 0

      !===============================================================
      ! NODE-BASED FIELDS: z, T, S
      !===============================================================
      do k = 1, nnkn

         ! Compute spatial boundary weight
         if (file_exists) then
            call bc_correction('node', k, nbc, bcid, bcrho, w)
         else
            w = 0.0_dp
         end if

         ! Apply boundary relaxation weight BEFORE physical check
         Aan(ne)%z(k) = w * Abk(ne)%z(k) + (1.0_dp - w) * Aan(ne)%z(k)

         ! Validate physical thresholds and metric increments for SSH
         call check_one_val( &
            vtype = 'Z', &
            va = Aan(ne)%z(k), &
            vb = Abk(ne)%z(k), &
            vmax = SSH_MAX, &
            vmin = SSH_MIN, &
            vnan = znan, vout = zout, vbig = zbig )

         !-------------------------
         ! Temperature & salinity
         !-------------------------
         do nl = 1, nnlv

            ! Apply boundary relaxation to Salinity
            Aan(ne)%s(nl,k) = w * Abk(ne)%s(nl,k) + (1.0_dp - w) * Aan(ne)%s(nl,k)

            ! Validate Salinity
            call check_one_val( &
               vtype = 'S', &
               va = Aan(ne)%s(nl,k), &
               vb = Abk(ne)%s(nl,k), &
               vmax = SAL_MAX, vmin = SAL_MIN, &
               vnan = snan, vout = sout, vbig = sbig)

            ! Apply boundary relaxation to Temperature
            Aan(ne)%t(nl,k) = w * Abk(ne)%t(nl,k) + (1.0_dp - w) * Aan(ne)%t(nl,k)

            ! Validate Temperature
            call check_one_val( &
               vtype = 'T', &
               va = Aan(ne)%t(nl,k), &
               vb = Abk(ne)%t(nl,k), &
               vmax = TEM_MAX, vmin = TEM_MIN, &
               vnan = tnan, vout = tout, vbig = tbig)
         end do

      end do

      !===============================================================
      ! ELEMENT-BASED FIELDS: u, v (Velocity Vectors)
      !===============================================================
      do ie = 1, nnel

         if (file_exists) then
            call bc_correction('elem', ie, nbc, bcid, bcrho, w)
         else
            w = 0.0_dp
         end if

         do nl = 1, nnlv
            ! Apply boundary relaxation to U and V currents
            Aan(ne)%u(nl,ie) = w * Abk(ne)%u(nl,ie) + (1.0_dp - w) * Aan(ne)%u(nl,ie)
            Aan(ne)%v(nl,ie) = w * Abk(ne)%v(nl,ie) + (1.0_dp - w) * Aan(ne)%v(nl,ie)

            ! PHYSICAL FIX: Validate currents to suppress non-physical kinetic energy spikes
            call check_one_val( &
               vtype = 'V', &
               va = Aan(ne)%u(nl,ie), &
               vb = Abk(ne)%u(nl,ie), &
               vmax = UV_MAX, vmin = UV_MIN, &
               vnan = uvnan, vout = uvout, vbig = uvbig)

            call check_one_val( &
               vtype = 'V', &
               va = Aan(ne)%v(nl,ie), &
               vb = Abk(ne)%v(nl,ie), &
               vmax = UV_MAX, vmin = UV_MIN, &
               vnan = uvnan, vout = uvout, vbig = uvbig)
         end do
      end do

      ! Safe reduction accumulation at the end of each member loop
      otot = otot + zout + sout + tout + uvout
      ntot = ntot + znan + snan + tnan + uvnan
      btot = btot + zbig + sbig + tbig + uvbig

   end do
!$OMP END PARALLEL DO

   if (file_exists) deallocate(bcid, bcrho)

   if (otot > 0) write(*,'(a,i8)') 'Total out-of-range corrections:', otot
   if (ntot > 0) write(*,'(a,i8)') 'Total NaN corrections: ', ntot
   if (btot > 0) write(*,'(a,i8)') 'Total excessive increment corrections:', btot

end subroutine bc_val_check_correct
!=======================================================================

!=======================================================================
! Read boundary condition file (two-stage read)
! icall==0 -> read number of BC nodes
! icall==1 -> read list of BC nodes and radii
!=======================================================================
subroutine read_bc_file(icall, bcfile, nbc, bcid, bcrho)
   implicit none
   character(len=*), intent(in) :: bcfile
   integer, intent(in) :: icall
   integer, intent(inout) :: nbc
   integer, intent(out) :: bcid(nbc)
   real(dp), intent(out) :: bcrho(nbc)

   integer :: i, ios

   if (icall == 0) then
      open(28, file=trim(bcfile), status='old', iostat=ios)
      if (ios /= 0) error stop "read_bc_file: Cannot open boundary file"
      read(28, *, iostat=ios) nbc
      if (ios /= 0) error stop "read_bc_file: Cannot read nbc"
      close(28)
      return
   end if

   open(28, file=trim(bcfile), status='old', iostat=ios)
   if (ios /= 0) error stop "read_bc_file: Cannot open boundary file"
   read(28, *, iostat=ios) nbc
   if (ios /= 0) error stop "read_bc_file: Cannot read nbc"
   
   do i = 1, nbc
      read(28, *, iostat=ios) bcid(i), bcrho(i)
      if (ios /= 0) then
         write(*,'(a,i5)') "read_bc_file: Error reading line", i
         error stop
      end if
   end do
   close(28)

end subroutine read_bc_file
!=======================================================================

!=======================================================================
! Compute weight for boundary damping based on distance from BC node.
! Thread-safe: accesses global arrays in a read-only manner.
!=======================================================================
subroutine bc_correction(stype, id, nbc, bcid, bcrho, w)
   implicit none

   character(len=*), intent(in) :: stype
   integer, intent(in) :: id
   integer, intent(in) :: nbc
   integer, intent(in) :: bcid(nbc)
   real(dp), intent(in) :: bcrho(nbc)
   real(dp), intent(out) :: w

   integer :: i, k
   real(dp) :: x, y, bcx, bcy, d, dmin, rho

   x = 0.0_dp ; y = 0.0_dp
   dmin = 1.0e15_dp
   rho = 0.0_dp

   if (stype == 'node') then
      x = xgv(id)
      y = ygv(id)
   else if (stype == 'elem') then
      ! element -> average vertices
      do i = 1, 3
         k = nen3v(i, id)
         x = x + xgv(k)
         y = y + ygv(k)
      end do
      x = x / 3.0_dp
      y = y / 3.0_dp
   end if

   ! Find closest boundary condition node
   do i = 1, nbc
      if (bcid(i) <= 0) cycle  ! Skip invalid BC IDs
      
      k = bcid(i)
      ! SAFETY: Check bounds
      if (k < 1 .or. k > size(xgv)) then
         write(*,'(a,i8)') 'WARNING: bc_correction: Invalid BC node ID:', k
         cycle
      end if
      
      bcx = xgv(k)
      bcy = ygv(k)
      d = sqrt((x-bcx)**2 + (y-bcy)**2)

      if (d < dmin) then
         dmin = d
         rho = bcrho(i)
      end if
   end do

   call find_weight_GC(rho, dmin, w)

end subroutine bc_correction
!=======================================================================

!=======================================================================
! Check a single scalar value (Thread-safe, no shared variables)
!
! IMPROVEMENTS (v2.0):
!   - Better relative increment logic for tracers
!   - Clearer distinction between SSH, velocity, and tracer bounds
!   - Safer handling of small background values
!=======================================================================
subroutine check_one_val(vtype, va, vb, vmax, vmin, vnan, vout, vbig)
   implicit none
   character(len=1), intent(in) :: vtype
   real(dp), intent(inout) :: va
   real(dp), intent(in) :: vb, vmin, vmax
   integer, intent(inout) :: vnan, vout, vbig

   real(dp), parameter :: max_rel_inc_tracer = 0.8_dp  ! 80% for T and S
   real(dp), parameter :: max_abs_ssh = 0.4_dp         ! Max 40 cm for SSH
   real(dp), parameter :: max_abs_uv  = 0.5_dp         ! Max 50 cm/s for currents
   real(dp), parameter :: small_ref = 1.0e-6_dp        ! Reference for small values
   
   real(dp) :: inc, rel_inc, scale_factor

   ! 1. Remove NaNs
   if (va /= va) then
      va = vb
      vnan = vnan + 1
      return
   end if

   ! 2. Bound increments using variable-specific physics
   inc = va - vb

   select case (vtype)
   case ('Z')
      ! SSH: Metric bounds to avoid gravity wave triggers
      if (abs(inc) > max_abs_ssh) then
         va = vb + sign(max_abs_ssh, inc)
         vbig = vbig + 1
      end if
      
   case ('V')
      ! VELOCITY: Metric absolute bounds to secure CFL stability
      if (abs(inc) > max_abs_uv) then
         va = vb + sign(max_abs_uv, inc)
         vbig = vbig + 1
      end if
      
   case default
      ! T & S: Proportional scaling for tracer consistency
      ! IMPROVED: Handle small background values more safely
      if (abs(vb) > small_ref) then
         rel_inc = abs(inc) / abs(vb)
         if (rel_inc > max_rel_inc_tracer) then
            scale_factor = max_rel_inc_tracer / rel_inc
            va = vb + inc * scale_factor
            vbig = vbig + 1
         end if
      else if (abs(inc) > 0.1_dp * abs(vb) + small_ref) then
         ! For very small background, use absolute bound
         va = vb + sign(small_ref + 0.1_dp * abs(vb), inc)
         vbig = vbig + 1
      end if
   end select

   ! 3. Enforce hard physical grid constraints
   if (va > vmax .or. va < vmin) then
      va = vb
      vout = vout + 1
   end if

end subroutine check_one_val
!=======================================================================

!=======================================================================
! Build mean and std of the ensemble (background or analysis)
! tflag = 'a' → analysis; otherwise background
!=======================================================================
subroutine make_mean_std(tflag)
   implicit none
   character(len=1), intent(in) :: tflag

   if (tflag == 'a') then
      call mean_state(Aan_m, Aan, nrens)
      call std_state(Aan_std, Aan, Aan_m, nrens)
   else
      call mean_state(Abk_m, Abk, nrens)
      call std_state(Abk_std, Abk, Abk_m, nrens)
   end if
end subroutine make_mean_std
!=======================================================================

!=======================================================================
subroutine allocate_all()
! Allocate all ensemble states and their mean/std containers
!=======================================================================
   implicit none
   integer :: ne

   if (.not. allocated(Abk)) allocate(Abk(nrens))
   if (.not. allocated(Aan)) allocate(Aan(nrens))

   call allocate_states(Abk_m, nnkn, nnel, nnlv)
   call allocate_states(Aan_m, nnkn, nnel, nnlv)
   call allocate_states(Abk_std, nnkn, nnel, nnlv)
   call allocate_states(Aan_std, nnkn, nnel, nnlv)

   do ne = 1, nrens
      call allocate_states(Abk(ne), nnkn, nnel, nnlv)
      call allocate_states(Aan(ne), nnkn, nnel, nnlv)
   end do

end subroutine allocate_all

!=======================================================================
! Write one state to a restart file (double precision → single precision)
!=======================================================================
subroutine write_state(Astate, filename)
   use mod_hydro
   use mod_ts
   use mod_restart
   implicit none

   type(states), intent(in) :: Astate
   character(len=*), intent(in) :: filename
   type(states4), allocatable :: A4

   allocate(A4)
   call allocate_states4(A4, nnkn, nnel, nnlv)

   ! Convert double→single precision container
   call states8to4(A4, Astate)

   ! Move data into hydrodynamic module variables
   call pull_state(A4)

   ! Write restart
   call rst_write(trim(filename), atime_an)

   deallocate(A4)
end subroutine write_state
!=======================================================================

!=======================================================================
! Read one state from restart file (single precision → double precision)
!=======================================================================
subroutine read_state(Astate, filename)
   use mod_restart
   implicit none

   type(states), intent(inout) :: Astate
   character(len=*), intent(in) :: filename
   type(states4), allocatable :: A4

   allocate(A4)
   call allocate_states4(A4, nnkn, nnel, nnlv)

   ! Read restart into SHYFEM global variables
   call rst_read(filename, atime_an)

   ! Transfer SHYFEM fields into A4
   call push_state(A4)

   ! Convert single→double and store into Astate
   call states4to8(Astate, A4)

   deallocate(A4)
end subroutine read_state
!=======================================================================

!=======================================================================
subroutine push_state(A4)
   use mod_hydro
   use mod_hydro_vel
   use mod_ts
   use mod_conz
   use mod_restart
   implicit none

   type(states4), intent(inout) :: A4
   integer :: nbc, i, k
   integer, allocatable :: bcid(:)
   real(dp), allocatable :: bcrho(:)
   logical :: file_exists

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

   !-----------------------------
   ! (2) Copy values to A4 (single precision container)
   !-----------------------------
   A4%u = utlnv
   A4%v = vtlnv
   A4%z = znv

   if (ibarcl_rst /= 0) then
      A4%t = tempv
      A4%s = saltv
   else
      A4%t = 0.0
      A4%s = 0.0
   end if

end subroutine push_state
!=======================================================================

!=======================================================================
! pull_state(A4)
!
!  Moves values from the A4 container into SHYFEM global arrays.
!  This is the symmetric counterpart of push_state.
!=======================================================================
subroutine pull_state(A4)
   use mod_hydro
   use mod_hydro_vel
   use mod_ts
   use mod_conz
   use mod_restart
   implicit none

   type(states4), intent(in) :: A4

   utlnv = A4%u
   vtlnv = A4%v
   znv   = A4%z

   if (ibarcl_rst /= 0) then
      tempv = A4%t
      saltv = A4%s
   end if

end subroutine pull_state
!=======================================================================

!=======================================================================
! Conversion routines between ensemble of states and matrices, adding
! parameter estimation support
!=======================================================================
subroutine tystate_to_matrix(ibrcl, nens, ndim, A, ndim_pe, pe_mat, Amat)
   implicit none
   integer, intent(in) :: ibrcl, nens, ndim, ndim_pe
   type(states), intent(in) :: A(nens)
   real(dp), intent(in), optional :: pe_mat(ndim_pe, nens)
   real(dp), intent(out) :: Amat(ndim, nens)
   integer :: i, dimuv, dimts, dimz
   integer :: offset

   dimz  = nnkn
   dimuv = nnlv * nnel
   dimts = nnlv * nnkn

   ! Sanity check
   if (2*dimuv + dimz > ndim) then
      write(*,'(a,i10,a,i10)') 'ERROR: tystate_to_matrix dimension mismatch: ', &
                               2*dimuv+dimz, ' > ndim=', ndim
      error stop
   end if

   ! 1. Copy base hydrodynamic components (U, V, Z) for all ensembles
   do i = 1, nens         
      Amat(1:dimuv, i) = reshape(A(i)%u, [dimuv])        
      Amat(dimuv+1:2*dimuv, i) = reshape(A(i)%v, [dimuv])
      Amat(2*dimuv+1:2*dimuv+dimz, i) = A(i)%z    
   end do             
                      
   ! 2. Handle baroclinic components (T, S)
   if (ibrcl > 0) then
      if (2*dimuv+dimz+2*dimts > ndim) then
         write(*,'(a,i10,a,i10)') 'ERROR: tystate_to_matrix (baroclinic) dimension mismatch: ', &
                                  2*dimuv+dimz+2*dimts, ' > ndim=', ndim
         error stop
      end if
      
      do i = 1, nens
         Amat(2*dimuv+dimz+1 : 2*dimuv+dimz+dimts, i) = reshape(A(i)%t, [dimts])
         Amat(2*dimuv+dimz+dimts+1 : 2*dimuv+dimz+2*dimts, i) = reshape(A(i)%s, [dimts])
      end do
      
      offset = 2*dimuv + dimz + 2*dimts
   else
      offset = 2*dimuv + dimz
   end if

   ! 3. Append PE data if present
   if (ndim_pe > 0 .and. present(pe_mat)) then
      if (offset + ndim_pe > ndim) then
         write(*,'(a,i10,a,i10)') 'ERROR: tystate_to_matrix (PE) dimension mismatch: ', &
                                  offset+ndim_pe, ' > ndim=', ndim
         error stop
      end if
      
      do i = 1, nens
         Amat(offset+1 : offset+ndim_pe, i) = pe_mat(:, i)
      end do
   end if

end subroutine tystate_to_matrix

!=======================================================================
subroutine matrix_to_tystate(ibrcl, nens, ndim, Amat, ndim_pe, pe_mat, A)
   implicit none
   integer, intent(in) :: ibrcl, nens, ndim, ndim_pe
   real(dp), intent(in) :: Amat(ndim, nens)
   real(dp), intent(out), optional :: pe_mat(ndim_pe, nens)
   type(states), intent(inout) :: A(nens)
   integer :: i, dimuv, dimts, dimz
   integer :: offset

   dimz  = nnkn
   dimuv = nnlv * nnel
   dimts = nnlv * nnkn

   ! Sanity check
   if (2*dimuv + dimz > ndim) then
      write(*,'(a,i10,a,i10)') 'ERROR: matrix_to_tystate dimension mismatch: ', &
                               2*dimuv+dimz, ' > ndim=', ndim
      error stop
   end if

   ! 1. Extract base hydrodynamic components (U, V, Z) for all ensembles
   do i = 1, nens
      A(i)%u = reshape(Amat(1:dimuv, i), [nnlv, nnel])
      A(i)%v = reshape(Amat(dimuv+1:2*dimuv, i), [nnlv, nnel])
      A(i)%z = Amat(2*dimuv+1:2*dimuv+dimz, i)
   end do

   ! 2. Handle baroclinic components (T, S)
   if (ibrcl > 0) then
      if (2*dimuv+dimz+2*dimts > ndim) then
         write(*,'(a,i10,a,i10)') 'ERROR: matrix_to_tystate (baroclinic) dimension mismatch: ', &
                                  2*dimuv+dimz+2*dimts, ' > ndim=', ndim
         error stop
      end if
      
      do i = 1, nens
         A(i)%t = reshape(Amat(2*dimuv+dimz+1 : 2*dimuv+dimz+dimts, i), [nnlv, nnkn])
         A(i)%s = reshape(Amat(2*dimuv+dimz+dimts+1 : 2*dimuv+dimz+2*dimts, i), [nnlv, nnkn])
      end do
      
      offset = 2*dimuv + dimz + 2*dimts
   else
      offset = 2*dimuv + dimz
   end if

   ! 3. Extract PE data if present
   if (ndim_pe > 0 .and. present(pe_mat)) then
      if (offset + ndim_pe > ndim) then
         write(*,'(a,i10,a,i10)') 'ERROR: matrix_to_tystate (PE) dimension mismatch: ', &
                                  offset+ndim_pe, ' > ndim=', ndim
         error stop
      end if
      
      do i = 1, nens
         pe_mat(:, i) = Amat(offset+1 : offset+ndim_pe, i)
      end do
   end if

end subroutine matrix_to_tystate

!=======================================================================
! Convert ensemble of qstates to matrix (2*ndim x nens)
! The first block holds q-fields, the second block holds state fields.
!=======================================================================
subroutine tyqstate_to_matrix(ibrcl, nens, ndim, A, Amat)
   implicit none
   integer,  intent(in)  :: ibrcl, nens, ndim
   type(qstates), intent(in) :: A(nens)
   real(dp), intent(out) :: Amat(2*ndim, nens)

   integer :: i, dimuv, dimts, dimz

   dimz  = nnkn
   dimuv = nnlv * nnel
   dimts = nnlv * nnkn

   ! Sanity check
   if (2*dimuv + dimz > ndim) then
      write(*,'(a,i10,a,i10)') 'ERROR: tyqstate_to_matrix dimension mismatch: ', &
                               2*dimuv+dimz, ' > ndim=', ndim
      error stop
   end if

   ! --- q fields (first block: 1 to ndim)
   do i = 1, nens
      Amat(1:dimuv, i) = reshape(A(i)%qu, [dimuv])
      Amat(dimuv+1:2*dimuv, i) = reshape(A(i)%qv, [dimuv])
      Amat(2*dimuv+1:2*dimuv+dimz, i) = A(i)%qz
   end do

   if (ibrcl > 0) then
      do i = 1, nens
         Amat(2*dimuv+dimz+1 : 2*dimuv+dimz+dimts, i) = reshape(A(i)%qt, [dimts])
         Amat(2*dimuv+dimz+dimts+1 : 2*dimuv+dimz+2*dimts, i) = reshape(A(i)%qs, [dimts])
      end do
   end if

   ! --- state fields (second block: ndim+1 to 2*ndim)
   do i = 1, nens
      Amat(ndim+1:ndim+dimuv, i) = reshape(A(i)%u, [dimuv])
      Amat(ndim+dimuv+1:ndim+2*dimuv, i) = reshape(A(i)%v, [dimuv])
      Amat(ndim+2*dimuv+1:ndim+2*dimuv+dimz, i) = A(i)%z
   end do

   if (ibrcl > 0) then
      do i = 1, nens
         Amat(ndim+2*dimuv+dimz+1 : ndim+2*dimuv+dimz+dimts, i) = reshape(A(i)%t, [dimts])
         Amat(ndim+2*dimuv+dimz+dimts+1 : ndim+2*dimuv+dimz+2*dimts, i) = reshape(A(i)%s, [dimts])
      end do
   end if
end subroutine tyqstate_to_matrix

!=======================================================================
! Convert matrix (2*ndim x nens) to ensemble of qstates
!=======================================================================
subroutine matrix_to_tyqstate(ibrcl, nens, ndim, Amat, A)
   implicit none
   integer,  intent(in)  :: ibrcl, nens, ndim
   real(dp), intent(in)  :: Amat(2*ndim, nens)
   type(qstates), intent(out) :: A(nens)

   integer :: i, dimuv, dimts, dimz

   dimz  = nnkn
   dimuv = nnlv * nnel
   dimts = nnlv * nnkn

   ! Sanity check
   if (2*dimuv + dimz > ndim) then
      write(*,'(a,i10,a,i10)') 'ERROR: matrix_to_tyqstate dimension mismatch: ', &
                               2*dimuv+dimz, ' > ndim=', ndim
      error stop
   end if

   ! --- q fields
   do i = 1, nens
      A(i)%qu = reshape(Amat(1:dimuv, i), [nnlv, nnel])
      A(i)%qv = reshape(Amat(dimuv+1:2*dimuv, i), [nnlv, nnel])
      A(i)%qz = Amat(2*dimuv+1:2*dimuv+dimz, i)
   end do

   if (ibrcl > 0) then
      do i = 1, nens
         A(i)%qt = reshape(Amat(2*dimuv+dimz+1 : 2*dimuv+dimz+dimts, i), [nnlv, nnkn])
         A(i)%qs = reshape(Amat(2*dimuv+dimz+dimts+1 : 2*dimuv+dimz+2*dimts, i), [nnlv, nnkn])
      end do
   end if

   ! --- state fields
   do i = 1, nens
      A(i)%u = reshape(Amat(ndim+1:ndim+dimuv, i), [nnlv, nnel])
      A(i)%v = reshape(Amat(ndim+dimuv+1:ndim+2*dimuv, i), [nnlv, nnel])
      A(i)%z = Amat(ndim+2*dimuv+1:ndim+2*dimuv+dimz, i)
   end do

   if (ibrcl > 0) then
      do i = 1, nens
         A(i)%t = reshape(Amat(ndim+2*dimuv+dimz+1 : ndim+2*dimuv+dimz+dimts, i), [nnlv, nnkn])
         A(i)%s = reshape(Amat(ndim+2*dimuv+dimz+dimts+1 : ndim+2*dimuv+dimz+2*dimts, i), [nnlv, nnkn])
      end do
   end if
end subroutine matrix_to_tyqstate

end module mod_ens_state
