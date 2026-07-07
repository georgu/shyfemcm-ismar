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
   use mod_layer_thickness
   use sigma
   use mod_area
   use evgeom
   use mod_depth
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
   call mod_layer_thickness_init(nnkn, nnel, nnlv)
   call init_sigma_info(nnlv,hlv)
   call mod_area_init(nnkn,nnlv)
   call ev_init(nnel)
   call mod_depth_init(nnkn,nnel)

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
! SUBROUTINE: val_check_correct
!
! PURPOSE:
!   Scans the analysed ensemble states (Aan) against background fields (Abk)
!   to enforce physical thresholds for Sea Surface Height (SSH), Salinity (S),
!   Temperature (T), and Velocity Vectors (u, v). Replaces anomalous values 
!   with prior background coefficients upon violation.
!=======================================================================
subroutine val_check_correct()
   implicit none

   ! Legitimate local loop counters and grid pointers
   integer :: ne, k, nl, ie
   integer :: zout, sout, tout, uout, vout
   integer :: znan, snan, tnan, unan, vnan
   integer :: zbig, sbig, tbig, ubig, vbig
   integer :: otot, ntot, btot

   ! Initialize global diagnostics accumulators
   otot = 0 ; ntot = 0 ; btot = 0

   !===================================================================
   ! Loop over ensemble members (Executed sequentially for data safety)
   !===================================================================
   do ne = 1, nrens

      ! Initialize member-specific quality control counters
      zout = 0 ; sout = 0 ; tout = 0; uout = 0; vout = 0
      znan = 0 ; snan = 0 ; tnan = 0; unan = 0; vnan = 0
      zbig = 0 ; sbig = 0 ; tbig = 0; ubig = 0; vbig = 0

      !===============================================================
      ! NODE-BASED FIELDS: Sea Surface Height (z), Temperature (T), Salinity (S)
      !===============================================================
      do k = 1, nnkn

         ! Validate physical thresholds and metric increments for SSH (Z field)
         call check_one_val( &
            vtype = 'Z', &
            va = Aan(ne)%z(k), &
            vb = Abk(ne)%z(k), &
            vmax = SSH_MAX, &
            vmin = SSH_MIN, &
            vnan = znan, vout = zout, vbig = zbig )

         !------------------------------------------------------------
         ! 3D Hydrographic Fields: Temperature & Salinity Profiles
         !------------------------------------------------------------
         do nl = 1, nnlv

            ! Validate Salinity profiles (S field)
            call check_one_val( &
               vtype = 'S', &
               va = Aan(ne)%s(nl,k), &
               vb = Abk(ne)%s(nl,k), &
               vmax = SAL_MAX, vmin = SAL_MIN, &
               vnan = snan, vout = sout, vbig = sbig)

            ! Validate Temperature profiles (T field)
            call check_one_val( &
               vtype = 'T', &
               va = Aan(ne)%t(nl,k), &
               vb = Abk(ne)%t(nl,k), &
               vmax = TEM_MAX, vmin = TEM_MIN, &
               vnan = tnan, vout = tout, vbig = tbig)
         end do

      end do

      !===============================================================
      ! ELEMENT-BASED FIELDS: Velocity Components (u, v Vectors)
      !===============================================================
      do ie = 1, nnel
         do nl = 1, nnlv
            
            ! Validate horizontal zonal velocity component (u field)
            call check_one_val( &
               vtype = 'V', &
               va = Aan(ne)%u(nl,ie), &
               vb = Abk(ne)%u(nl,ie), &
               vmax = VEL_MAX, vmin = VEL_MIN, &
               vnan = unan, vout = uout, vbig = ubig)

            ! Validate horizontal meridional velocity component (v field)
            ! MATHEMATICAL FIX: Integrated isolated 'v' counters to track components reliably
            call check_one_val( &
               vtype = 'V', &
               va = Aan(ne)%v(nl,ie), &
               vb = Abk(ne)%v(nl,ie), &
               vmax = VEL_MAX, vmin = VEL_MIN, &
               vnan = vnan, vout = vout, vbig = vbig)
         end do
      end do

      ! Concrete tracking reduction accumulation across active members
      otot = otot + zout + sout + tout + uout + vout
      ntot = ntot + znan + snan + tnan + unan + vnan
      btot = btot + zbig + sbig + tbig + ubig + vbig

   end do

   !===================================================================
   ! SERIAL DIAGNOSTICS LOGGING
   !===================================================================
   if (otot > 0) write(*,'(A,I8)') 'Total out-of-range corrections per member: ', otot / nrens
   if (ntot > 0) write(*,'(A,I8)') 'Total NaN corrections per member:          ', ntot / nrens
   if (btot > 0) write(*,'(A,I8)') 'Total excessive increment corrections per member: ', btot / nrens

end subroutine val_check_correct

!=======================================================================
! SUBROUTINE: check_one_val
!
! PURPOSE:
!   Validates a single scalar ensemble trajectory entry against physical bounds 
!   and structural increments. Eliminates NaNs, limits runaway filter shocks, 
!   and forces hard grid boundary thresholds.
!=======================================================================
subroutine check_one_val(vtype, va, vb, vmax, vmin, vnan, vout, vbig)
   implicit none
   
   character(len=1), intent(in)    :: vtype
   real(dp),         intent(inout) :: va
   real(dp),         intent(in)    :: vb, vmin, vmax
   integer,          intent(inout) :: vnan, vout, vbig

   ! Tuning thresholds for numerical filters
   real(dp), parameter :: max_rel_inc_tracer = 0.8_dp  ! 80% maximum relative shift for tracers
   real(dp), parameter :: max_abs_tracer     = 2.0_dp  ! 2.0 maximum absolute shift for T/S near zero
   real(dp), parameter :: max_abs_ssh        = 0.4_dp  ! Max 40 cm increment for SSH to prevent gravity waves
   real(dp), parameter :: max_abs_uv         = 0.5_dp  ! Max 50 cm/s increment for currents (CFL safety)
   real(dp), parameter :: small_ref          = 1.0e-3_dp ! Reference threshold for small background targets

   real(dp) :: inc, rel_inc, scale_factor

   ! 1. Immediate removal of NaNs via non-equality validation property
   if (va /= va) then
      va = vb
      vnan = vnan + 1
      return
   end if

   ! 2. Dynamically bound structural increments using specific physics
   inc = va - vb

   select case (vtype)
   case ('Z')
      ! Sea Surface Height (SSH) absolute update clipping
      if (abs(inc) > max_abs_ssh) then
         va = vb + sign(max_abs_ssh, inc)
         vbig = vbig + 1
      end if

   case ('V')
      ! Velocity vectors absolute update clipping for CFL stability
      if (abs(inc) > max_abs_uv) then
         va = vb + sign(max_abs_uv, inc)
         vbig = vbig + 1
      end if

   case default
      ! Thermodynamics Scalars (Temperature & Salinity)
      if (abs(vb) > small_ref) then
         ! Standard proportional scaling for robust background fields
         rel_inc = abs(inc) / abs(vb)
         if (rel_inc > max_rel_inc_tracer) then
            scale_factor = max_rel_inc_tracer / rel_inc
            va = vb + inc * scale_factor
            vbig = vbig + 1
         end if
      else
         ! For low-background configurations, switch to a reasonable physical
         ! absolute bound instead of choking the filter updates near zero
         if (abs(inc) > max_abs_tracer) then
            va = vb + sign(max_abs_tracer, inc)
            vbig = vbig + 1
         end if
      end if
   end select

   ! 3. Enforce hard physical grid limits (Absolute validation check)
   if (va > vmax .or. va < vmin) then
      va = vb
      vout = vout + 1
   end if

end subroutine check_one_val

!=======================================================================
! Build mean and std of the ensemble (background or analysis)
! tflag = 'a' analysis; otherwise background
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
   use mod_layer_thickness
   use mod_hydro
   use mod_hydro_vel
   use levels
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
   ! (2) Compute the velocities from the transports
   !-----------------------------
   ulnv = 0.
   vlnv = 0.
   where( hdenv > 0. )
       ulnv = utlnv / hdenv
       vlnv = vtlnv / hdenv
   end where

   !-----------------------------
   ! (3) Copy values to A4 (single precision container)
   !-----------------------------
   A4%u = ulnv
   A4%v = vlnv
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
   use mod_layer_thickness
   use mod_hydro
   use mod_hydro_vel
   use levels
   use mod_ts
   use mod_conz
   use mod_restart
   implicit none

   type(states4), intent(in) :: A4

   ulnv = A4%u
   vlnv = A4%v
   znv   = A4%z

   if (ibarcl_rst /= 0) then
      tempv = A4%t
      saltv = A4%s
   end if

   !-----------------------------
   ! Compute the the transports
   !-----------------------------
   utlnv = ulnv * hdenv
   vtlnv = vlnv * hdenv
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
