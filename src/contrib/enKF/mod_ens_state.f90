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

   ! Interface to external single-precision function ipint (from subnsu.f)
   !interface
   !   integer function ipint(i)
   !      integer, intent(in) :: i
   !   end function ipint
   !end interface

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
   if (( nnkn /= nkn ).or.( nnel /= nel )) error stop 'Horizontal dimensions of restart and basin differ.'

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

   if (bnew_ens == 0 .or. nanal > 1) then
      write(*,*) 'Loading existing ensemble...'
      do ne = 1, nrens
         call num2str(ne-1, nrel)
         rstname = 'an'//nal//'_'//'en'//nrel//'b.rst'
	 inquire(file=rstname,exist=bexist)
	 if( .not. bexist ) then
	   write(6,*) 'restart file does not exists: ',trim(rstname)
	   error stop 'error stop read_ensemble: no restart file'
	 end if
         call read_state(Abk(ne), rstname)
      end do

   else if (bnew_ens == 1 .and. nanal == 1) then
      write(*,*) 'Creating a new ensemble...'
      call num2str(0, nrel)
      rstname = 'an'//nal//'_en'//nrel//'b.rst'
      call read_state(Abk(1), rstname)
      call make_init_ens(Abk(1))

      do ne = 1, nrens
         call num2str(ne-1, nrel)
         rstname='an'//nal//'_'//'en'//nrel//'b.rst'
         call write_state(Abk(ne), rstname)
      end do

   else
      write(*,*) 'Invalid bnew_ens option'
      error stop
   end if

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

!=======================================================================
! Create the initial ensemble using random 2D perturbations
!=======================================================================
subroutine make_init_ens(Ain)
   use mod_restart
   use mod_para
   implicit none

   type(states), intent(in) :: Ain
   real(dp) :: kvec1(nnkn, nrens-1), kvec2(nnkn, nrens-1)
   integer :: ne, nl

   ! Generate perturbations for the free-surface level
   call make_2Dpert(kvec1, nnkn, nrens-1)

   Abk(1) = Ain

   do ne = 2, nrens
      Abk(ne) = Ain
      Abk(ne)%z = Ain%z + kvec1(:, ne-1) * sigma_init_z
   end do

   ! Temperature–salinity perturbation if baroclinic mode active
   if (ibarcl_rst /= 0) then
      call make_2Dpert(kvec1, nnkn, nrens-1)
      call make_2Dpert(kvec2, nnkn, nrens-1)

      do ne = 2, nrens
         do nl = 1, nnlv
            Abk(ne)%t(nl,:) = Ain%t(nl,:) + kvec1(:, ne-1) * sigma_init_t
            Abk(ne)%s(nl,:) = Ain%s(nl,:) + kvec2(:, ne-1) * sigma_init_s
         end do
      end do
   end if
end subroutine make_init_ens
!=======================================================================



!=======================================================================
! Check and correct ensemble fields near boundaries and ensure values
! remain physical (no NaN, no large increments, no out-of-range).
!
! This routine is OpenMP-parallel and now thread-safe.
!=======================================================================
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
      write(*,*) 'Warning: lbound.dat not found - no boundary correction applied.'
   end if

   otot = 0 ; ntot = 0 ; btot = 0

   !===================================================================
   ! Loop over ensemble members
   !===================================================================
   do ne = 1, nrens

      zout = 0 ; uvout = 0 ; sout = 0 ; tout = 0
      znan = 0 ; uvnan = 0 ; snan = 0 ; tnan = 0
      zbig = 0 ; uvbig = 0 ; sbig = 0 ; tbig = 0

      !===============================================================
      ! NODE-BASED FIELDS: z, T, S
      !===============================================================
!$OMP PARALLEL DO PRIVATE(k,nl,w)                  &
!$OMP PRIVATE(znan,zout,snan,sout,tnan,tout,zbig,tbig,sbig,uvbig) &
!$OMP SHARED(ne,Abk,Aan,file_exists,nbc,bcid,bcrho)
      do k = 1, nnkn

         ! Boundary weight w(k)
         if (file_exists) then
            call bc_correction('node', k, nbc, bcid, bcrho, w)
         else
            w = 0.0_dp
         end if

         !-------------------------
         ! Free-surface elevation z
         !-------------------------
         call check_one_val( &
            va = Aan(ne)%z(k), &
            vb = Abk(ne)%z(k), &
            vmax = SSH_MAX, &
            vmin = SSH_MIN, &
            vnan = znan, vout = zout, vbig = zbig )

         Aan(ne)%z(k) = w * Abk(ne)%z(k) + (1.0_dp - w) * Aan(ne)%z(k)

         !-------------------------
         ! Temperature & salinity
         !-------------------------
         do nl = 1, nnlv

            call check_one_val( &
               va = Aan(ne)%s(nl,k), &
               vb = Abk(ne)%s(nl,k), &
               vmax = SAL_MAX, vmin = SAL_MIN, &
               vnan = snan, vout = sout, vbig = sbig)

            Aan(ne)%s(nl,k) = w * Abk(ne)%s(nl,k) + (1.0_dp - w) * Aan(ne)%s(nl,k)

            call check_one_val( &
               va = Aan(ne)%t(nl,k), &
               vb = Abk(ne)%t(nl,k), &
               vmax = TEM_MAX, vmin = TEM_MIN, &
               vnan = tnan, vout = tout, vbig = tbig)

            Aan(ne)%t(nl,k) = w * Abk(ne)%t(nl,k) + (1.0_dp - w) * Aan(ne)%t(nl,k)
         end do

      end do
!$OMP END PARALLEL DO


      !===============================================================
      ! ELEMENT-BASED FIELDS: u, v
      !===============================================================
!$OMP PARALLEL DO PRIVATE(ie,nl,w) SHARED(ne,Abk,Aan,file_exists,nbc,bcid,bcrho)
      do ie = 1, nnel

         if (file_exists) then
            call bc_correction('elem', ie, nbc, bcid, bcrho, w)
         else
            w = 0.0_dp
         end if

         do nl = 1, nnlv
            Aan(ne)%u(nl,ie) = w * Abk(ne)%u(nl,ie) + (1.0_dp - w) * Aan(ne)%u(nl,ie)
            Aan(ne)%v(nl,ie) = w * Abk(ne)%v(nl,ie) + (1.0_dp - w) * Aan(ne)%v(nl,ie)
         end do
      end do
!$OMP END PARALLEL DO

      if (verbose) then
         if (zout>0 .or. sout>0 .or. tout>0 .or. uvout>0) then
            write(*,*) 'Member', ne, ' out of range counts = ', zout, sout, tout, uvout
         end if
      end if

      otot = otot + zout + sout + tout + uvout
      ntot = ntot + znan + snan + tnan + uvnan
      btot = btot + zbig + sbig + tbig + uvbig

   end do

   if (otot > 0) write(*,*) 'Total values out of range: ', otot
   if (ntot > 0) write(*,*) 'Total NaN values: ', ntot
   if (btot > 0) write(*,*) 'Total excessive increments: ', btot

end subroutine bc_val_check_correct
!=======================================================================



!=======================================================================
! Read boundary condition file (two-stage read)
! icall==0 → read number of BC nodes
! icall==1 → read list of BC nodes and radii
!=======================================================================
subroutine read_bc_file(icall, bcfile, nbc, bcid, bcrho)
   implicit none
   character(len=*), intent(in) :: bcfile
   integer, intent(in) :: icall
   integer, intent(inout) :: nbc
   integer, intent(out) :: bcid(nbc)
   real(dp), intent(out) :: bcrho(nbc)

   integer :: i

   if (icall == 0) then
      open(28, file=trim(bcfile), status='old')
      read(28,*) nbc
      close(28)
      return
   end if

   open(28, file=trim(bcfile), status='old')
   read(28,*) nbc   ! skip first line
   do i = 1, nbc
      read(28,*) bcid(i), bcrho(i)
   end do
   close(28)

end subroutine read_bc_file
!=======================================================================



!=======================================================================
! Compute weight for boundary damping based on distance from BC node.
!=======================================================================
subroutine bc_correction(stype, id, nbc, bcid, bcrho, w)
   implicit none

   character(len=*), intent(in) :: stype
   integer, intent(in) :: id
   integer, intent(in) :: nbc
   integer, intent(in) :: bcid(nbc)
   real(dp), intent(in) :: bcrho(nbc)
   real(dp), intent(out) :: w

   integer*4 :: kext
   integer :: i, k, kbc
   real(dp) :: x, y, bcx, bcy
   real(dp) :: d, dmin, rho

   integer*4 ipint

   x = 0.0_dp ; y = 0.0_dp
   dmin = 1.0e15_dp
   rho = 0.0_dp

   if (stype == 'node') then
      x = xgv(id)
      y = ygv(id)
   else
      ! element → average vertices
      do i = 1, 3
         k = nen3v(i, id)
         x = x + xgv(k)
         y = y + ygv(k)
      end do
      x = x / 3.0_dp
      y = y / 3.0_dp
   end if

   do i = 1, nbc
      kext = bcid(i)
      kbc = ipint(kext)
      bcx = xgv(kbc)
      bcy = ygv(kbc)
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
! Check a single scalar:
! - replace NaN with background
! - limit increment to 80%
! - ensure value stays within physical bounds
!=======================================================================
subroutine check_one_val(va, vb, vmax, vmin, vnan, vout, vbig)
   implicit none
   real(dp), intent(inout) :: va
   real(dp), intent(in) :: vb, vmin, vmax
   integer, intent(inout) :: vnan, vout, vbig
   real(dp), parameter :: max_inc = 0.8_dp
   real(dp) :: inc, rel_inc

   ! Replace NaN
   if (isnan(va) .and. .not. isnan(vb)) then
      va = vb
!$OMP ATOMIC
      vnan = vnan + 1
   end if

   ! Limit increment
   inc = va - vb
   if (abs(vb) > 0.0_dp) then
      rel_inc = abs(inc) / abs(vb)
   else
      rel_inc = 0.0_dp
   end if

   if (rel_inc > max_inc) then
      va = vb + inc * (max_inc / rel_inc)
!$OMP ATOMIC
      vbig = vbig + 1
   end if

   ! Enforce bounds
   if (va >= vmax .or. va <= vmin) then
      va = vb
!$OMP ATOMIC
      vout = vout + 1
   end if

end subroutine check_one_val
!=======================================================================

!=======================================================================
! Build mean and std of the ensemble (background or analysis)
! tflag = 'a' → analysis; background
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
   integer :: nbc, i, k, k_int
   integer, allocatable :: bcid(:)
   real(dp), allocatable :: bcrho(:)
   logical :: file_exists
   real(dp) :: zmean, zstd, dist
   integer, parameter :: fact = 8   ! z-threshold factor

   !-----------------------------
   ! (1) Check that SHYFEM arrays are allocated
   !-----------------------------
   if (.not. allocated(znv)) stop 'ERROR: znv not allocated in push_state'
   if (.not. allocated(utlnv)) stop 'ERROR: utlnv not allocated'
   if (.not. allocated(vtlnv)) stop 'ERROR: vtlnv not allocated'
   if (ibarcl_rst /= 0) then
      if (.not. allocated(tempv)) stop 'ERROR: tempv not allocated'
      if (.not. allocated(saltv)) stop 'ERROR: saltv not allocated'
   end if

   !-----------------------------
   ! (2) Boundary file reading
   !-----------------------------
   nbc = 1
   allocate(bcid(nbc), bcrho(nbc))
   inquire(file='lbound.dat', exist=file_exists)

   if (file_exists) then
      call read_bc_file(0, 'lbound.dat', nbc, bcid, bcrho)
      deallocate(bcid, bcrho)
      allocate(bcid(nbc), bcrho(nbc))
      call read_bc_file(1, 'lbound.dat', nbc, bcid, bcrho)
   else
      bcid = 1
      bcrho = 0.0_dp
   end if

   !-----------------------------
   ! (3) Validate free-surface elevation znv
   !-----------------------------
!   zmean = sum(znv) / nnkn
!   zstd  = sqrt( sum(znv*znv)/nnkn - zmean*zmean )
!
!   do k = 1, nnkn
!      do i = 1, nbc
!         k_int = ipint(bcid(i))   ! now safe due to interface
!         dist = sqrt( (xgv(k)-xgv(k_int))**2 + (ygv(k)-ygv(k_int))**2 )
!
!         if (dist > 1.5_dp * bcrho(i)) then
!            if (znv(k) > zmean + fact*zstd) then
!               write(*,*) 'Warning: large z-value at node ', k
!               znv(k) = zmean + fact*zstd
!            else if (znv(k) < zmean - fact*zstd) then
!               write(*,*) 'Warning: low z-value at node ', k
!               znv(k) = zmean - fact*zstd
!            end if
!         end if
!      end do
!   end do

   !-----------------------------
   ! (4) Copy values to A4
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
   else
      tempv = 0.0
      saltv = 0.0
   end if

   !call layer_thick(nnel)
end subroutine pull_state
!=======================================================================



!=======================================================================
! Conversion routines between ensemble of states and matrices
!=======================================================================
subroutine tystate_to_matrix(ibrcl, nens, ndim, A, Amat)
   implicit none
   integer, intent(in) :: ibrcl, nens, ndim
   type(states), intent(in) :: A(nens)
   real(dp), intent(out) :: Amat(ndim, nens)
   integer :: i, dimuv, dimts, dimz

   dimz  = nnkn
   dimuv = nnlv * nnel
   dimts = nnlv * nnkn

   do i = 1, nens
      Amat(1:dimuv, i) = reshape(A(i)%u, (/dimuv/))
      Amat(dimuv+1:2*dimuv, i) = reshape(A(i)%v, (/dimuv/))
      Amat(2*dimuv+1:2*dimuv+dimz, i) = A(i)%z
   end do

   if (ibrcl > 0) then
      do i = 1, nens
         Amat(2*dimuv+dimz+1 : 2*dimuv+dimz+dimts, i) = reshape(A(i)%t, (/dimts/))
         Amat(2*dimuv+dimz+dimts+1 : 2*dimuv+dimz+2*dimts, i) = reshape(A(i)%s, (/dimts/))
      end do
   end if
end subroutine tystate_to_matrix
!=======================================================================



!=======================================================================
subroutine matrix_to_tystate(ibrcl, nens, ndim, Amat, A)
   implicit none
   integer, intent(in) :: ibrcl, nens, ndim
   real(dp), intent(in) :: Amat(ndim, nens)
   type(states), intent(inout) :: A(nens)
   integer :: i, dimuv, dimts, dimz

   dimz  = nnkn
   dimuv = nnlv * nnel
   dimts = nnlv * nnkn

   do i = 1, nens
      A(i)%u = reshape(Amat(1:dimuv, i), (/nnlv, nnel/))
      A(i)%v = reshape(Amat(dimuv+1:2*dimuv, i), (/nnlv, nnel/))
      A(i)%z = Amat(2*dimuv+1:2*dimuv+dimz, i)
   end do

   if (ibrcl > 0) then
      do i = 1, nens
         A(i)%t = reshape(Amat(2*dimuv+dimz+1 : 2*dimuv+dimz+dimts, i), (/nnlv, nnkn/))
         A(i)%s = reshape(Amat(2*dimuv+dimz+dimts+1 : 2*dimuv+dimz+2*dimts, i), (/nnlv, nnkn/))
      end do
   end if
end subroutine matrix_to_tystate
!=======================================================================

!=======================================================================
! Convert ensemble of qstates → matrix (2*ndim x nens)
! The first block holds q-fields, the second block holds state fields.
!=======================================================================
subroutine tyqstate_to_matrix(ibrcl, nens, ndim, A, Amat)
   implicit none
   integer,  intent(in)  :: ibrcl
   integer,  intent(in)  :: nens
   integer,  intent(in)  :: ndim
   type(qstates), intent(in) :: A(nens)
   real(dp), intent(out) :: Amat(2*ndim, nens)

   integer :: i, dimuv, dimts, dimz

   dimz  = nnkn
   dimuv = nnlv * nnel
   dimts = nnlv * nnkn

   ! --- q fields
   do i = 1, nens
      Amat(1:dimuv, i) = reshape(A(i)%qu, (/dimuv/))
      Amat(dimuv+1:2*dimuv, i) = reshape(A(i)%qv, (/dimuv/))
      Amat(2*dimuv+1:2*dimuv+dimz, i) = A(i)%qz
   end do

   if (ibrcl > 0) then
      do i = 1, nens
         Amat(2*dimuv+dimz+1 : 2*dimuv+dimz+dimts, i)   = reshape(A(i)%qt, (/dimts/))
         Amat(2*dimuv+dimz+dimts+1 : 2*dimuv+dimz+2*dimts, i) = reshape(A(i)%qs, (/dimts/))
      end do
   end if

   ! --- state fields (after the first ndim block)
   do i = 1, nens
      Amat(ndim+1:ndim+dimuv, i) = reshape(A(i)%u, (/dimuv/))
      Amat(ndim+dimuv+1:ndim+2*dimuv, i) = reshape(A(i)%v, (/dimuv/))
      Amat(ndim+2*dimuv+1:ndim+2*dimuv+dimz, i) = A(i)%z
   end do

   if (ibrcl > 0) then
      do i = 1, nens
         Amat(ndim+2*dimuv+dimz+1 : ndim+2*dimuv+dimz+dimts, i) = reshape(A(i)%t, (/dimts/))
         Amat(ndim+2*dimuv+dimz+dimts+1 : ndim+2*dimuv+dimz+2*dimts, i) = reshape(A(i)%s, (/dimts/))
      end do
   end if
end subroutine tyqstate_to_matrix
!=======================================================================

!=======================================================================
! Convert matrix (2*ndim x nens) → ensemble of qstates
!=======================================================================
subroutine matrix_to_tyqstate(ibrcl, nens, ndim, Amat, A)
   implicit none
   integer,  intent(in)  :: ibrcl
   integer,  intent(in)  :: nens
   integer,  intent(in)  :: ndim
   real(dp), intent(in)  :: Amat(2*ndim, nens)
   type(qstates), intent(out) :: A(nens)

   integer :: i, dimuv, dimts, dimz

   dimz  = nnkn
   dimuv = nnlv * nnel
   dimts = nnlv * nnkn

   ! --- q fields
   do i = 1, nens
      A(i)%qu = reshape(Amat(1:dimuv, i), (/nnlv, nnel/))
      A(i)%qv = reshape(Amat(dimuv+1:2*dimuv, i), (/nnlv, nnel/))
      A(i)%qz = Amat(2*dimuv+1:2*dimuv+dimz, i)
   end do

   if (ibrcl > 0) then
      do i = 1, nens
         A(i)%qt = reshape(Amat(2*dimuv+dimz+1 : 2*dimuv+dimz+dimts, i), (/nnlv, nnkn/))
         A(i)%qs = reshape(Amat(2*dimuv+dimz+dimts+1 : 2*dimuv+dimz+2*dimts, i), (/nnlv, nnkn/))
      end do
   end if

   ! --- state fields
   do i = 1, nens
      A(i)%u = reshape(Amat(ndim+1:ndim+dimuv, i), (/nnlv, nnel/))
      A(i)%v = reshape(Amat(ndim+dimuv+1:ndim+2*dimuv, i), (/nnlv, nnel/))
      A(i)%z = Amat(ndim+2*dimuv+1:ndim+2*dimuv+dimz, i)
   end do

   if (ibrcl > 0) then
      do i = 1, nens
         A(i)%t = reshape(Amat(ndim+2*dimuv+dimz+1 : ndim+2*dimuv+dimz+dimts, i), (/nnlv, nnkn/))
         A(i)%s = reshape(Amat(ndim+2*dimuv+dimz+dimts+1 : ndim+2*dimuv+dimz+2*dimts, i), (/nnlv, nnkn/))
      end do
   end if
end subroutine matrix_to_tyqstate
!=======================================================================

end module mod_ens_state
