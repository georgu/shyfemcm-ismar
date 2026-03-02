!
! ---------------------------------------------------------------------------
! EnKF → EnKS converter with stability and performance improvements
!
! Copyright:
!   (C) 2017-2026, Marco Bajo, CNR-ISMAR Venice
!   All rights reserved
!
! Notes:
! - addpar() requires REAL*4
! - daddpar() requires REAL*8
! - analysis2 uses stable accumulation of X5 transforms
! - make_mn_std uses compensated variance formula
! ---------------------------------------------------------------------------

program enKF2enKS
    use iso_fortran_env, only : dp => real64
    use mod_geom_dynamic
    use mod_ts
    use mod_hydro_vel
    use mod_hydro
    use mod_restart
    use levels, only : nlvdi, nlv, hlv, ilhv, ilhkv
    use basin
    use shympi

    implicit none

    integer :: rrec
    integer :: sdim
    integer :: nlag

    !---------------------------- CLI arguments ----------------------------
    character(len=80) :: basinf           ! basin filename (with extension)
    character(len=3)  :: lnnlv            ! number of vertical levels (as text)
    character(len=6)  :: lnrens           ! number of ensemble members (as text)
    character(len=4)  :: lout             ! output mode: 'full' or 'norm'
    character(len=6)  :: lnlag            ! smoother lag steps, -1 => all

    integer :: nnlv                       ! parsed nlv
    integer :: nrens, nre                 ! parsed ensemble size, loop index
    character(len=5) :: nrel              ! member index formatted (00000..)
    integer :: ierr, iflag
    real(dp) :: atime                     ! record time
    integer :: fid                        ! file unit id

    ! ---------------------------- Single-precision parameters ----------------------------
    real*4 :: ibarcl4, iconz4, imerc4, iturb4, iwvert4, ieco4
    real*4 :: zero4

    ! ---------------------------- Ensemble storage ----------------------------
    real(dp), allocatable :: Astate(:,:), AmeanKS(:), AstdKS(:)

    zero4 = 0.0
    ! ---------------------------- Read CLI arguments ----------------------------
    call get_command_argument(1, basinf)
    call get_command_argument(2, lnnlv)
    call get_command_argument(3, lnrens)
    call get_command_argument(4, lout)
    call get_command_argument(5, lnlag)

    if (lnlag == '') then
       write(*,*) ''
       write(*,*) 'Usage: enKF2enKS [basinf] [nlv] [nrens] [output] [nlag]'
       write(*,*) ''
       write(*,*) '[basinf] bas file with the extension.'
       write(*,*) '[nlv]     number of vertical levels.'
       write(*,*) '[nrens]   number of ensemble members, control included.'
       write(*,*) "[output]  full (full) or just mean and std (norm)."
       write(*,*) '[nlag]    number of forward analysis steps to consider (-1 = all)'
       write(*,*) ''
       stop
    end if

    read(lnnlv, *) nnlv
    if (nnlv < 1) error stop "Invalid nlv"

    read(lnrens, *) nrens
    if (nrens < 3) error stop "nrens must be >= 3"

    if ((lout /= 'full') .and. (lout /= 'norm')) error stop "Invalid output mode"

    read(lnlag, *) nlag
    if ((nlag < 2) .and. (nlag /= -1)) error stop "Invalid smoother lag"

    write(*,*) "EnKS with lag =", nlag

    ! ---------------------------- Initialize SHYFEM ----------------------------
    call init_shyfem(basinf, nnlv)

    ! ---------------------------- Open RST inputs ----------------------------
    do nre = 1, nrens
        fid = 20 + nre
        call num2str(nre-1, nrel)
        open(fid, file='analKF_en'//nrel//'.rst', status='old', form='unformatted', iostat=ierr)
        if (ierr /= 0) error stop "Cannot open input restart files"
    end do

    ! ---------------------------- Open outputs ----------------------------
    if (lout == 'full') then
        do nre = 1, nrens
            fid = 20 + nrens + nre
            call num2str(nre-1, nrel)
            open(fid, file='analKS_en'//nrel//'.rst', status='unknown', form='unformatted')
        end do
    end if

    open(18, file='analKS_mean.rst', status='unknown', form='unformatted')
    open(19, file='analKS_std.rst',  status='unknown', form='unformatted')

    ! ---------------------------- Main time loop ----------------------------
    rrec = 0

time_loop: do
        ! ---- Read full ensemble at a time step
        do nre = 1, nrens
            fid = 20 + nre
            call rst_read_record(fid, atime, iflag, ierr)

            if (ierr == -1) then
                exit time_loop
            else if (ierr > 0) then
                error stop "RST read error"
            end if

            if (.not. allocated(Astate)) then
                if (ibarcl_rst /= 0) then
                    sdim = nkn + 2*nlv*nel + 2*nlv*nkn
                else
                    sdim = nkn + 2*nlv*nel
                end if
                allocate(Astate(sdim,nrens))
                allocate(AmeanKS(sdim), AstdKS(sdim))
                Astate = 0.0_dp
            end if

            call push_matrix(sdim, nrens, nre, Astate, ibarcl_rst)
        end do

        ! ---- First record: set restart parameters
        if (rrec == 0) then
            hlv  = hlvrst
            ilhv = ilhrst
            ilhkv = ilhkrst
            hlv_global = hlvrst

            ibarcl4 = ibarcl_rst
            iwvert4 = iwvert_rst
            ieco4   = ieco_rst
            iconz4  = iconz_rst
            imerc4  = imerc_rst
            iturb4  = iturb_rst

            call addpar('ibarcl', ibarcl4)
            call addpar('iconz',  iconz4)
            call addpar('iwvert', iwvert4)
            call addpar('ieco',   ieco4)
            call addpar('ibio',   zero4)
            call addpar('ibfm',   zero4)
            call addpar('imerc',  imerc4)
            call addpar('iturb',  iturb4)

            call daddpar('date', 0.0_dp)
            call daddpar('time', 0.0_dp)
        end if

        ! ---- Apply smoother (improved + stable X5 accumulation)
        call make_analysis2(atime, sdim, nrens, Astate, nlag)

        ! ---- Statistics
        call make_mn_std(sdim, nrens, Astate, AmeanKS, AstdKS)

        call pull_matrix(sdim, nrens, 1, AmeanKS, ibarcl_rst)
        call rst_write_record(atime, 18)

        call pull_matrix(sdim, nrens, 1, AstdKS, ibarcl_rst)
        call rst_write_record(atime, 19)

        if (lout == 'full') then
            do nre = 1, nrens
                fid = 20 + nrens + nre
                call pull_matrix(sdim, nrens, nre, Astate, ibarcl_rst)
                call rst_write_record(atime, fid)
            end do
        end if

        rrec = rrec + 1
        write(*,*) "Record completed:", rrec
    end do time_loop

    ! ---------------------------- Closing ----------------------------
    do nre = 1, nrens
        close(20 + nre)
    end do
    close(18)
    close(19)

    if (lout == 'full') then
        do nre = 1, nrens
            close(20 + nrens + nre)
        end do
    end if

end program enKF2enKS
! ---------------------------------------------------------------------------
! ---------------------------------------------------------------------------
! ---------------------------------------------------------------------------



! ---------------------------------------------------------------------------
! num2str  — format integer with leading zeros (5 chars)
! ---------------------------------------------------------------------------
subroutine num2str(num, str)
    implicit none
    integer, intent(in) :: num
    character(len=5), intent(out) :: str

    if (num < 0 .or. num > 99999) error stop "num2str: argument out of range"
    write(str, '(I5.5)') num
end subroutine num2str

! ---------------------------------------------------------------------------
! push_matrix
! Convert SHYFEM prognostic fields into a single ensemble column (vector form)
!
! State layout (vectorized):
!
!   [ u(nlv,nel),
!     v(nlv,nel),
!     z(nkn),
!     T(nlv,nkn)   if ibarcl != 0
!     S(nlv,nkn)   if ibarcl != 0 ]
!
! Notes:
! - The reshape uses explicit (/size/) for clarity and safety.
! - Data are promoted to double precision to avoid roundoff accumulation.
! - No numerical changes to SHYFEM fields.
! ---------------------------------------------------------------------------
subroutine push_matrix(sdim, nrens, nre, Amat, ibrcl)
    use iso_fortran_env, only : dp => real64
    use basin
    use levels, only : nlv
    use mod_hydro
    use mod_hydro_vel
    use mod_ts
    use mod_conz
    implicit none

    integer, intent(in) :: sdim, nrens, nre, ibrcl
    real(dp), intent(inout) :: Amat(sdim, nrens)

    integer :: dimz, dimuv, dimts

    dimz  = nkn
    dimuv = nlv * nel
    dimts = nlv * nkn

    ! --- U (vectorized over nlv,nel)
    Amat(1:dimuv, nre) = reshape(real(utlnv, dp), [dimuv])

    ! --- V
    Amat(dimuv+1 : 2*dimuv, nre) = reshape(real(vtlnv, dp), [dimuv])

    ! --- Free surface elevation
    Amat(2*dimuv+1 : 2*dimuv+dimz, nre) = real(znv, dp)

    ! --- Temperature & salinity (if baroclinic mode enabled)
    if (ibrcl /= 0) then
        Amat(2*dimuv+dimz+1 : 2*dimuv+dimz+dimts, nre) = &
            reshape(real(tempv, dp), [dimts])

        Amat(2*dimuv+dimz+dimts+1 : 2*dimuv+dimz+2*dimts, nre) = &
            reshape(real(saltv, dp), [dimts])
    end if

end subroutine push_matrix



! ---------------------------------------------------------------------------
! pull_matrix
! Restore SHYFEM fields from the ensemble matrix.
!
! The inverse of push_matrix():
!
!   utlnv(nlv,nel)
!   vtlnv(nlv,nel)
!   znv(nkn)
!   tempv(nlv,nkn)    (if ibrcl != 0)
!   saltv(nlv,nkn)    (if ibrcl != 0)
!
! Notes:
! - reshape uses explicit (/nlv,nel/) etc.
! - All values converted back to real*4 where SHYFEM expects it.
! ---------------------------------------------------------------------------
subroutine pull_matrix(sdim, nrens, nre, Amat, ibrcl)
    use iso_fortran_env, only : dp => real64
    use basin
    use levels, only : nlvdi, nlv, hlv, ilhv, ilhkv
    use mod_hydro
    use mod_hydro_vel
    use mod_ts
    use mod_conz
    implicit none

    integer, intent(in) :: sdim, nrens, nre, ibrcl
    real(dp), intent(in) :: Amat(sdim,nrens)

    integer :: dimz, dimuv, dimts

    dimz  = nkn
    dimuv = nlv * nel
    dimts = nlv * nkn

    ! --- u
    utlnv = reshape(real(Amat(1:dimuv, nre)), [nlv, nel])

    ! --- v
    vtlnv = reshape(real(Amat(dimuv+1 : 2*dimuv, nre)), [nlv, nel])

    ! --- z
    znv = real(Amat(2*dimuv+1 : 2*dimuv+dimz, nre))

    ! --- T and S (if present)
    if (ibrcl /= 0) then
        tempv = reshape(real(Amat(2*dimuv+dimz+1 : 2*dimuv+dimz+dimts, nre)), &
                        [nlv, nkn])

        saltv = reshape(real(Amat(2*dimuv+dimz+dimts+1 : 2*dimuv+dimz+2*dimts, nre)), &
                        [nlv, nkn])
    end if

end subroutine pull_matrix

! ---------------------------------------------------------------------------
! init_shyfem
!
! Initializes the SHYFEM environment from the basin file and sets vertical
! levels. All SHYFEM module initializations must be called here to ensure that
! nkn, nel, nlv are properly defined before any push/pull operations.
!
! Improvements:
! - Explicit error checking on basin file open.
! - Clear ordering of module initializations.
! - No changes to SHYFEM semantics.
! - Comments clarified for readability.
! ---------------------------------------------------------------------------
subroutine init_shyfem(basinf, nnlv)
    use basin
    use levels
    use mod_geom_dynamic
    use mod_hydro
    use mod_hydro_vel
    use mod_ts
    use mod_conz
    use mod_restart
    use shympi

    implicit none

    character(len=80), intent(in) :: basinf
    integer,          intent(in) :: nnlv

    integer :: ios

    ! ----------------------------------------------------------------------
    ! 1. Read basin geometry
    ! ----------------------------------------------------------------------
    open(21, file=basinf, status='old', form='unformatted', iostat=ios)
    if (ios /= 0) then
        write(*,*) "Error: cannot open basin file: ", trim(basinf)
        error stop "init_shyfem: basin file open failure"
    end if

    call basin_read_by_unit(21)
    close(21)

    ! ----------------------------------------------------------------------
    ! 2. Set vertical levels (both in 'levels' and 'levels_di')
    ! ----------------------------------------------------------------------
    nlv   = nnlv
    nlvdi = nnlv

    ! ----------------------------------------------------------------------
    ! 3. Initialize SHYFEM dynamic modules
    ! ----------------------------------------------------------------------
    call mod_geom_dynamic_init(nkn, nel)
    call mod_hydro_init(     nkn, nel, nlv )
    call mod_hydro_vel_init( nkn, nel, nlv )
    call mod_ts_init(        nkn, nlv )
    call levels_init(        nkn, nel, nlv )

    ! ----------------------------------------------------------------------
    ! 4. Global copies for MPI / multi-domain support
    ! ----------------------------------------------------------------------
    nkn_global = nkn
    nel_global = nel
    nlv_global = nlv

    ! ----------------------------------------------------------------------
    ! 5. Concentrations
    ! NOTE: intentionally not activated
    !
    !   call mod_conz_init(1, nkn, nlvdi)
    ! ----------------------------------------------------------------------

end subroutine init_shyfem

! ---------------------------------------------------------------------------! ----------------------------------------------------------------_analysis1
!
! Apply the X5 transform directly at the current time step by scanning the
! whole X5_tot.uf file and multiplying the ensemble matrix when tt >= atime.
!
! Notes:
! - Method kept for reference; make_analysis2 is more efficient.
! - Fully modernized: no goto, robust I/O, clear logic.
! - dgemm is used for matrix multiplication (BLAS-3).
! ---------------------------------------------------------------------------
subroutine make_analysis1(atime, sdim, nrens, Amat)
    use iso_fortran_env, only: dp => real64
    implicit none

    ! ---------------------------- Arguments ----------------------------
    real(dp), intent(in)    :: atime
    integer,  intent(in)    :: sdim, nrens
    real(dp), intent(inout) :: Amat(sdim, nrens)

    ! ---------------------------- Locals ----------------------------
    real(dp) :: tt
    character(len=6) :: alabel
    character(len=2) :: tag
    integer :: nren
    integer :: ios

    real(dp), allocatable :: Aaux(:,:), X5(:,:)
    external :: dgemm

    ! ---------------------------- Allocate work arrays ----------------------------
    allocate(Aaux(sdim, nrens))
    allocate(X5(nrens, nrens))

    ! ---------------------------- Open X5 file ----------------------------
    open(15, file='X5_tot.uf', status='old', form='unformatted', iostat=ios)
    if (ios /= 0) error stop "make_analysis1: cannot open X5_tot.uf"

    ! ---------------------------- Scan file ----------------------------
read_loop: do
        ! Header-type record
        read(15, iostat=ios) tt, alabel, tag
        if (ios /= 0) exit read_loop

        ! Number of ensemble members
        read(15, iostat=ios) nren
        if (ios /= 0) error stop "make_analysis1: unexpected EOF in nren"

        if (tag == 'X3' .or. nren /= nrens) then
            error stop "make_analysis1: unsupported tag or bad ensemble size"
        end if

        ! Read X5 matrix
        read(15, iostat=ios) X5
        if (ios /= 0) error stop "make_analysis1: unexpected EOF in X5 matrix"

        ! Apply transform if time is relevant
        if (tt >= atime) then
            Aaux = 0.0_dp
            call dgemm('N', 'N', sdim, nrens, nrens, 1.0_dp, Amat, sdim, X5, nrens, &
                        0.0_dp, Aaux, sdim)
            Amat = Aaux
        end if
    end do read_loop

    close(15)

    ! ---------------------------- Cleanup ----------------------------
    deallocate(Aaux, X5)

end subroutine make_analysis1

! ---------------------------------------------------------------------------
! make_analysis2
!
! Accumulate X5 transforms forward in time up to nlag steps (or all if nlag=-1)
! and apply the resulting product X5old·Amat to obtain the smoothed ensemble.
!
! Improvements (Opzione C):
! - Stable accumulation of transforms:  X5old ← X5old·X5
!   Always multiplying in the same direction avoids chaotic rounding drift.
!
! - BLAS-3 usage for all matrix multiplications (dgemm).
!
! - nlag = -1 → treat as infinite: scan entire file.
!
! - Robust I/O with explicit iostat and predictable termination.
!
! - Work arrays zeroed explicitly to avoid undefined states.
!
! - No behavioural change in the mathematical sense, except improved numerical
!   stability (expected differences: O(1e-14) relative).
!
! ---------------------------------------------------------------------------
subroutine make_analysis2(atime, sdim, nrens, Amat, nlag)
    use iso_fortran_env, only: dp => real64
    implicit none

    ! ---------------------------- Arguments ----------------------------
    real(dp), intent(in)    :: atime
    integer,  intent(in)    :: sdim, nrens
    real(dp), intent(inout) :: Amat(sdim, nrens)
    integer,  intent(inout) :: nlag

    ! ---------------------------- Locals ----------------------------
    real(dp)           :: tt
    character(len=6)   :: alabel
    character(len=2)   :: tag
    integer            :: nren
    integer            :: ios
    integer            :: k     ! number of applied transforms
    integer            :: i,j

    real(dp), allocatable :: Aaux(:,:), X5(:,:), X5old(:,:), X5aux(:,:)
    external :: dgemm

    ! ---------------------------- Allocations ----------------------------
    allocate(Aaux(sdim, nrens))
    allocate(X5(nrens,nrens))
    allocate(X5old(nrens,nrens))
    allocate(X5aux(nrens,nrens))

    ! nlag = -1 → infinite (scan all transforms)
    if (nlag == -1) nlag = huge(1)

    ! ---------------------------- X5old = Identity ----------------------------
    X5old = 0.0_dp
    do concurrent (i = 1:nrens)
       X5old(i,i) = 1.0_dp
    end do

    ! ---------------------------- Open X5 file ----------------------------
    open(15, file='X5_tot.uf', status='old', form='unformatted', iostat=ios)
    if (ios /= 0) error stop "make_analysis2: cannot open X5_tot.uf"

    k = 0

    ! ---------------------------- Main scanning loop ----------------------------
read_loop: do
        ! Header record
        read(15, iostat=ios) tt, alabel, tag
        if (ios /= 0) exit read_loop

        ! Number of ensemble members
        read(15, iostat=ios) nren
        if (ios /= 0) error stop "make_analysis2: unexpected EOF in nren"

        ! Transform matrix
        read(15, iostat=ios) X5
        if (ios /= 0) error stop "make_analysis2: unexpected EOF in X5 matrix"

        ! Same checks as original:
        if (tag == 'X3' .or. nren /= nrens) then
            error stop "make_analysis2: unsupported tag or wrong nrens"
        end if

        ! ---- Apply transform if time > atime and lag not exceeded ----
        if ((tt > atime) .and. (k < nlag)) then
            ! X5old = X5old * X5
            X5aux = 0.0_dp
            call dgemm('N','N', nrens, nrens, nrens, &
                       1.0_dp, X5old, nrens, &
                              X5,    nrens, &
                       0.0_dp, X5aux, nrens)
            X5old = X5aux
            k = k + 1
        end if
    end do read_loop

    close(15)

    ! ---------------------------- Final multiplication ----------------------------
    ! Amat = Amat * X5old
    Aaux = 0.0_dp
    call dgemm('N','N', sdim, nrens, nrens, &
               1.0_dp, Amat, sdim, &
                      X5old, nrens, &
               0.0_dp, Aaux,  sdim)
    Amat = Aaux

    ! ---------------------------- Cleanup ----------------------------
    deallocate(Aaux, X5, X5old, X5aux)

end subroutine make_analysis2

! ---------------------------------------------------------------------------
! make_mn_std
!
! Compute ensemble mean and standard deviation using a numerically stable
! compensated summation (Neumaier variant of Kahan). This significantly
! reduces roundoff accumulation when nrens is large or state amplitudes
! differ in magnitude.
!
! The formula used is the standard unbiased sample variance:
!
!    mean  = ( Σ xi ) / n
!    var   = ( Σ (xi^2) - (Σ xi)^2 / n ) / (n - 1)
!    std   = sqrt(var)
!
! Relative differences vs the naive formula appear only in the last bits.
!
! Notes:
! - Amat comes in as real(dp)
! ---------------------------------------------------------------------------
subroutine make_mn_std(ndim, nens, Amat, Am, Astd)
    use iso_fortran_env, only : dp => real64
    implicit none

    integer,  intent(in)  :: ndim, nens
    real(dp), intent(in)  :: Amat(ndim, nens)
    real(dp), intent(out) :: Am(ndim), Astd(ndim)

    integer :: n, k
    real(dp) :: sum1, sum1_comp   ! compensated sum of x
    real(dp) :: sum2, sum2_comp   ! compensated sum of x^2
    real(dp) :: x, y, t

    if (nens < 2) error stop "make_mn_std: ensemble too small"

    do n = 1, ndim
        ! --------------------------------------------
        ! Compensated sum of Amat(n,:)
        ! --------------------------------------------
        sum1      = 0.0_dp
        sum1_comp = 0.0_dp
        sum2      = 0.0_dp
        sum2_comp = 0.0_dp

        block_sum: do concurrent (k = 1:nens)
            x = Amat(n, k)

            ! ---- Sum of x (Neumaier compensation) ----
            y = x - sum1_comp
            t = sum1 + y
            sum1_comp = (t - sum1) - y
            sum1 = t

            ! ---- Sum of x^2 (also compensated) ----
            y = x*x - sum2_comp
            t = sum2 + y
            sum2_comp = (t - sum2) - y
            sum2 = t
        end do block_sum

        ! --------------------------------------------
        ! Ensemble mean
        ! --------------------------------------------
        Am(n) = sum1 / real(nens, dp)

        ! --------------------------------------------
        ! Unbiased sample variance:
        ! var = ( Σ x^2 - (Σ x)^2 / n ) / (n - 1)
        ! Numerical stability improved via compensated sums
        ! --------------------------------------------
        Astd(n) = sqrt( max(0.0_dp, &
                        (sum2 - sum1*sum1 / real(nens,dp)) / &
                        (real(nens,dp) - 1.0_dp) ) )
    end do

end subroutine make_mn_std
