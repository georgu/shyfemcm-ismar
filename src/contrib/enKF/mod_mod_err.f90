!======================================================================
!  Module: mod_mod_err
!
!  Purpose:
!    Build and manage the model-error component for an augmented EnKF
!    state. A temporally correlated (AR(1)-like) model error q is used:
!        q_{k} = alpha * q_{k-1} + sqrt(1 - alpha^2) * w_k
!    where w_k is spatially correlated white noise generated on a
!    supersampled grid, rotated by theta_er.
!
!  Precision:
!    - All time/coefficients in double precision (consistent with EnKF).
!    - Spatial fields retain their kinds as defined in the state types.
!
!  Files:
!    - Reads parameters from   'mod_err.info'
!    - (Optionally) reads old  'an<step>_moderr.bin'
!    - Writes model error file 'an<step>_moderr.bin'
!
!  Status codes / errors:
!    - Uses ERROR STOP with explicit messages on I/O/logic errors.
!
!======================================================================
module mod_mod_err
  use mod_init_enkf
  use mod_mod_states
  use mod_ens_state
  implicit none

  !--------------------------------------------------------------------
  ! Parameters for the computation of the model error
  !--------------------------------------------------------------------
  integer            :: nx_er, ny_er       ! number of x and y grid points
  integer            :: fmult_er           ! multiplier for supersampling (noise generation)
  real               :: theta_er           ! rotation of random fields (0 = East, anticlockwise)
  real               :: rsigma             ! relative (spatial) error amplitude
  double precision   :: dt_er              ! time between two analysis steps
  double precision   :: tau_er             ! time decorrelation (>= dt_er): dq/dt = -(1/tau)*q

  ! Augmented state containers
  type(qstates), allocatable :: Abk_aug(:)      ! augmented state (Abk, qA) after adding model error
  type(states),  allocatable, private :: qA(:)  ! model error fields (same geometry as Abk)

contains

!======================================================================
!  info_moderr
!  -----------
!  Read model-error parameters from 'mod_err.info'.
!======================================================================
  subroutine info_moderr
    implicit none
    integer :: ios

    open(22, file='mod_err.info', status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'info_moderr: cannot open mod_err.info'

    read(22, *, iostat=ios) nx_er, ny_er, fmult_er, theta_er, rsigma, dt_er, tau_er
    if (ios /= 0) error stop 'info_moderr: malformed mod_err.info'
    close(22)

    if (tau_er < dt_er) error stop 'info_moderr: tau_er must be >= dt_er'
  end subroutine info_moderr

!======================================================================
!  push_aug
!  --------
!  Create the augmented ensemble state by:
!   1) Generating new spatial noise (qA) on a supersampled grid,
!   2) Loading old error (if available) and blending with AR(1) alpha,
!   3) Adding model error to background Abk: Abk <- Abk + sqrt(dt)*rsigma*rho*qA,
!   4) Packing (Abk, qA) into Abk_aug.
!======================================================================
  subroutine push_aug
    implicit none

    integer           :: nst
    double precision  :: alpha, rho
    type(states)      :: Aaux             ! temporary container for scaled fields
    real              :: kvec(nnkn, nrens)
    integer           :: nx, ny, fmult
    real              :: theta
    real              :: mfact_r          ! real scalar for interface compatibility
    double precision  :: mfact            ! double precision factor
    integer           :: ne, ie, n, k

    write(*,*) ''
    write(*,*) 'Creating an augmented state with model errors'
    write(*,*) ''

    !---------------------------------------
    ! Define temporal parameters for AR(1)
    !---------------------------------------
    nst = nanal
    if (tau_er < dt_er) error stop 'push_aug: parameter error (tau_er < dt_er)'

    alpha = 1.0d0 - (dt_er / tau_er)
    if (alpha < 0.0d0) alpha = 0.0d0
    if (alpha > 1.0d0) alpha = 1.0d0

    ! rho scales the cumulative effect over nst steps 
    rho = sqrt( (1.0d0 - alpha)**2 / &
                ( dt_er * ( dble(nst) - 2.0d0*alpha - dble(nst)*alpha**2 + 2.0d0*alpha*(dble(nst)+1.0d0) ) ) )

    !---------------------------------------
    ! Make new white-noise field (spatial)
    !---------------------------------------
    nx     = nx_er
    ny     = ny_er
    fmult  = fmult_er
    theta  = theta_er

    !call make_2Dpert(kvec, nnkn, nrens, fmult, theta, nx, ny)	!ggu bug
    call make_2Dpert(kvec, nnkn, nrens)

    allocate(qA(nrens))
    do ne = 1, nrens
      call allocate_states(qA(ne), nnkn, nnel, nnlv)
      qA(ne) = 0.0            ! zero the state (component-wise)
      qA(ne)%z = kvec(:, ne)  ! store the generated noise into, e.g., free-surface component
    end do

    !---------------------------------------
    ! Blend with old error if it exists:
    !   qA <- alpha*qA_old + sqrt(1 - alpha^2)*qA
    !---------------------------------------
    call load_error(alpha)

    !---------------------------------------
    ! Add model error to background:
    !   Abk <- Abk + (sqrt(dt)*rsigma*rho) * ( qA ⊙ Abk )  (component-wise scaling)
    ! Original code uses Aaux = Abk(ne)*mfact; then Aaux = qA(ne)*Aaux; Abk += Aaux
    !---------------------------------------
    mfact = sqrt(dt_er) * dble(rsigma) * rho
    mfact_r = real(mfact)  ! if overloaded operators expect REAL scalar in your code base

    do ne = 1, nrens
      ! Scale background by mfact
      Aaux = Abk(ne) * mfact_r
      ! Apply spatially varying error (component-wise product)
      Aaux = qA(ne) * Aaux
      ! Add to background
      Abk(ne) = Abk(ne) + Aaux
    end do

    !---------------------------------------
    ! Build augmented states Abk_aug = (Abk, qA)
    !---------------------------------------
    allocate(Abk_aug(nrens))
    do ne = 1, nrens
      call allocate_qstates(Abk_aug(ne), nnkn, nnel, nnlv)
      call push_qstate(Abk(ne), qA(ne), Abk_aug(ne))
    end do

    deallocate(qA)
    deallocate(Abk)
  end subroutine push_aug

!======================================================================
!  pull_aug
!  --------
!  Split augmented state into (Abk, qA) and write qA to file
!  'an<nanal>_moderr.bin' (unformatted).
!======================================================================
  subroutine pull_aug
    implicit none

    type(states), allocatable :: qA(:)
    character(len=5)  :: nal
    character(len=18) :: fname
    integer :: ne, ios

    write(*,*) 'Saving model errors'

    allocate(Abk(nrens), qA(nrens))
    do ne = 1, nrens
      call allocate_states(qA(ne), nnkn, nnel, nnlv)
      call pull_qstate(Abk(ne), qA(ne), Abk_aug(ne))
    end do

    call num2str(nanal, nal)
    fname = 'an'//nal//'_moderr.bin'

    open(33, file=fname, form='unformatted', status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'pull_aug: cannot open output moderr file'

    do ne = 1, nrens
      write(33) qA(ne)%u
      write(33) qA(ne)%v
      write(33) qA(ne)%z
      write(33) qA(ne)%t
      write(33) qA(ne)%s
    end do
    close(33)

    deallocate(Abk_aug, qA)
  end subroutine pull_aug

!======================================================================
!  load_error
!  ----------
!  If old model-error file exists, read qA_old and blend:
!      qA <- alpha*qA_old + sqrt(1 - alpha^2) * qA
!  No goto; structured flow with early return when file is absent.
!======================================================================
  subroutine load_error(alpha)
    implicit none
    double precision, intent(in) :: alpha

    logical :: bfile
    integer :: ne, ios
    character(len=5)  :: nal
    character(len=18) :: fname

    type(states), allocatable :: qA1(:)
    type(states)              :: A2

    double precision :: s1

    ! Old analysis step file name
    call num2str(nanal-1, nal)
    fname = 'an'//nal//'_moderr.bin'

    inquire(file=fname, exist=bfile)
    if (.not. bfile) then
      write(*,*) 'Model error files not found (skipping old-error blending).'
      return
    end if

    write(*,*) 'Loading model error from files: ', trim(fname)

    open(22, file=fname, status='old', form='unformatted', action='read', iostat=ios)
    if (ios /= 0) error stop 'load_error: cannot open old moderr file'

    allocate(qA1(nrens))
    do ne = 1, nrens
      call allocate_states(qA1(ne), nnkn, nnel, nnlv)
      read(22, iostat=ios) qA1(ne)%u; if (ios /= 0) error stop 'load_error: read u failed'
      read(22, iostat=ios) qA1(ne)%v; if (ios /= 0) error stop 'load_error: read v failed'
      read(22, iostat=ios) qA1(ne)%z; if (ios /= 0) error stop 'load_error: read z failed'
      read(22, iostat=ios) qA1(ne)%t; if (ios /= 0) error stop 'load_error: read t failed'
      read(22, iostat=ios) qA1(ne)%s; if (ios /= 0) error stop 'load_error: read s failed'
    end do
    close(22)

    s1 = sqrt(max(0.0d0, 1.0d0 - alpha*alpha))

    do ne = 1, nrens
      A2     = qA1(ne) * real(alpha)     ! alpha * q_old
      qA(ne) = A2 + ( qA(ne) * real(s1) ) ! + sqrt(1-alpha^2) * q_new
    end do

    deallocate(qA1)
  end subroutine load_error

end module mod_mod_err
