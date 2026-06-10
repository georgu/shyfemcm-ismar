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
  integer            :: nx_er, ny_er
  integer            :: fmult_er
  double precision   :: theta_er
  double precision   :: rsigma
  double precision   :: dt_er
  double precision   :: tau_er

  ! Augmented state containers
  type(qstates), allocatable :: Abk_aug(:)
  type(states),  allocatable, private :: qA(:)

contains

!======================================================================
!  info_moderr
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
!======================================================================
  subroutine push_aug
    implicit none

    integer           :: nst, ne
    double precision  :: alpha, rho, mfact
    type(states)      :: Aaux
    double precision  :: kvec(nnkn, nrens)

    write(*,*) ''
    write(*,*) 'Creating an augmented state with model errors'
    write(*,*) ''

    !---------------------------------------
    ! AR(1) parameters
    !---------------------------------------
    nst = nanal
    if (tau_er < dt_er) error stop 'push_aug: tau_er < dt_er'

    alpha = 1.0d0 - dt_er / tau_er
    alpha = max(0.0d0, min(1.0d0, alpha))

    rho = sqrt( (1.0d0 - alpha)**2 / &
                ( dt_er * ( dble(nst)                         &
                          - 2.0d0*alpha                        &
                          - dble(nst)*alpha**2                 &
                          + 2.0d0*alpha*(dble(nst) + 1.0d0) ) ) )

    !---------------------------------------
    ! Generate spatial white noise
    !---------------------------------------
    call make_2Dpert(kvec, nnkn, nrens)

    allocate(qA(nrens))
    do ne = 1, nrens
      call allocate_states(qA(ne), nnkn, nnel, nnlv)
      qA(ne) = 0.0d0
      qA(ne)%z = kvec(:, ne)
    end do

    !---------------------------------------
    ! Blend with previous model error
    !---------------------------------------
    call load_error(alpha)

    !---------------------------------------
    ! Add model error to background
    !---------------------------------------
    mfact = sqrt(dt_er) * rsigma * rho

    do ne = 1, nrens
      Aaux    = Abk(ne) * mfact
      Aaux    = qA(ne) * Aaux
      Abk(ne) = Abk(ne) + Aaux
    end do

    !---------------------------------------
    ! Build augmented state
    !---------------------------------------
    allocate(Abk_aug(nrens))
    do ne = 1, nrens
      call allocate_qstates(Abk_aug(ne), nnkn, nnel, nnlv)
      call push_qstate(Abk(ne), qA(ne), Abk_aug(ne))
    end do

    deallocate(qA, Abk)
  end subroutine push_aug

!======================================================================
!  pull_aug
!======================================================================
  subroutine pull_aug
    implicit none

    type(states), allocatable :: qA(:)
    character(len=5)  :: nal
    character(len=18) :: fname
    integer           :: ne, ios

    write(*,*) 'Saving model errors'

    allocate(Abk(nrens), qA(nrens))
    do ne = 1, nrens
      call allocate_states(qA(ne), nnkn, nnel, nnlv)
      call pull_qstate(Abk(ne), qA(ne), Abk_aug(ne))
    end do

    call num2str(nanal, nal)
    fname = 'an'//nal//'_moderr.bin'

    open(33, file=fname, form='unformatted', status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'pull_aug: cannot open moderr file'

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
!======================================================================
  subroutine load_error(alpha)
    implicit none
    double precision, intent(in) :: alpha

    logical                       :: bfile
    integer                       :: ne, ios
    character(len=5)              :: nal
    character(len=18)             :: fname
    type(states), allocatable     :: qA1(:)
    type(states)                  :: A2
    double precision              :: s1

    call num2str(nanal-1, nal)
    fname = 'an'//nal//'_moderr.bin'

    inquire(file=fname, exist=bfile)
    if (.not. bfile) then
      write(*,*) 'Model error file not found; skipping blending'
      return
    end if

    open(22, file=fname, form='unformatted', action='read', iostat=ios)
    if (ios /= 0) error stop 'load_error: cannot open old moderr file'

    allocate(qA1(nrens))
    do ne = 1, nrens
      call allocate_states(qA1(ne), nnkn, nnel, nnlv)
      read(22) qA1(ne)%u
      read(22) qA1(ne)%v
      read(22) qA1(ne)%z
      read(22) qA1(ne)%t
      read(22) qA1(ne)%s
    end do
    close(22)

    s1 = sqrt(max(0.0d0, 1.0d0 - alpha*alpha))

    do ne = 1, nrens
      A2     = qA1(ne) * alpha
      qA(ne) = A2 + qA(ne) * s1
    end do

    deallocate(qA1)
  end subroutine load_error

end module mod_mod_err
