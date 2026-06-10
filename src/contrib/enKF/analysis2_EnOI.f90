subroutine analysis2_EnOI(A, psi, D, R, S, ndim, nrens, nrobs, verbose)
!=======================================================================
!  PURPOSE:
!    Computes the analysed state vector psi using the EnOI algorithm
!    with a static ensemble A and innovations D.
!
!  DESCRIPTION:
!    - Builds observation-space covariance SS' + N*R
!    - Performs SVD-based pseudo inversion
!    - Computes representer or full-space update depending on dims
!    - Ensures double precision and safe allocations
!=======================================================================

  use iso_fortran_env, only: dp => real64
  use m_multa
  implicit none

  !======================
  ! Arguments
  !======================
  integer, intent(in) :: ndim, nrens, nrobs
  real(dp), intent(in)    :: A(ndim,nrens)     ! static ensemble
  real(dp), intent(inout) :: psi(ndim)         ! state vector updated
  real(dp), intent(in)    :: D(nrobs)          ! innovations
  real(dp), intent(inout) :: R(nrobs,nrobs)    ! observation covariance (modified)
  real(dp), intent(in)    :: S(nrobs,nrens)    ! HA'
  logical, intent(in)     :: verbose

  !======================
  ! Local allocatables
  !======================
  real(dp), allocatable :: U(:,:), V(:,:), sig(:), work(:)
  real(dp), allocatable :: X1(:,:), X2(:,:), X3(:,:), Reps(:,:), I_N(:,:)
  real(dp), allocatable :: X4(:), X5(:)
  real(dp) :: sigsum, sigsum1, oneobs(1,1)
  integer  :: ierr, nrsigma, i, j, lwork

  integer, parameter :: target = 3*22+1   ! as in original code

  external dgemm, dgemv

  !=====================================================================
  ! 1) BUILD R := SS' + N*R
  !=====================================================================
  if (nrobs > 1) then

     call dgemm('N','T', nrobs, nrobs, nrens, 1.0_dp, &
                 S, nrobs, S, nrobs, real(nrens,dp), R, nrobs)

     ! Allocate SVD arrays
     allocate(U(nrobs,nrobs), V(nrobs,nrobs), sig(nrobs), stat=ierr)
     lwork = 2*max(3*nrobs+nrobs, 5*nrobs)
     allocate(work(lwork), stat=ierr)
     sig = 0.0_dp

     !================================================================
     ! SVD: R = U * diag(sig) * V^T
     !================================================================
     call dgesvd('A','A', nrobs, nrobs, R, nrobs, sig, U, nrobs, V, nrobs, &
                 work, lwork, ierr)
     deallocate(work)

     if (ierr /= 0) then
        print *, 'ERROR: dgesvd returned ierr = ', ierr
        stop
     endif

     !===========================================================
     ! Select dominant singular values (cumulative > 99.9%)
     ! Invert them for pseudo-inverse
     !===========================================================
     sigsum  = sum(sig)
     sigsum1 = 0.0_dp
     nrsigma = 0

     do i = 1, nrobs
        if (sigsum1/sigsum < 0.999_dp) then
           nrsigma = nrsigma + 1
           sigsum1 = sigsum1 + sig(i)
           sig(i)  = 1.0_dp / sig(i)
        else
           sig(i:nrobs) = 0.0_dp
           exit
        endif
     end do

     if (verbose) then
        print *, 'Dominant singular values: ', nrsigma, ' share=', sigsum1/sigsum
     endif

     !===========================================================
     ! Build X1 = U * diag(sig)
     !===========================================================
     allocate(X1(nrobs,nrobs))
     do i = 1, nrobs
        do j = 1, nrobs
           X1(i,j) = sig(i) * U(j,i)
        end do
     end do
     deallocate(sig)

     !===========================================================
     ! X2 = X1 * D     (nrobs x 1)
     !===========================================================
     allocate(X2(nrobs,1))
     call dgemm('N','N', nrobs, 1, nrobs, 1.0_dp, X1, nrobs, D, nrobs, 0.0_dp, X2, nrobs)
     deallocate(X1)

     !===========================================================
     ! X3 = V^T * X2   (nrobs x 1)
     !===========================================================
     allocate(X3(nrobs,1))
     call dgemm('T','N', nrobs, 1, nrobs, 1.0_dp, V, nrobs, X2, nrobs, 0.0_dp, X3, nrobs)
     deallocate(V)
     deallocate(X2)

  else
     !=================================================================
     ! Special case: only 1 observation
     !=================================================================
     allocate(X3(nrobs,1))
     oneobs = matmul(S,transpose(S)) + real(nrens,dp)*R
     X3(:,1) = D / oneobs(1,1)
  endif


  !=====================================================================
  ! 2) REPRESENTERS or FULL UPDATE depending on size rule
  !=====================================================================
  if (2_8*ndim*nrobs < 1_8*nrens*(nrobs+ndim)) then
     !==============================================================
     ! REPRESENTERS APPROACH (few observations)
     !==============================================================
     if (verbose) print *, 'analysis: Representer approach is used'

     allocate(Reps(ndim,nrobs))
     call dgemm('N','T', ndim, nrobs, nrens, 1.0_dp, A, ndim, S, nrobs, 0.0_dp, Reps, ndim)

     ! Remove mean in obs space: I_N = I - (1/nrobs)*11^T
     allocate(I_N(nrobs,nrobs))
     I_N = -1.0_dp/real(nrobs,dp)
     do i=1,nrobs
        I_N(i,i) = I_N(i,i) + 1.0_dp
     end do

     call dgemm('N','N', ndim, nrobs, nrobs, 1.0_dp, Reps, ndim, I_N, nrobs, 0.0_dp, Reps, ndim)
     deallocate(I_N)

     ! psi ← psi + Reps * X3
     call dgemv('N', ndim, nrobs, 1.0_dp, Reps, ndim, X3, 1, 1.0_dp, psi, 1)

     deallocate(Reps)
     deallocate(X3)

  else
     !==============================================================
     ! FULL-SPACE UPDATE (many observations)
     !==============================================================
     allocate(X4(nrens))
     call dgemm('T','N', nrens, 1, nrobs, 1.0_dp, S, nrobs, X3, nrobs, 0.0_dp, X4, nrens)

     ! Remove ensemble mean: X5 = I_N * X4
     allocate(I_N(nrens,nrens))
     I_N = -1.0_dp / real(nrens,dp)
     do i = 1, nrens
       I_N(i,i) = I_N(i,i) + 1.0_dp
     end do

     allocate(X5(nrens))
     X5 = matmul(I_N, X4)

     deallocate(I_N)
     deallocate(X4)
     deallocate(X3)

     ! psi ← psi + A * X5
     call dgemv('N', ndim, nrens, 1.0_dp, A, ndim, X5, 1, 1.0_dp, psi, 1)
     deallocate(X5)

     if (verbose) print *, 'Analysis, TEM(1) =', psi(target)
  endif

end subroutine analysis2_EnOI
