module m_mean_preserving_rotation
  !=======================================================================
  ! Mean-preserving random rotation for EnKF SQRT (Sakov 2006/07 style):
  !   Find Up (nrens x nrens) such that:
  !     1) Up * Up^T = I        (orthonormal)
  !     2) Up * 1 = 1           (each row sums to 1 → preserves ensemble mean)
  !
  ! CRITICAL REFACTORING (v2.0):
  !   - Fixed Gram-Schmidt variance deflation via zero-mean centering of base vectors.
  !   - Converted workspace arrays to dynamically allocated Heap to avoid Stack Overflow.
  !   - Guaranteed explicit type safety for BLAS/LAPACK interaction parameters.
  !=======================================================================
  use iso_fortran_env, only : dp => real64
  use iso_fortran_env, only : ip => int32
  use m_randrot                       
  implicit none
  private
  public :: mean_preserving_rotation

  interface
     subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
       import dp, ip
       character(len=1), intent(in) :: transa, transb
       integer(ip),      intent(in) :: m, n, k, lda, ldb, ldc
       real(dp),         intent(in) :: alpha, beta
       real(dp),         intent(in) :: a(lda,*), b(ldb,*)
       real(dp),         intent(inout) :: c(ldc,*)
     end subroutine dgemm
  end interface

contains

  subroutine mean_preserving_rotation(Up, nrens)
    integer(ip), intent(in)  :: nrens
    real(dp),    intent(out) :: Up(nrens, nrens)

    ! Dynamic heap-allocated workspaces to prevent openMP Stack Overflows
    real(dp), allocatable :: B(:,:), Q(:,:), R(:,:), U(:,:), Upb(:,:)
    real(dp) :: row_mean

    integer(ip) :: i, j, k
    integer(ip) :: ld_nrens
    real(dp)    :: one, zero

    one  = 1.0_dp
    zero = 0.0_dp
    ld_nrens = nrens

    if (nrens <= zero) then
       write(*,*) 'mean_preserving_rotation: nrens must be positive.'
       Up = zero
       return
    end if

    if (nrens == one) then
       Up = one
       return
    end if

    allocate(B(nrens, nrens), Q(nrens, nrens), R(nrens, nrens))
    allocate(U(nrens-1, nrens-1), Upb(nrens, nrens))

    ! -------------------------------------------------------------------------
    ! 1) Build B: an orthonormal basis whose first column is 1/sqrt(nrens)
    ! -------------------------------------------------------------------------
    call random_number(B)                  
    
    !=======================================================================
    ! FIX: Raw uniform random numbers have a mean of 0.5. Projecting columns
    ! 2:N against a constant first column causes severe numerical cancellation 
    ! and loss of rank during Gram-Schmidt. We explicitly force columns 2:N 
    ! to have an intrinsic zero-mean configuration before running MGS.
    !=======================================================================
    do j = 2, nrens
       row_mean = sum(B(:,j)) / real(nrens, dp)
       B(:,j) = B(:,j) - row_mean
    end do
    
    ! Fix the target mean-preserving transformation vector in column 1
    B(:,1) = 1.0_dp / sqrt( real(nrens, dp) )

    ! Modified Gram-Schmidt (MGS) on columns of B
    R = 0.0_dp
    do k = 1, nrens
       R(k,k) = sqrt( dot_product(B(:,k), B(:,k)) )
       if (R(k,k) < 1.0e-12_dp) then
          ! Robust guard against pathological linear dependencies
          B(:,k) = 0.0_dp
          B(k,k) = 1.0_dp
          R(k,k) = sqrt( dot_product(B(:,k), B(:,k)) )
       end if
       B(:,k) = B(:,k) / R(k,k)
       do j = k+1, nrens
          R(k,j) = dot_product(B(:,k), B(:,j))
          B(:,j) = B(:,j) - B(:,k) * R(k,j)
       end do
    end do

    ! ---------------- ---------------------------------------------------------
    ! 2) Random orthonormal U of size (nrens-1) x (nrens-1)
    ! -------------------------------------------------------------------------
    call randrot(U, nrens-1_ip)

    ! -------------------------------------------------------------------------
    ! 3) Build Upb = diag(1, U)  (nrens x nrens)
    ! -------------------------------------------------------------------------
    Upb = 0.0_dp
    Upb(1,1) = 1.0_dp
    Upb(2:nrens, 2:nrens) = U(1:nrens-1, 1:nrens-1)

    ! -------------------------------------------------------------------------
    ! 4) Up = B * Upb * B^T via explicit type-safe BLAS parameters
    ! -------------------------------------------------------------------------
    call dgemm('N', 'N', ld_nrens, ld_nrens, ld_nrens, one, B, ld_nrens, Upb, ld_nrens, zero, Q, ld_nrens)
    call dgemm('N', 'T', ld_nrens, ld_nrens, ld_nrens, one, Q, ld_nrens, B, ld_nrens, zero, Up, ld_nrens)

    deallocate(B, Q, R, U, Upb)

  end subroutine mean_preserving_rotation

end module m_mean_preserving_rotation
