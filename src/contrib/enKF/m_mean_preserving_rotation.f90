module m_mean_preserving_rotation
  !!
  !! Mean-preserving random rotation for EnKF SQRT (Sakov 2006/07 style):
  !!   Find Up (nrens x nrens) such that:
  !!     1) Up * Up^T = I        (orthonormal)
  !!     2) Up * 1 = 1           (each row sums to 1 → preserves ensemble mean)
  !!
  !! Construction:
  !!   - Build an orthonormal basis B with the first column equal to 1/sqrt(nrens)
  !!   - Draw an arbitrary orthonormal U of size (nrens-1) x (nrens-1)
  !!   - Form Upb = diag(1, U)
  !!   - Set Up = B * Upb * B^T
  !!
  use iso_fortran_env, only : dp => real64
  use iso_fortran_env, only : ip => int32
  use m_randrot                       ! Assumed to provide: randrot(U, n) in double precision
  implicit none

  ! Provide an explicit interface for BLAS dgemm (double precision).
  ! Adjust name mangling/linking flags as needed when linking BLAS/LAPACK.
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
    !!
    !! Arguments (unchanged from your original):
    !!   Up    [out] : real(dp) :: Up(nrens, nrens)  ← mean-preserving orthonormal rotation
    !!   nrens [in]  : integer(ip)                   ← ensemble size
    !!
    integer(ip), intent(in)  :: nrens
    real(dp),    intent(out) :: Up(nrens, nrens)

    ! Workspace
    real(dp) :: B(nrens, nrens)     ! Orthonormal basis with e1 = 1/sqrt(nrens)
    real(dp) :: Q(nrens, nrens)     ! GEMM workspace
    real(dp) :: R(nrens, nrens)     ! R factors used during MGS (only upper triangle is touched)
    real(dp) :: U(nrens-1, nrens-1) ! Random orthonormal (nrens-1)x(nrens-1)
    real(dp) :: Upb(nrens, nrens)   ! Block-diagonal [1, U]

    integer(ip) :: i, j, k
    real(dp)    :: one, zero

    one  = 1.0_dp
    zero = 0.0_dp

    ! Basic checks
    if (nrens <= 0) then
       write(*,*) 'mean_preserving_rotation: nrens must be positive.'
       Up = 0.0_dp
       return
    end if

    if (nrens == 1_ip) then
       ! Trivial 1x1 case: Up = [1]
       Up = 1.0_dp
       return
    end if

    ! -------------------------------------------------------------------------
    ! 1) Build B: an orthonormal basis whose first column is 1/sqrt(nrens)
    ! -------------------------------------------------------------------------
    call random_number(B)                  ! Uniform(0,1) entries
    B(:,1) = 1.0_dp / sqrt( real(nrens, dp) )

    ! Modified Gram-Schmidt on columns of B, overwriting in place.
    ! Only R(k,k) and R(k,j) for j>k are used during the process.
    R = 0.0_dp
    do k = 1, nrens
       R(k,k) = sqrt( dot_product(B(:,k), B(:,k)) )
       if (R(k,k) == 0.0_dp) then
          ! Extremely unlikely, but guard against pathological random draw.
          B(:,k) = 0.0_dp
          B(1,k) = 1.0_dp
          R(k,k) = sqrt( dot_product(B(:,k), B(:,k)) )
       end if
       B(:,k) = B(:,k) / R(k,k)
       do j = k+1, nrens
          R(k,j) = dot_product(B(:,k), B(:,j))
          B(:,j) = B(:,j) - B(:,k) * R(k,j)
       end do
    end do
    ! (Optional) Orthonormality checks — kept but commented out.
    ! do k = 1, nrens
    !    do j = k, min(k+14, nrens)
    !       write(*,'(15f10.4)', advance='no') dot_product(B(:,j), B(:,k))
    !    end do
    !    write(*,*) ''
    ! end do

    ! -------------------------------------------------------------------------
    ! 2) Random orthonormal U of size (nrens-1) x (nrens-1)
    ! -------------------------------------------------------------------------
    call randrot(U, nrens-1)
    ! (Optional) Orthonormality checks — kept but commented out.
    ! do k = 1, nrens-1
    !    do j = k, min(k+14, nrens-1)
    !       write(*,'(15f10.4)', advance='no') dot_product(U(:,j), U(:,k))
    !    end do
    !    write(*,*) ''
    ! end do

    ! -------------------------------------------------------------------------
    ! 3) Build Upb = diag(1, U)  (nrens x nrens)
    ! -------------------------------------------------------------------------
    Upb = 0.0_dp
    Upb(1,1) = 1.0_dp
    Upb(2:nrens, 2:nrens) = U(1:nrens-1, 1:nrens-1)

    ! (Optional) Orthonormality checks — kept but commented out.
    ! do k = 1, nrens
    !    do j = k, min(k+14, nrens)
    !       write(*,'(15f10.4)', advance='no') dot_product(Upb(:,j), Upb(:,k))
    !    end do
    !    write(*,*) ''
    ! end do

    ! -------------------------------------------------------------------------
    ! 4) Up = B * Upb * B^T
    !    Use BLAS dgemm for efficiency and numerical robustness.
    ! -------------------------------------------------------------------------
    call dgemm('N', 'N', nrens, nrens, nrens, one,  B,  nrens, Upb, nrens, zero, Q,  nrens)
    call dgemm('N', 'T', nrens, nrens, nrens, one,  Q,  nrens, B,   nrens, zero, Up, nrens)

    ! (Optional) Diagnostics — kept but commented out.
    ! do k = 1, nrens
    !    write(*,*) 'Up row sum:', k, sum(Up(k,1:nrens))
    ! end do
    ! do k = 1, nrens
    !    do j = k, min(k+14, nrens)
    !       write(*,'(15f10.4)', advance='no') dot_product(Up(:,j), Up(:,k))
    !    end do
    !    write(*,*) ''
    ! end do

  end subroutine mean_preserving_rotation

end module m_mean_preserving_rotation
