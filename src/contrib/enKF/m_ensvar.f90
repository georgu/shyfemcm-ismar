module m_ensvar
  !===============================================================================
  !  PURPOSE
  !    Provides the subroutine `ensvar` to compute the ensemble variance
  !    along the ensemble dimension (columns) of matrix A.
  !
  !  FORMULA
  !      var(i) = 1/(nrens-1) * Σ_j  ( A(i,j) - ave(i) )²
  !
  !    - The routine expects that the ensemble mean `ave(:)` has already
  !      been computed by the caller.
  !===============================================================================
  use iso_fortran_env, only : dp => real64
  implicit none
  private
  public :: ensvar
contains

  !=============================================================================
  !  SUBROUTINE: ensvar
  !
  !  INPUTS:
  !    nx        - length of the state vector  
  !    nrens     - number of ensemble members  
  !    A(nx,nrens)  - ensemble matrix  
  !    ave(nx)      - ensemble mean (precomputed)
  !
  !  OUTPUTS:
  !    var(nx)   - ensemble variance
  !
  !  REMARKS:
  !    - Requires nrens >= 2 because it uses the unbiased estimator:
  !          var = 1/(nrens-1) * sum (...)
  !=============================================================================
  subroutine ensvar(A, ave, var, nx, nrens)
    integer,  intent(in)  :: nx
    integer,  intent(in)  :: nrens
    real(dp), intent(in)  :: A(nx, nrens)
    real(dp), intent(in)  :: ave(nx)
    real(dp), intent(out) :: var(nx)

    integer :: j

    if (nrens < 2) then
       error stop 'ensvar: nrens must be >= 2'
    end if

    ! Initialize accumulator
    var = 0.0_dp

    ! Sum squared deviations across all ensemble members
    do j = 1, nrens
       var(:) = var(:) + (A(:,j) - ave(:)) * (A(:,j) - ave(:))
    end do

    ! Apply the unbiased normalization factor
    var = var / real(nrens - 1, dp)

  end subroutine ensvar

end module m_ensvar
