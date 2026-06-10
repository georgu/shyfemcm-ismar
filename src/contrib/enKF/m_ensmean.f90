module m_ensmean
  !===============================================================================
  !  PURPOSE
  !    Provides the subroutine `ensmean` to compute the ensemble mean
  !    along the ensemble dimension (columns) of matrix A.
  !
  !===============================================================================
  use iso_fortran_env, only : dp => real64
  implicit none
  private
  public :: ensmean
contains

  !=============================================================================
  !  SUBROUTINE: ensmean
  !
  !  INPUTS:
  !    nx      - size of the state vector (number of rows)
  !    nrens   - number of ensemble members (columns of A)
  !    A(nx,nrens) - ensemble matrix; each column = one ensemble member
  !
  !  OUTPUT:
  !    ave(nx) - ensemble mean: ave = (1/nrens) * sum_j A(:,j)
  !
  !  REMARKS:
  !    - Simple, explicit loop over ensemble members.
  !    - Requires nrens >= 1.
  !=============================================================================
  subroutine ensmean(A, ave, nx, nrens)
    integer,  intent(in)  :: nx
    integer,  intent(in)  :: nrens
    real(dp), intent(in)  :: A(nx, nrens)
    real(dp), intent(out) :: ave(nx)

    integer :: j

    ! Basic input check
    if (nrens < 1) then
      error stop 'ensmean: nrens must be >= 1'
    end if

    ! Initialize the accumulator with the first member
    ave(:) = A(:,1)

    ! Accumulate the remaining ensemble members
    do j = 2, nrens
      ave(:) = ave(:) + A(:,j)
    end do

    ! Finalize with normalization
    ave = (1.0_dp / real(nrens, dp)) * ave

  end subroutine ensmean

end module m_ensmean
