! Read the initial settings
!
module mod_init_enkf
  implicit none

  integer :: nrens, nanal
  character (len=80) :: basfile, obsfile
  double precision :: atime
  integer :: bnew_ens
  integer :: bmod_err


contains


  subroutine read_info

  implicit none

  open(20, file='analysis.info', status='old')

  read(20,*) nrens	! number of ens members
  read(20,*) nanal		! analysis step
  read(20,*) basfile	! name of bas file (no extension)
  read(20,*) atime	! current time of the observations
  read(20,*) obsfile	! name of obs file list
  read(20,*) bnew_ens	! 1 to create a new initial ens of states
  read(20,*) bmod_err	! 1 to use an augmented state with mod err

  close(20)

  if (mod(nrens,2) == 0)&
     error stop 'read_info: n of ens members must be odd, with the control as first.'

  write(*,*) 'time: ',atime
  write(*,*) 'n. of ens members: ',nrens
  write(*,*) 'bnew_ens: ',bnew_ens
  write(*,*) 'bmod_err: ',bmod_err
  
  end subroutine read_info

!***************************************************************

  subroutine set_model_params
  use mod_dimensions
  use basin
  implicit none

  open(21, file=basfile, status='old', form='unformatted')
  call basin_read_by_unit(21)
  close(21)

  if ((nkn /= nnkn).or.(nel /= nnel)) error stop "read_basin: dim error"

  ! init modules
  !
  call init_shyfem_vars(nnkn,nnel,nnlv)

  ! addpar for restart
  !
  call add_rst_params

  end subroutine set_model_params


end module mod_init_enkf
