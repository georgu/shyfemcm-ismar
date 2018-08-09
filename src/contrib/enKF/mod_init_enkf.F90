!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
!
! Read the initial settings
!
module mod_init_enkf
  implicit none

  integer :: nrens, nanal
  character (len=80) :: basfile, obsfile
  character (len=80) :: ostring
  double precision :: atime
  integer :: bnew_ens
  integer :: mode_an


contains


  subroutine read_info

  use iso8601
  implicit none
  integer ierr
  integer date,time

  open(20, file='analysis.info', status='old')

  read(20,*) nrens	! number of ens members
  read(20,*) nanal		! analysis step
  read(20,*) basfile	! name of bas file (no extension)
  read(20,*) ostring	! current time of the observations, string format
  read(20,*) obsfile	! name of obs file list
  read(20,*) bnew_ens	! 1 to create a new initial ens of states
  read(20,*) mode_an	! 0 normal; 1 mod err; 2 mod params

  close(20)

  if (mod(nrens,2) == 0)&
     error stop 'read_info: n of ens members must be odd, with the control as first.'

  call string2date(trim(ostring),date,time,ierr)
  if (ierr /= 0) error stop 'read_info: invalid date string'
  call dts_to_abs_time(date,time,atime)

  write(*,*) 'time of the analysis step: ',trim(ostring)
  write(*,*) 'absolute time: ',atime
  write(*,*) 'n. of ens members: ',nrens
  write(*,*) 'ensemble state creation: ',bnew_ens
  write(*,*) 'mode of simulation: ',mode_an
  
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
