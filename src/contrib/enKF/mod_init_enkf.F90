!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
!
! Read the initial settings
!
module mod_init_enkf
  implicit none

  integer :: nnlv, nrens, nanal
  character (len=80) :: basfile, obsfile
  character (len=80) :: ostring
  double precision :: atime_an
  integer :: bnew_ens


contains


  subroutine read_info

  use iso8601
  use mod_para
  implicit none
  integer ierr
  integer date,time

  open(20, file='analysis.info', status='old')

  read(20,*) nnlv	! number of vertical levels
  read(20,*) nrens	! number of ens members
  read(20,*) nanal	! analysis step
  read(20,*) basfile	! name of bas file (no extension)
  read(20,*) ostring	! current time of the observations, string format
  read(20,*) obsfile	! name of obs file list
  read(20,*) bnew_ens	! 1 to create a new initial ens of states
  read(20,*) rmode	! analysis method
  read(20,*) is_local	! localisation or not

  close(20)

  if ((rmode /= 11).and.(rmode /= 12).and.(rmode /= 13).and.(rmode /= 21).and.(rmode /= 22).and.(rmode /= 23)) &
	  error stop 'Wrong analysis method.'

  if ((is_local /= 0).and.(is_local /= 1)) error stop 'Wrong localisation option.'

  if (mod(nrens,2) == 0)&
     error stop 'read_info: n of ens members must be odd, with the control as first.'

  call string2date(trim(ostring),date,time,ierr)
  if (ierr /= 0) error stop 'read_info: invalid date string'
  call dts_to_abs_time(date,time,atime_an)

  write(*,*) 'time of the analysis step: ',trim(ostring)
  write(*,*) 'n. of ens members: ',nrens
  
  end subroutine read_info

end module mod_init_enkf
