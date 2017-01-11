module m_read_ensemble
! Reads the modelled forecast from nrens restarts and stores it in the matrix A
   use mod_restart
contains

  subroutine read_ensemble
   use mod_enkf
   use mod_hydro
   use mod_ts

   implicit none

   integer ne

   type(states4) :: A4

   ! reads ens
   do ne = 1,nrens
      call rst_read(ne,na,tobs)
      A4%u = utlnv
      A4%v = vtlnv
      A4%z = znv
      if( ibarcl_rst.le.0 ) then
        tempv = 0.
        saltv = 0.
      end if 
      A4%t = tempv
      A4%s = saltv
      A(ne)=A4
   end do

   return
  end subroutine read_ensemble

!********************************************************

  subroutine rst_read(ne,nan,tt)

  use basin
  use mod_dimensions
  use mod_geom_dynamic
  use mod_hydro
  use mod_hydro_vel
  use mod_ts
  use mod_conz
  implicit none

  integer, intent(in) :: ne, nan
  real, intent(in) :: tt

  character(len=3) :: nrel,nal
  integer io
  integer it,nvers,nrec,iflag,ierr
  real atime
  integer date,time
  character(len=16) rstname

  integer, save :: nkn0,nel0,nlv0
  integer, save :: icall = 0

  call num2str(ne,nrel)
  call num2str(nan,nal)

  rstname='an'//nal//'_'//'en'//nrel//'b.rst'

!---- first call
  if( icall.eq.0 ) then

    ! checks dimensions and time
    open(25,file=rstname,status='old',form='unformatted',iostat=io)
    if( io.ne.0 ) stop 'rst_read: Error opening file'
    call rst_skip_record(25,atime,it,nvers,nrec,nkn0,nel0,nlv0,iflag,ierr)
    close(25)

    if( ( nkn0.ne.nkn ).or.( nel0.ne.nel ) ) stop "rst_read: dim bas error"
    if( ( nkn0.ne.nnkn ).or.( nel0.ne.nnel ).or.( nlv0.ne.nnlv ) ) stop "rst_read: dim ens error"

    ! init shyfem variables
    call mod_geom_dynamic_init(nkn0,nel0)
    call mod_hydro_init(nkn0,nel0,nlv0)
    call mod_hydro_vel_init(nkn0,nel0,nlv0)
    call mod_ts_init(nkn0,nlv0)  
    if( iconz_rst.gt.0 ) call mod_conz_init(iconz_rst,nkn0,nlv0)

    icall = icall + 1

  end if
!---- end first call

  open(24,file=trim(rstname),status='old',form='unformatted',iostat=io)
  if( io.ne.0 ) stop 'rst_read: Error opening file'

 89  call rst_read_record(atime,it,24,ierr)
     if( it.ne.nint(tt) ) goto 89

  close(24)

  if( it.ne.nint(tt) ) stop 'Error in rst file time'

  end subroutine rst_read

end module m_read_ensemble
