module m_write_ensemble
! Writes the ensemble held in A to ensembleA.uf
contains

  subroutine write_ensemble
   use mod_enkf
   use mod_hydro
   use mod_ts
   use m_read_ensemble

   implicit none

   integer ne

   type(states4) A4

   do ne = 1,nrens
      call rst_read(ne,na,tobs) !This is to load var not present in the ens state. To be removed
      A4 = A(ne)
      utlnv = A4%u
      vtlnv = A4%v
      znv = A4%z
      tempv = A4%t
      saltv = A4%s
      call rst_write(ne,na,tobs)
   end do
  end subroutine write_ensemble


!********************************************************

  subroutine rst_write(ne,nan,tt)

  use mod_restart

  implicit none

  integer, intent(in) :: ne, nan
  double precision, intent(in) :: tt

  character(len=3) :: nrel,nal
  character(len=1) :: stype
  integer it
  character(len=16) rstname
  double precision ddate,dtime

  stype = 'a'

  if( ne.eq.-1 ) then
    write(*,*) 'Writing average background state...'
    nrel='avr'
    stype='b'
  elseif( ne.eq.-2 ) then
    write(*,*) 'Writing average analysis state...'
    nrel='avr'
  else
    call num2str(ne,nrel)
  end if

  call num2str(nan,nal)

  rstname = 'an'//nal//'_'//'en'//nrel//stype//'.rst'

  ! adds parameters
  call addpar('ibarcl',float(ibarcl_rst))
  call addpar('iconz',float(iconz_rst))
  call addpar('ibfm',0.)
  ddate = date_rst
  dtime = time_rst
  call daddpar('date',ddate)
  call daddpar('time',dtime)

  open(34,file=rstname,form='unformatted')
  it = nint(tt)
  call rst_write_record(it,34)
  close(34)

  end subroutine rst_write

end module m_write_ensemble

