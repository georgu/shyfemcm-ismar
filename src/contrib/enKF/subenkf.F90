!********************************************************

  subroutine rst_read(nnkn,nnel,nnlv,rstname,tt)

  use basin
  use mod_restart
  use mod_geom_dynamic
  use levels, only : nlvdi,nlv
  use mod_hydro
  use mod_hydro_vel
  use mod_ts
  use mod_conz
  implicit none

  integer, intent(in) :: nnkn,nnel,nnlv
  character(len=*), intent(in) :: rstname
  double precision, intent(in) :: tt

  integer io
  integer it,nvers,nrec,iflag,ierr
  double precision atime
  integer date,time

  integer, save :: nkn0,nel0,nlv0
  integer, save :: icall = 0

!---- first call
  if( icall.eq.0 ) then

    ! checks dimensions and time
    open(25,file=rstname,status='old',form='unformatted',iostat=io)
    if( io.ne.0 ) stop 'rst_read: Error opening file'
    call rst_skip_record(25,atime,it,nvers,nrec,nkn0,nel0,nlv0,iflag,ierr)
    close(25)

    if( ( nkn0.ne.nkn ).or.( nel0.ne.nel ) ) stop "rst_read: dim bas error"
    if( ( nkn0.ne.nnkn ).or.( nel0.ne.nnel ).or.( nlv0.ne.nnlv ) ) stop "rst_read: dim ens error"

    ! setting the value of nlv
    nlv = nlv0

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

!********************************************************

  subroutine rst_write(rstname,tt)

  use mod_restart

  implicit none

  character(len=*), intent(in) :: rstname
  double precision, intent(in) :: tt

  integer it
  double precision ddate,dtime

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

!********************************************************

  subroutine find_node(x,y,ie,ik)
  use shyfile
  use basin
  implicit none

  real, intent(in) :: x,y
  integer, intent(in) :: ie
  integer, intent(out) :: ik
  
  integer i,ikmin
  real d,dmin

  dmin=10e10
  do i = 1,3
     ik = nen3v(i,ie)
     d = ( x - xgv(ik) )**2 + ( y - ygv(ik) )**2
     if( d.le.dmin ) then
       dmin = d
       ikmin = ik
     end if
  end do
  ik = ikmin

  end subroutine find_node

!********************************************************

  subroutine random2(work1,n)
! Returns a vector of random values N(variance=1,mean=0)
! From Evensen's code 
   implicit none
   integer, intent(in) :: n
   real,   intent(out) :: work1(n)
   real,   allocatable :: work2(:)
   real, parameter :: pi=3.14159253589

   allocate (work2(n))

   call random_number(work1)
   call random_number(work2)
   work1= sqrt(-2.0*log(work1))*cos(2.0*pi*work2)

   deallocate(work2)

  end subroutine random2

!********************************************************

  subroutine check_innov_val(inn,typeo)
  implicit none
  double precision, intent(inout) :: inn
  character(len=*), intent(in) :: typeo

  if( trim(typeo).eq.'level' ) then
    if ( abs(inn) .gt. 1. ) then
       write(*,*) 'Warning: level innovation too large: ',inn
       write(*,*) '         Setting it to zero...'
       inn = 0.
    end if
  else
    write(*,*) 'check_innov_val: warning observation not implemented'
  end if

  end subroutine check_innov_val

!********************************************************

  subroutine num2str(num,str)
  implicit none
  integer, intent(in) :: num
  character(len=3), intent(out) :: str

  if( (num.ge.0).and.(num.lt.10) ) then
    write(str,'(a2,i1)') '00',num
  elseif( (num.ge.10).and.(num.lt.100) ) then
    write(str,'(a1,i2)') '0',num
  elseif( (num.ge.100).and.(num.lt.1000) ) then
    write(str,'(i3)') num
  else
    stop 'num2str: num out of range'
  end if

  end subroutine num2str

!********************************************************

  subroutine make_pert(Ap,n)
  use mod_dimensions
  use mod_states
  implicit none

  integer, intent(in) :: n
  type(states), intent(out) :: Ap(n)

  
  end subroutine make_pert
