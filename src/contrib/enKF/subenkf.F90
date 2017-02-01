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
  integer it,ierr
  double precision atime

  nlv = nnlv
 
  ! init shyfem variables
  call mod_geom_dynamic_init(nnkn,nnel)
  call mod_hydro_init(nnkn,nnel,nnlv)
  call mod_hydro_vel_init(nnkn,nnel,nnlv)
  call mod_ts_init(nnkn,nnlv)  
  !call mod_conz_init(iconz_rst,nnkn,nnlv)

  open(24,file=trim(rstname),status='old',form='unformatted',iostat=io)
  if( io.ne.0 ) stop 'rst_read: Error opening file'

 89  call rst_read_record(atime,it,24,ierr)
     if( ierr.ne.0 ) goto 90
     if( it.ne.nint(tt) ) goto 89

  close(24)

  return

 90 stop 'Error in the restart file'

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

  subroutine make_pert(vec,n,nens,fmult,theta,nx,ny)
  use basin
  use m_sample2D
  implicit none

  integer, intent(in) :: n,nens
  real, intent(out) :: vec(n,nens)
  integer, intent(in) :: fmult	!Mult factor for the super-sampling
  real, intent(in) :: theta	!Rotation of the random fields (theta=0 is east, rotation anticlocwise)
  integer, intent(in) :: nx,ny	!number of grid points in x and y direction for the regular grid

  real x1,x2,y1,y2,x0,y0,xlength,ylength
  real dx,dy,rx,ry
  real*4 sdx,sdy,sx0,sy0
  logical samp_fix,verbose

  real, allocatable :: mat(:,:,:)
  real*4, allocatable :: mat4(:,:),vec4fem(:)

  integer ne

  !----------------------------------------------------
  ! Computes some geometric parameters from the grid
  !----------------------------------------------------
  x1 = minval(xgv)
  x2 = maxval(xgv)
  y1 = minval(ygv)
  y2 = maxval(ygv)
  x0 = floor(x1)
  y0 = floor(y1)
  xlength = ceiling(x2) - floor(x1)
  ylength = ceiling(y2) - floor(y1)

  dx = xlength / float(nx - 1)
  dy = ylength / float(ny - 1)

  rx = 100. * dx
  ry = 100. * dy
  rx = rx/sqrt(3.0) !?
  ry = ry/sqrt(3.0)

  verbose = .false.
  samp_fix = .true.	!keep true

  if( verbose ) then
    write(*,'(a20,2f8.4,i5,f8.4)') 'x0,xlength,nx,dx: ',x0,xlength,nx,dx
    write(*,'(a20,2f8.4,i5,f8.4)') 'y0,ylength,ny,dy: ',y0,ylength,ny,dy
    write(*,'(a14,2f10.4,1x,f5.1)') 'rx,ry,theta: ',rx,ry,theta
  end if

  !----------------------------------------------------
  ! creates the sample
  !----------------------------------------------------
  allocate(mat(nx,ny,nens))
  call sample2D(mat,nx,ny,nens,fmult,dx,dy,rx,ry,theta,samp_fix,verbose)

  !----------------------------------------------------
  ! Interpolates over the FEM grid
  !----------------------------------------------------
  write(*,*) 'Interpolating 2D field over the FEM grid...'
  sdx = dx
  sdy = dy
  sx0 = x0
  sy0 = y0
  call setgeo(sx0,sy0,sdx,sdy,-999.)

  allocate(mat4(nx,ny),vec4fem(n))
  do ne = 1,nens
    mat4 = mat(:,:,ne)
    call am2av(mat4,vec4fem,nx,ny)
    vec(:,ne) = vec4fem
  end do
  deallocate(mat,mat4,vec4fem)

  end subroutine make_pert


