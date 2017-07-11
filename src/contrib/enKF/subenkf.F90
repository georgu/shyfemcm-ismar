  subroutine init_shyfem_vars(nk,ne,nl)
  use basin
  use mod_restart
  use mod_geom_dynamic
  use levels
  use mod_hydro
  use mod_hydro_vel
  use mod_ts
  use mod_conz

  implicit none
  integer, intent(in) :: nk,ne,nl

  call mod_geom_dynamic_init(nk,ne)
  call mod_hydro_init(nk,ne,nl)
  call mod_hydro_vel_init(nk,ne,nl)
  call mod_ts_init(nk,nl)  
  !call mod_conz_init(iconz_rst,nk,nl)
  call levels_init(nk,ne,nl)

  nlvdi = nl
  nlv = nl
  end subroutine init_shyfem_vars

!********************************************************

  subroutine add_rst_pars
  use mod_restart
  implicit none
  
  call addpar('ibarcl',0.)
  call addpar('iconz',0.)
  call addpar('ibfm',0.)
  call daddpar('date',0.)
  call daddpar('time',0.)

  end subroutine add_rst_pars

!********************************************************

  subroutine rst_read(rstname,tt,date,time)

  implicit none

  character(len=*), intent(in) :: rstname
  double precision, intent(in) :: tt
  integer,intent(out) :: date,time

  integer io
  integer it,ierr
  double precision atime


  open(24,file=trim(rstname),status='old',form='unformatted',iostat=io)
  if (io /= 0) error stop 'rst_read: Error opening file'

 89  call rst_read_record(atime,it,24,ierr)
     if (it /= nint(tt)) goto 89
     if (ierr /= 0) goto 90

  close(24)

  call dts_from_abs_time(date,time,atime-tt)

  return

 90 write(*,*) 'Error in the restart file. Are you sure that the first observation'
    write(*,*) 'has the same time of the restart files?'
    error stop


  end subroutine rst_read

!********************************************************

  subroutine rst_write(rstname,tt,date,time)

  use mod_restart
  use levels, only : hlv

  implicit none

  character(len=*), intent(in) :: rstname
  double precision, intent(in) :: tt
  integer, intent(in) :: date,time

  integer it

  ! adds parameters
  !
  call putpar('ibarcl',float(ibarcl_rst))
  call putpar('iconz',float(iconz_rst))
  call putpar('ibfm',0.)
  call dputpar('date',dfloat(date))
  call dputpar('time',dfloat(time))

  ! In 2D barotropic hlv is set to 10000.
  !
  if (size(hlv) == 1) hlv = 10000.

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
    error stop 'num2str: num out of range'
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

!********************************************************

  subroutine find_el_node(x,y,ie,ik)
  use basin
  implicit none

  real, intent(in) :: x,y
  integer, intent(out) :: ie,ik
  real*4 x4,y4
  real dst,dstmax
  integer iik,ii

  !-----------
  ! Finds the grid element and the node 
  ! nearest to the observation (the subroutine is in real4)
  !-----------
  x4 = x
  y4 = y
  call find_element(x4,y4,ie)
  if( ie.eq.0 ) then
     write(*,*) 'find_el_node: observations must be inside the grid.'
     write(*,*) 'x, y: ',x4,y4
     error stop
  end if

  dstmax = 1e15
  do ii = 1,3
     iik = nen3v(ii,ie)
     dst = sqrt( (xgv(iik)-x4)**2 + (ygv(iik)-y4)**2 )
     if( dst.lt.dstmax ) then
       dstmax = dst
       ik = iik
     end if
  end do

  end subroutine find_el_node

!********************************************************

  subroutine make_0Dpert(vflag,n,na,id,vec,t,tau)
  implicit none
  character(len=1), intent(in) :: vflag
  integer, intent(in) :: n,na,id
  real, intent(out) :: vec(n)
  double precision, intent(in) :: t,tau
  character(len=3) :: nal,idl
  character(len=17) :: pfile
  logical bfile
  integer nf
  real, allocatable :: vec_old(:)
  double precision t_old
  double precision alpha,dt

  ! make a new perturbation
  !
  call random_vec(vec,n)

  ! if exist load an old perturbation and merge
  !
  call num2str(na-1,nal)
  call num2str(id,idl)
  pfile = vflag // 'pert_' // nal // '_' // idl // '.bin'
  inquire(file=pfile,exist=bfile)
  if (bfile) then
     open(22,file=pfile,status='old',form='unformatted')
     read(22) nf
     if (nf /= n) error stop 'make_level_pert: dimension mismatch'
     read(22) t_old
     allocate(vec_old(nf))
     read(22) vec_old
     close(22)

     dt = t - t_old
     alpha = 1. - (dt/tau) 
     vec = alpha * vec_old + sqrt(1 - alpha**2) * vec
  else
     !write(*,*) 'No old file with perturbations: ',pfile
     continue
  end if

  ! save the last perturbation
  !
  call num2str(na,nal)
  pfile = vflag // 'pert_' // nal // '_' // idl // '.bin'
  open(32,file=pfile,form='unformatted')
  write(32) n
  write(32) t
  write(32) vec
  close(32)

  end subroutine make_0Dpert

!********************************************************

  subroutine random_vec(v,vdim)
  use m_random
  implicit none
  integer vdim
  real v(vdim),vaux(vdim-1)
  real aaux,ave
  integer n

  call random(vaux,vdim-1)
  ! remove outlayers
  do n = 1,vdim-1
     aaux = vaux(n)
     if( abs(aaux).ge.3. ) then
       aaux = aaux/abs(aaux) * (abs(aaux)-floor(abs(aaux)) + 1.) 
     end if
     vaux(n) = aaux
  end do

  ! set mean eq to zero
  ave = sum(vaux)/float(vdim-1)
  vaux = vaux - ave

  v(1) = 0.
  v(2:vdim) = vaux

  end subroutine random_vec
