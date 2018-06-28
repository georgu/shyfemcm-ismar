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

  subroutine add_rst_params

  use mod_restart

  implicit none
  
  call addpar('ibarcl',0.)
  call addpar('iconz',0.)
  call addpar('ibfm',0.)
  call daddpar('date',0.)
  call daddpar('time',0.)

  end subroutine add_rst_params

!********************************************************

  subroutine rst_read(rstname,aatime)

  implicit none

  character(len=*), intent(in) :: rstname
  double precision, intent(in) :: aatime

  integer io
  integer ierr
  double precision atimef


  open(24,file=trim(rstname),status='old',form='unformatted',iostat=io)
  if (io /= 0) error stop 'rst_read: Error opening file'

 89  call rst_read_record(24,atimef,ierr)
     if (ierr /= 0) goto 90
     if (atimef /= aatime) goto 89

  close(24)

  return

 90 write(*,*) 'Error in the restart file. Are you sure that the analysis step'
    write(*,*) 'has a time present in the restart records?'
    error stop


  end subroutine rst_read

!********************************************************

  subroutine rst_write(rstname,aatime)

  use mod_restart
  use levels, only : hlv

  implicit none

  character(len=*), intent(in) :: rstname
  double precision, intent(in) :: aatime

  ! adds parameters
  !
  call putpar('ibarcl',float(ibarcl_rst))
  call putpar('iconz',float(iconz_rst))
  call putpar('ibfm',0.)
  !call dputpar('date',dfloat(date))
  !call dputpar('time',dfloat(time))

  ! In 2D barotropic hlv is set to 10000.
  !
  if (size(hlv) == 1) hlv = 10000.

  open(34,file=rstname,form='unformatted')
  call rst_write_record(aatime,34)
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

  ikmin = -999
  dmin=10e10
  do i = 1,3
     ik = nen3v(i,ie)
     d = ( x - xgv(ik) )**2 + ( y - ygv(ik) )**2
     if (d <= dmin) then
       dmin = d
       ikmin = ik
     end if
  end do
  ik = ikmin

  end subroutine find_node

!********************************************************

   subroutine layer_thick
   ! Set zenv from znv
   ! A value bigger than the depth hm3v is set
   use mod_dimensions
   use mod_hydro
   use basin

   implicit none
   integer ie,ii,k
   real*4 z,h
   real*4 hmin

   hmin = 0.03
   
   do ie = 1,nnel 
     do ii = 1,3
        k = nen3v(ii,ie)
        zenv(ii,ie) = znv(k)

        z = zenv(ii,ie)
        h = hm3v(ii,ie)
        if (z + h < hmin) then
           z = -h + hmin
           zenv(ii,ie) = z
        end if

     end do
   end do
   
   end subroutine layer_thick

!********************************************************

  subroutine num2str(num,str)

  implicit none

  integer, intent(in) :: num
  character(len=3), intent(out) :: str

  if ((num >= 0).and.(num < 10)) then
    write(str,'(a2,i1)') '00',num
  elseif ((num >= 10).and.(num < 100)) then
    write(str,'(a1,i2)') '0',num
  elseif ((num >= 100).and.(num < 1000)) then
    write(str,'(i3)') num
  else
    error stop 'num2str: num out of range'
  end if

  end subroutine num2str

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
     if (abs(aaux) >= 3.) then
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

  ! white noise if tau is lower than 0
  !
  if (tau <= 0.) return

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

  subroutine make_2Dpert(vec,n,nens,fmult,theta,nx,ny)

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

  !dx = xlength / float(nx - 1)
  !dy = ylength / float(ny - 1)
  dx = 0.25
  dy = 0.25

  rx = 5
  ry = 5
  !rx = rx/sqrt(3.0) !?
  !ry = ry/sqrt(3.0)

  verbose = .false.
  samp_fix = .true.	!keep true

  if (verbose) then
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

  end subroutine make_2Dpert

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
  if (ie == 0) then
     write(*,*) 'find_el_node: observations must be inside the grid.'
     write(*,*) 'x, y: ',x4,y4
     error stop
  end if

  dstmax = 1e15
  do ii = 1,3
     iik = nen3v(ii,ie)
     dst = sqrt( (xgv(iik)-x4)**2 + (ygv(iik)-y4)**2 )
     if (dst < dstmax) then
       dstmax = dst
       ik = iik
     end if
  end do

  end subroutine find_el_node

!********************************************************
  
  subroutine superobs_horiz_el(no,x,y,ostatus,val1,val2)

  use basin

  implicit none

  integer,intent(in) :: no
  real,intent(in) :: x(no),y(no)
  integer,intent(inout) :: ostatus(no)
  real,intent(inout) :: val1(no),val2(no)
  integer nobs(nel)
  integer, allocatable :: ieobs(:)
  integer, allocatable :: oindex(:,:)
  integer omax
  real*4 x4,y4
  integer n,ie,nn
  real avval

  ! find the index of the elements containing obs
  !
  allocate(ieobs(no))
  ieobs = -999
  nobs = 0
  do n = 1,no
     if (ostatus(n) > 0) cycle

     x4 = x(no)
     y4 = y(no)
     call find_element(x4,y4,ie)
     ieobs(n) = ie
     nobs(ie) = nobs(ie) + 1
  end do

  ! find the maximum number of obs inside an element
  !
  omax = 0
  do ie = 1,nel
     omax = max(omax,nobs(ie))
  end do

  ! find the final obs indexes
  !
  allocate(oindex(0:omax,nel))
  oindex = 0
  do n = 1,no
     if (ostatus(n) > 0) cycle

     ie = ieobs(n)
     nn = oindex(0,ie)
     nn = nn + 1
     oindex(nn,ie) = n
     oindex(0,ie) = nn
  end do

  ! average the values and assign the status
  !
  do ie = 1,nel
     nn = oindex(0,ie)
     if (nn == 0) cycle
     write(*,*) 'making super-observation...'
     avval = sum(val1(oindex(1:nn,ie)))/nn
     val1(oindex(1,ie)) = avval
     avval = sum(val2(oindex(1:nn,ie)))/nn
     val2(oindex(1,ie)) = avval
     ostatus(oindex(1,ie)) = 1
     ostatus(oindex(2:nn,ie)) = 2
  end do
  
  end subroutine superobs_horiz_el

!********************************************************

  subroutine save_X5(alabel,tt)

  implicit none
  character, intent(in) :: alabel*6
  double precision, intent(in) :: tt

  integer :: nrens
  integer :: nrobs
  real, allocatable :: X5(:,:)
  real, allocatable :: X3(:,:)
  real, allocatable :: S(:,:)
  character(len=2) :: tag2

  ! read X5.uf depending on the mult method (X5 or X3,S)
  open(unit=25,file='X5.uf',form='unformatted',status='old')
   read(25) tag2
   rewind 25
   if (tag2 == 'X3') then
      read(25) tag2,nrens,nrobs
      allocate(X3(nrens,nrens),S(nrobs,nrens))
      rewind 25
      read(25) tag2,nrens,nrobs,X3,S
   else
      read(25) tag2,nrens
      allocate(X5(nrens,nrens))
      rewind 25
      read(25) tag2,nrens,X5
   end if
  close(25)
  
  ! write the total X5_tot.uf adding time and type of analysis (global or local)
  write(*,*) 'saving X5 matrix for enKS'
  open(unit=35,file='X5_tot.uf',form='unformatted',status='unknown',access='append')
   write(35) tt,alabel,tag2
   if (tag2 == 'X3') then
      write(35) nrens,nrobs
      write(35) X3,S
   else
      write(35) nrens
      write(35) X5
   end if
  close(35)

  end subroutine save_X5
