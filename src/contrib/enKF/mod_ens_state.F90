!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
!
module mod_ens_state

  use mod_dimensions
  use mod_mod_states
  use mod_para
  use mod_init_enkf

  implicit none

  type(states), dimension(:), allocatable :: Abk,Aan      ! ensemble states
  type(states) :: Abk_m, Aan_m                    ! mean state
  type(states) :: Abk_std, Aan_std     ! standard deviation old and new

contains

!********************************************************

  subroutine read_ensemble

   use basin
   use levels
   use mod_geom_dynamic
   use mod_hydro
   use mod_hydro_vel
   use mod_ts
   use mod_conz
   use mod_gotm_aux

   use mod_restart

   implicit none

   character(len=5) :: nrel,nal
   character(len=80) rstname
   integer ne

   ! read basin and check dimensions
   open(21, file=basfile, status='old', form='unformatted')
   call basin_read_by_unit(21)
   close(21)
   if ((nkn /= nnkn).or.(nel /= nnel)) error stop "read_basin: dim error"

   ! set vertical levels
   nlv = nnlv
   nlvdi = nnlv

   ! init some shyfem vars
   call mod_geom_dynamic_init(nkn,nel)
   call mod_hydro_init(nkn,nel,nlv)
   call mod_hydro_vel_init(nkn,nel,nlv)
   call mod_ts_init(nkn,nlv)  
   call levels_init(nkn,nel,nlv)
   call mod_gotm_aux_init(nkn,nlv)
   ! init concentration, this is a issue
   !call mod_conz_init(1,nkn,nlvdi)


   ! Allocates the state Abk to store the ens states
   if (.not. allocated(Abk)) allocate(Abk(nrens))

   ! init to zero
   do ne = 1,nrens
      Abk(ne) = 0.
   end do

   call num2str(nanal,nal)

   if ((bnew_ens == 0) .or. (nanal.gt.1)) then

     write(*,*) 'loading an ensemble of initial states...'
     do ne = 1,nrens
        call num2str(ne-1,nrel)
        rstname='an'//nal//'_'//'en'//nrel//'b.rst'
        call read_state(Abk(ne),rstname)
     end do

   else if ((bnew_ens == 1) .and. (nanal == 1)) then

     write(*,*) 'creating a new ensemble of initial states...'
     !read an input restart file
     call num2str(0,nrel)
     rstname = 'an'//nal//'_'//'en'//nrel//'b.rst'
     call read_state(Abk(1),rstname)

     call make_init_ens(Abk(1))
     
     !save the initial ens in new restart files
     call num2str(nanal,nal)
     do ne = 1,nrens
        call num2str(ne-1,nrel)
        rstname='an'//nal//'_'//'en'//nrel//'b.rst'
        call write_state(Abk(ne),rstname)
     end do

   else

     write(*,*) 'not a valid option for bnew_ens'
     error stop

   end if

   return
  end subroutine read_ensemble

!********************************************************

  subroutine write_ensemble
   implicit none

   character(len=5) :: nrel,nal
   character(len=80) rstname
   integer ne

   write(*,*) 'writing the restart files...'
   call num2str(nanal,nal)
   do ne = 1,nrens
      call num2str(ne-1,nrel)
      rstname='an'//nal//'_'//'en'//nrel//'a.rst'
      call write_state(Aan(ne),rstname)
   end do

   rstname='an'//nal//'_mean_a.rst'
   call write_state(Aan_m,rstname)
   rstname='an'//nal//'_std_a.rst'
   call write_state(Aan_std,rstname)

  end subroutine write_ensemble

!********************************************************

  subroutine make_init_ens(Ain)
   use basin
   use mod_restart
   use mod_para

   implicit none

   type(states),intent(in) :: Ain

   real kvec1(nnkn,nrens-1),kvec2(nnkn,nrens-1)

   integer ne,nl

   ! this id for the level
   call make_2Dpert(kvec1,nnkn,nrens-1)

   Abk(1) = Ain
   do ne = 2,nrens
      Abk(ne) = Ain
      Abk(ne)%z = Ain%z + (kvec1(:,ne-1) * sigma_init_z)
   end do

   if (ibarcl_rst /= 0) then
      call make_2Dpert(kvec1,nnkn,nrens-1)
      call make_2Dpert(kvec2,nnkn,nrens-1)

      do ne = 2,nrens
        do nl = 1,nnlv
         Abk(ne)%t(nl,:) = Ain%t(nl,:) + (kvec1(:,ne-1) * sigma_init_t)
         Abk(ne)%s(nl,:) = Ain%s(nl,:) + (kvec2(:,ne-1) * sigma_init_s)
	end do
      end do
   end if

  end subroutine make_init_ens

!********************************************************

  subroutine bc_val_check_correct(Abkg,Aanl)
   implicit none
   type(states),intent(in) :: Abkg(nrens)
   type(states),intent(inout) :: Aanl(nrens)

   real vala,valb

   integer ne,k,ie,nl

   integer zout,uvout,sout,tout,otot
   integer znan,uvnan,snan,tnan,ntot

   character(len=80),parameter :: bcfile = 'lbound.dat'
   logical :: file_exists
   integer :: nbc
   integer,allocatable :: bcid(:)
   real,allocatable :: bcrho(:)
   real :: w
   
   nbc = 1
   w = 0.
   inquire(file=bcfile,exist=file_exists)
   if (file_exists) then
      ! read file
      write(*,*) 'File to correct values near the boundaries found.'
      allocate(bcid(nbc),bcrho(nbc))
      call read_bc_file(0,bcfile,nbc,bcid,bcrho)
      deallocate(bcid,bcrho)
      allocate(bcid(nbc),bcrho(nbc))
      call read_bc_file(1,bcfile,nbc,bcid,bcrho)
   else
      write(*,*) 'Warning: no file to correct values near the boundaries.'
      write(*,*) 'Make a file lbound.dat. In the first line the n. of BC nodes.'
      write(*,*) 'In the other lines, in each line, the external number and a damping radius.'
   end if

   otot = 0
   ntot = 0
   do ne = 1,nrens

    zout = 0
    uvout = 0
    sout = 0
    tout = 0

    znan = 0
    uvnan = 0
    snan = 0
    tnan = 0

!$OMP PARALLEL PRIVATE(k,nl,valb,vala,znan,zout,snan,sout,tnan,tout,w),SHARED(Abkg,Aanl,nbc,bcid,bcrho)
!$OMP DO
    do k = 1,nnkn

      ! BC correction
      if (file_exists) call bc_correction('node',k,nbc,bcid,bcrho,w)

      ! level
      valb = Abkg(ne)%z(k)
      vala = Aanl(ne)%z(k)
      vala = w * valb + (1. - w) * vala

      ! check val
      call check_one_val(vala,valb,SSH_MAX,SSH_MIN,znan,zout)
      Aanl(ne)%z(k) = vala

      do nl = 1,nnlv
        ! salinity
	valb = Abkg(ne)%s(nl,k)
	vala = Aanl(ne)%s(nl,k)
        vala = w * valb + (1. - w) * vala
        ! check val
	call check_one_val(vala,valb,SAL_MAX,SAL_MIN,snan,sout)
	Aanl(ne)%s(nl,k) = vala

	! temperature
	valb = Abkg(ne)%t(nl,k)
	vala = Aanl(ne)%t(nl,k)
        vala = w * valb + (1. - w) * vala
        ! check val
	call check_one_val(vala,valb,TEM_MAX,TEM_MIN,tnan,tout)
	Aanl(ne)%t(nl,k) = vala
      end do
    end do
!$OMP ENDDO
!$OMP END PARALLEL


!$OMP PARALLEL PRIVATE(ie,nl,valb,vala,uvnan,uvout,w),SHARED(Abkg,Aanl,nbc,bcid,bcrho)
!$OMP DO
    do ie = 1,nnel
      ! BC correction
      if (file_exists) call bc_correction('elem',ie,nbc,bcid,bcrho,w)
      do nl = 1,nnlv
        ! current
	valb = Abkg(ne)%u(nl,ie)
	vala = Aanl(ne)%u(nl,ie)
        vala = w * valb + (1. - w) * vala
        ! check val
	call check_one_val(vala,valb,VEL_MAX,VEL_MIN,uvnan,uvout)
	Aanl(ne)%u(nl,ie) = vala

	valb = Abkg(ne)%v(nl,ie)
	vala = Aanl(ne)%v(nl,ie)
        vala = w * valb + (1. - w) * vala
        ! check val
	call check_one_val(vala,valb,VEL_MAX,VEL_MIN,uvnan,uvout)
	Aanl(ne)%v(nl,ie) = vala
      end do
    end do
!$OMP ENDDO
!$OMP END PARALLEL

    otot = zout + tout + sout + uvout
    ntot = znan + tnan + snan + uvnan

    if ((verbose).and.(otot > 0)) write(*,*) 'Ensemble member: ',ne
    if ((verbose).and.(zout > 0)) write(*,*) 'Number of levels out of range: ', zout
    if ((verbose).and.(tout > 0)) write(*,*) 'Number of temperatures out of range: ', tout
    if ((verbose).and.(sout > 0)) write(*,*) 'Number of salinities out of range: ', sout
    if ((verbose).and.(uvout > 0)) write(*,*) 'Number of velocities out of range: ', uvout

    if ((verbose).and.(ntot > 0)) write(*,*) 'Ensemble member: ',ne
    if ((verbose).and.(znan > 0)) write(*,*) 'Number of nan levels: ', znan
    if ((verbose).and.(tnan > 0)) write(*,*) 'Number of nan temperatures: ', tnan
    if ((verbose).and.(snan > 0)) write(*,*) 'Number of nan salinities: ', snan
    if ((verbose).and.(uvnan > 0)) write(*,*) 'Number of nan velocities: ', uvnan

   end do

   if (otot > 0) write(*,*) 'Number of variables out of range: ', otot
   if (ntot > 0) write(*,*) 'Number of nan variables: ', ntot

  end subroutine bc_val_check_correct

!********************************************************
  subroutine read_bc_file(icall,bcfile,nbc,bcid,bcrho)
  implicit none
  character(len=*) :: bcfile
  integer,intent(in) :: icall
  integer,intent(inout) :: nbc
  integer,intent(out) :: bcid(nbc)
  real,intent(out) :: bcrho(nbc)
  integer i

  if (icall == 0) then
     open(28,file=trim(bcfile),status='old',form='formatted')
     read(28,*) nbc
     return
  end if

  do i=1,nbc
     read(28,*) bcid(i),bcrho(i)
  end do
  close(28)

  end subroutine read_bc_file

!********************************************************
  subroutine bc_correction(stype,id,nbc,bcid,bcrho,w)
  use basin
  implicit none
  character(len=4),intent(in) :: stype
  integer,intent(in) :: id
  integer,intent(in) :: nbc
  integer,intent(in) :: bcid(nbc)
  real,intent(in) :: bcrho(nbc)
  real,intent(out) :: w
  integer :: i,k,kbc
  real :: bcx,bcy
  real :: x,y,d,dmin,rho
  integer :: ipint

  x = 0.
  y = 0.
  d = 1.e15
  dmin = 1.e15
  rho = 0.
  if (stype == 'node') then
     x = xgv(id)
     y = ygv(id)
  else
     do i = 1,3
        k = nen3v(i,id)
        x = x + xgv(k)
        y = y + ygv(k)
     end do
     x = x/3.
     y = y/3.
  end if

  do i = 1,nbc
     kbc = ipint(bcid(i))
     bcx = xgv(kbc)
     bcy = ygv(kbc)
     d = sqrt((x-bcx)**2 + (y-bcy)**2)
     if (d < dmin) then
	dmin = d
	rho = bcrho(i)
     end if
  end do

  call find_weight_GC(rho,dmin,w)


  !write(*,*) 'bc_correction: ',rho,dmin,w

  end subroutine bc_correction


!********************************************************
  subroutine check_one_val(va,vb,vmax,vmin,vnan,vout)
  implicit none
  real,intent(inout) :: va
  real,intent(in) :: vb
  real,intent(in) :: vmin,vmax
  integer,intent(inout) :: vnan,vout

  ! isnan and bk is a number
  if ( isnan(va) .and. (.not. isnan(vb) ) ) then
     va = vb
!$OMP CRITICAL
     vnan = vnan + 1
!$OMP END CRITICAL
   end if

   ! out of range
   if (( va >= vmax ) .or. ( va <= vmin )) then
      va = vb
!$OMP CRITICAL
      vout = vout + 1
!$OMP END CRITICAL
   end if

  end subroutine check_one_val

!********************************************************

  subroutine make_mean_std(tflag)

  implicit none
  character(len=1), intent(in) :: tflag

  if (tflag == 'a') then
     call mean_state(nrens,Aan,Aan_m)
     call std_state(nrens,Aan,Aan_std)
  else
     call mean_state(nrens,Abk,Abk_m)
     call std_state(nrens,Abk,Abk_std)
  end if

  end subroutine make_mean_std

!********************************************************

  subroutine write_state(Astate,filename)
  use mod_hydro
  use mod_ts
  implicit none
  type(states),intent(in) :: Astate
  character(len=*),intent(in) :: filename

  type(states4),allocatable :: A4

  allocate(A4)
  call states8to4(A4,Astate)
  call pull_state(A4)
  call rst_write(trim(filename),atime_an)
  deallocate(A4)

  end subroutine write_state

!********************************************************

  subroutine read_state(Astate,filename)

  implicit none

  type(states),intent(out) :: Astate
  character(len=*),intent(in) :: filename

  type(states4),allocatable :: A4

  allocate(A4)
  Astate = 0.
  call states8to4(A4,Astate)
  call rst_read(filename,atime_an)
  call push_state(A4)
  call states4to8(Astate,A4)
  deallocate(A4)

  end subroutine read_state

!********************************************************

   subroutine push_state(A4)
    use mod_hydro
    use mod_hydro_vel
    use mod_ts
    use mod_conz
    use mod_restart
    implicit none
 
    type(states4),intent(inout) :: A4

    ! no significant differences by using currents rather than transports
    A4%u = utlnv
    A4%v = vtlnv
    A4%z = znv
    if (ibarcl_rst /= 0) then
       A4%t = tempv
       A4%s = saltv
    else
       A4%t = 0.
       A4%s = 0.
    end if
   
   end subroutine push_state

!********************************************************

   subroutine pull_state(A4)
    use mod_hydro
    use mod_hydro_vel
    use mod_ts
    use mod_conz
    use mod_restart
    implicit none
 
    type(states4),intent(in) :: A4

    ! no significant differences by using currents rather than transports
    utlnv = A4%u
    vtlnv = A4%v
    znv = A4%z
    if (ibarcl_rst /= 0) then
       tempv = A4%t
       saltv = A4%s
    else
       tempv = 0.
       saltv = 0.
    end if

    ! make zenv
    call layer_thick
   end subroutine pull_state

end module mod_ens_state
