!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
!
module mod_ens_state

  use mod_init_enkf
  use mod_mod_states
  use mod_para

  implicit none

  type(states), dimension(:), allocatable :: Abk,Aan      ! ensemble states
  type(states) :: Abk_m, Aan_m                    ! mean state
  type(states) :: Abk_std, Aan_std     ! standard deviation old and new

contains


!********************************************************

  subroutine allocate_all
  implicit none
  integer ne

   ! Allocates the state Abk to store the ens states
   if (.not. allocated(Abk)) allocate(Abk(nrens))
   if (.not. allocated(Aan)) allocate(Aan(nrens))

   ! allocate and init to zero
   call allocate_states(Abk_m,nnkn,nnel,nnlv)
   call allocate_states(Aan_m,nnkn,nnel,nnlv)
   call allocate_states(Abk_std,nnkn,nnel,nnlv)
   call allocate_states(Aan_std,nnkn,nnel,nnlv)
   do ne = 1,nrens
      call allocate_states(Abk(ne),nnkn,nnel,nnlv)
      call allocate_states(Aan(ne),nnkn,nnel,nnlv)
   end do

  end subroutine allocate_all


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
   use shympi

   use mod_restart

   implicit none

   character(len=5) :: nrel,nal
   character(len=80) rstname
   integer ne

   ! read basin and check dimensions
   open(21, file=basfile, status='old', form='unformatted')
   call basin_read_by_unit(21)
   close(21)

   ! set dimensions
   nnkn = nkn
   nnel = nel
   nlv = nnlv
   !nlvdi = nnlv

   ! init some shyfem vars
   call mod_geom_dynamic_init(nkn,nel)
   call mod_hydro_init(nkn,nel,nlv)
   call mod_hydro_vel_init(nkn,nel,nlv)
   call mod_ts_init(nkn,nlv)  
   call levels_init(nkn,nel,nlv)
   call mod_gotm_aux_init(nkn,nlv)
   call shympi_set_hlv(nlv,hlv)
   call shympi_init(.false.)
   ! init concentration, this is a issue
   !call mod_conz_init(1,nkn,nlvdi)

   call allocate_all

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

  subroutine bc_val_check_correct
   implicit none

   real vala,valb,vala_old

   integer ne,k,ie,nl

   integer zout,uvout,sout,tout,otot
   integer znan,uvnan,snan,tnan,ntot
   integer zbig,uvbig,sbig,tbig,btot

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
   btot = 0
   do ne = 1,nrens

    zout = 0
    uvout = 0
    sout = 0
    tout = 0

    znan = 0
    uvnan = 0
    snan = 0
    tnan = 0

    zbig = 0
    uvbig = 0
    sbig = 0
    tbig = 0


!$OMP PARALLEL PRIVATE(k,nl,valb,vala,znan,zout,snan,sout,tnan,tout,w),SHARED(ne,Abk,Aan,nbc,bcid,bcrho)
!$OMP DO
    do k = 1,nnkn

      ! BC correction
      if (file_exists) call bc_correction('node',k,nbc,bcid,bcrho,w)

      ! level
      ! Correct only with the non perturbed value
      valb = Abk(1)%z(k)
      !valb = Abk(ne)%z(k)
      vala = Aan(ne)%z(k)
      vala = w * valb + (1. - w) * vala

      ! check val
      call check_one_val(vala,valb,SSH_MAX,SSH_MIN,znan,zout,zbig)
      Aan(ne)%z(k) = vala

      do nl = 1,nnlv
        ! salinity
	valb = Abk(ne)%s(nl,k)
	vala = Aan(ne)%s(nl,k)
	vala = w * valb + (1. - w) * vala
        ! check val
	call check_one_val(vala,valb,SAL_MAX,SAL_MIN,snan,sout,sbig)
	Aan(ne)%s(nl,k) = vala

	! temperature
	valb = Abk(ne)%t(nl,k)
	vala = Aan(ne)%t(nl,k)
        vala = w * valb + (1. - w) * vala
        ! check val
	call check_one_val(vala,valb,TEM_MAX,TEM_MIN,tnan,tout,tbig)
	Aan(ne)%t(nl,k) = vala
      end do
    end do
!$OMP ENDDO
!$OMP END PARALLEL


!$OMP PARALLEL PRIVATE(ie,nl,valb,vala,uvnan,uvout,w),SHARED(ne,Abk,Aan,nbc,bcid,bcrho)
!$OMP DO
    do ie = 1,nnel
      ! BC correction
      if (file_exists) call bc_correction('elem',ie,nbc,bcid,bcrho,w)
      do nl = 1,nnlv
        ! Current, check only the speed
	! BC correction
        Aan(ne)%u(nl,ie) = w * Abk(ne)%u(nl,ie) + (1. - w) * Aan(ne)%u(nl,ie)
        Aan(ne)%v(nl,ie) = w * Abk(ne)%v(nl,ie) + (1. - w) * Aan(ne)%v(nl,ie)

        ! Increment correction
	valb = sqrt(Abk(ne)%u(nl,ie)**2 + Abk(ne)%v(nl,ie)**2)
	vala = sqrt(Aan(ne)%u(nl,ie)**2 + Aan(ne)%v(nl,ie)**2)
	vala_old = vala
	call check_one_val(vala,valb,VEL_MAX,-0.001,uvnan,uvout,uvbig)

	Aan(ne)%u(nl,ie) = Aan(ne)%u(nl,ie)*(vala/vala_old)
	Aan(ne)%v(nl,ie) = Aan(ne)%v(nl,ie)*(vala/vala_old)

      end do
    end do
!$OMP ENDDO
!$OMP END PARALLEL

    otot = zout + tout + sout + uvout
    ntot = znan + tnan + snan + uvnan
    btot = zbig + tbig + sbig + uvbig

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

    if ((verbose).and.(btot > 0)) write(*,*) 'Ensemble member: ',ne
    if ((verbose).and.(zbig > 0)) write(*,*) 'Number of level increments too high: ', zbig
    if ((verbose).and.(tbig > 0)) write(*,*) 'Number of temperature increments too high: ', tbig
    if ((verbose).and.(sbig > 0)) write(*,*) 'Number of salinity increments too high: ', sbig
    if ((verbose).and.(uvbig > 0)) write(*,*) 'Number of velocity increments too high: ', uvbig

   end do

   if (otot > 0) write(*,*) 'Number of variables out of range: ', otot
   if (ntot > 0) write(*,*) 'Number of nan variables: ', ntot
   if (btot > 0) write(*,*) 'Number of increments too high: ', btot

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
  subroutine check_one_val(va,vb,vmax,vmin,vnan,vout,vbig)
  implicit none
  real,intent(inout) :: va
  real,intent(in) :: vb
  real,intent(in) :: vmin,vmax
  integer,intent(inout) :: vnan,vout,vbig
  real, parameter :: max_inc = 0.8
  real rel_inc,inc

   ! isnan and bk is a number
   if ( isnan(va) .and. (.not. isnan(vb) ) ) then
     va = vb
!$OMP CRITICAL
     vnan = vnan + 1
!$OMP END CRITICAL
   end if

   ! set a maximum increment
   inc = va-vb
   rel_inc = abs(inc)/abs(vb)
   if (rel_inc > max_inc) then
      va = vb + inc * (max_inc/rel_inc)
      !write(*,*) 'Increment correction: ',max_inc,rel_inc,inc,vb,va,va-vb
!$OMP CRITICAL
      vbig = vbig + 1
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
     call std_state(nrens,nnkn,nnel,nnlv,Aan,Aan_std)
  else
     call mean_state(nrens,Abk,Abk_m)
     call std_state(nrens,nnkn,nnel,nnlv,Abk,Abk_std)
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
  call allocate_states4(A4,nnkn,nnel,nnlv)

  call states8to4(A4,Astate)
  call pull_state(A4)
  call rst_write(trim(filename),atime_an)
  deallocate(A4)

  end subroutine write_state

!********************************************************

  subroutine read_state(Astate,filename)

  implicit none

  type(states),intent(inout) :: Astate
  character(len=*),intent(in) :: filename

  type(states4),allocatable :: A4

  allocate(A4)
  call allocate_states4(A4,nnkn,nnel,nnlv)
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
    call layer_thick(nnel)
   end subroutine pull_state

!-------------------------------------------------------------------
! subroutines to switch between type and matrix formats
!-------------------------------------------------------------------

   subroutine tystate_to_matrix(ibrcl,nens,ndim,A,Amat)
      implicit none
      integer, intent(in) :: ibrcl
      integer, intent(in) :: nens
      integer, intent(in) :: ndim
      type(states), intent(in) :: A(nens)
      real, intent(out) :: Amat(ndim,nens)

      integer i,dimuv,dimts,dimz

      dimz = nnkn
      dimuv = nnlv*nnel
      dimts = nnlv*nnkn
      do i = 1,nens
         Amat(1:dimuv,i) = reshape(A(i)%u,(/dimuv/))
         Amat(dimuv+1:2*dimuv,i) = reshape(A(i)%v,(/dimuv/))
         Amat(2*dimuv+1:2*dimuv+dimz,i) = A(i)%z
      end do
      if (ibrcl > 0) then
         do i = 1,nens
            Amat(2*dimuv+dimz+1:2*dimuv+dimz+dimts,i) = reshape(A(i)%t,(/dimts/))
            Amat(2*dimuv+dimz+dimts+1:2*dimuv+dimz+2*dimts,i) = reshape(A(i)%s,(/dimts/))
	 end do
      end if
   end subroutine tystate_to_matrix

   subroutine matrix_to_tystate(ibrcl,nens,ndim,Amat,A)
      implicit none
      integer, intent(in) :: ibrcl
      integer, intent(in) :: nens
      integer, intent(in) :: ndim
      real, intent(in) :: Amat(ndim,nens)
      type(states), intent(inout) :: A(nens)

      integer i,dimuv,dimts,dimz

      dimz = nnkn
      dimuv = nnlv*nnel
      dimts = nnlv*nnkn
      do i = 1,nens
         A(i)%u = reshape(Amat(1:dimuv,i),(/nnlv,nnel/))
         A(i)%v = reshape(Amat(dimuv+1:2*dimuv,i),(/nnlv,nnel/))
         A(i)%z = Amat(2*dimuv+1:2*dimuv+dimz,i)
      end do
      if (ibrcl > 0) then
         do i = 1,nens
            A(i)%t = reshape(Amat(2*dimuv+dimz+1:2*dimuv+dimz+dimts,i),(/nnlv,nnkn/))
            A(i)%s = reshape(Amat(2*dimuv+dimz+dimts+1:2*dimuv+dimz+2*dimts,i),(/nnlv,nnkn/))
	 end do
      end if
   end subroutine matrix_to_tystate

   subroutine tyqstate_to_matrix(ibrcl,nens,ndim,A,Amat)
      implicit none
      integer, intent(in) :: ibrcl
      integer, intent(in) :: nens
      integer, intent(in) :: ndim
      type(qstates), intent(in) :: A(nens)
      real, intent(out) :: Amat(2*ndim,nens)

      integer i,dimuv,dimts,dimz

      dimz = nnkn
      dimuv = nnlv*nnel
      dimts = nnlv*nnkn
      do i = 1,nens
         Amat(1:dimuv,i) = reshape(A(i)%qu,(/dimuv/))
         Amat(dimuv+1:2*dimuv,i) = reshape(A(i)%qv,(/dimuv/))
         Amat(2*dimuv+1:2*dimuv+dimz,i) = A(i)%qz
      end do
      if (ibrcl > 0) then
         do i = 1,nens
            Amat(2*dimuv+dimz+1:2*dimuv+dimz+dimts,i) = reshape(A(i)%qt,(/dimts/))
            Amat(2*dimuv+dimz+dimts+1:2*dimuv+dimz+2*dimts,i) = reshape(A(i)%qs,(/dimts/))
	 end do
      end if

      do i = 1,nens
         Amat(ndim+1:ndim+dimuv,i) = reshape(A(i)%u,(/dimuv/))
         Amat(ndim+dimuv+1:ndim+2*dimuv,i) = reshape(A(i)%v,(/dimuv/))
         Amat(ndim+2*dimuv+1:ndim+2*dimuv+dimz,i) = A(i)%z
      end do

      if (ibrcl > 0) then
         do i = 1,nens
            Amat(ndim+2*dimuv+dimz+1:ndim+2*dimuv+dimz+dimts,i) = reshape(A(i)%t,(/dimts/))
            Amat(ndim+2*dimuv+dimz+dimts+1:ndim+2*dimuv+dimz+2*dimts,i) = reshape(A(i)%s,(/dimts/))
	 end do
      end if
   end subroutine tyqstate_to_matrix

   subroutine matrix_to_tyqstate(ibrcl,nens,ndim,Amat,A)
      implicit none
      integer, intent(in) :: ibrcl
      integer, intent(in) :: nens
      integer, intent(in) :: ndim
      real, intent(in) :: Amat(2*ndim,nens)
      type(qstates), intent(out) :: A(nens)

      integer i,dimuv,dimts,dimz

      dimz = nnkn
      dimuv = nnlv*nnel
      dimts = nnlv*nnkn
      do i = 1,nens
         A(i)%qu = reshape(Amat(1:dimuv,i),(/nnlv,nnel/))
         A(i)%qv = reshape(Amat(dimuv+1:2*dimuv,i),(/nnlv,nnel/))
         A(i)%qz = Amat(2*dimuv+1:2*dimuv+dimz,i)
      end do
      if (ibrcl > 0) then
         do i = 1,nens
            A(i)%qt = reshape(Amat(2*dimuv+dimz+1:2*dimuv+dimz+dimts,i),(/nnlv,nnkn/))
            A(i)%qs = reshape(Amat(2*dimuv+dimz+dimts+1:2*dimuv+dimz+2*dimts,i),(/nnlv,nnkn/))
	 end do
      end if

      do i = 1,nens
         A(i)%u = reshape(Amat(ndim+1:ndim+dimuv,i),(/nnlv,nnel/))
         A(i)%v = reshape(Amat(ndim+dimuv+1:ndim+2*dimuv,i),(/nnlv,nnel/))
         A(i)%z = Amat(ndim+2*dimuv+1:ndim+2*dimuv+dimz,i)
      end do
      if (ibrcl > 0) then
         do i = 1,nens
            A(i)%t = reshape(Amat(ndim+2*dimuv+dimz+1:ndim+2*dimuv+dimz+dimts,i),(/nnlv,nnkn/))
            A(i)%s = reshape(Amat(ndim+2*dimuv+dimz+dimts+1:ndim+2*dimuv+dimz+2*dimts,i),(/nnlv,nnkn/))
         end do
      end if
   end subroutine matrix_to_tyqstate


end module mod_ens_state
