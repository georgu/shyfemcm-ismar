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

   implicit none

   type(states),intent(in) :: Ain

   ! parameters for the initial ensemble of states
   !
   integer nx_in,ny_in           !number of x and y grid points
   integer fmult_in              !mult factor to determine the supersampling
   real theta_in                 !rotation of the random fields (0 East, anticlockwise)
   real sigma_in                 !standard deviation of the fields (level)

   type(states), save :: Apert
   real kvec(nnkn,nrens-1),evec(nnel,nrens-1)

   integer ne,n,ie,k

   open(21, file='init_ens.info', status='old')
   read(21,*) nx_in,ny_in,fmult_in,theta_in,sigma_in
   close(21)

   ! perturbation for z
   call make_2Dpert(kvec,nnkn,nrens-1,fmult_in,theta_in,nx_in,ny_in)

   do ne = 1,nrens

     if (ne == 1) then
       Abk(ne) = Ain
     else
       Apert = 0.
       do k = 1,nnkn
          Apert%z(k) = kvec(k,ne-1) * sigma_in
       end do
       Abk(ne) = Ain + Apert
     end if

   end do

  end subroutine make_init_ens

!********************************************************

  subroutine check_and_correct(Abkg,Aanl)
   implicit none
   type(states),intent(in) :: Abkg(nrens)
   type(states),intent(inout) :: Aanl(nrens)

   real vala,valb,vala2,valb2

   integer ne,k,ie,nl

   integer zout,uvout,sout,tout,otot
   integer znan,uvnan,snan,tnan,ntot

   write(*,*) 'Checking the analysis values of the variables...'

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

    do k = 1,nnkn

      !---
      ! level
      valb = Abkg(ne)%z(k)
      vala = Aanl(ne)%z(k)

      ! isnan and bk is a number
      if ( isnan(vala) .and. (.not. isnan(valb) ) ) then
        !if (verbose) write(*,*) 'Warning: nan in water level. Node: ',k
        vala = valb
        znan = znan + 1
      end if

      ! out of range
      if (( vala > SSH_MAX ) .or. ( vala < SSH_MIN )) then
	vala = valb
	zout = zout + 1
      end if

      Aanl(ne)%z(k) = vala
      !---

      do nl = 1,nnlv

        !---
        ! salinity
	valb = Abkg(ne)%s(nl,k)
	vala = Aanl(ne)%s(nl,k)

        ! isnan and bk is a number
        if ( isnan(vala) .and. (.not. isnan(valb) ) ) then
           !if (verbose) write(*,*) 'Warning: nan in salility. Node, level: ',k,nl
	   vala = valb
	   snan = snan + 1
	end if
	! out of range
        if (( vala > SAL_MAX ) .or. ( vala < SAL_MIN )) then
          vala = valb
	  sout = sout + 1
	end if

	Aanl(ne)%s(nl,k) = vala
        !---

        !---
	! temperature
	valb = Abkg(ne)%t(nl,k)
	vala = Aanl(ne)%t(nl,k)
	
        ! isnan and bk is a number
        if ( isnan(vala) .and. (.not. isnan(valb)) ) then
           !if (verbose) write(*,*) 'Warning: nan in temperature. Node, level: ',k,nl
	   vala = valb
	   tnan = tnan + 1
	end if
	! out of range
        if (( vala > TEM_MAX ) .or. ( vala < TEM_MIN )) then
          vala = valb
	  tout = tout + 1
	end if

	Aanl(ne)%t(nl,k) = vala
        !---

      end do

    end do


    do ie = 1,nnel
      do nl = 1,nnlv
	
        !---
        ! current
	valb = Abkg(ne)%u(nl,ie)
	vala = Aanl(ne)%u(nl,ie)

	valb2 = Abkg(ne)%v(nl,ie)
	vala2 = Aanl(ne)%v(nl,ie)

	! isnan and bk is a number
        if ( isnan(vala) .and. (.not. isnan(valb)) ) then
	   !if (verbose) write(*,*) 'Warning: nan in velocity. El., level: ',ie,nl
	   vala = valb
	   uvnan = uvnan + 1
	end if
        if ( isnan(vala2) .and. (.not. isnan(valb2)) ) then
	   !if (verbose) write(*,*) 'Warning: nan in velocity. El., level: ',ie,nl
	   vala2 = valb2
	   uvnan = uvnan + 1
	end if

	if (( vala > VEL_MAX ) .or. ( vala < VEL_MIN )) then
	   vala = valb
	   uvout = uvout + 1
	end if
	if (( vala2 > VEL_MAX ) .or. ( vala2 < VEL_MIN )) then
	   vala2 = valb2
	   uvout = uvout + 1
	end if

	Aanl(ne)%u(nl,ie) = vala
	Aanl(ne)%v(nl,ie) = vala2
        !---

      end do
    end do

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

  end subroutine check_and_correct


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
