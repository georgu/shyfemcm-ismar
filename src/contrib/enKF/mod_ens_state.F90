module mod_ens_state

  use mod_init_enkf
  use mod_mod_states
  use mod_para

  implicit none

  type(states), allocatable :: A(:)      ! ensemble states
  type(states) :: Am                     ! mean state
  type(states) :: Astd_old, Astd_new     ! standard deviation old and new


contains

!********************************************************

  subroutine read_ensemble(date,time)

   implicit none

   integer,intent(out) :: date,time

   type(states4) :: A4
   type(states) :: Ap(nrens-1)
   character(len=3) :: nrel,nal
   character(len=16) rstname
   integer ne

   ! Allocates the state A to store the ens states
   allocate(A(nrens))

   ! init to zero
   do ne = 1,nrens
      A(ne) = 0.
   end do
   A4 = A(1)

   call num2str(nanal,nal)

   if ((bnew_ens == 0) .or. (nanal.gt.1)) then

     write(*,*) 'loading an ensemble of initial states...'
     do ne = 1,nrens
        call num2str(ne-1,nrel)
        rstname='an'//nal//'_'//'en'//nrel//'b.rst'
        call rst_read(rstname,atime,date,time)
        call push_state(A4)
        A(ne)=A4
     end do

   else if ((bnew_ens == 1) .and. (nanal == 1)) then

     write(*,*) 'creating a new ensemble of initial states...'
     !read an input restart file
     call num2str(0,nrel)
     rstname = 'an'//nal//'_'//'en'//nrel//'b.rst'
     call rst_read(rstname,atime,date,time)

     !push the vars into the state and makes the ens
     call push_state(A4)
     call make_init_ens(A4)
     
     !save the initial ens in new restart files
     call num2str(nanal,nal)
     do ne = 1,nrens
        call num2str(ne-1,nrel)
        A4 = A(ne)
        call pull_state(A4)
        rstname='an'//nal//'_'//'en'//nrel//'b.rst'
        call rst_write(rstname,atime,date,time)
     end do

   else

     write(*,*) 'not a valid option for bnew_ens'
     error stop

   end if

   return
  end subroutine read_ensemble

!********************************************************

  subroutine write_ensemble(date,time)
   implicit none
   integer,intent(in) :: date,time

   type(states4) :: A4
   character(len=3) :: nrel,nal
   character(len=16) rstname
   integer ne
   integer dt,tm

   write(*,*) 'writing the analysis restart files...'
   call num2str(nanal,nal)
   do ne = 1,nrens
      call num2str(ne-1,nrel)
      rstname='an'//nal//'_'//'en'//nrel//'b.rst'
      call rst_read(rstname,atime,dt,tm) !This is to load var not present 
                                        ! in the ens state. It should be removed.
      A4 = A(ne)
      call pull_state(A4)
      rstname='an'//nal//'_'//'en'//nrel//'a.rst'
      call rst_write(rstname,atime,date,time)
   end do
  end subroutine write_ensemble

!********************************************************

  subroutine make_init_ens(A4)
   use basin

   implicit none

   type(states4),intent(in) :: A4

   ! parameters for the initial ensemble of states
   !
   integer nx_in,ny_in           !number of x and y grid points
   integer fmult_in              !mult factor to determine the supersampling
   real theta_in                 !rotation of the random fields (0 East, anticlockwise)
   real sigma_in                 !standard deviation of the fields (level)

   type(states), allocatable, save :: Apert(:)
   type(states) :: Aaux
   real kvec(nnkn,nrens-1),evec(nnel,nrens-1)

   integer ne,n,ie,k

   open(21, file='init_ens.info', status='old')
   read(21,*) nx_in,ny_in,fmult_in,theta_in,sigma_in
   close(21)

   Aaux = A4

   allocate(Apert(nrens-1))

   ! perturbation for z
   call make_2Dpert(kvec,nnkn,nrens-1,fmult_in,theta_in,nx_in,ny_in)

   do ne = 1,nrens-1
     call assign_states(Apert(ne),0.)

     do k = 1,nnkn
        Apert(ne)%z(k) = kvec(k,ne) * sigma_in
     end do

   end do

   ! The first state is unperturbed
   A(1) = Aaux
   do ne = 2,nrens
      A(ne) = add_states(Aaux,Apert(ne-1))
   end do

  end subroutine make_init_ens

!********************************************************

  subroutine mean_state(Amean)
  ! make the mean state
  implicit none
  type(states),intent(out) :: Amean
  integer ne
  real inrens
  real, parameter :: eps = 1.e-15

  Amean = eps
  do ne = 1,nrens
    Amean = Amean + A(ne)
  end do
  inrens = 1./float(nrens)
  Amean = Amean * inrens

  end subroutine mean_state

!********************************************************

  subroutine std_state(Amean,Astd)
  ! make the standard deviation of the states
  implicit none
  type(states),intent(in) :: Amean
  type(states),intent(out) :: Astd
  integer ne
  real inrens
  real, parameter :: eps = 1.e-15

  Astd = eps
  do ne = 1,nrens
     Astd = Astd + ((A(ne) - Amean) * (A(ne) - Amean))
  end do

  inrens = 1./float(nrens-1)
  Astd = Astd * inrens
  Astd = root_state(Astd)

  end subroutine std_state

!********************************************************

  subroutine inflate_state(Astdo,Astdn,Amean)
  ! Multiplicative state inflation according to Whitaker J. S. et al. 2012
  ! Relaxation-to-prior-spread (RTPS) method
  ! A = A * (alpha * (Astdo - Astdn)/Astdn + 1)
  ! alpha ~ 0.5 
  ! see mod_para
  implicit none
  type(states),intent(in) :: Astdo,Astdn,Amean
  type(states) :: Aaux, Apert
  integer ne

  write(*,*) 'RTPS inflation, alpha = ',alpha_infl

  Aaux = Astdo - Astdn
  Aaux = Aaux / Astdn
  Aaux = alpha_infl * Aaux 
  Aaux = Aaux + 1.

  do ne = 1,nrens
     Apert = (A(ne) - Amean) * Aaux
     A(ne) = Amean + Apert 
  enddo
  
  end subroutine inflate_state
  

!********************************************************

  subroutine write_state(date,time,Astate,filename)
  use mod_hydro
  use mod_ts
  implicit none
  integer,intent(in) :: date,time
  type(states),intent(in) :: Astate
  character(len=*),intent(in) :: filename

  type(states4) :: A4

  A4 = Astate
  call pull_state(A4)
  write(*,*) 'writing file: ',trim(filename)
  call rst_write(trim(filename),atime,date,time)
  end subroutine write_state


!********************************************************

   subroutine push_state(AA)
    use mod_hydro
    use mod_hydro_vel
    use mod_ts
    use mod_conz
    implicit none
 
    type(states4) :: AA

    AA%z = znv
    ! no significant differences by using currents rather than transports
    AA%u = utlnv
    AA%v = vtlnv
    AA%t = tempv
    AA%s = saltv
   
   end subroutine push_state

!********************************************************

   subroutine pull_state(AA)
    use mod_hydro
    use mod_hydro_vel
    use mod_ts
    use mod_conz
    implicit none
 
    type(states4) :: AA

    znv = AA%z
    ! no significant differences by using currents rather than transports
    utlnv = AA%u
    vtlnv = AA%v
    tempv = AA%t
    saltv = AA%s

    ! make zenv
    call layer_thick
   end subroutine pull_state

end module mod_ens_state
