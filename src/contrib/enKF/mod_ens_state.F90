module mod_ens_state

  use mod_init_enkf
  use mod_mod_states

  implicit none

  type(states), allocatable :: A(:)      ! ensemble states
  type(states) :: Am                     ! average state


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

   ! perturbation for ze
   call make_2Dpert(kvec,nnkn,nrens-1,fmult_in,theta_in,nx_in,ny_in)

   do ne = 1,nrens-1
     call assign_states(Apert(ne),0.)

     do ie = 1,nnel
        do n = 1,3
           k = nen3v(n,ie)
           Apert(ne)%ze(n,ie) = kvec(k,ne) * sigma_in
        end do
     end do

   end do

   ! The first state is unperturbed
   A(1) = Aaux
   do ne = 2,nrens
      A(ne) = add_states(Aaux,Apert(ne-1))
   end do

  end subroutine make_init_ens

!********************************************************

  subroutine average_mat(date,time,rstw)
  use mod_hydro
  use mod_ts
  implicit none
  integer,intent(in) :: date,time
  integer,intent(in) :: rstw
  integer ne
  character(len=16) :: rstname
  character(len=3) :: nal

  type(states4) :: A4

  ! make the average
  !
  Am = 0.
  do ne = 1,nrens
    Am = Am + A(ne)
  end do
  Am = states_real_mult(Am,1./float(nrens))

  ! write in a file
  !
  A4 = Am
  call pull_state(A4)
  call num2str(nanal,nal)
  if (rstw == -1) then
        write(*,*) 'writing the ens mean background restart'
        rstname = 'an'//nal//'_enavrb.rst'
        call rst_write(rstname,atime,date,time)
  elseif (rstw == -2) then
        write(*,*) 'writing the ens mean analysis restart'
        rstname = 'an'//nal//'_enavra.rst'
        call rst_write(rstname,atime,date,time)
  end if

  end subroutine average_mat

!********************************************************

   subroutine push_state(AA)
    use mod_hydro
    use mod_hydro_vel
    use mod_ts
    use mod_conz
    implicit none
    !real*4 hly(nnlv,nnel)
 
    type(states4) :: AA

    AA%z = znv
    !AA%ze = zenv
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
    !use basutil
    implicit none
 
    type(states4) :: AA

    znv = AA%z
    ! make zenv
    !AA%ze = zenv
    call layer_thick

    ! no significant differences by using currents rather than transports
    utlnv = AA%u
    vtlnv = AA%v
    tempv = AA%t
    saltv = AA%s

   end subroutine pull_state

!********************************************************

   subroutine layer_thick
   ! Set zenv bigger than the depth hm3v
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

end module mod_ens_state
