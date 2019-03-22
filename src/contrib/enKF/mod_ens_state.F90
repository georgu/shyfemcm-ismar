!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
!
module mod_ens_state

  use mod_dimensions
  use mod_mod_states
  use mod_para
  use mod_init_enkf

  implicit none

  type(states), dimension(:), allocatable :: Ashy      ! ensemble states
  type(states) :: Ashy_m                     ! mean state
  type(states) :: Ashy_stdb, Ashy_stda     ! standard deviation old and new

contains

!********************************************************

  subroutine read_ensemble

   implicit none

   character(len=3) :: nrel,nal
   character(len=16) rstname
   integer ne

   ! Allocates the state Ashy to store the ens states
   allocate(Ashy(nrens))

   ! init to zero
   do ne = 1,nrens
      Ashy(ne) = 0.
   end do

   call num2str(nanal,nal)

   if ((bnew_ens == 0) .or. (nanal.gt.1)) then

     write(*,*) 'loading an ensemble of initial states...'
     do ne = 1,nrens
        call num2str(ne-1,nrel)
        rstname='an'//nal//'_'//'en'//nrel//'b.rst'
        call read_state(Ashy(ne),rstname)
     end do

   else if ((bnew_ens == 1) .and. (nanal == 1)) then

     write(*,*) 'creating a new ensemble of initial states...'
     !read an input restart file
     call num2str(0,nrel)
     rstname = 'an'//nal//'_'//'en'//nrel//'b.rst'
     call read_state(Ashy(1),rstname)

     call make_init_ens(Ashy(1))
     
     !save the initial ens in new restart files
     call num2str(nanal,nal)
     do ne = 1,nrens
        call num2str(ne-1,nrel)
        rstname='an'//nal//'_'//'en'//nrel//'b.rst'
        call write_state(Ashy(ne),rstname)
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

   character(len=3) :: nrel,nal
   character(len=16) rstname
   integer ne

   write(*,*) 'writing the analysis restart files...'
   call num2str(nanal,nal)
   do ne = 1,nrens
      call num2str(ne-1,nrel)
      rstname='an'//nal//'_'//'en'//nrel//'b.rst'
      call rst_read(rstname,atime) !This is to load var not present 
                                   ! in the ens state. It should be removed.

      rstname='an'//nal//'_'//'en'//nrel//'a.rst'
      call write_state(Ashy(ne),rstname)
   end do
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
       Ashy(ne) = Ain
     else
       Apert = 0.
       do k = 1,nnkn
          Apert%z(k) = kvec(k,ne-1) * sigma_in
       end do
       Ashy(ne) = Ain + Apert
     end if

   end do

  end subroutine make_init_ens

!********************************************************

  subroutine make_mean_std(tflag,na)

  implicit none
  integer, intent(in) :: na
  character(len=1), intent(in) :: tflag
  character(len=3) :: nal
  character(len=80) :: filinm,filins

  call num2str(na,nal)
  filinm = 'an'//nal//'_mean_'//tflag//'.rst'
  filins = 'an'//nal//'_std_'//tflag//'.rst'

  Ashy_m = mean_state(nrens,Ashy)
  call write_state(Ashy_m,filinm)

  if (tflag == 'a') then
     Ashy_stda = std_state(nrens,Ashy)
     call write_state(Ashy_stda,filins)
  else
     Ashy_stdb = std_state(nrens,Ashy)
     call write_state(Ashy_stdb,filins)
  end if

  end subroutine make_mean_std

!********************************************************

  subroutine inflate_state
  ! Multiplicative state inflation  Whitaker J. S. et al. 2012
  ! 1- Relaxation-to-prior-spread (RTPS) method (Whitaker J. S. et al. 2012)
  !    Ashy' = Ashy' * (alpha * (Astdo - Astdn)/Astdn + 1)
  !    alpha ~ 0.1, see mod_para
  ! 2- Simple multiplication: Ashy' = Ashy' (1 + alpha)
  !
  implicit none

  type(states) :: Aaux, Apert
  integer ne

  if (type_infl == 1) then

     write(*,*) 'RTPS inflation, alpha = ',alpha_infl
     call rtps_inflation(alpha_infl,nrens,Ashy,Ashy_m,Ashy_stdb,Ashy_stda)

  else if (type_infl == 2) then

     write(*,*) 'Multiplication inflation, alpha = ',alpha_infl
     call mult_inflation(alpha_infl,nrens,Ashy,Ashy_m)


  else
  
     write(*,*) 'No inflation'
     return

  end if
  
  end subroutine inflate_state
  

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
  call rst_write(trim(filename),atime)
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
  call rst_read(filename,atime)
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
    implicit none
 
    type(states4),intent(inout) :: A4

    A4%z = znv
    ! no significant differences by using currents rather than transports
    A4%u = utlnv
    A4%v = vtlnv
    A4%t = tempv
    A4%s = saltv
   
   end subroutine push_state

!********************************************************

   subroutine pull_state(A4)
    use mod_hydro
    use mod_hydro_vel
    use mod_ts
    use mod_conz
    implicit none
 
    type(states4),intent(in) :: A4

    znv = A4%z
    ! no significant differences by using currents rather than transports
    utlnv = A4%u
    vtlnv = A4%v
    tempv = A4%t
    saltv = A4%s

    ! make zenv
    call layer_thick
   end subroutine pull_state

end module mod_ens_state
