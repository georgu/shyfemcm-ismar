!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
!
module mod_ens_state

  use mod_dimensions
  use mod_mod_states
  use mod_para
  use mod_init_enkf

  implicit none

  type(states), dimension(:), allocatable :: A      ! ensemble states
  type(states) :: Am                     ! mean state
  type(states) :: Astd_b, Astd_a     ! standard deviation old and new

contains

!********************************************************

  subroutine read_ensemble

   implicit none

   character(len=3) :: nrel,nal
   character(len=16) rstname
   integer ne

   ! Allocates the state A to store the ens states
   allocate(A(nrens))

   ! init to zero
   do ne = 1,nrens
      A(ne) = 0.
   end do

   call num2str(nanal,nal)

   if ((bnew_ens == 0) .or. (nanal.gt.1)) then

     write(*,*) 'loading an ensemble of initial states...'
     do ne = 1,nrens
        call num2str(ne-1,nrel)
        rstname='an'//nal//'_'//'en'//nrel//'b.rst'
        call read_state(A(ne),rstname)
     end do

   else if ((bnew_ens == 1) .and. (nanal == 1)) then

     write(*,*) 'creating a new ensemble of initial states...'
     !read an input restart file
     call num2str(0,nrel)
     rstname = 'an'//nal//'_'//'en'//nrel//'b.rst'
     call read_state(A(1),rstname)

     call make_init_ens(A(1))
     
     !save the initial ens in new restart files
     call num2str(nanal,nal)
     do ne = 1,nrens
        call num2str(ne-1,nrel)
        rstname='an'//nal//'_'//'en'//nrel//'b.rst'
        call write_state(A(ne),rstname)
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
      call write_state(A(ne),rstname)
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

   type(states), allocatable, save :: Apert
   real kvec(nnkn,nrens-1),evec(nnel,nrens-1)

   integer ne,n,ie,k

   open(21, file='init_ens.info', status='old')
   read(21,*) nx_in,ny_in,fmult_in,theta_in,sigma_in
   close(21)

   allocate(Apert)

   ! perturbation for z
   call make_2Dpert(kvec,nnkn,nrens-1,fmult_in,theta_in,nx_in,ny_in)

   do ne = 1,nrens

     if (ne == 1) then
       A(ne) = Ain
     else
       Apert = 0.
       do k = 1,nnkn
          Apert%z(k) = kvec(k,ne-1) * sigma_in
       end do
       A(ne) = Ain + Apert
     end if

   end do

   deallocate(Apert)

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

  call mean_state
  call std_state(tflag)

  call write_state(Am,filinm)
  if (tflag == 'a') then
     call write_state(Astd_a,filins)
  else
     call write_state(Astd_b,filins)
  end if

  end subroutine make_mean_std


!********************************************************

  subroutine mean_state
  ! make the mean state
  implicit none
  integer ne
  real inrens
  real, parameter :: eps = 1.e-15

  Am = eps
  do ne = 1,nrens
    Am = Am + A(ne)
  end do
  inrens = 1./float(nrens)
  Am = Am * inrens

  end subroutine mean_state

!********************************************************

  subroutine std_state(label)
  ! make the standard deviation of the states
  implicit none
  character(len=1), intent(in) :: label

  type(states), allocatable :: Astd,Apert
  integer ne
  real inrens
  real, parameter :: eps = 1.e-15

  allocate(Astd,Apert)
  Astd = eps
  do ne = 1,nrens
     Apert =  A(ne) - Am
     Apert =  Apert * Apert
     Astd =  Astd + Apert
  end do
  deallocate(Apert)

  inrens = 1./float(nrens-1)
  Astd = Astd * inrens
  Astd = root_state(Astd)

  if (label == 'a') then
     Astd_a = Astd
  else
     Astd_b = Astd
  end if
  deallocate(Astd)

  end subroutine std_state

!********************************************************

  subroutine inflate_state
  ! Multiplicative state inflation  Whitaker J. S. et al. 2012
  ! 1- Relaxation-to-prior-spread (RTPS) method (Whitaker J. S. et al. 2012)
  !    A' = A' * (alpha * (Astdo - Astdn)/Astdn + 1)
  !    alpha ~ 0.1, see mod_para
  ! 2- Simple multiplication: A' = A' (1 + alpha)
  !
  implicit none

  type(states), allocatable :: Aaux, Apert
  integer ne

  if (type_infl == 1) then

     write(*,*) 'RTPS inflation, alpha = ',alpha_infl

     allocate(Aaux,Apert)

     Aaux = Astd_b - Astd_a
     Aaux = Aaux / Astd_a
     Aaux = alpha_infl * Aaux 
     Aaux = Aaux + 1.

     do ne = 1,nrens
        Apert = A(ne) - Am
        A(ne) = Am + (Apert * Aaux)
     enddo
     deallocate(Aaux,Apert)

  else if (type_infl == 2) then

     write(*,*) 'Multiplication inflation, alpha = ',alpha_infl

     allocate(Apert)
     do ne = 1,nrens
        Apert = A(ne) - Am
        A(ne) = Am + (Apert * alpha_infl)
     enddo
     deallocate(Apert)

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
  A4 = Astate
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
  A4 = Astate
  call rst_read(filename,atime)
  call push_state(A4)
  Astate = A4
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
