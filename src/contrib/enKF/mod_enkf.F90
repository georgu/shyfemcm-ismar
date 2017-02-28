  module mod_enkf

  use mod_states
  use mod_observations

  implicit none

  character (len=80), save :: basfile, obsfile, rstfile

  ! parameters for the initial ensemble of states
  integer, save :: is_new_ens
  integer nx_in,ny_in		!number of x and y grid points
  integer fmult_in		!mult factor to determine the supersampling
  real theta_in			!rotation of the random fields (0 East, anticlockwise)
  real sigma_in			!standard deviation of the fields (level)

  ! parameters for the computation of the model error
  integer, save :: is_mod_err
  integer nx_er,ny_er		!number of x and y grid points
  integer fmult_er		!mult factor to determine the supersampling 
  real theta_er			!rotation of the random fields (0 East, anticlockwise)
  real sigma_er			!standard deviation of the fields (level)
  double precision dt_er	!time between 2 analysis steps
  double precision tau_er	!time decorrelation (>=dt) of the old error: dq/dt = -(1/tau)*q

  integer, save :: nrens, na

  ! Observations
  double precision, save :: tobs
  integer, save :: nobs_lev,nobs_tot
  type(levels), save :: olev

  type(states), allocatable, save  :: A(:) 	! ensemble states
  type(states), save  :: Am			! average state
  type(states), allocatable :: qA(:)		! model error
  type(dstates), allocatable, save  :: Aaug(:) 	! double state with model error

  real, save, allocatable :: D(:,:)		! matrix holding perturbed measurments
  real, save, allocatable :: E(:,:)		! matrix holding perturbations (mode=?3)
  real, save, allocatable :: R(:,:)		! Obs error cov matrix

  real, save, allocatable :: S(:,:)		! matrix holding HA`
  real, save, allocatable :: innov(:)		! innovation vector holding d-H*mean(A)


  contains


!********************************************************
!********************************************************
!********************************************************

  subroutine read_info

  implicit none

  integer n

  open(20, file='analysis.info', status='old')

  read(20,*) nrens	! number of ens members
  read(20,*) na		! analysis step
  read(20,*) basfile	! name of bas file (no extension)
  read(20,*) tobs	! current time of the observations
  read(20,*) obsfile	! name of obs file list
  read(20,*) is_new_ens	! 1 to create a new initial ens of states
  read(20,*) is_mod_err	! 1 to use an augmented state with mod err

  close(20)

  if( mod(nrens,2).eq.0 ) stop 'read_info: n of ens members must be odd, with the control as first.'

  if( is_new_ens.eq.1 ) then
    open(21, file='init_ens.info', status='old')
    read(21,*) nx_in,ny_in,fmult_in,theta_in,sigma_in	
    close(21)
  end if
  if( is_mod_err.eq.1 ) then
    open(22, file='mod_err.info', status='old')
    read(22,*) nx_er,ny_er,fmult_er,theta_er,sigma_er,dt_er,tau_er
    close(22)
  end if
  
  end subroutine read_info

!********************************************************

  subroutine read_basin

  use basin
  implicit none

  open(21, file=basfile, status='old', form='unformatted')
  call basin_read_by_unit(21)
  close(21)

  if( ( nkn.ne.nnkn ).or.( nel.ne.nnel) ) stop "read_basin: dim error"

  end subroutine read_basin


!********************************************************

  subroutine read_obs
  implicit none

  integer ios
  character(len=80) :: line
  character(len=80), allocatable :: ofile(:)
  integer, allocatable :: nrec(:)

  double precision :: eps = 300.	! 300 seconds
  double precision tt
  character(len=6) :: ty
  real x, y, z, v, stdv
  integer n,nfile,klev

     write(*,*) 'Observation file list: ',trim(obsfile)

!    Reads the obs list
     n = 1
     open(25,file = obsfile, status = 'old', form = 'formatted', iostat = ios)
     if( ios.ne.0 ) stop 'read_obs: error opening file'
 88   read(25,*,end=98) line
      n = n + 1
     goto 88
 98  continue
     nfile = n - 1
     rewind(unit=25)

     allocate(ofile(nfile),nrec(nfile))

     do n = 1,nfile
        read(25,*,err=95) ofile(n)
     end do
     close(25)

!-------------------------------
! Reads every file one time to find the number and the type of each obs
!-------------------------------
  klev = 0
  do n = 1,nfile

     open(26,file=ofile(n), status = 'old', form = 'formatted', iostat = ios)
     if( ios.ne.0 ) stop 'read_obs: error opening file'

 89  read(26,*,end=99) tt, ty

     ! Takes only records with times near tobs
     if( abs(tt - tobs) .lt. eps ) then
        if( trim(ty).eq.'level' ) then
          klev = klev + 1
        else
          stop 'Observation type not still implemented.'
        end if
     end if

     goto 89

 99  close(26)
  end do
  nobs_lev = klev

  ! allocate
  if( nobs_lev.gt.0 ) then
      allocate(olev%t(nobs_lev))
      allocate(olev%x(nobs_lev))
      allocate(olev%y(nobs_lev))
      allocate(olev%val(nobs_lev))
      allocate(olev%std(nobs_lev))
  end if

!-------------------------------
! reads the second time and store
!-------------------------------
  klev = 0
  do n = 1,nfile

     open(26,file=ofile(n), status = 'old', form = 'formatted', iostat = ios)
     if( ios.ne.0 ) stop 'read_obs: error opening file'

 90  read(26,*,end=100) tt, ty, x, y, z, v, stdv

     ! Takes only records with times near tobs
     if( abs(tt - tobs) .lt. eps ) then
        if( trim(ty).eq.'level' ) then
          klev = klev + 1

          olev%t(klev) = tt
          olev%x(klev) = x
          olev%y(klev) = y
          olev%val(klev) = v
          olev%std(klev) = stdv
        else
          stop 'Observation type not still implemented.'
        end if
     end if

     goto 90

 100  close(26)
  end do

  nobs_tot = nobs_lev
  return

 95 write(*,*) 'read_obs error reading file: ',trim(ofile(n))
    stop

  end subroutine read_obs

!********************************************************

  subroutine average_mat(rstw)

  use mod_hydro
  use mod_ts
  implicit none
  integer ne
  integer rstw
  character(len=16) :: rstname
  character(len=3) :: nal

  type(states4) :: A4

  ! makes the average
  Am = 0.
  do ne = 1,nrens
    Am = Am + A(ne)
  end do
  Am = states_real_mult(Am,1./float(nrens))

  ! writes in a file
  A4 = Am
  call pull_state(A4)
  call num2str(na,nal)
  if( rstw.eq.-1 ) then
        write(*,*) 'Writing average background state...'
        rstname = 'an'//nal//'_enavrb.rst'
        call rst_write(rstname,tobs)
  elseif( rstw.eq.-2 ) then
        write(*,*) 'Writing average analysis state...'
        rstname = 'an'//nal//'_enavra.rst'
        call rst_write(rstname,tobs)
  end if

  end subroutine average_mat


!********************************************************

  subroutine make_matrices

! R (only used if mode=?1 or ?2) (no for low-rank sq root)

  use m_random
  implicit none

  integer iel
  real*4 x4,y4

  real rand_v(nrens)
  real inn,ave,var
  integer ne, i

  integer n

  ! Reads observations and defines nobs_tot
  call read_obs

  allocate(D(nobs_tot,nrens),E(nobs_tot,nrens),R(nobs_tot,nobs_tot))
  allocate (S(nobs_tot,nrens),innov(nobs_tot))

  R(:,:) = 0.

  do n = 1,nobs_tot

     if( nobs_tot.ne.nobs_lev ) stop 'Other obs not implemented'

     ! create ensemble of observation perturbations
     call random(rand_v,nrens)
     ave = sum(rand_v)/float(nrens)
     rand_v = rand_v - ave
     var = dot_product(rand_v,rand_v)/float(nrens-1)
     rand_v = sqrt(1.0/var)*rand_v
     E(n,:) = olev%std(n) * rand_v(:)

     ! Finds the finite element nearest to the obs
     x4 = olev%x(n)
     y4 = olev%y(n)
     call find_element(x4,y4,iel)

     ! compute S matrix (HA') and the ensemble of innovations D
     do ne = 1,nrens
        S(n,ne) = sum( A(ne)%ze(:,iel) - Am%ze(:,iel) )/3.
        D(n,ne) = olev%val(n) + E(n,ne) - sum(A(ne)%ze(:,iel))/3.
     end do

     ! find innovation
     ave = sum(Am%ze(:,iel))/3.
     inn = olev%val(n) - ave
     call check_innov_val(inn,'level')
     innov(n) = inn
     write(*,'(a,i5,3f8.4)') ' nobs, vobs, vmod, innov: ',n,olev%val(n),ave,inn
  
     ! optional R
     !R(n,n) = olev%std(n)**2
  end do

  end subroutine make_matrices

!********************************************************

  subroutine read_ensemble

   implicit none

   integer ne

   type(states4) :: A4
   type(states) :: Ap(nrens-1)

   character(len=3) :: nrel,nal
   character(len=16) rstname

   ! Allocates the state A to store the ens states
   allocate(A(nrens))

   call num2str(na,nal)

   if( (is_new_ens.eq.0).or.(na.gt.1) ) then

     write(*,*) 'Loading an ensemble of initial states'
     do ne = 1,nrens
        call num2str(ne-1,nrel)
        rstname='an'//nal//'_'//'en'//nrel//'b.rst'
        call rst_read(nnkn,nnel,nnlv,rstname,tobs)
        call push_state(A4)
        A(ne)=A4
     end do

   else if( (is_new_ens.eq.1).and.(na.eq.1) ) then

     write(*,*) '*****************************************'
     write(*,*) 'Creating a new ensemble of initial states'
     write(*,*) '*****************************************'
     !read an input restart file
     call num2str(0,nrel)
     rstname = 'an'//nal//'_'//'en'//nrel//'b.rst'
     call rst_read(nnkn,nnel,nnlv,rstname,tobs)

     !push the vars into the state and makes the ens
     call push_state(A4)
     call make_init_ens(A4)
     
     !save the initial ens in new restart files
     call num2str(na,nal)
     do ne = 1,nrens
        call num2str(ne-1,nrel)
        A4 = A(ne)
        call pull_state(A4)
        rstname='an'//nal//'_'//'en'//nrel//'b.rst'
        call rst_write(rstname,tobs)
     end do

   else

     write(*,*) 'Not a valid option for is_new_ens'
     stop

   end if

   return
  end subroutine read_ensemble


!********************************************************

  subroutine write_ensemble
   implicit none

   integer ne

   type(states4) :: A4

   character(len=3) :: nrel,nal
   character(len=16) rstname

   call num2str(na,nal)
   do ne = 1,nrens
      call num2str(ne-1,nrel)
      rstname='an'//nal//'_'//'en'//nrel//'b.rst'
      call rst_read(nnkn,nnel,nnlv,rstname,tobs) !This is to load var not present 
                                                 ! in the ens state. It must be removed.
      A4 = A(ne)
      call pull_state(A4)
      rstname='an'//nal//'_'//'en'//nrel//'a.rst'
      call rst_write(rstname,tobs)
   end do
  end subroutine write_ensemble

!********************************************************

   subroutine push_state(AA)
    use mod_hydro
    use mod_hydro_vel
    use mod_ts
    use mod_conz
    implicit none
 
    type(states4) :: AA

    AA%u = utlnv
    AA%v = vtlnv
    AA%ze = zenv
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

    utlnv = AA%u
    vtlnv = AA%v
    zenv = AA%ze
    tempv = AA%t
    saltv = AA%s
   
   end subroutine pull_state


!********************************************************

  subroutine make_init_ens(A4)
   use basin
   implicit none
   type(states4),intent(in) :: A4

   type(states), allocatable, save :: Apert(:)
   type(states) :: Aaux
   real kvec(nnkn,nrens-1),evec(nnel,nrens-1)

   integer ne,n,ie,k

   Aaux = A4

   allocate(Apert(nrens-1))

   ! perturbation for ze
   call make_pert(kvec,nnkn,nrens-1,fmult_in,theta_in,nx_in,ny_in)

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

  subroutine push_aug
   use basin
   implicit none
   integer nst 		!= na-1	number of time steps from the begin of the assimilation
   double precision alpha,rho

   type(states) :: A1
   real kvec(nnkn,nrens)

   ! Parameters for the sample fields
   integer nx,ny 	!grid dimension
   integer fmult	!the start ensemble is fmult*nrens
   real theta		!rotation of the fields (0 East, anticlockwise)

   real mfact
   integer ne,ie,n,k

   write(*,*) '*********************************************'
   write(*,*) 'Creating an augmented state with model errors'
   write(*,*) '*********************************************'

   !---------------------------------------
   ! defines parameters for the model error
   !---------------------------------------
   nst = na 
 
   if( tau_er.lt.dt_er ) stop 'make_aug: parameter error'
 
   alpha = 1. - (dt_er/tau_er)
   rho=sqrt( (1.0-alpha)**2 / (dt_er*(float(nst) - 2.0*alpha - float(nst)*alpha**2 + 2.0*alpha**(nst+1))) )
 
   !---------------------------------------
   !makes the new white noise field
   !---------------------------------------
   call make_pert(kvec,nnkn,nrens,fmult_er,theta_er,nx_er,ny_er)
 
   allocate(qA(nrens))
   do ne = 1,nrens
      call assign_states(qA(ne),0.)
 
      do ie = 1,nnel
        do n = 1,3
          k = nen3v(n,ie)
          qA(ne)%ze(n,ie) = kvec(k,ne)
        end do
      end do

   end do
 
   !---------------------------------------
   !if exist old error q0 load it and add to qA (q1 = alpha*q0 + sqrt(1-alpha**2)*w)
   !---------------------------------------
   call load_error(alpha,tobs-dt_er)
 
   !---------------------------------------
   !compute the new state qA = A + sqrt(dt)*sigma*rho*q1
   !---------------------------------------
   ! change this if you want errors not only in zeta
   mfact = sqrt(dt_er) * sigma_er * rho
   do ne = 1,nrens
      A1 = states_real_mult(qA(ne),mfact)
      A(ne) = A(ne) + A1
   end do
    
   !---------------------------------------
   !make the augmented state Aaug = (A,qA)
   !---------------------------------------
   allocate(Aaug(nrens))
   do ne = 1,nrens
      call push_dstate(A(ne),qA(ne),Aaug(ne))
   end do
   deallocate(A,qA)
 
  end subroutine push_aug

!********************************************************

  subroutine load_error(alpha,tt)

   implicit none
   real, intent(in) :: alpha
   double precision, intent(in) :: tt

   logical :: file_exist
   integer ne
   character(len=3) :: nrel,nal
   character(len=19) rstname
   type(states4) :: A4
   type(states) :: A1,A2

   real mfact

   ! Old analysis step
   call num2str(na-1,nal)

   ! Check if error files exist
   do ne = 1,nrens
      call num2str(ne-1,nrel)
      rstname='an'//nal//'_'//'en'//nrel//'_err.rst'
      inquire(file=rstname, exist=file_exist)
      if(.not. file_exist) goto 777
   end do

   ! Add the old error to the new one
   write(*,*) '********'
   write(*,*) 'Loading model error from files'
   do ne = 1,nrens
      call num2str(ne-1,nrel)
      rstname='an'//nal//'_'//'en'//nrel//'_err.rst'
      call rst_read(nnkn,nnel,nnlv,rstname,tt)
      call push_state(A4)

      ! q1 = alpha*q0 + sqrt(1-alpha**2)*w
      A1 = A4
      A1 = states_real_mult(A1,alpha)

      mfact = sqrt(1 - alpha**2)
      A2 = states_real_mult(qA(ne),mfact)

      qA(ne) = A1 + A2

   end do

   return

 777 continue

   write(*,*) '********'
   write(*,*) 'Model error files not found'

  end subroutine load_error

!********************************************************

  subroutine pull_aug
   implicit none

   integer ne
   character(len=3) :: nrel,nal
   character(len=19) rstname
   type(states4) :: A4

   write(*,*) '********'
   write(*,*) 'Saving model errors'

   allocate(A(nrens),qA(nrens))
   do ne = 1,nrens
      call pull_dstate(A(ne),qA(ne),Aaug(ne))
   end do
   deallocate(Aaug)

   call num2str(na,nal)

   do ne = 1,nrens
      A4 = qA(ne)
      call pull_state(A4)

      call num2str(ne-1,nrel)
      rstname='an'//nal//'_'//'en'//nrel//'_err.rst'
      call rst_write(rstname,tobs)
   end do

  end subroutine pull_aug

  end module mod_enkf
