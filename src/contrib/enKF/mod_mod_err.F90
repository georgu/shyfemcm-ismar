module mod_mod_err

 use mod_init_enkf
 use mod_mod_states
 
 implicit none

 ! parameters for the computation of the model error
 !
 integer nx_er,ny_er		!number of x and y grid points
 integer fmult_er		!mult factor to determine the supersampling 
 real theta_er			!rotation of the random fields (0 East, anticlockwise)
 real sigma_er			!standard deviation of the fields (level)
 double precision dt_er	!time between 2 analysis steps
 double precision tau_er	!time decorrelation (>=dt) of the old error: dq/dt = -(1/tau)*q

 type(qstate), allocatable :: Aaug(:) 	! double state with model error

 type(states), allocatable, private :: qA(:)		! model error

contains

!********************************************************

  subroutine init_moderr
  
  implicit none

   open(22, file='mod_err.info', status='old')
   read(22,*) nx_er,ny_er,fmult_er,theta_er,sigma_er,dt_er,tau_er
   close(22)

  end subroutine init_moderr

!********************************************************

  subroutine push_aug(Ain)
  
   use basin
   
   implicit none
   type(states), intent(inout)  :: Ain(nrens)
   integer nst 		!= nanal-1	number of time steps from the begin of the assimilation
   double precision alpha,rho

   type(states) :: A1
   real kvec(nnkn,nrens)

   ! Parameters for the sample fields
   integer nx,ny 	!grid dimension
   integer fmult	!the start ensemble is fmult*nrens
   real theta		!rotation of the fields (0 East, anticlockwise)

   real mfact
   integer ne,ie,n,k

   write(*,*) ''
   write(*,*) 'Creating an augmented state with model errors'
   write(*,*) ''

   !---------------------------------------
   ! defines parameters for the model error
   !---------------------------------------
   nst = nanal 
 
   if (tau_er < dt_er) error stop 'make_aug: parameter error'
 
   alpha = 1. - (dt_er/tau_er)
   rho=sqrt( (1.0-alpha)**2 /&
     (dt_er*(float(nst) - 2.0*alpha - float(nst)*alpha**2 + 2.0*alpha**(nst+1))) )
 
   !---------------------------------------
   !makes the new white noise field
   !---------------------------------------
   call make_2Dpert(kvec,nnkn,nrens,fmult_er,theta_er,nx_er,ny_er)
 
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
   call load_error(alpha,atime-dt_er)
 
   !---------------------------------------
   !compute the new state qA = Ain + sqrt(dt)*sigma*rho*q1
   !---------------------------------------
   ! change this if you want errors not only in zeta
   mfact = sqrt(dt_er) * sigma_er * rho
   do ne = 1,nrens
      A1 = states_real_mult(qA(ne),mfact)
      Ain(ne) = Ain(ne) + A1
   end do
    
   !---------------------------------------
   !make the augmented state Aaug = (Ain,qA)
   !---------------------------------------
   allocate(Aaug(nrens))
   do ne = 1,nrens
      call push_qstate(Ain(ne),qA(ne),Aaug(ne))
   end do
   deallocate(qA)
 
  end subroutine push_aug

!********************************************************

  subroutine pull_aug(date,time,Aout)
  
  use mod_ens_state
  
  implicit none
  integer,intent(in) :: date,time
  type(states), intent(out)  :: Aout(nrens)

  type(states),allocatable :: qA(:)
  character(len=3) :: nrel,nal
  character(len=19) rstname
  type(states4) :: A4
  integer ne

  write(*,*) 'Saving model errors'

  allocate(qA(nrens))
  do ne = 1,nrens
     call pull_qstate(Aout(ne),qA(ne),Aaug(ne))
  end do

  call num2str(nanal,nal)

  do ne = 1,nrens
     A4 = qA(ne)
     call pull_state(A4)

     call num2str(ne-1,nrel)
     rstname='an'//nal//'_'//'en'//nrel//'_err.rst'
     call rst_write(rstname,atime,date,time)
  end do

  end subroutine pull_aug

!********************************************************

  subroutine load_error(alpha,tt)
  
   use mod_ens_state
   
   implicit none
   double precision, intent(in) :: alpha
   double precision, intent(in) :: tt

   logical :: bfile
   integer ne
   character(len=3) :: nrel,nal
   character(len=19) rstname
   type(states4) :: A4
   type(states) :: A1,A2
   integer date,time

   real mfact,aalpha

   aalpha = alpha

   ! Old analysis step
   call num2str(nanal-1,nal)

   ! Check if error files exist
   do ne = 1,nrens
      call num2str(ne-1,nrel)
      rstname='an'//nal//'_'//'en'//nrel//'_err.rst'
      inquire(file=rstname, exist=bfile)
      if (.not. bfile) goto 777
   end do

   ! Add the old error to the new one
   write(*,*) 'Loading model error from files'
   do ne = 1,nrens
      call num2str(ne-1,nrel)
      rstname='an'//nal//'_'//'en'//nrel//'_err.rst'
      call rst_read(rstname,tt,date,time)
      call push_state(A4)

      ! q1 = alpha*q0 + sqrt(1-alpha**2)*w
      A1 = A4
      A1 = states_real_mult(A1,aalpha)

      mfact = sqrt(1 - alpha**2)
      A2 = states_real_mult(qA(ne),mfact)

      qA(ne) = A1 + A2

   end do

   return

 777 continue

   write(*,*) 'Model error files not found'

  end subroutine load_error


end module mod_mod_err
