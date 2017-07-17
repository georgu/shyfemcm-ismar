module mod_mod_err

 use mod_init_enkf
 use mod_mod_states
 use mod_ens_state
 
 implicit none

 ! parameters for the computation of the model error
 !
 integer nx_er,ny_er		!number of x and y grid points
 integer fmult_er		!mult factor to determine the supersampling 
 real theta_er			!rotation of the random fields (0 East, anticlockwise)
 real rsigma			!relative error
 double precision dt_er	!time between 2 analysis steps
 double precision tau_er	!time decorrelation (>=dt) of the old error: dq/dt = -(1/tau)*q

 type(qstate), allocatable :: Aaug(:) 	! double state with model error

 type(states), allocatable, private :: qA(:)		! model error

contains

!********************************************************

  subroutine init_moderr
  
  implicit none

   open(22, file='mod_err.info', status='old')
   read(22,*) nx_er,ny_er,fmult_er,theta_er,rsigma,dt_er,tau_er
   close(22)

  end subroutine init_moderr

!********************************************************

  subroutine push_aug(Ain)
  
   use basin
   
   implicit none
   type(states), intent(inout)  :: Ain(nrens)
   integer nst 		!= nanal-1	number of time steps from the begin of the assimilation
   double precision alpha,rho

   type(states) :: A1,A2
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
 
      qA(ne)%z = kvec(:,ne)

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
   call load_error(alpha)
 
   !---------------------------------------
   !compute the new state qA = Ain + sqrt(dt)*sigma*rho*q1
   !---------------------------------------
   mfact = sqrt(dt_er) * rsigma * rho
   do ne = 1,nrens
      ! find the spatial factor from the relative error (rsigma)
      A1 = states_real_mult(A(ne),mfact)
      ! find the error
      A2 = states_states_mult(qA(ne),A1)
      ! add the error
      Ain(ne) = Ain(ne) + A2
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

  subroutine pull_aug(Aout)
  
  implicit none

  type(states), intent(out)  :: Aout(nrens)

  type(states),allocatable :: qA(:)
  character(len=3) :: nrel,nal
  character(len=16) fname
  type(states4) :: A4
  integer ne

  write(*,*) 'Saving model errors'

  allocate(qA(nrens))
  do ne = 1,nrens
     call pull_qstate(Aout(ne),qA(ne),Aaug(ne))
  end do

  call num2str(nanal,nal)
  fname='an'//nal//'_moderr.bin'
  open(33,file=fname,form='unformatted')
  write(33) qA
  close(33)

  deallocate(qA)
  deallocate(Aaug)

  end subroutine pull_aug

!********************************************************

  subroutine load_error(alpha)
  
   implicit none
   double precision, intent(in) :: alpha

   logical :: bfile
   integer ne
   character(len=3) :: nrel,nal
   character(len=16) fname
   type(states),allocatable :: A1(:),A2(:)

   real mfact,aalpha

   aalpha = alpha

   allocate(A1(nrens),A2(nrens))

   ! Old analysis step
   !
   call num2str(nanal-1,nal)
   fname = 'an'//nal//'_moderr.bin'

   ! Check if error files exist
   !
   inquire(file=fname, exist=bfile)
   if (.not. bfile) goto 777

   ! Add the old error to the new one
   write(*,*) 'Loading model error from files'
   open(22,file=fname,status='old',form='unformatted')
   read(22) A1
   close(22)

   do ne = 1,nrens
      A1(ne) = states_real_mult(A1(ne),aalpha)
      mfact = sqrt(1 - alpha**2)
      A2(ne) = states_real_mult(qA(ne),mfact)
      qA(ne) = A1(ne) + A2(ne)
   end do

   deallocate(A1,A2)

   return

 777 continue

   write(*,*) 'Model error files not found'

  end subroutine load_error


end module mod_mod_err
