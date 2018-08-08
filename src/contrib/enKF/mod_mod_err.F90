!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
!
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

 type(qstates), allocatable :: Ashy_aug(:) 	! double state with model error

 type(states), allocatable, private :: qA(:)		! model error

contains

!********************************************************

  subroutine info_moderr
  
  implicit none

   open(22, file='mod_err.info', status='old')
   read(22,*) nx_er,ny_er,fmult_er,theta_er,rsigma,dt_er,tau_er
   close(22)

  end subroutine info_moderr

!********************************************************

  subroutine push_aug
  
   !use basin
   
   implicit none
   integer nst 		!= nanal-1	number of time steps from the begin of the assimilation
   double precision alpha,rho

   type(states), allocatable :: Aaux
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
      qA(ne) = 0.
      qA(ne)%z = kvec(:,ne)
   end do
 
   !---------------------------------------
   !if exist old error q0 load it and add to qA (q1 = alpha*q0 + sqrt(1-alpha**2)*w)
   !---------------------------------------
   call load_error(alpha)
 
   !---------------------------------------
   !compute the new state qA = Ashy + sqrt(dt)*sigma*rho*q1
   !---------------------------------------
   mfact = sqrt(dt_er) * rsigma * rho
   allocate(Aaux)
   do ne = 1,nrens
      ! find the spatial factor from the relative error (rsigma)
      Aaux = Ashy(ne) * mfact
      ! find the error
      Aaux = qA(ne) * Aaux
      ! add the error
      Ashy(ne) = Ashy(ne) + Aaux
   end do
   deallocate(Aaux)
    
   !---------------------------------------
   !make the augmented state Ashy_aug = (Ashy,qA)
   !---------------------------------------
   allocate(Ashy_aug(nrens))
   do ne = 1,nrens
      call push_qstate(Ashy(ne),qA(ne),Ashy_aug(ne))
   end do
   deallocate(qA)
   deallocate(Ashy)
 
  end subroutine push_aug

!********************************************************

  subroutine pull_aug
  
  implicit none

  type(states),allocatable :: qA(:)
  character(len=3) :: nrel,nal
  character(len=16) fname
  integer ne

  write(*,*) 'Saving model errors'

  allocate(Ashy(nrens),qA(nrens))
  do ne = 1,nrens
     call pull_qstate(Ashy(ne),qA(ne),Ashy_aug(ne))
  end do

  call num2str(nanal,nal)
  fname='an'//nal//'_moderr.bin'
  open(33,file=fname,form='unformatted')
  write(33) qA
  close(33)

  deallocate(Ashy_aug,qA)

  end subroutine pull_aug

!********************************************************

  subroutine load_error(alpha)
  
   implicit none
   double precision, intent(in) :: alpha

   logical :: bfile
   integer ne
   character(len=3) :: nrel,nal
   character(len=16) fname
   type(states),allocatable :: qA1(:),A2

   real aalpha

   aalpha = alpha

   allocate(qA1(nrens),A2)

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
   read(22) qA1
   close(22)

   do ne = 1,nrens
      A2 = qA1(ne) * aalpha
      qA(ne) = A2 + (qA(ne) * sqrt(1 - alpha**2))
   end do

   deallocate(qA1,A2)

   return

 777 continue

   write(*,*) 'Model error files not found'

  end subroutine load_error


end module mod_mod_err
