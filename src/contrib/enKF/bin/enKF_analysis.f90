!------------------------------------------------------------------------------
   program enKF_analysis
!------------------------------------------------------------------------------
! Load the restart files Make the enKF analysis and save the analysis states in
! output restart files
!------------------------------------------------------------------------------
   implicit none

!---------------------
! include for restart
   include 'param.h'
   include 'nlevel.h'
   include 'hydro.h'
   include 'geom_dynamic.h'
   include 'basin.h'
   include 'ts.h'
   include 'hydro_vel.h'
   include 'conz.h'

   double precision tt

!---------------------
! asspar parameters
   character assdir*40
   integer nrens
   integer date, time
   double precision tt_anf,tt_eanf,tt_oanf,tt_oend,tt_eend,tt_end

!--- Restart
   integer itrst
   character rstname*70
   double precision atime

!--- Observation vars
   integer, parameter :: nrmax = 1000	! Maximum number of observations
   integer nrobs			! Number of observations
   real val_obs(nrmax), std_obs(nrmax)
   real x_obs(nrmax), y_obs(nrmax), z_obs(nrmax)
   character (len = 10) type_obs(nrmax)

!--- Dimensions
   integer ndim             ! Dimension of model states

!---------------------
! Vars for the analysis routine
   real, allocatable :: A(:,:)!ndim,nrens)    ! ensemble matrix
   real, allocatable :: R(:,:)!nrobs,nrobs)   ! matrix holding R (only used if mode=?1 or ?2)
   real, allocatable :: D(:,:)!nrobs,nrens)   ! matrix holding perturbed measurments
   real, allocatable :: E(:,:)!nrobs,nrens)   ! matrix holding perturbations (only used if mode=?3)
   real, allocatable :: S(:,:)!nrobs,nrens)   ! matrix holding HA` 
   real, allocatable :: innov(:)!nrobs)     ! vector holding d-H*mean(A)

   logical verbose          ! Printing some diagnostic output

   real truncation       ! The ratio of variaince retained in pseudo inversion (0.99)

   integer mode             ! first integer means (EnKF=1, SQRT=2)
                                           ! Second integer is pseudo inversion
                                           !  1=eigen value pseudo inversion of SS'+(N-1)R
                                           !  2=SVD subspace pseudo inversion of SS'+(N-1)R
                                           !  3=SVD subspace pseudo inversion of SS'+EE'

   logical update_randrot   ! Normally true; false for all but first grid point
                                           ! updates when using local analysis since all grid
                                           ! points need to use the same rotation.
!---------------------
   real, allocatable :: Am(:) !ndim    ! average state

   integer ne
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

!----------------------------------    
! Reads parameters
!----------------------------------    
    read(5,*) tt
    call read_asspar(23,assdir,nrens,date,time,tt_anf,tt_end,tt_eanf,tt_eend,tt_oanf,tt_oend)
    nlv = 8 !restart read needs it

    ! Analysis parameters
    mode = 23
    verbose = .true.
    truncation = 0.99
    update_randrot = .true.

!--------------------------------
! Read basin
!--------------------------------
    open(21,file=trim(assdir)//'/nador2006_3321.bas',status='old',form='unformatted')
    call sp13rr(21,nkndim,neldim)
    close(21)

!-----------------------------------
! Load restart files and fill A(ndim,nrens)
!-----------------------------------
! rstfile_list.txt : list of restart files
  open(22,file=trim(assdir)//'/rstfile_list.txt', status='old')
     
  do ne = 1,nrens

     read(22,'(a70)') rstname
     call rst_read(23,rstname,date,time,tt)

     !----------
     ! Here we defined the model state to be analysed and store it
     ! into a matrix for each ensemble member
     if (ne.eq.1) then
	ndim = nkn + 2 * nel * nlv
	allocate(A(ndim,nrens))
     end if

     ! Store A
     A(:,ne) = [ znv(1:nkn), utlnv(1:nlv,1:nel), vtlnv(1:nlv,1:nel) ]
     !----------
       
  end do
  close(22)

!--------------------------------
! Make mean(A) matrix holding the average state
!--------------------------------
   write(*,*) 'Computing the mean value Am of the ensemble A'
   allocate(Am(ndim))
   call average_mat(A,Am,ndim,nrens)
   
!--------------------------------
! Make D(nrobs,nrens)   matrix holding perturbed measurments
!--------------------------------
   ! load original observations
   call read_obs(tt,22,trim(assdir)//'/observations.bin',nrobs,type_obs,x_obs,y_obs,z_obs,val_obs,std_obs)
   allocate(D(nrobs,nrens))
   call make_random_D(nrobs,nrens,val_obs,std_obs,D)

!--------------------------------
! Make S(nrobs,nrens), matrix holding HA`, and innov(nrobs), vector holding d-H*mean(A)
! TODO: add more model vars
!--------------------------------
   allocate(S(nrobs,nrens),innov(nrobs))
   call make_S_innov(nrens,nrobs,ndim,A,Am,val_obs,x_obs,y_obs,z_obs,type_obs,S,innov)

!--------------------------------
! Make R(nrobs,nrobs)  matrix holding R (only used if mode=?1 or ?2) (no for low-rank sq root)
!--------------------------------
   allocate(R(nrobs,nrobs))
   R = 0

!--------------------------------
! Make E(nrobs,nrens)  matrix holding perturbations (only used if mode=?3) (no for low-rank sq root)
!--------------------------------
   allocate(E(nrobs,nrens))
   E = 0

!--------------------------------
! Call the analysis routine
!--------------------------------
   call analysis(A, R, E, S, D, innov, ndim, nrens, nrobs, verbose, truncation,mode,update_randrot)
   write(*,*) 'Analysis done'


!--------------------------------
! Save the output in different restart files
!--------------------------------

!! To write a restart
!    open(30,file=trim(assdir)//'test.rst',form='unformatted')
!    call addpar('ibarcl',1)
!    call addpar('iconz',1)
!    call addpar('ibfm',0)
!    call addpar('date','20040101')
!    call addpar('time','000000')
!    call wrrst(itrst,30)
!    close(30)

end



!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!					SUBROUTINES
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------

!------------------------------------------------------
   subroutine average_mat(M,Mmean,ni,nj)
!------------------------------------------------------
   implicit none
   integer, intent(in) :: ni, nj
   real, intent(in) :: M(ni,nj)
   real, intent(out) :: Mmean(ni)
   integer i

   do i = 1,ni
      Mmean(i) = sum(M(i,1:nj)) / nj
   end do

   end subroutine

!------------------------------------------------------
   subroutine make_random_D(nrobs,nrens,obs,ostd,D)
!------------------------------------------------------
   implicit none
   integer, intent(in) :: nrobs,nrens
   real, intent(in) :: obs(nrobs),ostd(nrobs)
   real, intent(inout) :: D(nrobs,nrens)
   real work1(nrens)
   integer n

   do n = 1,nrobs
      call random2(work1,nrens)
      D(n,:) = obs(n) + ( ostd(n) * work1 )
   end do

   end subroutine

!------------------------------------------------------
   subroutine random2(work1,n)
!------------------------------------------------------
!  From Evenseen, slightly modified
!  Returns a vector of random values N(variance=1,mean=0)
   !use mod_dimensions
   implicit none
   integer, intent(in) :: n
   real, intent(out) :: work1(n)
   real, allocatable :: work2(:)
   real, parameter :: pi = 4. * atan(1.)

   allocate (work2(n))

   call random_number(work1)
   call random_number(work2)
   work1= sqrt(-2.0*log(work1))*cos(2.0*pi*work2)

   deallocate(work2)
   end subroutine random2


!------------------------------------------------------
   subroutine make_S_innov(nrens,nrobs,ndim,A,Am,vobs,x,y,z,otype,S,innov)
!------------------------------------------------------
! Finds S = HA` : the model perturbations interpolated in the observation
! locations
! Remember: C = SS^t + EE^t 
! with SS^t the cov matrix of model in the obs space
   implicit none

!---------------------
! include
   include 'param.h'
   include 'nlevel.h'
   include 'hydro.h'
   include 'geom_dynamic.h'
   include 'basin.h'
   include 'ts.h'
   include 'hydro_vel.h'
   include 'conz.h'

   integer, intent(in) :: nrens,nrobs,ndim
   real, intent(in) :: A(ndim,nrens)
   real, intent(in) :: Am(ndim)
   real, intent(in) :: vobs(nrobs)
   real, intent(in) :: x(nrobs),y(nrobs),z(nrobs)    !position of the observations
   character (len=*), intent(in) :: otype(nrobs)     !Obs type
   real, intent(out) :: S(nrobs,nrens)
   real, intent(out) :: innov(nrobs)
   
   real xe(3),ye(3) 
   integer kn(3)

   integer intri,inval
   logical istri	! Model variable corresponding to obs is defined in element (vel)
   logical isvert	! Model variable corresponding to obs has a vertical dimension

   integer no,ne
   integer ie,ii
   real val
   real kw(3)	! spatial weights for linear interpolation

   ! Obs in the timestep
   do no = 1,nrobs 

      inval = 0
      do ie = 1,nel
         do ii = 1,3
                kn(ii) = nen3v(ii,ie)
                xe(ii) = xgv(kn(ii))
                ye(ii) = ygv(kn(ii))
         end do
         !write(*,*) xe,ye,x(no),y(no) 	!obs coords are wrong TODO
         inval = intri(xe,ye,x(no),y(no))
         if( inval.eq.1 ) exit
      end do
      if( inval.eq.0 ) stop 'Observation outside the domain' 
      

      if( otype(no)(7:9).eq.'vel' ) then 
          istri = .true. !value on the elements
	  write(*,*) 'TODO, not yet implemented'
          stop
      else
          istri = .false.
          call int_weight(xe,ye,x(no),y(no),kw)
      end if
      
      ! Assign value depending on model var
      if( otype(no).eq.'times_lev' ) then
          isvert = .false.
          do ne = 1,nrens
             val = kw(1)*A(kn(1),ne) + kw(2)*A(kn(2),ne) + kw(3)*A(kn(3),ne)
             S(no,ne) = val
          end do
          val = kw(1)*Am(kn(1)) + kw(2)*Am(kn(2)) + kw(3)*Am(kn(3))
          innov(no) = vobs(no) - val
      else
          isvert = .true.
	  write(*,*) 'TODO, not yet implemented'
          stop
      end if

   end do

   end subroutine

!------------------------------------------------------
   subroutine int_weight(xe,ye,xo,yo,d)
!------------------------------------------------------
   implicit none
   real, intent(in) :: xe(3),ye(3)
   real, intent(in) :: xo,yo
   real, intent(out) :: d(3)

   integer ii
   real dtot

   do ii = 1,3
      d(ii) = (xe(ii) - xo)**2. + (ye(ii) - yo)**2.
   end do

   dtot = sum(d)

   d = d / dtot

   end subroutine

!------------------------------------------------------
   subroutine rst_read(iunit,filin,date,time,tt)
!------------------------------------------------------
   implicit none
!---------------------
! include for restart
   include 'param.h'
   include 'nlevel.h'
   include 'hydro.h'
   include 'geom_dynamic.h'
   include 'basin.h'
   include 'ts.h'
   include 'hydro_vel.h'
   include 'conz.h'
!---------------------
   integer iunit
   character(len=*), intent(in) :: filin
   double precision tt, atime
   integer date, time, it
   integer ierr

   it = nint(tt)

   open(iunit,file=filin,status='old',form='unformatted')
   call dts_to_abs_time(date,time,atime)
   atime = atime + it
   write(*,*) 'Reading restart file: ',filin
   call rdrst(atime,iunit,ierr)
   close(iunit)

   end subroutine
