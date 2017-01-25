  module mod_enkf

  use mod_states
  use mod_observations

  implicit none

  character (len=80), save :: basfile, obsfile, rstfile
  integer, save :: is_new_ens, is_mod_err

  integer, save :: nrens, na

  ! Observations
  double precision, save :: tobs
  integer, save :: nobs_lev,nobs_tot
  type(levels), save :: olev

  type(states), allocatable, save  :: A(:) 		! ens states
  type(states), save  :: Am		! mean state

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
  
  if( mod(nrens,2).eq.0 ) stop 'read_info: n of ens members must be odd, with the control as first.'

  ! Allocates the type A to store the ens states
  allocate(A(nrens))

  close(20)

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
  Am = states_real_mult(Am,1./nrens)

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

  subroutine make_D_E_R

! R (only used if mode=?1 or ?2) (no for low-rank sq root)

  implicit none
  real rand_v(nrens)
  integer n,ne

  allocate(D(nobs_tot,nrens),E(nobs_tot,nrens),R(nobs_tot,nobs_tot))

  call read_obs

  R = 0.	!Observations are indipendent
  do n = 1,nobs_lev

     ! Makes a random vector
     call random2(rand_v,nrens)
     rand_v = rand_v - (sum(rand_v)/nrens)

     E(n,:) = olev%std(n) * rand_v(:)
     D(n,:) = E(n,:) + olev%val(n)

     R(n,n) = olev%std(n)**2
  end do

  end subroutine make_D_E_R

!********************************************************

  subroutine make_S_innov

  implicit none

  integer iel
  integer n
  double precision av_mod,inn
  real*4 x4,y4
  double precision sk
  integer ne, i

  allocate (S(nobs_tot,nrens),innov(nobs_tot))
  
  if( nobs_lev.gt.0 ) then
    do n = 1,nobs_lev

       ! Finds the nearest element
       x4 = olev%x(n)
       y4 = olev%y(n)
       call find_element(x4,y4,iel)

         do ne = 1,nrens
          sk = 0.
          do i = 1,3
             sk = sk + ( A(ne)%ze(i,iel) - Am%ze(i,iel) )
          end do
          S(n,ne) = sk/3.
         end do

         av_mod = sum(Am%ze(:,iel))/3.
         inn = olev%val(n) - av_mod
     
         call check_innov_val(inn,'level')

         innov(n) = inn
         write(*,'(a,i5,3f8.4)') ' nobs, vobs, vmod, innov: ',n,olev%val(n),av_mod,inn

    end do
  end if

  end subroutine make_S_innov

!********************************************************

  subroutine read_ensemble

   implicit none

   integer ne

   type(states4) :: A4
   type(states) :: Ap(nrens-1)

   character(len=3) :: nrel,nal
   character(len=16) rstname

   call num2str(na,nal)

   ! reads the ensemble or makes a new one
   if( (is_new_ens.eq.0).or.(na.ne.0) ) then

     do ne = 1,nrens
        call num2str(ne-1,nrel)
        rstname='an'//nal//'_'//'en'//nrel//'b.rst'
        call rst_read(nnkn,nnel,nnlv,rstname,tobs)
        call push_state(A4)
        A(ne)=A4
     end do

   else if( (is_new_ens.eq.1).and.(na.eq.0) ) then

     ! reads just the control member
     call num2str(0,nrel)
     rstname='an'//nal//'_'//'en'//nrel//'b.rst'
     call push_state(A4)
     A(1)=A4

     ! makes the perturbed members
     call make_pert(Ap,nrens-1) !TODO
     do ne = 2,nrens
        A(ne) = A(1) + Ap(ne-1)
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

   type(states4) A4

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


  end module mod_enKF
