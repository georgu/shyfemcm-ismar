!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
!
!------------------------------------------------------------------------------
! Convert a enKF simulation saved in rst format into a enKS analysis.
!------------------------------------------------------------------------------
! This program needs a rst file with the background of the KF run and
! the "X5_tot.uf" file, which contains the X5 matrix at each timestep
!------------------------------------------------------------------------------
program enKF2enKS

  use basin
  use levels, only : nlv
  implicit none

  character(len=4) :: arg1
  character(len=80) :: arg2
  character(len=9) :: arg3
  character(len=4) :: arg4,arg5
  character(len=80) :: enfile
  integer nrens,nre
  character(len=5) :: nrel
  integer ierr
  double precision atime
  integer fid
  real, allocatable :: Astate(:,:),AmeanKF(:),AstdKF(:),AmeanKS(:),AstdKS(:)
  integer rrec,sdim
  integer nlag

  ! read number of ens member from input
  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)
  call get_command_argument(4, arg4)
  call get_command_argument(5, arg5)
  if ((trim(arg2) /= 'norm').and.(trim(arg2) /= 'full')) then
	  write(*,*) ''
	  write(*,*) 'Usage: enKF2enKS [nrens] [output] [nlag] [file-type] [ks-option]'
	  write(*,*) ''
	  write(*,*) '[nrens] is the n. of ens members, control included.'
	  write(*,*) '[output] full (full) or just mean and std (norm).'
	  write(*,*) '[nlag] number of forward analysis steps to consider (-1 to consider all)'
	  write(*,*) '[file_type] can be back or anal'
	  write(*,*) '[ks-option] can be noks or ks, to make the KS analysis'
	  write(*,*) ''
	  stop
  end if
  read(arg1,*) nrens
  if (nrens <= 0) error stop 'enKF2enKS: not a valid nrens number'

  read(arg3,*) nlag
  write(*,*) 'Lagged smoother. N. of analysis steps: ',nlag

  !--- init variables
  call init_shyfem
  call add_rst_params

  !--- open rst files  
  do nre = 1,nrens
     fid = 20 + nre
     call num2str(nre-1,nrel)
     if (arg4 == 'back') then
        enfile = 'backKF_en'//nrel//'.rst'
     else
        enfile = 'analKF_en'//nrel//'.rst'
     end if
     open(fid,file=trim(enfile),status='old',form='unformatted',iostat=ierr)
     if (ierr /= 0) error stop 'enKF2enKS: error opening file'
  end do

  if ((trim(arg2) == 'full') .and. (trim(arg5) == 'ks')) then
    write(*,*) 'Writing the full output'
    do nre = 1,nrens
       fid = 20 + nrens + nre
       call num2str(nre-1,nrel)
       enfile = arg4//'KS_en'//nrel//'.rst'
       open(fid,file=trim(enfile),status='unknown',form='unformatted')
    end do
  end if

  open(16,file=arg4//'KF_mean.rst',status='unknown',form='unformatted')
  open(17,file=arg4//'KF_std.rst',status='unknown',form='unformatted')


  if (trim(arg5) == 'ks') then
    open(18,file=arg4//'KS_mean.rst',status='unknown',form='unformatted')
    open(19,file=arg4//'KS_std.rst',status='unknown',form='unformatted')
  end if


  !--------------------------- loop on time records
  rrec = 0
 89 continue  

  !--- read ensemble
  do nre = 1,nrens
     fid = 20 + nre
     call rst_read_record(fid,atime,ierr)
     if (ierr == -1) goto 90
     if (ierr > 0) goto 91

     if (.not. allocated(Astate)) then
	 sdim = nkn + 2*nlv*nel + 2*nlv*nkn
	 allocate(Astate(sdim,nrens))
	 allocate(AmeanKF(sdim),AmeanKS(sdim))
	 allocate(AstdKF(sdim),AstdKS(sdim))
	 Astate = 0.0
     end if

     call push_matrix(sdim,nrens,nre,Astate)
  end do



  !--- make mean and std of the ensemble Kalman Filter
  call make_mn_std(sdim,nrens,Astate,AmeanKF,AstdKF)

  call pull_state(sdim,AmeanKF)
  call rst_write_rec(atime,16)
  call pull_state(sdim,AstdKF)
  call rst_write_rec(atime,17)
 
  if (trim(arg5) == 'ks') then

    !--- make analysis
    ! Results differ for a very small epsilon. 2nd method is much faster.
    !call make_analysis1(atime,sdim,nrens,Astate)
    call make_analysis2(atime,sdim,nrens,Astate,nlag,arg4)

    !--- make mean and std of the ensemble Kalman Smoother
    call make_mn_std(sdim,nrens,Astate,AmeanKS,AstdKS)

    call pull_state(sdim,AmeanKS)
    call rst_write_rec(atime,18)
    call pull_state(sdim,AstdKS)
    call rst_write_rec(atime,19)

    !--- write ensemble 
    if (trim(arg2) == 'full') then
      do nre = 1,nrens
         fid = 20 + nrens + nre
         call pull_matrix(sdim,nrens,nre,Astate)
         call rst_write_rec(atime,fid)
      end do
    end if

  end if

  rrec = rrec + 1
  write(*,*) 'N. of record: ',rrec

 goto 89
  !--------------------------- end loop on time records

  !--- close rst files
 90 continue
  do nre = 1,nrens
     fid = 20 + nre
     close(fid)
  end do

  close(16)
  close(17)

  if (trim(arg5) == 'ks') then
    close(18)
    close(19)
    if (trim(arg2) == 'full') then
      do nre = 1,nrens
         fid = 20 + nrens + nre
         close(fid)
      end do
    end if
  end if

  stop

 91 error stop 'enKF2enKS error reading file'

end program enKF2enKS


!*******************************************************************
!*******************************************************************

!********************************************************

subroutine num2str(num,str)

  implicit none

  integer, intent(in) :: num
  character(len=5), intent(out) :: str

  if ((num >= 0).and.(num < 10)) then
    write(str,'(a4,i1)') '0000',num
  elseif ((num >= 10).and.(num < 100)) then
    write(str,'(a3,i2)') '000',num
  elseif ((num >= 100).and.(num < 1000)) then
    write(str,'(a2,i3)') '00',num
  elseif ((num >= 1000).and.(num < 10000)) then
    write(str,'(a1,i4)') '0',num
  elseif ((num >= 10000).and.(num < 100000)) then
    write(str,'(i5)') num
  else
    error stop 'num2str: num out of range'
  end if

end subroutine num2str

!*******************************************************************

subroutine push_matrix(sdim,nrens,nre,Amat)
  ! the sequence is: u,v,z,t,s
  use basin
  use levels, only : nlv
  use mod_hydro
  use mod_hydro_vel
  use mod_ts
  use mod_conz

  implicit none

  integer, intent(in) :: sdim,nrens,nre
  real, intent(inout) :: Amat(sdim,nrens)

  integer dimz,dimuv,dimts

  dimz = nkn
  dimuv = nlv*nel
  dimts = nlv*nkn

  Amat(1:dimuv,nre) = reshape(dble(utlnv),(/dimuv/))
  Amat(dimuv+1:2*dimuv,nre) = reshape(dble(vtlnv),(/dimuv/))
  Amat(2*dimuv+1:2*dimuv+dimz,nre) = dble(znv)
  Amat(2*dimuv+dimz+1:2*dimuv+dimz+dimts,nre) = reshape(dble(tempv),(/dimts/))
  Amat(2*dimuv+dimz+dimts+1:2*dimuv+dimz+2*dimts,nre) = reshape(dble(saltv),(/dimts/))

end subroutine push_matrix

!*******************************************************************

subroutine pull_matrix(sdim,nrens,nre,Amat)
  ! the sequence is: u,v,z,t,s
  use basin
  use levels, only : nlv
  use mod_hydro
  use mod_hydro_vel
  use mod_ts
  use mod_conz

  implicit none

  integer, intent(in) :: sdim,nrens,nre
  real, intent(in) :: Amat(sdim,nrens)

  integer dimz,dimuv,dimts

  dimz = nkn
  dimuv = nlv*nel
  dimts = nlv*nkn

  utlnv = reshape(real(Amat(1:dimuv,nre)),(/nlv,nel/))
  vtlnv = reshape(real(Amat(dimuv+1:2*dimuv,nre)),(/nlv,nel/))
  znv = real(Amat(2*dimuv+1:2*dimuv+dimz,nre))
  tempv = reshape(real(Amat(2*dimuv+dimz+1:2*dimuv+dimz+dimts,nre)),(/nlv,nkn/))
  saltv = reshape(real(Amat(2*dimuv+dimz+dimts+1:2*dimuv+dimz+2*dimts,nre)),(/nlv,nkn/))

end subroutine pull_matrix

!*******************************************************************

subroutine pull_state(nvec,Avec)
  ! the sequence is: u,v,z,t,s
  use basin
  use levels, only : nlv
  use mod_hydro
  use mod_hydro_vel
  use mod_ts
  use mod_conz

  implicit none

  integer, intent(in) :: nvec
  real, intent(in) :: Avec(nvec)

  integer dimz,dimuv,dimts

  dimz = nkn
  dimuv = nlv*nel
  dimts = nlv*nkn

  utlnv = reshape(real(Avec(1:dimuv)),(/nlv,nel/))
  vtlnv = reshape(real(Avec(dimuv+1:2*dimuv)),(/nlv,nel/))
  znv = real(Avec(2*dimuv+1:2*dimuv+dimz))
  tempv = reshape(real(Avec(2*dimuv+dimz+1:2*dimuv+dimz+dimts)),(/nlv,nkn/))
  saltv = reshape(real(Avec(2*dimuv+dimz+dimts+1:2*dimuv+dimz+2*dimts)),(/nlv,nkn/))

end subroutine pull_state

!********************************************************

subroutine rst_write_rec(atimea,fid)

  use mod_restart
  use levels, only : hlv

  implicit none

  integer, intent(in) :: fid
  double precision, intent(in) :: atimea
  real*4 :: svar

  ! adds parameters
  !
  svar = 1.
  call putpar('ibarcl',svar)
  svar = 0.
  call putpar('iconz',svar)
  call putpar('ibfm',svar)
  call putpar('imerc',svar)
  call putpar('ibio',svar)

  ! In 2D barotropic hlv is set to 10000.
  !
  if (size(hlv) == 1) hlv = 10000.

  call rst_write_record(atimea,fid)

end subroutine rst_write_rec

!********************************************************

subroutine init_shyfem

  use basin
  use mod_restart
  use mod_geom_dynamic
  use levels
  use mod_hydro
  use mod_hydro_vel
  use mod_ts
  use mod_conz

  implicit none

  double precision atime
  integer nvers,nrec,iflag,ierr,iconzrst
  character(len=5) :: nrel

  call num2str(0,nrel)
  open(14,file='backKF_en'//nrel//'.rst',status='old',form='unformatted',iostat=ierr)
  if (ierr /= 0) error stop 'enKF2enKS: error opening file'
  call rst_skip_record(14,atime,nvers,nrec,nkn,nel,nlv,iconzrst,iflag,ierr)
  close(14)

  if (.not. allocated(hm3v)) allocate(hm3v(3,nel))

  call mod_geom_dynamic_init(nkn,nel)
  call mod_hydro_init(nkn,nel,nlv)
  call mod_hydro_vel_init(nkn,nel,nlv)
  call mod_ts_init(nkn,nlv)
  if (iconzrst > 0) call mod_conz_init(iconzrst,nkn,nlv)
  call levels_init(nkn,nel,nlv)

end subroutine init_shyfem

!********************************************************

  subroutine add_rst_params

  use mod_restart

  implicit none
  real*4 :: svar
  double precision :: dvar

  svar = 1.
  call addpar('ibarcl',svar)
  svar = 0.
  call addpar('iconz',svar)
  call addpar('ibfm',svar)
  call addpar('imerc',svar)
  call addpar('ibio',svar)

  dvar = 0.
  call daddpar('date',dvar)
  call daddpar('time',dvar)

  end subroutine add_rst_params

!********************************************************

  subroutine make_analysis1(atime,sdim,nrens,Amat)
  implicit none
  double precision, intent(in) :: atime
  integer, intent(in) :: sdim,nrens
  real, intent(inout) :: Amat(sdim,nrens)

  double precision tt
  character(len=6) :: alabel
  character(len=2) :: tag
  integer nren
  real, allocatable :: Aaux(:,:),X5(:,:)

    allocate(Aaux(sdim,nrens))
    allocate(X5(nrens,nrens))


    open(15,status='old',form='unformatted',file='X5_tot.uf')
 45 read(15,end=44) tt,alabel,tag

     read(15) nren

     if ((tag == 'X3').or.(nren /= nrens)) then
        error stop 'Error reading X5_tot.uf.'
     end if

     read(15) X5

     Aaux = 0.0
     if (tt >= atime) then
        !Aaux = matmul(Amat,X5)
	! This way should be fine (as used in Geir's code)
	call dgemm('n','n',sdim,nrens,nrens,1.0,Amat,sdim,X5,nrens,0.0,Aaux,sdim)
	Amat = Aaux
     end if

     goto 45

 44 close(15)

  deallocate(Aaux,X5)

  end subroutine make_analysis1

!********************************************************

  subroutine make_analysis2(atime,sdim,nrens,Amat,nlag,ftype)
  implicit none
  double precision, intent(in) :: atime
  integer, intent(in) :: sdim,nrens
  real, intent(inout) :: Amat(sdim,nrens)
  integer, intent(inout) :: nlag
  character(len=4), intent(in) :: ftype

  double precision tt
  character(len=6) :: alabel
  character(len=2) :: tag
  integer nren
  real, allocatable :: Aaux(:,:),X5(:,:),X5old(:,:),X5aux(:,:)
  integer n,k

    allocate(Aaux(sdim,nrens))
    allocate(X5(nrens,nrens),X5old(nrens,nrens),X5aux(nrens,nrens))

    if (nlag == -1) nlag = 1000000000

    ! initial set the identity
    X5old = 0.
    forall(n = 1:nrens) X5old(n,n) = 1. 

    open(15,status='old',form='unformatted',file='X5_tot.uf')
    k = 0
 45 read(15,end=44) tt,alabel,tag

     read(15) nren

     if ((tag == 'X3').or.(nren /= nrens)) then
        error stop 'X3 not implemented or bad n. of ens members'
     end if

     read(15) X5

     if (ftype == 'anal') then
      if ((tt > atime).and.(k <= nlag)) then	!This should be correct using the analysis rst
	X5aux = 0.
	! This way should be fine (as used in Geir's code)
 	call dgemm('n','n',nrens,nrens,nrens,1.0,X5old,nrens,X5,nrens,0.0,X5aux,nrens)
	X5old = X5aux
      end if
     else
      if ((tt >= atime).and.(k <= nlag)) then	!This should be correct using the back rst
	X5aux = 0.
	! This way should be fine (as used in Geir's code)
 	call dgemm('n','n',nrens,nrens,nrens,1.0,X5old,nrens,X5,nrens,0.0,X5aux,nrens)
	X5old = X5aux
      end if
     end if

     k = k + 1

     goto 45

 44 close(15)

    Aaux = 0.0
    call dgemm('n','n',sdim,nrens,nrens,1.0,Amat,sdim,X5old,nrens,0.0,Aaux,sdim)
    Amat = Aaux

    deallocate(Aaux,X5,X5old,X5aux)

  end subroutine make_analysis2


!********************************************************

  subroutine make_mn_std(ndim,nens,Amat,Am,Astd)
  implicit none

  integer, intent(in) :: ndim,nens
  real, intent(in) :: Amat(ndim,nens)
  real, intent(out) :: Am(ndim),Astd(ndim)
  real esum,esumsq
  integer n

  if (nens < 2) error stop 'Ensemble too small'

  do n = 1,ndim
     esum = sum(Amat(n,:))
     esumsq = sum(Amat(n,:)*Amat(n,:))
     Am(n) = esum / dble(nens)
     Astd(n) = sqrt( (esumsq - esum*esum/dble(nens)) / (dble(nens) - 1.) )
  end do

  end subroutine make_mn_std

