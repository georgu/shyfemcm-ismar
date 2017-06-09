	program make_ens_meteo
	implicit none

	integer iunit, iformat
	integer ounit
	character(len=80) :: filein
	double precision dtime  !time stamp
	integer nvers           !version of file format
	integer np              !size of data (horizontal, nodes or elements)
	integer nvar            !number of variables to write
        integer ntype           !type of information contained
        integer datetime(2)     !date and time information
        integer ierr		!return error code
	real*4 regpar(7)          !regular array params
	integer lmax		!vertical values
	real*4,allocatable :: hlv(:)          !vertical structure
	character*50 string
	integer,allocatable :: ilhkv(:)
	real*4,allocatable :: hd(:)
	integer nlvddi
	real*4,allocatable :: data(:,:),dataens(:,:)
	real rerr
	integer nrens
	integer irec
	integer i,nr

	integer nx,ny
	real dx,dy
	real,allocatable :: pmat1(:,:,:), pmat2(:,:,:)

	nrens = 30
	rerr = .3	!30% of relative error

	filein = 'wind.fem'
	np = 0
	call fem_file_read_open(filein,np,iformat,iunit)

        !--------------------------------------------------
	! Reads and writes records
        !--------------------------------------------------
	irec = 0
	do

	  irec = irec + 1
          !--------------------------------------------------
	  ! Reads headers
          !--------------------------------------------------
 	  call fem_file_read_params(iformat,iunit,dtime
     +                          ,nvers,np,lmax,nvar,ntype,datetime,ierr)
	  if( ierr .lt. 0 ) exit
	  write(*,*) 'Time: ',dtime
	  allocate(hlv(lmax))
	  nlvddi = lmax
	  call fem_file_read_2header(iformat,iunit,ntype,lmax
     +                  ,hlv,regpar,ierr)
	  nx = nint(regpar(1))
	  ny = nint(regpar(2))
	  allocate(ilhkv(np),hd(np),data(nx,ny))
	  allocate(dataens(nx,ny))

          !--------------------------------------------------
	  ! Makes 2D perturbed normalised fields and update them
          !--------------------------------------------------
	  dx = regpar(5)
	  dy = regpar(6)
	  if(.not.allocated(pmat1)) allocate(pmat1(nx,ny,nrens))
	  if(.not.allocated(pmat2)) allocate(pmat2(nx,ny,nrens))
	  call make_pert(irec,dtime,nrens,nx,ny,dx,dy,pmat1,pmat2)
          !--------------------------------------------------

          !--------------------------------------------------
	  ! Writes headers
          !--------------------------------------------------
	  do nr = 1,nrens

	     ounit = iunit + 10 + nr
	     call fem_file_write_header(iformat,ounit,dtime
     +                          ,nvers,np,lmax
     +                          ,nvar,ntype
     +                          ,nlvddi,hlv,datetime,regpar)

	  end do


          !--------------------------------------------------
	  ! Loop on variables
          !--------------------------------------------------
 	  do i = 1,nvar

             !--------------------------------------------------
	     ! Reads variable
             !--------------------------------------------------
	     call fem_file_read_data(iformat,iunit
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvddi,data
     +                          ,ierr)
	     write(*,*) 'Reading: ',string
	  
	     do nr = 1,nrens
                !--------------------------------------------------
	        ! Makes 2D perturbed real wind fields
                !--------------------------------------------------
		select case (i)
		 case (1)
		   dataens = data + (rerr * data) * pmat1(:,:,nr)
		   !dataens = pmat1(:,:,nr)	!check rand field
		 case (2)
		   dataens = data + (rerr * data) * pmat2(:,:,nr)
		   !dataens = pmat2(:,:,nr)	!check rand field
		 case (3)
		   dataens = data
		end select
                !--------------------------------------------------

                !--------------------------------------------------
	        ! Writes record
                !--------------------------------------------------
	        ounit = iunit + 10 + nr
	        call fem_file_write_data(iformat,ounit
     +                          ,nvers,np,lmax
     +                          ,string
     +                          ,ilhkv,hd
     +                          ,nlvddi,dataens)
                !--------------------------------------------------
	     end do

	  end do	! loop on variables
          !--------------------------------------------------

	  deallocate(ilhkv,hd,data,dataens)
	  deallocate(hlv)

	end do		!loop on records
        !--------------------------------------------------

	! closes files
	close(iunit)
	do nr = 1,nrens
	   ounit = iunit + 10 + nr
	   close(ounit)
        end do

	end program make_ens_meteo

!***********************************************************
!***********************************************************
!***********************************************************
        !--------------------------------------------------
	subroutine make_pert(irec,tt,nrens,nx,ny,dx,dy,pmat1,pmat2)
        !--------------------------------------------------
	use m_sample2D
	implicit none
	integer, intent(in) :: irec,nrens,nx,ny
	double precision, intent(in) :: tt
	real, intent(in) :: dx,dy
	real, intent(out) :: pmat1(nx,ny,nrens), pmat2(nx,ny,nrens)
	real, allocatable, save :: pmat1_old(:,:,:),pmat2_old(:,:,:)
	double precision, save :: tt_old
	real rx,ry
	integer fmult	!Mult factor for the super-sampling
	real theta	!Rotation of the random fields (theta=0 is east, rotation anticlocwise)
	logical verbose,samp_fix
	double precision dt_er        !time between 2 analysis steps
	double precision tau_er       !time decorrelation (>=dt) of the old error: dq/dt = -(1/tau)*q
	double precision alpha

	tau_er = 86400*3	! 3 days
	theta = 0
	fmult = 8
	rx = 30. * dx
	ry = 30. * dy
	rx = rx/sqrt(3.0)
	ry = ry/sqrt(3.0)
	verbose = .true.
	samp_fix = .true.     !keep true

	! Allocate old matrices
	if(.not.allocated(pmat1_old)) allocate(pmat1_old(nx,ny,nrens))
	if(.not.allocated(pmat2_old)) allocate(pmat2_old(nx,ny,nrens))

	! Make the random fields
	call sample2D(pmat1,nx,ny,nrens,fmult,dx,dy,rx,ry,theta
     +               ,samp_fix,verbose)

	call sample2D(pmat2,nx,ny,nrens,fmult,dx,dy,rx,ry,theta
     +               ,samp_fix,verbose)

	! Merge with the old fields in order to have time correlation
        if( irec.ne.1 ) then
	  
	  dt_er = tt - tt_old
	  alpha = 1. - (dt_er/tau_er)
	  pmat1 = alpha * pmat1_old + sqrt(1 - alpha**2) * pmat1
	  pmat2 = alpha * pmat2_old + sqrt(1 - alpha**2) * pmat2

	endif
         
	! Update old fields
        pmat1_old = pmat1
        pmat2_old = pmat2
	tt_old = tt

	end subroutine make_pert
