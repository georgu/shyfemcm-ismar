!------------------------------------
! make a list of times with at least 1 observation
!------------------------------------
program make_antime_list

 use iso8601
 implicit none
 character(len=20) :: date1,date2
 character(len=80) :: filelist
 character(len=8) :: dtstr
 double precision :: dt
 integer :: d,t
 integer :: ierr
 double precision :: atime1,atime2,atime
 integer :: n,nsteps,nn
 integer :: k,ktot
 character(len=80) :: line
 character(len=5) :: oflag
 character(len=80) :: ofile
 character(len=20) :: dstr

 integer :: nr
 double precision,allocatable :: atval(:),val(:)
 integer, allocatable :: nrec(:), isfile(:)
 double precision,allocatable :: atval_tot(:,:),val_tot(:,:)
 double precision, parameter :: flag = -999.
 integer :: iok,iflag

 call get_command_argument(1, date1)
 call get_command_argument(2, date2)
 call get_command_argument(3, dtstr)
 call get_command_argument(4, filelist)

 if (trim(filelist) .eq. '') then
      write(*,*) ''
      write(*,*) 'Make a list of times with at least one observation'
      write(*,*) ''
      write(*,*) 'Usage:'
      write(*,*) ''
      write(*,*) './make_antime_list [date1] [date2] [dt] [filelist]'
      write(*,*) ''
      write(*,*) 'date1: initial date (yyyy-mm-dd::HH:MM:SS)'
      write(*,*) 'date2: final date (yyyy-mm-dd::HH:MM:SS)'
      write(*,*) 'dt: minimum timestep between observations (sec)'
      write(*,*) 'filelist: list of observation files'
      write(*,*) ''
      stop
 end if

 ! date1
 call string2date(trim(date1),d,t,ierr)
 if ( ierr /= 0 ) error stop 'Invalid date'
 call dts_to_abs_time(d,t,atime1)

 ! date2
 call string2date(trim(date2),d,t,ierr)
 if ( ierr /= 0 ) error stop 'Invalid date'
 call dts_to_abs_time(d,t,atime2)

 ! time step
 read(dtstr,*,iostat=ierr) dt
 if ( ierr /= 0 ) error stop 'Invalid timestep'
 if (( dt <= 0 ) .or. ( dt > 864000 )) error stop 'Bad timestep'

 if ( atime2 - atime1 < 2*dt ) error stop 'Bad times'

 nsteps = nint((atime2-atime1)/dt + 1)

! Read a list of files and store observations
 allocate(atval(nsteps),val(nsteps))
 k = 0
 open(20,file=trim(filelist),iostat=ierr,form='formatted',status='old')
 if ( ierr /= 0 ) error stop 'Error opening the file list'
 85  read(20,*,end=95) line
 k = k + 1
 goto 85
 95  rewind(20)
 ktot = k

 allocate(nrec(ktot),atval_tot(ktot,nsteps),val_tot(ktot,nsteps))
 allocate(isfile(ktot))
 atval_tot = -999.
 val_tot = -999.
 do k = 1,ktot

    read(20,*) oflag,ofile

    write(*,*) 'Reading: ',trim(ofile)
    call read_file_obs(ofile,atime1,atime2,nsteps,nr,atval,val,dt)
    nrec(k) = nr
    atval_tot(k,:) = atval
    val_tot(k,:) = val

 end do
 close(20)

 open(30,file='antime_list.txt',form='formatted')
 ! loop on all times
 do n = 1,nsteps

    atime = atime1 + (n-1) * dt
    !call dts_from_abs_time(d,t,atime)
    !call date2string(d,t,dstr)
    !write(*,*) 'Time: ',dstr

    isfile = 0
    iflag = 0
    iok = 0
    ! loop on files
    do k = 1,ktot

     ! loop on records in the file
     do nn = 1,nrec(k)

       ! if times are the same
       if (nint(atime-atval_tot(k,nn)) == 0) then
	  if (nint(val_tot(k,nn)-flag) == 0) then
		  iflag = iflag + 1
	  else
		  iok = iok + 1
		  isfile(k) = 1
	  end if
       end if

     end do
    end do

    if (iok > 0) then
	    !write(30,*) iok,iflag,atime
	    call dts_from_abs_time(d,t,atime)
	    call date2string(d,t,dstr)

	    write(30,'(a20,1x,i5)') dstr,ktot
	    write(30,'(9999i2)') isfile(1:ktot)
    end if

 end do

 close(30)

end program make_antime_list


subroutine read_file_obs(ofile,atime1,atime2,nsteps,nrec,atval,val,dt)
  use iso8601
  implicit none
  character(len=80),intent(in) :: ofile
  integer,intent(in) :: nsteps
  double precision,intent(in) :: atime1,atime2,dt
  integer,intent(out) :: nrec
  double precision,intent(inout) :: atval(nsteps),val(nsteps)
  integer ierr
  character(len=20) :: dstring
  integer :: d,t
  double precision :: atime,v,ctime
  integer :: k,j

  open(21,file=trim(ofile),status='old',iostat=ierr)
  if ( ierr /= 0 ) error stop 'Error opening file'

  k = 0
 93 read(21,*,end=92) dstring,v

  call string2date(trim(dstring),d,t,ierr)
  if ( ierr /= 0 ) error stop 'Invalid date'
  call dts_to_abs_time(d,t,atime)

  do j = 0,nsteps-1
     ctime = atime1 + j * dt
     if (atime == ctime) then
        k = k + 1
        atval(k) = atime
        val(k) = v
     end if
  end do

  goto 93

 92 close(21)

  nrec = k

  if ( nrec > nsteps ) error stop 'wrong nrec'

end subroutine read_file_obs
