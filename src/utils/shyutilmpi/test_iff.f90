
	program test_iff

	use intp_fem_file

	implicit none

	integer nc
	integer nintp,nvar,id
	integer date,time
	integer ierr
	double precision dtime,atime,atime0
	double precision dstart,dend,dt
	real values(10)
	real, parameter :: flag = -999.
	character*80 file
	character*80 aend,astart
	character*20 aline

	date = 20130101
	time = 0

	call dts_to_abs_time(date,time,atime0)

        nc = command_argument_count()
	if( nc <= 0 ) stop 'no file given'
        call get_command_argument(nc,file)

	call iff_init_global_simplified(date)

	nintp = 2
	dtime = -1
	nvar = 5
	call iff_ts_init(dtime,file,nintp,nvar,id)

	write(6,*) 'file has been initialized'

	astart = '2013-01-01::00:00:00'
	call dts_string2time(astart,atime,ierr)
	if( ierr /= 0 ) stop 'error stop converting date'
	dstart = atime - atime0
	aend = '2014-01-01::00:00:00'
	call dts_string2time(aend,atime,ierr)
	if( ierr /= 0 ) stop 'error stop converting date'
	dend = atime - atime0
	dt = 86400.
	dt = 3600.
	dt = 50000.

	write(6,*) 'starting time loop'

	dtime = dstart
	do
	  values(1:nvar) = flag
	  if( iff_file_has_data(id,dtime) ) then
	    call iff_ts_intp(id,dtime,values)
	    atime = dtime + atime0
	    call dts_format_abs_time(atime,aline)
	    write(6,'(a20,f10.0,5f9.3)') aline,dtime,values(1:nvar)
	  end if
	  dtime = dtime + dt
	  if( dtime > dend ) exit
	end do

	call iff_close_file(id)

	end


