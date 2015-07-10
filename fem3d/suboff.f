c
c $Id$
c
c routines for offline data handling
c
c revision log :
c
c 13.06.2013    ggu     new routines written from scratch
c 17.06.2013    ggu     eliminated compiler warnings
c 25.03.2014    ggu     new offline (for T/S)
c 06.05.2015    ccf     write offline to .off file
c 06.05.2015    ccf     read offline from offlin file in section name
c
c****************************************************************

	subroutine offline(mode)

c handles offline version

c-----------------------------------------------------
c
c parameters:
c
c	mode - input parameter
c
c	mode = 1	write to file
c	mode = 2	read from file
c
c	idtoff - parameter set in STR file
c
c	idtoff = 0	nothing (no offline routines called)
c	idtoff > 0	write offline data file with time step idtoff
c	idtoff < 0	reads offline data from file
c	idtoff = -1	uses offline hydro results
c	idtoff = -2	uses offline T/S results
c	idtoff = -4	uses offline turbulence results
c
c combinations are possible: -3,-7
c
c-----------------------------------------------------

	implicit none

	include 'param.h'

	integer mode

	integer nintp			!2 (linear) or 4 (cubic) are possible
	!parameter (nintp=2)		!grade of interpolation
	parameter (nintp=4)		!grade of interpolation

	include 'femtime.h'

	integer ioffline
	common /ioffline/ioffline
	save /ioffline/

	double precision dtr
	integer time(nintp)
	save time
	double precision ut(nlvdim,neldim,nintp)
	double precision vt(nlvdim,neldim,nintp)
	double precision ze(3,neldim,nintp)
	double precision wn(0:nlvdim,nkndim,nintp)
	double precision zn(nkndim,nintp)
	double precision sn(nlvdim,nkndim,nintp)
	double precision tn(nlvdim,nkndim,nintp)
	save dtr
	save ut,vt,ze,wn,zn,sn,tn

	integer iwhat,iread
	integer idtoff,itmoff,itstart
	integer ierr,ig
	real dt
	character*60 name,status
	save idtoff, itmoff, iwhat

	real getpar

        integer ifemop, ifileo
	integer iu,it1,it2,itoff
	save iu,it1,it2,itoff

	integer icall
	save icall
	data icall /0/

	if( icall .lt. 0 ) return

c-------------------------------------------------------------
c initialize
c-------------------------------------------------------------

	if( icall .eq. 0 ) then
	  ioffline = 0

          call convert_date('itmoff',itmoff)
          call convert_time('idtoff',idtoff)

	  if( it .lt. itmoff ) return

  	  if( idtoff .eq. 0 ) iwhat = 0		!nothing
	  if( idtoff .gt. 0 ) iwhat = 1		!write
	  if( idtoff .lt. 0 ) iwhat = 2		!read

	  if( iwhat .le. 0 ) icall = -1
	  if( idtoff .eq. 0 ) icall = -1
	  if( icall .lt. 0 ) return
	  if( iwhat .eq. 1 ) then
            if (iu .le. 0 ) then
              iu = ifemop('.off','unform','new') !unit for writing offline
              if( iu .le. 0 ) then
                write(6,*) 'iu = ',iu
                stop 'error stop offline: cannot open output file'
              end if
	      write(6,*) 'Start writing offline file'
            end if
	  else
            call getfnm('offlin',name)
            iu = ifileo(1,name,'unformatted','old')
            if( iu .le. 0 ) then
              write(6,*) '*** Cannot find offline file: '
              write(6,*) name
              stop
            end if
            write(6,*) '---------------------------------------------'
            write(6,*) '... performing offline from file: '
            write(6,*) name
            write(6,*) '---------------------------------------------'
	  end if
	  call off_init(dtr,ut,vt,ze,wn,zn,sn,tn)
	  itoff = itmoff + idtoff
	  ioffline = -idtoff
	end if

c-------------------------------------------------------------
c do different modes
c-------------------------------------------------------------

	if( mode .ne. iwhat ) return

	if( mode .eq. 1 ) then

c	  -------------------------------------------------------------
c	  accumulate and write data
c	  -------------------------------------------------------------

	  call get_timestep(dt)
	  call off_accum(dt,dtr,ut,vt,ze,wn,zn,sn,tn)

	  if( icall .eq. 0 ) then
	    call off_aver(dtr,ut,vt,ze,wn,zn,sn,tn)
	    call off_write(iu,itmoff,ut,vt,ze,wn,zn,sn,tn)
	    call off_init(dtr,ut,vt,ze,wn,zn,sn,tn)
	    icall = 1
	  end if

	  if( it .lt. itoff ) return

	  call off_aver(dtr,ut,vt,ze,wn,zn,sn,tn)
	  call off_write(iu,it,ut,vt,ze,wn,zn,sn,tn)
	  call off_init(dtr,ut,vt,ze,wn,zn,sn,tn)
	  itoff = itoff + idtoff

	else if( mode .eq. 2 ) then

c	  -------------------------------------------------------------
c	  read data and put into hydro structures
c	  -------------------------------------------------------------

	  if( icall .eq. 0 ) then
	    do ig=1,nintp
	      call off_read(iu,ig,time,ut,vt,ze,wn,zn,sn,tn,ierr,iread)
	      if( ierr .ne. 0 ) goto 97
	    end do
	    call can_do_offline(iread)
	    if( it .lt. time(1) ) goto 99
	    call get_timestep(dt)
	    if( it .eq. itmoff ) then
	      itstart = it
	    else
	      itstart = max(it-nint(dt),itmoff)
	    end if
	    call off_intp_all(iu,nintp,itstart,time,ut,vt,ze,wn,zn,sn,tn)
	    icall = 1
	  end if

	  call off_intp_all(iu,nintp,it,time,ut,vt,ze,wn,zn,sn,tn)

	  !call off_check(1,time,ut,vt,ze,wn,zn,sn,tn)
	  !call off_check(2,time,ut,vt,ze,wn,zn,sn,tn)

	else

c	  -------------------------------------------------------------
c	  error in mode
c	  -------------------------------------------------------------

	  write(6,*) 'mode = ',mode,'  iwhat = ',iwhat
	  stop 'error stop offline: value for mode not allowed'

	end if

c-------------------------------------------------------------
c end of routine
c-------------------------------------------------------------

	return
   97	continue
	write(6,*) time
	write(6,*) nintp,ig
	stop 'error stop offline: read error at start'
   99	continue
	write(6,*) it,time
	stop 'error stop offline: no time available'
	end

c****************************************************************

	subroutine is_offline(type,boff)

c type: 1 hydro, 2 T/S, 4 turb, combinations are possible: 3,7
c type == 0 -> any offline

	implicit none

	integer type	!should we use this offline data?
	logical boff	!data is available and should be used (return)

	integer ioffline
	common /ioffline/ioffline
	save /ioffline/

	integer iwhat

	iwhat = ioffline		!this is what we want (from idtoff)

	if( iwhat .le. 0 ) then		!no offline
	  boff = .false.
	else if( type .eq. 0 ) then	!general
	  boff = .true.
	  !boff = iwhat .gt. 0
	else if( type .eq. 1 ) then	!hydro
	  boff = mod(iwhat/1,2) .ne. 0
	else if( type .eq. 2 ) then	!T/S
	  boff = mod(iwhat/2,2) .ne. 0
	else if( type .eq. 4 ) then	!turbulence
	  boff = mod(iwhat/4,2) .ne. 0
	else
	  write(6,*) 'value for type not allowed: ',type
	  stop 'error stop is_offline: type'
	end if
	  
	end

c****************************************************************

	subroutine can_do_offline(iread)

	implicit none

	integer iread		!this is what we get

	integer ioffline
	common /ioffline/ioffline
	save /ioffline/

	logical bwhat,bread
	integer iwhat,i

	iwhat = ioffline	!this is what we want

	i = 1
	do while( i .le. 4 )
	  bwhat = mod(iwhat/i,2) .ne. 0
	  bread = mod(iread/i,2) .ne. 0
	  if( bwhat .and. .not. bread ) goto 99
	  i = i * 2
	end do

	return
   99	continue
	write(6,*) 'iread = ',iread,'  iwhat = ',iwhat
	write(6,*) 'type = ',i
	write(6,*) 'offline data requested has not been read'
	stop 'error stop can_do_offline: no such data'
	end

c****************************************************************
c****************************************************************
c****************************************************************

	subroutine off_intp_all(iu,nintp,it,time,ut,vt,ze,wn,zn,sn,tn)

	implicit none

	include 'param.h'

	integer iu
	integer nintp
	integer it
	integer time(4)
	double precision ut(nlvdim,neldim,nintp)
	double precision vt(nlvdim,neldim,nintp)
	double precision ze(3,neldim,nintp)
	double precision wn(0:nlvdim,nkndim,nintp)
	double precision zn(nkndim,nintp)
	double precision sn(nlvdim,nkndim,nintp)
	double precision tn(nlvdim,nkndim,nintp)

	include 'nbasin.h'

	include 'levels.h'

	include 'hydro_vel.h'
	include 'hydro.h'
	include 'ts.h'

	logical boff,bhydro,bts
	integer ierr,iread
	integer ip,i,itnext

	integer ieof
	save ieof
	data ieof / 0 /

c	---------------------------------------------------------
c	initialize
c	---------------------------------------------------------

	ip = 2
	if( nintp .eq. 4 ) ip = 3

	call is_offline(1,bhydro)		!hydro
	call is_offline(2,bts)			!T/S

c	---------------------------------------------------------
c	find new records for time
c	---------------------------------------------------------

	do while( ieof .eq. 0 .and. it .gt. time(ip) )
	  call off_next_record(iu,itnext,ieof)
	  if( ieof .ne. 0 ) exit
	  call off_copy(nintp,time,ut,vt,ze,wn,zn,sn,tn)
	  call off_read(iu,nintp,time,ut,vt,ze,wn,zn,sn,tn,ierr,iread)
	end do

	if( it .gt. time(nintp) ) goto 99

	!write(67,*) it,(time(i),i=1,nintp)

c	---------------------------------------------------------
c	pre processing
c	---------------------------------------------------------

	if( bhydro ) then
	  call copy_uvz
	  call copy_depth
	end if

c	---------------------------------------------------------
c	interpolation
c	---------------------------------------------------------

	!if( nintp .eq. 2 ) then
	!  call off_intp2(it,time,ut,vt,ze,wn,zn,sn,tn)
	!else if( nintp .eq. 4 ) then
	!  call off_intp4(it,time,ut,vt,ze,wn,zn,sn,tn)
	!else
	!  write(6,*) 'nintp = ',nintp
	!  stop 'error stop off_intp_all: nintp not possible'
	!end if

	if( bhydro ) then
	  call off_intp(nintp,it,time,nlvdim,neldim,ilhv,nel,ut,utlnv)
	  call off_intp(nintp,it,time,nlvdim,neldim,ilhv,nel,vt,vtlnv)
	  call off_intp(nintp,it,time,1,3*neldim,ilhv,3*nel,ze,zenv)
	  call off_intp(nintp,it,time,nlvdim+1,nkndim,ilhkv,nkn,wn,wlnv)
	  call off_intp(nintp,it,time,1,nkndim,ilhkv,nkn,zn,znv)
	end if

	if( bts ) then
	  call off_intp(nintp,it,time,nlvdim,nkndim,ilhkv,nkn,sn,saltv)
	  call off_intp(nintp,it,time,nlvdim,nkndim,ilhkv,nkn,tn,tempv)
	end if

c	---------------------------------------------------------
c	post processing
c	---------------------------------------------------------

	if( bhydro ) then
	  call make_new_depth
	  call uvint
          call ttov
          call make_prvel
	end if

	if( bts ) then
	  call rhoset_shell
	end if

c	---------------------------------------------------------
c	end of routine
c	---------------------------------------------------------

	return
   99	continue
	write(6,*) 'time to interpolate: it = ',it
	write(6,*) 'time values available in time(): '
	write(6,*) (time(i),i=1,nintp)
	stop 'error stop off_intp_all: no such time'
	end

c****************************************************************

	subroutine off_intp(nintp,it,time,nlvddi,ndim,il,n,dval,rval)

	implicit none

	integer nintp
	integer it
	integer time(nintp)
	integer nlvddi,ndim
	integer il(ndim)
	integer n
	double precision dval(nlvddi,ndim,nintp)
	real rval(nlvddi,ndim)

	integer l,lmax,i,j
	real x(4),y(4),t

	real intp_neville

	if( nintp .lt. 2 .or. nintp .gt. 4 ) then
	  write(6,*) 'nintp = ',nintp
	  stop 'error stop off_intp: nintp not possible'
	end if

	t = it
	do j=1,nintp
	  x(j) = time(j)
	end do

	do i=1,n
	  lmax = 1
	  if( nlvddi .gt. 1 ) lmax = il(i)
	  do l=1,lmax
	    do j=1,nintp
	      y(j) = dval(l,i,j)
	    end do
	    rval(l,i) = intp_neville(nintp,x,y,t)
	  end do
	end do

	end

c****************************************************************

	subroutine off_intp4(it,time,ut,vt,ze,wn,zn)

	implicit none

	include 'param.h'

	integer it
	integer time(4)
	double precision ut(nlvdim,neldim,4)
	double precision vt(nlvdim,neldim,4)
	double precision ze(3,neldim,4)
	double precision wn(0:nlvdim,nkndim,4)
	double precision zn(nkndim,4)

	include 'nbasin.h'
	include 'levels.h'
	include 'hydro_vel.h'
	include 'hydro.h'

	integer ie,ii,k,l,lmax,i,nintp
	real x(4),y(4),t

	real intp_neville

	nintp = 4

	t = it
	do i=1,nintp
	  x(i) = time(i)
	end do
	
	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
	    do i=1,nintp
	      y(i) = ut(l,ie,i)
	    end do
	    utlnv(l,ie) = intp_neville(nintp,x,y,t)
	    do i=1,nintp
	      y(i) = vt(l,ie,i)
	    end do
	    vtlnv(l,ie) = intp_neville(nintp,x,y,t)
	  end do
	  do ii=1,3
	    do i=1,nintp
	      y(i) = ze(ii,ie,i)
	    end do
	    zenv(ii,ie) = intp_neville(nintp,x,y,t)
	  end do
	end do

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax-1
	    do i=1,nintp
	      y(i) = wn(l,k,i)
	    end do
	    wlnv(l,k) = intp_neville(nintp,x,y,t)
	  end do
	  do i=1,nintp
	    y(i) = zn(k,i)
	  end do
	  znv(k) = intp_neville(nintp,x,y,t)
	end do

	end

c****************************************************************

	subroutine off_intp2(it,time,ut,vt,ze,wn,zn)

	implicit none

	include 'param.h'

	integer it
	integer time(2)
	double precision ut(nlvdim,neldim,2)
	double precision vt(nlvdim,neldim,2)
	double precision ze(3,neldim,2)
	double precision wn(0:nlvdim,nkndim,2)
	double precision zn(nkndim,2)

	include 'nbasin.h'

	include 'levels.h'

	include 'hydro_vel.h'
	include 'hydro.h'

	integer ie,ii,k,l,lmax
	integer it1,it2
	double precision rr,rt

	it1 = time(1)
	it2 = time(2)

	rr = 0.
	if( it2 .gt. it1 ) rr = float(it-it1)/float(it2-it1)
	rt = 1. - rr
	
	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
	    utlnv(l,ie) = rt*ut(l,ie,1) + rr*ut(l,ie,2)
	    vtlnv(l,ie) = rt*vt(l,ie,1) + rr*vt(l,ie,2)
	  end do
	  do ii=1,3
	    zenv(ii,ie) = rt*ze(ii,ie,1) + rr*ze(ii,ie,2)
	  end do
	end do

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax-1
	    wlnv(l,k) = rt*wn(l,k,1) + rr*wn(l,k,2)
	  end do
	  znv(k) = rt*zn(k,1) + rr*zn(k,2)
	end do

	end

c****************************************************************
c****************************************************************
c****************************************************************

	subroutine off_copy(nintp,time,ut,vt,ze,wn,zn,sn,tn)

	implicit none

	include 'param.h'

	integer nintp
	integer time(1)
	double precision ut(nlvdim,neldim,1)
	double precision vt(nlvdim,neldim,1)
	double precision ze(3,neldim,1)
	double precision wn(0:nlvdim,nkndim,1)
	double precision zn(nkndim,1)
	double precision sn(nlvdim,nkndim,1)
	double precision tn(nlvdim,nkndim,1)

	include 'nbasin.h'

	include 'levels.h'

	integer ie,ii,k,l,lmax
	integer ito,ifrom

	do ito=1,nintp-1

	  ifrom = ito + 1

	  time(ito) = time(ifrom)

	  do ie=1,nel
	    lmax = ilhv(ie)
	    do l=1,lmax
	      ut(l,ie,ito) = ut(l,ie,ifrom)
	      vt(l,ie,ito) = vt(l,ie,ifrom)
	    end do
	    do ii=1,3
	      ze(ii,ie,ito) = ze(ii,ie,ifrom)
	    end do
	  end do

	  do k=1,nkn
	    lmax = ilhkv(k)
	    do l=0,lmax
	      wn(l,k,ito) = wn(l,k,ifrom)
	      sn(l,k,ito) = sn(l,k,ifrom)
	      tn(l,k,ito) = tn(l,k,ifrom)
	    end do
	    zn(k,ito) = zn(k,ifrom)
	    sn(lmax,k,ito) = sn(lmax,k,ifrom)
	    tn(lmax,k,ito) = tn(lmax,k,ifrom)
	  end do

	end do

	end

c****************************************************************
	
	subroutine off_init(dtr,ut,vt,ze,wn,zn,sn,tn)

	implicit none

	include 'param.h'

	double precision dtr
	double precision ut(nlvdim,neldim)
	double precision vt(nlvdim,neldim)
	double precision ze(3,neldim)
	double precision wn(0:nlvdim,nkndim)
	double precision zn(nkndim)
	double precision sn(nlvdim,nkndim)
	double precision tn(nlvdim,nkndim)

	include 'nbasin.h'

	include 'levels.h'

	integer ie,ii,k,l,lmax

	dtr = 0.

	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
	    ut(l,ie) = 0.
	    vt(l,ie) = 0.
	  end do
	  do ii=1,3
	    ze(ii,ie) = 0.
	  end do
	end do

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=0,lmax
	    wn(l,k) = 0.
	    sn(l,k) = 0.
	    tn(l,k) = 0.
	  end do
	  zn(k) = 0.
	  sn(lmax,k) = 0.
	  tn(lmax,k) = 0.
	end do

	end

c****************************************************************
	
	subroutine off_accum(dt,dtr,ut,vt,ze,wn,zn,sn,tn)

	implicit none

	include 'param.h'

	real dt
	double precision dtr
	double precision ut(nlvdim,neldim)
	double precision vt(nlvdim,neldim)
	double precision ze(3,neldim)
	double precision wn(0:nlvdim,nkndim)
	double precision zn(nkndim)
	double precision sn(nlvdim,nkndim)
	double precision tn(nlvdim,nkndim)

	include 'nbasin.h'

	include 'levels.h'
	include 'hydro_vel.h'
	include 'hydro.h'
	include 'ts.h'

	integer ie,ii,k,l,lmax
	double precision dtt

	dtt = dt
	dtr = dtr + dtt
	dtr = 1.

	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
	    ut(l,ie) = ut(l,ie) + utlnv(l,ie) * dtt
	    vt(l,ie) = vt(l,ie) + vtlnv(l,ie) * dtt
	    ut(l,ie) = utlnv(l,ie)
	    vt(l,ie) = vtlnv(l,ie)
	  end do
	  do ii=1,3
	    ze(ii,ie) = ze(ii,ie) + zenv(ii,ie) * dtt
	    ze(ii,ie) = zenv(ii,ie)
	  end do
	end do

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax-1
	    wn(l,k) = wn(l,k) + wlnv(l,k) * dtt
	    wn(l,k) = wlnv(l,k)
	    sn(l,k) = sn(l,k) + saltv(l,k) * dtt
	    sn(l,k) = saltv(l,k)
	    tn(l,k) = tn(l,k) + tempv(l,k) * dtt
	    tn(l,k) = tempv(l,k)
	  end do
	  zn(k) = zn(k) + znv(k) * dtt
	  zn(k) = znv(k)
	  sn(lmax,k) = sn(lmax,k) + saltv(lmax,k) * dtt
	  sn(lmax,k) = saltv(lmax,k)
	  tn(lmax,k) = tn(lmax,k) + tempv(lmax,k) * dtt
	  tn(lmax,k) = tempv(lmax,k)
	end do

	end

c****************************************************************
	
	subroutine off_aver(dtr,ut,vt,ze,wn,zn,sn,tn)

	implicit none

	include 'param.h'

	double precision dtr
	double precision ut(nlvdim,neldim)
	double precision vt(nlvdim,neldim)
	double precision ze(3,neldim)
	double precision wn(0:nlvdim,nkndim)
	double precision zn(nkndim)
	double precision sn(nlvdim,nkndim)
	double precision tn(nlvdim,nkndim)

	include 'nbasin.h'

	include 'levels.h'

	integer ie,ii,k,l,lmax
	double precision rr

	rr = 0.
	if( dtr .gt. 0. ) rr = 1. / dtr

	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
	    ut(l,ie) = ut(l,ie) * rr
	    vt(l,ie) = vt(l,ie) * rr
	  end do
	  do ii=1,3
	    ze(ii,ie) = ze(ii,ie) * rr
	  end do
	end do

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax-1
	    wn(l,k) = wn(l,k) * rr
	    sn(l,k) = sn(l,k) * rr
	    tn(l,k) = tn(l,k) * rr
	  end do
	  zn(k) = zn(k) * rr
	  sn(lmax,k) = sn(lmax,k) * rr
	  tn(lmax,k) = tn(lmax,k) * rr
	end do

	end

c****************************************************************
c****************************************************************
c****************************************************************

	subroutine off_check(ig,time,ut,vt,ze,wn,zn,sn,tn)

	implicit none

	include 'param.h'

	integer ig
	integer time(1)
	double precision ut(nlvdim,neldim,1)
	double precision vt(nlvdim,neldim,1)
	double precision ze(3,neldim,1)
	double precision wn(0:nlvdim,nkndim,1)
	double precision zn(nkndim,1)
	double precision sn(nlvdim,nkndim,1)
	double precision tn(nlvdim,nkndim,1)

	include 'nbasin.h'

	include 'levels.h'


	include 'hydro_vel.h'

	include 'hydro_print.h'


	include 'hydro.h'

	integer ie,ii,k,l,lmax
	integer ierr
	real utmax,umax,zmax,wmax,smax,tmax

	ierr = 0
	utmax = 10000.
	umax = 10.
	zmax = 10.
	wmax = 10.
	smax = 100.
	tmax = 100.

	do ie=1,nel
	  lmax = ilhv(ie)
	  do l=1,lmax
	    call off_check_val('ut',ie,l,real(ut(l,ie,ig)),utmax,ierr)
	    call off_check_val('vt',ie,l,real(vt(l,ie,ig)),utmax,ierr)
	    call off_check_val('utlnv',ie,l,utlnv(l,ie),utmax,ierr)
	    call off_check_val('vtlnv',ie,l,vtlnv(l,ie),utmax,ierr)
	    call off_check_val('utlov',ie,l,utlov(l,ie),utmax,ierr)
	    call off_check_val('vtlov',ie,l,vtlov(l,ie),utmax,ierr)
	    call off_check_val('ulnv',ie,l,ulnv(l,ie),umax,ierr)
	    call off_check_val('vlnv',ie,l,vlnv(l,ie),umax,ierr)
	    call off_check_val('ulov',ie,l,ulov(l,ie),umax,ierr)
	    call off_check_val('vlov',ie,l,vlov(l,ie),umax,ierr)
	  end do
	  do ii=1,3
	    call off_check_val('ze',ie,ii,real(ze(ii,ie,ig)),zmax,ierr)
	  end do
	end do

	do k=1,nkn
	  lmax = ilhkv(k)
	  do l=1,lmax-1
	    call off_check_val('wn',k,l,real(wn(l,k,ig)),wmax,ierr)
	  end do
	  do l=1,lmax
	    call off_check_val('sn',k,l,real(sn(l,k,ig)),smax,ierr)
	    call off_check_val('tn',k,l,real(tn(l,k,ig)),tmax,ierr)
	  end do
	  call off_check_val('zn',k,0,real(zn(k,ig)),zmax,ierr)
	end do

	if( ierr .gt. 0 ) then
	  write(6,*) 'errors checking variables read from file'
	  write(6,*) time(ig),ig,ierr
	  stop 'error stop off_check: out of range'
	else
	  !write(6,*) 'finished offline error check... ok... ',it
	end if

	end

c****************************************************************

	subroutine off_check_val(what,iek,l,val,vmax,ierr)

	implicit none

	character*(*) what
	integer iek,l
	real val,vmax
	integer ierr

	if( abs(val) .gt. vmax ) then
	  write(6,*) what,iek,l,val
	  ierr = ierr + 1
	end if

	end

c****************************************************************
c****************************************************************
c****************************************************************
	
	subroutine off_write(iu,it,ut,vt,ze,wn,zn,sn,tn)

	implicit none

	include 'param.h'

	integer iu,it
	double precision ut(nlvdim,neldim)
	double precision vt(nlvdim,neldim)
	double precision ze(3,neldim)
	double precision wn(0:nlvdim,nkndim)
	double precision zn(nkndim)
	double precision sn(nlvdim,nkndim)
	double precision tn(nlvdim,nkndim)

	include 'nbasin.h'

	include 'levels.h'

	integer ie,ii,k,l,lmax

	write(iu) it,nkn,nel,3
	write(iu) (ilhv(ie),ie=1,nel)
	write(iu) (ilhkv(k),k=1,nkn)
	write(iu) ((ut(l,ie),l=1,ilhv(ie)),ie=1,nel)
	write(iu) ((vt(l,ie),l=1,ilhv(ie)),ie=1,nel)
	write(iu) ((ze(ii,ie),ii=1,3),ie=1,nel)
	write(iu) ((wn(l,k),l=1,ilhkv(k)-1),k=1,nkn)
	write(iu) (zn(k),k=1,nkn)
	write(iu) ((sn(l,k),l=1,ilhkv(k)),k=1,nkn)
	write(iu) ((tn(l,k),l=1,ilhkv(k)),k=1,nkn)

	end

c****************************************************************

	subroutine off_read(iu,ig,time,ut,vt,ze,wn,zn,sn,tn,ierr,iread)

	implicit none

	include 'param.h'

	integer iu,ig
	integer time(1)
	double precision ut(nlvdim,neldim,1)
	double precision vt(nlvdim,neldim,1)
	double precision ze(3,neldim,1)
	double precision wn(0:nlvdim,nkndim,1)
	double precision zn(nkndim,1)
	double precision sn(nlvdim,nkndim,1)
	double precision tn(nlvdim,nkndim,1)
	integer ierr
	integer iread

	include 'nbasin.h'

	include 'levels.h'

	integer ilhaux(neldim)
	integer ilhkaux(nkndim)

	integer ie,ii,k,l,lmax,it
	integer nknaux,nelaux
	integer type

	read(iu,err=99,end=98) it,nknaux,nelaux,iread
	if( nkn .ne. nknaux .or. nel .ne. nelaux ) goto 97
	if( iread .ne. 3 ) goto 96
	!write(6,*) 'offline record read: ',it,ig
	time(ig) = it
	read(iu) (ilhaux(ie),ie=1,nel)
	read(iu) (ilhkaux(k),k=1,nkn)
	call off_check_vertical(nel,ilhaux,ilhv)
	call off_check_vertical(nkn,ilhkaux,ilhkv)
	read(iu) ((ut(l,ie,ig),l=1,ilhv(ie)),ie=1,nel)
	read(iu) ((vt(l,ie,ig),l=1,ilhv(ie)),ie=1,nel)
	read(iu) ((ze(ii,ie,ig),ii=1,3),ie=1,nel)
	read(iu) ((wn(l,k,ig),l=1,ilhkv(k)-1),k=1,nkn)
	read(iu) (zn(k,ig),k=1,nkn)
	read(iu) ((sn(l,k,ig),l=1,ilhkv(k)),k=1,nkn)
	read(iu) ((tn(l,k,ig),l=1,ilhkv(k)),k=1,nkn)

	ierr = 0

	return
   96	continue
	write(6,*) 'type: ',type
	stop 'error stop off_read: we must have type == 3'
   97	continue
	write(6,*) 'nkn,nknaux: ',nkn,nknaux
	write(6,*) 'nel,nelaux: ',nel,nelaux
	stop 'error stop off_read: parameter mismatch'
   98	continue
	write(6,*) 'EOF encountered: ',iu,ig
	ierr = -1
	return
	!stop 'error stop off_read: EOF encountered'
   99	continue
	write(6,*) iu,ig
	stop 'error stop off_read: error reading record'
	end

c****************************************************************

	subroutine off_check_vertical(n,ilaux,il)

	implicit none

	integer n
	integer ilaux(n)
	integer il(n)

	integer i

	do i=1,n
	  if( il(i) .le. 0 ) il(i) = ilaux(i)
	  if( il(i) .ne. ilaux(i) ) then
	    write(6,*) i,il(i),ilaux(i)
	    stop 'error stop off_check_vertical: not compatible'
	  end if
	end do

	end 

c****************************************************************
c****************************************************************
c****************************************************************

	subroutine off_next_record(iu,it,ierr)

	implicit none

	integer iu,it,ierr

	integer nknaux,nelaux

	read(iu,err=99,end=98) it,nknaux,nelaux
	backspace(iu)
	ierr = 0

	return
   98	continue
	it = 0
	ierr = -1
	return
   99	continue
	write(6,*) iu
	stop 'error stop off_next_record: error reading record'
	end

c****************************************************************
c****************************************************************
c****************************************************************
c
c	subroutine get_timestep(dt)
c	end
c	subroutine make_new_depth
c	end
c	subroutine uvint
c	end
c	subroutine ttov
c	end
c	subroutine make_prvel
c	end
c	subroutine copy_uvz
c	end
c	subroutine off_test
c	call offline(1,iwhat)
c	end
c	program off_main
c	call off_test
c	end
c
c****************************************************************

