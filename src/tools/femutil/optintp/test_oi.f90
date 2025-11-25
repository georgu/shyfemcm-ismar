
	program test_oi

	implicit none

	integer i
	real pi,rad
	real t,dt,tt,tmax
	real zo,zb,za,zn
	real alpha,w,so,sb,tau,r,rr
	character*70 header
	character*70 text,text1

	pi = 4.*atan(1.)
	rad = 180./pi
	tt = 43200.
	tmax = 3*tt
	dt = 300.
	dt = 3600.
	dt = 300
	tau = 0.
	tau = 10800.
	tau = 300.
	tau = 0
	sb = 0.3
	sb = 0.1
	sb = 2.0
	sb = 2.0
	so = 0.1
	alpha = so**2 / sb**2
	w = 1/(1.+alpha)
	zb = 2.

	write(text,2000) 'tau=',nint(tau),' sb=',nint(10.*sb),' dt=',nint(dt)
 2000	format(a,i5,a,i3,a,i4)
	text = adjustl(text)
	text1 = text
	do i=1,len(trim(text1))
	  if( text1(i:i) == ' ' ) text1(i:i) = '_'
	end do
	write(6,'(a)') text
	write(77,'(a)') trim(text)
	write(6,'(a)') text1
	write(88,'(a)') trim(text1)

	rr = 0.
	if( tau > 0 ) then
	  r = dt/tau
	  rr = 1./(1.+r)
	end if

	header = '#       time         obs  background' &
     &				// '    analysis         new'
	write(6,'(a)') header
	write(66,'(a)') header

	t = 0.
	do while( t < tmax )

	  t = t + dt
	  zo = 3.+sin((t/tt)*pi)
	  za = zb + w*(zo-zb)
	  zn = rr*zb + (1.-rr)*za
	  write(6,1000) t,zo,zb,za,zn
	  write(66,1000) t,zo,zb,za,zn
	  zb = zn

	end do

 1000	format(f12.2,4f12.4)
	end
