!**************************************************************
! Subroutines to read and write files
!**************************************************************

!----------------------------------------
 subroutine read_asspar(iunit,assdir,nrens,date,time,tt_anf,tt_end,tt_eanf,tt_eend,tt_oanf,tt_oend)
 implicit none
 integer, intent(in) :: iunit

 character, intent(out) :: assdir*40
 integer, intent(out) :: nrens
 integer date, time
 double precision, intent(out) :: tt_anf,tt_eanf,tt_oanf,tt_oend,tt_eend,tt_end
 
 open(iunit,file='asspar.txt',status='old')

 read(iunit,'(a40)') assdir	! working dir
 read(iunit,*) nrens		! nr ens member
 read(iunit,*) date		! date at time 0: yyyymmdd
 read(iunit,*) time		! time
 read(iunit,*) tt_anf		! initial time
 read(iunit,*) tt_end		! final time
 read(iunit,*) tt_eanf		! time first ens creation
 read(iunit,*) tt_eend		! last time with ensemble
 read(iunit,*) tt_oanf		! time first observation
 read(iunit,*) tt_oend		! time last observation
 
 close(iunit)

 end subroutine


!----------------------------------------
 subroutine read_obs_list(otype,filin,nfile)
 implicit none
 character (len = *), dimension (1:100), intent(out) :: filin
 character (len = *), dimension (1:100), intent(out) :: otype
 integer n,nfile
 
 n = 1
 open(20,file='obs_list.txt',status='old')
900 read(20,'(a10,1x,a60)',end = 901) otype(n), filin(n)
    n = n + 1
    goto 900
901 close(20)
 nfile = n - 1

 end subroutine

!----------------------------------------
 subroutine open_times_lev(iunit,tsfile,kend,x,y,z,std)
 implicit none
 integer, intent(in) :: iunit
 character (len = *), intent(in) :: tsfile
 integer, intent(out) :: kend
 real, intent(out) :: x,y,z,std

 open(iunit,file=tsfile,status='old')
 read(iunit,*) kend, x, y, z, std

 end subroutine

!----------------------------------------
 subroutine read_times_lev(iunit,kend,timets,valts)
 implicit none
 integer, intent(in) :: iunit
 integer, intent(in) :: kend
 double precision timets(*), valts(*)
 integer k, kend_real

 k = 1
900 read(iunit,*,end = 901) timets(k), valts(k)
    k = k + 1
    goto 900
901 close(iunit)
 kend_real = k - 1

 if (kend_real.ne.kend) stop 'read_times_lev: bad file length'

 end subroutine

!----------------------------------------
 subroutine read_obs(tt,iunit,filin,nobs,otype,ox,oy,oz,oval,ostd)
 implicit none
 double precision, intent(in) :: tt
 integer, intent(in) :: iunit
 character (len = *), intent(in) :: filin
 real, intent(out) ::  ox(*),oy(*),oz(*),oval(*),ostd(*)
 character (len=10), intent(out) :: otype(*)
 integer, intent(out) :: nobs

 double precision ttv
 real x,y,z,val,stdv
 character (len=10) :: typev

 integer k

 k = 0
 open(iunit,file=filin,status='old',form='unformatted')
 55 read(iunit,end=56) ttv,typev,x,y,z,val,stdv
    if( ttv.ne.tt ) goto 55
      k = k + 1
      write(*,*) 'Observation found at t = ',ttv,'  type: ',typev
      otype(k) = typev
      ox(k) = x
      oy(k) = y
      oz(k) = z
      oval(k) = val
      ostd(k) = stdv
    goto 55
 56 close(iunit)

 if( k.eq.0 ) write(*,*) 'Warning! No observation found'

 ! Number of observation at timestep tt
 nobs = k

 end subroutine


