
!*********************************************************************

	program parse_bfm

! handle description for bfm variables

	implicit none

	integer, parameter :: ndim = 100

	integer ios
	integer in_description
	integer ivar,iv
	integer ivars(ndim)
	character*256 string
	character*256 long
	character*20 short
	character*256 longs(ndim)
	character*20 shorts(ndim)

	in_description = 0

	do
	  read(5,'(a)',iostat=ios) string
	  if( ios /= 0 ) exit
	  if( string(1:6) == 'bfm={' ) then
	    in_description = 1
	    write(6,*) 'start of bfm section found'
	    call clean_first_entry(string)
	  end if
	  if( string(1:2) == '}' .and. in_description == 1 ) exit
	  if( in_description == 1 ) then
	    call scan_bfm(string,ivar,short,long)
	    iv = iv + 1
	    if( iv > ndim ) stop 'error stop parse_bfm: ndim'
	    call adjust_string(long)
	    ivars(iv) = ivar
	    shorts(iv) = short
	    longs(iv) = long
	    call check_entries(iv,ivars,shorts,longs)
	  end if
	end do

	write(6,*) 'variables read: ',iv

	call write_code(iv,ivars,shorts,longs)

	end

!*********************************************************************

	subroutine clean_first_entry(string)

! cleans first entry from intro string

	implicit none

	character*(*) string

	integer l,i

	l = len(string)

	do i=1,l
	  if( string(i:i) == '"' ) exit
	  string(i:i) = ' '
	end do

	end

!*********************************************************************

	subroutine scan_bfm(string,ivar,short,long)

! scans string and returns short and long description

	implicit none

	character*(*) string
	integer ivar
	character*(*) short,long

	integer l,i,ioff,j
	real f
	character*20 strings(3),s
	character*20 svar,saux
	character*1 c

	integer iscans,istos,iscanf

	! get first three parts of string

	strings = ' '
	i = iscans(string,strings,3)
	if( i /= 3 ) then
	  write(6,*) trim(string)
	  stop 'error stop scan_bfm: cannot parse'
	end if
	s = strings(1)
	saux = strings(3)

!	parse first part and get svar and short

	l = len(s)
	do i=1,l
	  c = s(i:i)
	  if( c == '_' ) then
	    !write(6,*) '_ found: ',i
	    svar = s(i+1:i+3)
	    short = s(i+7:i+9)
	    if( short /= saux ) then
	      write(6,*) short,' ',saux
	      stop 'error stop scan_bfm: incompatibility'
	    end if
	    j = iscanf(svar,f,1)
	    if( j /= 1 ) stop 'error stop scan_bfm: converting number'
	    ivar = nint(f)
	  end if
	end do

!	skip first 3 strings to get to description

	l = len(string)
	do i=1,l
	  if( string(i:i+4) == '# '//short(1:3) ) exit
	end do
	if( i > l ) then
	  stop 'error stop scan_bfm: string not found'
	end if
	long = string(i+5:)

	write(6,*) ivar,'  ',trim(svar),'  ',trim(short),'  ',trim(long)

	end

!*********************************************************************

	subroutine check_entries(iv,ivars,shorts,longs)

! checks entries for double entries

	implicit none

	integer iv
	integer ivars(iv)
	character*(*) shorts(iv)
	character*(*) longs(iv)

	integer i
	logical bdouble

	do i=1,iv-1
	  bdouble = .false.
	  if( ivars(i) == ivars(iv) ) bdouble = .true.
	  if( shorts(i) == shorts(iv) ) bdouble = .true.
	  if( longs(i) == longs(iv) ) bdouble = .true.
	  if( bdouble ) then
	    write(6,*) 'double entry: '
	    write(6,*) ivars(i),ivars(iv)
	    write(6,*) shorts(i),shorts(iv)
	    write(6,*) longs(i),longs(iv)
	    stop 'error stop check_entries: double entry'
	  end if
	end do

	end

!*********************************************************************

	subroutine squeeze(string)

! condenses multiple white space into one

	implicit none

	character*(*) string

	integer l,i,is

	l = len(string)
	is = 1

	do i=2,l
	  if( string(i:i) == ' ' ) then
	    if( string(is:is) == ' ' ) cycle
	  end if
	  is = is + 1
	  string(is:is) = string(i:i)
	end do

	string(is+1:) = ' '

	end

!*********************************************************************

	subroutine adjust_string(long)

! adjusts long description

	implicit none

	character*(*) long

	integer j

	call squeeze(long)		!squeeze white space

	! separate units with [ from description

	j = index(long,' mmol')
	if( j > 0 ) long = long(1:j) // '[' // long(j+1:)
	j = index(long,' mg')
	if( j > 0 ) long = long(1:j) // '[' // long(j+1:)

	long = trim(long) // ']'	!end unit with ]

	end

!*********************************************************************

	subroutine write_code(iv,ivars,shorts,longs)

! writes out code to be included into shyfem

	implicit none

	integer iv
	integer ivars(iv)
	character*(*) shorts(iv)
	character*(*) longs(iv)

	integer i,lmax
	character*80 short,long
	character*80 line
	character*4, parameter :: prefix = 'bfm_'

	open(1,file='bfm_strings.f90',status='unknown',form='formatted')

	write(1,*)
	write(1,*) '	subroutine populate_bfm_strings'
	write(1,*)
	write(1,*) '	! set up strings for BFM'
	write(1,*)
	write(1,*) '	use shyfem_strings'
	write(1,*)
	write(1,*) '	implicit none'
	write(1,*)

	lmax = 0
	do i=1,iv
	  long = adjustl(longs(i))
	  lmax = max(lmax,len_trim(long))
	  
	  line = '''' // trim(long) // ''','
	  write(1,'(a,i4,a)') '	call strings_add_new(' // trim(line)  &
     &			, ivars(i) , ')'
	end do

	write(1,*)

	do i=1,iv
	  short = adjustl(shorts(i))
	  if( prefix /= ' ' ) short = prefix // short
	  
	  line = ',''' // trim(short) // ''''
	  write(1,'(a,i4,a)') '	call strings_set_short('  &
     &			, ivars(i) , trim(line) // ')'
	end do

	write(1,*)
	write(1,*) '	end subroutine'
	write(1,*)

	close (1)

	write(6,*) 'max length of description: ',lmax

	end

!*********************************************************************

