
!*********************************************************************
!
! please run as: a.out [-check|-nocheck] what < shyfem_vars.py
!
! where what is bfm, afm, or similar
!
!*********************************************************************

	program parse_bfm

! handle description for bfm variables

	use clo

	implicit none

	integer, parameter :: ndim = 500

	logical bcheck
	integer ios
	integer in_description
	integer ivar,iv
	integer nolong,ls
	integer ivars(ndim)
	character*256 string
	character*256 long
	character*20 short
	character*256 longs(ndim)
	character*20 shorts(ndim)
	character*20 what

	logical is_comment

	bcheck = .true.
	bcheck = .false.
	what = 'bfm'
	what = 'afm'

	nolong = 0
	in_description = 0
	iv = 0

	call init_parse_bfm(what,bcheck)
	write(6,*) 'processing strings for ',trim(what)
	ls = len_trim(what) + 2

	do
	  read(5,'(a)',iostat=ios) string
	  if( ios /= 0 ) exit
	  if( string(1:ls) == trim(what)//'={' ) then
	    in_description = 1
	    write(6,*) 'start of bfm section found'
	    call clean_first_entry(string)
	  end if
	  if( string(1:1) == '}' .and. in_description == 1 ) exit
	  if( in_description == 1 ) then
	    !if( string == ' ' ) cycle
	    if( is_comment(string) ) cycle
	    call scan_bfm(string,ivar,short,long)
	    if( short == ' ' ) cycle
	    iv = iv + 1
	    if( iv > ndim ) stop 'error stop parse_bfm: ndim'
	    call adjust_string(long)
	    ivars(iv) = ivar
	    shorts(iv) = short
	    longs(iv) = long
	    if( long == ' ' ) nolong = nolong + 1
	    if( bcheck ) call check_entries(iv,ivars,shorts,longs)
	  end if
	end do

	if( in_description == 0 ) then
	  write(6,*) '*** could not find descriptor ',trim(what)
	  stop 'error stop parse_bfm: no descriptor'
	end if

	write(6,*) 'variables read: ',iv
	if( nolong > 0 ) then
	  write(6,*) 'variables with no description: ',nolong
	end if

	call write_code(what,iv,ivars,shorts,longs)

	call exit(99)

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
	  if( string(i:i) == '#' ) exit
	  string(i:i) = ' '
	end do

	end

!*********************************************************************

	function is_comment(string)

! checks if line is a comment

	implicit none

	logical is_comment
	character*(*) string

	integer l,i

	l = len(string)
	is_comment = .true.

	do i=1,l
	  if( string(i:i) /= ' ' ) exit
	end do

	if( i > l ) return			! empty - treat as comment
	if( string(i:i) == '#' ) return
	is_comment = .false.

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
	character*80 strings(3),s
	character*80 strings2(3),ss
	character*80 strings3(3),sss
	character*20 svar,saux
	character*1 c

	integer iscans,istos,iscanf

	! get first three parts of string

	strings = ' '
	i = iscans(string,strings,3)
	if( i < 1 ) then
	  write(6,*) trim(string)
	  stop 'error stop scan_bfm: cannot parse'
	end if
	s = strings(1)
	saux = strings(3)	! this is start of description
	long = ' '

!	parse first part and get svar and short

	call delete_chars(s,'":,''')

	strings2 = ' '
	i = iscans(s,strings2,3)
	if( i /= 2 ) then
	  write(6,*) trim(s)
	  stop 'error stop scan_bfm: cannot parse first part of string'
	end if

	ss = strings2(1)
	if( ss(1:4) /= "var_" ) then
	  write(6,*) trim(ss)
	  stop 'error stop scan_bfm: cannot find variable number'
	end if
	svar = ss(5:)
	j = iscanf(svar,f,1)
	if( j /= 1 ) stop 'error stop scan_bfm: converting number'
	ivar = nint(f)
	short = strings2(2)
	if( short == 'NULL' ) short = ' '
	if( short == ' ' ) return		!end of list

!	find # and insert rest into long

	i = index(string,'#')
	if( i == 0 ) then		!no description
	  !write(6,*) 'cannot find #'
	  !write(6,*) trim(string)
	  !stop 'error stop scan_bfm: no #'
	else
	  long = string(i+1:)
	end if
	call skip_leading_white_space(long)
	if( long == ' ' ) long = short

	if( short == saux ) then !first part of description is equal to short
	  call skip_next_word(long)
	end if

	write(6,*) ivar,'  ',trim(svar),'  ',trim(short),'  ',trim(long)

	return
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
	    write(6,*) i,ivars(i),' ',trim(shorts(i)),' ',trim(longs(i))
	    write(6,*) iv,ivars(iv),' ',trim(shorts(iv)),' ',trim(longs(iv))
	    !stop 'error stop check_entries: double entry'
	    call exit(77)
	  end if
	end do

	end

!*********************************************************************

	subroutine skip_next_word(string)

! skips next word in string

	implicit none

	character*(*) string

	integer i,ioff

	ioff = 1
	call skipwh(string,ioff)
	if( ioff > 1 ) string = string(ioff:)

	i = index(string,' ')
	if( i == 0 ) then
	  string = ' '
	else
	  string = string(i+1:)
	end if

	ioff = 1
	call skipwh(string,ioff)
	if( ioff > 1 ) string = string(ioff:)

	end

!*********************************************************************

	subroutine delete_chars(string,chars)

! deletes single chars from string

	implicit none

	character*(*) string,chars

	integer l,lc
	integer i,j

	l = len(string)
	lc = len(chars)

	do i=1,l
	  do j=1,lc
	    if( string(i:i) == chars(j:j) ) string(i:i) = ' '
	  end do
	end do

	end

!*********************************************************************

	subroutine skip_leading_white_space(string)

	implicit none

	character*(*) string

	integer ioff

	ioff = 1
	call skipwh(string,ioff)
	if( ioff > 1 ) string = string(ioff:)

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

! adjusts long description - separate units with [] from description

	implicit none

	character*(*) long

	integer j

	call squeeze(long)		!squeeze white space

	j = index(long,' mmol')
	if( j > 0 ) long = long(1:j) // '[' // trim(long(j+1:)) // ']'
	j = index(long,' mg')
	if( j > 0 ) long = long(1:j) // '[' // trim(long(j+1:)) // ']'
	j = index(long,' adimensional')
	if( j > 0 ) long = long(1:j) // '[' // trim(long(j+1:)) // ']'
	j = index(long,' deg Celsius')
	if( j > 0 ) long = long(1:j) // '[' // trim(long(j+1:)) // ']'

	end

!*********************************************************************

	subroutine write_code(what,iv,ivars,shorts,longs)

! writes out code to be included into shyfem

	implicit none

	character*(*) what
	integer iv
	integer ivars(iv)
	character*(*) shorts(iv)
	character*(*) longs(iv)

	integer i,lmax
	character*80 short,long
	character*80 line,file
	character*20 :: prefix

	prefix=trim(what) // '_'
	file='new_strings.f90'
	file=trim(prefix) // 'strings.f90'
	write(6,*) 'writing file ' // trim(file)
	open(1,file=file,status='unknown',form='formatted')

	write(1,*)
	line = '	subroutine populate_strings_' // trim(what)
	write(1,*) trim(line)
	write(1,*)
	line = '	! set up strings for ' // trim(what)
	write(1,*) trim(line)
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
	  if( prefix /= ' ' ) short = trim(prefix) // short
	  
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

	subroutine init_parse_bfm(what,bcheck)

	use clo

	implicit none

	character*(*) what
	logical bcheck

	call clo_init('parse_bfm','what','1.0')
	call clo_add_info('parses strings for BFM routines')
        call clo_add_option('check',.false.,'check if entries are unique' &
     &				,'do not check if entries are unique')
        call clo_parse_options
        call clo_get_option('check',bcheck)
	!write(6,*) 'bcheck = ',bcheck

        call clo_get_file(1,what)
        if( what == ' ' ) call clo_usage

	end

!*********************************************************************


