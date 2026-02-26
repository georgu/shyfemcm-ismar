
	program test_random

	use mod_random

	implicit none

	integer i
	real r

	!call init_random_seed
	call set_random_seed(6789)

	do i=1,10
	  r = grand()
	  write(6,'(i5,f10.7)') i,r
	  write(66,'(i5,f10.7)') i,r
	end do

	end

