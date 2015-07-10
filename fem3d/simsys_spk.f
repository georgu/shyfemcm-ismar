c
c revision log :
c
c 12.01.2009	ggu	new file for system routines
c 31.03.2009	ggu	call renamed to spk_*
c 25.05.2015	ggu	some calls changed (pass array in)
c
c******************************************************************

        subroutine system_initialize

        implicit none

        include 'param.h'
        include 'common.h'

        write(6,*) '----------------------------------------'
        write(6,*) 'initializing matrix inversion routines'
        write(6,*) 'using Sparskit routines'
        write(6,*) '----------------------------------------'

        end

c******************************************************************

	subroutine system_init

	implicit none

        include 'param.h'
	include 'common.h'

	call spk_init_system

	end

c******************************************************************

	subroutine system_solve_z(n,z)

	implicit none

	integer n
	real z(n)

        include 'param.h'
	include 'common.h'

	call spk_solve_system(n,z)

	end

c******************************************************************

	subroutine system_assemble(n,m,kn,mass,rhs)

	implicit none

	integer n,m
	integer kn(3)
	real mass(3,3)
	real rhs(3)

        include 'param.h'
	include 'common.h'

	integer i,j,kk

	integer loclp,loccoo
	external loclp,loccoo

        do i=1,3
          do j=1,3
            kk=loccoo(kn(i),kn(j),n,m)
            if(kk.gt.0) coo(kk) = coo(kk) + mass(i,j)
          end do
          rvec(kn(i)) = rvec(kn(i)) + rhs(i)
        end do

	end

c******************************************************************

        subroutine system_adjust_z(n,z)

        implicit none

	integer n
	real z(n)

        include 'param.h'
	include 'common.h'

        integer k

        do k=1,n
          z(k) = rvec(k)
        end do

        end

c******************************************************************

        subroutine system_add_rhs(dt,n,array)

        implicit none

        real dt
	integer n
        real array(n)

        include 'param.h'
	include 'common.h'

        integer k

        do k=1,n
          rvec(k) = rvec(k) + dt * array(k)
        end do

        end

c******************************************************************

