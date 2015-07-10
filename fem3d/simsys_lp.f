c
c $Id: simsys_lp.f,v 1.3 2009-03-31 15:37:07 georg Exp $
c
c matrix inversion administration (Gaussian elimination - traditional)
c
c revision log :
c
c 12.01.2009	ggu	new file for system routines
c 31.01.2009	ggu	prepared for double precision, new aux vect vs1v,...
c
c******************************************************************
c
c to change from real to double precision
c change amat here and vs*v in common.h
c
c******************************************************************

	subroutine system_initialize

	implicit none

	include 'param.h'
	include 'common.h'

	write(6,*) '----------------------------------------'
	write(6,*) 'initializing matrix inversion routines'
	write(6,*) 'using Gaussian elimination'
	write(6,*) '----------------------------------------'

	end

c******************************************************************

	subroutine system_init

	implicit none

	include 'param.h'
	include 'common.h'
	include 'nbasin.h'

        double precision amat(matdim)
        common /amat/amat
        save /amat/

	!call lp_init_system(nkn,mbw,amat,vs1v)
	call dlp_init_system(nkn,mbw,amat,vs1v)

	end

c******************************************************************

	subroutine system_solve_z

	implicit none

	include 'param.h'
	include 'common.h'
	include 'nbasin.h'

        double precision amat(matdim)
        common /amat/amat
        save /amat/

	!call lp_solve_system(nkn,mbw,amat,vs1v,is2v,vs3v)
	call dlp_solve_system(nkn,mbw,amat,vs1v,is2v,vs3v)

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

        double precision amat(matdim)
        common /amat/amat
        save /amat/

	integer i,j,kk

	integer loclp,loccoo
	external loclp,loccoo

        do i=1,3
          do j=1,3
            kk=loclp(kn(i),kn(j),n,m)
            if(kk.gt.0) amat(kk) = amat(kk) + mass(i,j)
          end do
          vs1v(kn(i)) = vs1v(kn(i)) + rhs(i)
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
          z(k) = vs1v(k)
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
          vs1v(k) = vs1v(k) + dt * array(k)
        end do

        end

c******************************************************************

