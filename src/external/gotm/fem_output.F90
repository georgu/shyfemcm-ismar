
!=================================================================
        module fem_gotm_output
!=================================================================

        implicit none

        logical, save :: boutput = .true.

!=================================================================
        end module fem_gotm_output
!=================================================================

	subroutine set_gotm_output(bw)

	use fem_gotm_output

	implicit none

	logical bw

	boutput = bw

	end subroutine set_gotm_output

!*****************************************************************

