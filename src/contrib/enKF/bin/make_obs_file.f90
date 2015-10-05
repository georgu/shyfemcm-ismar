!---------------------------------
 program make_obs_files
!---------------------------------
!---------------------------------
!---------------------------------
!! TODO a generic progam working with different observations!!
!---------------------------------
!---------------------------------
! Reads a text file with the names of the obs files and
! makes one obs file with all the measurements, each row
! is:  
! time type_of_observation x_coord y_coord z_coord obs_value obs_error
 implicit none

 character (len = 60), dimension (1:100) :: filin
 character (len = 10), dimension (1:100) :: otype

 integer n, nend
 integer k, kend

 ! Timeseries
 real x_ts, y_ts, z_ts, std_ts
 double precision, allocatable :: t_ts(:)
 real, allocatable :: val_ts(:)
 double precision tt

! asspar parameters
 character assdir*40
 integer nrens
 integer date, time
 double precision tt_anf,tt_eanf,tt_oanf,tt_oend,tt_eend,tt_end

 ! Read the list of observation files
 call read_obs_list(otype,filin,nend)

 ! Reads the parameters of the simulation
 call read_asspar(23,assdir,nrens,date,time,tt_anf,tt_end,tt_eanf,tt_eend,tt_oanf,tt_oend)

 open(30,file=trim(assdir)//'/observations.bin',form='unformatted')
 open(31,file=trim(assdir)//'/obs_times.dat',form='formatted')

 ! Loop on obs files
 do n = 1,nend

    if (trim(otype(n)).eq.'times_lev') then

       ! Reads a time-series of level
       call open_times_lev(25,filin(n),kend,x_ts,y_ts,z_ts,std_ts)
       allocate(t_ts(kend),val_ts(kend))
       call read_times_lev(25,kend,t_ts,val_ts)

       ! write in a binary file
       do k =  1,kend
          tt = t_ts(k)
          if(( tt.ge.tt_oanf ).and.( tt.le.tt_oend )) then
            write(30) tt,otype(n),x_ts,y_ts,z_ts,val_ts(k),std_ts
            write(31,*) tt
          end if
       end do

    else if (trim(otype(n)).eq.'times_sal') then
       print*, 'sal'
    else
       write(*,*) 'Warning: obs type not found  ', trim(otype(n))
    end if

 end do
 close(30)
 close(31)

 end
