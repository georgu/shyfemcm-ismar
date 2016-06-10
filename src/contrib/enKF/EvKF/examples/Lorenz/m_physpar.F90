module m_physpar
contains
subroutine physpar(tfin,var,dt,mode,inobs)
   use mod_dimensions  
   implicit none
   real,            intent(out) :: tfin
   type(variances), intent(out) :: var
   real,            intent(out) :: dt
   integer,         intent(out) :: mode
   integer,         intent(out) :: inobs

! The weigths should here be specified as follows:
   open(20,file='infile')
      read(20,*)tfin              ! Final time 
      read(20,*)var%dyn        ! Model error variance
      read(20,*)var%ini        ! Initial error variance
      read(20,*)var%mes        ! Measurement error variance
      read(20,*)inobs          ! Read observations from file mesA.dat
      read(20,*)mode           ! Operation: 0 - pure ensemble int
                               !            1 - ensemble Kalman filter
                               !            2 - ensemble Kalman smoother
                               !            3 - Lagged ensemble Kalman smoother
   close(20)

   dt=(tfin-0.0)/float(ndim-1)

   print '(a,f10.6)','Time step dt=',dt
   print '(a,f10.2)','Final time T=',tfin
   print '(a,g13.4)','var%dyn=',var%dyn
   print '(a,g13.4)','var%ini=',var%ini
   print '(a,g13.4)','var%mes=',var%mes
   
end subroutine physpar
end module m_physpar
