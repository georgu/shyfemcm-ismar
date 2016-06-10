program assim 
   use mod_dimensions       ! Parameter and type  declarations
   use m_eksRfact
   use m_eksRfact_lag
   use m_eksreps
   use m_eksreps_lag
   use m_eksscalc
   use m_eksscalc_lag
   use m_kfilt
   use m_lorenz
   use m_obs
   use m_physpar
   use m_random
   use m_random2
   use m_ranmean
   use m_ranvar

   implicit none
   type(variances) vars

   integer itsol             ! start point for iteration
   integer, parameter :: lag=75
   real A(3,ndim,nrsamp)
   real AA(3,nrsamp)
   real ave(3,ndim)          ! predicted state variable
   real var(3,ndim)          ! predicted state variable
   real cov(3,ndim)          ! predicted state variable
   real corr(3,ndim)         ! predicted state variable
   real sol(3,ndim)          ! current accepted state variable
   real ref(3,ndim)          ! Reference solution
   real cfc(3,ndim)          ! Central forecast solution
   real d(3,nrmes)           ! vector of measurements
  
   real bcoeff(nrtobs,0:nrsamp)   ! data misfit
   real S(nrtobs,nrsamp)     ! H(A-M) matrix
   real RR(nrtobs,nrtobs)    ! Smoother representer matrix
   real Reps(nrtobs,3)    ! Matrix with influence functions
   real h(1,1)            ! h(nrmes,ndim) Linear functional L[w]
   integer ih(nrmes)      ! Measurement pointers
   real v(nrmes,nrmes)    ! Measurement error covariance matrix
   real v2(nrmes,nrmes)   ! work array
   real t,dt,tfin,t1,t2
   integer na,nb
   real zini(3)
   real deltaobs
   real covar(3,3)    
   integer i,j,k,m,n,nc,ican,p,ipp,jgm,nc2,mm
   integer irejec,iave
   integer nsing
   integer mode
   integer inobs
   character*1 tag

   real bcoeff_lag(nrobs,0:nrsamp)   ! data misfit
   real S_lag(nrobs,nrsamp)     ! H(A-M) matrix
   real RR_lag(nrobs,nrobs)    ! Smoother representer matrix
   real Reps_lag(nrobs,3)    ! Matrix with influence functions

   call random(sol)
   do i=1,ndim
      write(20,*)i,sol(1:3,i)
   enddo

! Reading initial parameters
   call physpar(tfin,vars,dt,mode,inobs)

! Simulating reference case
   ref(1,1)= 1.508870
   ref(2,1)=-1.531271
   ref(3,1)= 25.46091
   na=1
   nb=ndim
   call lorenz(na,nb,dt,ref,0.0)


! Defining initial cond for central forecast
   call random2(cfc(1,1),3)
   cfc(:,1)=sqrt(vars%ini)*cfc(:,1)+ref(:,1)


! Generating initial ensemble
   do i=1,nrsamp
      call random2(A(1,1,i),3)
      A(:,1,i)=cfc(:,1)+sqrt(vars%ini)*A(:,1,i)
   enddo


! Generating measurements
   call obs(ref,d,ih,v,dt,vars,deltaobs,inobs)



! Time stepping
   do m=1,nrmes

      na=nint(deltaobs*float(m-1)/dt+1.0)
      nb=nint(deltaobs*float(m)/dt+1.0)
      print '(2(a,i4))','na=',na,' nb=',nb

      call lorenz(na,nb,dt,cfc,0.0) 
      do j=1,nrsamp
         call lorenz(na,nb,dt,A(1,1,j),vars%dyn) 
      enddo

      if (mode == 1) then
         ! Kalman filter analysis
         AA(:,:)=A(:,nb,:)
         call kfilt(AA,d(1,m),cfc(1,nb),vars)
         A(:,nb,:)=AA(:,:)

      elseif (mode == 2) then
         ! Accumulating S and residuals
         nb=nint(deltaobs*float(m)/dt+1.0)
         AA(:,:)=A(:,nb,:)
         call eksscalc(cfc,AA,S,bcoeff,d,m,vars)

      elseif (mode == 3) then
         ! Accumulating S_lag and residuals
         nb=nint(deltaobs*float(m)/dt+1.0)
         AA(:,:)=A(:,nb,:)
         call eksscalc_lag(cfc,AA,S_lag,bcoeff_lag,d,m,vars)
         print *,'eksscalc_lag ok'

         print *,'call eksRfact'
         call eksRfact_lag(S_lag,RR_lag,bcoeff_lag,nsing,vars%mes)
         print *,'eksRfact_lag ok'

         do i=max(1,nb-lag),nb
            AA(:,:)=A(:,i,:)
            call eksreps_lag(S_lag,AA,RR_lag,Reps_lag)
            do mm=1,nrobs
               cfc(:,i)=cfc(:,i)+bcoeff_lag(mm,0)*Reps_lag(mm,:)
               do j=1,nrsamp
                  A(:,i,j)=A(:,i,j)+bcoeff_lag(mm,j)*Reps_lag(mm,:)
               enddo
            enddo
         enddo

      endif

   enddo

   if (mode == 2) then   ! Smoother analysis
      print *,'call eksRfact'
      call eksRfact(S,RR,bcoeff,nsing,vars%mes)
      print *,'eksRfact ok'

      do i=1,ndim
         print *,'i=',i,ndim
         AA(:,:)=A(:,i,:)
         call eksreps(S,AA,RR,Reps)
         do m=nrtobs,nrtobs-nsing+1,-1
            cfc(:,i)=cfc(:,i)+bcoeff(m,0)*Reps(m,:)
            do j=1,nrsamp
               A(:,i,j)=A(:,i,j)+bcoeff(m,j)*Reps(m,:)
            enddo
         enddo
      enddo
   endif

   call ranmean(A,ave)
   call ranvar(A,ave,var)
   write(tag,'(i1.1)')mode

   open(10,file='sol'//tag//'.dat')
      do i=1,ndim
         write(10,'(16g12.4)')float(i-1)*dt,ref(1:3,i),&
                                            cfc(1:3,i),&
                                            ave(1:3,i),&
                                            var(1:3,i),&
                                            (ave(1:3,i)-ref(1:3,i))**2
      enddo
   close(10)

   open(10,file='mes'//tag//'.dat')
      do i=1,nrmes
         write(10,*)float(i)*deltaobs, d(1:3,i), ih(i)
      enddo
   close(10)

end
