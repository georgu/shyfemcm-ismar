program brown
! EnKF and EnKS program for a scalar red noise model
   implicit none
   character(len=3) caselabel
   integer nrens ! number of particles
   real modvar   ! variance per time unit
   real inivar   ! initial variance
   real obsvar   ! observation variance
   real timeint  ! length of time interval
   integer nstep ! number of time steps
   integer nrdt  ! number of assimilations


   real, parameter :: pi=3.1415927

   real gamma   ! Model bias term
   real bias    ! Estimate model bias or not (1.0 or 0.0)
   real biasvar ! Initial variance in model bias ensemble


   real ave,var,aveq,varq,dt,t,obs,sigma,alpha,rho,new,aveb,varb,inicond
   integer i,j,m,n,nrass

   real, allocatable, dimension(:,:)   ::  D,E,S,X3,X4,X5,ES
   real, allocatable, dimension(:,:)   ::  Astore,Qstore,Bstore

   real, allocatable, dimension(:) :: A,Q,dbeta,workx,worky,B




   open(10, file='infile.in')
      read (10,*)caselabel
      read (10,*)nrens
      read (10,*)timeint
      read (10,*)nstep
      read (10,*)inicond
      read (10,*)inivar
      read (10,*)modvar
      read (10,*)alpha
      read (10,*)obsvar
      read (10,*)nrdt
      read (10,*)gamma
      read (10,*)bias
      read (10,*)biasvar
   close(10)


   print *,'CASE          = ',caselabel
   print *,'Ensemble size = ',nrens
   print *,'Time interval = ',timeint
   print *,'Nr time steps = ',nstep
   print *,'Initial value = ',inicond
   print *,'Initial var   = ',inivar
   print *,'Model var     = ',modvar
   print *,'Alphas        = ',alpha
   print *,'Meas var      = ',obsvar
   print *,'Meas times    = ',nrdt
   print *,'Gamma         = ',gamma
   print *,'bias          = ',bias
   print *,'bias var      = ',biasvar
   nrass=nint(float(nstep)/float(nrdt+1))
   print *,'Ass step      = ',nrass

   dt=timeint/float(nstep)
   print *,'dt            = ',dt

   print *,'model error time scale [tau=dt/(1-alpha)] = ',dt/(1.0-alpha)

   allocate(A(nrens), Q(nrens), dbeta(nrens),workx(nrens),worky(nrens), B(nrens))
   allocate(D(1,nrens), E(1,nrens), S(1,nrens), ES(1,nrens))
   allocate(X3(1,nrens))
   allocate(X4(nrens,nrens))
   allocate(X5(nrens,nrens))

   allocate(Astore(nrens,0:nstep), Qstore(nrens,0:nstep), Bstore(nrens,0:nstep))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Intitialize ensemble to N(inicond,inivar)
   call random_number(workx)
   call random_number(worky)
   A = sqrt(-2.0*log(workx))*cos(2.0*pi*worky)
   ave=sum(A)/float(nrens)
   A=A-ave
   var=dot_product(A,A)/float(nrens-1)
   A=sqrt(1.0/var)*A
   A=sqrt(inivar)*A
   A=A+inicond

! Q
   call random_number(workx)
   call random_number(worky)
   Q = sqrt(-2.0*log(workx))*cos(2.0*pi*worky)

! B
   call random_number(workx)
   call random_number(worky)
   B = sqrt(-2.0*log(workx))*cos(2.0*pi*worky)
   aveb=sum(B)/float(nrens)
   B=B-aveb
   B=sqrt(biasvar)*B

! Store each time step for smoother
   Astore(:,0)=A(:)
   Qstore(:,0)=Q(:)
   Bstore(:,0)=B(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   m=0
   n=nint(float(nstep)/timeint)
   t=0.0
   obs=99.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diagnostics for EnKF
   open(10, file='sol'//caselabel//'.dat')
   ave=sum(A)/float(nrens)
   var=dot_product(A-ave,A-ave)/float(nrens-1)
   aveq=sum(Q)/float(nrens)
   varq=dot_product(Q-aveq,Q-aveq)/float(nrens-1)
   aveb=sum(B)/float(nrens)
   varb=dot_product(B-aveb,B-aveb)/float(nrens-1)
   write(10,'(99f10.3)') t,obs,                                 &
                         ave,   ave-sqrt(var),   ave+sqrt(var), &
                         aveq, aveq-sqrt(varq), aveq+sqrt(varq), &
                         aveb, aveb-sqrt(varb), aveb+sqrt(varb)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computing rho
   rho=sqrt( (1.0-alpha)**2 / (dt*(float(n) - 2.0*alpha -float(n)*alpha**2 + 2.0*alpha**(n+1))) )
   print *, 'rho/n/dt/alpha= ',rho,n,dt,alpha

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Time step loop 
   do j=1,nstep
! Model equations
      call random_number(workx)
      call random_number(worky)
      dbeta = sqrt(-2.0*log(workx))*cos(2.0*pi*worky)
      Q=alpha*Q+sqrt(1.0-alpha**2)*dbeta

      call random_number(workx)
      call random_number(worky)
      dbeta = sqrt(-2.0*log(workx))*cos(2.0*pi*worky)
      B=B+sqrt(0.1*modvar)*sqrt(dt)*dbeta

      A=A+(gamma+B*bias)*dt+sqrt(modvar)*sqrt(dt)*rho*Q 

! Store each time step for smoother
      Astore(:,j)=A(:)
      Qstore(:,j)=Q(:)
      Bstore(:,j)=B(:)

! Diagnostics for EnKF
      ave=sum(A)/float(nrens)
      var=dot_product(A-ave,A-ave)/float(nrens-1)
      aveq=sum(Q)/float(nrens)
      varq=dot_product(Q-aveq,Q-aveq)/float(nrens-1)
      aveb=sum(B)/float(nrens)
      varb=dot_product(B-aveb,B-aveb)/float(nrens-1)
      t=t+dt
      obs=10.0
      write(10,'(99f10.3)') t,obs,                                 &
                            ave,   ave-sqrt(var),   ave+sqrt(var), &
                            aveq, aveq-sqrt(varq), aveq+sqrt(varq),&
                            aveb, aveb-sqrt(varb), aveb+sqrt(varb)


! Assimilation step
      if (mod(j,nrass) == 0) then
         m=m+1
         print *,'Analysis',t, j, m

! create observation with std dev = sqrt(obsvar) and outliers eliminated
         do i=1,10
            call random_number(workx(1))
            call random_number(worky(1))
            dbeta(1) = sqrt(-2.0*log(workx(1)))*cos(2.0*pi*worky(1))
            if (abs(dbeta(1)) < 2.0) exit
         enddo
         obs=inicond+sqrt(obsvar)*dbeta(1)

! create ensemble of observation perturbations
         call random_number(workx)
         call random_number(worky)
         dbeta = sqrt(-2.0*log(workx))*cos(2.0*pi*worky)
         ave=sum(dbeta)/float(nrens)
         dbeta=dbeta-ave
         var=dot_product(dbeta,dbeta)/float(nrens-1)
         dbeta=sqrt(1.0/var)*dbeta
         E(1,:)=sqrt(obsvar)*dbeta(:)

! compute ensemble of innovations
         D(1,:)=obs+E(1,:)-A(:)

! compute S matrix (HA')
         ave=sum(A)/float(nrens)
         S(1,:)=A(:)-ave

! standard EnKF analysis for one single observation
         ES=E+S
!         sigma=dot_product(ES(1,:),ES(1,:))  ! old approximate algorithm from Evensen 2003
         sigma=dot_product(E(1,:),E(1,:))+dot_product(S(1,:),S(1,:))

         X3=D/sigma

         X4=matmul(transpose(S),X3)

         X5=X4
         do i=1,nrens
            X5(i,i)=X5(i,i)+1.0
         enddo

! update A, Q and B
         A=matmul(A,X5)
         Q=matmul(Q,X5)
         B=matmul(B,X5)

! compute EnKS analysis from time 0 to present time.
         do i=0,j
            Astore(:,i)=matmul(Astore(:,i),X5)
            Qstore(:,i)=matmul(Qstore(:,i),X5)
            Bstore(:,i)=matmul(Bstore(:,i),X5)
         enddo

! Diagnostics for EnKF
         ave=sum(A)/float(nrens)
         var=dot_product(A-ave,A-ave)/float(nrens-1)
         aveq=sum(Q)/float(nrens)
         varq=dot_product(Q-aveq,Q-aveq)/float(nrens-1)
         aveb=sum(B)/float(nrens)
         varb=dot_product(B-aveb,B-aveb)/float(nrens-1)
         write(10,'(99f10.3)') t,obs,                                 &
                               ave,   ave-sqrt(var),   ave+sqrt(var), &
                               aveq, aveq-sqrt(varq), aveq+sqrt(varq),&
                               aveb, aveb-sqrt(varb), aveb+sqrt(varb)
      endif
   enddo
   close(10)   


! Diagnostics for smoother solution
   t=0.0
   new=sum(Astore(:,0))/float(nrens)-0.1
   open(10, file='smo'//caselabel//'.dat')
   do j=0,nstep
      ave=sum(Astore(:,j))/float(nrens)
      var=dot_product(Astore(:,j)-ave,Astore(:,j)-ave)/float(nrens-1)
      aveq=sum(Qstore(:,j))/float(nrens)
      varq=dot_product(Qstore(:,j)-aveq,Qstore(:,j)-aveq)/float(nrens-1)
      aveb=sum(Bstore(:,j))/float(nrens)
      varb=dot_product(Bstore(:,j)-aveb,Bstore(:,j)-aveb)/float(nrens-1)

      write(10,'(99f10.3)') t,new,                                 &
                            ave,   ave-sqrt(var),   ave+sqrt(var), &
                            aveq, aveq-sqrt(varq), aveq+sqrt(varq),&
                            aveb, aveb-sqrt(varb), aveb+sqrt(varb)
      t=t+dt
      if (j<nstep) then
         aveq=sum(Qstore(:,j+1))/float(nrens)
         aveb=sum(Bstore(:,j+1))/float(nrens)
      endif
      new=new+(gamma+bias*aveb)*dt+sqrt(dt)*sqrt(modvar)*rho*aveq
   enddo
   close(10)



end program

