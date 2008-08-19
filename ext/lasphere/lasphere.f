      program lahyc
c*******************************************************
c                                                      *
c  This is a Lagrangian 1D hydrodynamics code.         *
c                                                      *
c  this is the main driver of the integration.         *
c                                                      *
c*******************************************************
c
      implicit double precision (a-h,o-z)
c
      common /steps/ maxstep,ifreq
      common /tim / time, dt, dtold, tmax
c
c--read material data
c
      call readmatdata
c
c--read initial condition 
c
      call initial
c
c--main loop 
c
      nstep=0
   10 continue
      nstep=nstep + 1
c
c--find time step using stability criterion
c
      call tstep
c
c--advance the system of hydro equation by one time step
c
      call hydro
c
c--produce output if desired
c
      if(mod(nstep,ifreq).eq.0)then
         call printout(nstep)
      end if
c
c--check if simulation is finished
c
      if(time.gt.tmax.or.nstep.eq.maxstep) go to 90
c
      go to 10
c
   90 continue
      stop
      end
      subroutine hydro
c****************************************************************
c                                                               *
c  this subroutine advances the system of hydro equations by    *
c  one time step.                                               *
c                                                               *
c  Description of the variables :                               *
c  ------------------------------                               *
c                                                               *
c  ncell        :  number of cells                              *
c  ncell1       :  number of cell edges                         *
c  rho          :  density (cell centered)                      *
c  p            :  pressure (cell centered)                     *
c  u            :  specific internal energy (cell centered)     *
c  q            :  artificial viscous stress (cell centered)    *
c  deltam       :  mass of the cell (cell centered)             *
c  pold,qold,rhold  : pressure artificial viscous stress and    *
c                     density at previous time step             *
c  x            :  cell boundaries (edge centered)              *
c  v            :  velocity (edge centered)                     *
c  vold         :  velocity at previous time step               *
c  gamma        :  ratio of specific heat                       *
c  cq, cl       :  quadratic and linear artificial viscous      *
c                  stress coefficients                          *
c  tkin, uint   :  kinetic and specific internal energy         *
c  time         :  current time                                 *
c  dt, dtold    :  current and previous time step               *
c  tmax         :  maximum time for the simulation              *
c                                                               *
c****************************************************************
c
      implicit double precision (a-h,o-z)
c
      parameter (ncell=1000)
      parameter (ncell1=ncell + 1)
c
      common /cellc/ rho(0:ncell1), p(0:ncell1), u(0:ncell1), 
     1               q(0:ncell1), deltam(0:ncell1), pold(0:ncell1),
     2               qold(0:ncell1), rhold(0:ncell1), mat(0:ncell1)
      common /celle/ x(0:ncell1), xold(0:ncell1), v(0:ncell1), 
     1               vold(0:ncell1)
      common /frict/ gam
      common /shock/ cq,cl
      common /energ/ tkin, uint
      common /tim / time, dt, dtold, tmax
c
      logical first
c
      data first/.true./
      data pi4/12.56637/
      data pi43/4.1887902/
      data gg/6.67d-8/
c
c--update time
c
      time= time + dt
      if(first) then
         dtn=dt
         first=.false.
      else
         dtn=0.5*(dt+dtold) 
      end if
c
c--update velocity using momentum equation
c
      v(0)=0.d0
      xmsum=0.d0
      do k=1,ncell
         km05=k - 1
         kp05=k
         ak=pi4*x(k)**2
         akp1=pi4*x(k+1)**2
         akm1=pi4*x(k-1)**2
         akp05=0.5d0*(akp1+ak)
         akm05=0.5d0*(ak+akm1)
         vold(k)=v(k)
         deltamk=0.5d0*(deltam(km05)+deltam(kp05))
         gradp=ak*(p(kp05) - p(km05))
         gradq=0.5d0*q(kp05)*(3*akp05-ak) - 0.5*q(km05)*(3*akm05-ak)
         xmsum=xmsum + deltam(km05)
         grav=gg*xmsum/(x(k)*x(k))
         v(k)=v(k) - dtn*(gradp + gradq)/deltamk - dtn*grav - 
     1        dtn*gam*v(k)
      enddo
      v(ncell1)=v(ncell)
c
c--detemine new coordinate using new velocity
c
      do k=0,ncell1
         xold(k)=x(k)
         x(k)=x(k) + v(k)*dt
      enddo
c
c--update density using continuity equation
c
      do kp05=0,ncell
         k1=kp05 + 1
         k=kp05 
         rhold(kp05)=rho(kp05)
         rho(kp05)=3.*deltam(kp05)/(pi4*(x(k1)**3-x(k)**3)) 
      enddo
c
c--update q value
c
      do kp05=0,ncell
         k=kp05
         k1=kp05 + 1
         qold(kp05)=q(kp05)
         q(kp05)=0.
         akp1=pi4*x(k1)**2
         ak=pi4*x(k)**2
         akp05=0.5*(akp1+ak)
         gradv=(v(k1)*(akp05-akp1/3.) - v(k)*(akp05-ak/3.))
c
c--compute gradient of velocity
c
         if(gradv.lt.0.)then
            dv=v(k1) - v(k)
            rhonp05=0.5*(rho(kp05)+rhold(kp05))
            alpha=1.5*cq**2*rhonp05*abs(dv)/akp05
            q(kp05)=-alpha*gradv
c
c--linear term
c
            pk=0.5*(p(kp05)+pold(kp05))
            cs=sqrt(1.6666*abs(pk)/rhonp05)
            alpha=1.5*cl*rhonp05*cs/akp05
            q(kp05)=q(kp05) - alpha*gradv
         end if
      enddo
c
c--update energy value by computing temporary pressure first
c
      do kp05=0,ncell
         k=kp05
         k1=kp05 + 1
         akp1=pi4*x(k1)**2
         ak=pi4*x(k)**2
         akp05=0.5*(akp1+ak)
         dvol=(akp1*v(k1)-ak*v(k))*dt/deltam(kp05)
c        dvol=pi43*((x(k1)**3-x(k)**3)-(xold(k1)**3-xold(k)**3))
         uprime=u(kp05) - p(kp05)*dvol
         matk=mat(kp05)
         rhok=rho(kp05)
         call eos(matk,rhok,uprime,pprime)
         pnp05=0.5*(p(kp05)+pprime)
         pdv=pnp05*(akp1*v(k1) - ak*v(k))
         qnp05=0.5*(q(kp05)+qold(kp05))
         dq=0.5*qnp05*(v(k1)*(3.*akp05-akp1) - v(k)*(3.*akp05-ak))
         u(kp05)=u(kp05) - dt*(pdv+dq)/deltam(kp05)
c        u(kp05)=u(kp05) - dt*(pdv)/deltam(kp05)
      enddo
c
c--update pressure
c
      do kp05=0,ncell
         pold(kp05)=p(kp05)
         matk=mat(kp05)
         uk=u(kp05)
         rhok=rho(kp05)
         call eos(matk,rhok,uk,press)
         p(kp05)=press
      enddo
c
c--compute energies
c
      tkin=0.
      uint=0.
      do kp05=1,ncell
         k=kp05
         k1=kp05 + 1
         deltamk1=0.5*(deltam(kp05)+deltam(kp05+1))
         deltamk=0.5*(deltam(kp05)+deltam(kp05-1))
         tnk05=0.25*(deltamk1*v(k1)*vold(k1)+deltamk*v(k)*vold(k))
         tkin=tkin + tnk05
         uint=uint + u(kp05)*deltam(kp05)
      enddo
c
      return
      end
      subroutine eos (imat,rho,u,pr)
c************************************************
c                                               *
c           material constants                  *
c                                               *
c   imat  material                              *
c   ----  --------                              *
c                                               *
c    1   granite        11   dry tuff           *
c    2   basalt         12   alluvium           *
c    3   aluminum       13   anorthosite 1pp    *
c    4   copper         14   anorthosite hpp    *
c    5   iron 130pt     15   andesite           *
c    6   lucite         16   water              *
c    7   limestone      17   pure ice           *
c    8   halite         18   5% silicate ice    *
c    9   oil shale      19   30% silicate ice   *
c   10   wet tuff       20   special            *
c                                               *
c************************************************
c
      implicit double precision (a-h,o-z)
c
      parameter (nmatmax=20)
c
      common /till / rozero(nmatmax), smalla(nmatmax), uzero(nmatmax), 
     1               smallb(nmatmax), capa(nmatmax), es(nmatmax), 
     2               esp(nmatmax), alpha(nmatmax), beta(nmatmax), 
     3               capb(nmatmax)
c
      data prlim/-1.d15/, tiny/1.d-30/
c
c--define constants
c
      rozn=rozero(imat)
      ezn=uzero(imat)
      capan=capa(imat)
      esn=es(imat)
      espn=esp(imat)
      difesn=espn-esn
      capbn=capb(imat)
      smallaj=smalla(imat)
      smallbj=smallb(imat)
      alphaj=alpha(imat)
      betaj=beta(imat)
c
      rho2=rho*rho
      eta=rho/rozn
      xmu=eta-1.d0
      rap=rozn/rho
c
c--condensed phase, compute pressure and sound speed
c
      p1=smallbj/(u/(ezn*eta*eta)+1.d0)
      pressc=(smallaj+p1)*u*rho+capan*xmu+capbn*xmu*xmu
      p1=u/(ezn*eta*eta)+1.d0
      c2=p1*p1
      c3=p1-1.d0
      c4=smallbj*u*(3.d0*c3+1.d0)/c2
      dpdi=smallaj*rho+smallbj*rho/c2
      dpdrho=smallaj*u+capan/rozn+2.*capbn*xmu/rozn+c4
      cvelc=dpdrho+dpdi*pressc/rho2
c
c--expanded phase, compute pressure and sound speed
c
      p1=u/ezn/eta**2+1.d0
      p2=smallaj*u*rho
      p3=smallbj*u*rho/p1
      p4=capan*xmu
      vow=1./eta
      p6=betaj*(vow-1.d0)
      p6=min(p6,70.d0)
      p6=exp(-p6)
      p7=alphaj*(vow-1.d0)**2
      p7=min(p7,70.d0)
      p7=exp(-p7)
      pressv=p2+(p3+p4*p6)*p7
      c1=smallaj*u
      c2=p1**2
      c3=p1-1.d0
      c4=smallbj*u*(3.*c3+1.d0)/c2
      c5=p6*p7*capan
      c6=2.*alphaj*(vow-1.d0)
      dpdrho=c1+p7*c4+p7*p3*rozn*c6/rho2+c5*(1.d0/rozn+xmu*rozn/
     1       rho2*(c6+betaj))
      dpdi=smallaj*rho+p7/c2*smallbj*rho
      cvelv=dpdrho+dpdi*pressv/rho2
      if (cvelv.lt.0.d0) cvelv=0.
c
c--pick up the right state
c
      pr=pressc
      vsound=cvelc
      if(rap.ge.1.d0.and.u.gt.espn) then
         pr=pressv
         vsound=cvelv
      end if
      if(rap.ge.1.d0.and.u.gt.esn.and.u.lt.espn) then
         pr=(pressv*(u-esn) + pressc*(espn-u))/difesn
         vsound=(cvelv*(u-esn) + cvelc*(espn-u))/difesn
      end if
      cmin=.25d0*capan/rozn
      vsound=max(cmin,vsound)
      pr=max(pr,prlim)
c
      return
      end
      subroutine tstep
c****************************************************************
c                                                               *
c  This subroutine computes the time step according to the      *
c  Courant condition and using the criterion on the artificial  *
c  viscosity.                                                   *
c                                                               *
c****************************************************************
c
      implicit double precision (a-h,o-z)
c
      parameter (ncell=1000)
      parameter (ncell1=ncell + 1)
c
      common /cellc/ rho(0:ncell1), p(0:ncell1), u(0:ncell1), 
     1               q(0:ncell1), deltam(0:ncell1), pold(0:ncell1),
     2               qold(0:ncell1), rhold(0:ncell1), mat(0:ncell1)
      common /celle/ x(0:ncell1), xold(0:ncell1), v(0:ncell1), 
     1               vold(0:ncell1)
      common /shock/ cq,cl
      common /tim / time, dt, dtold, tmax
c
      data c1/0.2/
c
c--compute the largest inverse time step over the whole grid 
c
      dt1=0.d0
      do kp05=0,ncell
         k=kp05
         k1=kp05 + 1
         dx=x(k1) - x(k)
c
c--viscous time step
c
         dtq1=4.d0*cq**2*abs((v(k1)-v(k))/dx)
c
c--courant condition
c
         cs=sqrt(1.66666*abs(p(kp05))/rho(kp05))
         dtc1=(cs+0.5d0*abs(v(k1)+v(k)))/dx
c
c--find maximum
c
         dt1=max(dt1,dtq1,dtc1)
       enddo
c
c--define time step as being the smallest
c
      dtold=dt
      dt=c1/dt1
c
      return
      end
      subroutine initial
c***************************************************************
c                                                              *
c  this subroutine reads the initial values for all quantities *
c                                                              *
c***************************************************************
c
      implicit double precision (a-h,o-z)
c
      parameter (ncell=1000)
      parameter (ncell1=ncell + 1)
c
      common /cellc/ rho(0:ncell1), p(0:ncell1), u(0:ncell1), 
     1               q(0:ncell1), deltam(0:ncell1), pold(0:ncell1),
     2               qold(0:ncell1), rhold(0:ncell1), mat(0:ncell1)
      common /celle/ x(0:ncell1), xold(0:ncell1), v(0:ncell1), 
     1               vold(0:ncell1)
      common /frict/ gam
      common /rec  / irec
      common /shock/ cq,cl
      common /steps/ maxstep,ifreq
      common /tim / time, dt, dtold, tmax
c
      data pi43/4.1887902/
c
      irec=0
c
c--set artificial viscosity constants
c
      cq=2.d0
      cl=1.d0
c
c--read initial conditions
c
      open(15,file='inlasphere')
      read(15,*)xmtot
      read(15,*)imat1,f1,rho1,u1
      read(15,*)imat2,f2,rho2,u2
      read(15,*)tmax,maxstep,ifreq
      read(15,*)gam
c
c--setup initial conditions
c
c  a) masses, densities, energies and material type
c
      n1=f1*ncell
      write(*,*)n1,imat1
      if(n1.ne.0)then
         deltam1=f1*xmtot/n1
         do i=0,n1
            deltam(i)=deltam1
            rho(i)=rho1
            u(i)=u1
            mat(i)=imat1
            q(i)=0.0d0
         enddo
      end if
      n2=ncell-n1
      write(*,*)n2,imat2
      if(n2.ne.0)then
         deltam2=f2*xmtot/n2
         if(n1.eq.0)then
            ndep=0
         else
            ndep=n1+1
         end if
         do i=ndep,ncell1
            deltam(i)=deltam2
            rho(i)=rho2
            u(i)=u2
            mat(i)=imat2
            q(i)=0.0d0
         enddo
      end if
c
c  b) position, velocities
c
      x(0)=0.d0
      v(0)=0.d0
      do i=1,n1
         x(i)=(deltam(i-1)/(pi43*rho1) + x(i-1)**3)**0.3333333333
         v(i)=0.d0
      enddo
      do i=n1+1,ncell1
         x(i)=(deltam(i-1)/(pi43*rho2) + x(i-1)**3)**0.3333333333
         v(i)=0.d0
      enddo
      write(*,*)x(ncell1)
c
c  c) pressure
c
      do i=0,ncell
         mati=mat(i)
         rhoi=rho(i)
         ui=u(i)
         call eos(mati,rhoi,ui,press)
         p(i)=press
      enddo
c
      return
      end
      subroutine printout (nstep)
c**************************************************************
c                                                             *
c  This subroutine prints out all the results.                *
c                                                             *
c**************************************************************
c 
      implicit double precision (a-h,o-z)
c
      parameter (ncell=1000)
      parameter (ncell1=ncell + 1)
c
      common /cellc/ rho(0:ncell1), p(0:ncell1), u(0:ncell1), 
     1               q(0:ncell1), deltam(0:ncell1), pold(0:ncell1),
     2               qold(0:ncell1), rhold(0:ncell1), mat(0:ncell1)
      common /celle/ x(0:ncell1), xold(0:ncell1), v(0:ncell1), 
     1               vold(0:ncell1)
      common /rec  / irec
      common /shock/ cq,cl
      common /energ/ tkin, uint
      common /tim / time, dt, dtold, tmax
c
      character*8 filename
      character*5 file
      character*3 number
c
      data file/'body-'/
c
c--open appropriate file
c
      irec=irec + 1
      write(number,200)irec
  200 format(i3.3)
      filename=file // number
      open(21,file=filename)
c
c--print general quantities first
c
      etot=uint + tkin
      write(*,100)nstep,time,dt,etot,tkin,uint
  100 format(1x,i5,' time steps. Current time : ',1pd12.5,
     1' current dt : ',1pd12.5,/,
     2' total energy : ',1pd12.5,' kinetic energy : ',1pd12.5,
     3' internal energy : ',1pd12.5,/)
c
c--print all cell quantities at cell centers
c
c     write(21,101)
c 101 format(' zone     x          v            p          q',
c    1'           rho          u')
      do kp05=0,ncell
         k=kp05
         k1=kp05 + 1
         xc=0.5*(x(k1)+x(k))
         vc=0.5*(v(k1)+v(k))
         write(21,102)kp05,xc,vc,p(kp05),q(kp05),rho(kp05),u(kp05),
     1                mat(kp05)
  102    format(1x,i4,1x,6(1pd11.4,1x),i3)
      enddo
c
      return
      end
      subroutine readmatdata
c************************************************
c                                               *
c  subroutine reading material data             * 
c                                               *
c************************************************
c
      implicit double precision (a-h,o-z)
c
c      parameter (nmatmax=20)
      parameter (nmatmax=21)
c
      common /till / rozero(nmatmax), smalla(nmatmax), uzero(nmatmax), 
     1               smallb(nmatmax), capa(nmatmax), es(nmatmax), 
     2               esp(nmatmax), alpha(nmatmax), beta(nmatmax), 
     3               capb(nmatmax)
      common /shear/ xmu(nmatmax)
      common /yield/ umelt(nmatmax), yie(nmatmax)
      common /fracm/ pweib(nmatmax), cweib(nmatmax)
c
      character*50 dumy
c
c--open data file
c
c      open(30,file='../material-constants.dat')
      open(30,file='tillotson.input')
c
c--skip headers
c
      do i=1,12
         read(30,100)dumy
  100    format(a50)
      enddo
c
c--read number of materials and lines
c
      read(30,*)nmat
      im=mod(nmat,5)
      nline=nmat/5 
      if(im.ne.0)nline=nline + 1
c
c--read data
c
c  i) reference density rozero
c
      read(30,100)dumy
      read(30,100)dumy
      read(30,100)dumy
      imat=0
      do i=1,nline
         imatmin=(i-1)*5 + 1
         imatmax=min(nmat,imatmin+4)
         read(30,*)(rozero(j),j=imatmin,imatmax)
      enddo
c
c  ii) bulk modulus capa
c
      read(30,100)dumy
      read(30,100)dumy
      read(30,100)dumy
      imat=0
      do i=1,nline
         imatmin=(i-1)*5 + 1
         imatmax=min(nmat,imatmin+4)
         read(30,*)(capa(j),j=imatmin,imatmax)
      enddo
c
c  iii) nonlinear compressive term capb
c
      read(30,100)dumy
      read(30,100)dumy
      read(30,100)dumy
      imat=0
      do i=1,nline
         imatmin=(i-1)*5 + 1
         imatmax=min(nmat,imatmin+4)
         read(30,*)(capb(j),j=imatmin,imatmax)
      enddo
c
c  iv) smalla
c
      read(30,100)dumy
      read(30,100)dumy
      read(30,100)dumy
      imat=0
      do i=1,nline
         imatmin=(i-1)*5 + 1
         imatmax=min(nmat,imatmin+4)
         read(30,*)(smalla(j),j=imatmin,imatmax)
      enddo
c
c  v) smallb
c
      read(30,100)dumy
      read(30,100)dumy
      read(30,100)dumy
      imat=0
      do i=1,nline
         imatmin=(i-1)*5 + 1
         imatmax=min(nmat,imatmin+4)
         read(30,*)(smallb(j),j=imatmin,imatmax)
      enddo
c
c  vi) alpha
c
      read(30,100)dumy
      read(30,100)dumy
      read(30,100)dumy
      imat=0
      do i=1,nline
         imatmin=(i-1)*5 + 1
         imatmax=min(nmat,imatmin+4)
         read(30,*)(alpha(j),j=imatmin,imatmax)
      enddo
c
c  vii) beta
c
      read(30,100)dumy
      read(30,100)dumy
      read(30,100)dumy
      imat=0
      do i=1,nline
         imatmin=(i-1)*5 + 1
         imatmax=min(nmat,imatmin+4)
         read(30,*)(beta(j),j=imatmin,imatmax)
      enddo
c
c  viii) uzero
c
      read(30,100)dumy
      read(30,100)dumy
      read(30,100)dumy
      imat=0
      do i=1,nline
         imatmin=(i-1)*5 + 1
         imatmax=min(nmat,imatmin+4)
         read(30,*)(uzero(j),j=imatmin,imatmax)
      enddo
c
c  ix) es (sublimation point)
c
      read(30,100)dumy
      read(30,100)dumy
      read(30,100)dumy
      imat=0
      do i=1,nline
         imatmin=(i-1)*5 + 1
         imatmax=min(nmat,imatmin+4)
         read(30,*)(es(j),j=imatmin,imatmax)
      enddo
c
c  x) esp
c
      read(30,100)dumy
      read(30,100)dumy
      read(30,100)dumy
      imat=0
      do i=1,nline
         imatmin=(i-1)*5 + 1
         imatmax=min(nmat,imatmin+4)
         read(30,*)(esp(j),j=imatmin,imatmax)
      enddo
c
c  xi) shear modulus xmu
c
c      read(30,100)dumy
c      read(30,100)dumy
c      read(30,100)dumy
c      imat=0
c      do i=1,nline
c         imatmin=(i-1)*5 + 1
c         imatmax=min(nmat,imatmin+4)
c         read(30,*)(xmu(j),j=imatmin,imatmax)
c      enddo
c
c  xii) melt energy umelt
c
c      read(30,100)dumy
c      read(30,100)dumy
c      read(30,100)dumy
c      imat=0
c      do i=1,nline
c         imatmin=(i-1)*5 + 1
c         imatmax=min(nmat,imatmin+4)
c         read(30,*)(umelt(j),j=imatmin,imatmax)
c      enddo
c
c  xiii) yielding yie
c
c      read(30,100)dumy
c      read(30,100)dumy
c      read(30,100)dumy
c      imat=0
c      do i=1,nline
c         imatmin=(i-1)*5 + 1
c         imatmax=min(nmat,imatmin+4)
c         read(30,*)(yie(j),j=imatmin,imatmax)
c      enddo
c
c  xiv) fracture model pweib
c
c      read(30,100)dumy
c      read(30,100)dumy
c      read(30,100)dumy
c      imat=0
c      do i=1,nline
c         imatmin=(i-1)*5 + 1
c         imatmax=min(nmat,imatmin+4)
c         read(30,*)(pweib(j),j=imatmin,imatmax)
c      enddo
c
c  xv) fracture model cweib
c
c      read(30,100)dumy
c      read(30,100)dumy
c      read(30,100)dumy
c    imat=0
c      do i=1,nline
c         imatmin=(i-1)*5 + 1
c         imatmax=min(nmat,imatmin+4)
c         read(30,*)(cweib(j),j=imatmin,imatmax)
c      enddo
c
      close(30)
      return
      end 
