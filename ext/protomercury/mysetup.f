      program setparticles
c*****************************************************************
c                                                                *
c  Program to set up particles using the close-pack algorithm.   *
c                                                                *
c*****************************************************************
c
      implicit double precision (a-h,o-z)
      parameter (idim=200000)
c
      common /part0/ npart, x(idim), y(idim), z(idim), h(idim),
     1               alpha(idim),matter(idim)
      common /part1/ u(idim),temp(idim), rho(idim),pr(idim)
      common /part2/ pmass(idim), vx(idim), vy(idim), vz(idim)
c      common /part3/ sxx(idim), syy(idim), sxy(idim), sxz(idim),
c     1               syz(idim), dm(idim) 
      common /units/ umass, udist, udens, utime, uergg, uergcc    
c
c--read input options
c
      open(40,file='set-particles.in')
      read(40,*)rr,napprox
      close(40)
      write(*,*)'input read: ',rr, napprox
c
      call unit
c
      rmax=rr/udist
      write(*,*) 'rmax:',rmax
c
c--compute particle spacing
c
      call setdist(rmax,napprox,hsep,xmin,xmax,
     1             ymin,ymax,zmin,zmax)
c
c--read the computed tables
c
      call getdata
c
c--set particles on lattice
c
      call setonlattice(rmax,hsep,xmin,xmax,ymin,ymax,
     1                  zmin,zmax)
c
c--set velocities
c
      call setvelocities 
c
c--set deviatoric stress tensor
c
c      call setdeviatoric
c
c--set damage
c
c      call setdamage
c
c--set remaining quantities
c
      n1=npart
      n2=0
      t=0.d0/utime
      trot=0.d0
      tkin=0.d0
      tgrav=0.d0
      tterm=0.d0
      do i=1,npart
         alpha(i)=0.1d0
      enddo
c
c--write output
c
      open(11,file='protomr000',form='unformatted')
      write(11)npart,n1,n2,t,trot,tkin,tgrav,tterm,
     &         (x(i),i=1,npart),
     &         (y(i),i=1,npart),
     &         (z(i),i=1,npart),
     &         (vx(i),i=1,npart),
     &         (vy(i),i=1,npart),
     &         (vz(i),i=1,npart),
     &         (u(i),i=1,npart),
     &         (h(i),i=1,npart),
     &         (pmass(i),i=1,npart),
     &         (rho(i),i=1,npart),
     &         (temp(i),i=1,npart),
     &         (pr(i),i=1,npart),
     &         (alpha(i),i=1,npart),
     &         (matter(i),i=1,npart)
      close(11)
c 
c      do i=1,npart
c      if (rho(i) .lt. 0.d0) then
c      write(*,*) 'rho(i):',rho(i)
c      endif
c      enddo
c
      stop
      end
      subroutine setdist(rmax,napprox,hsep,xmin,xmax,
     1                   ymin,ymax,zmin,zmax)
c*******************************************************************
c                                                                  *
c  computes the spacing between the particles                      *
c                                                                  *
c*******************************************************************
c
      implicit double precision (a-h,o-z)
      pi=acos(-1.d0)
      third=1.d0/3.d0
c
      volume=4.d0*third*pi*rmax**3.d0
      hsep=(volume/napprox)**third
      xmin=-rmax
      xmax=rmax
      ymin=-rmax
      ymax=rmax
      zmin=-rmax
      zmax=rmax
      write(*,*) 'hsep', hsep
c
      return
      end
      subroutine inside(ri,xi,yi,zi,rmax,iflg)
c************************************************************
c                                                           *
c  iflg=0 if particle is outside object                     *
c  iflg=1 if particle is inside  object                     *
c                                                           *
c************************************************************
c
      implicit double precision (a-h,o-z)
      iflg=0
c
      ri=sqrt(xi*xi + yi*yi + zi*zi)
      ratio=ri*ri/(rmax*rmax)
      if(ratio.lt.1.d0)then
         iflg=1
      endif
c
      return
      end
      subroutine setonlattice(rmax,hsep,xmin,xmax,ymin,ymax,
     1                        zmin,zmax)
c****************************************************************
c                                                               *
c  set particles on a close-packed lattice                      *
c                                                               *
c****************************************************************
c
      implicit double precision (a-h,o-z)
      parameter (idim=200000)
c
      common /part0/ npart, x(idim), y(idim), z(idim), h(idim),
     1               alpha(idim), matter(idim)
      common /part1/ u(idim),temp(idim),rho(idim),pr(idim)
      common /part2/ pmass(idim), vx(idim), vy(idim), vz(idim)
c
c--set particles on lattice
c
c      hsep=hsep*2.d0**(1.d0/6.d0)
      dx=0.5d0*hsep
      dy=hsep/sqrt(12.d0)
      dz=hsep*sqrt(2.d0/3.d0)
c
      xstep=0.5d0*hsep
      ystep=xstep*sqrt(3.d0)
      zstep=2.d0*dz
c
c--first layer
c
      npart=0
      factor=0.5d0
      bfct=1.d0/sqrt(2.d0)
      do y0=ymin,ymax,ystep
         xmin=xmin + factor*hsep
         factor=-factor
         do x0=xmin,xmax,hsep
            do z0=zmin,zmax,zstep
               xi=x0
               yi=y0
               zi=z0
               call inside(ri,xi,yi,zi,rmax,iflg)
               if(iflg.eq.1)then
                  npart=npart + 1
                  x(npart)=xi
                  y(npart)=yi
                  z(npart)=zi
                  h(npart)=hsep
                  call setthermo(ri,dens,druck,toru,energy,imat)
                  rho(npart)=dens
                  pmass(npart)=bfct*rho(npart)*h(npart)**3.d0 
                  pr(npart)=druck
                  temp(npart)=toru
                  u(npart)=energy
                  matter(npart)=imat
               end if
c
c--second layer
c
               xi=x0+dx
               yi=y0+dy
               zi=z0+dz
               call inside(ri,xi,yi,zi,rmax,iflg)
               if(iflg.eq.1)then
                  npart=npart + 1
                  x(npart)=xi
                  y(npart)=yi
                  z(npart)=zi
                  h(npart)=hsep
                  call setthermo(ri,dens,druck,toru,energy,imat)
                  rho(npart)=dens 
                  pmass(npart)=bfct*rho(npart)*h(npart)**3.d0
                  pr(npart)=druck
                  temp(npart)=toru
                  u(npart)=energy
                  matter(npart)=imat  
               end if
            enddo
         enddo
      enddo
c
      write(*,*)'total number of particles set: ', npart
c
      return
      end
c      
      subroutine setvelocities
c*****************************************************************
c                                                                *
c  set velocities                                                *
c                                                                *
c*****************************************************************
c
      implicit double precision (a-h,o-z) 
      parameter (idim=200000)
c
      common /part0/ npart, x(idim), y(idim), z(idim), h(idim),
     1               alpha(idim),matter(idim)
      common /part2/ pmass(idim), vx(idim), vy(idim), vz(idim)
c
c--set velocities
c
      do i=1,npart
         vx(i)=0.d0
         vy(i)=0.d0
         vz(i)=0.d0
      enddo
c
      return
      end
c
c
      subroutine setthermo(ri,dens,druck,toru,energy,imat)
c*****************************************************************
c                                                                *
c  set density,pressures,temperatures and inner energy           *
c                                                                *
c*****************************************************************
c
      implicit double precision (a-h,o-z)
      common /tab/ nlines,rt(10000),rhot(10000),presst(10000),
     1             tempt(10000),ut(10000),ut1(10000)
      common /units/ umass, udist, udens, utime, uergg, uergcc
c
      open (31,file='interpolated',status='old')
      call interp(rt,rhot,nlines,ri,dens)
      call interp(rt,presst,nlines,ri,druck)
      call interp(rt,tempt,nlines,ri,toru)
      call interp(rt,ut,nlines,ri,energy)
      write(31,*) toru,ri*udist,dens*udens,druck*uergcc,energy*uergg
c         
      if (dens .lt. 5.d0/udens) then
         imat=4
         else
         imat=5
      endif
c
      return
      end
c
      subroutine getdata
c*********************************************************
c Read the data from the corresponding file              *
c*********************************************************
      implicit double precision (a-h,o-z)
      integer nlines
c      
      common /tab/nlines,rt(10000),rhot(10000),presst(10000),
     1           tempt(10000),ut(10000),ut1(10000)
      common /units/ umass, udist, udens, utime, uergg, uergcc
      dimension rtt(10000),rhott(10000),presstt(10000),utt(10000)
c
      open(7,file='protomercury',status='old')
      do 15 i=1,500000 !nlines
         read(7,*,end=16) tempti,rti,rhoti,pressti,uti
         tempt(i)=tempti
         rtt(i)=rti
         rt(i)=rtt(i)/udist
         rhott(i)=rhoti
         rhot(i)=rhott(i)/udens
         presstt(i)=pressti
         presst(i)=presstt(i)/uergcc
         utt(i)=uti
         ut(i)=utt(i)/uergg      
 15   continue
 16   continue
      nlines=i-1
c      write(*,*) 'getdata:',nlines
c      write(*,*) (rt(i), i=1,500)
c      write(*,*) (rhot(i)*udens,i=1,468)
c      pause
c      write(*,*) (presst(i),i=1,500)
      return
      end
c                     
      subroutine interp(rt,valt,nlines,rad,val)
c*********************************************************
c                                                        *
c  linear interpolation in a table                       *
c                                                        *
c*********************************************************
c
      implicit double precision (a-h,o-z)
      double precision rt,valt
      parameter (ilines=10000)
c
      dimension rt(ilines), valt(ilines)
c
c--find location in table
c
 1    do i=1,ilines
         it=i
         if(rt(it).le. rad) go to 10
      enddo
   10 continue
c      write(*,*) it, rt(it),rad
c
c--beginning of table
c
      if(it.eq.1)then
         val=valt(1)
         return
      end if
c
c--end of table
c     
c      if (it .ge. nlines) then
c          val=valt(nlines)
c          Write(*,*) 'here 10:',val
c          return
c      endif
c
c--linear interpolation
c
      i2=it
      i1=it-1
      dr=rt(i1)-rt(i2)
      dval=valt(i1)-valt(i2)
c
c      call  LinIntPol(rad, valt(i1),rt(i1),valt(i2),rt(i2),valip)
c      if (valip .lt. 0.d0) then
c      write(*,*) 'rad,valt(i1),rt(i1),valt(i2),rt(i2),valip:'
c      write(*,*) rad,valt(i1),rt(i1),valt(i2),rt(i2),valip 
c      endif
c         
      val=valt(i1) + (rad-rt(i1))*dval/dr
c
      return
      end     
c
c      SUBROUTINE LinIntPol(x,y1,x1,y2,x2,y)
c--makes linear interpolation for f(x) (Dim n) at the point x, if 
c--(x1,f(x1)=y1) and (x2,f(x2)=y2)
c--are given (dimension y1, y2 also n of course) using Lagrange's formula
c       DOUBLE PRECISION  y1,y2,y
c       DOUBLE PRECISION x,x1,x2
c
c
c       y=((x-x2)/(x1-x2))*y1+((x-x1)/(x2-x1))*y2
c
c       END SUBROUTINE LinIntPol
c
      include 'funit.f'
      

























































