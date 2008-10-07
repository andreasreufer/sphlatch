      program setupcollision
c************************************************************
c                                                           *
c  Setup initial conditions for the simulation of a         *
c  collision between two bodies                             *
c                                                           *  
c  Parameters to specify:                                   *
c  ----------------------                                   *
c     iphys:  1 sph only                                    *
c             2 sph + solid                                 *
c     ieos:   1 tillotson                                   *
c             2 aneos                                       *
c     target: radius                                        *
c             nb. of particles                              *
c             fraction of core mass                         *
c             core:                                         * 
c                  material type                            *
c                  specific energy or temperature           *
c                  nb. of fragments                         *
c                  void fraction                            *
c             envelope:                                     *
c                  material type                            *
c                  specific energy or temperature           *
c                  nb. of fragments                         *
c                  void fraction                            *
c     projectile: radius                                    *               
c                 nb of particles                           *
c                 fraction of core mass                     *
c                 core:                                     * 
c                      material type                        *
c                      specific energy or temperature       *
c                      nb. of fragments                     *
c                      void fraction                        *
c                 envelope:                                 *
c                      material type                        *
c                      specific energy or temperature       *
c                      nb. of fragments                     *
c                      void fraction                        *
c                                                           *
c  Output                                                   *
c  ------                                                   *
c  - xdr file in ParaSPH format                             *
c  - cracks file                                            *
c  - xxx.init file                                          * 
c  - times file with suggested output times                 *
c                                                           *
c  Author: W. Benz                                          *
c  Original version: 8.12.2003                              *
c                                                           *
c************************************************************
c
      parameter (idim=2300000)
c
      common /sph   / x(idim)     , y(idim) , z(idim)    ,
     1                vx(idim)    , vy(idim), vz(idim)   ,
     2                h(idim)     , id(idim), pmass(idim),
     3                matter(idim), p(idim) , T(idim)    ,
     4                rho(idim)   , u(idim)
      common /solids/ dm(idim)    , rft(idim), sxx(idim) ,
     1                sxy(idim)   , sxz(idim), syy(idim) , 
     2                syz(idim)
      common /work  / xtt(idim), ytt(idim), ztt(idim), htt(idim),
     1                dmtt(idim)
      common /rseed/ iseed
      common /oldc / ioldcode
      common /integ/ ipind(idim), ipdam(idim), iprem(idim), 
     1               iremove(idim)
c
      character*70 title
      character*15 filename
c
c--initialize quantities
c
      pi=acos(-1.)
      vfac=1./sqrt(2.)
c
      do i=1,idim
         id(i)=0
         ipind(i)=0
         ipdam(i)=0
         iprem(i)=0
      enddo
      iseed=-765431
      call initrand2(iseed)
c
c--read run parameters
c  ioldcode=1 old SPH code
c           0 ParaSPH code
c
c  igeom=1 for spheres
c  igeom=2 for cylinders
c
      open(30,file='paramrun.in')
      read(30,100)title
  100 format(a70)
      read(30,101)filename
  101 format(a15)
      read(30,*)iphys, ieos, igeom, ioldcode
      read(30,*)rt, zt, nt, ftc, mattc, torutc, nfragtc, voidtc, 
     1          ftm, mattm, torutm, nfragtm, voidtm
      read(30,*)rp, zp, np, fpc, matpc, torupc, nfragpc, voidpc, 
     1          fpm, matpm, torupm, nfragpm, voidpm
      read(30,*)angle, vrel
      close(30)
c
      write(*,*)'---> output type: ',ioldcode
c
      if(iphys.eq.1)numvars=14
      if(iphys.eq.2)numvars=21
c
c--initialize equation of state 
c  ----------------------------
c
c  a) Tillotson eos
c
      if(ieos.eq.1)then
         call tillinit
      endif
c
c  b) ANEOS eos
c
      if(ieos.eq.2)then
         call aneosinit(mattc,mattm,matpc,matpm)
      end if
c
c--make target
c  ===========
c
      call setparticles(iphys,ieos,igeom,rt,zt,nt,ftc,mattc,torutc,ftm,
     1                  mattm,torutm,rcoret,hsept,nsett)
c
c--pre-fracture target
c  ===================
c
c  a) core
c
      if(nfragtc.gt.0)then
         rcoret2=rcoret**2
         np=0
         ntot=nsett
         do i=1,ntot
            r2=x(i)**2+y(i)**2+z(i)**2
            if(r2.lt.rcoret2)then
               np=np + 1
               ipind(np)=i
               xtt(np)=x(i)
               ytt(np)=y(i)
               ztt(np)=z(i)
               htt(np)=h(i)
               dmtt(np)=dm(i)
            endif
         enddo
         call fracture(np,nfragtc,rt,voidtc)
         call cleanup(ntot,nsett,np)
      endif
c
c  b) mantle
c
      if(nfragtm.gt.0)then
         rcoret2=rcoret**2
         np=0
         ntot=nsett
         do i=1,ntot
            r2=x(i)**2+y(i)**2+z(i)**2
            if(r2.gt.rcoret2)then
               np=np + 1
               ipind(np)=i
               xtt(np)=x(i)
               ytt(np)=y(i)
               ztt(np)=z(i)
               htt(np)=h(i)
               dmtt(np)=dm(i)
            endif
         enddo
         call fracture(np,nfragtm,rt,voidtm)
         call cleanup(ntot,nsett,np)
      endif
c
c--make projectile
c  ===============
c
      call setparticles(iphys,ieos,igeom,rp,zp,np,fpc,matpc,torupc,fpm,
     1                  matpm,torupm,rcorep,hsepp,nsetp)
c
c--pre-fracture projectile
c  ========================
c
c  a) core
c
      if(nfragtc.gt.0)then
         rcorep2=rcorep**2
         np=0
         ntot=nsett+nsetp
         do i=nsett+1,ntot
            r2=x(i)**2+y(i)**2+z(i)**2
            if(r2.lt.rcorep2)then
               np=np + 1
               ipind(np)=i
               xtt(np)=x(i)
               ytt(np)=y(i)
               ztt(np)=z(i)
               htt(np)=h(i)
               dmtt(np)=dm(i)
            endif
         enddo
         call fracture(np,nfragpc,rp,voidpc)
         call cleanup(ntot,nsetp,np)
      endif
c
c  b) mantle
c
      if(nfragpm.gt.0)then
         rcorep2=rcorep**2
         np=0
         ntot=nsett+nsetp
         do i=nsett+ 1,ntot
            r2=x(i)**2+y(i)**2+z(i)**2
            if(r2.gt.rcorep2)then
               np=np + 1
               ipind(np)=i
               xtt(np)=x(i)
               ytt(np)=y(i)
               ztt(np)=z(i)
               htt(np)=h(i)
               dmtt(np)=dm(i)
            endif
         enddo
         call fracture(np,nfragpm,rp,voidpm)
         call cleanup(ntot,nsetp,np)
      endif
c
c--assign flaws
c  ============
c
      call setflaws(igeom,nsett,rt,zt,nsetp,rp,zp)
c
c--setup collision
c  ===============
c
      call setcol(igeom,nsett,nsetp,rt,zt,rp,zp,angle,vrel,x1,y1,z1,
     1            x2,y2,z2)
c
c--damage projectile to avoid short timesteps
c  ==========================================
c
      do i=nsett+1,nsett+nsetp
         dm(i)=1.0
      enddo
c
c--compute collisional characteristics
c  ===================================
c
      xmt=0.
      xmtd=0.
      xmtc=0.
      xmtcd=0.
      tkint=0.
      volumet=0.
      rcoret2=rcoret**2
      do i=1,nsett
         xmt=xmt + pmass(i)
         r2=(x(i)-x1)**2 + (y(i)-y1)**2 + (z(i)-z1)**2
         if(r2.le.rcoret2)then
            xmtc=xmtc + pmass(i)
            if(dm(i).gt.0.)xmtcd=xmtcd + pmass(i)
         endif
         tkint=tkint + 0.5*pmass(i)*(vx(i)**2+vy(i)**2+vz(i)**2)
         if(dm(i).gt.0.)xmtd=xmtd + pmass(i)
         volumet=volumet + vfac*h(i)**3
      enddo
      rhot=xmt/volumet
      xmtc=xmtc/xmt
      xmtcd=xmtcd/xmt
      xmte=1.0 - xmtc
      xmtd=xmtd/xmt
      xmted=xmtd - xmtcd
c
      xmp=0.
      xmpc=0.
      tkinp=0.
      xmpd=0.
      xmpcd=0.
      rcorep2=rcorep**2
      volumep=0.
      ncp=0
      do i=nsett+1,nsett+nsetp
         xmp=xmp + pmass(i)
         r2=(x(i)-x2)**2 + (y(i)-y2)**2 + (z(i)-z2)**2
         if(r2.le.rcorep2)then
            xmpc=xmpc + pmass(i)
            ncp=ncp + 1
            if(dm(i).gt.0.)xmpcd=xmpcd + pmass(i)
         endif
         tkinp=tkinp + 0.5*pmass(i)*(vx(i)**2+vy(i)**2+vz(i)**2)
         if(dm(i).gt.0.)xmpd=xmpd + pmass(i)
         volumep=volumep + vfac*h(i)**3
      enddo
      rhop=xmp/volumep
      xmpc=xmpc/xmp
      xmpcd=xmpcd/xmp
      xmpe=1.0 - xmpc
      xmpd=xmpd/xmp
      xmped=xmpd - xmpcd
c
      q=tkinp/xmt
c
c--write init file
c  ===============
c
      open(30,file='impact.init')
      write(30,110)ieos
  110 format('INITIAL CONDITIONS:  (eos: ',i3,')',///)
c
      if(ieos.eq.1)then
         write(30,120)nsett, rt, xmt, rhot, xmtc, mattc, torutc, 
     1                xmtcd, xmte, mattm, torutm, xmted
  120    format('- target characteristics    : ',/,
     1          '     - nb of particles      : ',i9,/,
     1          '     - radius (cm)          : ',1pe12.4,/,
     2          '     - mass (g)             : ',1pe12.4,/,
     3          '     - mean density         : ',1pe12.4,/,
     3          '     - core fraction        : ',1pe12.4,/,
     4          '            mat. type       : ',i3,/,
     5          '            u (erg/g)       : ',1pe12.4,/,
     6          '            frac. damaged   : ',1pe12.4,/,
     6          '     - mantle fraction      : ',1pe12.4,/,
     7          '            mat. type       : ',i3,/,
     8          '            u (erg/g)       : ',1pe12.4,/,
     9          '            frac. damaged   : ',1pe12.4)
      endif
      if(ieos.eq.2)then
         write(30,121)nsett, rt, xmt, rhot, xmtc, mattc, torutc, 
     1                xmtcd, xmte, mattm, torutm, xmted
  121    format('- target characteristics    : ',/,
     1          '     - nb of particles      : ',i9,/,
     1          '     - radius (cm)          : ',1pe12.4,/,
     2          '     - mass (g)             : ',1pe12.4,/,
     3          '     - mean density         : ',1pe12.4,/,
     3          '     - core fraction        : ',1pe12.4,/,
     4          '            mat. type       : ',i3,/,
     5          '            T (K)           : ',1pe12.4,/,
     6          '            frac. damaged   : ',1pe12.4,/,
     6          '     - mantle fraction      : ',1pe12.4,/,
     7          '            mat. type       : ',i3,/,
     8          '            T (K)           : ',1pe12.4,/,
     9          '            frac. damaged   : ',1pe12.4)
      endif
c
      if(ieos.eq.1)then
         write(30,130)nsetp, rp, xmp, rhop, xmpc, matpc, torupc, 
     1                xmpcd, xmpe, matpm, torupm, xmped
  130    format('- projectile characteristics: ',/,
     1          '     - nb of particles      : ',i9,/,
     1          '     - radius (cm)          : ',1pe12.4,/,
     2          '     - mass (g)             : ',1pe12.4,/,
     3          '     - mean density         : ',1pe12.4,/,
     3          '     - core fraction        : ',1pe12.4,/,
     4          '            mat. type       : ',i3,/,
     5          '            u (erg/g)       : ',1pe12.4,/,
     6          '            frac. damaged   : ',1pe12.4,/,
     6          '     - mantle fraction      : ',1pe12.4,/,
     7          '            mat. type       : ',i3,/,
     8          '            u (erg/g)       : ',1pe12.4,/,
     9          '            frac. damaged   : ',1pe12.4)
      endif
      if(ieos.eq.2)then
         write(30,131)nsetp, rp, xmp, rhop, xmpc, matpc, torupc, 
     1                xmpcd, xmpe, matpm, torupm, xmped
  131    format('- projectile characteristics: ',/,
     1          '     - nb of particles      : ',i9,/,
     1          '     - radius (cm)          : ',1pe12.4,/,
     2          '     - mass (g)             : ',1pe12.4,/,
     3          '     - mean density         : ',1pe12.4,/,
     3          '     - core fraction        : ',1pe12.4,/,
     4          '            mat. type       : ',i3,/,
     5          '            T (K)           : ',1pe12.4,/,
     6          '            frac. damaged   : ',1pe12.4,/,
     6          '     - mantle fraction      : ',1pe12.4,/,
     7          '            mat. type       : ',i3,/,
     8          '            T (K)           : ',1pe12.4,/,
     9          '            frac. damaged   : ',1pe12.4)
      endif
c
      write(30,140)angle, vrel, q
  140 format('- collision characteristics :',/,
     1       '     - angle of incidence   : ',1pe12.4,/,
     2       '     - relative velocity    : ',1pe12.4,/,
     3       '     - Q                    : ',1pe12.4)
c
c--provide suggested times for output
c  ==================================
c
      open(25,file='times')
      ntime=10
      nrad=8
      tau=nrad*rt/abs(vrel)
      dtau=tau/ntime
      do i=1,ntime
         write(25,150)i*dtau
  150    format('output.saveTime        = ',e12.4) 
      enddo
      close(25)
c
c--set material type consistent with ParaSPH
c  =========================================
c
      if(ioldcode.eq.0)then
         do i=nsett+1,nsett+nsetp
            matter(i)=matter(i) + 256
         enddo

c
c--write xdr output file
c  =====================
c
         syz(1)=1.0
         npart=nsett+nsetp
         time=0.0
         iflag=2
         call xdrheader(iflag,filename,title,numvars,npart,time)
         do ipart=1,npart
            call xdrparticle(iflag,ipart,numvars)
         enddo
         call closexdrfile
      endif
c
c--this part for the old SPH program
c
      if(ioldcode.eq.1)then
         npart=nsett+nsetp
         write(*,*)nsett,nsetp
         time=0.
         gama=1.6
         tkin=0.
         tterm=0.
         open(55,file='sphold-001',form='unformatted')
         write(55)npart,nsett,nsetp,
     1      time,gama,(h(i),i=1,npart),tkin,
     2      tterm,(x(i),i=1,npart),(y(i),i=1,npart),(z(i),i=1,npart),
     3      (vx(i),i=1,npart),(vy(i),i=1,npart),(vz(i),i=1,npart),
     4      (u(i),i=1,npart),(T(i),i=1,npart),(pmass(i),i=1,npart),
     5      (matter(i),i=1,npart),(rho(i),i=1,npart),(sxx(i),i=1,npart),
     6      (syy(i),i=1,npart),(sxy(i),i=1,npart),(sxz(i),i=1,npart),
     7      (syz(i),i=1,npart),(dm(i),i=1,npart)
         close(55)
       endif
c
      stop
      end
      subroutine setparticles(iphys,ieos,igeom,rmax,zlen,napprox,fc,
     1                        imatc,toruc,fm,imatm,torum,rcore,hsep,
     2                        nset)
c*****************************************************************
c                                                                *
c  set particles inside a sphere of radius rmax using a close    *
c  packed setting                                                *
c                                                                *
c*****************************************************************
c
      parameter (idim=2300000)
c
      common /sph   / x(idim)     , y(idim) , z(idim)    ,
     1                vx(idim)    , vy(idim), vz(idim)   ,
     2                h(idim)     , id(idim), pmass(idim),
     3                matter(idim), p(idim) , T(idim)    ,
     4                rho(idim)   , u(idim)
c
c--compute particle spacing
c
      call setdist(igeom,rmax,zlen,napprox,hsep,xmin,xmax,
     1             ymin,ymax,zmin,zmax)
c
c--set particles on lattice
c
c  a) find where to start filling up particles
c
      istart=1
      do while (id(istart).ne.0)
         istart=istart + 1
      enddo
      call setonlattice(istart,igeom,rmax,zlen,hsep,xmin,xmax,ymin,ymax,
     1                  zmin,zmax,nset)
c
c--set velocities
c
      call setvelocities(istart,nset)
c
c--only for solid bodies
c
      if(iphys.eq.2)then
c
c  a) set deviatoric stress tensor
c
         call setdeviatoric(istart,nset)
c
c--set damage
c
         call setdamage(istart,nset)
c
      endif
c
c--find equilibrium densities corresponding to u/T
c
      call findrho(ieos,imatc,toruc,rhoc,uc,pc)
      write(*,*)'core density found: ',imatc,rhoc
      call findrho(ieos,imatm,torum,rhom,um,pm)
      write(*,*)'mantle density found: ',imatm,rhom
c
c--compute internal structure related quantities
c
      if(fc.gt.0.)then
c        frac=(fm/fc)*(rhoc/rhom) + 1.0
c        rcore=rmax/frac**0.333333
         frac=(1./fc -1.0)*(rhoc/rhom) + 1.0
         rcore=rmax/(frac**0.3333333)
      else
         rcore=0.0
      endif
      write(*,*)'core radius found: ',rcore
c
c--set density and masses
c
      call setrhom(istart,nset,rcore,imatc,rhoc,imatm,rhom)
c
c--set thermodynamics
c
      call setthermo(ieos,istart,nset,rcore,toruc,uc,pc,torum,um,pm)
c
      return
      end
      subroutine setdist(igeom,rmax,zlen,napprox,hsep,xmin,xmax,
     1                   ymin,ymax,zmin,zmax)
c*******************************************************************
c                                                                  *
c  computes the spacing between the particles                      *
c                                                                  *
c*******************************************************************
c
      pi=acos(-1.)
      third=1./3.
c
      if(igeom.eq.1)then
         volume=4.*third*pi*rmax**3
      elseif(igeom.eq.2)then
         volume=pi*rmax*rmax*zlen
      endif
      hsep=(volume/napprox)**third
      xmin=-rmax
      xmax=rmax
      ymin=-rmax
      ymax=rmax
      if(igeom.eq.1)then
         zmin=-rmax
         zmax=rmax
      elseif(igeom.eq.2)then
         zmin=-0.5*zlen
         zmax=0.5*zlen
      endif
c
      return
      end
      subroutine inside(igeom,xi,yi,zi,rmax,zlen,iflg)
c************************************************************
c                                                           *
c  iflg=0 if particle is outside object                     *
c  iflg=1 if particle is inside  object                     *
c                                                           *
c************************************************************
c
      iflg=0
c
      if(igeom.eq.1)then
         r2=xi*xi + yi*yi + zi*zi
         ratio=r2/(rmax*rmax)
         if(ratio.le.1.)then
            iflg=1
         endif
      elseif(igeom.eq.2)then
         rcyl2=xi*xi + yi*yi
         ratio=rcyl2/(rmax*rmax)
         if(ratio.le.1)then
            iflg=1
         endif
      endif
c
      return
      end
      subroutine setonlattice(istart,igeom,rmax,zlen,hsep,xmin,xmax,
     1                        ymin,ymax,zmin,zmax,nset)
c****************************************************************
c                                                               *
c  set particles on a close-packed lattice                      *
c                                                               *
c****************************************************************
c
      parameter (idim=2300000)
c
      common /sph   / x(idim)     , y(idim) , z(idim)    ,
     1                vx(idim)    , vy(idim), vz(idim)   ,
     2                h(idim)     , id(idim), pmass(idim),
     3                matter(idim), p(idim) , T(idim)    ,
     4                rho(idim)   , u(idim)
c
c--set particles on lattice
c
      dx=0.5*hsep
      dy=hsep/sqrt(12.)
      dz=hsep*sqrt(2./3.)
c
      xstep=0.5*hsep
      ystep=xstep*sqrt(3.)
      zstep=2.*dz
c
c--first layer
c
      npart=istart - 1
      nset=0
      factor=0.5
      do y0=ymin,ymax,ystep
         xmin=xmin + factor*hsep
         factor=-factor
         do x0=xmin,xmax,hsep
            do z0=zmin,zmax,zstep
               xi=x0
               yi=y0
               zi=z0
               call inside(igeom,xi,yi,zi,rmax,zlen,iflg)
               if(iflg.eq.1)then
                  npart=npart + 1
                  nset=nset + 1
                  x(npart)=xi
                  y(npart)=yi
                  z(npart)=zi
                  id(npart)=npart
                  h(npart)=hsep
               end if
c
c--second layer
c
               xi=x0+dx
               yi=y0+dy
               zi=z0+dz
               call inside(igeom,xi,yi,zi,rmax,zlen,iflg)
               if(iflg.eq.1)then
                  npart=npart + 1
                  nset=nset + 1
                  x(npart)=xi
                  y(npart)=yi
                  z(npart)=zi
                  id(npart)=npart
                  h(npart)=hsep
               end if
            enddo
         enddo
      enddo
c
c--reset center of mass at origin
c
      cmx=0.
      cmy=0.
      cmz=0.
      do i=istart,istart + nset - 1 
         cmx=cmx + x(i)
         cmy=cmy + y(i)
         cmz=cmz + z(i)
      enddo
      cmx=cmx/nset
      cmy=cmy/nset
      cmz=cmz/nset
      do i=istart,istart + nset - 1 
         x(i)=x(i) - cmx
         y(i)=y(i) - cmy
         z(i)=z(i) - cmz
      enddo
c
      write(*,*)'total number of particles set: ', npart - istart + 1
c
      return
      end
      subroutine setvelocities(istart,nset)
c*****************************************************************
c                                                                *
c  set velocities                                                *
c                                                                *
c*****************************************************************
c
      parameter (idim=2300000)
c
      common /sph   / x(idim)     , y(idim) , z(idim)    ,
     1                vx(idim)    , vy(idim), vz(idim)   ,
     2                h(idim)     , id(idim), pmass(idim),
     3                matter(idim), p(idim) , T(idim)    ,
     4                rho(idim)   , u(idim)
c
c--set velocities
c
      do i=istart,istart + nset - 1
         vx(i)=0.
         vy(i)=0.
         vz(i)=0.
      enddo
c
      return
      end
      subroutine setdeviatoric(istart,nset)
c*****************************************************************
c                                                                *
c  set deviatoric stress tensor                                  *
c                                                                *
c*****************************************************************
c
      parameter (idim=2300000)
c
      common /solids/ dm(idim)    , rft(idim), sxx(idim) ,
     1                sxy(idim)   , sxz(idim), syy(idim) , 
     2                syz(idim)
c
c--set deviatoric stresss tensor
c
      do i=istart,istart + nset - 1
         sxx(i)=0.
         syy(i)=0.
         sxy(i)=0.
         sxz(i)=0.
         syz(i)=0.
      enddo
c
      return
      end
      subroutine setdamage(istart,nset)
c*****************************************************************
c                                                                *
c  set damage                                                    *
c                                                                *
c*****************************************************************
c
      parameter (idim=2300000)
c
      common /solids/ dm(idim)    , rft(idim), sxx(idim) ,
     1                sxy(idim)   , sxz(idim), syy(idim) , 
     2                syz(idim)
c
c
c--set deviatoric stresss tensor
c
      do i=istart,istart + nset - 1
         dm(i)=0.
         rft(i)=0.
      enddo
c
      return
      end
      subroutine setrhom(istart,nset,rcore,imatc,rhoc,imatm,rhom)
c*****************************************************************
c                                                                *
c  set matter, density and masses                                *
c                                                                *
c*****************************************************************
c
      parameter (idim=2300000)
c
      common /sph   / x(idim)     , y(idim) , z(idim)    ,
     1                vx(idim)    , vy(idim), vz(idim)   ,
     2                h(idim)     , id(idim), pmass(idim),
     3                matter(idim), p(idim) , T(idim)    ,
     4                rho(idim)   , u(idim)
c
c--set density and masses
c
      rcore2=rcore*rcore
      factor=1./sqrt(2.)
      nmatc=0
      nmate=0
      do i=istart,istart + nset - 1
         r2=x(i)**2 + y(i)**2 + z(i)**2
         if(r2.le.rcore2)then
            matter(i)=imatc
            rho(i)=rhoc
            nmatc=nmatc + 1
         else
            matter(i)=imatm
            rho(i)=rhom
            nmate=nmate + 1
         endif
         pmass(i)=factor*rho(i)*h(i)**3
      enddo
      write(*,*)'nb. of core particles:  ', nmatc
      write(*,*)'nb. of mantleparticles: ', nmate
c
      return
      end
      subroutine setthermo(ieos,istart,nset,rcore,toruc,uc,pc,torum,
     1                     um,pm)
c*****************************************************************
c                                                                *
c  set matter, density and masses                                *
c                                                                *
c*****************************************************************
c
      parameter (idim=2300000)
c
      common /sph   / x(idim)     , y(idim) , z(idim)    ,
     1                vx(idim)    , vy(idim), vz(idim)   ,
     2                h(idim)     , id(idim), pmass(idim),
     3                matter(idim), p(idim) , T(idim)    ,
     4                rho(idim)   , u(idim)
c
c--set internal energy and temperature
c
      rcore2=rcore*rcore
c
c  a) Tillotson eos
c
      if(ieos.eq.1) then
         do i=istart,istart + nset - 1
            r2=x(i)**2+y(i)**2+z(i)**2 
            if(r2.le.rcore2)then
               u(i)=uc
               T(i)=0.0
               p(i)=pc
            else
               u(i)=um
               T(i)=0.0
               p(i)=pm
            endif
         enddo
         return
      endif
c
c  b) Aneos eos
c
      if(ieos.eq.2) then
         do i=istart,istart + nset - 1
            r2=x(i)**2+y(i)**2+z(i)**2 
            if(r2.le.rcore2)then
               u(i)=uc
               T(i)=toruc
               p(i)=pc
            else
               u(i)=um
               T(i)=torum
               p(i)=pm
            endif
         enddo
         return
      endif
c
      return
      end
      subroutine findrho(ieos,imat,toru,rho0,u0,p0)
c**************************************************************
c                                                             *
c  find the equilibrium density corresponding to a given u/T  *
c  assuming a cold state                                      *
c                                                             *
c**************************************************************
c
      data itermax/50/
      data tinyp/1.e5/
c
c--iterate on density until minimum pressure is achieved
c  -----------------------------------------------------
c
c  a) set bounds
c
      rhomax=10.0
      if(ieos.eq.1)then
         u0=toru
         call eost(imat,u0,rhomax,csi,pmax)
      endif
      if(ieos.eq.2)then
         temp0=toru
         call eosa(imat,temp0,rhomax,pmax,ui,csi,kpai)
      endif
      if(pmax.lt.0.)then
         write(*,*)'upper limit not correct!'
         stop
      endif
c
      rhomin=1.0
      if(ieos.eq.1)then
         u0=toru
         call eost(imat,u0,rhomin,csi,pmin)
      endif
      if(ieos.eq.2)then
   10    continue
         temp0=toru
         call eosa(imat,temp0,rhomin,pmin,ui,csi,kpai)
         if(abs(pmin).lt.1.e-5)then
            rhomin=1.1*rhomin
            go to 10
         endif
      endif
      if(pmin.gt.0.)then
         write(*,*)'lower limit not correct!'
         stop
      endif
c
c  b) iterate using midpoint method
c
      do i=1,itermax
         rho=0.5*(rhomax+rhomin)
         if(ieos.eq.1)then
            u=toru
            call eost(imat,u,rho,cs,p)
         endif
         if(ieos.eq.2)then
            temp0=toru
            call eosa(imat,temp0,rho,p,u,cs,kpa)
         endif
c
         if(abs(p).lt.tinyp)then
            rho0=rho
            u0=u
            p0=p
            return
         endif
c
         if(p.gt.0.)then
            rhomax=rho
         else
            rhomin=rho
         endif
      enddo
c
      write(*,*)'max iteration steps exceeded!'
      stop
      end
      subroutine setcol(igeom,n1,n2,radt,zt,radp,zp,angle,vrel,x1,y1,z1,
     1                  x2,y2,z2)
c********************************************************************
c                                                                   *
c  set up initial conditions for collisions                         *
c                                                                   *
c********************************************************************
c
      parameter (idim=2300000)
c
      common /sph   / x(idim)     , y(idim) , z(idim)    ,
     1                vx(idim)    , vy(idim), vz(idim)   ,
     2                h(idim)     , id(idim), pmass(idim),
     3                matter(idim), p(idim) , T(idim)    ,
     4                rho(idim)   , u(idim)
c
      ntot=n1+n2
      pi=acos(-1.0)
c
c--compute cm positions of bodies
c  (in the frame of the CM of the target)
c
      if(igeom.eq.1)then
         dtot=radt + radp + 2.005*(h(1)+h(n1+1))
         offset=radt*sin(angle)
         dist=sqrt(dtot**2 - offset**2)
      elseif(igeom.eq.2)then
         dtot=0.5*zt + 0.5*zp + 2.005*(h(1)+h(n1+1)) 
         offset=0.0
         dist=dtot
         if(angle.ne.0.)then
            write(*,*)'setup not yet defined for this angle',angle
            stop
         endif
      endif
c
      x1=0.0
      y1=0.0
      z1=0.0
      vz1=0.0
      vy1=0.0
      vx1=0.0
c
      x2=offset
      y2=0.0
      z2=dist
      vx2=0.
      vy2=0.
      vz2=vrel
c
c--set bodies
c
      do i=1,n1
         x(i)=x(i) + x1
         y(i)=y(i) + y1
         z(i)=z(i) + z1
         vx(i)=vx(i) + vx1
         vy(i)=vy(i) + vy1
         vz(i)=vz(i) + vz1
      enddo
c
      do i=n1+1,ntot
         x(i)=x(i) + x2
         y(i)=y(i) + y2
         z(i)=z(i) + z2
         vx(i)=vx(i) + vx2
         vy(i)=vy(i) + vy2
         vz(i)=vz(i) + vz2
      enddo
c
      return
      end
      subroutine setflaws(igeom,n1,rt,zt,n2,rp,zp)
c**************************************************************
c                                                             *
c  set up a crack distribution following a Weibull            *
c  distribution.                                              *
c                                                             *
c**************************************************************
c
      parameter (nmatmax=21)
c
      double precision cweib
      common /fracm/ pweib(nmatmax), cweib(nmatmax)
      common /oldc / ioldcode
c
c--set material dependent constants
c
      npart=n1+n2
      call fractmodel(npart)
c
c--for fold ortran programm
c
      if(ioldcode.eq.1)then
         open(41,file='cracks')
         do i=1,nmatmax
            write(41,*)i,pweib(i),cweib(i)
         enddo
      endif
c
c--first object
c  ------------
c
      iallow=1
      ishield=0
      i1=1
      i2=n1
      write(*,*)'setting flaws in body 1: ',i1,i2
      call setf(igeom,iallow,ishield,i1,i2,rt,zt)
c
c--second object
c  -------------
c
      iallow=0
      ishield=0
      i1=n1+1
      i2=n1+n2
      write(*,*)'setting flaws in body 2: ',i1,i2
      call setf(igeom,iallow,ishield,i1,i2,rp,zp)
c
c     close(41)
c
      return
      end
      subroutine setf(igeom,iallow,ishield,i1,i2,rad,zlen)
c**********************************************************
c                                                         *
c  set flaws in objects                                   *
c                                                         *
c**********************************************************
c
      parameter (idim=2300000)
      parameter (nmatmax=21)
c
      common /sph   / x(idim)     , y(idim) , z(idim)    ,
     1                vx(idim)    , vy(idim), vz(idim)   ,
     2                h(idim)     , id(idim), pmass(idim),
     3                matter(idim), p(idim) , T(idim)    ,
     4                rho(idim)   , u(idim)
      common /fracm/ pweib(nmatmax), cweib(nmatmax)
      common /fract/ epsmin(idim), epsmax(idim), 
     1               grav(idim), nflaws(idim), xm(idim), 
     2               acoef(idim), young(idim)
      common /till / rozero(nmatmax), smalla(nmatmax), 
     1               uzero(nmatmax),  smallb(nmatmax), 
     2               capa(nmatmax), es(nmatmax), 
     3               esp(nmatmax), alpha(nmatmax), 
     4               beta(nmatmax), capb(nmatmax)
      common /sorte/ indxe(idim)
      common /rseed/ iseed
      common /oldc / ioldcode
c
      double precision const, cweib, pw, vol, epsmeantrue,
     1                 xflawmax, xflaw, randm, rand2
      dimension ebin(200), nebin(200)
c
c--initialize
c
      onesixth=1./6.
      pi=acos(-1.)
      nused=i2 - i1 + 1 
      do i=i1,i2
         nflaws(i)=0
         epsmin(i)=1e30
         epsmax(i)=0.
         grav(i)=0.0
      enddo
c
c--compute volume
c
      if(igeom.eq.1)then
         vol=4.d0*pi*rad*rad*rad/3.
      elseif(igeom.eq.2)then
         vol=pi*rad*rad*zlen
      endif
c
c--compute gravity shielding if necessary
c
      if(ishield.eq.1)then 
         rad2=rad**2
         do i=i1,i2
            rho0=rozero(matter(i))
            ri=sqrt(x(i)**2+y(i)**2+z(i)**2)
            pressg=1.39696e-7*rho0**2*(rad2-ri**2)
            grav(i)=max(pressg,0.)
         enddo
      end if
c
c--no flaws
c
      if(iallow.eq.0)then
         ei=1.e20
         xi=1.0
         ni=0
         ai=0.
         if(ioldcode.eq.0)then
            do i=i1,i2
               call writecracks(ei,xi,ni,ai,young(i),grav(i))
            enddo
         endif
         if(ioldcode.eq.1)then
            do i=i1,i2
               write(41)ei,xi,ni,ai,young(i), grav(i)
            enddo
         endif
         return
      end if
c
c--set flaws
c
      nflawset=15000000
      xfstmax=300
      epsmeantrue=0.0
      do iflaw=1,nflawset
c
c--pick particle
c
         randm=rand2(iseed)
         ip=int(randm*nused) + i1
         imat=matter(ip)
         pw=1.d0/pweib(imat)
         const=1.d0/((dble(cweib(imat))*vol)**pw)
c
c--set flaw
c
         randm=rand2(iseed)
         xflaw=randm*xfstmax
         eps=const*xflaw**pw
         epsmeantrue=epsmeantrue + eps
         epsmin(ip)=min(epsmin(ip),eps)
         epsmax(ip)=max(epsmax(ip),eps)
         nflaws(ip)=nflaws(ip) + 1
      enddo
      epsmeantrue=epsmeantrue/nflawset
c
c--check for holes
c
      do ip=i1,i2
         if(nflaws(ip).eq.0)then
            imat=matter(ip)
            pw=1.d0/pweib(imat)
            const=1.d0/((dble(cweib(imat))*vol)**pw)
            randm=rand2(iseed)
            xflaw=randm*xfstmax
            eps=const*xflaw**pw
            epsmin(ip)=min(epsmin(ip),eps)
            epsmax(ip)=max(epsmax(ip),eps)
            nflaws(ip)=nflaws(ip) + 1
         endif
      enddo
c
c--compute average epsmin and normalize
c
      call indexx(nused,epsmin,indxe)
      epsminbar=0.
      nave=min(nused,1000)
      nave=nused
      do i=nused,nused-nave+1,-1
         ip=indxe(i)
         epsminbar=epsminbar + epsmin(ip)
      enddo
      epsminbar=epsminbar/nave
      ratio=epsmeantrue/epsminbar
c      do ip=i1,i2
c         epsmin(ip)=ratio*epsmin(ip)
c         epsmax(ip)=ratio*epsmax(ip)
c      enddo
c
c--compute local Weibull exponent
c
      do ip=i1,i2
         if(nflaws(ip).eq.1) then
            xm(ip)=1.
         else 
            xm(ip)=log(float(nflaws(ip)))/log(epsmax(ip)/epsmin(ip))
         endif
      enddo
c
c--compute flaw statistics
c
      epsminbar=0.
      epsminmin=1e30
      epsminmax=0.
      nfmin=300000000
      nfmax=0
      nftot=0
      do ip=i1,i2
         epsminbar=epsminbar + epsmin(ip)
         epsminmin=min(epsminmin,epsmin(ip))
         epsminmax=max(epsminmax,epsmin(ip))
         nfmin=min(nfmin,nflaws(ip))
         nfmax=max(nfmax,nflaws(ip))
         nftot=nftot + nflaws(ip)
      enddo
      epsminbar=epsminbar/nused
c
c--build distribution
c
      nbin=100
      do i=1,nbin
         ebin(i)=0.
         nebin(i)=0
      enddo
      de=(1.1*epsminmax-0.9*epsminmin)/nbin
      do i=1,nbin
         emin=0.9*epsminmin + (i-1)*de
         emax=emin + de
         do ip=i1,i2
            if(epsmin(ip).ge.emin.and.epsmin(ip).lt.emax)then
               nebin(i)=nebin(i) + 1
               ebin(i)=0.5*(emin+emax)
            endif
         enddo
      enddo
      nsum=0.
      do i=1,nbin
         if(nebin(i).ne.0)then
            xbin=float(nebin(i))/float(nused)
            write(97,*)ebin(i),xbin
         endif
         nsum=nsum + nebin(i) 
      enddo
      write(*,*)'total found: ',nsum
c
c--write flaw file
c
c  a) new code
c
      if(ioldcode.eq.0)then
         xmbar=0.
         do i=i1,i2
            call writecracks(epsmin(i),xm(i),nflaws(i),acoef(i),
     1                       young(i),grav(i))
            xmbar=xmbar + xm(i)
         enddo
         xmbar=xmbar/(i2-i1+1)
      endif
c
c  b) old code
c
      if(ioldcode.eq.1)then
         xmbar=0.
         do i=i1,i2
            write(41,100)epsmin(i), xm(i), nflaws(i),
     1                   acoef(i), young(i), grav(i)
  100       format(3(1pe13.6,1x),i6,2x,3(1pe13.6,1x))
            xmbar=xmbar + xm(i)
         enddo
         xmbar=xmbar/(i2-i1+1)
      endif
c
c--write statistics
c
      write(*,*)'tot. nb. of flaws set in body    : ', nflawset,nftot
      write(*,*)'min and max nb. of flaws in part.: ', nfmin,nfmax
      write(*,*)'average strain (theory, actual)  : ', epsmeantrue,
     1                                                 epsminbar 
      write(*,*)'min and max activation strains   : ', epsminmin,
     1                                                 epsminmax
      write(*,*)'avergae Weibull m                : ', xmbar
c
      return
      end
      subroutine fractmodel(npart)
c************************************************
c                                               *
c  subroutine to initialize fracture model      *
c                                               *
c************************************************
c
      parameter (idim=2300000)
      parameter (nmatmax=21)
c
      common /sph   / x(idim)     , y(idim) , z(idim)    ,
     1                vx(idim)    , vy(idim), vz(idim)   ,
     2                h(idim)     , id(idim), pmass(idim),
     3                matter(idim), p(idim) , T(idim)    ,
     4                rho(idim)   , u(idim)
      common /fract/ epsmin(idim), epsmax(idim), 
     1               grav(idim), nflaws(idim), xm(idim), 
     2               acoef(idim), young(idim)
      common /shear/ xmu(nmatmax)
      common /till / rozero(nmatmax), smalla(nmatmax), 
     1               uzero(nmatmax),  smallb(nmatmax), 
     2               capa(nmatmax), es(nmatmax), 
     3               esp(nmatmax), alpha(nmatmax), 
     4               beta(nmatmax), capb(nmatmax)
c
c--read material data file
c
      call matdata 
c
c--obtain rest of material quantities from tillotson eos
c
      call tillinit
c
c--compute single-flaw slope of damage integral and young's modulus
c
      do i=1,npart
         imat=matter(i)
         bulk=capa(imat)
         xmui=xmu(imat)
c
c--longitudital sound speed (Rayleigh sound speed)
c
         cg=0.4d0*sqrt((bulk+1.33333333d0*xmui)/rozero(imat))
c
c--explicit growth
c
         acoef(i)=cg/(2.d0*h(i))
c
c--young's modulus
c
         yg=9.d0*bulk/(3.d0*bulk+xmui)
         young(i)=xmui*yg
      enddo
c
      return
      end
      subroutine matdata
c*************************************************
c                                                *
c  subroutine reading material data for elastic, * 
c  plastic and fracture model                    *
c                                                *
c  --> read has been changed to accomodate C     *
c      arrays in ParaSPH                         *
c                                                *
c*************************************************
c
      parameter (nmatmax=21)
c
      double precision cweib
c
      common /shear/ xmu(nmatmax)
      common /yield/ umelt(nmatmax), yie(nmatmax)
      common /fracm/ pweib(nmatmax), cweib(nmatmax)
c
      character*70 filen
      character*50 dumy
      character*2 path
c
      data path /'./'/
c
c--open data file
c
      filen=path // 'matdata.input'
      open(30,file=filen)
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
c  i) shear modulus xmu
c
      read(30,100)dumy
      read(30,100)dumy
      read(30,100)dumy
      do i=1,nline
         imatmin=(i-1)*5 + 1 
         imatmax=min(nmat,imatmin+4)
         read(30,*)(xmu(j),j=imatmin,imatmax)
      enddo
c
c  ii) melt energy umelt
c
      read(30,100)dumy
      read(30,100)dumy
      read(30,100)dumy
      do i=1,nline
         imatmin=(i-1)*5 + 1 
         imatmax=min(nmat,imatmin+4)
         read(30,*)(umelt(j),j=imatmin,imatmax)
      enddo
c
c  iii) yielding yie
c
      read(30,100)dumy
      read(30,100)dumy
      read(30,100)dumy
      do i=1,nline
         imatmin=(i-1)*5 + 1 
         imatmax=min(nmat,imatmin+4)
         read(30,*)(yie(j),j=imatmin,imatmax)
      enddo
c
c  iv) fracture model pweib
c
      read(30,100)dumy
      read(30,100)dumy
      read(30,100)dumy
      do i=1,nline
         imatmin=(i-1)*5 + 1 
         imatmax=min(nmat,imatmin+4)
         read(30,*)(pweib(j),j=imatmin,imatmax)
      enddo
c
c  v) fracture model cweib
c
      read(30,100)dumy
      read(30,100)dumy
      read(30,100)dumy
      do i=1,nline
         imatmin=(i-1)*5 + 1 
         imatmax=min(nmat,imatmin+4)
         read(30,*)(cweib(j),j=imatmin,imatmax)
      enddo
      close(30)
c
      return
      end
      subroutine fracture(np,nseed,rad,voidf)
c**************************************************************
c                                                             *
c  pre-fracture a target                                      *
c                                                             *
c**************************************************************
c
      parameter (idim=2300000)
      parameter (nbseed=1)
      parameter (ndam=2)
c
      common /work  / xtt(idim), ytt(idim), ztt(idim), htt(idim),
     1                dmtt(idim)
      common /rseed/ iseed
      common /list / nlistpart(nbseed), listpart(nbseed,idim)
      common /integ/ iflag(idim), ipdam(idim), iprem(idim), 
     1               listdam(idim)
c
      dimension neilist(idim), listseed(nbseed) 
c
      double precision rand2, random 
c
c--set parameters
c
      expdist=3.0
c
c--initialize flag
c
      do i=1,idim
         iflag(i)=0
         do j=1,nseed
            listpart(j,i)=0
         enddo
      enddo
      do i=1,nbseed
         nlistpart(i)=0
      enddo
c
c--set seed according to suitable probability distribution
c
      nused=0
      do i=1,nseed
    5    call randombody(np,ipart)
         if(iflag(ipart).ne.0)goto 5
         r=sqrt(xtt(ipart)**2+ytt(ipart)**2+ztt(ipart)**2)
         frac=(r/rad)**expdist
         random=rand2(iseed)
         if(frac.lt.random)goto 5
         nused=nused + 1
         listseed(nused)=ipart
         iflag(ipart)=nused
         nlistpart(nused)=nlistpart(nused) + 1
         listpart(nused,nlistpart(nused))=ipart
      enddo
      nused=nseed
c
c--build linked list
c
      call mllist (np)
c
c--pick a fragment at random
c
      nstep=0
      do while (nused.lt.np)
         nstep=nstep + 1
         call randombody(nseed,ifrag)
         icent=listseed(ifrag)
         if(mod(nstep,10000).eq.0)then
            write(*,*)'nseed, ifrag, icent',nseed,ifrag,icent
            write(*,*)'seed done: ',nstep
            write(*,*)'particles done: ',nused
            write(*,*)'particles to go: ',np-nused
         endif
c
c  a) find all neighbors and attribute them to fragment 
c
         call findneigh(icent,nlist,neilist)
         do j=1,nlist
            jpart=neilist(j)
            if(iflag(jpart).eq.0)then
               nused=nused + 1
               iflag(jpart)=ifrag
               nlistpart(ifrag)=nlistpart(ifrag) + 1
               listpart(ifrag,nlistpart(ifrag))=jpart
            endif
         enddo
c
c  b) pick start particle for fragment
c
         nfi=nlistpart(ifrag)
         call randombody(nfi,ifn)
         listseed(ifrag)=listpart(ifrag,ifn)
      enddo
c
c--set damage
c  ----------
c
      do ipart=1,np
         ifrag=iflag(ipart)
         call findneigh(ipart,nlist,neilist)
         notsame=0
         do j=1,nlist
            jpart=neilist(j)
            if(iflag(jpart).ne.ifrag)then
               notsame=notsame + 1
            endif
         enddo
         if(notsame.gt.ndam)then
            dmtt(ipart)=1.0
         endif
      enddo
      do i=1,np
         if(dmtt(i).eq.1.)then
            iflag(i)=-abs(iflag(i))
            ipdam(i)=1
         endif
      enddo
c
c--set voids
c
      nd=0 
      do i=1,np
         if(nseed.gt.1)then
            if(dmtt(i).eq.1.)then
               nd=nd + 1
               listdam(nd)=i
            endif
         else
            nd=nd + 1
            listdam(nd)=i
         endif
      enddo
      nrem=voidf*np
      if(nrem.gt.nd)then
         write(*,*)'void fraction unreacheable!'
         stop
      endif
c
      do i=1,nrem
   30    call randombody(nd,ir)
         jpart=listdam(ir)
         if(iprem(jpart).eq.0)then
            iprem(jpart)=1
         else
            goto 30
         endif
      enddo
c
      return
      end
      subroutine randombody(nb,i)
c***************************************************
c                                                  *
c returns a particle number chosen at random       *
c between 1 and nb                                 *
c                                                  *
c***************************************************
c
      common/rseed / iseed
c
      double precision random, rand2
c
      random=rand2(iseed)
c
      xnum=nb*random
      i=int(xnum) + 1
      i=min(nb,i)
c
      return
      end
      subroutine findneigh(ipart,nlist,neilist)
c************************************************************
c                                                           *
c  This subroutine provides a list of particles within a    *
c  radius 2h of a particle.                                 *
c  This subroutine uses the linked list algorithm.          *
c                                                           *
c************************************************************
c
      parameter (idim=2300000)
      parameter (ncell=55)
c
      common /sph   / x(idim)     , y(idim) , z(idim)    ,
     1                vx(idim)    , vy(idim), vz(idim)   ,
     2                h(idim)     , id(idim), pmass(idim),
     3                matter(idim), p(idim) , T(idim)    ,
     4                rho(idim)   , u(idim)
      common /llisi/ ncellused, ll(idim),
     1               icell(0:ncell,0:ncell,0:ncell)
      common /llisr/ xcell(0:ncell),ycell(0:ncell),zcell(0:ncell),
     2               xmapl(0:ncell), xmapu(0:ncell), xlow, xup,
     3               ymapl(0:ncell), ymapu(0:ncell), ylow, yup,
     4               zmapl(0:ncell), zmapu(0:ncell), zlow, zup
      common /sortc/ irankx(idim), iranky(idim), irankz(idim),
     1               nslice, nx, ny, nz
c
      dimension neilist(idim)
c
      nlist=0
      xi=x(ipart)
      yi=y(ipart)
      zi=z(ipart)
      rcut=1.05*h(ipart)
      rcut2=rcut*rcut
c
c--find cells to serch for neighbors
c  ---------------------------------
c
c  1) X axis
c  
      ratioi=irankx(ipart)/nslice
      ic=nint(ratioi)
c
c  a) minimum x index
c
      xmin=max(xi-rcut,xlow)
      icxmin=ic
    1 continue
      if(xmin.lt.xmapl(icxmin))then
         icxmin=icxmin - 1
         go to 1
      end if
c
c  b) maximum x index
c
      xmax=min(xi+rcut,xup)
      icxmax=ic
    2 continue
      if(xmax.gt.xmapu(icxmax))then
         icxmax=icxmax + 1
         go to 2
      end if
c
c--2) Y axis
c
      ratioj=iranky(ipart)/nslice
      jc=nint(ratioj)
c
c  a) minimum y index
c
      ymin=max(yi-rcut,ylow)
      icymin=jc
    3 continue
      if(ymin.lt.ymapl(icymin))then
         icymin=icymin - 1
         go to 3
      end if
c
c  b) maximum y index
c
      ymax=min(yi+rcut,yup)
      icymax=jc
    4 continue
      if(ymax.gt.ymapu(icymax))then
         icymax=icymax + 1
         go to 4
      end if
c
c--Z axis
c
      ratiok=irankz(ipart)/nslice
      kc=nint(ratiok)
c
c  e) minimum z index
c
      zmin=max(zi-rcut,zlow)
      iczmin=kc
    5 continue
      if(zmin.lt.zmapl(iczmin))then
         iczmin=iczmin - 1
         go to 5
      end if
c
c  f) maximum z index
c
      zmax=min(zi+rcut,zup)
      iczmax=kc
    6 continue
      if(zmax.gt.zmapu(iczmax))then
         iczmax=iczmax + 1
         go to 6
      end if
c
c--find potential neighbors
c  ------------------------
c
      do ic=icxmin,icxmax
         do jc=icymin,icymax
            do kc=iczmin,iczmax
               j=icell(ic,jc,kc)
               do while (j.ne.0)
                  dx=xi-x(j)
                  dy=yi-y(j)
                  dz=zi-z(j)
                  d2=dx*dx + dy*dy + dz*dz
                  if(d2.le.rcut2)then
                     if(j.ne.ipart) then
                        nlist=nlist + 1
                        neilist(nlist)=j
                     endif
                  endif
                  j=ll(j)
               enddo
            enddo
         enddo
      enddo

c
      return
      end
      double precision function rand2(idum)
c*******************************************************c
c                                                       c
c     random number generator (w. press)                c
c                                                       c
c*******************************************************c
c
      implicit double precision (a-h,o-z)
c
      parameter (ntab=32)
c
      common /jrand/ im1, im2, imm1, ia1, ia2, iq1, iq2, ir1,
     1               ir2, ndiv, iv(ntab), iy, idum2
      common /rrand/ am, eps, rnmx
c
c--compute random number
c  ---------------------
c
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if(idum.lt.0) idum=idum + im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if(idum2.lt.0) idum2=idum2 + im2
      j=1 + iy/ndiv
      iy=iv(j) - idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy + imm1
      rand2=min(am*iy,rnmx)
c
      return
      end
      subroutine initrand2(idum)
c*******************************************************c
c                                                       c
c     initialize random number generator                c
c                                                       c
c*******************************************************c
c
      implicit double precision (a-h,o-z)
c
      parameter (ntab=32)
c
      common /jrand/ im1, im2, imm1, ia1, ia2, iq1, iq2, ir1,
     1               ir2, ndiv, iv(ntab), iy, idum2
      common /rrand/ am, eps, rnmx
c
c--set generator constants
c  -----------------------
c
      im1=2147483563
      im2=2147483399
      am=1.0/im1
      imm1=im1-1
      ia1=40014
      ia2=40692
      iq1=53668
      iq2=52774
      ir1=12211
      ir2=3791
      ndiv=1+imm1/ntab
      eps=1.2e-14
      rnmx=1.0e0-eps
      idum2=123456789
      iy=0
      do i=1,ntab
         iv(i)=0
      enddo
c
c--initialize quantities
c  ---------------------
c
      idum=max(-idum,1)
      idum2=idum
      do j=ntab+8,1,-1
         k=idum/iq1
         idum=ia1*(idum-k*iq1)-k*ir1
         if(idum.lt.0) idum=idum+im1
         if(j.le.ntab) iv(j)=idum
      enddo
      iy=iv(1)
c
      return
      end
      subroutine mllist (npart)
c************************************************************
c                                                           *
c  this subroutine builds the linked list by superposing a  *
c  cubic lattice on the computational box. The linked list  *
c  is obtained using the rank of the particles sorted       *
c  along increasing coordinates.                            *
c                                                           *
c************************************************************
c
      parameter (idim=2300000)
      parameter (nocc=10)
      parameter (ncell=55)
c
      common /llisi/ ncellused, ll(idim),
     1               icell(0:ncell,0:ncell,0:ncell)
      common /llisr/ xcell(0:ncell),ycell(0:ncell),zcell(0:ncell),
     2               xmapl(0:ncell), xmapu(0:ncell), xlow, xup,
     3               ymapl(0:ncell), ymapu(0:ncell), ylow, yup,
     4               zmapl(0:ncell), zmapu(0:ncell), zlow, zup
      common /sorth/ indxh(idim), irankh(idim)
      common /sortc/ irankx(idim), iranky(idim), irankz(idim),
     1               nslice, nx, ny, nz
c
      common /work  / x(idim), y(idim), z(idim), h(idim),
     1                dm(idim)
c
c--check if enough space in list
c
      ratio=float(npart)/float(nocc)
      ncellused=nint(ratio**0.3333333333+1)
      nslice=nocc*ncellused*ncellused
      nx=npart/nslice
      ny=nx
      nz=nx
      if(ncellused.gt.ncell) then
         write(*,*)'not enough space in linked list!'
         stop
      endif
c
c--rank particles in increasing coordinate
c  ---------------------------------------
c
c  a) x coordinate
c
      call indexx(npart,x,indxh)
      xlow=x(indxh(1))
      xup=x(indxh(npart))
      call rank(npart,indxh,irankx)
c
c  b) y coordinate
c
      call indexx(npart,y,indxh)
      ylow=y(indxh(1))
      yup=y(indxh(npart))
      call rank(npart,indxh,iranky)
c
c  c) z coordinate
c
      call indexx(npart,z,indxh)
      zlow=z(indxh(1))
      zup=z(indxh(npart))
      call rank(npart,indxh,irankz)
c
c--set all cells to zero
c  ---------------------
c
      do i=1,idim
         ll(i)=0
      enddo
      do i=0,ncellused
         do j=0,ncellused
            do k=0,ncellused
               icell(i,j,k)=0
            enddo
         enddo
      enddo
      do i=0,ncellused
         xmapu(i)=-1e20
         xmapl(i)=1e20
      enddo
      do i=0,ncellused
         ymapu(i)=-1e20
         ymapl(i)=1e20
      enddo
      do i=0,ncellused
         zmapu(i)=-1e20
         zmapl(i)=1e20
      enddo
c
c--build linked list and cell boundary map
c  ---------------------------------------
c
      do ipart=1,npart
         ratioi=float(irankx(ipart)-1)/float(nslice)
         i=int(ratioi)
         ratioj=float(iranky(ipart)-1)/float(nslice)
         j=int(ratioj)
         ratiok=float(irankz(ipart)-1)/float(nslice)
         k=int(ratiok)
         ll(ipart)=icell(i,j,k)
         icell(i,j,k)=ipart
c
         xmapu(i)=max(xmapu(i),x(ipart))
         xmapl(i)=min(xmapl(i),x(ipart))
         ymapu(j)=max(ymapu(j),y(ipart))
         ymapl(j)=min(ymapl(j),y(ipart))
         zmapu(k)=max(zmapu(k),z(ipart))
         zmapl(k)=min(zmapl(k),z(ipart))
      enddo
c
c--sort particles according to smoothing length
c  --------------------------------------------
c
      call indexx(npart,h,indxh)
      call rank(npart,indxh,irankh)
c
      return
      end
      subroutine indexx(n,arr,indx)
c************************************************************
c                                                           *
c  Sort arrays using the Quicksort algorithm (Press 2nd)    *
c                                                           *
c************************************************************
c
      parameter (m=7,nstack=200)
c
      dimension arr(n), indx(n)
      dimension istack(nstack)
c
      do j=1,n
         indx(j)=j
      enddo
c
      jstack=0
      l=1
      ir=n
c
    1 continue
      if(ir-l.lt.m)then
         do j=l+1,ir
            indxt=indx(j)
            a=arr(indxt)
            do i=j-1,1,-1
               if(arr(indx(i)).le.a)go to 2
               indx(i+1)=indx(i)
            enddo
            i=0
    2       indx(i+1)=indxt
         enddo
         if(jstack.eq.0)return
         ir=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack - 2
      else
         k=(l+ir)/2
         itemp=indx(k)
         indx(k)=indx(l+1)
         indx(l+1)=itemp
         if(arr(indx(l+1)).gt.arr(indx(ir)))then
            itemp=indx(l+1)
            indx(l+1)=indx(ir)
            indx(ir)=itemp
         end if
         if(arr(indx(l)).gt.arr(indx(ir)))then
            itemp=indx(l)
            indx(l)=indx(ir)
            indx(ir)=itemp
         end if
         if(arr(indx(l+1)).gt.arr(indx(l)))then
            itemp=indx(l+1)
            indx(l+1)=indx(l)
            indx(l)=itemp
         end if
c
         i=l+1
         j=ir
         indxt=indx(l)
         a=arr(indxt)
c
    3    continue
         i=i+1
         if(arr(indx(i)).lt.a)go to 3
    4    continue
         j=j-1
         if(arr(indx(j)).gt.a)go to 4
         if(j.lt.i)go to 5
c
         itemp=indx(i)
         indx(i)=indx(j)
         indx(j)=itemp
         go to 3
c
    5    continue
         indx(l)=indx(j)
         indx(j)=indxt
         jstack=jstack + 2
         if(jstack.gt.nstack)then
            write(*,*)'nstack too small in indexx!'
            stop
         end if
c
         if(ir-i+1.ge.j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
         else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
         end if
      end if
      go to 1
c
      end
      subroutine rank(n,indx,irank)
c*************************************************************
c                                                            *
c  subroutine ranking sorted particles.                      *
c                                                            *
c*************************************************************
c
      dimension indx(n), irank(n)
c
      do j=1,n
         irank(indx(j))=j
      enddo
c
      return
      end
      subroutine eost (imat,ui,rhoi,vsoundi,pri)
c************************************************
c                                               *
c   compute the tillotson equation              *
c   of state                                    *
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

      parameter (nmatmax=21)
c
      common /till / rozero(nmatmax), smalla(nmatmax), 
     1               uzero(nmatmax),  smallb(nmatmax), 
     2               capa(nmatmax), es(nmatmax), 
     3               esp(nmatmax), alpha(nmatmax), 
     4               beta(nmatmax), capb(nmatmax)
c
      rhoi2=rhoi*rhoi
c
c--normalisation of the constants
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
      eta=rhoi/rozn
      xmu=eta-1.
      rap=rozn/rhoi
c
c--condensed phase, compute pressure and sound speed
c
      p1=smallbj/(ui/(ezn*eta*eta)+1.)
      pressc=(smallaj+p1)*ui*rhoi+capan*xmu+capbn*xmu*xmu
c        if(abs(pressc).lt.capan*roundoff)pressc=0.
      p1=ui/(ezn*eta*eta)+1.
      c2=p1*p1
      c3=p1-1.
      c4=smallbj*ui*(3.*c3+1.)/c2
      dpdi=smallaj*rhoi+smallbj*rhoi/c2
      dpdrho=smallaj*ui+capan/rozn+2.*capbn*xmu/rozn+c4
      cvelc=dpdrho+dpdi*pressc/rhoi2
c
c--expanded phase, compute pressure and sound speed
c
      p1=ui/ezn/eta**2+1.
      p2=smallaj*ui*rhoi
      p3=smallbj*ui*rhoi/p1
      p4=capan*xmu
      vow=1./eta
      p6=betaj*(vow-1.)
      p6=min(p6,70.)
      p6=exp(-p6)
      p7=alphaj*(vow-1.)**2
      p7=min(p7,70.)
      p7=exp(-p7)
      pressv=p2+(p3+p4*p6)*p7
      c1=smallaj*ui
      c2=p1**2
      c3=p1-1.
      c4=smallbj*ui*(3.*c3+1.)/c2
      c5=p6*p7*capan
      c6=2.*alphaj*(vow-1.)
      dpdrho=c1+p7*c4+p7*p3*rozn*c6/rhoi2+c5*(1./rozn+xmu*rozn/
     1          rhoi2*(c6+betaj))
      dpdi=smallaj*rhoi+p7/c2*smallbj*rhoi
      cvelv=dpdrho+dpdi*pressv/rhoi2
      if (cvelv.lt.0.) cvelv=0.
c
c--pick up the right state
c
      pri=pressc
      vsoundi=cvelc
      if(rap.ge.1..and.ui.gt.espn) then
         pri=pressv
         vsoundi=cvelv
      end if
      if(rap.ge.1.and.ui.gt.esn.and.ui.lt.espn) then
         pri=(pressv*(ui-esn) + pressc*(espn-ui))/difesn
         vsoundi=(cvelv*(ui-esn) + cvelc*(espn-ui))/difesn
      end if
      cmin=.25*capan/rozn
      vsoundi=sqrt(max(cmin,vsoundi))
c
      return
      end
      subroutine eosa(imat,tempi,rhoi,pri,ui,csi,kpai)
c************************************************
c                                               *
c  ANEOS equation of state                      *
c                                               *
c   imat  material                              *
c   ----  --------                              *
c                                               *
c    1   granite        11   dry tuff           *
c    2   dunite         12   alluvium           *
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
      double precision tti,rrhoi,ppi,uui,ssi,ccvi,ddpdti,
     1                 ddpdri,ffkrosi,ccsi,ffve,ffva
c
      data evtemp/11604.8/
c
c--double precision and transform T in ev
c  --------------------------------------
c
      rrhoi=rhoi
      tti=tempi/evtemp
c
c--compute pressure and all other quantities
c  -----------------------------------------
c
      call aneos (tti,rrhoi,ppi,uui,ssi,ccvi,ddpdti,
     1            ddpdri,ffkrosi,ccsi,kpai,imat,ffve,ffva)
      write(*,*)'calling aneos: ',imat,tti,rrhoi,ppi
c
c--transform in single precision
c  -----------------------------
c
      ui=uui
      pri=ppi
      csi=ccsi
c
      return
      end
      subroutine cleanup(npart,npset,np)
c***************************************************************************
c                                                                          *
c  clean up particle list                                                  *
c                                                                          *
c***************************************************************************
c
      parameter (idim=2300000)
c
      common /sph   / x(idim)     , y(idim) , z(idim)    ,
     1                vx(idim)    , vy(idim), vz(idim)   ,
     2                h(idim)     , id(idim), pmass(idim),
     3                matter(idim), p(idim) , T(idim)    ,
     4                rho(idim)   , u(idim)
      common /solids/ dm(idim)    , rft(idim), sxx(idim) ,
     1                sxy(idim)   , sxz(idim), syy(idim) , 
     2                syz(idim)
      common /integ/ ipind(idim), ipdam(idim), iprem(idim), 
     1               iremove(idim)
c
      do i=1,idim
         iremove(i)=0
      enddo
c
c--set quantities
c
      do i=1,np
c
c  a) damage
c
         if(ipdam(i).eq.1)then
            dm(ipind(i))=1.0
         endif
c
c  b) removal flag
c
         if(iprem(i).eq.1)then
            iremove(ipind(i))=1
         endif
c
      enddo
c
c--remove particles and plug holes 
c
      ishift=0
      do i=1,npart
         if(iremove(i).eq.0)then
            x(i-ishift)=x(i)
            y(i-ishift)=y(i)
            z(i-ishift)=z(i)
            vx(i-ishift)=vx(i)
            vy(i-ishift)=vy(i)
            vz(i-ishift)=vz(i)
            h(i-ishift)=h(i)
            id(i-ishift)=id(i)
            pmass(i-ishift)=pmass(i)
            matter(i-ishift)=matter(i)
            p(i-ishift)=p(i)
            T(i-ishift)=T(i)
            rho(i-ishift)=rho(i)
            u(i-ishift)=u(i)
            dm(i-ishift)=dm(i)
            rft(i-ishift)=rft(i)
            sxx(i-ishift)=sxx(i)
            sxy(i-ishift)=sxy(i)
            sxz(i-ishift)=sxz(i)
            syy(i-ishift)=syy(i)
            syz(i-ishift)=syz(i)
         else
            ishift=ishift + 1
         end if
      enddo
c
      npset=npset - ishift
c
      return
      end
      subroutine xdrheader(iflag,filename,title,numvars,npart,time)
c************************************************************
c                                                           *
c  read/write header of xdr file                            *
c  iflag=1  : read                                          *
c  ifalg=2  : write                                         *
c                                                           *
c************************************************************
c
      parameter (maxvars=21)
c
      character*70 title
      character*15 filename
      character*56 varnam(maxvars), varnam2(maxvars)
c
      data varnam /'x ','y ','z ','vx ','vy ','vz ','h ','id ',
     1             'mass ','mat ','p ','T ','rho ','u ','dm ',
     2             'rft ','sxx ','sxy ','sxz ','syy ','syz '/
c
c--header
c
c  a) read
c 
      if(iflag.eq.1)then
         j=index(filename,' ')
         filename(j:j)=char(0)
         call openreadxdrfile(filename,varnam2,numvars,npart,time,
     1                        title)
         write(*,*)varnam2
      endif
c
c  b) write (add a zero to all fortran strings because of C)
c
      if(iflag.eq.2)then
         j=index(title,' ')
         title(j:j)=char(0)
         do i=1,numvars
            j=index(varnam(i),' ')
            varnam(i)(j:j)=char(0)
         enddo
         j=index(filename,' ')
         filename(j:j)=char(0)
         call openwritexdrfile(filename,varnam,numvars,npart,time,
     1                         title)
      endif
c
      return
      end
      subroutine tillinit
c************************************************
c                                               *
c  subroutine reading data for the Tillotson    * 
c  equation of state                            *
c                                               *
c************************************************
c
      parameter (nmatmax=21)
c
      common /till / rozero(nmatmax), smalla(nmatmax), 
     1               uzero(nmatmax),  smallb(nmatmax), 
     2               capa(nmatmax), es(nmatmax), 
     3               esp(nmatmax), alpha(nmatmax), 
     4               beta(nmatmax), capb(nmatmax)
c
      character*70 filen
      character*50 dumy
      character*2 path
c
      data path /'./'/
c
c--open data file
c
      filen=path // 'tillotson.input'
      open(30,file=filen)
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
      imat=-1
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
      do i=1,nline
         imatmin=(i-1)*5 + 1 
         imatmax=min(nmat,imatmin+4)
         read(30,*)(esp(j),j=imatmin,imatmax)
      enddo
c
      close(30)
c
      return
      end
      subroutine xdrparticle(iflag,ipart,numvars)
c************************************************************
c                                                           *
c  read/write all variables for particle ipart              *
c  iflag=1    : read                                        *
c  iflag=2    : write                                       *
c                                                           *
c************************************************************
c
      parameter (idim=2300000)
      parameter (maxvars=21)
c
      common /sph   / x(idim)     , y(idim) , z(idim)    ,
     1                vx(idim)    , vy(idim), vz(idim)   ,
     2                h(idim)     , id(idim), pmass(idim),
     3                matter(idim), p(idim) , T(idim)    ,
     4                rho(idim)   , u(idim)
      common /solids/ dm(idim)    , rft(idim), sxx(idim) ,
     1                sxy(idim)   , sxz(idim), syy(idim) ,
     2                syz(idim)
c
      dimension particle(maxvars)
c
c--read attributes for one particle
c
      if(iflag.eq.1)then
         call readparticle(particle,numvars)
         x(ipart)=particle(1)
         y(ipart)=particle(2)
         z(ipart)=particle(3)
         vx(ipart)=particle(4)
         vy(ipart)=particle(5)
         vz(ipart)=particle(6)
         h(ipart)=particle(7)
         id(ipart)=particle(8)
         pmass(ipart)=particle(9)
         matter(ipart)=particle(10)
         p(ipart)=particle(11)
         T(ipart)=particle(12)
         rho(ipart)=particle(13)
         u(ipart)=particle(14)
         dm(ipart)=particle(15)
         rft(ipart)=particle(16)
         sxx(ipart)=particle(17)
         sxy(ipart)=particle(18)
         sxz(ipart)=particle(19)
         syy(ipart)=particle(20)
         syz(ipart)=particle(21)
      endif
c
c--write attributes for one particle
c
      if(iflag.eq.2)then
         particle(1)=x(ipart)
         particle(2)=y(ipart)
         particle(3)=z(ipart)
         particle(4)=vx(ipart)
         particle(5)=vy(ipart)
         particle(6)=vz(ipart)
         particle(7)=h(ipart)
         particle(8)=id(ipart)
         particle(9)=pmass(ipart)
         particle(10)=matter(ipart)
         particle(11)=p(ipart)
         particle(12)=T(ipart)
         particle(13)=rho(ipart)
         particle(14)=u(ipart)
         particle(15)=dm(ipart)
         particle(16)=rft(ipart)
         particle(17)=sxx(ipart)
         particle(18)=sxy(ipart)
         particle(19)=sxz(ipart)
         particle(20)=syy(ipart)
         particle(21)=syz(ipart)
         call writeparticle(particle,numvars)
      endif
c
      return
      end
