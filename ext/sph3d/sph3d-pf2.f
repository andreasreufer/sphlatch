      subroutine tmpoffaneosinit
c************************************************************
c                                                           *
c  initialize the ANEOS equation of state                   *
c                                                           *
c************************************************************
c
      parameter (nmatmax=21)
c
      common /fileos/ klst, kinp
      common /numat/ nmat, ntype(nmatmax)
c
      dimension nntype(nmatmax)
c
c--open aneos related files
c  ------------------------
c
      klst=20
      kinp=21
      open(klst,file='aneos.output') 
      open(kinp,file='aneos.input')
c
c--set material types
c  ------------------
c
      do i=1,nmat
         nntype(i)=-ntype(i)
      enddo
c
c--initialize
c  ----------
c
      call aneos2 (1,nmat,0,nntype)
c
c--close open files
c  ----------------
c
      close(20)
      close(21)
c
      return
      end
      subroutine angmom
c************************************************************
c                                                           *
c  compute the total angular momentum of the system         *
c                                                           *
c************************************************************
c
      double precision dangto, dangx, dangy, dangz, pmi, xi, 
     1                 yi, zi, dvxi, dvyi, dvzi, cmxtot, cmytot,
     2                 cmztot, cmtotvx, cmtotvy, cmtotvz
c
      parameter (idim=300000)
c
      common /parti/ npart
      common /partr/ x(idim),y(idim),z(idim),vx(idim),
     1               vy(idim),vz(idim),u(idim),rho(idim),h(idim)
      common /carai/ matter(idim), itype(idim)
      common /carar/ pmass(idim)
      common /out1 / cmx1, cmy1, cmz1, vcmx1, vcmy1, vcmz1, romean1,
     1               romax1, dmax1, zmax1, tmean1, tmax1, hmi1, hma1,
     2               damax1, damean1
      common /out2 / cmx2, cmy2, cmz2, vcmx2, vcmy2, vcmz2, romean2,
     1               romax2, dmax2, zmax2, tmean2, tmax2, hmi2, hma2,
     2               damax2, damean2
      common /bodyi/ n1, n2
      common /bodyr/ fmas1, fmas2
      common /angm / angx, angy, angz, angto
      common /logun/ iprint, iterm, idisk1, idisk2, idisk3
c
c--compute center of mass position and velocity
c  --------------------------------------------
c
      xmtot=fmas1+fmas2
      cmtotx=(fmas1*cmx1+fmas2*cmx2)/xmtot
      cmtoty=(fmas1*cmy1+fmas2*cmy2)/xmtot
      cmtotz=(fmas1*cmz1+fmas2*cmz2)/xmtot
      cmtotvx=(fmas1*vcmx1+fmas2*vcmx2)/xmtot
      cmtotvy=(fmas1*vcmy1+fmas2*vcmy2)/xmtot
      cmtotvz=(fmas1*vcmz1+fmas2*vcmz2)/xmtot
c
c--compute total angular momentum
c  ------------------------------
c
      dangx=0.d0
      dangy=0.d0
      dangz=0.d0
      do i=1,npart
         pmi=pmass(i)
         xi=x(i)-cmtotx
         yi=y(i)-cmtoty
         zi=z(i)-cmtotz
         vxi=vx(i)-cmtotvx
         vyi=vy(i)-cmtotvy
         vzi=vz(i)-cmtotvz
         dangx=dangx + pmi*(yi*vzi - vyi*zi)
         dangy=dangy + pmi*(vxi*zi - xi*vzi)
         dangz=dangz + pmi*(xi*vyi - vxi*yi)
      enddo
c
c--total angular momentum
c  ----------------------
c
      dangto=sqrt(dangx**2+dangy**2+dangz**2)
      angx=dangx
      angy=dangy
      angz=dangz
      angto=dangto
c
      return
      end
      subroutine damping(npart,iupdat,vx,vy,vz,fx,fy,fz)
c************************************************************
c                                                           *
c  subroutine adding an overall damping to the momentum     *
c  equation.                                                *
c                                                           *
c************************************************************
c
      parameter (idim=300000)
c
      common /dissi/ damp
c
      dimension iupdat(idim)
      dimension vx(idim), vy(idim), vz(idim)
      dimension fx(idim), fy(idim), fz(idim)
c
c--add general damping
c  -------------------
c
      do i=1,npart
         ipart=iupdat(i)
         fx(ipart)=fx(ipart) - damp*vx(ipart)
         fy(ipart)=fy(ipart) - damp*vy(ipart)
         fz(ipart)=fz(ipart) - damp*vz(ipart)
      enddo
c
      return
      end
      subroutine density(npart,iupdat,rho,drho)
c************************************************************
c                                                           *
c  this routine computes the change in density solving the  *
c  continuity equation                                      *
c                                                           *
c************************************************************
c
      parameter (idim=300000)
c
      dimension iupdat(idim),rho(idim),drho(idim)
c
      common /srate/ epsxx(idim), epsyy(idim), epszz(idim),
     1               epsxy(idim), epsxz(idim), epsyz(idim)
      common /divve/ divv(idim), divmax1, divmax2
      common /typef/ istrength, ifrac, ieos, iker
      common /carai/ matter(idim), itype(idim)
      common /scode/ istr, i1d0, i1d1, i2d0, i2d1
c
      do i=1,npart
         ipart=iupdat(i)
         drho(ipart)=-rho(ipart)*divv(ipart)
      enddo
c
      return
      end
      subroutine deriv(t,npart,istep,x,y,z,vx,vy,vz,u,rho,h,
     1                 sxx,syy,sxy,sxz,syz,dm,fx,fy,fz,du,drho,dh,
     2                 dsxx,dsyy,dsxy,dsxz,dsyz,ddm)
c******************************************************************
c                                                                 *
c  This subroutine is the driver for all right-hand-side terms    *
c  for all differential equations.                                *
c                                                                 *
c******************************************************************
c
      double precision umass
      double precision cmx1, cmy1, cmz1, cmvx1, cmvy1, cmvz1,
     1                 cmx2, cmy2, cmz2, cmvx2, cmvy2, cmvz2,
     2                 totm, pmi, vrelob, distc
c
      parameter (idim=300000)
      parameter (ilist=idim)
c
      dimension x(idim),y(idim),z(idim),vx(idim),vy(idim),vz(idim),
     1          u(idim),rho(idim),h(idim),dm(idim)
      dimension sxx(idim),syy(idim),sxy(idim),sxz(idim),syz(idim)
      dimension dsxx(idim),dsyy(idim),dsxy(idim),dsxz(idim),
     1          dsyz(idim)
      dimension fx(idim),fy(idim),fz(idim),du(idim),drho(idim),
     1          dh(idim),ddm(idim)
      dimension iupdat(idim)
c
      common /neigh/ neilist(ilist), ilen(idim)
      common /numpa/ avalpha, avbeta
      common /tempe/ temp(idim)
      common /carai/ matter(idim), itype(idim)
      common /carar/ pmass(idim)
      common /bodyi/ n1, n2
      common /bodyr/ fmas1, fmas2
      common /divve/ divv(idim), divmax1, divmax2
      common /eosq / pr(idim), vsound(idim), vsmax
      common /ener1/ dq(idim)
      common /reduc/ vonmises(idim), reduce(idim)
      common /damag/ ratioft(idim), ratiofs(idim)
      common /units/ umass, udist, udens, utime, uergg, uergcc
      common /typef/ istrength, ifrac, ieos, iker
      common /dissi/ damp
      common /bndr / xmin,xmax,ymin,ymax,zmin,zmax
c
c--compute thermodynamical properties
c  ----------------------------------
c
      if(ieos.eq.1)call eost(npart,rho,u)
      if(ieos.eq.2)call eosa(npart,rho,u)
c
c--plastic yielding
c  ----------------
c
      call plastic(npart,u,dm,sxx,syy,sxy,sxz,syz)
c
c--compute reduction factor for stresses
c  -------------------------------------
c
      call reduces(npart,dm)
c
c--build linked list
c  -----------------
c
      call mllist (npart,iupdat,x,y,z,h)
c
c--initialize arrays
c  -----------------
c
      call zeros(npart,rho,sxx,syy,sxy,sxz,syz,fx,fy,fz)
c
c--flag particles
c  --------------
c
      call flagp(npart)
c
c--compute interactions
c  --------------------
c
      call forces(npart,iupdat,x,y,z,vx,vy,vz,
     1            rho,h,sxx,syy,sxy,sxz,syz,fx,fy,fz)
c
c--compute all local quantities
c  ----------------------------
c
c  a)compute density derivative
c
      call density(npart,iupdat,rho,drho)
c
c--compute energy derivative 
c
      call energ(npart,iupdat,du)
c
c--compute derivatives of deviatoric stress tensor
c
      call deviator(npart,iupdat,sxx,syy,sxy,
     1              sxz,syz,dsxx,dsyy,dsxy,dsxz,dsyz)
c
c--integrate fracture model
c
      call fracture(npart,iupdat,dm,sxx,syy,sxy,sxz,syz,ddm)
c
c--add general damping
c
      if(damp.ne.0.) then
         call damping(npart,iupdat,vx,vy,vz,fx,fy,fz)
      end if
c
c--compute density averaged divergence of the velocity
c
      divmax1=0.
      divmax2=0.
      robar1=0.
      robar2=0.
      do i=1,npart
         ipart=iupdat(i)
         if(ipart.le.n1)then
c           divmax1=divmax1+rho(ipart)*divv(ipart)**2
c           robar1=robar1 + rho(ipart)
            divmax1=max(divmax1,abs(divv(ipart)))
         else
c           divmax2=divmax2+rho(ipart)*divv(ipart)**2
c           robar2=robar2 + rho(ipart)
            divmax2=max(divmax2,abs(divv(ipart)))
         end if
      enddo
c     divmax1=sqrt(divmax1/robar1)
c     if(n2.ne.0)then
c        divmax2=sqrt(divmax2/robar2)
c     end if
c
c--compute time derivative of smoothing length
c  -------------------------------------------
c
      call hdot(npart,iupdat,h,dh)
c
c--compute some diagnostics
c
c a) objects are approaching or not
c
      if(n2.ne.0)then
         cmx1=0.d0
         cmy1=0.d0
         cmz1=0.d0
         cmvx1=0.d0
         cmvy1=0.d0
         cmvz1=0.d0
         totm=0.d0
         do i=1,n1
            pmi=pmass(i)
            totm=totm + pmi
            cmx1=cmx1 + pmi*x(i)
            cmy1=cmy1 + pmi*y(i)
            cmz1=cmz1 + pmi*z(i)
            cmvx1=cmvx1 + pmi*vx(i)
            cmvy1=cmvy1 + pmi*vy(i)
            cmvz1=cmvz1 + pmi*vz(i)
         enddo
         cmx1=cmx1/totm
         cmy1=cmy1/totm
         cmz1=cmz1/totm
         cmvx1=cmvx1/totm
         cmvy1=cmvy1/totm
         cmvz1=cmvz1/totm
c
         cmx2=0.d0
         cmy2=0.d0
         cmz2=0.d0
         cmvx2=0.d0
         cmvy2=0.d0
         cmvz2=0.d0
         totm=0.d0
         do i=n1+1,npart
            pmi=pmass(i)
            totm=totm + pmi
            cmx2=cmx2 + pmi*x(i)
            cmy2=cmy2 + pmi*y(i)
            cmz2=cmz2 + pmi*z(i)
            cmvx2=cmvx2 + pmi*vx(i)
            cmvy2=cmvy2 + pmi*vy(i)
            cmvz2=cmvz2 + pmi*vz(i)
         enddo
         cmx2=cmx2/totm
         cmy2=cmy2/totm
         cmz2=cmz2/totm
         cmvx2=cmvx2/totm
         cmvy2=cmvy2/totm
         cmvz2=cmvz2/totm
c
         distc=sqrt((cmx1-cmx2)**2+(cmy1-cmy2)**2+
     1             (cmz1-cmz2)**2)
         vrelob=((cmvx1-cmvx2)*(cmx1-cmx2) + 
     1           (cmvy1-cmvy2)*(cmy1-cmy2) +
     2           (cmvz1-cmvz2)*(cmz1-cmz2))/distc
      else
         vrelob=0.0d0
      endif
c
c b) other statistics
c
      hmin=1.e30
      hmax=-hmin
      ilenmax=0
      ilenmin=1000000
      prmax=-1e30
      prmin=1.e30
      rhomin=1.e30
      rhomax=-rhomin
      umin=1.e30
      umax=-umin
      vonmean=0.
      vonmin=1.e30
      vonmax=-vonmin
      nfrac=0
      dmmin=1e30
      dmmax=-dmmin 
      dmmean=0.
      ratioftmax=-1e30
      ratiofsmax=-1e30
      do i=1,npart
         ipart=iupdat(i)
         if(ipart.le.npart)then
            if(ilen(ipart).gt.ilenmax)then
               ilenmax=ilen(ipart)
               jlenmax=ipart
            end if
            if(ilen(ipart).lt.ilenmin)then
               ilenmin=ilen(ipart)
               jlenmin=ipart
            end if
            prmax=max(prmax,pr(ipart))
            prmin=min(prmin,pr(ipart))
            if(h(ipart).lt.hmin)then
               hmin=h(ipart)
            end if
            if(h(ipart).gt.hmax)then
               hmax=h(ipart)
            end if
            rhomin=min(rhomin,rho(ipart))
            rhomax=max(rhomax,rho(ipart))
            umin=min(umin,u(ipart))
            umax=max(umax,u(ipart))
            vonmin=min(vonmin,vonmises(ipart))
            vonmax=max(vonmax,vonmises(ipart))
            vonmean=vonmean + vonmises(ipart)
            if(dm(ipart).gt.1.e-5) then
               nfrac=nfrac + 1
            endif
            dmmin=min(dmmin,dm(ipart)**3)
            dmmax=max(dmmax,dm(ipart)**3)
            ratioftmax=max(ratioft(ipart),ratioftmax)
            ratiofsmax=max(ratiofs(ipart),ratiofsmax)
         end if
      enddo
      write(*,*)'deriv: nb. of dt  ',istep
      write(*,*)'deriv: neighbors  ',jlenmin,ilenmin,jlenmax,ilenmax
      write(*,*)'deriv: smooth. l. ',hmin,hmax
      write(*,*)'deriv: pressure   ',prmax,prmin
      write(*,*)'deriv: density    ',rhomin,rhomax
      write(*,*)'deriv: energy     ',umin,umax
      write(*,*)'deriv: divergence ',divmax1, divmax2
      write(*,*)'deriv: plastic    ',vonmin,vonmax,vonmean/npart
      write(*,*)'deriv: fracture   ',nfrac,dmmin,dmmax
      write(*,*)'deriv: rt,rs      ',ratioftmax, ratiofsmax
      write(*,*)'deriv: v rel.     ',vrelob
c
      return
      end
      subroutine deviator(npart,iupdat,sxx,syy,sxy,sxz,
     1                    syz,dsxx,dsyy,dsxy,dsxz,dsyz)
c***************************************************************
c                                                              *
c  subroutine to compute the derivative of the deviatoric      *
c  stress tensor.                                              *
c                                                              *
c***************************************************************
c
      parameter(idim=300000)
      parameter (nmatmax=21)
c
      dimension iupdat(idim)
      dimension sxx(idim),syy(idim),sxy(idim),sxz(idim),
     1          syz(idim)
      dimension dsxx(idim),dsyy(idim),dsxy(idim),dsxz(idim),
     1          dsyz(idim)
c
      common /bodyi/ n1, n2
      common /carai/ matter(idim), itype(idim)
      common /carar/ pmass(idim)
      common /srate/ epsxx(idim), epsyy(idim), epszz(idim),
     1               epsxy(idim), epsxz(idim), epsyz(idim)
      common /reduc/ vonmises(idim), reduce(idim)
      common /rrate/ rxy(idim), rxz(idim), ryz(idim)
      common /shear/ xmu(nmatmax)
      common /divve/ divv(idim), divmax1, divmax2
      common /typef/ istrength, ifrac, ieos, iker
c
c--if no material strength, zero derivatives
c  -----------------------------------------
c
      if(istrength.eq.0)then
         do i=1,npart
            ipart=iupdat(i)
            dsxx(ipart)=0.
            dsyy(ipart)=0.
            dsxy(ipart)=0.
            dsxz(ipart)=0.
            dsyz(ipart)=0.
         enddo
         return
      end if
c
c--compute derivative of deviatoric stress tensor
c  ----------------------------------------------
c
      do i=1,npart
         ipart=iupdat(i)
         div3=(epsxx(ipart)+epsyy(ipart)+epszz(ipart))/3.
         imat=matter(ipart)
         xmu2=2.*xmu(imat)
c
         sxxi=sxx(ipart)
         syyi=syy(ipart)
         szzi=-sxxi-syyi
         sxyi=sxy(ipart)
         sxzi=sxz(ipart)
         syzi=syz(ipart)
c
         epsxxi=epsxx(ipart)
         epsyyi=epsyy(ipart)
         epsxyi=epsxy(ipart)
         epsxzi=epsxz(ipart)
         epsyzi=epsyz(ipart)
c
         rxyi=rxy(ipart)
         rxzi=rxz(ipart)
         ryzi=ryz(ipart)
c
         dsxx(ipart)=xmu2*(epsxxi-div3)+2.*sxyi*rxyi+
     1                  2.*sxzi*rxzi
         dsyy(ipart)=xmu2*(epsyyi-div3)-2.*sxyi*rxyi+
     1                  2.*syzi*ryzi
         dsxy(ipart)=xmu2*epsxyi+(syyi-sxxi)*rxyi+
     1                  sxzi*ryzi+syzi*rxzi
         dsxz(ipart)=xmu2*epsxzi+(szzi-sxxi)*rxzi-
     1                  sxyi*ryzi+syzi*rxyi
         dsyz(ipart)=xmu2*epsyzi+(szzi-syyi)*ryzi-
     1                  sxyi*rxzi-sxzi*rxyi
      enddo
c
c--zero-out derivative in impactor if required
c
      if(istrength.eq.2)then
         do i=1,npart
            ipart=iupdat(i)
            if(ipart.gt.n1)then
               dsxx(ipart)=0.
               dsyy(ipart)=0.
               dsxy(ipart)=0.
               dsxz(ipart)=0.
               dsyz(ipart)=0.
            endif
         enddo
      endif
c
      return
      end
      subroutine endrun
c************************************************************
c                                                           *
c  this subroutine ends the runstream.                      *
c                                                           *
c************************************************************
c
      common /logun/ iprint, iterm, idisk1, idisk2, idisk3
      common /actio/ namerun
c
      character*20 namerun
c
c--get time and date
c
      call getime(ih,im,is,heure)
      call getdat(j,m,ia) 
      call getused(tused)
c
c--write end page
c
      iw=1
      do i=1,iw
         write(iprint,100)namerun,j,m,ia,ih,im,is
  100    format(/////,1x,'SPH run ',a20,' ended normally on  : ',i2,'/' 
     1          ,i2,'/',i4,'  at ',i2,' h. ',i2,' min. ',i2,' sec.')
         write(iprint,101)tused/60.
  101    format(1x,'cpu time used for this run :',f9.3,' min.')
      enddo
      close(iprint)
c
c--write ok file
c
      open(27,file='ok')
      write(27,*)'termination ok!'
      close(27)
c
      stop
      end
      subroutine energ(npart,iupdat,du)
c************************************************************
c                                                           *
c  this routine computes the change in internal energy      *
c                                                           *
c************************************************************
c
      parameter(idim=300000)
c
      dimension iupdat(idim), du(idim)
c 
      common /ener1/ dq(idim)
c
c--set specific internal energy derivative 
c
      do i=1,npart
         ipart=iupdat(i)
         du(ipart)=0.5*dq(ipart)
      enddo
c
      return
      end
      subroutine eosa(npart,rho,u)
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
      parameter (idim=300000)
      parameter (nmatmax=21)
c
      double precision pi,csi,ti,rhoi,ui,prlim
c
      dimension rho(idim), u(idim)
c
      common /bodyi/ n1, n2
      common /carai/ matter(idim), itype(idim)
      common /carar/ pmass(idim)
      common /cgas / gamma
      common /eosq / pr(idim), vsound(idim), vsmax
      common /fileos/ klst, kinp
      common /gasfr/ fgas1, fgas2
      common /logun/ iprint, iterm, idisk1, idisk2, idisk3
      common /numat/ nmat, ntype(nmatmax)
      common /tempe/ temp(idim)
      common /wban2/ nntype(nmatmax)
c
      data evtemp/11604.8/
      data prlim/-1.d15/, tini/300./
c
c--initialisation
c
      fgas1=0.
      fgas2=0.
      vsmax=0.
      tempmax=0.
      prmax=-1e30
c
c--compute pressure and sound speed for all particles
c
      imax=1
      do ipart=1,npart
c
c--transform density and specific internal energy to cgs
c
         rhoi=rho(ipart)
         ui=u(ipart)
         mati=matter(ipart)
         ti=temp(ipart)/evtemp
c
c--compute pressure by iterating on temperature until energy
c  given by aneos equals energy given by the code.
c  Use the Brent method (Numerical recipes p 251)
c
         call rooten (ipart,rhoi,ui,mati,ti,pi,csi,kpai)
c
c--store values 
c
         pr(ipart)=max(pi,prlim)
         vsound(ipart)=csi
         temp(ipart)=ti*evtemp
         vsmax=max(vsmax,vsound(ipart))
         if(abs(pr(ipart)).gt.prmax)then
            imax=ipart
            prmax=abs(pr(ipart))
         end if
         if(kpai.eq.6) then
            if(ipart.le.n1)then
               fgas1=fgas1 + pmass(ipart)
            else
               fgas2=fgas2 + pmass(ipart)
            endif
         endif
      enddo
c
      return
      end
      subroutine eosmg(npart,rho)
c************************************************************
c                                                           *
c  subroutine implementing a Murnaghan eos                  *
c                                                           *
c************************************************************
c
      parameter (idim=300000)
c
      common /eosq / pr(idim), vsound(idim), vsmax
c
      dimension rho(idim)
c
c--define material constants here ice
c
      rho0=0.94
      n=1
      bulk=6.5e10/n
      vsound0=sqrt(bulk/rho0)
c
      do ipart=1,npart
         r1=rho(ipart)/rho0
         pr(ipart)=bulk*(r1**n - 1.)
         vsound(ipart)=vsound0
      enddo
c
      return
      end
      subroutine eost (npart,rho,u)
c************************************************
c                                               *
c  Tillotson equation of state                  *
c                                               *
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
c************************************************
c
      parameter (idim=300000)
      parameter (nmatmax=21)
c
      dimension rho(idim), u(idim)
c
      common /till / rozero(nmatmax), smalla(nmatmax), uzero(nmatmax), 
     1               smallb(nmatmax), capa(nmatmax), es(nmatmax), 
     2               esp(nmatmax), alpha(nmatmax), beta(nmatmax), 
     3               capb(nmatmax)
      common /carai/ matter(idim), itype(idim)
      common /carar/ pmass(idim)
      common /eosq / pr(idim), vsound(idim), vsmax
      common /cgas / gamma
      common /numat/ nmat, ntype(nmatmax)
c
      double precision cvelv, dpdrho, dpdi, pressc, rhoi2
c
      data prlim/-1.e15/, tiny/1.e-30/
c
      vsmax=0.
c
      do ipart=1,npart
         rhoi=rho(ipart)
         rhoi2=rhoi*rhoi
         ui=u(ipart)
         imat=matter(ipart)
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
         vsoundi=max(cmin,vsoundi)
         pri=max(pri,prlim)
c
         pr(ipart)=pri
         vsound(ipart)=sqrt(vsoundi+tiny)
         vsmax=max(vsmax,vsound(ipart))
      enddo
c
      return
      end

      subroutine eospg (npart,rho,u)
c************************************************************
c                                                           *
c  This subroutine computes the pressure and sound speed    *
c  for all particles on a list assuming a perfect gas       *
c  equation of state.                                       *
c                                                           *
c************************************************************
c
      parameter (idim=300000)
c
      dimension rho(1), u(1)
c
      common /eosq / pr(idim), vsound(idim), vsmax
      common /carai/ matter(idim), itype(idim)
      common /carar/ pmass(idim)
      common /cgas / gamma
c
c--initialize quantities
c
      vsmax=0.
      gama1=gamma - 1.
c
c--isothermal equation of state
c
      if(gamma.eq.1.)then
         do ipart=1,npart
            pr(ipart)=u(ipart)*rho(ipart)
            vsound(ipart)=sqrt(pr(ipart)/rho(ipart))
            vsmax=max(vsmax,vsound(ipart))
         enddo
         return
      end if
c
c--perfect gas equation of state
c
      do ipart=1,npart
         pr(ipart)=u(ipart)*gama1*rho(ipart)
         vsound(ipart)=sqrt(gamma*pr(ipart)/rho(ipart))
         vsmax=max(vsound(ipart),vsmax)
      enddo
      return
c
      end
      subroutine evol
c************************************************************
c                                                           *
c  this routine drives the time evolution of the system.    *
c                                                           *
c************************************************************
c
c--set all quantities needed for the run
c
      call preset
c
c--write all quantities on listings
c
      call header
c
c--evolve the system
c  -----------------
c
      istep=0
      do i=1,1000000
         call integs(istep)
      enddo
c
      return
      end
      subroutine findneigh(ipart,xi,yi,zi,rcut,nlist)
c************************************************************
c                                                           *
c  This subroutine provides a list of particles within a    *
c  radius 2h of a particle.                                 *
c  This subroutine uses the linked list algorithm.          *
c                                                           *
c************************************************************
c
      parameter (idim=300000)
      parameter (ncell=83)
      parameter (ilist=idim)
c
      common /neigh/ neilist(ilist), ilen(idim)
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
      irankhi=irankh(ipart)
      nlist=0
c
c--find cells to search for neighbors
c  ---------------------------------
c
c  1) X axis
c  
      ratioi=float(irankx(ipart)-1)/float(nslice)
      ic=int(ratioi)
c
c  a) minimum x index
c
      xmin=xi-rcut
      icxmin=ic
    1 continue
      if(xmin.le.xmapl(icxmin))then
         if(icxmin.eq.0)go to 2
         icxmin=icxmin - 1
         go to 1
      end if
c
c  b) maximum x index
c
    2 continue
      xmax=xi+rcut
      icxmax=ic
    3 continue
      if(xmax.ge.xmapu(icxmax))then
         if(icxmax.eq.nx)go to 4
         icxmax=icxmax + 1
         go to 3
      end if
c
c--2) Y axis
c
    4 continue
      ratioj=float(iranky(ipart)-1)/float(nslice)
      jc=int(ratioj)
c
c  a) minimum y index
c
      ymin=yi-rcut
      icymin=jc
    5 continue
      if(ymin.le.ymapl(icymin))then
         if(icymin.eq.0)go to 6
         icymin=icymin - 1
         go to 5
      end if
c
c  b) maximum y index
c
    6 continue
      ymax=yi+rcut
      icymax=jc
    7 continue
      if(ymax.ge.ymapu(icymax))then
         if(icymax.eq.ny)go to 8
         icymax=icymax + 1
         go to 7
      end if
c
c--Z axis
c
    8 continue
      ratiok=float(irankz(ipart)-1)/float(nslice)
      kc=int(ratiok)
c
c  e) minimum z index
c
      zmin=zi-rcut
      iczmin=kc
    9 continue
      if(zmin.le.zmapl(iczmin))then
         if(iczmin.eq.0)go to 10
         iczmin=iczmin - 1
         go to 9
      end if
c
c  f) maximum z index
c
   10 continue
      zmax=zi+rcut
      iczmax=kc
   11 continue
      if(zmax.ge.zmapu(iczmax))then
         if(iczmax.eq.nz)go to 12
         iczmax=iczmax + 1
         go to 11
      end if
c
c--find potential neighbors
c  ------------------------
c
   12 continue
      do ic=icxmin,icxmax
         do jc=icymin,icymax
            do kc=iczmin,iczmax
               j=icell(ic,jc,kc)
               do while (j.ne.0)
                  if(irankhi.gt.irankh(j))then
                     nlist=nlist + 1
                     neilist(nlist)=j
                  end if
                  j=ll(j)
               enddo
            enddo
         enddo
      enddo
c
      return
      end
      subroutine flagp(npart)
c*********************************************************
c                                                        *
c  Subroutine to flag particles that will require spring *
c  interaction                                           *
c                                                        *
c*********************************************************
c
      parameter(idim=300000)
c
      common /bodyi/ n1, n2
      common /carai/ matter(idim), itype(idim)
      common /reduc/ vonmises(idim), reduce(idim)
      common /typef/ istrength, ifrac, ieos, iker
      common /scode/ istr, i1d0, i1d1, i2d0, i2d1
c
      data tiny/1.e-15/
c
c--set code
c
      i1d0=10
      i1d1=11
      i2d0=-10
      i2d1=-11 
c
c--if not solid, dont flag
c
      if(istrength.eq.0)then
         istr=0
         do ipart=1,npart
            if(ipart.lt.n1)then
               itype(ipart)=i1d0
            else
               itype(ipart)=i2d0
            end if
         enddo
         return
      end if
c
c--flag particles
c
      istr=i1d0*i1d0
      do ipart=1,npart
         redi=reduce(ipart)
         if(ipart.le.n1)then
            itype(ipart)=i1d0
            if(redi.lt.tiny)itype(ipart)=i1d1
         else
            itype(ipart)=i2d0
            if(redi.lt.tiny)itype(ipart)=i2d1
         end if
      enddo
c
      return
      end
      subroutine forces(npart,iupdat,x,y,z,vx,vy,vz,
     1                  rho,h,sxx,syy,sxy,sxz,syz,fx,fy,fz)
c*********************************************************
c                                                        *
c  This subroutine computes the force on the particles   *
c  that need to have their forces evaluated.             *
c  These particles (np) are stored in the array iupdat.  *
c                                                        *
c*********************************************************
c
      parameter (idim=300000)
      parameter (ilist=idim)
      parameter (itable=40010)
c
      dimension iupdat(idim)
      dimension x(idim),y(idim),z(idim),vx(idim),vy(idim),
     1          vz(idim),rho(idim),h(idim)
      dimension sxx(idim), syy(idim), sxy(idim), sxz(idim),
     1          syz(idim), u(idim)
      dimension fx(idim), fy(idim), fz(idim)
c
      common /neigh/ neilist(ilist), ilen(idim)
      common /reduc/ vonmises(idim), reduce(idim)
      common /table/ wij(0:itable), grwij(0:itable), dvtable,
     1               v2max
      common /srate/ epsxx(idim), epsyy(idim), epszz(idim),
     1               epsxy(idim), epsxz(idim), epsyz(idim)
      common /avis / xmum(idim)
      common /rrate/ rxy(idim), rxz(idim), ryz(idim)
      common /divve/ divv(idim), divmax1, divmax2
      common /ener1/ dq(idim)
      common /carai/ matter(idim), itype(idim)
      common /carar/ pmass(idim)
      common /eosq / pr(idim), vsound(idim), vsmax
      common /numpa/ avalpha, avbeta
      common /typef/ istrength, ifrac, ieos, iker
      common /scode/ istr, i1d0, i1d1, i2d0, i2d1
c
      data epsil/1.e-2/
c
c--compute forces
c  --------------
c
      nbnei=0
      nbneimax=0
      vmax=sqrt(v2max)
      do i=1,npart
         ipart=iupdat(i)
         xi=x(ipart)
         yi=y(ipart)
         zi=z(ipart)
         vxi=vx(ipart)
         vyi=vy(ipart)
         vzi=vz(ipart)
c
         pmassi=pmass(ipart)
         rhoi=rho(ipart)
         pri=pr(ipart)
         vsoundi=vsound(ipart)
         hi=h(ipart)
         itypei=itype(ipart)
         redi=reduce(ipart)
c
         sxxi=redi*sxx(ipart)
         syyi=redi*syy(ipart)
         szzi=-sxxi - syyi
         sxyi=redi*sxy(ipart)
         sxzi=redi*sxz(ipart)
         syzi=redi*syz(ipart)
c
c--initialize auxilliary quantities
c
         fxi=0.
         fyi=0.
         fzi=0.
         dqi=0.
         divvi=0.
         epsxxi=0.
         epsyyi=0.
         epszzi=0.
         epsxyi=0.
         epsxzi=0.
         epsyzi=0.
         rxyi=0.
         rxzi=0.
         ryzi=0.
         ileni=0
         xmumaxi=0.
c
c--find interacting neighbors
c
         rcut=vmax*hi 
         call findneigh(ipart,xi,yi,zi,rcut,nlist)
c
c--loop over neighbors
c
         do k=1,nlist
            j=neilist(k)
c
c--define mean h
c
            hmean=0.5*(hi+h(j))
            hmean21=1./(hmean*hmean)
c
c--check if real neighbors
c
            dx=xi - x(j) 
            dy=yi - y(j)
            dz=zi - z(j)
            rij2=dx*dx + dy*dy + dz*dz 
            v2=rij2*hmean21
            if(v2.ge.v2max)go to 50
c
c--compute interactions
c  --------------------
c
c  Use pressure and artificial viscosity between damaged and 
c  undamaged particles as well as between particles of a different
c  object.
c
            rhoj=rho(j)
            rhoij=rhoi*rhoj
            itypeij=itypei*itype(j)
            redj=reduce(j)
c
            hmean41=hmean21*hmean21
            index=v2/dvtable
            dxx=v2-index*dvtable
            index1=index + 1
            dgrwdx=(grwij(index1)-grwij(index))/dvtable
            grwtij=(grwij(index) + dgrwdx*dxx)*hmean41/hmean
            grpmj=pmass(j)*grwtij
            grpmrj=grpmj/rhoj
            grpmi=pmassi*grwtij
            grpmri=grpmi/rhoi
c
            dvx=vxi - vx(j) 
            dvy=vyi - vy(j)
            dvz=vzi - vz(j)
            exx=dvx*dx
            eyy=dvy*dy
            ezz=dvz*dz
            projv=exx + eyy + ezz
c
c  1) artificial viscosity 
c
            tij=0.
            if(projv.lt.0.)then
               vsbar=0.5*(vsoundi + vsound(j))
               robar=0.5*(rhoi+rhoj)
               f=projv*hmean/(rij2+epsil/hmean21)
               xmumaxi=max(xmumaxi,abs(f))
               xmum(j)=max(xmum(j),xmumaxi)
               tij=(avbeta*f*f - f*avalpha*vsbar)/robar
            end if
c
c  2) pressure
c
            diag=(pri+pr(j))/rhoij + tij
            pj=grpmj*diag
            pi=grpmi*diag
            fxi=fxi - pj*dx
            fyi=fyi - pj*dy
            fzi=fzi - pj*dz
            fx(j)=fx(j) + pi*dx
            fy(j)=fy(j) + pi*dy
            fz(j)=fz(j) + pi*dz
c
c  3) divergence velocity
c
            divvi=divvi - projv*grpmrj
            divv(j)=divv(j) - projv*grpmri
c
c  4) internal energy (pressure + artificial viscosity)
c
            dqi=dqi + pj*projv
            dq(j)=dq(j) + pi*projv
c
c  5) material strength if pairs belong to same object
c
            if(itypeij.eq.istr)then
               sigxxij=(sxxi + redj*sxx(j))/rhoij
               sigyyij=(syyi + redj*syy(j))/rhoij
               sigzzij=(szzi - redj*sxx(j) - redj*syy(j))/rhoij
               sigxyij=(sxyi + redj*sxy(j))/rhoij
               sigxzij=(sxzi + redj*sxz(j))/rhoij
               sigyzij=(syzi + redj*syz(j))/rhoij
               tx=sigxxij*dx + sigxyij*dy + sigxzij*dz
               ty=sigxyij*dx + sigyyij*dy + sigyzij*dz
               tz=sigxzij*dx + sigyzij*dy + sigzzij*dz
               fxi=fxi + grpmj*tx
               fyi=fyi + grpmj*ty
               fzi=fzi + grpmj*tz
               fx(j)=fx(j) - grpmi*tx
               fy(j)=fy(j) - grpmi*ty
               fz(j)=fz(j) - grpmi*tz
c
c  6) elastic energy
c
               pscal=(tx*dvx + ty*dvy + dz*dvz) 
               dqi=dqi - pscal*grpmj
               dq(j)=dq(j) - pscal*grpmi
c
c  5) strain rate tensor 
c
               epsxxi=epsxxi - grpmrj*exx
               epsyyi=epsyyi - grpmrj*eyy
               epszzi=epszzi - grpmrj*ezz
               epsxx(j)=epsxx(j) - grpmri*exx
               epsyy(j)=epsyy(j) - grpmri*eyy
               epszz(j)=epszz(j) - grpmri*ezz
               exy=0.5*(dvx*dy+dvy*dx)
               exz=0.5*(dvx*dz+dvz*dx)
               eyz=0.5*(dvy*dz+dvz*dy)
               epsxyi=epsxyi - grpmrj*exy
               epsxzi=epsxzi - grpmrj*exz
               epsyzi=epsyzi - grpmrj*eyz
               epsxy(j)=epsxy(j) - grpmri*exy
               epsxz(j)=epsxz(j) - grpmri*exz
               epsyz(j)=epsyz(j) - grpmri*eyz
c
c  6) rotation rate tensor
c
               rrxy=0.5*(dvx*dy-dvy*dx)
               rrxz=0.5*(dvx*dz-dvz*dx)
               rryz=0.5*(dvy*dz-dvz*dy)
               rxyi=rxyi - grpmrj*rrxy
               rxzi=rxzi - grpmrj*rrxz
               ryzi=ryzi - grpmrj*rryz
               rxy(j)=rxy(j) - grpmri*rrxy
               rxz(j)=rxz(j) - grpmri*rrxz
               ryz(j)=ryz(j) - grpmri*rryz
            end if
c
c--count neighbors
c
            ileni=ileni + 1
            ilen(j)=ilen(j) + 1
c
   50       continue
         enddo
c
c--store final quantities
c
         fx(ipart)=fx(ipart) + fxi
         fy(ipart)=fy(ipart) + fyi
         fz(ipart)=fz(ipart) + fzi
         dq(ipart)=dq(ipart) + dqi
         divv(ipart)=divv(ipart) + divvi
c
         epsxx(ipart)=epsxx(ipart) + epsxxi
         epsyy(ipart)=epsyy(ipart) + epsyyi
         epszz(ipart)=epszz(ipart) + epszzi
         epsxy(ipart)=epsxy(ipart) + epsxyi
         epsxz(ipart)=epsxz(ipart) + epsxzi
         epsyz(ipart)=epsyz(ipart) + epsyzi
c
         rxy(ipart)=rxy(ipart) + rxyi
         rxz(ipart)=rxz(ipart) + rxzi
         ryz(ipart)=ryz(ipart) + ryzi
c
         ilen(ipart)=ilen(ipart) + ileni
c
         xmum(ipart)=max(xmum(ipart),xmumaxi)
c
      enddo
c
      return
      end
      subroutine fracdata
c************************************************
c                                               *
c  read fracture model quantities               *
c                                               *
c************************************************
c
      parameter (idim=300000)
      parameter (nmatmax=21)
c
      double precision cweib
c
      common /parti/ npart
      common /typef/ istrength, ifrac, ieos, iker
      common /fract/ acoef(idim), epsmin(idim), grav(idim),
     1               xm(idim), nflaws(idim), young(idim)
      common /mohrc/ sigma0(idim), sigm0(idim), friction(idim), 
     1               fkmin(idim), fkmax(idim)
c
c--read brittle fracture model variables and constants
c  ---------------------------------------------------
c
      if(ifrac.ge.1) then
         open(30,file='cracks')
         do i=1,nmatmax
            read(30,*)dum
         enddo
         do i=1,npart
            read(30,*)epsmin(i), xm(i), nflaws(i), 
     1                acoef(i), young(i), grav(i)
         enddo
         close(30)
      end if
c
c--read Mohr-Coulomb related fracture constants
c  --------------------------------------------
c
      if(ifrac.eq.2.or.ifrac.eq.3)then
         open(30,file='slips')
         do i=1,20
            read(30,*)j,sigm0(i),fkmin(i),fkmax(i)
         enddo
         do i=1,npart
            read(30,*)sigma0(i),friction(i)
         enddo
         close (30)
      end if
c
      return
      end
      subroutine fracture(npart,iupdat,dm,sxx,syy,sxy,sxz,
     1                    syz,ddm)
c************************************************************
c                                                           *
c  Subroutine computing crack growth.                       *
c                                                           *
c  ifrac=0   : no fracture                                  *
c  ifrac=1   : tensile failure only using Weibull flaws     *
c  ifrac=2   : shear failure using Mohr-Coulomb             *
c  ifrac=3   " both shear and tensile filure                *
c                                                           *              
c************************************************************
c
      parameter(idim=300000)
c
      dimension sxx(idim), syy(idim), sxy(idim), sxz(idim), 
     1          syz(idim), dm(idim), ddm(idim), iupdat(idim)
c
      common /fract/ acoef(idim), epsmin(idim), grav(idim),
     1               xm(idim), nflaws(idim), young(idim)
      common /mohrc/ sigma0(idim), sigm0(idim), friction(idim), 
     1               fkmin(idim), fkmax(idim)
      common /reduc/ vonmises(idim), reduce(idim)
      common /dmx  / dnflmx(idim)
      common /carai/ matter(idim), itype(idim)
      common /carar/ pmass(idim)
      common /eosq / pr(idim), vsound(idim), vsmax
      common /damag/ ratioft(idim), ratiofs(idim)
      common /bodyi/ n1, n2
      common /bodyr/ fmas1, fmas2
      common /typef/ istrength, ifrac, ieos, iker
c
      data pi2/1.570796327/
      data tiny/1.e-20/
      data epsbig/1.e10/
c
c--no fracture, return 
c  -------------------
c
      if(ifrac.eq.0)then
         do i=1,npart
            ipart=iupdat(i)
            ddm(ipart)=0.
         enddo
         return
      end if
c
c--apply fracture model
c  --------------------
c
      do i=1,npart
         ipart=iupdat(i)
         ddm(ipart)=0.
         ratioft(ipart)=0.
         ratiofs(ipart)=0.
         dnflmx(ipart)=1.0
c
c--if not already totally damaged, and if the particle is 
c  fracturable, compute time derivative of damage
c
         if(dm(ipart).lt.1.0)then
c
c--define total stress
c
            pri=pr(ipart)
            redi=reduce(ipart)
            gi=grav(ipart)
            sxxi=redi*sxx(ipart) 
            syyi=redi*syy(ipart) 
            szzi=-sxxi - syyi 
            sxxi=sxxi - pri - gi
            syyi=syyi - pri - gi
            szzi=szzi - pri - gi
            sxyi=redi*sxy(ipart)
            sxzi=redi*sxz(ipart)
            syzi=redi*syz(ipart)
c
c--make principal axis transformation on the total stress, including
c  self-gravity term
c
            call paxis(sxxi,syyi,szzi,sxyi,sxzi,syzi,
     1                 sig1,sig2,sig3)
c
            stmax=max(sig1,sig2,sig3)
            stmin=min(sig1,sig2,sig3)
c
c--tensile failure (Weibull)
c  -------------------------
c
            if(ifrac.eq.1.or.ifrac.eq.3)then
               youngi=max(young(ipart)*(1.-dm(ipart)**3),tiny)
               strain=stmax/youngi
               epsmini=epsmin(ipart)
               ratioft(ipart)=strain/epsmini
               rati=ratioft(ipart)
               if(rati.gt.1.0.and.(epsmini.lt.epsbig))then
                  xmi=xm(ipart)
                  xfi=nflaws(ipart)
                  pop=min(rati**xmi,xfi)
                  pop3=pop**0.3333333333333
                  fpop=pop/xfi 
                  dnflmx(ipart)=fpop**0.3333333333333
                  ddm(ipart)=pop3*acoef(ipart)
               end if
            end if
c
c--shear failure (Mohr Coulomb)
c  ----------------------------
c
            if(ifrac.eq.2.or.ifrac.eq.3)then
               if(ratioft(ipart).le.1.) then
                  s3=stmax
                  s1=stmin
                  ds31=0.5*(s3 - s1)
                  xki=friction(ipart)
                  sig0i=sigma0(ipart)*redi
                  angle=pi2 - atan(1./xki)
                  si=sin(angle)
                  co=cos(angle)
                  sigman=0.5*(s1+s3) + ds31*co
                  sigmas= ds31*si
                  sigtot=sigmas + xki*sigmas
                  ratiofs(ipart)=sigtot/(sig0i+tiny)
                  if(ratiofs(ipart).ge.1.)then
                     dnflmx(ipart)=1.0
                     ddm(ipart)=acoef(ipart) 
                  end if
               end if
            endif
         end if
      enddo
c
      return
      end
      subroutine hdot (npart,iupdat,h,dh)
c************************************************************
c                                                           *
c  this subroutine computes the derivative of the smoothing *
c  length.                                                  *
c                                                           *
c************************************************************
c
      parameter (idim=300000)
      parameter (ilist=idim)
c
      common /bodyi/ n1, n2
      common /bodyr/ fmas1, fmas2
      common /neigh/ neilist(ilist), ilen(idim)
      common /nlim / neimin, neimax
      common /divve/ divv(idim), divmax1, divmax2
      common /reduc/ vonmises(idim), reduce(idim)
      common /srate/ epsxx(idim), epsyy(idim), epszz(idim),
     1               epsxy(idim), epsxz(idim), epsyz(idim)
      common /typef/ istrength, ifrac, ieos, iker
      common /carai/ matter(idim), itype(idim)
      common /scode/ istr, i1d0, i1d1, i2d0, i2d1
c
      dimension h(idim), dh(idim), iupdat(idim)
c
c--compute derivative of h, try to enforce finite range
c  of neighbors
c
      rangen=10
c
      do i=1,npart
         ipart=iupdat(i)
         if(istrength.ge.1)then
            itypei=itype(ipart)
            if((itypei.eq.i1d0) .or. (itypei.eq.i2d0))then
               divvi=epsxx(ipart)+epsyy(ipart)+epszz(ipart)
            else
               divvi=divv(ipart)
            end if
         else
            divvi=divv(ipart)
         end if
c
c--use formula to compute hdot
c
         dhi=h(ipart)*divvi/3.
c
c--determine neighbor statistics
c
         if(ipart.le.n1)then
            dhs=-divmax1*h(ipart)
         else
            dhs=-divmax2*h(ipart)
         end if
         dhl=-dhs
         neireal=ilen(ipart)
         dnsup=max(neimax-neireal,-50)
         dninf=max(neireal-neimin,-50)
         if(dnsup.lt.rangen) then
c
c--number of neighbors too close to upper limit
c
            wsex1=exp(dnsup/5.)
            wsex2=1./wsex1
            dhi=(wsex1*dhi+wsex2*dhs)/(wsex1+wsex2)
         else if(dninf.lt.rangen) then
c
c--number of neighbors too close to lower limit
c
            wiex1=exp(dninf/5.)
            wiex2=1./wiex1
            dhi=(wiex1*dhi+wiex2*dhl)/(wiex1+wiex2)
         end if
         dh(ipart)=dhi
      enddo
c
      return
      end
      subroutine header 
c************************************************************
c                                                           *
c  this routine writes on first page of listing the value   *
c  of all variable defined at the start of the run.         *
c                                                           *
c************************************************************
c
      parameter (idim=300000)
      parameter (nmatmax=21)
c
      double precision angmom, umass
      double precision cweib
c
      common /parti/ npart
      common /partr/ x(idim),y(idim),z(idim),vx(idim),
     1               vy(idim),vz(idim),u(idim),rho(idim),h(idim)
      common /units/ umass, udist, udens, utime, uergg, uergcc
      common /typef/ istrength, ifrac, ieos, iker
      common /bndr / xmin,xmax,ymin,ymax,zmin,zmax
      common /times/ t, dt
      common /dissi/ damp
      common /numpa/ avalpha, avbeta
      common /nlim / neimin, neimax
      common /bodyi/ n1, n2
      common /bodyr/ fmas1, fmas2
      common /cgas / gamma
      common /numat/ nmat, ntype(nmatmax)
      common /integ/ courant, dtmax
      common /recor/ irec
      common /fracm/ pweib(nmatmax), cweib(nmatmax)
      common /mohrc/ sigma0(idim), sigm0(idim), friction(idim), 
     1               fkmin(idim), fkmax(idim)
      common /logun/ iprint, iterm, idisk1, idisk2, idisk3
      common /files/ file1, file2, file3
c
      character*35 var
      character*7 file1, file2, file3
c
c--write units
c
      angmom=umass*dble(udist**2)/dble(utime)
      velo=udist/utime
      write(iprint,100)umass,udist,udens,utime,velo,uergg,angmom
  100 format(//,' the computations are done in the following units',
     1       /,
     2' units of :  mass       :',1pe12.4,'   distance    :',1pe12.4,/,
     3'             density    :',1pe12.4,'   time        :',1pe12.4,/,
     4'             velocity   :',1pe12.4,'   energy/mass :',1pe12.4,/,
     5'             ang. mom.  :',1pe12.4,///)
c
c--write options
c
      var='specific internal energy'
      write(iprint,101)var
  101 format(' variable of state used : density and ',a35,/)
      write(iprint,102)nmat,(ntype(i),i=1,nmat)
  102 format(' number of different type of material used : ',i2,/,
     1       ' type(s) : ',11(i2,1x))
c
c--fracture constants
c
      if(ifrac.gt.0)then
         write(iprint,115)
  115    format(/,' fracture constants: ')
         do j=1,nmat
            im=ntype(j)
            write(iprint,116)im,pweib(im),cweib(im),sigm0(im),
     1                       fkmin(im),fkmax(im)
  116       format(' material type: ',i2,/,
     1             '    m= ',1pe12.4,'    k= ',1pe12.4,/,
     2             ' sig0= ',1pe12.4,2x,1pe12.4,' < k < ',
     3             1pe12.4)
         enddo
      end if
c
c--mass fractions
c
      write(iprint,103)npart,(fmas1+fmas2),n1,fmas1,n2,fmas2
  103 format(/,
     1' total number of particles used : ',i6,/,
     2' fraction of total mass         :',1pe12.4,/,
     3' distribution : object 1 :  number of particles : ',i6,/, 
     4'                            mass fraction       : ',1pe12.4,/,
     5'                object 2 :  number of particles : ',i6,/,
     6'                            mass fraction       : ',1pe12.4)
c
c--options of the code
c
      write(iprint,105)ieos,istrength,ifrac,damp
  105 format(/,' Options set during this run:',/,
     1         ' ----------------------------',/,
     3' eos type               : ',i2,' mat. strength  : ',i2,/
     4' fracture               : ',i2,' general damp.  : ',f3.1,/)
c
c--numerical constants
c
      write(iprint,109)avalpha,avbeta
  109 format(' Numerical constants used during this run:',/,
     1       ' -----------------------------------------',/,
     2       ' artificial viscosity  alpha   : ',1pe12.3,/,
     3       '                       beta    : ',1pe12.3)
      write(iprint,110)courant,dtmax
  110 format(' integration parameters        :',/,
     1       ' courant number                : ',1pe12.3,/,
     2       ' maximum time step authorized  : ',1pe12.3)
      write(iprint,111)neimin, neimax
  111 format(' algorithm for neighbor search :   linked list '/,
     3          ' limits on number of neighbors :  ',i3,' to ',i3,/)
c
c--kernel information
c
      write(iprint,112)iker
  112 format(' Kernel used for this run: ',i2,/,
     1       ' -------------------------',/)
c
c--write name of file used
c
      write(iprint,114)file1
  114 format(' name of input file : ',a7,//)
c
      call flush(iprint)
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
      subroutine inform 
c************************************************************
c                                                           *
c  this routine computes relevant quantities for print out  *
c                                                           *
c************************************************************
c
      parameter (idim=300000)
c
      common /parti/ npart
      common /partr/ x(idim),y(idim),z(idim),vx(idim),
     1               vy(idim),vz(idim),u(idim),rho(idim),h(idim)
      common /damm / dm(idim)
      common /damag/ ratioft(idim), ratiofs(idim)
      common /tempe/ temp(idim)
      common /typef/ istrength, ifrac, ieos, iker
      common /carai/ matter(idim), itype(idim)
      common /carar/ pmass(idim)
      common /out1 / cmx1, cmy1, cmz1, vcmx1, vcmy1, vcmz1, romean1,
     1               romax1, dmax1, zmax1, tmean1, tmax1, hmi1, hma1,
     2               damax1, damean1
      common /out2 / cmx2, cmy2, cmz2, vcmx2, vcmy2, vcmz2, romean2,
     1               romax2, dmax2, zmax2, tmean2, tmax2, hmi2, hma2,
     2               damax2, damean2
      common /bodyi/ n1, n2
      common /bodyr/ fmas1, fmas2
      common /logun/ iprint, iterm, idisk1, idisk2, idisk3
c
c--compute : center of mass, velocity of cm, mean density and dispersion
c  for first object
c
      cmx1=0.
      cmy1=0.
      cmz1=0.
      vcmx1=0.
      vcmy1=0.
      vcmz1=0.
      cmx2=0.
      cmy2=0.
      cmz2=0.
      vcmx2=0.
      vcmy2=0.
      vcmz2=0.
c
      romean1=0.
      romax1=0.
      romean2=0.
      romax2=0.
      sigro1=0.
      hmi1=1.e30
      hma1=0.
      hmi2=1.e30
      hma2=0.
      damax1=0.
      damean1=0.
      damax2=0.
      damean2=0.
c
      do i=1,n1
         pmassi=pmass(i)
         cmx1=cmx1 + pmassi*x(i)
         cmy1=cmy1 + pmassi*y(i)
         cmz1=cmz1 + pmassi*z(i)
c
         vcmx1=vcmx1 + pmassi*vx(i)
         vcmy1=vcmy1 + pmassi*vy(i)
         vcmz1=vcmz1 + pmassi*vz(i)
c
         hmi1=min(hmi1,h(i))
         hma1=max(hma1,h(i))
c
         romean1=romean1 + rho(i)
         romax1=max(romax1,rho(i))
c
         dmi=dm(i)**3
         damax1=max(damax1,dmi)
         damean1=damean1 + dmi
      enddo
c
      cmx1=cmx1/fmas1
      cmy1=cmy1/fmas1
      cmz1=cmz1/fmas1
c
      vcmx1=vcmx1/fmas1
      vcmy1=vcmy1/fmas1
      vcmz1=vcmz1/fmas1
c
      romean1=romean1/n1
      damean1=damean1/n1
c
c--if available compute mean temperature
c
      if(ieos.ge.2) then
         tmean1=0.
         tmax1=0.
         do i=1,n1
            tmean1=tmean1 + temp(i)
            tmax1=max(tmax1,temp(i))
         enddo
         tmean1=tmean1/n1
      end if
c
c--compute : center of mass, velocity of cm, mean density and dispersion
c  for second object (if existing)
c
      if(n2.eq.0) go to 30
c
      do i=n1+1,npart
         pmassi=pmass(i)
         cmx2=cmx2 + pmassi*x(i)
         cmy2=cmy2 + pmassi*y(i)
         cmz2=cmz2 + pmassi*z(i)
c
         vcmx2=vcmx2 + pmassi*vx(i)
         vcmy2=vcmy2 + pmassi*vy(i)
         vcmz2=vcmz2 + pmassi*vz(i)
c
         hmi2=min(hmi2,h(i))
         hma2=max(hma2,h(i))
c
         romean2=romean2 + rho(i)
         romax2=max(romax2,rho(i))
c
         dmi=dm(i)**3
         damax2=max(damax2,dmi)
         damean2=damean2 + dmi
      enddo
c
      cmx2=cmx2/fmas2
      cmy2=cmy2/fmas2
      cmz2=cmz2/fmas2
c
      vcmx2=vcmx2/fmas2
      vcmy2=vcmy2/fmas2
      vcmz2=vcmz2/fmas2
c
      romean2=romean2/n2
      damean2=damean2/n2
c
c--if available compute mean temperature
c
      if(ieos.ge.2) then
         tmean2=0.
         tmax2=0.
         do i=n1+1,npart
            tmean2=tmean2 + temp(i)
            tmax2=max(tmax2,temp(i))
         enddo
         tmean2=tmean2/n2
      end if
c
c--compute maximum distance
c
   30 continue
      dmax1=0.
      zmax1=0.
      do i=1,n1
         dz=z(i) - cmz1
         d2=(x(i)-cmx1)**2 + (y(i)-cmy1)**2 + dz*dz
         dmax1=max(dmax1,d2)
         zmax1=max(zmax1,abs(dz))
      enddo
      dmax1=sqrt(dmax1)
c
      dmax2=0.
      zmax2=0.
      do i=n1+1,npart
         dz=z(i) - cmz2
         d2=(x(i)-cmx2)**2 + (y(i)-cmy2)**2 + dz*dz
         dmax2=max(dmax2,d2)
         zmax2=max(zmax2,abs(dz))
      enddo
      dmax2=sqrt(dmax2)
c
c--compute energies
c
      call toten
c
c--compute total angular momentum
c
      call angmom
c
c--write dump on disk
c
      iflg=1
      call wdump (iflg)
c
c--write global results on listing
c
      call prout
   60 continue
c
      return
      end
      subroutine integs(istep)
c************************************************************
c                                                           *
c  this subroutine integrates the system of differential    *
c  equations over one dump step.                            *
c                                                           *
c************************************************************
c
      parameter (idim=300000)
c
      common /parti/ npart
      common /partr/ x(idim),y(idim),z(idim),vx(idim),
     1               vy(idim),vz(idim),u(idim),rho(idim),h(idim)
      common /s    / sxx(idim), syy(idim), sxy(idim), sxz(idim), 
     1               syz(idim)
      common /damm / dm(idim)
      common /tpoui/ ntwout
      common /tpour/ twout(100)
      common /typef/ istrength, ifrac, ieos, iker
      common /tming/ tstep, irange
      common /integ/ courant, dtmax
      common /times/ t, dt
      common /logun/ iprint, iterm, idisk1, idisk2, idisk3
c
c--get clock time before starting dump step
c
      call getused (tbefor)
c
c--advance system one dump step
c
      call step (istep,npart,t,dt,x,y,z,vx,vy,vz,u,rho,h,sxx,syy,
     1           sxy,sxz,syz,dm)
c
c--get clock time and compute time needed for this timestep
c
      call getused (tafter)
      tstep=(tafter - tbefor)
c
c--write time used for dump step
c
      write(iprint,100)tstep/60.,dt
  100 format(1x,'time used: ',1pe12.3,' min., dt: ',1pe12.3)
c
c--compute and print out the state of the system
c
      call inform
c
c--check for run termination
c
      if(t.ge.twout(ntwout)) call endrun

      return
      end
      subroutine ktable
c***********************************************************
c                                                          *
c  This subroutine builds a table for the various values   *
c  of the kernel, the gradient of the kernel, the mass     *
c  fraction and the potential energy.                      *
c  The entry is v**2.                                      *
c                                                          *
c***********************************************************
c
      parameter (itable=40010, itab=itable-10)
c
      common /kerne/ cnormk
      common /table/ wij(0:itable), grwij(0:itable), dvtable,
     1               v2max
      common /typef/ istrength, ifrac, ieos, iker
      common /logun/ iprint, iterm, idisk1, idisk2, idisk3
c
      data pi/3.141592654/, tiny/5.e-5/
c
c--build kernel tables
c  -------------------
c
c  a) regular W3 kernel
c
      if(iker.eq.1) then
         v2max=4.00
         dvtable=v2max/itab
         i1=1.0/dvtable
         cnormk=3./(2.*pi)
c
         wij(0)=cnormk*2./3.
         grwij(0)=0.
         do i=1,i1
            v2=i*dvtable
            v=sqrt(v2)
            w=((2.-v)**3 - 4.*(1.-v)**3)/6.
            dw=-0.5*(2.-v)**2 + 2.*(1.-v)**2
            wij(i)=cnormk*w
            grwij(i)=cnormk*dw/v
         enddo
c
         do i=i1+1,itab
            v2=i*dvtable
            v=sqrt(v2)
            w=((2.-v)**3)/6.
            dw=-0.5*(2.-v)**2
            wij(i)=cnormk*w
            grwij(i)=cnormk*dw/v
         enddo
c
         do i=itab+1,itable
            wij(i)=0.
            grwij(i)=0.
         enddo
      end if
c
c  b) new W4 kernel
c
      if(iker.eq.2)then
         v2max=6.25
         dvtable=v2max/itab
         i1=0.25/dvtable
         i2=2.25/dvtable
         cnormk=6./(5.*24.*pi)
c
         wij(0)=cnormk*230./16.
         grwij(0)=0.
         do i=1,i1
            v2=i*dvtable
            v=sqrt(v2)
            w=(2.5-v)**4 - 5.*(1.5-v)**4 + 10.*(0.5-v)**4
            dw=-4.*(2.5-v)**3 + 20.*(1.5-v)**3 - 40.*(0.5-v)**3
            wij(i)=cnormk*w
            grwij(i)=cnormk*dw/v
         enddo
c
         do i=i1+1,i2
            v2=i*dvtable
            v=sqrt(v2)
            w=(2.5-v)**4 - 5.*(1.5-v)**4
            dw=-4.*(2.5-v)**3 + 20.*(1.5-v)**3
            wij(i)=cnormk*w
            grwij(i)=cnormk*dw/v
         enddo
c
         do i=i2+1,itab
            v2=i*dvtable
            v=sqrt(v2)
            w=(2.5-v)**4
            dw=-4.*(2.5-v)**3
            wij(i)=cnormk*w
            grwij(i)=cnormk*dw/v
         enddo
c
         do i=itab+1,itable
            wij(i)=0.
            grwij(i)=0.
         enddo
      end if
c
c--check that kernel is properly normalized
c
      wtot=0.
      npoin=10000
      dv=sqrt(v2max)/(npoin-1)
      do i=1,npoin
         v=(i-1)*dv
         v2=v*v
         index=min(v2/dvtable,float(itab))
         dxx=v2-index*dvtable
         index1=min(index + 1,itab)
         dwdx=(wij(index1)-wij(index))/dvtable
         w=(wij(index) + dwdx*dxx)
         wtot=wtot + v*v*w
      enddo
      wtot=4.*pi*wtot*dv
      error=wtot - 1.
      if(abs(error).gt.tiny)then
         write(*,*)'error in kernel normalization: ',error
         stop
      end if
c
      return
      end
      subroutine lunit 
c************************************************************
c                                                           *
c  this routine attributes logical unit numbers.            *
c                                                           *
c************************************************************
c
      common /logun/ iprint, iterm, idisk1, idisk2, idisk3
c
      iprint=60
      iterm=50
      idisk1=11
      idisk2=12
      idisk3=13
c
      return
      end
      subroutine matdata
c*************************************************
c                                                *
c  subroutine reading material data for elastic, * 
c  plastic and fracture model                    *
c                                                *
c*************************************************
c
      parameter (idim=300000)
      parameter (nmatmax=21)
c
      double precision cweib
c
      common /shear/ xmu(nmatmax)
      common /yield/ umelt(nmatmax), yie(nmatmax)
      common /fracm/ pweib(nmatmax), cweib(nmatmax)
c
      character*50 dumy
c
c--open data file
c
      open(30,file='matdata.input')
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
      imat=0
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
      imat=0
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
      imat=0
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
      imat=0
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
      imat=0
      do i=1,nline
         imatmin=(i-1)*5 + 1
         imatmax=min(nmat,imatmin+4)
         read(30,*)(cweib(j),j=imatmin,imatmax)
      enddo
      close(30)
c
      return
      end
      subroutine matype
c************************************************************
c                                                           *
c  find type of material present and masses                 *
c                                                           *
c************************************************************
c
      double precision pmi,fm1,fm2
c
      parameter (idim=300000)
      parameter (nmatmax=21)
c
      common /parti/ npart
      common /partr/ x(idim),y(idim),z(idim),vx(idim),
     1               vy(idim),vz(idim),u(idim),rho(idim),h(idim)
      common /numat/ nmat, ntype(nmatmax)
      common /carai/ matter(idim), itype(idim)
      common /carar/ pmass(idim)
      common /bodyi/ n1, n2
      common /bodyr/ fmas1, fmas2
c
c--find how much different type of material there is and what are
c  they
c
      do i=1,21
         ntype(i)=9999
      enddo
c
      nmat=0
      do 30 i=1,npart
         j=matter(i)
           do k=1,11
              if(j.eq.ntype(k)) go to 30
           enddo
      nmat=nmat + 1
      ntype(nmat)=j
   30 continue
c
c--find partial masses
c
      fm1=0.d0
      do i=1,n1
         pmi=pmass(i)
         fm1=fm1 + pmi
      enddo
      fm2=0.d0
      do i=n1+1,npart
         pmi=pmass(i)
         fm2=fm2 + pmi
      enddo
      fmas1=fm1
      fmas2=fm2
c
      return
      end
      subroutine mesop
c******************************************************************
c                                                                 *
c  checks for messages                                            *
c                                                                 *
c******************************************************************
c
      logical iex
c
c--check for termination
c
      inquire(file='stop',exist=iex)
      if(iex) then
         write(*,*)'job stopped!'
         iflg=0
         call wdump(iflg)
         open(43,file='ok')
         write(43,*)'termination ok!'
         close(43)
         stop
      end if
c
c--check for abort
c
      inquire(file='abort',exist=iex)
      if(iex) then
         write(*,*)'job aborted!'
         stop
      end if
c
      return
      end
      subroutine mllist (npart,iupdat,x,y,z,h)
c************************************************************
c                                                           *
c  this subroutine builds the linked list by superposing a  *
c  cubic lattice on the computational box. The linked list  *
c  is obtained using the rank of the particles sorted       *
c  along increasing coordinates.                            *
c                                                           *
c************************************************************
c
      parameter (idim=300000)
      parameter (nocc=1)
      parameter (ncell=83)
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
      dimension x(idim), y(idim), z(idim), h(idim)
      dimension iupdat(idim)
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
c--build updat list
c  ----------------
c
      inew=0
      do i=npart,1,-1
         inew=inew + 1
         ipart=indxh(j)
         iupdat(inew)=i
      enddo
c
      return
      end  
      subroutine options
c************************************************************
c                                                           *
c  this subroutine defines all options desired for the run  *
c                                                           *
c************************************************************
c
      double precision umass
c
      common /typef/ istrength, ifrac, ieos, iker
      common /units/ umass, udist, udens, utime, uergg, uergcc
      common /bndr / xmin,xmax,ymin,ymax,zmin,zmax
      common /recor/ irec
      common /dissi/ damp
      common /tpoui/ ntwout
      common /tpour/ twout(100)
      common /tming/ tstep, irange
      common /integ/ courant, dtmax
      common /times/ t, dt
      common /files/ file1, file2, file3
      common /actio/ namerun
      common /logun/ iprint, iterm, idisk1, idisk2, idisk3
c
      character*20 namerun
      character*7 file1, file2, file3
c
c--open input file
c
      open(iterm,file='insph')
c
c--determine options 
c  -----------------
c
c a) name of run
c
      read(iterm,100)namerun
  100 format(a20)
c
c b) basename of binary input data file and dump number
c
      read(iterm,101) file1
  101 format(a7)
      read(iterm,*)irec
c
c c) units for the code
c
      read(iterm,*)umass, udist
c
c d) read physical options
c
      read(iterm,*)istrength,ifrac,damp,ieos
c
c e) read integration options
c
      read(iterm,*)courant,dtmax,ntwout
      do i=1,ntwout
         read(iterm,*)twout(i)
      enddo
c
      return
      end
      subroutine paxis(sxxi,syyi,szzi,sxyi,sxzi,syzi,sig1,sig2,
     1                 sig3)
c****************************************************************
c                                                               *
c  subroutine performing the principal axis transformation of a *
c  symmetric tensor. The 3 eigenvalues are returned.            *
c                                                               *
c****************************************************************
c
      data pi2/6.283185307/, pi4/12.56637061/ 
c
      third=1./3.
      sig1=0.
      sig2=0.
      sig3=0.
c
c--find temporary norm to avoid overflow
c
      smax=max(abs(sxxi),abs(syyi),abs(szzi),abs(sxyi),abs(sxzi),
     1         abs(syzi))
      if(smax.ne.0.)then
         sxxi=sxxi/smax
         syyi=syyi/smax
         szzi=szzi/smax
         sxyi=sxyi/smax
         sxzi=sxzi/smax
         syzi=syzi/smax
c
c--compute 3 invariants
c
         xi1=sxxi+syyi+szzi
         xi2=-(syyi*szzi+szzi*sxxi+sxxi*syyi) +
     1         syzi*syzi + sxzi*sxzi + sxyi*sxyi
         xi3=sxxi*syyi*szzi+2.*syzi*sxzi*sxyi - 
     1       sxxi*syzi*syzi - syyi*sxzi*sxzi -
     2       szzi*sxyi*sxyi
c
c--find eigenvalue by soving cubic equation
c  y**3 + p*y**2 + q*y + r =0
c  (the stress tensor being symmetric, roots are real)
c   
c
         p=-xi1
         q=-xi2
         r=-xi3
         a=third*(3.*q-p*p)
         b=(2.*p*p*p-9.*p*q+27.*r)/27.
c
c--make sure there are 3 real roots (round offs)
c
         a1=a*a*a/27.
         rroots=0.25*b*b + a1
         if(rroots.lt.0.) then
            t1=2.*sqrt(-third*a)
            p3=third*p
            phi=acos(-0.5*b/sqrt(-a1))
            phi1=third*phi
            phi2=third*(phi+pi2)
            phi3=third*(phi+pi4)
            sig1=(t1*cos(phi1)-p3)*smax
            sig2=(t1*cos(phi2)-p3)*smax
            sig3=(t1*cos(phi3)-p3)*smax
         end if
      end if
c
      return
      end
      subroutine plastic(npart,u,dm,sxx,syy,sxy,sxz,syz)
c************************************************************
c                                                           *
c  subroutine applying von Mises criterion to limit the     *
c  deviatoric stress.                                       *
c                                                           *
c************************************************************
c
      parameter(idim=300000)
      parameter (nmatmax=21)
c
      double precision umass
c
      common /bodyi/ n1, n2
      common /bodyr/ fmas1, fmas2
      common /typef/ istrength, ifrac, ieos, iker
      common /carai/ matter(idim), itype(idim)
      common /yield/ umelt(nmatmax), yie(nmatmax)
      common /reduc/ vonmises(idim), reduce(idim)
      common /units/ umass, udist, udens, utime, uergg, uergcc
c
      dimension sxx(idim), syy(idim), sxy(idim), sxz(idim),
     1          syz(idim), u(idim), dm(idim)
c
      data tiny/1.e-15/, roundoff/1.e-5/
c
c--if no strength, return
c
      if(istrength.eq.0)then
         do ipart=1,npart
            vonmises(ipart)=1.0
         enddo
         return
      endif
c
c--compute von Mises criterion
c  ---------------------------
c
      do ipart=1,npart
         imat=matter(ipart)
         ratio=u(ipart)/umelt(imat)
         if(ratio.lt.roundoff)then
            yst=yie(imat)
         else
            yst=max(yie(imat)*(1.-ratio),0.)
            if(yst.le.tiny)then
               sxx(ipart)=0.
               syy(ipart)=0.
               sxy(ipart)=0.
               sxz(ipart)=0.
               syz(ipart)=0.
               vonmises(ipart)=0.
               go to 10
            end if
         end if 
c
         dami=1. - dm(ipart)**3
         sxxi=(dami*sxx(ipart)+tiny)/yst + tiny
         syyi=(dami*syy(ipart)+tiny)/yst + tiny
         szzi=-sxxi - syyi
         sxyi=(dami*sxy(ipart)+tiny)/yst + tiny
         sxzi=(dami*sxz(ipart)+tiny)/yst + tiny
         syzi=(dami*syz(ipart)+tiny)/yst + tiny
c
c--compute invariant and von Mises criterion
c
         vonmises(ipart)=1.0
         xj2=0.5*(sxxi*sxxi+syyi*syyi+szzi*szzi) +
     1       sxyi*sxyi+sxzi*sxzi+syzi*syzi+tiny
         vonmises(ipart)=min(sqrt(1./(3.*xj2)),1.)
c
         sxx(ipart)=vonmises(ipart)*sxx(ipart)
         syy(ipart)=vonmises(ipart)*syy(ipart)
         sxy(ipart)=vonmises(ipart)*sxy(ipart)
         sxz(ipart)=vonmises(ipart)*sxz(ipart)
         syz(ipart)=vonmises(ipart)*syz(ipart)
c
   10    continue
      enddo
c
      return
      end
      subroutine preset
c************************************************************
c                                                           *
c  This routine is setting up all quantities needed before  *
c  starting a simulation.                                   *
c                                                           *
c************************************************************
c
      common /files/ file1, file2, file3
      common /numpa/ avalpha, avbeta
      common /cgas / gamma
      common /nlim / neimin, neimax
      common /typef/ istrength, ifrac, ieos, iker
      common /logun/ iprint, iterm, idisk1, idisk2, idisk3
c
      character*7 file1,file2,file3
c
c--read initial conditions
c  -----------------------
c
      call rdump
c
c--make an inventory of material types used
c  ----------------------------------------
c
      call matype
c
c--determine material properties
c  -----------------------------
c
      call matdata
c
c--initialize equation of state 
c  ----------------------------
c
c  a) ANEOS eos
c
      if(ieos.eq.2)then
         call aneosinit
      end if
c
c  b) Tillotson eos
c
      if(ieos.eq.1)then
         call tillinit
      endif
c
c--initialize fracture model
c  -------------------------
c
      call fracdata
c
c--set kernel and compute tables
c  -----------------------------
c
      iker=1
      call ktable
c
c--set constant for artifial viscosity
c  -----------------------------------
c
      avbeta=3.0
      avalpha=1.5
c
c--set minimum and maximum limit of neighbors 
c  ------------------------------------------
c
      neimin=25
      neimax=100
c
      return
      end
      subroutine prout
c************************************************************
c                                                           *
c  this routine prints out all interesting quantities at    *
c  the present time                                         *
c                                                           *
c************************************************************
c
      double precision umass
c
      common /times/ t, dt
      common /bodyi/ n1, n2
      common /bodyr/ fmas1, fmas2
      common /ener2/ tkin, tterm
      common /angm / angx, angy, angz, angto
      common /cgas / gamma
      common /recor/ irec
      common /gasfr/ fgas1, fgas2
      common /typef/ istrength, ifrac, ieos, iker
      common /out1 / cmx1, cmy1, cmz1, vcmx1, vcmy1, vcmz1, romean1,
     1               romax1, dmax1, zmax1, tmean1, tmax1, hmi1, hma1,
     2               damax1, damean1 
      common /out2 / cmx2, cmy2, cmz2, vcmx2, vcmy2, vcmz2, romean2,
     1               romax2, dmax2, zmax2, tmean2, tmax2, hmi2, hma2,
     2               damax2, damean2
      common /files/ file1, file2, file3
      common /units/ umass, udist, udens, utime, uergg, uergcc
      common /logun/ iprint, iterm, idisk1, idisk2, idisk3
c
      character*7 what, file1, file2, file3
c
c--write characteristics of system
c  -------------------------------
c
c  a) dump number
c
      write(iprint,100)irec
  100 format(//,' ------> DUMP #: ',i4,'  W R I T T E N ',
     1'  O N    D I S K <------',/)
c
c  b) current time
c
      write(iprint,101)t
  101 format(1x,'TIME  : ',1pe16.10)
c
c  c) energies + total angular momentum
c
      total=tkin + tterm 
      write(iprint,103)total,tkin,tterm,angto
  103 format(' general properties of system : ',/,
     1   ' total energy                       : ',1pe14.5,/,
     2   ' kinetic energy                     : ',1pe14.5,/,
     4   ' internal energy                    : ',1pe14.5,/,
     6   ' total angular momentum             : ',1pe14.5)
c
c  d) object no 1
c
    5 continue
      if(ieos.eq.1.or.ieos.eq.3) then
         write(iprint,107)n1,cmx1,cmy1,cmz1,vcmx1,vcmy1,
     1                    vcmz1,hmi1,hma1,dmax1,zmax1,romean1,
     2                    romax1,damean1,damax1
  107    format(/,' object number 1 (',i6,' particles ) : ',/,
     1 ' center of mass  x :',1pe14.5,'  y :',1pe14.5,'  z :',1pe14.5,/,
     2 ' velocity cm    vx :',1pe14.5,' vy :',1pe14.5,' vz :',1pe14.5,/,
     3 ' smoothing l.  min :',1pe14.5,' max:',1pe14.5,/,
     4 ' max. dist. cm   r :',1pe14.5,'  z :',1pe14.5,/,
     5 ' density      mean :',1pe14.5,' max:',1pe14.5,/,
     6 ' damage       mean :',1pe14.5,' max:',1pe14.5,/)
      else
         write(iprint,108)n1,cmx1,cmy1,cmz1,vcmx1,vcmy1,vcmz1,
     1                    hmi1,hma1,dmax1,zmax1,romean1,romax1,tmean1,
     2                    tmax1,fgas1/fmas1,damean1,damax1
  108    format(/,' object number 1 (',i6,' particles ) : ',/,
     1 ' center of mass  x :',1pe14.5,'  y :',1pe14.5,'  z :',1pe14.5,/,
     2 ' velocity cm    vx :',1pe14.5,' vy :',1pe14.5,' vz :',1pe14.5,/,
     3 ' smoothing l.  min :',1pe14.5,' max:',1pe14.5,/,
     4 ' max. dist. cm   r :',1pe14.5,'  z :',1pe14.5,/,
     5 ' density      mean :',1pe14.5,' max:',1pe14.5,/,
     6 ' mean temperature  :',1pe14.5,6x,'max. temperature :',1pe14.5,/,
     7 ' fraction molten   :',1pe14.5,/,
     8 ' damage       mean :',1pe14.5,' max:',1pe14.5,/)
      end if
c
      if(n2.eq.0) go to 10
c
c  e) object no 2
c
      if(ieos.eq.1.or.ieos.eq.3) then
         write(iprint,109)n2,cmx2,cmy2,cmz2,vcmx2,vcmy2,vcmz2,
     1                    hmi2,hma2,dmax2,zmax2,romean2,romax2,
     2                    damean2,damax2
  109    format(/,' object number 2 (',i6,' particles ) : ',/,
     1 ' center of mass  x :',1pe14.5,'  y :',1pe14.5,'  z :',1pe14.5,/,
     2 ' velocity cm    vx :',1pe14.5,' vy :',1pe14.5,' vz :',1pe14.5,/,
     3 ' smoothing l.  min :',1pe14.5,' max:',1pe14.5,/,
     4 ' max. dist. cm   r :',1pe14.5,'  z :',1pe14.5,/,
     5 ' density      mean :',1pe14.5,' max:',1pe14.5,/,
     6 ' damage       mean :',1pe14.5,' max:',1pe14.5,/)
      else
         write(iprint,110)n2,cmx2,cmy2,cmz2,vcmx2,vcmy2,vcmz2,
     1                    hmi2,hma2,dmax2,zmax2,romean2,romax2,tmean2,
     2                    tmax2,fgas2/fmas2,damean2,damax2
  110    format(/,' object number 2 (',i6,' particles ) : ',/,
     1 ' center of mass  x :',1pe14.5,'  y :',1pe14.5,'  z :',1pe14.5,/,
     2 ' velocity cm    vx :',1pe14.5,' vy :',1pe14.5,' vz :',1pe14.5,/,
     3 ' smoothing l.  min :',1pe14.5,' max:',1pe14.5,/,
     4 ' max. dist. cm   r :',1pe14.5,'  z :',1pe14.5,/,
     5 ' density      mean :',1pe14.5,' max:',1pe14.5,/,
     6 ' mean temperature  :',1pe14.5,6x,'max. temperature :',1pe14.5,/,
     7 ' fraction molten   :',1pe14.5,/,
     7 ' damage       mean :',1pe14.5,' max:',1pe14.5,/)
      end if
c
   10 continue
c
      call flush(iprint)
      return
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
      subroutine reduces(npart,dm)
c************************************************************
c                                                           *
c  subroutine reducing deviatoric stress tensor in case of  *
c  plastic or damaged flows. The pressure is modified only  *
c  if negative. Also, find maximum flow divergence.         *
c                                                           *
c************************************************************
c
      parameter(idim=300000)
c
      dimension dm(idim)
c
      common /bodyi/ n1, n2
      common /bodyr/ fmas1, fmas2
      common /eosq / pr(idim), vsound(idim), vsmax
      common /reduc/ vonmises(idim), reduce(idim)
      common /typef/ istrength, ifrac, ieos, iker
c
c--no strength, no reduction
c
      if(istrength.eq.0)then
         do ipart=1,npart
            reduce(ipart)=1.
         enddo
         return
      end if
c
c--compute reduction factor and reduce pressure if negative
c
      do ipart=1,npart
         dmi=1.-dm(ipart)**3
         reduce(ipart)=dmi
         if(pr(ipart).lt.0.)then
            pr(ipart)=dmi*pr(ipart)
         end if
      enddo
c
c--zero out tension in impactor if required
c
      if(istrength.eq.2)then
         do ipart=1,npart
            pr(ipart)=max(pr(ipart),0.)
         enddo 
      endif
c
      return
      end
      subroutine rdump 
c************************************************************
c                                                           *
c  this routine reads a dump into memory                    *
c                                                           *
c************************************************************
c
      parameter (idim=300000)
c
      common /parti/ npart
      common /partr/ x(idim),y(idim),z(idim),vx(idim),
     1               vy(idim),vz(idim),u(idim),rho(idim),h(idim)
      common /s    / sxx(idim), syy(idim), sxy(idim), sxz(idim), 
     1               syz(idim)
      common /damm / dm(idim)
      common /tempe/ temp(idim)
      common /typef/ istrength, ifrac, ieos, iker
      common /carai/ matter(idim), itype(idim)
      common /carar/ pmass(idim)
      common /cgas / gamma
      common /times/ t, dt
      common /bodyi/ n1, n2
      common /bodyr/ fmas1, fmas2
      common /recor/ irec
      common /ener2/ tkin, tterm
      common /files/ file1, file2, file3
      common /logun/ iprint, iterm, idisk1, idisk2, idisk3
c
      character*10 filename
      character*7 file1,file2,file3
      character*3 number
c
c--read inner particles
c  --------------------
c
      write(number,100)irec
  100 format(i3.3)
      filename=file1 // number
      open(idisk1,file=filename,form='unformatted')
      read(idisk1)npart,n1,n2,
     1     t,gamma,(h(i),i=1,npart),tkin,
     2     tterm,(x(i),i=1,npart),(y(i),i=1,npart),(z(i),i=1,npart),
     3     (vx(i),i=1,npart),(vy(i),i=1,npart),(vz(i),i=1,npart),
     4     (u(i),i=1,npart),(temp(i),i=1,npart),(pmass(i),i=1,npart),
     5     (matter(i),i=1,npart),(rho(i),i=1,npart),(sxx(i),i=1,npart),
     6     (syy(i),i=1,npart),(sxy(i),i=1,npart),(sxz(i),i=1,npart),
     7     (syz(i),i=1,npart),(dm(i),i=1,npart)
      close(idisk1)
c
      return
      end
      subroutine rooten (ipart,rhoi,ui,mati,ti,pi,csi,kpai)
c************************************************************
c                                                           *
c  this subroutine finds the temperature and pressure given *
c  density and specific internal energy. The subroutine     *
c  iterates on temperature until energies are equal. The    *
c  root is found using the Brent method (Numerical Recipes  *
c  p 251)                                                   *
c                                                           *
c************************************************************
c
      implicit double precision (a-h,o-z)
c
      data itmax/30/ , eps/1.d-6/
      data prec/1.d-4/
      data evtemp/11604.8d0/
      data tzero/1.d-6/
c
c--initial braket of temperature (in ev)
c
      t1=ti/5.d0
      t2=5*ti
c
c--this part for isothermal relaxation only
c
c     a=300.0d0/evtemp
c     call aneos (a,rhoi,pi,ei,si,cvi,dpdti,dpdri,fkrosi,
c    1            csi,kpai,mati,fve,fva)
c     ti=a
c     ui=ei
c     return
c
c--make sure of lower boundary, if not return energy for zero
c  temperature
c
      do iter=1,6
         a=t1
         call aneos (a,rhoi,pi,ei,si,cvi,dpdti,dpdri,fkrosi,
     1               csi,kpai,mati,fve,fva)
         fa=ei-ui
         if(fa.lt.0.d0) go to 2
         t1=t1/10.d0 
      enddo
      ui=ei
      ti=tzero
      return
c
c--find upper bound for temperature
c
    2 do iter=1,itmax
         b=t2
         call aneos (b,rhoi,pi,ei,si,cvi,dpdti,dpdri,fkrosi,
     1               csi,kpai,mati,fve,fva)
         fb=ei-ui
         if(fb.gt.0.d0) go to 10
         t2=3.d0*t2
      enddo
      write(*,*)'temperature out of bounds!'
      stop
c
c--root is between braket, start iterations
c
   10 continue
      fc=fb
      do i=1,itmax
         if(fb*fc.gt.0.d0) then
            c=a
            fc=fa
            d=b-a
            e=d
         end if
         if(dabs(fc).lt.dabs(fb)) then
            a=b
            b=c
            c=a
            fa=fb
            fb=fc
            fc=fa
         end if
         tm=0.5d0*(c-b)
         tol1=2.d0*eps*abs(b)
         if(abs(tm).lt.tol1.or.abs(fb/ui).lt.prec) then
            ti=b
            return
         end if
         if(abs(e).ge.tol1.and.abs(fa).gt.abs(fb)) then
            s=fb/fa
            if(a.eq.c) then
               p=2.d0*tm*s
               q=1.d0-s
            else
               q=fa/fc
               r=fb/fc
               p=s*(2.d0*tm*q*(q-r)-(b-a)*(r-1.d0))
               q=(q-1.d0)*(r-1.d0)*(s-1.d0)
            end if
            if(p.gt.0.d0) q=-q
            p=abs(p)
            if(2.d0*p.lt.dmin1(3.d0*tm*q-abs(tol1*q),abs(e*q))) then
               e=d
               d=p/q
            else
               d=tm
               e=d
            end if
         else
            d=tm
            e=d
         end if
         a=b
         fa=fb
         if(abs(d).gt.tol1) then
            b=b + d
         else
            b=b + sign(tol1,tm)
         end if
         call aneos (b,rhoi,pi,ei,si,cvi,dpdti,dpdri,fkrosi,
     1               csi,kpai,mati,fve,fva)
         fb=ei-ui
      enddo
c
c--did not converge write(iprint,error message
c
      write(*,*)'rooten: T iteration did not converge!'
      stop
      end
      program sph
c************************************************************
c************************************************************
c************************************************************
c******************                    **********************
c******************    S    P    H     **********************
c******************                    **********************
c************************************************************
c************************************************************
c************************************************************
c                                                           *
c  This is a three-dimensional hydro code based on the      *
c  so-called smoothed particles hydrodynamics method.       *
c                                                           *
c************************************************************
c                                                           *
c  W. Benz        Los Alamos National Laboratory 08/31/86   *
c                 Harvard College Observatory    01/05/90   *
c                 Steward Observatory            30/09/91   *
c                                                           *
c************************************************************
c                                                           *
c  general description                                      *
c  ===================                                      *
c                                                           *
c  This code uses the 3D SPH formulation to solve the hydro *
c  conservation equation.                                   *
c                                                           *
c  This version of the code includes 1) gas pressure        *
c                                    2) self-gravity        *
c                                    3) dissipation (shocks)*
c                                                           *
c  material properties are described though the use of an   *
c  appropriate equation of state     1) perfect gas         *
c                                    2) aneos               *
c                                    3) tillotson           *
c                                                           *
c************************************************************
c
c
      common /infor/ version
      common /actio/ namerun
      common /logun/ iprint, iterm, idisk1, idisk2, idisk3
c
      character*40 version
      character*20 namerun
c
      logical iex
c
      data version /'SPH-10.00-UBE-02.01.03'/
c
c--set logical units
c
      call lunit
c
c--open output file, check first if it exists
c
      inquire(file='outsph',exist=iex)
      if(iex) then
         write(*,*)'OUTSPH already exists!'
         stop
      end if
      open(iprint,file='outsph')
c
c--define options of this run
c
      call options
c
c--define code units
c
      call unit
c
c--evolve system
c
      call evol
c
      stop
      end
      subroutine step (istep,npart,t,dt,x,y,z,vx,vy,vz,u,rho,h,
     1                 sxx,syy,sxy,sxz,syz,dm)
c************************************************************
c                                                           *
c  This subroutine integrates the system of equations using *
c  a Runge-Kutta-Fehlberg integrator of second order.       *
c  Particles are allowed to have individual time-steps.     *
c  All particles are synchronized every dtout at which time *
c  the subroutine is exited.                                *
c                                                           *
c************************************************************
c
      parameter (idim=300000)
c
      double precision umass
c
      logical ifirst
      logical flag(idim)
c
      dimension x(idim),y(idim),z(idim),vx(idim),vy(idim),vz(idim),
     1          u(idim),rho(idim),h(idim)
      dimension sxx(idim),syy(idim),sxy(idim),sxz(idim),syz(idim)
      dimension dm(idim)
c
      common /f1   / f1x(idim),f1y(idim),f1z(idim), 
     2               f1u(idim),f1rho(idim),f1h(idim),
     3               d1sxx(idim),d1syy(idim),d1sxy(idim),d1sxz(idim),
     4               d1syz(idim),d1dm(idim)
      common /f2   / f2x(idim),f2y(idim),f2z(idim), 
     2               f2u(idim),f2rho(idim),f2h(idim),
     3               d2sxx(idim),d2syy(idim),d2sxy(idim),d2sxz(idim),
     4               d2syz(idim),d2dm(idim)

      common /dissi/ damp
      common /bodyi/ n1, n2
      common /bodyr/ fmas1, fmas2
      common /eosq / pr(idim), vsound(idim), vsmax
      common /densi/ temp(idim)
      common /typef/ istrength, ifrac, ieos, iker
      common /bndr / xmin,xmax,ymin,ymax,zmin,zmax
      common /tpoui/ ntwout
      common /tpour/ twout(100)
      common /integ/ courant, dtmax
      common /reduc/ vonmises(idim), reduce(idim)
      common /dmx  / dnflmx(idim)
      common /srate/ epsxx(idim), epsyy(idim), epszz(idim),
     1               epsxy(idim), epsxz(idim), epsyz(idim)
      common /units/ umass, udist, udens, utime, uergg, uergcc
      common /minis/ umincgs, rhomincgs, smincgs, dmmin, hmin
      common /logun/ iprint, iterm, idisk1, idisk2, idisk3
c
      logical iex
c
      data idump/20/
      data ifirst/.true./
      data tiny/1.e-15/
c
c--get minimum values in code units
c
      umincgs=1.0e5
      rhomincgs=0.05
      hmin=0.001
      smincgs=1.e8
      dmmin=0.03
      umin=umincgs/uergg
      rhomin=rhomincgs/udens
      smin=smincgs/uergcc
c
c--Compute next dump time and initialise variables.
c
      do i=1,ntwout
         write(*,*)ntwout,twout(i)
         if(t.lt.twout(i))then
            tnext=twout(i)
            go to 1
         end if
      enddo 
    1 continue
      write(*,*)'t and tnext: ',t,tnext
c
c--define coefficients for integrator
c
      biga=1./3.
      bigb=0.5
      bigc=0.5
c
c--the first time, flag particles, set dt and compute forces
c
      if(ifirst) then
         ifirst=.false.
         dt=dtmax/32.
         call deriv(t,npart,istep,x,y,z,vx,vy,vz,u,rho,h,sxx,
     1              syy,sxy,sxz,syz,dm,f1x,f1y,f1z,f1u,
     2              f1rho,f1h,d1sxx,d1syy,d1sxy,d1sxz,d1syz,d1dm)
         call tstep(npart,t,dt,h,u,rho,sxx,syy,sxy,sxz,syz,dm,f1x,
     1              f1y,f1z,f1u,f1h,f1rho,d1sxx,d1syy,d1sxy,d1sxz,
     2              d1syz,d1dm)

      end if
c
c--integration
c  -----------
c
      idt=0 
   10 continue
      istep=istep + 1
      t=t + dt
      idt=idt + 1
c
c--use predictor to obtain temporary values at end of time step
c  ------------------------------------------------------------
c
      dt2=0.5*dt**2
      do i=1,npart
         x(i)=x(i) + vx(i)*dt + dt2*f1x(i)
         y(i)=y(i) + vy(i)*dt + dt2*f1y(i)
         z(i)=z(i) + vz(i)*dt + dt2*f1z(i)
         vx(i)=vx(i) + dt*f1x(i)
         vy(i)=vy(i) + dt*f1y(i)
         vz(i)=vz(i) + dt*f1z(i)
         rho(i)=rho(i) + dt*f1rho(i)
         rho(i)=max(rho(i),rhomin)
         u(i)=u(i) + dt*f1u(i)
         u(i)=max(u(i),umin)
         h(i)=h(i) + dt*f1h(i)
         h(i)=max(h(i),hmin)
         sxx(i)=sxx(i) + dt*d1sxx(i)
         syy(i)=syy(i) + dt*d1syy(i)
         sxy(i)=sxy(i) + dt*d1sxy(i)
         sxz(i)=sxz(i) + dt*d1sxz(i)
         syz(i)=syz(i) + dt*d1syz(i)
         dm(i)=dm(i) + dt*d1dm(i)
         dm(i)=min(dm(i),1.0)
      enddo

c
c--get forces at end of timestep using predicted values
c
      call deriv(t,npart,istep,x,y,z,vx,vy,vz,u,rho,h,sxx,
     1           syy,sxy,sxz,syz,dm,f2x,f2y,f2z,f2u,f2rho,
     2           f2h,d2sxx,d2syy,d2sxy,d2sxz,d2syz,d2dm)
c
c--use corrector to get final values
c
      dt2=0.5*dt**2
      do i=1,npart
         x(i)=x(i) + biga*dt2*(f2x(i)-f1x(i))
         y(i)=y(i) + biga*dt2*(f2y(i)-f1y(i))
         z(i)=z(i) + biga*dt2*(f2z(i)-f1z(i))
         vx(i)=vx(i) + bigb*dt*(f2x(i)-f1x(i))
         vy(i)=vy(i) + bigb*dt*(f2y(i)-f1y(i))
         vz(i)=vz(i) + bigb*dt*(f2z(i)-f1z(i))
         u(i)=u(i) + bigc*dt*(f2u(i)-f1u(i))
         u(i)=max(u(i),umin)
         rho(i)=rho(i) + bigc*dt*(f2rho(i)-f1rho(i))
         rho(i)=max(rho(i),rhomin)
         h(i)=h(i) + bigc*dt*(f2h(i)-f1h(i))
         h(i)=max(h(i),hmin)
         sxx(i)=sxx(i) + bigc*dt*(d2sxx(i)-d1sxx(i))
         syy(i)=syy(i) + bigc*dt*(d2syy(i)-d1syy(i))
         sxy(i)=sxy(i) + bigc*dt*(d2sxy(i)-d1sxy(i))
         sxz(i)=sxz(i) + bigc*dt*(d2sxz(i)-d1sxz(i))
         syz(i)=syz(i) + bigc*dt*(d2syz(i)-d1syz(i))
         dmoldi=dm(i)
         if(dm(i).lt.1.0) dm(i)=dm(i) + bigc*dt*(d2dm(i)-d1dm(i))
         if(dm(i).ge.1.0) dm(i)=dmoldi
         if(dm(i).ge.1.0) then
            dm(i)=1.0
            sxx(i)=0.
            syy(i)=0.
            sxy(i)=0.
            sxz(i)=0.
            syz(i)=0.
            dm(i)=1.0
         endif
      enddo
c
c--switch acceleration for next time step
c
      do i=1,npart
         f1x(i)=f2x(i)
         f1y(i)=f2y(i)
         f1z(i)=f2z(i)
         f1u(i)=f2u(i)
         f1rho(i)=f2rho(i)
         f1h(i)=f2h(i)
         d1sxx(i)=d2sxx(i)
         d1syy(i)=d2syy(i)
         d1sxy(i)=d2sxy(i)
         d1sxz(i)=d2sxz(i)
         d1syz(i)=d2syz(i)
         d1dm(i)=d2dm(i)
         if(dm(i).ge.1.)then
            sxx(i)=0.
            syy(i)=0.
            sxy(i)=0.
            sxz(i)=0.
            syz(i)=0.
            d1sxx(i)=0.
            d1syy(i)=0.
            d1sxy(i)=0.
            d1sxz(i)=0.
            d1syz(i)=0.
            d1dm(i)=0.
         end if
      enddo
c
c--compute new time step
c
      call tstep(npart,t,dt,h,u,rho,sxx,syy,sxy,sxz,syz,dm,f1x,
     1           f1y,f1z,f1u,f1h,f1rho,d1sxx,d1syy,d1sxy,d1sxz,
     2           d1syz,d1dm)
c
c--check for messages
c
      call mesop   
c
c--write temporary dump
c
      if(mod(idt,idump).eq.0)then
         iflg=0
         call wdump(iflg)
      endif
c
c--check if system has been integrated to dump time,
c  if yes exit, if no take another time step
c
      if(t.lt.tnext)go to 10
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
      parameter (idim=300000)
      parameter (nmatmax=21)
c
      common /till / rozero(nmatmax), smalla(nmatmax), uzero(nmatmax), 
     1               smallb(nmatmax), capa(nmatmax), es(nmatmax), 
     2               esp(nmatmax), alpha(nmatmax), beta(nmatmax), 
     3               capb(nmatmax)
c
      character*50 dumy
c
c--open data file
c
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
      return
      end
      subroutine tstep(npart,t,dt,h,u,rho,sxx,syy,sxy,sxz,syz,dm,fx,
     1                 fy,fz,fu,fh,frho,dsxx,dsyy,dsxy,dsxz,dsyz,ddm)
c************************************************************
c                                                           *
c  This subroutine computes the size of the time step given *
c  all the time derivatives.                                *
c                                                           *
c  1) Description of variables in calling sequence:         *
c  ================================================         *
c  a) INPUT                                                 *
c  --------                                                 *
c  npart   : number of particles                            *
c  t       : current time                                   *
c  h       : smoothing length                               *
c  u       : specific internal energy                       *
c  fx      : total acceleration in x direction              *
c  fy      : total acceleration in y direction              *
c  fz      : total acceleration in z direction              *
c  fu      : time derivative of specific energy             *
c  fh      : time derivative of smoothing length            *
c  frho    : time derivative of density                     *
c  ddm     : time derivative of damage                      *
c                                                           *
c  b) OUTPUT                                                *
c  ---------                                                *
c  dt      : time step                                      *
c                                                           *
c  2) Quantities computed in common blocks:                 *
c  ========================================                 *
c  none                                                     *
c                                                           *
c  3) Subroutine called                                     *
c  ====================                                     *
c  none                                                     *
c                                                           *
c************************************************************
c
      parameter (idim=300000)
c
      double precision umass
c
      common /integ/ courant, dtmax
      common /avis / xmum(idim)
      common /eosq / pr(idim), vsound(idim), vsmax
      common /bodyi/ n1, n2
      common /numpa/ avalpha, avbeta
      common /minis/ umincgs, rhomincgs, smincgs, dmmin, hmin
      common /units/ umass, udist, udens, utime, uergg, uergcc
c
      dimension u(idim), rho(idim), h(idim), sxx(idim), syy(idim), 
     1          sxy(idim), sxz(idim), syz(idim), dm(idim)
      dimension fx(idim), fy(idim), fz(idim), fu(idim), fh(idim),
     1          frho(idim), dsxx(idim), dsyy(idim), dsxy(idim),
     2          dsxz(idim), dsyz(idim), ddm(idim)
c  
      data factor/0.7/
c
      umin=umincgs/uergg
      rhomin=rhomincgs/udens
      smin=smincgs/uergcc
c
c--compute timestep
c
      dt1min=1e30
      dt2min=1e30
      dt3min=1e30
      dt4min=1e30
      dt5min=1e30
      dt6min=1e30
      dt7min=1e30
      dt8min=1e30
      dt9min=1e30
      dt10min=1e30
      dt11min=1e30
      iflag=0
c
      do i=1,npart
         hi=h(i)
c
c--total acceleration condition
c
         atot2=fx(i)**2+fy(i)**2+fz(i)**2
         r1=(hi*hi/atot2)**0.25
         if(r1.lt.dt1min)then
            dt1min=r1
            imin1=i
         end if
c
c--Courant condition + AV
c
         denom=vsound(i) + 1.2*avalpha*(vsound(i)+2.*xmum(i))
         r2=hi/denom
         if(r2.lt.dt2min)then
            dt2min=r2
            imin2=i
         endif
c
c--energy condition
c
         if(u(i).gt.2.*umin)then
            if(fu(i).ne.0.)then
               r3=factor*(u(i)+umin)/abs(fu(i))
               if(r3.lt.dt3min)then
                  dt3min=r3
                  imin3=i
               end if
            end if
         end if
c
c--smoothing length condition
c
         if(h(i).gt.2.*hmin)then
            if(fh(i).ne.0.)then
               r4=factor*h(i)/abs(fh(i))
               if(r4.lt.dt4min)then
                  dt4min=r4
                  imin4=i
               end if
            end if
         end if
c
c--density condition
c
         if(rho(i).gt.2.*rhomin)then
            if(frho(i).ne.0.)then
               r5=factor*(rho(i)+rhomin)/abs(frho(i))
               if(r5.lt.dt5min)then
                  dt5min=r5
                  imin5=i
               end if
            end if
         end if
c
c--damage condition
c
         if(dm(i).gt.2.*dmmin)then
            if(ddm(i).ne.0)then
               r6=factor*(dm(i)+dmmin)/abs(ddm(i))
               if(r6.lt.dt6min)then
                  dt6min=r6
                  imin6=i
               end if
            end if
         end if
c
c--deviatoric stress tensor conditions
c
c  a) sxx
c
         if(abs(sxx(i)).gt.2.*smin)then
            if(dsxx(i).ne.0.)then
               r7=factor*(abs(sxx(i))+smin)/abs(dsxx(i))
               if(r7.lt.dt7min)then
                  dt7min=r7
                  imin7=i
               end if
            end if
         end if
c
c  b) syy
c
         if(abs(syy(i)).gt.2.*smin)then
            if(dsyy(i).ne.0.)then
               r8=factor*(abs(syy(i))+smin)/abs(dsyy(i))
               if(r8.lt.dt7min)then
                  dt8min=r8
                  imin8=i
               end if
            end if
         end if
c
c  c) sxy
c
         if(abs(sxy(i)).gt.2.*smin)then
            if(dsxy(i).ne.0.)then
               r9=factor*(abs(sxy(i))+smin)/abs(dsxy(i))
               if(r9.lt.dt9min)then
                  dt9min=r9
                  imin9=i
               end if
            end if
         end if
c
c  d) sxz
c
         if(abs(sxz(i)).gt.2.*smin)then
            if(dsxz(i).ne.0.)then
               r10=factor*(abs(sxz(i))+smin)/abs(dsxz(i))
               if(r10.lt.dt10min)then
                  dt10min=r10
                  imin10=i
               end if
            end if
         end if
c
c  e) syz
c
         if(abs(syz(i)).gt.2.*smin)then
            if(dsyz(i).ne.0.)then
               r11=factor*(abs(syz(i))+smin)/abs(dsyz(i))
               if(r11.lt.dt11min)then
                  dt11min=r11
                  imin11=i
               end if
            end if
         end if
c
   10    continue
      enddo
c
c--choose final timestep
c
      dtmin=min(dt1min,dt2min,dt3min,dt4min,dt5min,dt6min,
     1          dt7min,dt8min,dt9min,dt10min,dt11min)
      dt=min(courant*dtmin,dtmax)
c
c--set flag
c
      if(dtmin.eq.dt1min)then
         iflag=1
         it=imin1
      end if
      if(dtmin.eq.dt2min)then
         iflag=2
         it=imin2
      end if
      if(dtmin.eq.dt3min)then
         iflag=3
         it=imin3
      end if
      if(dtmin.eq.dt4min)then
         iflag=4
         it=imin4
      end if
      if(dtmin.eq.dt5min)then
         iflag=5
         it=imin5
      end if
      if(dtmin.eq.dt6min)then
         iflag=6
         it=imin6
      end if
      if(dtmin.eq.dt7min)then
         iflag=7
         it=imin7
      end if
      if(dtmin.eq.dt8min)then
         iflag=8
         it=imin8
      end if
      if(dtmin.eq.dt9min)then
         iflag=9
         it=imin9
      end if
      if(dtmin.eq.dt10min)then
         iflag=10
         it=imin10
      end if
      if(dtmin.eq.dt11min)then
         iflag=11
         it=imin11
      end if
c
      if(iflag.eq.1)then
         write(*,*)'t, dt (acceleration): ',t,dt,it
      end if
      if(iflag.eq.2)then
         write(*,*)'t, dt (courant): ',t,dt,it
      end if
      if(iflag.eq.3)then
         write(*,*)'t, dt (energy): ',t,dt,it,u(it),fu(it)
      end if
      if(iflag.eq.4)then
         write(*,*)'t, dt (h): ',t,dt,it
      end if
      if(iflag.eq.5)then
         write(*,*)'t, dt (rho): ',t,dt,it
      end if
      if(iflag.eq.6)then
         write(*,*)'t, dt (dm): ',t,dt,it,dm(it),ddm(it)
      end if
      if(iflag.eq.7)then
         write(*,*)'t, dt (sxx): ',t,dt,it,sxx(it),dsxx(it)
      end if
      if(iflag.eq.8)then
         write(*,*)'t, dt (syy): ',t,dt,it,syy(it),dsyy(it)
      end if
      if(iflag.eq.9)then
         write(*,*)'t, dt (sxy): ',t,dt,it,sxz(it),dsxy(it)
      end if
      if(iflag.eq.10)then
         write(*,*)'t, dt (sxz): ',t,dt,it,sxz(it),dsxz(it)
      end if
      if(iflag.eq.11)then
         write(*,*)'t, dt (syz): ',t,dt,it,syz(it),dsyz(it)
      end if
c
      return
      end
      subroutine toten
c************************************************************
c                                                           *
c  this routine computes all energies per unit total mass   *
c  kinetic, rotational, potential and internal.             *
c                                                           *
c************************************************************
c
      parameter (idim=300000)
c
      common /parti/ npart
      common /partr/ x(idim),y(idim),z(idim),vx(idim),
     1               vy(idim),vz(idim),u(idim),rho(idim),h(idim)
      common /typef/ istrength, ifrac, ieos, iker
      common /out1 / cmx1, cmy1, cmz1, vcmx1, vcmy1, vcmz1, romean1,
     1               romax1, dmax1, zmax1, tmean1, tmax1, hmi1, hma1,
     2               damax1, damean1
      common /out2 / cmx2, cmy2, cmz2, vcmx2, vcmy2, vcmz2, romean2,
     1               romax2, dmax2, zmax2, tmean2, tmax2, hmi2, hma2,
     2               damax2, damean2
      common /tempe/ temp(idim)
      common /carai/ matter(idim), itype(idim)
      common /carar/ pmass(idim)
      common /cgas / gamma
      common /times/ t, dt
      common /ener1/ dq(idim)
      common /bodyi/ n1, n2
      common /bodyr/ fmas1, fmas2
      common /eosq / pr(idim), vsound(idim), vsmax
      common /ener2/ tkin, tterm
      common /logun/ iprint, iterm, idisk1, idisk2, idisk3
c
c--compute first mechanical energies 
c
      tkin=0.
      tterm=0.
c
c--compute energies (kinetic,rotational,potential,thermal) 
c
      do i=1,npart
         vxi=vx(i)
         vyi=vy(i)
         vzi=vz(i)
         vtot2=vxi**2 + vyi**2 + vzi**2
         tkin=tkin + pmass(i)*vtot2
         ttherm=u(i)*pmass(i)
         tterm=tterm + ttherm
      enddo
c
c--normalizations
c
      tkin=0.5*tkin
c
   99 return
      end
      subroutine unit
c************************************************************
c                                                           *
c  this routine computes the transformation between the     *
c  physical units (cgs) and the units used in the code.     *
c                                                           *
c************************************************************
c
      double precision umass
      double precision ud1,ut1,ue1,ueg1,uec1
c
      common /units/ umass, udist, udens, utime, uergg, uergcc
      common /norml/ fnorm
      common /logun/ iprint, iterm, idisk1, idisk2, idisk3
c
      data gg/6.67e-8/
c
c--transformation factor for :
c
c  a) density
c
      ud1=dble(umass)/dble(udist)**3
      udens=ud1
      udens=1.
c
c  b) time
c
      ut1=dsqrt(dble(udist)**3/(dble(gg)*dble(umass)))
      utime=ut1
      utime=1.
c
c  c) ergs
c
      ue1=dble(umass)*dble(udist)**2/dble(utime)**2
      uerg=ue1
      uerg=1.
c
c  c) ergs per gram
c
      ueg1=dble(udist)**2/dble(utime)**2
      uergg=ueg1
      uergg=1.
c
c  d) ergs per cc
c
      uec1=dble(umass)/(dble(udist)*dble(utime)**2)
      uergcc=uec1
      uergcc=1.
c
   99 return
      end
      subroutine wdump (iflag)
c************************************************************
c                                                           *
c  this routine writes a dump on disk                       *
c                                                           *
c************************************************************
c
      parameter (idim=300000)
c
      common /parti/ npart
      common /partr/ x(idim),y(idim),z(idim),vx(idim),
     1               vy(idim),vz(idim),u(idim),rho(idim),h(idim)
      common /s    / sxx(idim), syy(idim), sxy(idim), sxz(idim), 
     1               syz(idim)
      common /damm / dm(idim)
      common /damag/ ratioft(idim), ratiofs(idim)
      common /tempe/ temp(idim)
      common /carai/ matter(idim), itype(idim)
      common /carar/ pmass(idim)
      common /recor/ irec
      common /cgas / gamma
      common /times/ t, dt
      common /bodyi/ n1, n2
      common /bodyr/ fmas1, fmas2
      common /ener2/ tkin, tterm
      common /files/ file1, file2, file3
      common /logun/ iprint, iterm, idisk1, idisk2, idisk3
c
      character*10 filename
      character*7 file1,file2,file3
      character*3 number
c
      if(iflag.eq.1)then
         irec=irec + 1
         write(number,100)irec
  100    format(i3.3)
         filename=file1 // number
      else
         filename='sphtmp'
      end if
c
c--write
c
      open(idisk1,file=filename,form='unformatted')
      write(idisk1)npart,n1,n2,
     1      t,gamma,(h(i),i=1,npart),tkin,
     2      tterm,(x(i),i=1,npart),(y(i),i=1,npart),(z(i),i=1,
     3      npart),(vx(i),i=1,npart),(vy(i),i=1,npart),(vz(i),i=1,
     4      npart),(u(i),i=1,npart),(temp(i),i=1,npart),
     5      (pmass(i),i=1,npart),(matter(i),i=1,npart),(rho(i),
     6      i=1,npart),(sxx(i),i=1,npart),(syy(i),i=1,npart),
     7      (sxy(i),i=1,npart),(sxz(i),i=1,npart),(syz(i),i=1,
     8      npart),(dm(i),i=1,npart),(ratioft(i),i=1,npart),
     9      (ratiofs(i),i=1,npart)
      close(idisk1)
c
      return
      end
      subroutine zeros(npart,rho,sxx,syy,sxy,sxz,syz,fx,fy,fz)
c************************************************************
c                                                           *
c  subroutine to initialize all arrays.                     *
c                                                           *
c************************************************************
c
      parameter (idim=300000)
      parameter (ilist=idim)
c
      common /divve/ divv(idim), divmax1, divmax2
      common /ener1/ dq(idim)
      common /avis / xmum(idim)
      common /neigh/ neilist(ilist), ilen(idim)
      common /rrate/ rxy(idim), rxz(idim), ryz(idim)
      common /srate/ epsxx(idim), epsyy(idim), epszz(idim),
     1               epsxy(idim), epsxz(idim), epsyz(idim)
      common /eosq / pr(idim), vsound(idim), vsmax
      common /reduc/ vonmises(idim), reduce(idim)
c
      dimension fx(idim), fy(idim), fz(idim), rho(idim)
      dimension sxx(idim), syy(idim), sxy(idim), sxz(idim),
     1          syz(idim)
c
      do ipart=1,npart
c
c--set arrays to zero
c
         ilen(ipart)=0
         fx(ipart)=0.
         fy(ipart)=0.
         fz(ipart)=0.
         dq(ipart)=0.
         divv(ipart)=0.
         epsxx(ipart)=0.
         epsyy(ipart)=0.
         epszz(ipart)=0.
         epsxy(ipart)=0.
         epsxz(ipart)=0.
         epsyz(ipart)=0.
         rxy(ipart)=0.
         rxz(ipart)=0.
         ryz(ipart)=0.
         xmum(ipart)=0.
      enddo
c
      return
      end
      subroutine zzmachine
c************************************************************
c                                                           *
c  This subroutine contains all machine dependant coding.   *
c                                                           *
c  1) Description of variables in calling sequence:         *
c  ================================================         *
c  a) INPUT                                                 *
c  --------                                                 *
c  none                                                     *
c                                                           *
c  b) OUTPUT                                                *
c  ---------                                                *
c  none                                                     *
c                                                           *
c  2) Quantities computed in common blocks:                 *
c  ========================================                 *
c  none                                                     *
c                                                           *
c  2) Subroutine called                                     *
c  ====================                                     *
c  GETDAT    : get current calendar date                    *
c  GETIME    : get current clock time                       *
c                                                           *
c************************************************************
c
      dimension iarray(3),tarray(2),istatb(12)
c
c--get date
c
      entry getdat (ij,im,iy)
      call idate (ij,im,iy)
      return
c
c--get time
c
      entry getime(ih,im,is,fhour)
      call itime(iarray)
      ih=iarray(1)
      im=iarray(2)
      is=iarray(3)
      fhour=ih + im/60. + is/3600.
      return
c
c--get time used since begining
c
      entry getused (tused)
      tused=etime(tarray)
      return
c
c--check for file status
c
      entry statfile(file,ifsize)
c     ifile=stat(file,istatb)
c     ifsize=istatb(8)
      return
c
c--get argument on command line
c
c     entry getcom(job)
c     call getarg(1,job)
c     return
c
      end
