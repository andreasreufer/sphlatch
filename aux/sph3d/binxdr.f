      program binxdr
c************************************************************
c                                                           *
c  subroutine to transform sph binary files to xdr files    *
c                                                           *
c************************************************************
c
      parameter (maxvars=22)
      parameter (maxpart=600000)
c
      common /xdrvar/ store(maxvars)
      common /sphvar/ x(maxpart),y(maxpart),z(maxpart),vx(maxpart),
     1                vy(maxpart),vz(maxpart),u(maxpart),h(maxpart),
     2                rho(maxpart),pmass(maxpart), matter(maxpart),
     3                temp(maxpart)
      common /streng/ sxx(maxpart), syy(maxpart), sxz(maxpart), 
     1                sxy(maxpart), syz(maxpart), dm(maxpart), 
     2                ratioft(maxpart), ratiofs(maxpart)
c
      character*28 varnam(maxvars)
      character*80 title, outname
      character*3 dumpn
      character*4 suffix
      character*7 sphbasename
      character*11 filename
c
      data varnam /'x','y','z','vx','vy','vz','h','u','T',
     1             'p','phase','pmass','rho','dm','rt','rs','sxx',
     2             'syy','sxy','sxz','syz','ipart'/
      data suffix/'.xdr'/
c
c--set variable number
c     
      numvars=22
c
c--Make fortran strings 0 terminated like C wants it
c
      do i=1,numvars
         j=index(varnam(i),' ')
         varnam(i)(j:j)=char(0)
      enddo
c
c--read options
c
      write(*,*) 'Enter title-> '
      read(*,fmt='(a)') title
      j=index(title,' ')
      title(j:j)=char(0)
c
      write(*,*)'single file or set of files (1/2)'
      read(*,*)ising
      if(ising.eq.1)then
         write(*,*) 'Enter input SPH filename-> '
         read(*,101) filename
  101    format(a11)
         ndmp1=1
         ndmp2=1
      end if
      if(ising.eq.2)then
         write(*,*) 'Enter input SPH file basename-> '
         read(*,100) sphbasename
  100    format(a7)
         write(*,*) 'Enter start, end dump -> '
         read(*,*) ndmp1, ndmp2
      end if
      write(*,*)'equation of state:  1: Tillotson'
      write(*,*)'                    2: Aneos'
      read(*,*)ieos
c
c--loop over files to translate
c  ----------------------------
c
      do j=ndmp1, ndmp2
         if(ising.eq.2)then
            write(dumpn,110)j
  110       format(i3.3)
            filename=sphbasename(1:7) // dumpn
            jin=index(filename,' ') - 1
            outname=filename(1:jin)//suffix//char(0)
         else
            jin=index(filename,' ') - 1
            outname=filename(1:jin)//suffix//char(0) 
         end if
c
c--read data on sph file
c
      open(unit=11,file=filename,form='unformatted', status='old')
      read(11)npart,n1,n2,
     1      t,gamma,(h(i),i=1,npart),tkin,
     2      tterm,(x(i),i=1,npart),(y(i),i=1,npart),(z(i),i=1,
     3      npart),(vx(i),i=1,npart),(vy(i),i=1,npart),(vz(i),i=1,
     4      npart),(u(i),i=1,npart),(temp(i),i=1,npart),(pmass(i),
     5      i=1,npart),(matter(i),i=1,npart),(rho(i),i=1,npart),
     6      (sxx(i),i=1,npart),(syy(i),i=1,npart),(sxy(i),i=1,npart),
     7      (sxz(i),i=1,npart),(syz(i),i=1,npart),(dm(i),i=1,npart),
     9      (ratioft(i),i=1,npart),(ratiofs(i),i=1,npart)
       close(11)
c
c--initialize equation of state 
c  ----------------------------
c
c  a) inventory of material types
c
      call matype(npart)
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
         call aneosinit
      end if
c
c--open xdr file
c  -------------
c
      call openwritexdrfile( outname, varnam, numvars, npart,
     $                     t, title )
c
c--store quantities one particle at a time
c  ---------------------------------------
c
      do i=1,npart
         indx=1
         store(indx) = x(i)
         indx=indx+1
         store(indx) = y(i)
         indx=indx+1
         store(indx) = z(i)
         indx=indx+1
         store(indx) = vx(i)
         indx=indx+1
         store(indx) = vy(i)
         indx=indx+1
         store(indx) = vz(i)
         indx=indx+1
         store(indx) = h(i)
         indx=indx+1
         store(indx) = u(i)
         indx=indx+1
         store(indx) = temp(i)
         if(ieos.eq.1)then
            imat=matter(i)
            ui=u(i)
            rhoi=rho(i)
            call eost(imat,ui,rhoi,csvi,pi)
            kpai=0
         endif
         if(ieos.eq.2)then
            imat=matter(i)
            tempi=temp(i)
            rhoi=rho(i)
            call eosa(imat,tempi,rhoi,pi,ui,csi,kpai)
         endif
         indx=indx+1
         store(indx) = pi
         indx=indx+1
         store(indx) = kpai
         indx=indx+1
         store(indx) = pmass(i)
         indx=indx+1
         store(indx) = rho(i)
         indx=indx+1
         store(indx) = dm(i)**3
         indx=indx+1
         store(indx) = ratioft(i)
         indx=indx+1
         store(indx) = ratiofs(i)
         indx=indx+1
         store(indx) = sxx(i)
         indx=indx+1
         store(indx) = syy(i)
         indx=indx+1
         store(indx) = sxy(i)
         indx=indx+1
         store(indx) = sxz(i)
         indx=indx+1
         store(indx) = syz(i)
         indx=indx+1
         store(indx)=i
c
c--write particle #i data to xdrfile
c  ---------------------------------
c
         call writeparticle(store, numvars)
      enddo
c
c--close xdr file
c  --------------
c
         call closexdrfile
      enddo
c
      stop
      end
c
      subroutine indexx2(n,arrin,indx)
c**********************************************************
c                                                         *
c  subroutine for sorting see W. Press.                   *
c                                                         *
c**********************************************************
c
      parameter (idim=600000)
c
      dimension arrin(idim),indx(idim)
c
      do 10 i=1,n
      indx(i)=i
   10 continue
      l=n/2 + 1
      ir=n
   15 continue
      if(l.gt.1) then
         l=l-1
         indxt=indx(l)
         q=arrin(indxt)
      else
         indxt=indx(ir)
         q=arrin(indxt)
         indx(ir)=indx(1)
         ir=ir - 1
         if(ir.eq.1) then
            indx(1)=indxt
            return
         end if
      end if
      i=l
      j=l+l
   20 if(j.le.ir) then
         if(j.lt.ir) then
            if(arrin(indx(j)).lt.arrin(indx(j+1)))j=j+1
         end if
         if(q.lt.arrin(indx(j)))then
            indx(i)=indx(j)
            i=j
            j=j+j
         else
            j=ir+1
         end if
         go to 20
      end if
      indx(i)=indxt
      go to 15
      end
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
      close(30)
c
      return
      end
      subroutine eost (imat,ui,rhoi,csv,pri)
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
      common /till / rozero(nmatmax), smalla(nmatmax), uzero(nmatmax), 
     1               smallb(nmatmax), capa(nmatmax), es(nmatmax), 
     2               esp(nmatmax), alpha(nmatmax), beta(nmatmax), 
     3               capb(nmatmax)
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
c    2   basalt         12   alluvium           *
c    3   aluminum       13   anorthosite 1pp    *
c    4   dunite         14   anorthosite hpp    *
c    5   iron 130pt     15   andesite           *
c    6   lucite         16   water              *
c    7   limestone      17   pure ice           *
c    8   halite         18   5% silicate ice    *
c    9   oil shale      19   30% silicate ice   *
c   10   wet tuff       20   special            *
c                                               *
c************************************************
c
      parameter (nmatmax=21)
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
      subroutine matype(npart)
c************************************************************
c                                                           *
c  find type of material present                            *
c                                                           *
c************************************************************
c
      parameter (maxpart=600000)
      parameter (nmatmax=21)
c
      common /sphvar/ x(maxpart),y(maxpart),z(maxpart),vx(maxpart),
     1                vy(maxpart),vz(maxpart),u(maxpart),h(maxpart),
     2                rho(maxpart),pmass(maxpart), matter(maxpart),
     3                temp(maxpart)
      common /numat/ nmat, ntype(nmatmax)
c
c--find how much different type of material 
c
      do i=1,nmatmax
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
      return
      end
