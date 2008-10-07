      program binxdr
c************************************************************
c                                                           *
c  subroutine to transform sph binary files to xdr files    *
c                                                           *
c************************************************************
c
      implicit double precision(a-h,o-z)
      parameter (maxvars=15)
      parameter (maxpart=200000)
      real store(15)
c
c      common /xdrvar/ store(maxvars)
      common /sphvar/ x(maxpart),y(maxpart),z(maxpart),vx(maxpart),
     1                vy(maxpart),vz(maxpart),u(maxpart),h(maxpart),
     2                rho(maxpart),pmass(maxpart),matter(maxpart),
     3               temp(maxpart),alpha(maxpart),pr(maxpart)
      
c      common /units/ umass, udist, udens, utime, uergg, uergcc     
c
      character*28 varnam(maxvars)
      character*80 title, outname
      character*3 dumpn
      character*4 suffix
      character*7 sphbasename
      character*11 filename
c
      data varnam /'x','y','z','vx','vy','vz','u','h',    
     1           'pmass','rho','temp','pr','alpha','matter','ipart'/
      data suffix/'.xdr'/
c
c--set variable number
c     
      numvars=15
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
c      write(*,*)'equation of state:  1: Tillotson'
c      write(*,*)'                    2: Aneos'
c      read(*,*)ieos
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
      read(11)npart,nt,np,t,trot,tkin,tgrav,tterm,
     &       (x(i),i=1,npart),
     &       (y(i),i=1,npart),
     &       (z(i),i=1,npart),
     &       (vx(i),i=1,npart),
     &       (vy(i),i=1,npart),
     &       (vz(i),i=1,npart), 
     &       (u(i),i=1,npart),
     &       (h(i),i=1,npart), 
     &       (pmass(i),i=1,npart),
     &       (rho(i),i=1,npart),
     &       (temp(i),i=1,npart),
     &       (pr(i),i=1,npart),
     &       (alpha(i),i=1,npart),
     &       (matter(i),i=1,npart)
       close(11)
       write(*,*) npart
c
c--initialize equation of state 
c  ----------------------------
c
c  a) inventory of material types
c
c      call matype(npart)
c
c  a) Tillotson eos
c
c      if(ieos.eq.1)then
c         call tillinit
c      endif
c
c  b) ANEOS eos
c
c      if(ieos.eq.2)then
c         call aneosinit
c      end if
c
c--open xdr file
c  -------------
c
      call openwritexdrfile(outname,varnam,numvars,npart,t,title)
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
         store(indx) = u(i)
         indx=indx+1
         store(indx) = h(i)
         indx=indx+1
         store(indx) = pmass(i)
         indx=indx+1
         store(indx)= rho(i)
         indx=indx+1
         store(indx) = temp(i)
         indx=indx+1
         store(indx)=pr(i)
         indx=indx+1
         store(indx)=alpha(i)
         indx=indx+1         
         store(indx)=matter(i)
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
      parameter (idim=200000)
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
c
c--
c
c      subroutine unit
c************************************************************
c                                                           *
c  This subroutine computes the unit tranformation factors  *
c  from cgs to units defined by the specification of a mass *
c  (umass) and a distance (udist). Note that the G=1.       *
c                                                           *
c************************************************************
c
c      implicit double precision (a-h,o-z)
c
c-----Given in [gg]=cm**3/g*s**2
c
c      parameter (gg=6.6725985d-8)
c
c      double precision dudist, dutime
c      double precision ud1,ueg1
c
c      common /logun/ iprint, iterm, idisk1, idisk2, idisk3, idisk4
c      common /units/ umass, udist, udens, utime, uergg, uergcc
c
c
c      dudist=udist
c
c--In the next file we write down the transformation from the computer
c  to CGS-UNITS for the SPECIFIC DENSITY, TIME, INNER ENERGY AND
c  PRESSURE.
c
c      open(8,file='computer-cgs', status='old') 
c
c--transformation factor for :
c
c  a) density
c
c      ud1=umass/dudist**3
c      udens=ud1
cc      write(8,*) 'UDENS: ',udens
c
c  b) time
c
c      dutime=sqrt(dudist**3/(gg*umass))
c      utime=dutime
c      write(8,*) 'UTIME: ',utime
      
c
c  c) ergs per gram
c
c      ueg1=dudist**2/dutime**2
c      uergg=ueg1
c      write(8,*) 'UERGG (for energy units transformation): ', uergg
c
c  d) ergs per cubic centimeter
c
c      uec1=umass/(dudist*dutime**2)
c      uergcc=uec1
c      write(8,*) 'UERGCC(for pressure units transformation): ', uergcc
c      return
c
c      end


















