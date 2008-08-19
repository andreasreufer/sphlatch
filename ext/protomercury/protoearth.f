      Program interior 
      implicit double precision (a-h,o-z)
      integer kmax,kount,nok,nvar,nbad,maxmat
      parameter (nvar=3)
      parameter(maxmat=21)
      parameter(MAXLOOP=50)
      integer j,n
      integer loopcount,itdone
      integer errid,nstp
      double precision small,RadHigh,RadLow,diff,RhoUpperLimit
      double precision mc_vs_mtot,tsurf,rhosurf
      double precision ystart(nvar),dxsav
      double precision x1,x2,eps,h1,hmin,xp(200),yp(50,200)
      integer imat,nntype(21),klst,kinp,kpai
      common /path/ kmax,kount,dxsav,xp,yp
      common /fileos/ klst, kinp
      common /eqo/ imat,kpai
      common /ouput/ rhoi,ui,ti,pi 
      common /ins/ rhotry
c
      common /values/ radius,totrad,evtemp,totmass,ppi,uui,tti,zeta
      common /radIteration/ errid,nstp,loopcount,itdone
c
      external derivs,rkqs,aneos2,aneos,errorhand
c
c------------------------------------------------------------
c      Define variables and initial conditions.             !
c------------------------------------------------------------
c                                                           !
c      write(*,*) 'totrad [cm]:'                            !
c      read(*,*) totrad                                     !
      totrad        = 5.5d+8                         !
c                                                           !
c      write(*,*) 'total mass [g]:'                          !
c      read(*,*)   totmass                                   !
      totmass       = 5.33d+27                            !
c                                                           !
c      write(*,*) 'mass of core/total mass [0-1]:'          !
c      read (*,*) mc_vs_mtot                                !
      mc_vs_mtot    = 0.31d0                                !
c                                                           !
c      write(*,*) 'Upper limit for mean density [g/cm^3]:'  !
c      read(*,*) RhoUpperLimit                              !
      RhoUpperLimit = 19.d0                                 !
c                                                           !  
c      write(*,*) 'Lower limit for mean density [g/cm^3]:'  !
c      read(*,*) RhoLowerLimit                              !
      RhoLowerLimit = 1.8d0                                 !
c                                                           !
c      write(*,*) 'Surface temperature [K]:'                !
c      read(*,*) tsurf                                      !
      tsurf         = 1350.d0                               !
c                                                           !
c      write(*,*) 'Surface density [g/cm^3]:'               !
c      read(*,*) rhosurf                                    !
      rhosurf       = 7.8d0                                  !
c                                                           !
c-----------------------------------------------------------!
c
      RadHigh  = (3.d0*totmass/(12.566371d0*RhoLowerLimit))**(1.d0/3.d0)
      RadLow   = (3.d0*totmass/(12.566371d0*RhoUpperLimit))**(1.d0/3.d0)
      write(*,*)'RadHigh,RadLow',RadHigh,RadLow
c      
      if ((totrad.lt.RadLow).or.(totrad.gt.RadHigh)) then
         totrad = (RadHigh*RadLow)**(1.d0/2.d0)
c         write(*,*)'Unlikely radius. New value:',totrad
c         pause 
      endif
c
      loopcount = 0                             ! Number of loops for the radius iteration
      itdone    = 0                             ! itdone indicates radius iteration status (0:running,1:done)
c
      do 1203 loopcount = 1,MAXLOOP             ! Do at most MAXLOOP loops (RADIUS ITERATION)
c
          errid     = 0                         ! Error index (concerning rad.iteration)
          small     = 1.d-5                     ! Desired upper limit for radius(of center to last int.point)/totrad
          kmax      = 200                       ! Maximum number of values to be stored
          x1        = 0.d0                      ! Integrate from planets surface -
          x2        = 1.d0-(mc_vs_mtot)**(2.d0/3.d0)  ! - to core boundary. Int.variable is mu=1-(m/mtot)^(2/3)
          dxsav     = (x2-x1)/199.d0                  
          eps       = 1.d-16                    ! Desired accuracy for Runge-Kutta (see rkqs)
          h1        = 1.d-12                    ! Starting value for stepsize (Runge-Kutta)
          hmin      = 0.0d0                     ! Minimum stepsize
          ystart(1) = log(tsurf)                ! Functions are:  y(1): log(temperature)
          ystart(2) = 1.d0                      !                 y(2): (radius/totrad)^2
          ystart(3) = log(rhosurf)              !                 y(3): log(density)
          kinp      = 21                               
          klst      = 22                                
          imat      = 4                         ! Material defining variable. 4:dunite
          kpai      = 4
          nmat      = 1
          nntype(1) = -imat
c
          open(klst,file='aneosdunite.output')  
          open(kinp,file='aneosdunite.input')        
c
          if (loopcount.eq.1) then              
             call aneos2(1,nmat,0,nntype)       ! Initialize material just in the first loop of rad.integration
          endif
c     
          call odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,  ! Start integration of mantle
     &         derivs,rkqs)                     
c          
          write(*,*) nok,nbad                   ! How many accepted/rejected (by rkqs) steps ?
c
          if (itdone.eq.1) then                 ! If rad.iteration is done, list function values -
             do 99 j = 1,kount                  ! - for data storage points
                write(*,*) xp(j),yp(1,j),yp(2,j),yp(3,j)
 99          continue
          endif
c     
          close(21)                             ! Close the aneosdunite files
          close(22)
c
          write(*,*)'Radii ratio:',(yp(2,kount)/yp(2,1))**(1.d0/2.d0)  ! Ratio of r(core)/totrad
c     
          x1        = 1.d0-(mc_vs_mtot)**(2.d0/3.d0)     ! Now integrate from core boundary -
          x2        = 1.d0                               ! - to center
          dxsav     = (x2-x1)/199.d0
          eps       = 1.d-16                             
          h1        = 1.d-12                             
          hmin      = 0.0d0                              
          ystart(1) = yp(1,kount)               ! starting funct.value: last funct.value of mantle integration
          ystart(2) = yp(2,kount)               
c                                               ! ystart(3): density is found by iteration (see below)
c     
          klst      = 24
          kinp      = 23
          imat      = 5                                ! 5:iron
          kpai      = 4
          nmat      = 1
          nntype(1) = -imat
c     
          open(klst,file='aneosiron.output')
          open(kinp,file='aneosiron.input') 
c     
          if (loopcount.eq.1) then
             call aneos2(1,nmat,0,nntype)       ! Initialize material just in first loop of rad.integration
          endif   
c-------------------------------------------------------------------------------     
c                                                          ! With this LOOP the correnspoding pressure is to be found.
c                                                          ! Density varies and is set as boundary value (after found).
c                                                                              !
          rhoinf = 7.5d0                                   ! Assumed lower limit for density at core boundary          
          rhosup = 13.d0                                   ! Corresponding higher limit
c                                                                              !   
c                                                                              !
          do 18 n = 1,1.d+3                                                    !
c                                                                              !
             rhotry=(rhoinf+rhosup)/2                                          !
c                                                                              !
             call aneos(exp(ystart(1))/evtemp,rhotry,pptry,uui,ssi,ccvi,       !
     &                  ddpdti,ddpdri,ffkrosi,ccsi,kpai,imat,ffve,ffva)        !
c     write(*,*) 'rhotry:',rhotry                                              !
c     write(*,*) 'pptry:',pptry                                                !
c     write(*,*) '(pi-pptry)/pi:',dabs((pi-pptry)/pi)                          !
             if (pptry .lt. pi .or. pptry .lt. 0.d0) then                      !
                rhoinf = rhotry                                                !
                rhosup = rhosup                            ! If assumed pressure is too small, increase density
             endif                                                             !
             if (pptry .gt. pi) then                                           !
                rhoinf = rhoinf                                                !
                rhosup = rhotry                            ! If assumed pressure is too big, decrease density   
             endif                                                             !
             if (dabs((pi-pptry)/pi) .lt. 1.d-8) exit      ! If assumed pressure is close enough, exit                  
 18       continue                                                             !
c                                                                              !
          ystart(3) = log(rhotry)                          ! Set density at core boundary to last tried value
c                                                                              !
c-------------------------------------------------------------------------------   
c
          call odeint(ystart,nvar,x1,x2,eps,h1,hmin,       ! Now that all necessary boundary values are known -
     &                nok,nbad,derivs,rkqs)                ! - start integration of core

c     
c     write(*,*)'ystart', ystart(1),ystart(2),ystart(3)
c     write(*,*)'x1,x2,eps,h1,hmin:', x1,x2,eps,h1,hmin
c
          write(*,*) nok,nbad
c
          if (itdone.eq.1)then                  ! If rad.iteration is done, list function values - 
             do 100 j = 1,kount                 ! - for data storage points
                write(*,*) xp(j),yp(1,j),yp(2,j),yp(3,j) 
 100         continue
          endif
c
          close(23)                             ! Close the aneosiron files
          close(24)
c
          write(*,*) 'radius/totrad:', radius/totrad    ! How close did we get to the center ?
c     
          if (radius/totrad.gt.small) then      ! Not close enough ?
             errid = 4
          endif
c
          if (itdone.eq.1) exit                 ! After the radius iteration succeeded (and the program ran once more) exit
c     
          call errorhand(RadHigh,RadLow,small)  ! Subroutine ERRORHAND adjusts totrad, depending on encountered error
c   
 1203  continue
c      
      if ((loopcount.eq.MAXLOOP).and.(itdone.ne.1)) then     ! When rad.iteration was performed MAXLOOP times -
         write(*,*) 'No radius found'                        ! - and no good enough value for totrad was found -
      endif                                                  ! - exit the program
c
      stop
      end
c
c--Stepper program that takes one "quality-cotrolled"          
c--Runge-Kutta step.                                          
c
      subroutine rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
c*******************************************************************
c USES derivs,rkck                                                 * 
c Fifth-order Runge-Kutta step with monitoring of local truncation *
c error to ensure accuracy and adjust stepsize. Input are the      *
c dependent variable vector y(1:n) and ist derivative dydx(1:n) at *
c the starting value of the independent variable x. Also input are *
c the stepsize to be attempted "htry",the required accuracyeps, and*
c vector yscal(1:n) against which the vector is scaled. On output  *
c y and x are replaced by their new values, "hdid' is the stepsize *
c that was actually acomplished, and "hnext" is the estimated next *
c stepsize. "derivs" is the user-supplied subroutine that computes *
c the right-hand side derivatives.                                 *
c*******************************************************************
c
      implicit double precision(a-h,o-z)
      integer n,nmax
      double precision eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      external derivs
      parameter (nmax=50)
c--USES derivs,rkck
      integer i,errid,MAXSTP,nstp,loopcount,itdone
      double precision errmax,h,xnew,yerr(nmax),ytemp(nmax),SAFETY,
     & PGROW,PSHRNK,ERRCON
      parameter (SAFETY=0.9d0,PGROW=-.2d0,PSHRNK=-.25d0,ERRCON=1.89d-4)
c
      common /values/ radius,totrad,evtemp,totmass,ppi,uui,tti,zeta
      common /radIteration/ errid,nstp,loopcount,itdone
c
c
      h = htry
 1    call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
      errmax = 0.d0
      do 100 i = 1,n
        errmax = max(errmax,abs(yerr(i)/yscal(i)))
 100  continue
      errmax = errmax/eps
      if(errmax.gt.1.d0)then
        h = SAFETY*h*(errmax**PSHRNK)
        if(h.lt.0.1d0*h)then
          h = 0.1d0*h
        endif
        xnew = x+h
c        write(*,*) 'xnew,x,h:',xnew,x,h
c        if(xnew.eq.x)pause 'stepsize underflow in rkqs'
c
        if (xnew.eq.x) then
           errid = 1
           return
        endif
c
        goto 1
      else
        if(errmax.gt.ERRCON)then
          hnext = SAFETY*h*(errmax**PGROW)
        else
          hnext = 5.d0*h
        endif
        hdid = h
        x    = x+h
        do 120 i = 1,n
          y(i) = ytemp(i)
 120     continue
        return
      endif
      end
c
c--The routine where algorithm is implemented and the result
c--on its accuraccy is tested.
c
      subroutine rkck(y,dydx,n,x,h,yout,yerr,derivs)
c********************************************************************
c USES derivs.                                                      *
c Given values for n variables y and their derivatives dydx known   *
c at x, use the fifth-order Cash-Karp Runge-Kutta method to         *
c advance the solution over an interval h and return the incremented* 
c variables as yout. Also return an estimate of the local truncation*
c error in yout using the embedded fourth-order method. The user    *
c supplies the subroutine "derivs(x,y,dydx)", which returns         * 
c derivatives dydx at x.                                            *
c********************************************************************
c
      implicit double precision (a-h,o-z)
      integer n,nmax
      double precision h,x,dydx(n),y(n),yerr(n),yout(n)
      external derivs
      parameter (nmax=50)
      integer i
      double precision ak2(nmax),ak3(nmax),ak4(nmax),ak5(nmax),
     *ak6(nmax),
     *ytemp(nmax),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,B52,B53,
     *B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
      parameter (A2=.2d0,A3=.3d0,A4=.6d0,A5=1.d0,A6=.875d0,
     *B21=0.2d0,B31=3.d0/40.d0,
     *B32=9.d0/40.d0,B41=0.3d0,B42=-0.9d0,B43=1.2d0,B51=-11.d0/54.d0,
     *B52=2.5d0,
     *B53=-70.d0/27.d0,B54=35.d0/27.d0,B61=1631.d0/55296.d0,
     *B62=175.d0/512.d0,
     *B63=575.d0/13824.d0,B64=44275.d0/110592.d0,B65=253.d0/4096.d0,
     *C1=37.d0/378.d0,
     *C3=250.d0/621.d0,C4=125.d0/594.d0,C6=512.d0/1771.d0,
     *DC1=C1-2825.d0/27648.d0,
     *DC3=C3-18575.d0/48384.d0,DC4=C4-13525.d0/55296.d0,
     *DC5=-277.d0/14336.d0,
     *DC6=C6-0.25d0)
      do 11 i = 1,n
        ytemp(i) = y(i)+B21*h*dydx(i)
11    continue
      call derivs(x+A2*h,ytemp,ak2)
      do 12 i = 1,n
        ytemp(i) = y(i)+h*(B31*dydx(i)+B32*ak2(i))
12    continue
      call derivs(x+A3*h,ytemp,ak3)
      do 13 i = 1,n
        ytemp(i) = y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
13    continue
      call derivs(x+A4*h,ytemp,ak4)
      do 14 i = 1,n
        ytemp(i) = y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
14    continue
      call derivs(x+A5*h,ytemp,ak5)
      do 15 i = 1,n
        ytemp(i) = y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+
     *B65*ak5(i))
15    continue
      call derivs(x+A6*h,ytemp,ak6)
      do 16 i = 1,n
        yout(i) = y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
16    continue
      do 17 i = 1,n
        yerr(i) = h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*
     *ak6(i))
17    continue
      return
      end
c
c
c******************************************************************
c Full-fledged "driver" for Runge-Kutta with adaptive step-size   *
c controls.                                                       *
c******************************************************************
c
c
      subroutine odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,
     $                    rkqs)
      implicit double precision (a-h,o-z)
      integer nok,nbad,nvar,KMAXX,MAXSTP,nmax
      double precision eps,h1,hmin,x1,x2,ystart(nvar),TINY
      external derivs,rkqs
      parameter (KMAXX=200,MAXSTP=10000,nmax=50,TINY=1.d-30)
      integer i,kmax,kount,nstp,errid,loopcount,itdone
      double precision dxsav,h,hdid,hnext,x,xsav,dydx(nmax),xp(KMAXX)
      double precision y(nmax),yp(nmax,KMAXX),yscal(nmax)
      double precision radius,totrad,evtemp,totmass,ppi,uui,uuti
      common /eqo/imat,kpai
      common /path/ kmax,kount,dxsav,xp,yp
      common /values/radius,totrad,evtemp,totmass,ppi,uui,tti,zeta
      common /ouput/ rhoi,ui,ti,pi 
      common /radIteration/ errid,nstp,loopcount,itdone
c
c      pause
      x     = x1
      h     = sign(h1,x2-x1)
      nok   = 0
      nbad  = 0
      kount = 0
      do i = 1,nvar
	 y(i) = ystart(i)
	 enddo
         if (kmax.gt.0) xsav = x-2.d0*dxsav
         do nstp = 1,MAXSTP
	    call derivs(x,y,dydx)
	    do i = 1,nvar
	       yscal(i) = abs(y(i))+abs(h*dydx(i))+TINY
	       enddo
               if (kmax.gt.0) then
		 if (abs(x-xsav).gt.abs(dxsav)) then
	           if (kount.lt.kmax-1) then
		       kount = kount+1
		       xp(kount) = x
		       do i = 1,nvar
		          yp(i,kount) = y(i)
		       enddo
		       xsav = x
		   endif
		 endif
	       endif
	   if ((x+h-x2)*(x+h-x1).gt.0.d0) h = x2-x
c           

           call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
c
           if (errid.eq.1) then
              return
           endif
c
           if (itdone.eq.1) then
              open(20,file='protomercurynew', status='old')
              write(20,*)ti*evtemp,radius,rhoi,pi,ui
           endif
c
c           write(*,*) radius/totrad
           if (hdid.eq.h) then
		nok  = nok+1
	   else
		nbad = nbad+1
	   endif
c
	   if ((x-x2)*(x2-x1).ge.0.d0) then 	! are we done?
		do i = 1,nvar
		   ystart(i) = y(i)
		enddo
		if (kmax.ne.0) then
			kount = kount+1
			xp(kount) = x
			do i = 1,nvar
			   yp(i,kount) = y(i)
			enddo
		endif
		return
	   endif
c	   if(dabs(hnext).lt.hmin) pause 'stepsize smaller than minimum'
c
           if(dabs(hnext).lt.hmin) then
              errid = 2
c              hnext = hmin
              return
           endif
c
	   h = hnext
	enddo
        errid = 3
c        pause 'too many steps in odeint'
        return
	end
       
c
c************************************************************************
c Die folgende SUBROUTINE gibt das System der Diff.Gleichungen          *
c wieder.                                                               *
c************************************************************************
c 
      subroutine derivs(x,y,dydx)
      implicit double precision (a-h,o-z)
      double precision x,y(3),dydx(3)
      double precision radius,totmass
      double precision gg,pifour,evtemp,cp,alpha
      double precision tti,rrhoi,ppi,uui,ssi,ccvi,xi
      double precision ddpdti,ddpdri,ffkrosi,ccsi,ffve,ffva
      double precision totrad,plambda,zeta,eta,dxideta,dxidlam
      integer imat,kpai
      external aneos
      common /eqo/imat,kpai
      common /values/ radius,totrad,evtemp,totmass,ppi,uui,tti,zeta
      common /rest/ rrhoi,cp,alpha
      common /rest1/ dpdti,dpdri,dxideta,dxidlam,cvi
      common /ouput/ rhoi,ui,ti,pi
      common/rest2/ eta,xi,plambda
      data evtemp/11604.8d0/
c
c--
c
      data gg/6.6725985d-8/
      data pifour/12.56637d0/
c
      tti    = exp(y(1))/evtemp
      radius = totrad*(y(2))**(1.d0/2.d0)
      rrhoi  = exp(y(3))
      zeta   = (radius/totrad)**2.d0  
c    
      call aneos (tti,rrhoi,ppi,uui,ssi,ccvi,ddpdti,
     &            ddpdri,ffkrosi,ccsi,kpai,imat,ffve,ffva)            
c
      rhoi    = rrhoi
      ti      = tti
      pi      = ppi
      ui      = uui
      cvi     = ccvi
      dpdti   = ddpdti
      dpdri   = ddpdri
      
c            
      cp      = cvi+ti*dpdti*dpdti/(dpdri*rhoi*rhoi)
      alpha   = dpdti/(3.d0*rhoi*dpdri)
c 
      dxideta = dpdti*ti/pi
      dxidlam = dpdri*rhoi/pi
c
      zeta    = (radius/totrad)**2.d0  
      eta     = log(ti)
      xi      = log(pi)
      plambda = log(rhoi)
c
c      write(*,*)'pi',pi
c      write(*,*) 'LAMBDA,ZETA,ETA,XI',plambda,zeta,eta,xi
      dydx(1) = alpha/(cp*exp(plambda))*3.d0*gg*totmass**2.d0/
     &        (2.d0*pifour*totrad**4.d0)*
     &        ((1-x)/zeta)**2.d0
c
c      dydx(1) = 3.d0*gg*totmass**2.d0/(dpdti*exp(eta)*     ! Temp. Gradient for
c     &        2.d0*pifour*totrad**4.d0)*                   ! smaller or gas
c     &        ((1-x)/zeta)**2.d0                           ! planets
c
      dydx(2) = -3.d0*totmass*(1-x)**(1.d0/2.d0) / (pifour*totrad**3.d0*
     &           exp(plambda)*(zeta)**(1.d0/2.d0))
      dydx(3) = 3.d0*gg*totmass*totmass*(1-x)*(1-x)/
     &           (zeta*zeta*exp(xi)*2.d0*pifour*totrad**4.d0*dxidlam)
     &           - dxideta*dydx(1)/dxidlam     
c 
c-- Corresponding to the temp. gradient above this form of diff.eq. for
c-- density is to be used. Of course, under same conditions as for 
c-- temp. gradient.
c     
c      dydx(3) = 3.d0*gg*totmass*totmass*(1-x)*(1-x)/   
c     &        (zeta**2.d0*pifour*totrad**4.d0*exp(plambda)*dpdri) 
c     &        - dpdti*exp(eta)*dydx(1)/(exp(plambda)*dpdri)
      return
      end
c
c
c--------------------------------------------------------------------------
c
      subroutine errorhand(RadHigh,RadLow,small)
c
      implicit double precision (a-h,o-z)
      double precision RadHigh,RadLow,diff,small
      double precision radius,totrad,evtemp,totmass,ppi,uui,tti,zeta
      integer itdone,loopcount,nstp,errid
      character yes,no,reply
c
      common /values/radius,totrad,evtemp,totmass,ppi,uui,tti,zeta
      common /radIteration/ errid,nstp,loopcount,itdone

      diff = dabs(RadHigh-RadLow)
c
      if ((errid.ne.0).and.(itdone.ne.1)) then         
         if ((diff.gt.1.d3).or.(radius/totrad.gt.small)) then
c         
            if (errid.eq.1) then                ! stepsize underflow in rkqs
                RadLow = totrad
                write(*,*) 'RadHigh,RadLow:',RadHigh,RadLow
                write(*,*) 'Stepsize underflow; Next totrad:', totrad
                totrad = (RadLow+RadHigh)/2
                write(*,*) 'totrad,radius,small:',totrad,radius,small
                write(*,*) 'errid:',errid
                pause
             endif      
c
             if (errid.eq.2) then               ! stepsize smaller than hmin
                write(*,*) 'Stepsize too small; continue anyway? (y/n)'
                write(*,*) 'totrad,radius,small:',totrad,radius,small
                write(*,*) 'errid:',errid
                read(*,*) reply
                yes = 'y'
                no  = 'n'
               if(reply.eq.yes) return            
               if(reply.eq.no) then 
                  write(*,*)'hallo no'
               else 
                  write(*,*)'hallo '
               endif
             endif           
c
             if (errid.eq.3) then               ! too many steps in odeint 
                RadLow = totrad
                write(*,*) 'totrad,radius,small:',totrad,radius,small
                write(*,*) 'RadHigh,RadLow:',RadHigh,RadLow
c               totrad=(totrad**3.d0+(radius)**3.d0)**(1.d0/3.d0)            
                totrad=(RadHigh+RadLow)/2
                write(*,*) 'Too many steps; Next totrad:', totrad
                write(*,*) 'errid:',errid
c               pause
             endif
c
             if (errid.eq.4) then               ! too big hole in the center
                RadHigh = totrad
                write(*,*) 'RadHigh,RadLow:',RadHigh,RadLow
                write(*,*) 'totrad,radius,small:',totrad,radius,small
                totrad = max(((RadHigh+RadLow)/2),
     &            dabs((totrad**3.d0-(radius+4.d7)**3.d0))**(1.d0/3.d0))
                write(*,*) 'too big hole; Next totrad:', totrad  
                write(*,*) 'errid:', errid
c               pause
             endif
c
             if (errid.gt.4) then             
                write(*,*) 'Non specified error -> exit'
                pause
                stop
             endif 
c   
             write(*,*)'number of loops:',loopcount
             return  
c
          else
             itdone = 1
             write(*,*) 'RadHigh-RadLow:',diff
             write(*,*) ,itdone
             write(*,*) 'Final radius:  ',totrad
             write(*,*)
             pause
             return
          endif
      endif    
c      
      return
      end
c
c-----------------------------
      include 'aneosequ.f'   !
c-----------------------------
c
c                                                               :
c                         protoearth.f                          :                   aneosequ.f
c                                                               :
c                                       ________________        :
c                                      |                |       :
c        ----------------------------->| 1.main program |<---------------------------------------------
c        |                        ---->|________________|<-------------------------                   |
c        |                        |                             :                 |                   |
c        |                 _______|_______                      :                 |                   |
c        |                |               |                     :          _______|_______     _______|_______
c        |           ---->|  4. odeint    |<----                :         |               |   |               |
c        |           |    |_______________|    |                :         |     aneos2    |   |     aneos     |
c        |           |                         |                :         |_______________|   |_______________|
c        |     ______|_______           _______|_______         :             |   |   |           |   |   | 
c        |    |              |         |               |        :             |   |   |           |   |   |
c        |    |  5.derivs    |    ---->|    2.rkqs     |        :            ... ... ...          |  ... ...
c        |    |______________|    |    |_______________|        :                                 |
c        |           |            |                             :                                 |
c        |           |     _______|_______                      :             [ 5.derivs ] <-------
c        |           |    |               |                     : 
c        |           ---->|    3.rkck     |                     :
c  ______|_______         |_______________|                     :
c |              |                                              :
c |  6.errorhand |                                              :
c |______________|                                              :
c                                                               :    
c                                                               :
c                                                               :
c
c
