      subroutine unit
c************************************************************
c                                                           *
c  This subroutine computes the unit tranformation factors  *
c  from cgs to units defined by the specification of a mass *
c  (umass) and a distance (udist). Note that the G=1.       *
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
c  common     variables                                     *
c  ------     ---------                                     *
c  units   umass             : unit of mass                 *
c          udist             : unit of distance             *
c          udens             : unit of density              *
c          utime             : unit of time                 *
c          uergg             : unit of ergs/g               *
c                                                           *
c  3) Subroutine called                                     *
c  ====================                                     *
c  none                                                     *
c                                                           *
c************************************************************
c
      implicit double precision (a-h,o-z)
      parameter (gg=6.6725985d-8)
c
      double precision dudist, dutime
      double precision ud1,ueg1
c
      common /logun/ iprint, iterm, idisk1, idisk2, idisk3
      common /units/ umass, udist, udens, utime, uergg, uergcc
c
      udist=3.4171306022889d+8
      umass=7.425d+26
c
      dudist=udist
      open(28,file='transfactors',status='old')
c
c
c--transformation factor for :
c
c  a) density
c
      ud1=umass/dudist**3.d0
      udens=ud1
      write(*,*)'udens:',udens
      write(28,*)'udens:',udens
c
c  b) time
c
      dutime=sqrt(dudist**3.d0/(gg*umass))
      utime=dutime
      write(*,*) 'utime:',utime
      write(28,*) 'utime:',utime
c
c  c) ergs per gram
c
      ueg1=dudist**2.d0 / dutime**2.d0
      uergg=ueg1
      write(*,*) 'uergg:',uergg
      write(28,*) 'uergg:',uergg
c
c  d) ergs per cubic centimeter
c
      uec1=umass/(dudist*dutime**2.d0)
      uergcc=uec1
      write(*,*) 'uergcc:',uergcc
      write(28,*) 'uergcc:',uergcc
c
      return
      end








