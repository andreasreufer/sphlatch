#!/usr/bin/env ipython

import numpy as np
import sys

from numpy import sqrt, sin, arcsin, cos, arccos, \
    log, abs, tan, arctan2, pi, zeros, array

# constants to transform the initial parameters to physical values
Re  = 6.37814e8
Me  = 5.97360e27
lem = 3.40000e41 # Canup (2001) value
#lem = 2.88993e41 # actual (?) value

rad2deg = 360/(2.*pi)

# input: 

class GiantImpact(object):

  def __init__(self, mtar, mimp, Rtar, Rimp, G):
    nan = float('nan')
    self.mtar = mtar
    self.Rtar = Rtar
    
    self.mimp = mimp
    self.Rimp = Rimp

    self.G = G
    self.mtot = mtar + mimp
    self.rimp = Rtar + Rimp
    self.mu = G*self.mtot
    self.vesc = sqrt( 2.*self.mu / self.rimp )

  def getInitVinfL0(self, vinf, L0, relsep):
    rimp = self.rimp
    mimp = self.mimp
    vesc = self.vesc
    
    vimp = sqrt(vinf*vinf + vesc*vesc)
    Lgraz = vimp*mimp*rimp
    impa = arcsin( L0 / Lgraz )

    return self.getInitVimpAlpha(vimp, impa, relsep)


  def getInitVimpAlpha(self, vimp, impa, relsep):
    rimp = self.rimp
    mimp = self.mimp
    mtar = self.mtar
    vesc = self.vesc
    mu   = self.mu
    mtot = self.mtot

    k1    = rimp*(vimp*vimp)/mu
    bscal = sin(impa)
    Lgraz = vimp*mimp*rimp
    L0    = bscal*Lgraz
    betaimp = (pi/2.)-impa

    vinf = sqrt( vimp*vimp - vesc*vesc )
    
    e = sqrt((k1-1.)*(k1-1.)*cos(betaimp)*cos(betaimp) \
        + sin(betaimp)*sin(betaimp))
    thetaimp = arctan2(k1*sin(betaimp)*cos(betaimp), \
        k1*cos(betaimp)*cos(betaimp) - 1.)

    # store as common
    self.e = e
    self.thetaimp = thetaimp

    rperih = rimp *( 1. + e*cos(thetaimp) ) / ( 1. + e )
    vperih = sqrt( vinf*vinf + ( 2.*mu / rperih ) )

    # determine true anomaly theta0 for r0
    r0 = relsep*rimp
    v0 = sqrt( vinf*vinf + ( 2.*mu / r0 ) )
    theta0 = arccos( ( rperih*(1. + e ) - r0 ) / ( e * r0 ) )
    beta0  = arccos( L0 / ( r0*v0*mimp ) )

    if e > 1.:
      a = rperih / ( e - 1. )
      k2 = sqrt( mu / (a*a*a) )
      timp = ( (e*sqrt(e*e-1.)*sin(thetaimp))/(1. + e*cos(thetaimp) ) \
          - log( (sqrt(e*e-1.) + (e - 1.)*tan( thetaimp / 2. ))/      \
          (sqrt(e*e-1.) - (e - 1.)*tan( thetaimp / 2. )) ) ) / k2
  
      t0 = ( (e*sqrt(e*e-1.)*sin(theta0))/(1. + e*cos(theta0) ) \
          - log( (sqrt(e*e-1.) + (e - 1.)*tan( theta0 / 2. ))/  \
          (sqrt(e*e-1.) - (e - 1.)*tan( theta0 / 2. )) ) ) / k2
    else:
      k2 = sqrt( mu / ( 8. * rperih*rperih*rperih) )

      timp = abs(0.5*(tan(thetaimp / 2.) + \
          (1./3.)*pow(tan(thetaimp / 2.),3.)) / k2)
      t0   = abs(0.5*(tan(theta0   / 2.) + \
          (1./3.)*pow(tan(theta0   / 2.),3.)) / k2)


    r0vect = array([\
        -r0 * cos( pi/2. - thetaimp + theta0 ),\
        r0 * sin( pi/2. - thetaimp + theta0 ),\
        0.])

    v0vect = array([\
        -v0 * cos( theta0 - thetaimp - beta0 ),\
        v0 * sin( theta0 - thetaimp - beta0 ),\
        0.])

    r0tar = -( mimp / mtot )*r0vect
    r0imp =  ( mtar / mtot )*r0vect

    v0tar = -( mimp / mtot )*v0vect
    v0imp =  ( mtar / mtot )*v0vect

    logstr = "#  parameters:\n" +\
        "# vinf   = %12.6e  cm/s  " % vinf + "\n" +\
        "# vesc   = %12.6e  cm/s  " % vesc + "\n" +\
        "# L0     = %12.6e gcm^2/s" % L0   + "\n" +\
        "# Lgraz  = %12.6e gcm^2/s" % Lgraz + "\n" +\
        "# gamma  = %12.6f        " % (mimp / mtar) + "\n" +\
        "# bscal  = %12.6f        " % bscal + "\n" +\
        "# e      = %12.6f        " % e     + "\n" +\
        "#" + "\n" +\
        "#  perihelon:" + "\n" +\
        "# r      = %12.6e cm   " % rperih + "\n" +\
        "# v      = %12.6e cm/s " % vperih + "\n" +\
        "# theta  = %12.6f deg (true anomaly) " % 0. + "\n" +\
        "# t      = %12.6f s    " % ( timp ) + "\n" +\
        "#" + "\n" +\
        "#  impact:" + "\n" +\
        "# r      = %12.6e cm   " % rimp   + "\n" +\
        "# v      = %12.6e cm/s " % vimp   + "\n" +\
        "# theta  = %12.6f deg (true anomaly) " % (thetaimp*rad2deg) + "\n" +\
        "# beta   = %12.6f deg (impact angle) " % (betaimp*rad2deg) + "\n" +\
        "# t      = %12.6f s    " % 0. + "\n" +\
        "#" + "\n" +\
        "#  initial setup:" + "\n" +\
        "# r      = %12.6f rimp " % relsep + "\n" +\
        "# r      = %12.6e cm   " % r0 + "\n" +\
        "# v      = %12.6e cm/s " % v0 + "\n" +\
        "# theta  = %12.6f deg (true anomaly) " % (theta0*rad2deg) + "\n" +\
        "# beta   = %12.6f deg (impact angle) " % (beta0*rad2deg) + "\n" +\
        "# t      = %12.6f s    " % ( timp - t0 ) + "\n" +\
        "#" + "\n" +\
        "# r0tar   = " + str(r0tar) + " cm" + "\n" +\
        "# r0imp   = " + str(r0imp) + " cm" + "\n" +\
        "# v0tar   = " + str(v0tar) + " cm/s" + "\n" +\
        "# v0imp   = " + str(v0imp) + " cm/s" + "\n" +\
        "#" + "\n"

    return (r0tar, r0imp, v0tar, v0imp, timp - t0, logstr)

  def getOrbit(self):
    pass


