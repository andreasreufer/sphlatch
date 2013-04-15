#!/usr/bin/env ipython

import numpy as np

from numpy import sqrt, sin, cos, arcsin, arccos, pi, isnan

def interpolateRhoMat(rcur, rprof, rhoprof, matprof):
  if rcur <= rprof[0]:
    return (rhoprof[0], matprof[0])
  if rcur >= rprof[-1]:
    return (rhoprof[-1], matprof[-1])

  ilower = (rprof < rcur).nonzero()[0].max()
  iupper = (rprof > rcur).nonzero()[0].min()

  rlower = rprof[ilower]
  rupper = rprof[iupper]
  
  k = (rcur - rlower)/(rupper - rlower)
  rho = ( (1.-k)*rhoprof[ilower] + k*rhoprof[iupper] )
  mat = matprof[ilower]
  if k > 0.5:
    mat = matprof[iupper]

  return (rho, mat)


def getRelHitMass(Rimp, Rtar, bscal, rprof, rhoprof, matprof, res=200, maxmat=15):
  Rtot = Rimp + Rtar
  Y = bscal*Rtot

  dx = ( Rimp / res )
  dr = ( Rimp / res )
  
  M = np.zeros(maxmat)
  V = np.zeros(maxmat)
  
  for x in np.arange( -Rimp, Rimp, dx ):
    
    # calculate rmin and rmax for given x
    rmin = max(0., Y-Rtar + 1.e-6*dr)
    rmax = min( sqrt( Rimp*Rimp - x*x ), Rtar )

    for r in np.arange( rmin, rmax, dr):
      phi = 0.

      if ( r < (Rtar-Y) ):
        phi = 2*pi
      elif ( Y > (r+Rtar) ):
        phi = 0.
      else:
        phi = 2.*arccos( ( Y*Y + r*r - Rtar*Rtar ) / ( 2.*r*Y) ) 

      rprime = sqrt( x*x + r*r )
      (rho, mat) = interpolateRhoMat( rprime, rprof, rhoprof, matprof)

      dV = dr*phi*r*dx

      V[mat] += dV
      M[mat] += dV*rho

  return (M, V)



