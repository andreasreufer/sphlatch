from numpy import *

class Satellite(object):

  mu = 0.
  K = 0.
  r0 = zeros(3)
  v0 = zeros(3)

  beta0 = 0.
  
  theta0 = 0.
  e = 0.
  i = 0.
  argp = 0.
  longa = 0.

  rc = 0.

  def __init__(self, rm, vm, m, rM, vM, M, rc, G):
    self.mu = ( m + M ) # reduced mass of central body
    self.K = self.mu*G
    self.r0 = rm - rM
    self.v0 = vm - vM

    r0abs = sqrt(dot(self.r0,self.r0))
    v0abs = sqrt(dot(self.v0,self.v0))

    k1          = r0abs * ( v0abs * v0abs ) / self.K
    self.beta0  = arcsin( dot(self.r0, self.v0) / (r0abs * v0abs) )
    self.theta0 = arctan2(k1 * sin(self.beta0)*cos(self.beta0), \
        ( k1 * cos(self.beta0)*cos(self.beta0) - 1. ) )
    
    ez    = array( [0., 0., 1.] )
    h     = cross( self.r0, self.v0 )
    habs  = sqrt(dot(h,h))
    ev    = ( cross( self.v0, h ) / self.K ) - ( self.r0 / r0abs )

    nv    = ( cross( ez, h ) )
    nvabs = sqrt( dot( nv, nv ) )
    
    self.e = sqrt( dot( ev, ev ) )

    self.i     = arccos( h[2] / habs )
    self.argp  = arccos( dot( ev, nv ) / (self.e*nvabs) )
    if ( h[0] > 0. ):
      self.longa = arccos( -h[1] / nvabs )
    else:
      self.longa = 2.*pi - arccos( -h[1] / nvabs )

    self.rc = rc
    
class SatSystem(object):
  sats = []


  mu = 0.
  K = 0.
  r0 = zeros(3)
  v0 = zeros(3)

  beta0 = 0.
  
  theta0 = 0.
  e = 0.
  i = 0.
  argp = 0.
  longa = 0.

  rc = 0.


  def __init__(self, parts, G, mmin):
    cidx = where( parts.m[:] >= max( parts.m ) )[0][0]
    sidx = where( ( parts.m[:] >= mmin ) & ( parts.m < max( parts.m ) ) )[0]

    for i in sidx:
      self.sats.append( Satellite(parts.pos[i,:], parts.vel[i,:], parts.m[i],\
          parts.pos[cidx,:], parts.vel[cidx,:], parts.m[cidx], parts.rc[i], G) )

