from numpy import *


def rotmatr(theta, u):
  ux = u[0]
  uy = u[1]
  uz = u[2]

  c = cos(theta)
  s = sin(theta)

  R = array([[ux*ux + (1.-ux*ux)*c, ux*uy*(1.-c) - uz*s, ux*uz*(1.-c) + uy*s],\
             [ux*uy*(1.-c) + uz*s, uy*uy + (1.-uy*uy)*c, uy*uz*(1.-c) - ux*s],\
             [ux*uz*(1.-c) - uy*s, uy*uz*(1.-c) + ux*s, uz*uz + (1.-uz*uz)*c]])

  return R

class Satellite(object):

  m = 0.
  mu = 0.
  K = 0.
  r0 = zeros(3)
  v0 = zeros(3)
  
  rM0 = zeros(3)
  rm0 = zeros(3)

  beta0 = 0.
  
  theta0 = 0.
  e = 0.
  i = 0.
  argp = 0.
  longa = 0.

  rp = 0.
  rc = 0.

  def __init__(self, rm, vm, m, rM, vM, M, rc, G):
    self.m  = m
    self.mu = ( m + M ) # reduced mass of central body
    self.K = self.mu*G
    self.r0 = rm - rM
    self.v0 = vm - vM
    
    self.r0M = rM
    self.r0m = rm

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
    self.rp = r0abs * (1.+self.e*cos(self.theta0))/(1.+self.e)

    self.i     = arccos( h[2] / habs )
    self.argp  = arccos( dot( ev, nv ) / (self.e*nvabs) )
    if ( h[0] > 0. ):
      self.longa = arccos( -h[1] / nvabs )
    else:
      self.longa = 2.*pi - arccos( -h[1] / nvabs )
    
    if ( not nvabs > 0. ):
      self.argp  = 0.
      self.longa = 0.
      print "nvabs is zero!!! ",nvabs

    self.rc = rc
    print "r0 ",self.r0, r0abs
    print "theta0 ",self.theta0 / pi
    print "h ",h, habs
    print "ev ",ev
    print "nv ",nv, nvabs
    print "e, rp ",self.e, self.rp
    print "i, argp, longa ", self.i/pi, self.argp/pi, self.longa/pi

  def plotorbit(self, ax):

    xax = array( [1., 0., 0.] )
    zax = array( [0., 0., 1.] )

    R = dot( rotmatr(self.longa, zax), \
            dot( rotmatr(self.i, xax), rotmatr(self.argp, zax) ) )

    # get minimal & maximal theta
    theta = arange( -pi, pi, (pi/100.) )
    #theta = array( [-pi, 0., pi] )
    #theta = arange( -self.theta0, pi, (pi/100.) )

    r = zeros( (theta.size, 3), dtype=double)

    r[:,0] = self.rp*(1. + self.e)/(1.+self.e*cos(theta))*cos(theta)
    r[:,1] = self.rp*(1. + self.e)/(1.+self.e*cos(theta))*sin(theta)

    rcart = dot( r, R )
    
    rcart += self.r0M
    #rcart[:,0] += self.r0M[0]
    #rcart[:,1] += self.r0M[1]

    ax.plot( rcart[:,0], rcart[:,1], lw=0.3, color='red')
    #ax.plot( r[:,0], r[:,1], lw=0.3, color='red')
    return ( rcart, r, R )


class SatSystem(object):
  sats = []

  r0 = zeros(3)
  v0 = zeros(3)
  rc = 0.
  m  = 0.

  G    = 0.
  mmin = 0.


  def __init__(self, parts, G, mmin):
    cidx = where( parts.m[:] >= max( parts.m ) )[0][0]
    sidx = where( ( parts.m[:] >= mmin ) & ( parts.m < max( parts.m ) ) )[0]

    for i in sidx:
      self.sats.append( Satellite(parts.pos[i,:], parts.vel[i,:], \
          float(parts.m[i]), parts.pos[cidx,:], parts.vel[cidx,:], \
          float(parts.m[cidx]), float(parts.rc[i]), G) )
      

    self.r0 = parts.pos[cidx,:]
    self.v0 = parts.vel[cidx,:]
    self.m  = float( parts.m[cidx] )
    self.rc = float( parts.rc[cidx] )

    self.G = G
    self.mmin = mmin
