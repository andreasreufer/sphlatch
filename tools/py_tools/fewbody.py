import numpy as np

class FewBodies(object):
  pos = []
  vel = []
  m = []
  acc = []
  __oacc = []
  d = []
  rc = []

  Mm = 0.
  Mpos = []
  Mvel = []

  iM = []
  noc = 0

  G = 0.

  def __init__(self, parts, G, mmin):

    self.G = G

    csel = ( parts.m[:,0] > mmin )

    self.pos = (parts.pos[:,:]).compress(csel, axis=0)
    self.vel = (parts.vel[:,:]).compress(csel, axis=0)
    self.acc = np.zeros( self.pos.shape )
    self.m   = (parts.m[:,0]).compress(csel)
    self.rc  = (parts.rc[:,0]).compress(csel)
    self.dout = np.zeros( self.m.shape )
    self.__oacc = self.acc

    self.Mm = max( self.m )

    iM = (np.nonzero( self.m >= self.Mm ))[0][0]
    self.iM = iM
    self.noc = self.pos.shape[0]

    self.Mpos = self.pos[iM, :]
    self.Mvel = self.vel[iM, :]

  def forces(self):
    self.acc = np.zeros( self.pos.shape )
    for i in range(self.noc):
      for j in range(self.noc):
        if ( i != j ):
          rij = self.pos[i,:] - self.pos[j,:]
          r = np.sqrt( np.dot( rij, rij ) )
          rrr = r*r*r
          self.acc[i,:] -= self.G * rij * self.m[j] / rrr

  def advance(self, dt):
    dpos = self.vel*dt + self.acc*dt*dt
    self.pos += dpos
    #for i in range(self.noc):
    #  self.dout -= dpos*dt


    self.__oacc = self.acc
    self.forces()
    self.vel += 0.5*(self.__oacc + self.acc)*dt




