import numpy as np

class FewBodies(object):
  def __init__(self, parts, G, mmin):
    
    self.G = G
    csel = ( parts.m[:,0] > mmin )

    self.pos = (parts.pos[:,:]).compress(csel, axis=0)
    self.vel = (parts.vel[:,:]).compress(csel, axis=0)
    self.acc = np.zeros( self.pos.shape )
    self.m   = (parts.m[:,0]).compress(csel)
    self.rc  = (parts.rc[:,0]).compress(csel)
    
    self.s    = np.zeros( self.m.shape )
    self.__oacc = self.acc

    self.Mm = max( self.m )

    iM = (np.nonzero( self.m >= self.Mm ))[0][0]
    self.iM = iM
    self.noc = self.pos.shape[0]

    self.Mpos = self.pos[iM, :]
    self.Mvel = self.vel[iM, :]
  
    self.traj = []
    for i in range(self.noc):
      self.traj.append( np.array( [ self.pos[i,:] ] ) )
    
  def forces(self):
    self.acc = np.zeros( self.pos.shape )
    for i in range(self.noc):
      for j in range(self.noc):
        if ( i != j ):
          rij = self.pos[i,:] - self.pos[j,:]
          rree = np.dot( rij, rij ) + 1.e14
          re = np.sqrt( rree )
          #r = np.sqrt( np.dot( rij, rij ) )
          #rrr = r*r*r
          #self.acc[i,:] -= self.G * rij * self.m[j] / rrr
          self.acc[i,:] -= self.G * rij * self.m[j] / (rree*re)

  def advance(self, dt):
    dpos = self.vel*dt + self.acc*dt*dt
    self.pos += dpos
    for i in range(self.noc):
      self.s[i] += np.sqrt( np.dot( dpos[i,:], dpos[i,:] ) )

    self.__oacc = self.acc
    self.forces()
    self.vel += 0.5*(self.__oacc + self.acc)*dt
    self.t += dt

  def timestep(self, ds):
    dts = np.zeros( self.m.shape )
    dta = np.zeros( self.m.shape )
    
    # get timesteps until next dump
    for i in range(self.noc):
      dsi = ds - ( self.s[i] % ds )
      veli = np.sqrt( np.dot( self.vel[i,:], self.vel[i,:] ) )
      acci = np.sqrt( np.dot( self.acc[i,:], self.acc[i,:] ) )
      if ( acci > 0. ):
        dts[i] = ( -veli + np.sqrt( veli*veli + 4*acci*dsi ) ) / (2.*acci )
        dta[i] = 0.01*(veli / acci)
      else:
        print 'acci zero!'
        dts[i] = dsi / veli
        dta[i] = 1.e600
    #print 't = ',self.t,' dts = ', min(dts),' dta = ', min(dta)
    return min( min(dts), min(dta))

  def storepos(self, ds):
    for i in range(self.noc):
      dsi = min( ds - ( self.s[i] % ds ), self.s[i] % ds )
      if ( ( dsi / ds ) < 1.e-3 ):
        self.s[i] += dsi
        self.traj[i] = np.append(self.traj[i], [self.pos[i,:]], axis=0 )

  def integrate(self, t, ds):
    while (self.t < t):
      self.forces()
      self.advance(self.timestep(ds))
      self.storepos(ds)

