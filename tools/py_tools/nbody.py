import numpy as np

class Clumps(object):
  pos = []
  vel = []

  Mm = 0.
  Mpos = []
  Mvel = []
  iM = []

  def __init__(self, arg):
    self.pos = arg.pos
    self.vel = arg.vel
    self.Mm = float( max( arg.m ) )

    iM = int( (np.nonzero( arg.m >= max( arg.m ) ))[0] )
    self.iM = iM

    self.Mpos = self.pos[iM, :]
    self.Mvel = self.vel[iM, :]

   

