#!/usr/bin/env python
import numpy as np

class BodyFile(object):
  def __init__(self, file, m, mmat, r, T, h, nop, attrs, rmat=np.zeros([32])):
    self.file = file
    
    self.mmat = mmat
    self.T = T
    self.h = h
    self.m = m
    
    self.nop = nop
    self.r = r
    self.attrs = attrs
    self.rmat = rmat

  def isSimilar(self, bod, tol):
    # not logic, so that 'nan' acts like a wildcard
    fitMass = not abs( ( self.m - bod.m ) / self.m ) > tol
    fitRes  = not abs( ( self.h - bod.h ) / self.h ) > tol
    #fitTemp = not abs( ( self.T - bod.T ) / self.T ) > tol
    fitTemp = True
    return  fitMass and fitRes and fitTemp 


