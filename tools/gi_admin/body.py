#!/usr/bin/env python

class BodyFile(object):
  def __init__(self, file, m, r, T, h, attrs):
    self.file = file
    
    self.m = m
    self.T = T
    self.h = h
    
    self.r = r
    self.attrs = attrs

  def isSimilar(self, bod, tol):
    # not logic, so that 'nan' acts like a wildcard
    fitMass = not abs( ( self.m - bod.m ) / self.m ) > tol
    fitRes  = not abs( ( self.h - bod.h ) / self.h ) > tol
    fitTemp = not abs( ( self.T - bod.T ) / self.T ) > tol
    return fitMass and fitRes and fitTemp


