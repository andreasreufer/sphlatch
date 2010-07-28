#!/usr/bin/python

import shelve
import re
import os
import sys
import tables as pt
import numpy as np

from body import BodyFile

def getMetaData(dump):
  res = np.median(dump.h)
  m   = np.sum(dump.m)
  com = np.sum(dump.pos[:,:] * dump.m[:], axis=0) / m
  ri  = np.sqrt( np.sum( ( dump.pos[:,:] - com ) * ( dump.pos[:,:] - com ), \
      axis=1) )
  nop = dump.m.shape[0]
  r    = np.max(ri)
  attrs = {}
  for key in dump._v_attrs._f_list():
    attrs[key] = float( dump._v_attrs[key] )
  return (m, r, res, nop, attrs)

sfile = sys.argv[1]

wdir = "./"
if len(sys.argv) > 2:
  wdir = sys.argv[2]


bodydb = shelve.open(sfile)

for file in os.listdir(wdir):
  if re.search('\.h5part', file):
    ffile = os.path.abspath(wdir + "/" + file)
    print "match! "+ffile
    dumph = pt.openFile(ffile, "r")
    dump = dumph.root.current
    
    (m, r, h, nop, attrs) = getMetaData(dump)
    T = float(file[16:20])

    bodydb[ffile] = BodyFile(ffile, m, r, T, h, nop, attrs)
    bodydb.sync()
    dumph.close()

bodydb.sync()
bodydb.close()

