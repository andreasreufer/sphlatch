#!/usr/bin/env ipython

import shelve
import re
import os
import sys
import tables as pt
import numpy as np

from body import BodyFile
from h5part import H5PartDump

nan = float('nan')

def getMetaData(dump):
  filt = ( dump.mat[:,0] > 0 )
  h   = dump.h[:,0].compress(filt)
  m   = dump.m[:,0].compress(filt)
  pos = dump.pos[:,:].compress(filt, axis=0)
  mat = dump.mat[:,0].compress(filt)
  

  mtot = np.sum(m[:])
  com  = np.sum( ( dump.pos[:,:] * dump.m[:] ).compress(filt, axis=0) , \
      axis=0) / mtot
  nop  = m.shape[0]
  rhook = (dump.rho[:,0] > 0.5)
  posfilt = dump.pos[:,:].compress( rhook, axis=0)[:,:] - com
  res  = np.mean( dump.h[:,0].compress( rhook, axis=0) )
  rmax = np.sqrt( np.max( np.sum( posfilt * posfilt, axis = 1 ) ) )

  mmat = np.zeros([32])
  rmat = np.zeros([32])
  for i in range(0,32):
    filt = ( dump.mat[:,0]  == i )
    mmat[i] = np.sum( dump.m[:,0].compress(filt) )

    mmatfilt = ( rhook & ( mat == i ) )
    posfilt  = dump.pos[:,:].compress( mmatfilt, axis=0)[:,:] - com

    if len(posfilt) > 0:
      rmat[i]  = np.sqrt( np.max( np.sum( posfilt * posfilt, axis = 1 ) ) )
    else:
      rmat[i]  = nan
  
  attrs = {}
  for key in dump._v_attrs._f_list():
    attrs[key] = float( dump._v_attrs[key] )

  return (mtot, mmat, rmax, rmat, res, nop, attrs)


sfile = sys.argv[1]

wdir = "./"
if len(sys.argv) > 2:
  wdir = sys.argv[2]


bodydb = shelve.open(sfile)

for file in os.listdir(wdir):
  if re.search('\.h5part', file):
    ffile = os.path.abspath(wdir + "/" + file)
    print "match! "+ffile
    dumpf = H5PartDump(ffile)
    dump  = dumpf.getStepFirst()

    (m, mmat, r, rmat, h, nop, attrs) = getMetaData(dump)

    rnewstr = raw_input("proposed radius:  %e " % r)
    try:
      r = float(rnewstr)
      print (" using new radius of %e " % r)
    except:
      pass

    #T = float(file[16:20])
    #T = float(file[17:24])
    #T = float(file[19:26])
    T = float('nan')

    bodydb[ffile] = BodyFile(ffile, m, mmat, r, T, h, nop, attrs, rmat=rmat)
    print ffile, " r = ", r, " nop = ",nop 
    bodydb.sync()

bodydb.sync()
bodydb.close()

