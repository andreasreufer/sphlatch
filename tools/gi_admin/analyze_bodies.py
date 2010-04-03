#!/usr/bin/python

import shelve
import re
import os
import sys
import tables as pt
import numpy as np

from body import BodyFile

sfile = sys.argv[1]

wdir = "./"
if len(sys.argv) > 2:
  wdir = sys.argv[2]



def getMetaData(dump):
  res = np.median(dump.h)
  m   = np.sum(dump.m)
  com = np.sum(dump.pos[:,:] * dump.m[:], axis=0) / m
  ri  = np.sqrt( np.sum( ( dump.pos[:,:] - com ) * ( dump.pos[:,:] - com ), \
      axis=1) )
  r    = np.max(ri)

  attrs = {}
  for key in dump._v_attrs._f_list():
    attrs[key] = float( dump._v_attrs[key] )
  
  return (m, r, res, attrs)


bodydb = shelve.open(sfile)

for file in os.listdir(wdir):
  if re.search('\.h5part', file):
  #if file == "body_A_0.100me_T1500K_1.0dun_018k.h5part":

    ffile = os.path.abspath(wdir + "/" + file)
    print "match! "+ffile
    dumph = pt.openFile(ffile, "r")
    dump = dumph.root.current
    
    (m, r, res, attrs) = getMetaData(dump)
    T = file[16:20]

    bodydb[ffile] = BodyFile(ffile, m, r, T, res, attrs)
    bodydb.sync()
    
    dumph.close()


bodydb.close()


#bodies = [];
#bodies.append( (1.0, 2.e39) )

#dumpfile = open('bodies.pickle', 'w')
#P = pickle.Pickler(dumpfile)
#P.dump(bodies)
#dumpfile.close()

