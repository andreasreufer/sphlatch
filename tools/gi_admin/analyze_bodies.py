#!/usr/bin/python

import shelve
import re
import os
import sys
import tables

import body

sfile = sys.argv[1]

wdir = "./"
if len(sys.argv) > 2:
  wdir = sys.argv[2]



for file in os.listdir(wdir):
  if re.search('\.h5part', file):
    print "match! "+file


#bodies = [];
#bodies.append( (1.0, 2.e39) )

#dumpfile = open('bodies.pickle', 'w')
#P = pickle.Pickler(dumpfile)
#P.dump(bodies)
#dumpfile.close()

