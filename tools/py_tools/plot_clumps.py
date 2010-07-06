#!/usr/bin/env ipython

import numpy as np
import sys
from h5part import H5PartDump

import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt

MEg   = 5.9736e27

#fname = "/home0/areufer/pmc_colls/r0/sim_mtar001.000_mimp001.000_impa00.0_vimp01.0/dump0005778_T-1.0085e+04.h5part"
fname = "/Users/areufer/Documents/UniTAPS/pmc_colls/r0/sim_mtar001.000_mimp001.000_impa15.0_vimp01.5/clumps.h5part"

iname = sys.argv[1]
oname = sys.argv[2]
key   = sys.argv[3]

dump = H5PartDump(iname)


def getStepNum(sname):
  return int(sname.replace("/Step#",""))

snames = dump.getStepNames()
snames.sort(key=getStepNum)

nos = len(snames)
maxnoc = 0

for sname in snames:
  step = dump.getStep(sname)
  maxnoc = max(maxnoc,step.m.shape[0])

clmpmasses = np.zeros([nos,maxnoc])

cs = 0
for sname in snames:
  step = dump.getStep(sname)
  cnoc = step.m.shape[0]
  clmpmasses[cs,0:cnoc] = step.m[:,0]
  cs += 1

fig = plt.figure()
ax  = fig.add_subplot(111)


ax.semilogy(clmpmasses[:,0] / MEg, "r-")
for i in range(1,maxnoc):
  ax.semilogy(clmpmasses[:,i] / MEg, "b-")

plt.ylim(1.e-3,2.10)
plt.title(key)
plt.savefig(oname)



