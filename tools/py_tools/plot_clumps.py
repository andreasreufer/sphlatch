#!/usr/bin/env ipython

import numpy as np
import sys

from clumps import SimClumps

import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt

MEg   = 5.9736e27

iname = sys.argv[1]
oname = sys.argv[2]
key   = sys.argv[3]

clumps = SimClumps(iname)
clmpmasses = clumps.m
maxnoc     = clumps.maxnoc

fig = plt.figure()
ax  = fig.add_subplot(111)

ax.semilogy(clmpmasses[:,0] / MEg, "r-")
for i in range(1,maxnoc):
  ax.semilogy(clmpmasses[:,i] / MEg, "b-")

plt.ylim(1.e-3,2.10)
plt.title(key)
plt.savefig(oname)



