#!/usr/bin/env ipython

import matplotlib as mp
mp.use('Agg')
#mp.use('macosx')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

from sim_admin import SimAdmin, SimParam, SimSetConfig
from simulation import resolvePath
from gi_plot import *
from gi_viz  import GIviz, GIvizConfig
from clumps import ClumpsPlotConfig
from plot_helpers import *
import os
import stat
import time
import numpy as np
import pylab as pl

from plot_func import *

#rcParams['text.usetex']=True
mp.rc('text', usetex=True)
mp.rc('text.latex', preamble = '\usepackage{amssymb}, \usepackage{wasysym}')

nan = float('nan')

ssname = "c1"
if len(sys.argv) > 1:
  ssname = sys.argv[1]

sscfg = SimSetConfig()
sscfg.name = ssname
sscfg.dir = "/Users/areufer/Documents/UniTAPS/gi_sims"

simadm = SimAdmin(sscfg)

fig = plt.figure( figsize=(8, 6) )

x0 = 0.12
y0 = 0.10

dx = 0.21
dy = 0.27

params = []
if ssname == "c1":
  params = [ \
      ( SimParam(0.200, 1.000, nan, nan), [2,0] ),
      ( SimParam(0.700, 1.000, nan, nan), [0,0] ),
      ( SimParam(0.010, 0.100, nan, nan), [3,1] ),
      ( SimParam(0.020, 0.100, nan, nan), [2,1] ),
      ( SimParam(0.035, 0.100, nan, nan), [1,1] ),
      ( SimParam(0.070, 0.100, nan, nan), [0,1] ),
      ( SimParam(0.002, 0.010, nan, nan), [2,2] ),
      ( SimParam(0.007, 0.010, nan, nan), [0,2] ) ]

if ssname == "i1":
  params = [ \
      ( SimParam(0.200, 1.000, nan, nan), [2,0] ),
      ( SimParam(0.020, 0.100, nan, nan), [2,1] ),
      ( SimParam(0.002, 0.010, nan, nan), [2,2] ) ]

if ssname == "r3":
  params = [ \
      ( SimParam(0.100, 0.200, nan, nan), [0,2] ),
      ( SimParam(0.100, 0.500, nan, nan), [1,1] )]#,\
     #( SimParam(0.100, 1.000, nan, nan), [2,0] ) ]


bgax = plt.axes( [0.0, 0.0, 1.0, 1.0], frameon=False)
axs = []

zset = {}
for ( sparm, axv ) in params:
  sims = simadm.getSimsByFilter( sparm )
  ax   = plt.axes( [x0 + axv[0]*dx, y0 + axv[1]*dy, dx, dy] )

  sct2arr = ScatterToArray()
  valsims = 0
  for sim in sims:
    x = 0
    z = 0
    if xvar == "vimp":
      x = sim.params.vimprel
      z = sim.params.impa
    else:
      x = sim.params.impa
      z = sim.params.vimprel

    if filterFunc(sim):
      sct2arr.newPoint(z, x, plotFunc(sim) )
      zset[z] = []
      
      valsims += 1

  scta = sct2arr

  print len(sims), valsims, axv[0], axv[1]
  (zl, xll, resll) = sct2arr.getArray()

  for i in range(0, len(zl) ):
    (col, ls) = ("black", "+-") 
    if xvar == "vimp":
      (col, ls) = getAngleColor(zl[i])
    else:
      (col, ls) = getVimpColor(zl[i])

    if ylog:
      ax.semilogy( xll[i], resll[i], ls, color=col, markersize=6)
    else:
      ax.plot( xll[i], resll[i], ls, color=col, markersize=6)

  #Vhit = pl.loadtxt( "gamma_%3.2f0.txt" % ( sim.impb.m / sim.tarb.m ) )
  #ax.plot( Vhit[:,1], Vhit[:,2], 'k--', label="$V_{hit}$")
  
  if xvar == "vimp":
    ax.axis( [ 1.0, 4.0, yaxis[0], yaxis[1] ])
    ax.xaxis.set_ticks( (1.0, 2.0, 3.0, 4.0) )
    ax.xaxis.set_ticklabels( ("" , "" , "" , "") )
  else:
    ax.axis( [ 0., 90., yaxis[0], yaxis[1] ])
    ax.xaxis.set_ticks( (0.,15.,30.,45.,60.,75., 90.) )
    ax.xaxis.set_ticklabels( ("" , "" , "" , "" , "", "", "" ) )
  ax.grid(True)
  
  ax.yaxis.set_ticks( ytik )
  ylbll = []
  for i in range(0, len(ytik)):
    ylbll.append("")

  ax.yaxis.set_ticklabels( tuple( ylbll ) )

  axs.append(ax)

axselect = ()
if ssname == "c1" or ssname == "r3":
  axselect = (0,1)
if ssname == "i1":
  axselect = (0,)
for i in axselect:
  if xvar == "vimp":
    axs[i].set_xlabel(r"$v_{imp} / v_{esc}$")
    axs[i].xaxis.set_major_formatter(math_formatter)
  else:
    axs[i].set_xlabel(r"$\theta_{imp}  [^\circ]$")
    axs[i].xaxis.set_major_formatter(int_formatter)

axselect = ()
if ssname == "r3":
  axselect = (0,1)
if ssname == "c1":
  axselect = (1,5,7)
if ssname == "i1":
  axselect = (0,1,2)
for i in axselect:
  axs[i].yaxis.set_major_formatter(yfmt)
  axs[i].set_ylabel(ylbl)


if ssname == "c1":
  bgax.text( x0 + 3*dx + 0.07, 0.93, r"$\gamma = 0.10 $" )
  bgax.text( x0 + 2*dx + 0.07, 0.93, r"$\gamma = 0.20 $" )
  bgax.text( x0 + 1*dx + 0.07, 0.93, r"$\gamma = 0.35 $" )
  bgax.text( x0 + 0*dx + 0.07, 0.93, r"$\gamma = 0.70 $" )
  bgax.text( 0.03, y0 + (0+0.7)*dy, r"$M_{tar} = 1.0 M_{\oplus}$", rotation='vertical')
  bgax.text( 0.03, y0 + (1+0.7)*dy, r"$M_{tar} = 0.1 M_{\oplus}$", rotation='vertical')
  bgax.text( 0.03, y0 + (2+0.7)*dy, r"$M_{tar} = 0.01 M_{\oplus}$", rotation='vertical')
  
if ssname == "i1":
  bgax.text( x0 + 2*dx + 0.07, 0.93, r"$\gamma = 0.20 $" )
  bgax.text( 0.03, y0 + (0+0.7)*dy, r"$M_{tar} = 1.0 M_{\oplus}$", rotation='vertical')
  bgax.text( 0.03, y0 + (1+0.7)*dy, r"$M_{tar} = 0.1 M_{\oplus}$", rotation='vertical')
  bgax.text( 0.03, y0 + (2+0.7)*dy, r"$M_{tar} = 0.01 M_{\oplus}$", rotation='vertical')
  
if ssname == "r3":
  bgax.text( x0 + 2*dx + 0.07, 0.93, r"$\gamma = 0.10 $" )
  bgax.text( x0 + 1*dx + 0.07, 0.93, r"$\gamma = 0.20 $" )
  bgax.text( x0 + 0*dx + 0.07, 0.93, r"$\gamma = 0.50 $" )
  bgax.text( 0.03, y0 + (0+0.7)*dy, r"$M_{tar} = 1.0 M_{\oplus}$", rotation='vertical')
  bgax.text( 0.03, y0 + (1+0.7)*dy, r"$M_{tar} = 0.5 M_{\oplus}$", rotation='vertical')
  bgax.text( 0.03, y0 + (2+0.7)*dy, r"$M_{tar} = 0.2  M_{\oplus}$", rotation='vertical')
  


bgax.text( 0.03, 0.93, r"$\mathrm{" + sscfg.name + r"}$", size=20)


zvals = zset.keys()
zvals.sort()
for z in zvals:
  if xvar == "vimp":
    (col, ls) = getAngleColor(z)
    anglstr = "$%2.0f ^\circ$" % z
    bgax.plot( 99., 99., ls, label=anglstr, color=col, markersize=6)

  else:
    (col, ls) = getVimpColor(z)
    vimpstr = "$%3.2f$" % z
    bgax.plot( 99., 99., ls, label=vimpstr, color=col, markersize=6)


bgax.axis( [0., 1., 0., 1.] )

legtit = ""
if xvar == "vimp":
  legtit = r"$\theta_{imp} [^\circ]$"
else:
  legtit = r"$v_{imp} / v_{esc}$"

bgax.legend(loc=(0.80, 0.06), frameon=False, numpoints=1, ncol=2, columnspacing=0.5, prop=FontProperties(size=8), title=legtit)

bgax.xaxis.set_ticks( () )
bgax.yaxis.set_ticks( () )

plt.savefig("out_" + ssname + ".pdf")
exit()

