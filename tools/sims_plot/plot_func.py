#!/usr/bin/env ipython

from plot_helpers import *

nan = float('nan')

def plotFunc(sim):
  res = sim.results
  #ToW = abs(-res.Erotm[1] / res.Epotm[1])
  #ToW = abs(-res.Erotm[1] / sim.gi.Ered)
  #ToW = sim.results.Lm[1,2] / sim.gi.L
  #if ToW > 100.:
  #  ToW = nan
  #return (sim.results.mm[2:] / sim.impb.m).sum()
  #return sim.results.mm[0] / sim.impb.m
  return (res.mdiskmatm[1,:]).sum() / sim.mtot
  #return ToW

def filterFunc(sim):
  return sim.results.valid and ( sim.results.valtmax / sim.tcol > 5. )

xvar = "vimp"
#xvar = "angle"

yaxis = [ 1.e-3, 10.1 ]
ylog  = True
ylbl  = r"$E_{rot} / E_{pot}$"
yfmt  = tex_formatter
#ytik = (-2., -1., 0., 0.5, 1.0)
ytik = (1.e-2, 1.e-1, 1.0)


