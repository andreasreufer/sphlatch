#!/usr/bin/python

from simulation import SimSet
from sge_query import SGEquery

simset = SimSet("simset_config.sh")

#sim = simset.sims['mtar000.100_mimp000.100_impa75.0_vimp03.0']
sim = simset.sims['mtar001.000_mimp000.100_impa75.0_vimp01.0']

print sim.getState()
for i in range(3):
  print sim.next()

sim.setSGEjobid(42)
sim.setSGEstat("queued")
for i in range(3):
  print sim.next()

sim.setSGEstat("run")
for i in range(3):
  print sim.next()

sim.setSGEstat("finished")
for i in range(3):
  print sim.next()


