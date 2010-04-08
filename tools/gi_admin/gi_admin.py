#!/usr/bin/python

from simulation import SimSet
from sge_query import SGEquery

simset = SimSet("simset_config.sh")

#sim = simset.sims['mtar000.100_mimp000.100_impa75.0_vimp03.0']
#sim = simset.sims['mtar001.000_mimp000.100_impa75.0_vimp01.0']

for sim in simset.sims:
  for i in range(2):
    print sim.next()


