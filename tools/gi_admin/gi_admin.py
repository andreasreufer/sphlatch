#!/usr/bin/python

from simulation import SimSet
from sge_query import SGEquery

simset = SimSet("simset_config.sh")

sim = simset.sims['mtar000.100_mimp000.100_impa75.0_vimp03.0']
#pairs = sim._findBodies()
sim._prepare()

sim = simset.sims['mtar000.100_mimp000.100_impa75.0_vimp01.0']
sim._prepare()

sim = simset.sims['mtar001.000_mimp001.000_impa30.0_vimp01.0']
sim._prepare()

