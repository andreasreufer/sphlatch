#!/usr/bin/python

from simulation import Logger, SimSet
from sge_query import SGEquery

log = Logger("gi_admin.log")
simset = SimSet(log)
sim = simset.sims['mtar000.100_mimp000.100_impa75.0_vimp03.0']


pairs = sim._findBodies()

