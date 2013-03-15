#!/usr/bin/env ipython

import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt

from sim_admin import SimAdmin, SimParam, SimSetConfig
from simulation import resolvePath
from gi_plot import *
from gi_viz  import GIviz, GIvizConfig
from clumps import ClumpsPlotConfig
import os
import stat
import time
import numpy as np
import commands
import socket
import datetime

nan = float('nan')
eVinK = 11604.505

sscfg = SimSetConfig()
sscfg.name = "f3"

simadm = SimAdmin(sscfg)
sims = simadm._sims


