#!/usr/bin/python

import getpass

from simulation import SimSet, SimParam, SimSetConfig
from sge_query import SGEquery

class bcolors:
      black  = '\033[0;30m'
      red    = '\033[0;31m'
      green  = '\033[0;32m'
      brown  = '\033[0;33m'
      blue   = '\033[0;34m'
      purple = '\033[0;35m'
      cyan   = '\033[0;36m'
      lgrey  = '\033[0;37m'
      
      dgrey  = '\033[1;30m'
      lred   = '\033[1;31m'
      lgren  = '\033[1;32m'
      yellow = '\033[1;33m'
      lblue  = '\033[1;34m'
      lpurpl = '\033[1;35m'
      lcyan  = '\033[1;36m'
      white  = '\033[1;37m'

class SimAdmin(object):
  def __init__(self, ssetcfg):
    self._simset = SimSet(ssetcfg)
    self._sgeqry = SGEquery()

    self._sims = self._simset.sims
    self.user = getpass.getuser()
    
    for sim in self._sims.values():
      sim.getState()
    self._updateSGE()

  def __del__(self):
    pass

  def gatherTasks(self):
    pass

  def doTasks(self):
    pass

  def getJobs(self, state):
    sims = []
    for key in self._sims:
      sim = self._sims[key]
      if sim.state == state:
        sims.append(sim)
    return sims

  def getSims(self, filt=None):
    if filt==None:
      nan = float('nan')
      filt = SimParam(nan, nan, nan, nan)
    
    filtSims = []
    for sim in self._sims.values():
      if sim.params.isSimilar(filt, 0.01):
        filtSims.append(sim)
    return filtSims

  def showState(self, filt=None):
    showSims = self.getSims(filt)
    showSims.sort(key=lambda Simulation: Simulation.params.key)
    for sim in showSims:
      print sim.params.key + "   " + self._colorState(sim.state) + "   " +\
          "(%3i/" % sim.nodumps + "%3i)" % sim.totdumps


  def _colorState(self, state, width=10):
    col = bcolors.black
    if state == "error":
      col = bcolors.lred
    elif state == "run":
      col = bcolors.lgren
    elif state == "unprepared":
      col = bcolors.cyan
    elif state == "prepared":
      col = bcolors.purple
    elif state == "submitted":
      col = bcolors.brown
    elif state == "queued":
      col = bcolors.blue
    elif state == "finished":
      col = bcolors.lgrey
    elif state == "failed":
      col = bcolors.red
    elif state == "partial":
      col = bcolors.lpurpl
    return col + "{0:{1}}".format(state, width) + bcolors.black


  def _updateSGE(self):
    ssname = self._simset.cfg.name
    sgejobs = self._sgeqry.getJobs()

    for id in sgejobs:
      (sname, user, state) = sgejobs[id]
      if (user == self.user):
        fname = self._sgeqry.getFullname(id)
        prefix = fname.split("_")[0]
        if prefix == ssname:
          simname = fname.replace(prefix + "_", "")
          if self._simset.sims.has_key(simname):
            self._simset.sims[simname].setSGEjobid(id)
            self._simset.sims[simname].setSGEstat(state)



