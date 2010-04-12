#!/usr/bin/python

import getpass

from simulation import SimSet
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
  
  def __init__(self, cfgfile):
    self._simset = SimSet("t0_config.sh")
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

  def showAll(self):
    for key in self._sims:
      sim = self._sims[key]
      col = bcolors.black

      if sim.state == "error":
        col = bcolors.lred
      elif sim.state == "run":
        col = bcolors.lgren
      elif sim.state == "unprepared":
        col = bcolors.cyan
      elif sim.state == "prepared":
        col = bcolors.purple
      elif sim.state == "submitted":
        col = bcolors.brown
      elif sim.state == "queued":
        col = bcolors.blue
      elif sim.state == "finished":
        col = bcolors.lgrey

      print key + "    " + col + sim.state + bcolors.black

  def _updateSGE(self):
    ssname = self._simset.simsetcfg["SIMSETNAME"]
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



