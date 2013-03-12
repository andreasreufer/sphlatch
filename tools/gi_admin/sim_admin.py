#!/usr/bin/python

import getpass
import os
import tarfile
import os.path as path
import shelve

def resolvePath(_path):
    return path.abspath( path.expanduser(_path) )

from simulation import SimSet, SimParam, SimSetConfig, Simulation
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
    self.nosgejobs = 0
    
    self._updateSGE()
    
    for sim in self._sims.values():
      sim.getState()

    if path.exists(ssetcfg.resultsdb):
      resdb = shelve.open(resolvePath(self._simset.cfg.resultsdb))

      for key in resdb.keys():
        if self._sims.has_key(key) and resdb.has_key(key):
          self._sims[key].results = resdb[key]
      resdb.close()

  def __del__(self):
    pass

  def update(self):
    self._simset.discoverSims()
    self._updateSGE()
    for sim in self._sims.values():
      if sim.state == "run" or sim.state == "stuck":
        sim.getState()

  def storeResults(self):
    resdb = shelve.open(resolvePath(self._simset.cfg.resultsdb))
    for key in self._sims.keys():
      resdb[key] = self._sims[key].results
    resdb.close()

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

  def getSimsByFilter(self, filt=None):
    if filt==None:
      nan = float('nan')
      filt = SimParam(nan, nan, nan, nan)
    
    filtSims = []
    for sim in self._sims.values():
      if sim.params.isSimilar(filt, 0.01):
        filtSims.append(sim)
    return filtSims

  def getSimsByState(self, state):
    filtSims = []
    for sim in self._sims.values():
      if sim.state == state:
        filtSims.append(sim)
    return filtSims

  def lastLogStuck(self, ll):
    for sim in self.getSimsByState("stuck"):
      print sim.params.key
      sim._lastLog(ll)
      print "-----------------------------------------"

  def archiveSims(self, tfname, sims):
    archfiles = []
    tf = tarfile.open(tfname, mode="w:gz")

    for sim in sims:
      archfiles.extend( sim._getArchiveFiles() )

    cwd = os.getcwd() + "/"
    for f in archfiles:
       relf = f.replace(cwd, "")
       print relf
       tf.add(relf)
    #print archfiles
    tf.close()
    

  def newSim(self, params):
    if self._sims.has_key(params.key):
      print "sim: ",params.key," already exists!"
    else:
      csim = Simulation( params, self._simset.cfg )
      csim.next()
      self._sims[params.key] = csim
      return csim

  def newSimMatrix(self,mimpa,mtara,impaa,vimpa,cond="True"):
    for mimp in mimpa:
      for mtar in mtara:
        for impa in impaa:
          for vimp in vimpa:
            if ( cond ):
              self.newSim( SimParam(mimp, mtar, impa, vimp) )
  
  def newSimMatrixVimpImpaaComb(self,mimpa,mtara,vimpaimpaa,cond="True"):
    for mimp in mimpa:
      for mtar in mtara:
         for (impa, vimp) in vimpaimpaa:
            if ( cond ):
              self.newSim( SimParam(mimp, mtar, impa, vimp) )

  def forEachSim(self,func):
    for sim in self._sims.values():
      func(sim)

  def showState(self, filt=None):
    self.update()
    showSims = self.getSimsByFilter(filt)
    showSims.sort(key=lambda Simulation: Simulation.params.key)
    for sim in showSims:
      str = sim.params.key + "   " + self._colorState(sim.state) + "   " +\
          "(%3i/" % sim.nodumps + "%3i)" % sim.totdumps
      if sim.nodumps > 0:
         (lfile, ltime, lmtime) = sim.dumps[-1]
         str += "  %6.2f" % (ltime / sim.tcol)
      print str


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
    elif state == "stuck":
      col = bcolors.red
    elif state == "slow":
      col = bcolors.green
    return col + "{0:{1}}".format(state, width) + bcolors.black


  def _updateSGE(self):
    ssname = self._simset.cfg.name
    sgejobs = self._sgeqry.getJobs()
    self.nosgejobs = 0

    for id in sgejobs:
      (sname, user, state) = sgejobs[id]
      if (user == self.user):
        fname = self._sgeqry.getFullname(id)
	if len(fname) < 1:
	  continue
	if fname.count("_") < 1:
	  continue
        prefix = fname.split("_")[0]
        if prefix == ssname:
          simname = fname.replace(prefix + "_", "")
          if self._simset.sims.has_key(simname):
            self._simset.sims[simname].setSGEjobid(id)
            self._simset.sims[simname].setSGEstat(state)
	    if state == "queued" or state == "run":
	       self.nosgejobs += 1





