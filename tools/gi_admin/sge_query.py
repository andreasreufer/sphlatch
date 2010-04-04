import time
import commands

class SGEquery(object):
  def __init__(self):
    self.jobslist = {}
    self.timestamp = 0

  def getJobs(self,notolderthan=0):
    if self.timestamp < ( time.time() - notolderthan):
      self._refreshData()
    return self.jobslist

  def _refreshData(self):
    jobslist = {}
    (stat, runwaitraw) = commands.getstatusoutput("qstat -g d")
    if stat != 0:
      return
    for line in runwaitraw.splitlines()[2:]:
      lsplt = line.split()
      if len(lsplt) == 9:
        (id, prio, name, user, statestr, date, qtime, queue, slots) = line.split()
        jobslist[int(id)] = (name, user, "run")
      if len(lsplt) == 8:
        (id, prio, name, user, statestr, date, qtime, slots) = line.split()
        state = "unknown"
        if statestr == "qw":
          state = "queued"
        jobslist[int(id)] = (name, user, state)
  
    (stat, finishedraw) = commands.getstatusoutput("qstat -g d -s z")
    if stat != 0:
      return unknown
    for line in finishedraw.splitlines()[2:]:
      (id, prio, name, user, statestr, date, qtime, slots) = line.split()
      jobslist[int(id)] = (name, user, "finished")
  
    self.timestamp = time.time()
    self.jobslist = jobslist


class SGEdummy(object):
  def __init__(self):
    jobslist[18131] = ('pf-1', 'wbenz', 'run')
    jobslist[18132] = ('pf-1', 'wbenz', 'finished')

  def getJobs(self,notolderthan=0):
    return self.jobslist

  

