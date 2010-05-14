import time
import commands

class SGEquery(object):
  def __init__(self, notolderthan=10):
    self.jobslist = {}
    self.timestamp = 0
    self.notolderthan = notolderthan

  def getJobs(self,notolderthan=None):
    if notolderthan == None:
      notolderthan = self.notolderthan
    if self.timestamp < ( time.time() - notolderthan):
      self._refreshData()
    return self.jobslist
  
  def getFullname(self, jobid):
    (stat, out) = commands.getstatusoutput("qstat -j " + str(jobid))
    if stat == 0:
      for line in out.splitlines():
        if line.count("job_name:"):
          return line.split()[1]
    else:
      (stat, out) = commands.getstatusoutput("qacct -j " + str(jobid))
      if stat == 0:
        for line in out.splitlines():
          if line.count("jobname"):
            return line.split()[1]
    
    return None



  def delJobID(self, jobid):
    (stat, out) = commands.getstatusoutput("qdel -j " + str(jobid))
    return stat


  def delJobName(self, jobname):
    sgejobs = getJobs() 
    for sgejob in sgejobs:
      (name, user, state) = sgejob
      if name == jobname:
        (stat, out) = commands.getstatusoutput("qdel -j " + str(jobid))
        return


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
  
    #(stat, finishedraw) = commands.getstatusoutput("qstat -g d -s z")
    #if stat != 0:
    #  return unknown
    #for line in finishedraw.splitlines()[2:]:
    #  (id, prio, name, user, statestr, date, qtime, slots) = line.split()
    #  jobslist[int(id)] = (name, user, "finished")
  
    self.timestamp = time.time()
    self.jobslist = jobslist



class SGEdummy(object):
  def __init__(self):
    jobslist[18131] = ('pf-1', 'wbenz', 'run')
    jobslist[18132] = ('pf-1', 'wbenz', 'finished')

  def getJobs(self,notolderthan=0):
    return self.jobslist

  

