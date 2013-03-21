import time
import commands

class PBSquery(object):
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
    (stat, out) = commands.getstatusoutput("qstat -f " + str(jobid))
    if stat == 0:
      for line in out.splitlines():
        if line.count("Job_Name"):
          return line.split()[2]
    return None

  def delJobID(self, jobid):
    (stat, out) = commands.getstatusoutput("qdel " + str(jobid))
    return stat

  def delJobName(self, jobname):
    sgejobs = getJobs() 
    for sgejob in sgejobs:
      (name, user, state) = sgejob
      if name == jobname:
        (stat, out) = commands.getstatusoutput("qdel " + str(jobid))
        return

  def _refreshData(self):
    jobslist = {}
    # TODO: most probably the wrong parameter
    (stat, runwaitraw) = commands.getstatusoutput("qstat -f")
    if stat != 0:
      return
    for line in runwaitraw.splitlines()[2:]:
      (idstr, abbrname, user, timeused, statestr, queue) = line.split()

      state = "unknown"
      if statestr == "Q":
        state = "queued"
      if statestr == "R":
        state = "run"
      if statestr == "E":
        state = "finished"

      jobslist[ idstr.split('.')[0] ] = (abbrname, user, state)
  
    self.timestamp = time.time()
    self.jobslist = jobslist




  

