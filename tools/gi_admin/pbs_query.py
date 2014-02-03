import time
import commands
import xml.etree.ElementTree as ET

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
    if self.jobslist.has_key(jobid):
      (name, user, state) = self.jobslist[jobid]
      return name
    else:
      return ""
    
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
    (stat, runwaitraw) = commands.getstatusoutput("qstat -x")
    if stat != 0:
      return
    
    qroot = ET.fromstring(runwaitraw)

    for child in qroot:
      if child.tag == "Job":
        name = child.findtext('Job_Name')
        flid = child.findtext('Job_Id')
        sstr = child.findtext('job_state')
        user = child.findtext('Job_Owner').split('@')[0]
      
        state = "unknown"
        if sstr == "Q":
          state = "queued"
        if sstr == "R":
          state = "run"
        if sstr == "C":
          state = "partial"
        if sstr == "E":
          state = "partial"

        jobslist[ flid.split('.')[0] ] = (name, user, state)
  
    self.timestamp = time.time()
    self.jobslist = jobslist




  

