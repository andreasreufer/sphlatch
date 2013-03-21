class DUMMYquery(object):
  def __init__(self, notolderthan=10):
    self.jobslist = {}
    self.timestamp = 0
    self.notolderthan = notolderthan

  def getJobs(self,notolderthan=None):
    return self.jobslist
  
  def getFullname(self, jobid):
    return None

  def delJobID(self, jobid):
    return -1


  def delJobName(self, jobname):
    return -1

  def _refreshData(self):
    jobslist = {}
    self.jobslist = jobslist

