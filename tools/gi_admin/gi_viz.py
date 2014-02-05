import os
import os.path as path
import stat
import shelve
import commands
import sys
import shutil
import tempfile
import random
import sim_machine

from gi_plot import GIplotConfig, GIplot


class GIvizConfig(object):
  def __init__(self):
    self.scdir = "./givizscr"
    self.tasksperjob = 50
    self.subcmd = sim_machine.qsub_sgl
    if hasattr(sim_machine, "qsub_sglscr"):
      self.python = sim_machine.python
    self.axrel = [-32., 32., -18., 18.]

class GIvizTask(object):
  def __init__(self, pfile, cfile, ifile, plotcfg):
    self.pfile = pfile
    self.cfile = cfile
    self.ifile = ifile

    self.plotcfg = plotcfg

  def execute(self):
    gipl = GIplot(self.plotcfg)
    gipl.doPlot(self.pfile, self.cfile, self.ifile)
    
    if not path.exists(self.ifile):
      print self.ifile," does not exists!"
      print os.uname()[1]
    else:
      if path.getsize(self.ifile) == 0:
      	print self.ifile," is empty!"
      	print os.uname()[1]

class GIviz(object):
  def __init__(self, cfg):
    self.cfg     = cfg
    self.tasks = []
    
    self.scdir = os.path.abspath( self.cfg.scdir ) + "/"
    if not path.exists(self.scdir):
      os.mkdir(self.scdir)
    
  def animSim(self, sim, pcfgs):
    skey  = sim.params.key
    
    for pcfg in pcfgs:
      ibdir  = os.path.abspath( pcfg.imgdir ) + "/"
      isdir  = ibdir + skey + "/"
      vdir   = ibdir + "videos" + "/"
      iext   = pcfg.imgext

      if path.exists(isdir):
        if self.vizDone(sim, pcfg):
          if not path.exists(vdir):
            os.mkdir(vdir)
          vfile = vdir + skey + ".mp4"
          self.animDir(isdir, iext, vfile)

  def clearSim(self, sim, pcfgs):
    skey  = sim.params.key
    
    for pcfg in pcfgs:
      ibdir  = os.path.abspath( pcfg.imgdir ) + "/"
      isdir  = ibdir + skey + "/"
      iext   = pcfg.imgext
    
      def isimg(str):
        return os.path.splitext(str)[1] == iext

      if not path.exists(isdir):
        continue
      
      imglist = filter( isimg, os.listdir( isdir ) )
    
      for img in imglist:
        ifile = isdir + img
        if ( os.path.getsize(ifile) == 0 ):
	  (exstat, out) = commands.getstatusoutput("rm " + ifile)

    state = sim.state
    sim.state = "run"
    self.animSim(sim, pcfgs)
    
    for pcfg in pcfgs:
      ibdir  = os.path.abspath( pcfg.imgdir ) + "/"
      isdir  = ibdir + skey + "/"
      (exstat, out) = commands.getstatusoutput("rm -rf " + isdir)

    sim.state = state



  def vizDone(self, sim, pcfg):
    skey  = sim.params.key
      
    ibdir  = os.path.abspath( pcfg.imgdir ) + "/"
    isdir  = ibdir + skey + "/"
    vdir   = ibdir + "videos" + "/"
    vfile  = vdir + skey + ".mp4" 
    iext   = pcfg.imgext

    if not path.exists(isdir):
      return False

    if path.exists(vdir):
      if path.exists(vfile) and not sim.state == "run":
        return False

    def isimg(str):
      return os.path.splitext(str)[1] == iext
    
    imglist = filter( isimg, os.listdir( isdir ) )
    imglist.sort()

    necount = 0
    totcount = len(imglist)
    for img in imglist:
      if ( os.path.getsize(isdir + img) > 0 ):
        necount += 1

    print skey,' ',pcfg.imgdir,':    ',necount,'/',totcount
    return (necount == totcount)

  def doSims(self, sims, pcfgs):
    for sim in sims:
      if sim.nodumps > 0:
        self.gatherVizTasks(sim, pcfgs)
        self.runTasksSGE()

  def animSims(self, sims, pcfgs):
    for sim in sims:
      if sim.nodumps > 0:
        self.animSim(sim, pcfgs)
  
  def gatherVizTasks(self, sim, pcfgs):
    skey  = sim.params.key

    nt = 0

    ntid = len(self.tasks)
    for pcfg in pcfgs:
      for drec in sim.dumps:
        nt += 1
        self.addTask(sim, drec, pcfg)
    if nt > 0:
      print ( skey + ": " + str(nt))


  def addTask(self, sim, drec, pcfg):
    (dfile, dtime, dmtime) = drec
    dprfx = dfile.replace('.h5part','')
    pfile = sim.dir + dfile
    cfile = sim.dir + "clumps.h5part"
      
    ibdir = pcfg.imgdir
    if not path.exists(ibdir):
      os.mkdir(ibdir)
    
    skey  = sim.params.key
    
    isdir = ibdir + "/" + skey
    if not path.exists(isdir):
      os.mkdir(isdir)
      
    idir  = pcfg.imgdir + "/" + skey + "/"
    ifile = os.path.abspath( idir ) + "/" + dprfx + pcfg.imgext
    
    sparm = sim.params

    pcfg.mtar = sparm.mtar
    pcfg.mimp = sparm.mimp
    pcfg.impa = sparm.impa
    pcfg.vimp = sparm.vimprel
    pcfg.T    = sparm.temp
    pcfg.r    = sim.tarb.r
    
    ntid = len(self.tasks)
    if path.exists(pfile) and \
        path.exists(cfile) and \
        not path.exists(ifile):
          open(ifile,'w').close()
          ctask = GIvizTask(pfile, cfile, ifile, pcfg)
          ctask.id = ntid
          self.tasks.append(ctask)

  def clearTasks(self):
    self.tasks = []
  

  def runTasksSGE(self,jprfx="vizsim_"):
    jobs = []
    cjob = []
    tasks = self.tasks
    tasksperjob = self.cfg.tasksperjob

    if len(tasks) == 0:
      return

    jprfx += ( '%06i' % random.randint(0,1000000) ) + "_"
    slvname = self.scdir + jprfx + "shelve"
    drvname = self.scdir + jprfx + "driver.py"

    while len(tasks) > 0:
      cjob.append( tasks.pop() )
      if len(cjob) == tasksperjob:
        jobs.append( cjob )
        cjob = []

    if len(cjob) > 0:
      jobs.append( cjob )
    
    taskdb = shelve.open(slvname)
    taskdb.clear()

    jid  = 0
    for job in jobs:
      jobname = self.scdir + jprfx + ( "%03i" % jid ) + '.sh'
	  
      jobstr = "#!/bin/bash\n"
      if hasattr(sim_machine, "qsub_sglscr"):
        jobstr = sim_machine.qsub_sglscr
        jobstr = jobstr.replace('$JOBCMD', '')
        jobstr = jobstr.replace('$JOBNAME', jobname)

      for task in job:
        jobstr += ( self.cfg.python + " " + drvname + " " + str(task.id) + "\n" )
        taskdb[str(task.id)] = task

      jobfile = open(jobname,"w")
      print >>jobfile, jobstr
      jobfile.close()
      os.chmod(jobname, stat.S_IRWXU)
      print jobname

      jid += 1
    
    taskdb.close()
    
    if len(jobs) > 0:
      drvstr = "#!/usr/bin/env python\n" +\
          "from gi_viz import GIvizTask\n" +\
          "from gi_plot import GIplot\n" +\
          "import shelve, sys\n" +\
          "taskkey = str(sys.argv[1])\n" +\
          "tasks = shelve.open(\"" + slvname + "\",flag='r')\n" +\
          "tasks[taskkey].execute()\n" +\
          "tasks.close()\n"

      drvfile = open(drvname,"w")
      print >>drvfile, drvstr
      drvfile.close()
      os.chmod(drvname, stat.S_IRWXU)
  
    # submit the job
    oldwd = os.getcwd()
    os.chdir(self.scdir)

    jid = 0
    for job in jobs:
      jobname = jprfx + ( "%03i" % jid )
      jobsubstr = self.cfg.subcmd

      jobsubstr = jobsubstr.replace("$JOBNAME", jobname )
      jobsubstr = jobsubstr.replace("$JOBCMD",  jobname + ".sh")

      (exstat, out) = commands.getstatusoutput(jobsubstr)
      print out
      jid += 1
    os.chdir(oldwd)
    

  def runTasksLocal(self):
    while len(self.tasks) > 0:
      ( self.tasks.pop() ).execute

  def animDir(self, idir, iext, outf, ffmpegargs = '-b 20000k'):
    def isimg(str):
      return os.path.splitext(str)[1] == iext

    imglist = filter( isimg, os.listdir( idir ) )
    imglist.sort()
      
    #tmpdir = tempfile.mkdtemp( dir=idir )
    tmpdir = tempfile.mkdtemp()

    count = 0
    for img in imglist:
      os.symlink( idir + "/" + img, tmpdir + "/img" + '%06d' % count + iext )
      count += 1

    if path.exists(outf):
      os.remove(outf)

    cmd =  "ffmpeg -i " + tmpdir + "/img%06d" + iext + " " +\
        ffmpegargs + " " + outf
    (exstat, out) = commands.getstatusoutput(cmd)

    count = 0
    for img in imglist:
      os.remove(tmpdir + "/img" + '%06d' % count + iext)
      count += 1
    os.removedirs(tmpdir)

  

  def submitJob(self, scriptstr, scriptname, jobname):
    script = open(scriptname,"w")
    print >>script, scriptstr
    script.close()
    os.chmod(scriptname, stat.S_IRWXU)

    jobsubstr = self.cfg.subcmd

    jobsubstr = jobsubstr.replace("$JOBNAME", "viz_" + jobname)
    jobsubstr = jobsubstr.replace("$JOBCMD",  scriptname)
    
    oldwd = os.getcwd()
    os.chdir(self.scdir)
    (exstat, out) = commands.getstatusoutput(jobsubstr)
    print out
    os.chdir(oldwd)


  def clearScratch(self):
    filelist = os.listdir( self.scdir )
    filelist.sort()

    for file in filelist:
      os.remove(self.scdir + file)

  def clearImgDir(self,sims,pcfgs):
    idirs = []
    for pcfg in pcfgs:
      idirs.append(pcfg.imgdir)

    for sim in sims:
      key = sim.params.key
      for idir in idirs:
        cdir = idir + "/" + key
        filelist = os.listdir( cdir )
        for cfile in filelist:
          os.remove(cdir+"/"+cfile)
        os.removedirs(cdir)      
    


