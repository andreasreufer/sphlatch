import os
import os.path as path
import stat
import shelve
import commands
import sys
import shutil
import tempfile
import random

from gi_plot import GIplotConfig, GIplot

class GIvizConfig(object):
  def __init__(self):
    self.scdir = "./givizscr"
    self.tasksperjob = 3
    self.subcmd = "qsub -cwd -l qname=all.q -N $JOBNAME -b n $JOBCMD"

class GIvizTask(object):
  def __init__(self, pfile, cfile, ifile, plotcfg):
    self.pfile = pfile
    self.cfile = cfile
    self.ifile = ifile

    self.plotcfg = plotcfg

  def execute(self):
    gipl = GIplot(self.plotcfg)
    gipl.plotParticles(self.pfile, self.cfile, self.ifile)


class GIviz(object):
  def __init__(self, cfg):
    self.cfg     = cfg
    self.tasks = []
    
    self.scdir = os.path.abspath( self.cfg.scdir ) + "/"
    if not path.exists(self.scdir):
      os.mkdir(self.scdir)
    

  def gatherVizTasks(self, sim, plotcfgs):
    scdir = self.scdir
    skey  = sim.params.key
    sparm = sim.params

    for pcfg in plotcfgs:
      ibdir = pcfg.imgdir
      if not path.exists(ibdir):
        os.mkdir(ibdir)
        print "make ",ibdir

      isdir = ibdir + "/" + skey
      if not path.exists(isdir):
        os.mkdir(isdir)
        print "make ",isdir

      pcfg.mtar = sparm.mtar
      pcfg.mimp = sparm.mimp
      pcfg.impa = sparm.impa
      pcfg.vimp = sparm.vimprel
      pcfg.T    = sparm.temp

    ntid = len(self.tasks)
    for drec in sim.dumps:
      (dfile, dtime) = drec
      dprfx = dfile.replace('.h5part','')
      pfile = sim.dir + dfile
      cfile = sim.dir + "clumps.h5part"

      for pcfg in plotcfgs:
        idir  = pcfg.imgdir + "/" + skey + "/"
        ifile = os.path.abspath( idir ) + "/" + dprfx + pcfg.imgext

        if path.exists(pfile) and \
            path.exists(cfile) and \
            not path.exists(ifile):
              #open(ifile,'w').close()
              ctask = GIvizTask(pfile, cfile, ifile, pcfg)
              ctask.id = ntid
              self.tasks.append(ctask)
              ntid += 1


  def clearTasks(self):
    self.tasks = []

  def runTasksSGE(self,jprfx="vizsim_"):
    jobs = []
    cjob = []
    tasks = self.tasks
    tasksperjob = self.cfg.tasksperjob

    jprfx += ( '%06i' % random.randint(0,1000000) ) + "_"
    slvname = self.scdir + jprfx + "shelve"
    drvname = self.scdir + jprfx + "driver.py"

    print slvname, drvname

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

      for task in job:
        jobstr += ( drvname + " " + str(task.id) + "\n" )
        taskdb[str(task.id)] = task

      jobfile = open(jobname,"w")
      print >>jobfile, jobstr
      jobfile.close()
      os.chmod(jobname, stat.S_IRWXU)

      jid += 1
    
    taskdb.close()
    
    if len(jobs) > 0:
      drvstr = "#!/usr/bin/env ipython\n" +\
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
    

  def runTasksLocal(self):
    while len(tasks) > 0:
      ( tasks.pop() ).execute

  def animSim(self, sim):
    idir  = self.plotcfg.imgdir + sim.params.key + "/"
    iext  = self.plotcfg.imgext
    animext = ".mp4"
    adir = self.dir + "videos/"
    afile = adir + sim.params.key + animext
    ffmpegargs = '-b 20000k'

    def isimg(str):
      return os.path.splitext(str)[1] == iext

    imglist = filter( isimg, os.listdir( idir ) )
    imglist.sort()
      
    tmpdir = tempfile.mkdtemp( dir=idir )

    count = 0
    for img in imglist:
      os.symlink( idir + img, tmpdir + "/img" + '%06d' % count + iext )
      count += 1

    if not path.exists(adir):
      os.mkdir(adir)

    if path.exists(afile):
      os.remove(afile)

    cmd =  "ffmpeg -i " + tmpdir + "/img%06d" + iext + " " +\
        ffmpegargs + " " + afile
    print cmd
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
    os.chdir(oldwd)


  def clearScratch(self):
    filelist = os.listdir( self.scdir )
    filelist.sort()

    for file in filelist:
      os.remove(self.scdir + file)



