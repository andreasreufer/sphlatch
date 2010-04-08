import os
import os.path as path
import time
import commands
import shelve
import numpy as np
import shutil

from numpy import sqrt, sin, cos, arccos, log, abs, tan, arctan2, pi, zeros
from body import BodyFile
from setup_gi import GiantImpact
import pdb



rad2deg = 360./(2.*pi)
deg2rad = 1./rad2deg

# constants to transform the initial parameters to physical values
G   = 6.67429e-8
Re  = 6.37814e8
Me  = 5.97360e27
lem = 3.40000e41 # Canup (2001) value
#lem = 2.88993e41 # actual (?) value

def resolvePath(_path):
  return path.abspath( path.expanduser(_path) )

class Logger(object):
  logfile = []

  def __init__(self, filename):
    self.logfile = open(filename, 'a')

  def __del__(self):
    self.logfile.close()

  def write(self, str):
    if len(str) > 0:
      self.logfile.write(time.ctime() + ":   " + str + "\n")
      self.logfile.flush()


class SimParams(object):
  def __init__(self, mimp, mtar, impa, vimprel, cfg):
    self.mimp = mimp
    self.mtar = mtar
    self.impa = impa
    self.vimprel = vimprel

    self.cfg = cfg

    self.temp = float(cfg["TEMP"])

    self.key = 'mtar%07.3f' % mtar + '_' \
      + 'mimp%07.3f' % mimp + '_' \
      + 'impa%04.1f' % impa + '_' \
      + 'vimp%04.1f' % vimprel

class Task(object):
  def __init__(self, callfunc):
    self.callfunc = callfunc

  def __run__(self):
    self.callfunc()


class Simulation(object):
  # noCPUs, 
  def __init__(self,params):
    self.params = params
    ssbdir = resolvePath(params.cfg["SIMSETBDIR"]) + "/"
    self.dir = resolvePath(ssbdir + "sim_" + params.key) + "/"
    
    self.jobid = 0

    self.state = "unprepared"
    self.next = self._prepare

  def setError(self, str=None):
    if not str == None and hasattr(self, "Log"):
      self.Log.write("error: " + str)
    self.state = "error"
    self.next = self._doNothing
    return self.state

  def getState(self):
    if self.state == "error":
      return self.state

    if not path.exists(self.dir):
      self.state = "unprepared"
      self.next = self._prepare
      return self.state
    else:
      neededfiles = self.params.cfg["AUXFILES"].split()
      neededfiles.append(self.params.cfg["BINARY"])
      neededfiles.append("initial.h5part")

      for file in neededfiles:
        if not os.path.exists(self.dir + os.path.split(file)[1]):
          self.state = "unprepared"
          self.next = self._prepare
          return self.state
      # so everything's ready
      if self.jobid == 0:
        self.state = "prepared" 
        self.next = self._submit
      return self.state

    self.state = "unknown"
    self.next = self._doNothing
    return self.state

  def _doNothing(self):
    return self.state

  def _prepare(self):
    # if it does not yet exist, make dir
    if not path.exists(self.dir):
      os.mkdir(self.dir)

    # copy files
    auxf = self.params.cfg["AUXFILES"].split()
    
    self.Log = Logger(self.dir + "logfile.txt")
    
    for file in auxf:
      shutil.copy2(resolvePath(file), self.dir)
    
    self.Log.write("auxiliary files copied")
    
    if path.exists(self.dir + "initial.h5part"):
      os.remove(self.dir + "initial.h5part")

    # prepare bodies:
    #  find fitting bodies
    (self.tarb, self.impb) = self._findBodies()
    
    # get orbit
    mtar = self.tarb.m
    mimp = self.impb.m
    Rtar = self.tarb.r
    Rimp = self.impb.r

    G = float(self.params.cfg["GRAVCONST"])
    relsep = float(self.params.cfg["RELSEP"])

    gi = GiantImpact(mtar, mimp, Rtar, Rimp, G)
    vimp = self.params.vimprel * gi.vesc
    impa = self.params.impa * deg2rad

    (r0tar, r0imp, v0tar, v0imp, t0, logstr) = \
        gi.getInitVimpAlpha(vimp, impa, relsep)
    gilog = open(self.dir + "gi_setup.log", "w")
    print >>gilog, logstr
    gilog.close()

    # displace bodies
    cmds = []
    cmds.append("cp -Rpv " + self.tarb.file + " " + self.dir + "tarb.h5part")
    cmds.append("h5part_displace -i " + self.dir + "tarb.h5part " +\
        "--pos [%e,%e,%e] " % (r0tar[0], r0tar[1], r0tar[2]) +\
        "--vel [%e,%e,%e] " % (v0tar[0], v0tar[1], v0tar[2]) )
    cmds.append("cp -Rpv " + self.impb.file + " " + self.dir + "impb.h5part")
    cmds.append("h5part_displace -i " + self.dir + "impb.h5part " +\
        "--pos [%e,%e,%e] " % (r0imp[0], r0imp[1], r0imp[2]) +\
        "--vel [%e,%e,%e] " % (v0imp[0], v0imp[1], v0imp[2]) +\
        "--id 2000000" )
    cmds.append("h5part_combine -a " + self.dir + "tarb.h5part -b " +\
        self.dir + "impb.h5part -o " + self.dir + "initial.h5part")
    cmds.append("rm " + self.dir + "tarb.h5part " + self.dir + "impb.h5part")
    
    # write attributes (gravconst, attrs)
    attrkeys = self.params.cfg["ATTRKEYS"].split()
    attrvals = self.params.cfg["ATTRVALS"].split()

    attrkeys.extend(["time", "gravconst"])
    attrvals.extend([str(t0), str(G)])

    for key,val in zip(attrkeys, attrvals):
      cmds.append("h5part_writeattr -i " + self.dir + "initial.h5part " +\
          " -k " + key + " -v " + val)

    # prepare and copy binary
    if not self.params.cfg.has_key("BINFILE"):
      srcdir = resolvePath(self.params.cfg["SRCDIR"]) + "/"
      targ = self.params.cfg["MAKETARG"]
      oldwd = os.getcwd()
      os.chdir(srcdir)
      #(stat, out) = commands.getstatusoutput("make " + targ)
      os.chdir(oldwd)
      self.params.cfg["BINFILE"] = \
          resolvePath(srcdir + self.params.cfg["BINARY"])
    binf = self.params.cfg["BINFILE"]
    cmds.append("cp -Rpv " + binf + " " + self.dir)

    # now execute the commands
    for cmd in cmds:
      (stat, out) = commands.getstatusoutput(cmd)
      self.Log.write(out)
      if not stat == 0:
        self.setError(cmd + " failed")
        return self.state

    self.next = self._submit
    self.state = "prepared"
    return self.state

  def _submit(self):
    if not self.getState() == "prepared":
      self.setError("tried to submit from a non-prepared state")
      return self.state
    
    subcmd = self.params.cfg["SUBCMD"]
    nocpus = self.params.cfg["NOCPUS"]
    binary = self.params.cfg["BINARY"]
    runarg = self.params.cfg["RUNARGS"]
    sstnam = self.params.cfg["SIMSETNAME"]
    
    stoptime = 86400
    savetime = 3600

    subcmd = subcmd.replace('$NOCPUS' , str(nocpus))
    subcmd = subcmd.replace('$SIMNAME', sstnam + "_" + self.params.key)
    subcmd = subcmd.replace('$BINARY' , binary)
    
    subcmd = subcmd.replace('$SAVETIME' , str(savetime))
    subcmd = subcmd.replace('$STOPTIME' , str(stoptime))
    subcmd = subcmd.replace('$RUNARGS' ,  runarg)

    oldwd = os.getcwd()
    os.chdir(self.dir)
    (stat, out) = commands.getstatusoutput(subcmd)
    stat = 0
    print out
    os.chdir(oldwd)

    if not stat == 0:
      self.setError("submit command failed (" + subcmd + ")")
    else:
      self.next = self._doNothing
      self.state = "submitted"

    return self.state

    
  def setSGEstat(self, str):
    if str == "queued":
      self.next = self._doNothing
      self.state = str

    if str == "run":
      self.next = self._postproc
      self.state = str

    if str == "finished":
      self.next = self._finalize
      self.state = str

  def _postproc(self):
    if self.state == "run":
      print "postproc!"
      #check for new files
    else if self.state == "finished":
      self.next = self._doNothing
    else:
      self.setError("tried to posrproc from a non-submitted state")
    return self.state

  def _finalize(self):
    if not self.state == "finished":
      self.setError("tried to posrproc from a non-submitted state")
      return self.state



  def _findBodies(self):
    nan = float('nan')

    toler = float(self.params.cfg["BODTOLERANCE"])
    boddb = shelve.open(self.params.cfg["BODIESDB"])

    # find target candidates
    targtmpl = BodyFile("", Me*self.params.mtar, nan, \
        self.params.temp, nan, {})
    targcand = []
    for key in boddb.keys():
      if boddb[key].isSimilar(targtmpl, toler):
        targcand.append(boddb[key])

    # find impactor candidates
    impatmpl = BodyFile("", Me*self.params.mimp, nan, \
        self.params.temp, nan, {})
    impacand = []
    for key in boddb.keys():
      if boddb[key].isSimilar(impatmpl, toler):
        impacand.append(boddb[key])
    
    boddb.close()
    
    # find fitting pairs
    pairs = []
    for tarbd in targcand:
      for impbd in impacand:
        if abs((impbd.h - tarbd.h)/tarbd.h) < toler:
          pairs.append( (tarbd, impbd) )

    if len(pairs) == 0:
      pairs.append( (None, None) )
    else:
      # sort according to smoothing length
      pairs.sort(key=lambda bod: bod[0].h)
    
    return pairs[0]


class SimParams(object):
  def __init__(self, mimp, mtar, impa, vimprel, cfg):
    self.mimp = mimp
    self.mtar = mtar
    self.impa = impa
    self.vimprel = vimprel

    self.cfg = cfg

    self.temp = float(cfg["TEMP"])

    self.key = 'mtar%07.3f' % mtar + '_' \
      + 'mimp%07.3f' % mimp + '_' \
      + 'impa%04.1f' % impa + '_' \
      + 'vimp%04.1f' % vimprel


class SimSet(object):
  sims = {}
  def __init__(self, cfgfile):
    self.simsetcfg = {}
    
    if path.exists(cfgfile):
      execfile(cfgfile, self.simsetcfg)
      print cfgfile + " loaded"
    else:
      print cfgfile + " not found!   exiting ..."
    
    bdir = self.simsetcfg["SIMSETBDIR"]
    if not path.exists(bdir):
      os.mkdir(bdir)

    logfile = self.simsetcfg["LOGFILE"]
    self.Log = Logger(logfile)
    self.Log.write("new SimSet")
    
    mtara = np.array( self.simsetcfg["MTAR"].split(), dtype=float)
    mimpa = np.array( self.simsetcfg["MIMP"].split(), dtype=float)
    vimprela = np.array( self.simsetcfg["VIMPREL"].split(), dtype=float)
    impaa = np.array( self.simsetcfg["IMPA"].split(), dtype=float)

    cond = self.simsetcfg["SIMCOND"]
    if len(cond) < 1:
      cond = "True"

    for mtar in mtara:
      for mimp in mimpa:
        for vimprel in vimprela:
          for impa in impaa:
            if ( eval(cond) ):
              simpar = SimParams(mimp, mtar, impa, vimprel, self.simsetcfg)
              if not self.sims.has_key(simpar.key):
                self.sims[simpar.key] = Simulation( simpar )


  def __del__(self):
    pass

