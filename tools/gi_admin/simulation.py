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
from h5part import H5PartDump
from clumps import SimClumps, ClumpsPlotConfig

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


class Simulation(object):
  # noCPUs, 
  def __init__(self,params,cfg):
    self.params = params
    self.cfg    = cfg
    self.params.temp = cfg.T

    ssbdir = resolvePath(cfg.dir) + "/"
    self.dir = resolvePath(ssbdir + params.key) + "/"
    
    if path.exists(self.dir):
      self.Log = Logger(self.dir + "setup.log")
    else:
      self.Log = Logger("/dev/null")

    self.jobid = 0
    self.nodumps = 0

    self.state = "unknown"
    self._setOrbit()
    self.getState()

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
      neededfiles = self.cfg.auxfiles.split()
      neededfiles.append(self.cfg.binary)
      neededfiles.append("initial.h5part")

      for file in neededfiles:
        if not os.path.exists(self.dir + os.path.split(file)[1]):
          self.state = "unprepared"
          self.next = self._prepare
          return self.state
      
      print "find dumps for ",self.params.key
      self._findDumps()

      # so everything's ready
      if self.jobid == 0:
        if self.nodumps == 0:
          self.state = "prepared" 
          self.next = self._submit
        elif self.nodumps == self.totdumps:
          self.state = "finished"
          self.next = self._doNothing
        else:
          self.state = "partial"
          self.next = self._doNothing
      else:
        self.state = "run"
        self.next = self._postproc
      return self.state
    
    
    self.state = "unknown"
    self.next = self._doNothing
    return self.state

  def reset(self):
    self.state = "unprepared"
    self.next = self._prepare
    return self.state

  def _doNothing(self):
    return self.state
  
  def _setOrbit(self):
    # find fitting bodies
    (self.tarb, self.impb) = self._findBodies()

    mtar = self.tarb.m
    mimp = self.impb.m
    Rtar = self.tarb.r
    Rimp = self.impb.r
    
    G = self.cfg.G
    relsep = self.cfg.relsep
    
    # get giant impact
    gi = GiantImpact(mtar, mimp, Rtar, Rimp, G)
    vimp = self.params.vimprel * gi.vesc
    impa = self.params.impa * deg2rad

    (self.r0tar, self.r0imp, self.v0tar, self.v0imp, t0, self.gilogstr) = \
        gi.getInitVimpAlpha(vimp, impa, relsep)

    tcol = 2*(Rtar + Rimp)/gi.vesc
    self.gi = gi
    self.tcol = tcol
    self.vimp = vimp
    self.G    = G
    
    self.tstop = tcol * self.cfg.tcolmax
    self.tdump = tcol * self.cfg.tcoldump
    self.tanim = tcol * self.cfg.tcolanim

    self.t0 = t0
    
    self.totdumps = int( ( self.tstop - self.t0 ) / self.tdump )

    self.mtot = mtar + mimp
    self.nop  = self.tarb.nop + self.impb.nop


  def _findDumps(self):
    dumpdbfile = self.dir + "dumps"
    dumps = shelve.open(dumpdbfile)
    
    #for file in os.listdir(self.dir):
    #  if file[0:4] == "dump" and file[-6:] == "h5part":
    #    if not dumps.has_key(file):
    #      cmd = "h5part_readattr -k time -i " + self.dir + file
    #      (stat, out) = commands.getstatusoutput(cmd)
    #      time = float("nan")
    #      if stat == 0:
    #        time = float(out.split()[1])
    #        dumps[file] = time
  
    self.dumps = []
    for key in dumps.keys():
      self.dumps.append( (key, dumps[key]) )
    self.nodumps = len(self.dumps)
    
    dumps.sync()
    dumps.close()

    self.dumps.sort(key=lambda dmp: dmp[0])


  def _delAnimDumps(self):
    dumps = self.dumps
    tanim = self.tanim
    tdump = self.tdump
    delta = 5.e-3
    
    os.remove(self.dir + "dumps")
    self._findDumps()
    
    for (dfile, dtime) in self.dumps:
      ranim = dtime / tanim
      rdump = dtime / tdump

      deltaanim = abs( ranim - round(ranim) )
      deltadump = abs( rdump - round(rdump) )

      if deltaanim < delta and deltadump > delta:
        print "delete  ", dfile, dtime
        os.remove(self.dir + dfile)


  def _prepare(self):
    # if it does not yet exist, make dir
    if not path.exists(self.dir):
      os.mkdir(self.dir)
      self.Log = Logger(self.dir + "setup.log")

    # copy files
    auxf = self.cfg.auxfiles.split()
    
    for file in auxf:
      shutil.copy2(resolvePath(file), self.dir)
    
    self.Log.write("auxiliary files copied")
    
    if path.exists(self.dir + "initial.h5part"):
      os.remove(self.dir + "initial.h5part")

    # get orbit
    gilog = open(self.dir + "gi_setup.log", "w")
    print >>gilog, self.gilogstr
    gilog.close()
    self.Log.write("giant impact calculated")

    # displace bodies
    r0tar = self.r0tar
    r0imp = self.r0imp
    v0tar = self.v0tar
    v0imp = self.v0imp

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
    attrs = self.cfg.attr
    attrs.append( ("time",      self.t0   ) )
    attrs.append( ("gravconst", self.G    ) )
    attrs.append( ("tcol",     self.tcol ) )

    for (key,val) in attrs:
      cmds.append("h5part_writeattr -i " + self.dir + "initial.h5part " +\
          " -k " + key + " -v " + str(val))

    # prepare and copy binary
    binfile = resolvePath( self.cfg.srcdir ) + "/" + self.cfg.binary
    if not path.exists(binfile):
      oldwd = os.getcwd()
      os.chdir(self.cfg.srcdir)
      (stat, out) = commands.getstatusoutput("make " + self.cfg.maketarg)
      os.chdir(oldwd)
      self.Log.write("binary compiled")
    cmds.append("cp -Rpv " + binfile + " " + self.dir)

    # now execute the commands
    for cmd in cmds:
      (stat, out) = commands.getstatusoutput(cmd)
      self.Log.write(out)
      if not stat == 0:
        self.setError(cmd + " failed")
        return self.state
    
    self.next = self._submit
    self.Log.write("prepared state")
    self.state = "prepared"
    return self.state

  def _submit(self):
    if not self.getState() == "prepared":
      self.setError("tried to submit from a non-prepared state")
      return self.state
    
    subcmd = self.cfg.subcmd
    nocpus = self.cfg.nocpus
    binary = self.cfg.binary
    runarg = self.cfg.runargs
    sstnam = self.cfg.name

    self.jobname = sstnam + "_" + self.params.key
  
    self.Log.write("tcol = %9.3e" % self.tcol)
    self.Log.write("tstop = %9.3e" % self.tstop)

    subcmd = subcmd.replace('$NOCPUS' , str(nocpus))
    subcmd = subcmd.replace('$SIMNAME', self.jobname)
    subcmd = subcmd.replace('$BINARY' , binary)
    
    subcmd = subcmd.replace('$SAVETIME' , str(self.tanim))
    subcmd = subcmd.replace('$STOPTIME' , str(self.tstop))
    subcmd = subcmd.replace('$RUNARGS' ,  runarg)

    oldwd = os.getcwd()
    os.chdir(self.dir)
    (stat, out) = commands.getstatusoutput(subcmd)
    os.chdir(oldwd)

    if not stat == 0:
      self.setError("submit command failed (" + out + ")")
    else:
      self.next = self._doNothing
      self.state = "submitted"
      self.Log.write("job \"" + self.jobname + "\" submitted")
    return self.state

  def _abort(self):
    if not self.jobid == 0:
      delcmd = self.cfg.delcmd
      delcmd = delcmd.replace('$JOBID', str(self.jobid))
      (dstat, dout) = commands.getstatusoutput(delcmd)
      print dout

  def _del(self):
    if not self.state == "run":
      for root, dirs, files in os.walk(self.dir, topdown=False):
        for fname in files:
          os.remove(os.path.join(root, fname))
        for dname in dirs:
          os.rmdir(os.path.join(root, dname))
      os.rmdir(self.dir)
      self.state = "unprepared"
      self.next = self._prepare


  def _continue(self):
    pass

  def setSGEjobid(self, jobid):
    self.jobid = jobid

  def setSGEstat(self, str):
    if self.state == "error":
      return
      
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
    if self.state == "run" or self.state == "partial":
      self._findDumps()
      self._findClumps()

    elif self.state == "finished":
      self.Log.write("SGE reported finished state, was job " + str(self.id))
      self.next = self._doNothing
    
    else:
      self.setError("tried to postproc from a non-submitted state")
    
    return self.state

  def _finalize(self):
    if not self.state == "finished":
      self.setError("tried to finalize from non-finished state")
    else:
      self.next = self._doNothing
      self.state = "finalized"
      self.Log.write("finalized simulation")

    return self.state

  def _findBodies(self):
    nan = float('nan')
    toler = self.cfg.bodtol
    boddb = shelve.open(resolvePath(self.cfg.bodiesdb))

    # find target candidates
    targtmpl = BodyFile("", Me*self.params.mtar, nan, \
        self.params.temp, nan, nan, {})
    targcand = []
    for key in boddb.keys():
      if boddb[key].isSimilar(targtmpl, toler):
        targcand.append(boddb[key])

    # find impactor candidates
    impatmpl = BodyFile("", Me*self.params.mimp, nan, \
        self.params.temp, nan, nan, {})
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

  # change that
  def _clumpsExtent(self):
    import tables as pt
    clpdir  = self.dir + "clumps/"

    inf = float("inf")
    clpmin = np.array([ inf,  inf,  inf])
    clpmax = np.array([-inf, -inf, -inf])

    for file in os.listdir(clpdir):
      clpsh = pt.openFile(clpdir + file, "r")
      clps  = clpsh.root.current
      for clppos in clps.pos[:]:
        clpmin = np.minimum( clpmin, clppos )
        clpmax = np.maximum( clpmax, clppos )
      clpsh.close()

    self.clpmin = clpmin
    self.clpmax = clpmax

    return (clpmin, clpmax) 

  def _loadClumps(self):
    if path.exists(self.dir + "clumps.h5part"):
      self.clmps = SimClumps(self.dir + "clumps.h5part")
      self.clmps.plotMass(self.dir + "clumps.pdf", self.cfg.cpconf, self)




class SimParam(object):
  def __init__(self, mimp, mtar, impa, vimprel):
    self.mimp = mimp
    self.mtar = mtar
    self.impa = impa
    self.vimprel = vimprel

    self.key = 'mtar%07.3f' % mtar + '_' \
      + 'mimp%07.3f' % mimp + '_' \
      + 'impa%04.1f' % impa + '_' \
      + 'vimp%04.2f' % vimprel
  
  def isSimilar(self, par, tol):
    fitmimp = not abs( (self.mimp - par.mimp ) / self.mimp ) > tol
    fitmtar = not abs( (self.mtar - par.mtar ) / self.mtar ) > tol
    fitimpa = not abs( (self.impa - par.impa ) / self.impa ) > tol
    fitvimp = not abs( (self.vimprel - par.vimprel ) / self.vimprel ) > tol
    return fitmimp and fitmtar and fitimpa and fitvimp


class SimSet(object):
  sims = {}
  def __init__(self, cfg):
    self.cfg = cfg
    
    bdir = self.cfg.dir
    if not path.exists(bdir):
      os.mkdir(bdir)

    mtara = self.cfg.mtar
    mimpa = self.cfg.mimp
    vimprela = self.cfg.vimp
    impaa = self.cfg.impa

    cond = self.cfg.paramcond
    if len(cond) < 1:
      cond = "True"

    for mtar in mtara:
      for mimp in mimpa:
        for vimprel in vimprela:
          for impa in impaa:
            if ( eval(cond) ):
              simpar = SimParam(mimp, mtar, impa, vimprel)
              if not self.sims.has_key(simpar.key):
                self.sims[simpar.key] = Simulation( simpar, self.cfg)


class SimSetConfig(object):
  def __init__(self):
    self.dir = ""
    self.name = ""
    self.bodiesdb = ""

    self.G = G
    self.T = 0.
    self.mimp = np.array([])
    self.mtar = np.array([])
    self.vimp = np.array([])
    self.impa = np.array([])
    self.paramcond = "(mimp <= mtar)"

    self.bodtol = 0.1
    self.relsep = 3.0

    self.attr = []
    self.auxfiles = ""
    self.srcdir   = ""
    self.maketarg = ""
    self.binary   = ""

    self.subcmd  = ""
    self.delcmd  = ""
    self.nocpus  = ""
    self.runarg  = ""

    self.tcolmin  =   5.0
    self.tcolmax  =  50.
    self.tcoldump =   0.5
    self.tcolanim =   0.05
    
    self.stopnoclumps     =   2
    self.stopdnopdtcol    =  50
    self.stopdnopdtcolesc =  50

    self.cpconf = ClumpsPlotConfig()

