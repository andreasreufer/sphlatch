import os
import os.path as path
import stat
import time
import commands
import shelve
import numpy as np
import shutil
import re

from numpy import sqrt, sin, cos, arccos, log, abs, tan, arctan2, pi, zeros
from body import BodyFile
from setup_gi import GiantImpact
from h5part import H5PartDump
from clumps import SimClumps, ClumpsPlotConfig, MoonsPlotConfig
from disk import getDisk

import pdb
import sim_machine

rad2deg = 360./(2.*pi)
deg2rad = 1./rad2deg

# constants to transform the initial parameters to physical values
G   = 6.67429e-8
Re  = 6.37814e8
Me  = 5.97360e27
lem = 3.40000e41 # Canup (2001) value
#lem = 2.88993e41 # actual (?) value
nan = float('nan')

np.seterr(all='ignore')

def resolvePath(_path):
  return path.abspath( path.expanduser(_path) )

class Logger(object):

  def __init__(self, filename):
    self.filename = filename

  def __del__(self):
    pass

  def write(self, str):
    if len(str) > 0:
      logfile = open(self.filename, 'a')
      logfile.write(time.ctime() + ":   " + str + "\n")
      logfile.close()


class Simulation(object):
  # noCPUs, 
  def __init__(self,params,cfg):
    self.params = params
    self.cfg    = cfg
    self.params.temp = cfg.T
    self.idfilt = cfg.idfilt
    self.keepfiles = cfg.keepfiles.split()

    ssbdir = resolvePath(cfg.dir) + "/"
    self.dir = resolvePath(ssbdir + params.key) + "/"
    
    if path.exists(self.dir):
      self.Log = Logger(self.dir + "setup.log")
    else:
      self.Log = Logger("/dev/null")

    self.jobid = 0
    self.nodumps  = 0
    self.totdumps = 0

    self.state = "unknown"
    self._setOrbit()
    self.getState()
    self.results = SimResults()
    self.clumpjobid = 0

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
      
      print "find dumps for ",self.params.key
      self._findDumps()

      for file in neededfiles:
        if not os.path.exists(self.dir + os.path.split(file)[1]):
          self.state = "unprepared"
          self.next = self._prepare
          return self.state

      # so everything's ready
      if self.jobid == 0:
        if self.nodumps == 0:
          self.state = "prepared" 
          self.next = self._submit
        else:
          (dfile, dtime, mtime) = self.dumps[-1]
	  if dtime >= 0.99*self.tstop:
	    self.state = "finished"
	  else:
	    self.state = "partial"
      else:
        # there is a job id, so state was set by SimAdmin
	if self.state == "run":
          self.next = self._postproc
          if self._isSlow(5.):
            self.state = "slow"

	  if self._isLogStuck(1800.):
	    self.state = "stuck"
      
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
    if (self.tarb, self.impb) == (None, None):
      print "no suitable bodies found!"
      self.setError("no suitable bodies found!")
      return

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

    tcol = 2*(Rtar + Rimp)/vimp
    self.gi = gi
    self.tcol = tcol
    self.vimp = vimp
    self.G    = G
    self.mtot = mtar + mimp
    
    if hasattr(self.tarb, "rmat"):
      Rtarmat = self.tarb.rmat
      Rimpmat = self.impb.rmat

      ratmat = ( Rtar + Rimp ) / ( Rtarmat + Rimpmat )
      self.impamat = rad2deg*np.arcsin( ratmat*np.sin(impa) )
    
    self.impa    = rad2deg*impa
    
    self.tstop = tcol * self.cfg.tcolmax
    self.tdump = tcol * self.cfg.tcoldump
    self.tanim = tcol * self.cfg.tcolanim
    
    # allow a 1000fold maximum density compared to mean value
    self.hmin = 0.1*min(self.tarb.h, self.impb.h)
    self.pmin = 0.
    if self.mtot < (1.e-3 * Me):
      self.pmin = -1.e99

    self.t0 = t0
    
    self.totdumps = int( ( self.tstop - self.t0 ) / self.tdump )

    self.nop  = self.tarb.nop + self.impb.nop

  def _findDumps(self):
    dumps = {}
    for file in os.listdir(self.dir):
      if file[0:4] == "dump" and file[-6:] == "h5part" and len(file) == 31:
        dumps[file] = float(file[13:24])
  
    self.dumps = []
    for key in dumps.keys():
      self.dumps.append( (key, dumps[key], os.path.getmtime(self.dir + key)) )
    self.nodumps = len(self.dumps)
    self.dumps.sort(key=lambda dmp: dmp[0])

    if len(self.dumps) > 0:
      (lfile, self.tlast, lmtime) = self.dumps[-1]
    else:
      self.tlast = nan
    

  def _isSlow(self,k):
    nod = self.nodumps
    if nod > 1:
      mtimes = np.zeros(nod)
      for i in range(nod):
         (dfile, dtime, mtimes[i]) = self.dumps[i]
      
      self.dumpsdmtime = mtimes[1:] - mtimes[0:-1]
      return ( time.time() - mtimes[-1] ) > k*np.median( self.dumpsdmtime )
    return False

  def _isLogStuck(self,dt):
    return ( time.time() - os.path.getmtime(self.dir + "logDomain000") ) > dt

  def _delAnimDumps(self):
    dumps = self.dumps
    tanim = self.tanim
    tdump = self.tdump
    delta = 5.e-3
    
    self._findDumps()
    
    for (dfile, dtime, dmtime) in self.dumps:
      ranim = dtime / tanim
      rdump = dtime / tdump

      deltaanim = abs( ranim - round(ranim) )
      deltadump = abs( rdump - round(rdump) )

      if deltaanim < delta and deltadump > delta:
        print "delete  ", dfile, dtime
        os.remove(self.dir + dfile)

  def _delEscapeeOrphans(self):
    for file in os.listdir(self.dir):
      if "_esc.h5part" in file and \
      not os.path.exists(self.dir + file.replace("_esc","")):
	print "delete ",self.dir+file
        (exstat, out) = commands.getstatusoutput("rm "+self.dir+file)

  def _delDumpsSmaller(self,nob):
    for (dfile, dtime, dmtime) in self.dumps:
      if os.path.getsize(self.dir + dfile) < nob:
        print "delete  ", self.dir+dfile
	(exstat, out) = commands.getstatusoutput("rm "+self.dir+dfile)

    self._delEscapeeOrphans()

  def _delLastDumps(self,nod=5):
    for i in range(-1,-nod-1,-1):
      (dfile, dtime, dmtime) = self.dumps[i]
      dfile = self.dir + dfile
      efile = dfile.replace(".h5part","") + "_esc.h5part"
      
      print "delete  ", dfile
      (exstat, out) = commands.getstatusoutput("rm " + dfile)
      if os.path.exists(efile):
        print "delete  ", efile
        (exstat, out) = commands.getstatusoutput("rm " + efile)


  def _delAll(self):
    (exstat, out) = commands.getstatusoutput("rm -rf "+self.dir)
    self.state = "unprepared"
    self.next  = self._prepare
    self.nodumps = 0


  def _keepNthDumps(self,nth):
    dumps = self.dumps
    
    self._findDumps()
    
    for n in range( len(dumps) ):
       (dfile, dtime, dmtime) = dumps[n]
       if not ( n % nth == 0 ):
	  print "delete  ", dfile, dtime
	  os.remove(self.dir + dfile)
       else:
          print "keep    ", dfile, dtime

  def _prepare(self):
    # if it does not yet exist, make dir
    if not path.exists(self.dir):
      os.mkdir(self.dir)
      self.Log = Logger(self.dir + "setup.log")

    # copy files
    auxf = self.cfg.auxfiles.split()
    
    cmds = []
    for file in auxf:
      #shutil.copy2(resolvePath(file), self.dir)
      cmd = "cp -Rpv " + resolvePath(file) + " " + self.dir
      cmds.append("cp -Rpv " + resolvePath(file) + " " + self.dir )
    
    #self.Log.write("auxiliary files copied")
    
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

    cmds.append("cp -Rpv " + self.tarb.file + " " + self.dir + "tarb.h5part")
    cmds.append("h5part_displace -i " + self.dir + "tarb.h5part " +\
        "--pos [%e,%e,%e] " % (r0tar[0], r0tar[1], r0tar[2]) +\
        "--vel [%e,%e,%e] " % (v0tar[0], v0tar[1], v0tar[2]) )
    cmds.append("cp -Rpv " + self.impb.file + " " + self.dir + "impb.h5part")
    cmds.append("h5part_displace -i " + self.dir + "impb.h5part " +\
        "--pos [%e,%e,%e] " % (r0imp[0], r0imp[1], r0imp[2]) +\
        "--vel [%e,%e,%e] " % (v0imp[0], v0imp[1], v0imp[2]) +\
        "--id %d" % self.cfg.idfilt )
    cmds.append("h5part_combineA_ " + self.dir + "tarb.h5part " +\
        self.dir + "impb.h5part " + self.dir + "initial.h5part")
    cmds.append("rm " + self.dir + "tarb.h5part " + self.dir + "impb.h5part")
    
    # write attributes (gravconst, attrs)
    attrs = []
    attrs.append( ("time",      self.t0   ) )
    attrs.append( ("gravconst", self.G    ) )
    attrs.append( ("tcol",     self.tcol ) )
    attrs.append( ("idfilt",   self.cfg.idfilt ) )
    attrs.append( ("rmaxbound", self.vimp*self.tstop) )
    attrs.append( ("rmaxunbound", self.vimp*self.tstop) )
    attrs.append( ("hmin", self.hmin) )
    attrs.append( ("pmin", self.pmin) )

    for (key,val) in attrs:
      cmds.append("h5part_writeattr -i " + self.dir + "initial.h5part " +\
          " -k " + key + " -v " + str(val))
    
    for (key,val) in self.cfg.attr:
      cmds.append("h5part_writeattr -i " + self.dir + "initial.h5part " +\
          " -k " + key + " -v " + str(val))

    cmds.append("h5part_writeattr -i " + self.dir + "initial.h5part " +\
        " -k " + self.params.key + " -v 0. ")

    # prepare and copy binary
    binfile = resolvePath( self.cfg.srcdir ) + "/" + self.cfg.binary
    if not path.exists(binfile):
      oldwd = os.getcwd()
      os.chdir(self.cfg.srcdir)
      (exstat, out) = commands.getstatusoutput("make " + self.cfg.maketarg)
      os.chdir(oldwd)
      self.Log.write("binary compiled")
    cmds.append("cp -Rpv " + binfile + " " + self.dir)

    # now execute the commands
    for cmd in cmds:
      (exstat, out) = commands.getstatusoutput(cmd)
      self.Log.write(out)
      if not exstat == 0:
        self.setError(cmd + " failed")
        return self.state
    
    self.next = self._submit
    self.Log.write("prepared state")
    self.state = "prepared"
    return self.state
  
  def _setInitialCourant(self,courant):
    self.courant = courant
    (dfile,dtime,dmtime) = self.dumps[-1]
    cmd = "h5part_writeattr -i " + self.dir + "initial.h5part" +\
          " -k courant -v " + str(courant)
    (exstat, out) = commands.getstatusoutput(cmd)

  def _grabBinary(self):
    # prepare and copy binary
    binfile = resolvePath( self.cfg.srcdir ) + "/" + self.cfg.binary
    if not path.exists(binfile):
      return
    cmd = "cp -Rpv " + binfile + " " + self.dir
    (exstat, out) = commands.getstatusoutput(cmd)

  def _submit(self):
    if not ( self.state == "prepared" or self.state == "partial"):
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

    oldwd = os.getcwd()
    os.chdir(self.dir)
    
    if hasattr(sim_machine, "qsub_ompscr"):
      script = sim_machine.qsub_ompscr
      script = script.replace('$NOCPUS' , str(nocpus))
      script = script.replace('$SIMNAME', self.jobname)
      script = script.replace('$BINARY' , binary)
    
      script = script.replace('$SAVETIME' , str(self.tanim))
      script = script.replace('$STOPTIME' , str(self.tstop))
      script = script.replace('$RUNARGS' ,  runarg)

      scrfile = open("jobscript.sh", 'w')
      scrfile.write(script)
      scrfile.close()

      subcmd = subcmd.replace('$SCRIPT' , "jobscript.sh")
    else:
      subcmd = subcmd.replace('$NOCPUS' , str(nocpus))
      subcmd = subcmd.replace('$SIMNAME', self.jobname)
      subcmd = subcmd.replace('$BINARY' , binary)
    
      subcmd = subcmd.replace('$SAVETIME' , str(self.tanim))
      subcmd = subcmd.replace('$STOPTIME' , str(self.tstop))
      subcmd = subcmd.replace('$RUNARGS' ,  runarg)
      
    print "cd " + self.dir + ";" + subcmd
    (exstat, out) = commands.getstatusoutput(subcmd)
    os.chdir(oldwd)
      
    self.Log.write(subcmd)

    if not exstat == 0:
      self.setError("submit command failed (" + out + ")")
    else:
      self.next = self._doNothing
      self.state = "submitted"
      self.Log.write("job \"" + self.jobname + "\" submitted")
    return self.state
  
  def _extend(self, tcolmax):
    self.tstop = self.tcol * tcolmax
    self._resubmit()

  def _setCourant(self,courant):
    self.courant = courant
    dfile = "initial.h5part"

    if self.nodumps > 0:
      (dfile,dtime,dmtime) = self.dumps[-1]
    
    cmd = "h5part_writeattr -i " + self.dir + dfile +\
          " -k courant -v " + str(courant)
    if path.exists(self.dir + dfile):
      (exstat, out) = commands.getstatusoutput(cmd)
      print "set courant attribute in file " + dfile + " to ",str(courant)

  def _setAttribute(self,key,val):
    (dfile,dtime,dmtime) = self.dumps[-1]
    cmd = "h5part_writeattr -i " + self.dir + dfile +\
          " -k " + str(key) + " -v " + str(val)
    (exstat, out) = commands.getstatusoutput(cmd)
  
  def _resubmit(self):
    sdir = self.dir
    (dfile,dtime,dmtime) = self.dumps[-1]
    efile = dfile.replace(".h5part","_esc.h5part")
    
    cmd = "mv " + sdir + dfile + " " + sdir + "/initial.h5part"
    (exstat, out) = commands.getstatusoutput(cmd)
    
    cmd = "mv " + sdir + efile + " " + sdir + "/initial_esc.h5part"
    (exstat, out) = commands.getstatusoutput(cmd)
  
    self.state = "prepared"
    self._submit()
    self._findDumps()
  
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
      #os.rmdir(self.dir)
      self.state = "unprepared"
      self.next = self._prepare


  def _continue(self):
    pass

  def set_jobid(self, jobid):
    self.jobid = jobid

  def set_state(self, str):
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

  def _getArchiveFiles(self):
    archfiles = []
    self._findDumps()
    if self.nodumps > 0:
      (lfile, ltime, mtime) = self.dumps[-1]
      archfiles.append( self.dir + lfile )

      keepdumps = self.cfg.keepdumps
      for trel in keepdumps:
        archfiles.extend( self._getDump( trel ) )

      escfiles = []
      for afile in archfiles:
      	aesc = afile.replace(".h5part","_esc.h5part")
	if path.exists(aesc):
          escfiles.append( aesc )
      archfiles.extend( escfiles )
    
    for f in self.keepfiles:
      if path.exists(self.dir + f):
        archfiles.append( self.dir + f )
    
    return archfiles

  def _copyAllDumps(self, targ):
    sdir = self.dir
    tdir = targ + "/" + self.params.key + "/"
    
    if not path.exists(tdir):
      os.mkdir(tdir)
      print tdir

    for (dfile, dtime, mtime) in self.dumps:
      cmd = "cp -Rpv " + sdir + dfile + " " + tdir + dfile
      (exstat, out) = commands.getstatusoutput(cmd)
      print out

  def _getDump(self, trel):
    eps = 1.e-3
    mtchdumps = []
    for (dfile, dtime, dmtime) in self.dumps:
      if ( abs(dtime / self.tcol - trel) < eps ):
        mtchdumps.append(self.dir + dfile)
    return mtchdumps

  def _findBodies(self):
    toler = self.cfg.bodtol
    htoler= self.cfg.smotol
    boddb = shelve.open(resolvePath(self.cfg.bodiesdb))

    # find target candidates
    mmatdummy = np.zeros([32])
    targtmpl = BodyFile("", Me*self.params.mtar, mmatdummy, nan, \
        self.params.temp, nan, nan, {})
    targcand = []
    for key in boddb.keys():
      if boddb[key].isSimilar(targtmpl, toler):
        targcand.append(boddb[key])

    # find impactor candidates
    impatmpl = BodyFile("", Me*self.params.mimp, mmatdummy, nan, \
        self.params.temp, nan, nan, {})
    impacand = []
    for key in boddb.keys():
      if boddb[key].isSimilar(impatmpl, toler):
        impacand.append(boddb[key])

    boddb.close()
    
    # find fitting pairs by looking at mass ratio and smoothing length
    pairs = []
    gammaparm = (self.params.mimp / self.params.mtar )
    for tarbd in targcand:
      for impbd in impacand:
        hdev = abs((impbd.h - tarbd.h)/tarbd.h)
	gdev = abs((gammaparm - (impbd.m/tarbd.m))/gammaparm)
	if gdev < toler and hdev < htoler:
          pairs.append( (tarbd, impbd) )
    
    if len(pairs) == 0:
       return (None, None) 
    else:
       (btarb, bimpb) = pairs[0]
       desnop = self.cfg.desnop
       for (tarb, impb) in pairs:
         cfit = abs( float( desnop - ( tarb.nop +   impb.nop) )/desnop )
         bfit = abs( float( desnop - (btarb.nop +  bimpb.nop) )/desnop )
	 if ( cfit < bfit ):
           btarb = tarb
           bimpb = impb

       return (btarb, bimpb)

  def _plotMasses(self):
    if path.exists(self.dir + "clumps.h5part"):
      self.clmps = SimClumps(self)
      self.clmps.plotMass(self.dir + "clumps.pdf", self.cfg.cpconf)

  def _plotMoons(self):
    self.clmps = SimClumps(self, loadsat=True)
    self.clmps.plotLRmoons(self.dir + "moons.pdf", self.cfg.mpconf)

  def _plotDynamics(self):
    if path.exists(self.dir + "clumps.h5part"):
      self.clmps = SimClumps(self)
      self.clmps.plotDynamics(self.dir + "dynamics.pdf", self.cfg.cpconf)
      print self.dir + "dynamics.pdf"

  def _getDynamics(self):
    if path.exists(self.dir + "clumps.h5part"):
      self.clmps = SimClumps(self, loadeng=True)
      self.clmps.getDynamics2D()
      return True
    else:
      return False

  def _getResults(self,tmin=0.90,tmax=0.99):
    if path.exists(self.dir + "clumps.h5part"):
      try:
        self.clmps = SimClumps(self, loadeng=True)
        self.clmps.getResults(tmin, tmax)
        if self.nodumps > 0:
          (lfile, ltime, mtime) = self.dumps[-1]
          getDisk(self.dir + lfile, self)

      except:
        self.results.problem = True

  def _plotSummary(self):
    if path.exists(self.dir + "clumps.h5part"):
      try:
        self.clmps = SimClumps(self, loadeng=True)
        self.clmps.plotSummary(self.dir + "summary.pdf")
      
      except:
        self.results.problem = True
  
  def _getClumpEnergy(self,tmin=0.90,tmax=0.99):
    res = self.results
    try:
      clmps = SimClumps(self, loadclmp=True)
      res.Epotclmp = clmps.Epotclmp[-1,:]
      res.Erotclmp = clmps.Erotclmp[-1,:]
      res.Iclmp    = clmps.Iclmp[-1,:]
        
      print - res.Erotclmp[1] / res.Epotclmp[1], \
          res.mm[0] / clmps.m[-1,0], res.mm[1] / clmps.m[-1,1]
      print "Epotclmp and Erotclmp loaded!"
      
    except:
      noc = 0
      if len( res.mm ) > 0:
        noc = res.mm.shape[0]
          
      res.Epotclmp = nan*np.zeros(noc)
      res.Epotclmp = nan*np.zeros(noc)
      res.Iclmp    = nan*np.zeros(noc)
      print "Epotclmp and Erotclmp failed!"

  def _totEnergy(self):
    sdir  = self.dir
    for (dfile,dtime,dmtime) in self.dumps:
      cmd = "totEnergy__ " + sdir + dfile + " " + sdir + "clumps.h5part"
      (exstat, out) = commands.getstatusoutput(cmd)
      print out

  def _lastLog(self, n=10):
    sdir  = self.dir
    cmd = "tail -n " + str(n) + " " + sdir + "logDomain000"
    (exstat, out) = commands.getstatusoutput(cmd)
    return out

  def _lastDt(self):
    sdir  = self.dir
    cmd = "tail -n 100 " + sdir + "logDomain000 | grep dt:|tail -n1"
    (exstat, out) = commands.getstatusoutput(cmd)
    return out[24:]

  def _hdf5OK(self, file):
    (exstat, out) = commands.getstatusoutput("h5dump -A " + file)

    if exstat == 0:
      return True
    else:
      return False

  def _clumpsOK(self):
    cfile = self.dir + "clumps.h5part"
    if path.exists(cfile):
      return self._hdf5OK(self.dir + "clumps.h5part")
    else:
      return False
  
  def _redoClumpsLocally(self):
    print self.params.key
    sdir = self.dir
    
    orgclmpfile = sdir + "clumps.h5part"
    
    if not path.exists(orgclmpfile):
      cmd = "mv " + sdir + "clumps.h5part " + sdir + "clumps_sim.h5part"
      (exstat, out) = commands.getstatusoutput(cmd)
    else:
      cmd = "rm " + sdir + "clumps.h5part "
      (exstat, out) = commands.getstatusoutput(cmd)

    mmassclmp = 0.01*self.mtot
    mmassorb  = 0.01*self.mtot

    #for (dfile,dtime,dmtime) in self.dumps:
    if len(self.dumps) > 0:
      (dfile,dtime,dmtime) = self.dumps[-1]
      cmd = "find_clumps_A_ " + sdir + dfile + " " + sdir + \
            "clumps.h5part " + str(mmassclmp) + " " + str(mmassorb)
      (stat, out) = commands.getstatusoutput(cmd)
      print out

    
  def _redoClumps(self):
    if not len(self.cfg.subcmdsgl) > 0:
      self._redoClumpsLocally()
      return

    sdir = self.dir
    cmd = "rm " + sdir + "clumps.h5part"
    (exstat, out) = commands.getstatusoutput(cmd)
    
    jobscr  = self.dir + "redoclumps.sh"
    scrfile = open(jobscr,"w")
    
    mmassclmp = 1.e-4*self.mtot
    mmassorb  = 1.e-3*self.mtot

    for (dfile,dtime,dmtime) in self.dumps:
      #cmd = "find_clumps_A_ " + sdir + dfile + " " + sdir + \
      #      "clumps.h5part -1. -1."
      cmd = "find_clumps_A_ " + sdir + dfile + " " + sdir + \
            "clumps.h5part " + str(mmassclmp) + " " + str(mmassorb)
      print >>scrfile, cmd
      #(stat, out) = commands.getstatusoutput(cmd)
      #print out
    print >>scrfile, "rm " + jobscr
    scrfile.close()
    os.chmod(jobscr, stat.S_IRWXU)
    jobsub = self.cfg.subcmdsgl
    jobsub = jobsub.replace("$JOBNAME", "redoClp_"+self.params.key)
    jobsub = jobsub.replace("$JOBCMD" , "redoclumps.sh")

    oldwd = os.getcwd()
    os.chdir(self.dir)
    (exstat, out) = commands.getstatusoutput(jobsub)
    os.chdir(oldwd)
    print out

  def _findNoClumps(self):
    sdir = self.dir
    rcmd = "h5part_readattr -k noclumps -i "
    wcmd = "h5part_writeattr -k noclumps -v -2. -i "

    jobscr  = sdir + "doclumps.sh"
    (exstat, out) = commands.getstatusoutput("rm " + jobscr)
    scrfile = open(jobscr,"w")

    if not self.clumpjobid == 0:
      jobdel = self.cfg.delcmd
      jobdel = jobdel.replace("$JOBID", self.clumpjobid)
      (exstat, out) = commands.getstatusoutput(jobdel)
      print out

    for (dfile,dtime,dmtime) in self.dumps:
      (exstat, out) = commands.getstatusoutput(rcmd + sdir + dfile)
      if not exstat == 0:
        print sdir + dfile, " corrupted!"
	continue
      noc = int(float(out[10:]))
      if noc == -1:
        cmd = "find_clumps_A_ " + sdir + dfile + " " + sdir + \
            "clumps.h5part -1. -1."
        print >>scrfile, cmd
    print >>scrfile, "rm " + jobscr
    scrfile.close()
    
    os.chmod(jobscr, stat.S_IRWXU)
    jobsub = self.cfg.subcmdsgl
    jobsub = jobsub.replace("$JOBNAME", "doClp_"+self.params.key)
    jobsub = jobsub.replace("$JOBCMD" , "doclumps.sh")

    oldwd = os.getcwd()
    os.chdir(self.dir)
    (exstat, out) = commands.getstatusoutput(jobsub)
    os.chdir(oldwd)
    print out
    self.clumpjobid = re.split( '\W+', out )[2]

  def _redoClumpsIfNotOK(self):
    if not self._clumpsOK():
       self._redoClumps()

  def _delCorrupt(self):
    for (dfile, dtime, mtime) in self.dumps:
      if not self._hdf5OK(self.dir + dfile):
        (exstat, out) = commands.getstatusoutput("rm " + self.dir + dfile)
        print "rm " + self.dir + dfile
    if not self._clumpsOK():
      (exstat, out) = commands.getstatusoutput("rm " + self.dir + "clumps.h5part")
      print "rm " + self.dir + "clumps.h5part"

  def _redoEnergies(self):
    sdir = self.dir
    jobscr  = self.dir + "redoenergies.sh"
    scrfile = open(jobscr,"w")

    for (dfile,dtime,dmtime) in self.dumps:
      cmd = "get_energies__ " + sdir + dfile + " " + sdir + \
            "clumps.h5part "
      print >>scrfile, cmd
    print >>scrfile, "rm " + jobscr
    scrfile.close()
    os.chmod(jobscr, stat.S_IRWXU)
    jobsub = self.cfg.subcmdsgl
    jobsub = jobsub.replace("$JOBNAME", "redoEng_"+self.params.key)
    jobsub = jobsub.replace("$JOBCMD" , "redoenergies.sh")

    oldwd = os.getcwd()
    os.chdir(self.dir)
    (exstat, out) = commands.getstatusoutput(jobsub)
    print jobsub
    os.chdir(oldwd)
    print out

  def _getSecondaryImpact(self, rhomin):
    sdir = self.dir
    
    (dfile,dtime,dmtime) = self.dumps[-1]
    cmd = "find_fof " + sdir + dfile + " " + sdir + "last_fof.h5part 0. " + str(rhomin) + " 2.0 2 "
    (exstat, out) = commands.getstatusoutput(cmd)
    print out
    last_fof = SimClumps(self, "last_fof.h5part", loadsat=True, loadeng=True)
    (exstat, out) = commands.getstatusoutput("rm " + sdir + "last_fof.h5part")

    if last_fof.nofc[0] >= 2:
      r0vect = last_fof.pos[0,2,:] - last_fof.pos[0,1,:]
      v0vect = last_fof.vel[0,2,:] - last_fof.vel[0,1,:]
    
      mtar   = last_fof.m[0,1]
      mimp   = last_fof.m[0,2]
    
      Rtar   = last_fof.rclmp[0,1]
      Rimp   = last_fof.rclmp[0,2]

      G = self.cfg.G

      secgi = GiantImpact(mtar, mimp, Rtar, Rimp, G)
      secgi.getInitStateVecs(r0vect, v0vect)
      
      self.secgi = secgi
      self.last_fof = last_fof
      return True

    else:
      return False


class SimResults(object):
  def __init__(self):
    self.valid   = False
    self.problem = False
    self.blklst  = False
    self.noc = nan

    # mass
    self.mm = []
    self.mv = []
    
    # individual bodies, material
    self.mmatm = []
    self.mmatv = []
    
    # energies
    self.U0 = []
    self.Um = []
    self.Uv = []
    self.Ua = nan

    self.Epot0 = []
    self.Epotm = []
    self.Epotv = []

    self.Ekin0 = []
    self.Ekinm = []
    self.Ekinv = []
    
    self.Erot0 = []
    self.Erotm = []
    self.Erotv = []

    self.Etot0 = []
    self.dEtotm = []
    self.dEtotv = []
    
    self.Epotclmp = []
    self.Erotclmp = []
    self.Iclmp = []

    # dynamical parameters
    self.vinfm = []
    self.vinfv = []
    
    self.dblvarthetam = []
    self.dblvarthetav = []
    
    # mean densities
    self.rhom = []
    self.rhov = []

    # composition
    self.compm = []
    self.compv = []

    # rotational stuff
    self.Lm = []
    self.Lv = []
    
    self.Im = []
    self.Iv = []

    self.nopm = []
    self.nopv = []

    #self.trng = range(0., 0.)
    self.valtmin = nan
    self.valtmax = nan
    self.valnop  = 0

    # clump disk
    self.mdiskmatm = []
    self.mdiskmatv = []
    
    self.ediskmatm = []
    self.ediskmatv = []


  def resize(self,noc,nomat):
    self.noc = noc

    # mass
    self.mm = np.zeros([noc])
    self.mv = np.zeros([noc])
    
    # individual bodies, material
    self.mmatm = np.zeros([noc,nomat])
    self.mmatv = np.zeros([noc,nomat])
    
    # energies
    self.U0 = np.zeros([noc])
    self.Um = np.zeros([noc]) 
    self.Uv = np.zeros([noc]) 

    self.Epot0 = np.zeros([noc])
    self.Epotm = np.zeros([noc]) 
    self.Epotv = np.zeros([noc]) 

    self.EpotTm = nan
    self.EpotTv = nan

    self.Ekin0 = np.zeros([noc]) 
    self.Ekinm = np.zeros([noc]) 
    self.Ekinv = np.zeros([noc]) 
    
    self.Erot0 = np.zeros([noc]) 
    self.Erotm = np.zeros([noc]) 
    self.Erotv = np.zeros([noc]) 
    
    self.Epotclmp = np.zeros([noc])
    self.Erotclmp = np.zeros([noc])
    self.Iclmp = np.zeros([noc])

    # dynamical parameters
    self.vinfm = np.zeros([noc,3]) 
    self.vinfv = np.zeros([noc,3]) 
    
    self.dblvarthetam = np.zeros([noc]) 
    self.dblvarthetav = np.zeros([noc]) 
    
    # mean densities
    self.rhom = np.zeros([noc]) 
    self.rhov = np.zeros([noc]) 

    # composition
    self.compm = np.zeros([noc,nomat]) 
    self.compv = np.zeros([noc,nomat]) 

    # rotational stuff
    self.Lm = np.zeros([noc,3]) 
    self.Lv = np.zeros([noc,3])
    
    self.Im = np.zeros([noc]) 
    self.Iv = np.zeros([noc]) 

    self.nopm = np.zeros([noc]) 
    self.nopv = np.zeros([noc]) 
    
    # clump disk
    self.mdiskmatm = np.zeros([noc,nomat])
    self.mdiskmatv = np.zeros([noc,nomat])
    
    self.ediskmatm = np.zeros([noc,nomat])
    self.ediskmatv = np.zeros([noc,nomat])

  def blacklist(self):
    self.blklst = True
    self.valid  = False


class SimParam(object):
  def __init__(self, *args):
    self.valid = False
    argn = len(args)
    if ( argn == 4 ):
      self.mimp    = float(args[0])
      self.mtar    = float(args[1])
      self.impa    = float(args[2])
      self.vimprel = float(args[3])
      
      self.key = self._getKey()
      self.valid = True

    if ( argn == 1 ):
      key = args[0]
      
      valid = True
      (valid, self.mtar)    = self._attrmatch('mtar', valid, key)
      (valid, self.mimp)    = self._attrmatch('mimp', valid, key)
      (valid, self.impa)    = self._attrmatch('impa', valid, key)
      (valid, self.vimprel) = self._attrmatch('vimp', valid, key)

      #akey = self._getKey()
      #if akey == key:
      self.key   = key
      self.valid = valid

  def _attrmatch(self,attrkey,val,str):
    matcharr = re.findall(attrkey + '[0-9.e-]+', str)
    if ( not len(matcharr) == 1 ):
      return (False, float('nan') )
    else:
      return (val, float( matcharr[0].replace(attrkey, '') ) )

  def _getKey(self):
    retstr = ""
    if self.mtar < 0.001:
      retstr += ( 'mtar%07.1e' % self.mtar + '_' )
    else:
      retstr += ( 'mtar%07.3f' % self.mtar + '_' )
    
    if self.mimp < 0.001:
      retstr += ( 'mimp%07.1e' % self.mimp + '_' )
    else:
      retstr += ( 'mimp%07.3f' % self.mimp + '_' )

    retstr += 'impa%04.1f' % self.impa + '_' + \
        'vimp%04.2f' % self.vimprel
    #print retstr
    
    return retstr

    #return 'mtar%07.3f' % self.mtar + '_' \
    #    + 'mimp%07.3f' % self.mimp + '_' \
    #    + 'impa%04.1f' % self.impa + '_' \
    #    + 'vimp%04.2f' % self.vimprel

  def isSimilar(self, par, tol):
    fitmimp = not abs( (self.mimp - par.mimp ) / self.mimp ) > tol
    fitmtar = not abs( (self.mtar - par.mtar ) / self.mtar ) > tol
    fitimpa = not abs( (self.impa - par.impa ) / (self.impa+0.01) ) > tol
    fitvimp = not abs( (self.vimprel - par.vimprel ) / self.vimprel ) > tol
    return fitmimp and fitmtar and fitimpa and fitvimp



class SimSet(object):
  sims = {}
  def __init__(self, cfg):
    if cfg.dir == "":
      cfg.dir = sim_machine.basedir + cfg.name
    else:
      cfg.dir += "/" + cfg.name

    cfg.resultsdb = cfg.dir + "/results"
    cfg.bodiesdb  = cfg.dir + "/bodies"

    cfg.subcmd = sim_machine.qsub_omp
    cfg.delcmd = sim_machine.qdel
    cfg.subcmdsgl = sim_machine.qsub_sgl
    
    self.cfg = cfg
    
    bdir = self.cfg.dir
    if not path.exists(bdir):
      os.mkdir(bdir)
    
    for file in os.listdir(bdir):
      cparams = SimParam(file)
      if cparams.valid:
        self.sims[cparams.key] = Simulation( cparams, self.cfg )

  def discoverSims(self):
    self.sims = {}
    for file in os.listdir(self.cfg.dir):
      cparams = SimParam(file)
      if cparams.valid:
        self.sims[cparams.key] = Simulation( cparams, self.cfg )


class SimSetConfig(object):
  def __init__(self):
    self.dir = ""
    self.name = ""
    self.bodiesdb = ""
    self.resultdb = ""

    self.G = G
    self.T = 3.50e11

    self.bodtol = 0.05
    self.smotol = 0.50 
    self.desnop = 200000
    self.relsep = 5.0

    self.attr = []
    self.auxfiles = "./aux/ANEOS.QUARTZ ./aux/aneos_tables.hdf5"
    self.srcdir   = "./repos/sphlatch/apps/simple_sph"
    self.maketarg = "simple_sph_GCSHUmD_"
    self.binary   = "simple_sph_GCSHUmD_"

    self.subcmd  = ""
    self.delcmd  = ""
    self.nocpus  = ""
    self.runarg  = ""

    self.tcolmin  =   5.0
    self.tcolmax  =  50.
    self.tcoldump =   0.5
    self.tcolanim =   0.10
    
    self.stopnoclumps     =   2
    self.stopdnopdtcol    =  50
    self.stopdnopdtcolesc =  50

    self.idfilt = 2000000

    self.cpconf = ClumpsPlotConfig()
    self.mpconf = MoonsPlotConfig()

    self.keepdumps = [0., 1., 10., 50.]
    self.keepfiles = "clumps.h5part disk.hdf5 moons.pdf"


