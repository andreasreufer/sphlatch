#!/usr/bin/python

import pickle
import sys
import os
import os.path as path
import time
import commands
import shelve
import numpy as np
import shutil

from numpy import sqrt, sin, cos, arccos, log, abs, tan, arctan2, pi, zeros
from body import BodyFile
    
import pdb

(unprepared, prepared, queued, run, failed, finished, error) = range(7)
rad2deg = 360./(2.*pi)
deg2rad = 1./rad2deg

# constants to transform the initial parameters to physical values
G   = 6.67429e-8
Re  = 6.37814e8
Me  = 5.97360e27
lem = 3.40000e41 # Canup (2001) value
#lem = 2.88993e41 # actual (?) value



class Logger(object):
  logfile = []

  def __init__(self, filename):
    self.logfile = open(filename, 'a')

  def __del__(self):
    self.logfile.close()

  def write(self, str):
    self.logfile.write(time.ctime() + " " + str + "\n")
    self.logfile.flush()



class SimParams(object):
  def __init__(self, mimp, mtar, impa, vimp, cfg):
    self.mimp = mimp
    self.mtar = mtar
    self.impa = impa
    self.vimp = vimp

    self.cfg = cfg

    self.temp = float(cfg["TEMP"])

    self.key = 'mtar%07.3f' % mtar + '_' \
      + 'mimp%07.3f' % mimp + '_' \
      + 'impa%04.1f' % impa + '_' \
      + 'vimp%04.1f' % vimp

class Task(object):
  def __init__(self, callfunc):
    self.callfunc = callfunc

  def __run__(self):
    self.callfunc()


class Simulation(object):
  def __init__(self,params):
    self.params = params
    ssbdir = params.cfg["SIMSETBDIR"]
    self.dir = path.normpath(ssbdir + "/sim_" + params.key + "/")
    self.state = unprepared
    # fetch bodies


   #(unprepared, prepared, queued, run, failed, finished, error) = range(8)
  def getState(self):
    if not path.exists(self.dir):
      self.state = unprepared
      return self.state
    else:

      neededfiles = self.params.cfg["AUXFILES"].split()
      neededfiles.append(self.params.cfg["BINARY"])
      neededfiles.append("initial.h5part")

      for file in neededfiles:
        if not os.path.exists(self.dir + os.path.split(file)[1]):
          return self.state

      # so everything's ready
      if jobid == 0:
        self.state = prepared
        return self.state
      else:






  def openTasks(self):
    return Task(self.prepare)

  def _prepare(self):
    # if it does not yet exist, make dir
    if not path.exists(self.dir):
      os.mkdir(self.dir)

    # copy files
    ssbdir = self.params.cfg["SIMSETBDIR"]
    auxf = self.params.cfg["AUXFILES"].split()
    
    self.log = Logger(self.dir + "/logfile.txt")
    
    for file in auxf:
      shutil.copy2(ssbdir + file, self.dir)
    
    self.log.write("auxiliary files copied")
    
    #if not path.exists(self.dir + "initial.h5part"):
    
    # prepare bodies:
    #  find fitting bodies
    (self.tarb, self.impb) = self.findBodies()
    #  combine bodies

    #  write attributes
    #  
    # compile code
    #pass

  def _run(self):
    # check if everything's there
    pass

  def _archive(self):
    pass

  def _getKey(self):
    return params.key

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
    
    # find fitting pairs
    pairs = []
    for tarbd in targcand:
      for impbd in impacand:
        print 
        if abs((impbd.h - tarbd.h)/tarbd.h) < toler:
          pairs.append( (tarbd, impbd) )

    # sort according to smoothing length
    pairs.sort(key=lambda bod: bod[0].h)
    
    boddb.close()
    return pairs[0]


class SimSet(object):
  sims = {}
  #simsetcfg = {}
  #log = []

  def __init__(self, logger):
    self.log = logger
    self.simsetcfg = {}
    
    if path.exists("simset_config.sh"):
      execfile("simset_config.sh", self.simsetcfg)
      log.write("simset_config.sh loaded")
    else:
      log.write("simset_config.sh not found!")
      log.write("exiting")
      sys.exit(1)
    
    mtara = np.array( self.simsetcfg["MTAR"].split(), dtype=float)
    mimpa = np.array( self.simsetcfg["MIMP"].split(), dtype=float)
    vimpa = np.array( self.simsetcfg["VIMP"].split(), dtype=float)
    impaa = np.array( self.simsetcfg["IMPA"].split(), dtype=float)

    bdir = self.simsetcfg["SIMSETBDIR"]
    cond = self.simsetcfg["SIMCOND"]
    if len(cond) < 1:
      cond = "True"

    for mtar in mtara:
      for mimp in mimpa:
        for vimp in vimpa:
          for impa in impaa:
            if ( eval(cond) ):
              simpar = SimParams(mimp, mtar, impa, vimp, self.simsetcfg)
              if not self.sims.has_key(simpar.key):
                self.sims[simpar.key] = Simulation( simpar )


  def __del__(self):
    pass


log = Logger("gi_admin.log")
simset = SimSet(log)
sim = simset.sims['mtar000.100_mimp000.100_impa75.0_vimp03.0']
pairs = sim._findBodies()

