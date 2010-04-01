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

(unprepared, prepared, queued, run, failed, finished, error, unknown) = range(8)
rad2deg = 360./(2.*pi)
deg2rad = 1./rad2deg


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

    self.key = 'mtar%07.3f' % mtar + '_' \
      + 'mimp%07.3f' % mimp + '_' \
      + 'impa%04.1f' % impa + '_' \
      + 'vimp%04.1f' % vimp


class Body(object):
  def __init__(self, file, m, r, T, res):
    self.file = file
    self.m = m
    self.r = r
    self.T = T
    self.res = res



    
class Simulation(object):
  #state = unprepared
  #simdir = ""
  #jobid = -1

  def __init__(self,params):
    self.params = params
    ssbdir = params.cfg["SIMSETBDIR"]
    self.dir = path.normpath( ssbdir + "/" + params.key + "/" )

    # fetch bodies

  def prepare(self):
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
    #  combine bodies
    #  write attributes
    #  
    # compile code
    #pass

  def run(self):
    # check if everything's there
    pass

  def archive(self):
    pass

  def getKey(self):
    return params.key


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



def getJobStats():
  alljobs = {}

  (stat, runwaitraw) = commands.getstatusoutput("qstat -g d")
  if stat != 0:
    return
  for line in runwaitraw.splitlines()[2:]:
    lsplt = line.split()
    if len(lsplt) == 9:
      (id, prio, name, user, statestr, date, time, queue, slots) = line.split()
      alljobs[int(id)] = (name, user, run)
    if len(lsplt) == 8:
      (id, prio, name, user, statestr, date, time, slots) = line.split()
      state = unknown
      if statestr == "qw":
        state = queued
      alljobs[int(id)] = (name, user, state)
  
  (stat, finishedraw) = commands.getstatusoutput("qstat -g d -s z")
  if stat != 0:
    return unknown
  for line in finishedraw.splitlines()[2:]:
    print line.split()
    (id, prio, name, user, statestr, date, time, slots) = line.split()
    alljobs[int(id)] = (name, user, finished)

  alljobs['timestamp'] = time.time()

  return alljobs


log = Logger("gi_admin.log")
simset = SimSet(log)



