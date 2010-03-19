#!/usr/bin/python

import pickle
import sys
import os
import time
import commands

class Logger(object):
  logfile = []

  def __init__(self, filename):
    self.logfile = open(filename, 'a')

  def __del__(self):
    self.logfile.close()

  def write(self, str):
    self.logfile.write(time.ctime() + " " + str + "\n")


log = Logger("gi_admin.log")

if os.path.exists("simset_config.sh"):
  execfile("simset_config.sh")
  log.write("simset_config.sh loaded")
else:
  log.write("simset_config.sh not found!")
  sys.exit(1)


(unprepared, prepared, queued, run, failed, finished, error, unknown) = range(8)


#def genSimTable():
#class Simulation(object):
  #status
  #def 

class SimParams(object):
  def __init__(self, mimp, mtar, b, vimp):
    self.mimp = mimp
    self.mtar = mtar
    self.b    = b
    self.vimp = vimp

class Body(object):
  def __init__(self,file, m, rc, 


class Simulation(object):
  state = idle
  simdir = ""

  def __init__(self,params):
    self.params = params

  def prepare(self):
    pass

  def start(self):
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
