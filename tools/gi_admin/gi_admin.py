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
  print "haba"
else:
  log.write("simset_config.sh not found!")
  sys.exit(1)


#def genSimTable():




#class Simulation(object):
#  status
#
#  def 

(idle, queued, run, failed, finished, error, unknown) = range(5)
def getJobStatus(_id):
  (stat, out) = commands.getstatusoutput("qacct -j " + str(id) )
  if stat == 32512:
    return unknown
  if stat == 0:
    return finished

  (stat, out) = commands.getstatusoutput("qstat -g d")
  if stat != 0:
    return unknown

  lines = out.splitlines()
  for line in lines:
    (id, prio, name, user, state, date, time, queue, slots) = line.split()
    if id == _id:
      if state == 'r':
        return run
      else if state == 'qw':
        return queued
      else if state.find('E') >= 0:
        return error
  return unknown


