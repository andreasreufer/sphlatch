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

class SGEaccounter(object):


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


