#!/usr/bin/env python
import tables as pt

class H5PartDump(object):
  def __init__(self,fname):
    self.fname = fname
    self.pth = pt.openFile(fname, "r")

  def __del__(self):
    self.pth.close()

  def forEachStep(self,func):
    for grp in self.pth.walkGroups():
      if "Step" in grp._v_pathname:
        func(grp)

  def getStep(self,sname="/current"):
    pth = self.pth

    for node in pth.walkNodes():
      if node._v_pathname == sname:
        if type(node) == pt.group.Group:
          return node
        else:
          return self.getStep(node.target)
    return []

  def getStepNames(self):
    stpnam = []
    for node in self.pth.walkGroups():
      if "/Step" in node._v_pathname:
        stpnam.append(node._v_pathname)
    return stpnam


