#!/usr/bin/env python

import numpy as np
import matplotlib as mp

import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.mpl    as mpl


def latexExp10(str):
  if 'e+' in str:
    str = str.replace('e+0', 'e+')
    str = str.replace('e+', '*10^{')
    str += '}'
  if 'e-' in str:
    str = str.replace('e-0', 'e-')
    str = str.replace('e-', '*10^{')
    str += '}'
  return str

def expOnly10(xstr):
  return str(int(xstr.split("e")[-1]))

def retnullstring(x, pos):
    return ''

def retexpstring(x, pos):
  return '$' + latexExp10( '%4.1e' % (x/100.) ) + '$'

def retexponlystring(x, pos):
  return '$10^{' + expOnly10('%e' % x) + '}$' 

def retexpmstring(x, pos):
  return '$' + latexExp10( '%4.1e' % (x/100.) ) + 'm$'

def retmathstring(x, pos):
  return '$' + ( '%4.1f' % x ) + '$'

def rettexstring(x, pos):
  return '$' + str(x) + '$'
  
def retintstring(x, pos):
  return '$' + str(int(x)) + '$'


empty_formatter = plt.FuncFormatter(retnullstring)
exp_formatter = plt.FuncFormatter(retexpstring)
exponly_formatter = plt.FuncFormatter(retexponlystring)
expm_formatter = plt.FuncFormatter(retexpmstring)
math_formatter = plt.FuncFormatter(retmathstring)
tex_formatter = plt.FuncFormatter(rettexstring)
int_formatter = plt.FuncFormatter(retintstring)

def cnameToRGB(str):
  return np.array( col.hex2color( col.cnames[str] ) )

def cnameToRGBA(str, alpha=1.0):
  rgba = np.zeros(4)
  rgba[0:3] = np.array( col.hex2color( col.cnames[str] ) )
  rgba[3] = alpha
  return rgba

class ScatterToArray(object):
  def __init__(self):
    self.pts = {}

  def clear(self):
    self.pts = {}

  def newPoint(self, key, x, val):
    pts = self.pts
    if not pts.has_key(key):
      pts[key] = []
    pts[key].append( (x, val) )

  def getArray(self):
    pts = self.pts
    kl = pts.keys()
    kl.sort()

    xls = []
    valls = []

    for k in kl:
      cpts = pts[k]
      cpts.sort(key=lambda tpl: tpl[0])

      cxls = []
      cvalls = []
      for p in cpts:
        cxls.append(p[0])
        cvalls.append(p[1])

      xls.append(cxls)
      valls.append(cvalls)

    return (kl, xls, valls)

def getJetColor(x,xmin,xmax):
  xrel = (x - xmin) / (xmax-xmin)

  if xrel < 0.:
    xrel = 0.
  
  if xrel > 1.:
    xrel = 1.

  return ( mp.cm.jet(xrel) )

def getVimpColor(vimp):
  jmap = mp.cm.jet

  #vimps = [ (1.00,"red","x--"),
  #          (1.05,"orange","x--"),
  #          (1.10,"green","x--"),
  #          (1.15,"blue","x--"),
  #          (1.20,"cyan","x--"),
  #          (1.30,"brown","x--"),
  #          (1.40,"darkred","+--"),
  #          (1.50,"darkorange","+--"),
  #          (1.60,"darkgreen","+--"),
  #          (1.70,"darkblue","+--"),
  #          (1.80,"darkcyan","+--"),
  #          (2.00,"darkkhaki","+--"),
  #          (2.50,"pink","x--"),
  #          (3.00,"magenta","x--"),
  #          (3.50,"lightgreen","x--"),
  #          (4.00,"lightblue","x--")]
  
  vimps = [ (1.00,"red","x-"),
            (1.05,"orange","x-"),
            (1.10,"green","x-"),
            (1.15,"blue","x-"),
            (1.20,"cyan","x-"),
            (1.30,"brown","x-"),
            (1.40,"darkred","+-"),
            (1.50,"darkorange","+-"),
            (1.60,"darkgreen","+-"),
            (1.70,"darkblue","+-"),
            (1.80,"darkcyan","+-"),
            (2.00,"darkkhaki","+-"),
            (2.50,"pink","x-"),
            (3.00,"magenta","x-"),
            (3.50,"lightgreen","x-"),
            (4.00,"lightblue","x-")]


  for (vimpi,coli,ls) in vimps:
    if abs( vimpi - vimp ) < 1.e-2:
      return (coli, ls)
      #return (jmap(coli), ls)

  return ("grey","p--")


def getAngleColor(angl):
  angls = [ (0.10,"red","x-"),
            (15.0,"orange","x-"),
            (22.5,"green","x-"),
            (30.0,"blue","x-"),
            (37.5,"cyan","x-"),
            (45.0,"magenta","x-"),
            (52.5,"pink","+-"),
            (60.0,"darkorange","+-"),
            (75.0,"darkgreen","+-"),
            (82.5,"darkblue","+-"),
            (89.0,"darkcyan","+-"),
            (89.5,"darkcyan","+-")]

  for (angli,coli,ls) in angls:
    if abs( angli - angl ) < 1.e-2:
      return (coli, ls)

  return ("grey","p--")


