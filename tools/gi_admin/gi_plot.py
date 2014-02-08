#!/usr/bin/env python
import sys
import numpy as np
import commands

import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.mpl    as mpl

from fewbody import FewBodies
from h5part import H5PartDump
from plot_helpers import *
from const_cgs import *

import tables as pt

mp.rc('text.latex', preamble = '\usepackage{amssymb}, \usepackage{wasysym}')

def plotDummy(norma, physa, plt, pdump, cdump, cfg):
  pass

def plotDiskMass(norma, physa, plt, pdump, cdump, cfg):
  ML = cfg.ML
  
  # orbit: not set = 0, clump = 1, disk = 2, re-impact = 3, escape = 4
  disk = ( pdump.orbit[:,0] == 2 ) & ( pdump.clumpid[:,0] == 1)
  
  sio2 = ( pdump.mat[:,0] == 1 )
  ice  = ( pdump.mat[:,0] == 2 )
  iron = ( pdump.mat[:,0] == 5 )
  
  targ = ( pdump.id[:,0]  < cfg.idfilt )
  impa = ( pdump.id[:,0] >= cfg.idfilt )

  m = pdump.m[:,0]

  mdiskrocktarg = sum( m.compress( disk & sio2 & targ) )
  mdiskrockimpa = sum( m.compress( disk & sio2 & impa) )
  mdiskice      = sum( m.compress( disk & ice ) )
  mdiskiron     = sum( m.compress( disk & iron ) )

  mdiskrock   = mdiskrocktarg + mdiskrockimpa
  trgfracrock = 0.
  if mdiskrock > 0.:
     trgfracrock = mdiskrocktarg / mdiskrock
  
  txt = r"$\mathrm{Disk}:$" + "\n" + \
      r"$ m_{\mathrm{rock}} = " + ("%6.3f m_L" % ( mdiskrock / ML) ) + (" (%3.0f" % (trgfracrock*100.)) + r"\%_{\mathrm{target}} )$" + "\n" + \
      r"$ m_{\mathrm{ice}} = " + ("%6.3f m_L" % ( mdiskice / ML )) + "$\n" +\
      r"$ m_{\mathrm{iron}} = " + ("%6.3f m_L" % ( mdiskiron / ML )) + "$"
  norma.text( cfg.dmass_vc[0], cfg.dmass_vc[1], txt, \
      color=cfg.txtc, size=cfg.time_txts)

def plotGrid(norma, physa, plt, pdump, cdump, cfg):
  physa.grid(True, lw=0.1, color='darkgrey', alpha=0.5)
    
def plotSimParams(norma, physa, plt, pdump, cdump, cfg):
  paramstxt = '$m_T =' + ("%3.2f" % cfg.mtar) + r' M_E, ' +\
      'm_I = ' + ("%4.3f" % cfg.mimp) + 'M_E, ' +\
      r'v_{imp} = ' + ("%4.3f" % cfg.vimp) + r'v_{esc}, ' +\
      ("%2.0f" % cfg.impa) + r'^\circ$'
      #("%2.0f" % cfg.impa) + r'^\circ, '+\
      #'T = ' + ("%4.0f" % cfg.T ) + 'K$'

  norma.text( cfg.parm_vc[0], cfg.parm_vc[1], paramstxt, \
      color=cfg.txtc, size=cfg.parm_txts )

def plotSimParams(norma, physa, plt, pdump, cdump, cfg):
  paramstxt = '$M_{tar} =' + ("%3.2f" % cfg.mtar) + r' M_{\oplus}, ' +\
      'M_{imp} = ' + ("%4.3f" % cfg.mimp) + 'M_{\oplus}, ' +\
      r'v_{imp} = ' + ("%4.3f" % cfg.vimp) + r'v_{esc}, ' +\
      ("%2.0f" % cfg.impa) + r'^\circ$'
      #("%2.0f" % cfg.impa) + r'^\circ, '+\
      #'T = ' + ("%4.0f" % cfg.T ) + 'K$'

  norma.text( cfg.parm_vc[0], cfg.parm_vc[1], paramstxt, \
      color=cfg.txtc, size=cfg.parm_txts )

def colorBodyAndPhaseANEOS(pdump, cdump, filt, cfg):
  cmap = np.zeros((2,6,3))
  cmap[0,1] = cnameToRGB("darkgrey") # body 1 liqid/solid
  cmap[0,2] = cnameToRGB("red")      # body 1 vapour phase
  cmap[1,1] = cnameToRGB("grey")     # body 2 liqid/solid
  cmap[1,2] = cnameToRGB("orange")   # body 2 vapour phase

  phase = pdump.phase[:,0].compress(filt)[:]
  bod = np.int32( (pdump.id[:].compress(filt)[:] > cfg.idfilt) )

  color = cmap[bod, phase]
  pt2size = cfg.pt2size

  return (color, pt2size)




def colorBodyAndMat(pdump, cdump, filt, cfg):
  cmap = np.zeros((2,6,3))

  #cmap[0,0] = cnameToRGB("lawngreen")  # target SiO2
  cmap[0,0] = cnameToRGB("plum")       # target gas
  cmap[0,1] = cnameToRGB("red")        # target SiO2
  cmap[0,2] = cnameToRGB("darkcyan")   # target ice
  cmap[0,4] = cnameToRGB("red")        # target dunite
  cmap[0,5] = cnameToRGB("mediumblue") # target iron

  cmap[1,0] = cnameToRGB("lawngreen")  # impactor gas
  cmap[1,1] = cnameToRGB("orange")     # impactor SiO2
  cmap[1,2] = cnameToRGB("cyan")       # impactor ice
  cmap[1,4] = cnameToRGB("orange")     # impactor dunite
  cmap[1,5] = cnameToRGB("dodgerblue") # impactor iron

  mat = pdump.mat[:,0].compress(filt)[:]
  bod = np.int32( (pdump.id[:].compress(filt)[:] > cfg.idfilt) )

  color = cmap[bod, mat]
  pt2size = cfg.pt2size

  return (color, pt2size)

def colorBodyAndMatPaper(pdump, cdump, filt, cfg):
  cmap = np.zeros((2,6,3))

  #cmap[0,0] = cnameToRGB("lawngreen")  # target SiO2
  cmap[0,0] = cnameToRGB("plum")       # target gas
  cmap[0,1] = cnameToRGB("darkred")        # target SiO2
  cmap[0,2] = cnameToRGB("darkcyan")   # target ice
  cmap[0,4] = cnameToRGB("red")        # target dunite
  cmap[0,5] = cnameToRGB("mediumblue") # target iron

  cmap[1,0] = cnameToRGB("lawngreen")  # impactor gas
  cmap[1,1] = cnameToRGB("darkorange")     # impactor SiO2
  cmap[1,2] = cnameToRGB("cyan")       # impactor ice
  cmap[1,4] = cnameToRGB("orange")     # impactor dunite
  cmap[1,5] = cnameToRGB("dodgerblue") # impactor iron

  mat = pdump.mat[:,0].compress(filt)[:]
  bod = np.int32( (pdump.id[:].compress(filt)[:] > cfg.idfilt) )

  color = cmap[bod, mat]
  pt2size = cfg.pt2size

  return (color, pt2size)

def colorBodyAndOrbAndMat(pdump, cdump, filt, cfg):
  cmap = np.zeros((2,2,6,3))

  cmap[1,0,0] = cnameToRGB("lawngreen")  # target SiO2
  cmap[1,0,1] = cnameToRGB("red")        # target SiO2
  cmap[1,0,2] = cnameToRGB("darkcyan")   # target ice
  cmap[1,0,4] = cnameToRGB("red")        # target dunite
  cmap[1,0,5] = cnameToRGB("mediumblue") # target iron

  cmap[1,1,0] = cnameToRGB("lawngreen")  # target SiO2
  cmap[1,1,1] = cnameToRGB("orange")     # impactor SiO2
  cmap[1,1,2] = cnameToRGB("cyan")       # impactor ice
  cmap[1,1,4] = cnameToRGB("orange")     # impactor dunite
  cmap[1,1,5] = cnameToRGB("dodgerblue") # impactor iron

  cmap[0,:,:] = 0.5*cmap[1,:,:]

  mat = pdump.mat[:,0].compress(filt)[:]
  bod = np.int32( (pdump.id[:].compress(filt)[:] > cfg.idfilt) )
  orb = np.int32( (pdump.orbit[:].compress(filt)[:] == 2) & (pdump.clumpid[:].compress(filt)[:] == 1))
  
  color = cmap[orb, bod, mat]
  pt2size = cfg.pt2size

  return (color, pt2size)

def colorOrbit(pdump, cdump, filt, cfg):
  cmap = np.zeros((6,3))
  cmap[0] = cnameToRGB("darkgrey")
  cmap[1] = cnameToRGB("blue")
  cmap[2] = cnameToRGB("green")
  cmap[3] = cnameToRGB("orange")
  cmap[4] = cnameToRGB("red")

  orbit = pdump.orbit[:,0].compress(filt)[:]

  color = cmap[orbit]
  pt2size = cfg.pt2size

  return (color, pt2size)


def filtLaterOrb(pdump, cdump, filt, cfg):
  ldumpf = H5PartDump(cfg.lfile)
  sname  = (ldumpf.getStepNames())[0]
  ldump = ldumpf.getStep(sname)
  
  filt = (ldump.orbit[:,0] == cfg.orbtype ) & (ldump.clumpid[:,0] == 1 )
  return filt
  

def colorTemp(pdump, cdump, filt, cfg):
  T = (pdump.T[:,0].compress(filt)[:])*eVinK
  Tmin = cfg.Tmin
  Tmax = cfg.Tmax

  color = (T - Tmin) / (Tmax - Tmin)
  pt2size = cfg.pt2size
  return (color, pt2size)

def colorDens(pdump, cdump, filt, cfg):
  rho = (pdump.rho[:,0].compress(filt)[:])
  rhomin = cfg.rhomin
  rhomax = cfg.rhomax

  color = (rho - rhomin) / (rhomax - rhomin)
  pt2size = cfg.pt2size
  return (color, pt2size)

def colorPress(pdump, cdump, filt, cfg):
  p = (pdump.p[:,0].compress(filt)[:])
  pmin = cfg.pmin
  pmax = cfg.pmax

  color = (p - pmin) / (pmax - pmin)
  pt2size = cfg.pt2size
  return (color, pt2size)


def colorClumps(pdump, cdump, filt, cfg):
  cid = pdump.clumpid[:,0].compress(filt)[:]

  maxnoc = 10
  maxcid = maxnoc - 1
  cid[cid > maxcid] = maxcid
  
  cmap = np.zeros((maxnoc,4))

  cmap[0] = cnameToRGBA("lime") # escaping stuff
  cmap[1] = cnameToRGBA("red") #
  cmap[2] = cnameToRGBA("blue") # 
  cmap[3] = cnameToRGBA("yellow") # 
  cmap[4] = cnameToRGBA("orchid") # 
  cmap[5] = cnameToRGBA("aqua") # 
  cmap[6] = cnameToRGBA("orange") # 
  cmap[7] = cnameToRGBA("royalblue") # 
  cmap[8] = cnameToRGBA("fuchsia") # 
  cmap[9] = cnameToRGBA("darkgrey") # beyond maximum clump number

  color = cmap[cid]
  pt2size = cfg.pt2size
  return (color, pt2size)

def colorScalarLinear(pdump, cdump, filt, cfg):
  sname = cfg.scal_name
  for lv in pdump._v_leaves.values():
    if lv._v_name == sname:
      scal = lv[:,0].compress(filt)[:]

  scalmin = cfg.scal_min
  scalmax = cfg.scal_max
  
  if cfg.cbar_plot:
    cfg.cbar_ticksx = np.arange(0., 1.00001, 1./(cfg.cbar_noticks-1))

    cfg.cbar_tickstxt = []
    if cfg.cbar_noticks > 1:
      for x in  cfg.cbar_ticksx:
        scalsc = (x*(scalmax-scalmin)+scalmin)*cfg.cbar_sc
        str    = '$' + latexExp10( cfg.cbar_fmt % scalsc ) + cfg.cbar_ut + '$'
        cfg.cbar_tickstxt.append(str) 
  
  color = (scal - scalmin) / (scalmax - scalmin)
  pt2size = cfg.pt2size
  return (color, pt2size)

def colorScalarLog10(pdump, cdump, filt, cfg):
  sname = cfg.scal_name
  for lv in pdump._v_leaves.values():
    if lv._v_name == sname:
      scal = lv[:,0].compress(filt)[:]

  lscalmin = np.log10(cfg.scal_min)
  lscalmax = np.log10(cfg.scal_max)
  scalmin = cfg.scal_min
   
  if cfg.cbar_plot:
    cfg.cbar_ticksx = np.arange(0., 1.00001, 1./(cfg.cbar_noticks-1))
    cfg.cbar_tickstxt = []
    for x in  cfg.cbar_ticksx:
      scalsc = ( scalmin*pow(10., (lscalmax-lscalmin)*x) )*cfg.cbar_sc
      str    = '$' + latexExp10( cfg.cbar_fmt % scalsc ) + cfg.cbar_ut + '$'
      cfg.cbar_tickstxt.append(str) 
  
  color = (np.log10(scal) - lscalmin) / (lscalmax - lscalmin)
  pt2size = cfg.pt2size
  return (color, pt2size)



def interpolateScalar(norma, physa, plt, pdump, cdump, cfg):
  nx = cfg.dpi*cfg.xinch
  ny = cfg.dpi*cfg.yinch
  nx = 1280
  ny = 720

  r0 = [0.,0.,0.]
  rx = [0.,0.,0.]
  ry = [0.,0.,0.]

  dx = ( cfg.corrax[1] - cfg.corrax[0] ) / nx
  dy = ( cfg.corrax[3] - cfg.corrax[2] ) / ny

  rx[cfg.X] = dx
  ry[cfg.Y] = dy

  r0[cfg.X] = cfg.corrax[0] + dx*0.5
  r0[cfg.Y] = cfg.corrax[2] + dy*0.5
  

  r0str = str(r0).replace(" ","")
  rxstr = str(rx).replace(" ","")
  rystr = str(ry).replace(" ","")

  nxstr = str(int(nx))
  nystr = str(int(ny))

  
  tfile = cfg.pfile[:-7] + "_" + cfg.gridvarname + "_map.hdf5"

  cmd = cfg.gridcmd + " " + cfg.pfile + " " + tfile + " " + cfg.gridvarname + " " + nxstr + " " + nystr + " " + r0str + " " + rxstr + " " + rystr
  #cmd = "sph2grid_S " + cfg.pfile + " " + tfile + " " + cfg.gridvarname + " " + nxstr + " " + nystr + " " + r0str + " " + rxstr + " " + rystr
  (stat, out) = commands.getstatusoutput(cmd)
  if cfg.verbose:
    print out

  varmapf = pt.openFile(tfile)
  varmap = []
  if cfg.gridlog:
    varmap = np.log10( ( varmapf.root._g_loadChild(cfg.gridvarname) )[:,:] )
  else:
    varmap = ( varmapf.root._g_loadChild(cfg.gridvarname) )[:,:]
  varmapf.close()
  
  cfg.cf = physa.imshow(varmap.transpose(), alpha=1.0, origin='lower', aspect=9./16., cmap=cfg.gridcmap, vmin=cfg.gridvmin, vmax=cfg.gridvmax, extent=cfg.corrax)
  cfg.varmap = varmap
  
  (stat, out) = commands.getstatusoutput("rm " + tfile)
  if cfg.verbose:
    print out

def papercolors(cfg):
  cfg.axisbg = "white"
  cfg.clpc_ec = 'black'
  cfg.clpp_ec = 'black'
  cfg.clpr_ec = 'green'
  cfg.clpr_ls = "dotted"
  cfg.txtc    = 'black'


class GIplotConfig(object):
  def __init__(self):
    self.task = "partpos"

    # the plot dimensions
    self.xinch = 3.2
    self.yinch = 1.8
    self.dpi = 400

    self.ax = [-1.0e10, 1.0e10, -1.5e10, 1.5e10]
    
    self.axreluse = False
    self.axrel = [-32., 32., -18., 18.]
    self.axrelscl = 1.0

    # define the X/Y axes of the plot as dimension numbers
    self.X = 1
    self.Y = 0
    self.Z = 2

    self.filt = "clumpnegz"
    self.plotGas = True
    self.hmed = float('nan')

    self.ds = 5.e7
    self.pt2size = 4.00

    self.axisbg = "black"

    self.plotclmp = True
    self.plottraj = False
    self.clmptrajdt = 3600
    
    self.time = 0.
    self.G    = G

    self.ME = ME
    self.ML = ML
    self.Mmin    = -float('inf')
    self.MminPlt = 1.e-4*self.ME
    self.MminLbl = 5.e-3*self.ME
    self.sctsize = 0.02

    self.colorFunc = colorBodyAndMat
    self.cmap      = 'jet'
    self.txtc      = 'white'

    self.clpc_ec = 'white'
    self.clpc_fc = 'none'
    self.clpc_lw = 0.2
    self.clpc_al = 0.3

    self.clpp_ec = 'white'
    self.clpr_ec = 'green'
    self.clpp_fc = 'none'
    self.clpr_ls = "dotted"
    self.clpp_lw = 0.2
    self.clpp_al = 0.6
    self.clp_plotsurf = True
    self.clp_plotrche = True
    
    self.clpc_txts = 4

    self.scale_vc = [0.05, 0.06, 0.3, 0.06]
    self.scale_ut = 'm'
    self.scale_sc = 0.01
    self.scale_txts = 4
    
    self.time_vc = [0.05, 0.90]
    self.time_ut = 'h'
    self.time_sc = 0.000277777777
    self.time_txts = 4

    self.parm_vc = [0.5, 0.90]
    self.parm_txts = 4
    self.parm_txt = ""
    
    self.copy_vc = [0.65, 0.06]
    self.copy_txts = 3
    #self.copy_txt = "Andreas Reufer, University of Bern"
    self.copy_txt = "Andreas Reufer, Arizona State University"
  
    self.dmass_vc = [0.05, 0.2]

    self.cbar_plot        = False
    self.cbar_orientation = 'horizontal'
    self.cbar_extent      = [0.05, 0.85, 0.5, 0.01]
    self.cbar_noticks     = 4
    self.cbar_lw          = 0.1
    self.cbar_ut          = ''
    self.cbar_sc          = 1.
    self.cbar_txts        = 4
    self.cbar_fmt         = '%5.0f'

    self.idfilt = 2.e6

    self.rhosat = 2.65

    self.imgdir = "giplot"
    self.imgext = ".png"

    self.prePlot = [plotGrid]
    self.postPlot = [plotDummy]
    
    self.verbose = True 
    self.plotscale = True

    self.ORzoom  = False


cfg_ecc_m           = GIplotConfig()
cfg_ecc_m.parm_vc   = [0.6, 0.90]
cfg_ecc_m.dmass_vc  = [0.72, 0.20]
cfg_ecc_m.rmax      = 2.01e10
cfg_ecc_m.task      = "ecc"
cfg_ecc_m.colorFunc = colorBodyAndOrbAndMat
cfg_ecc_m.prePlot   = []
cfg_ecc_m.postPlot  = [plotSimParams, plotDiskMass]
cfg_ecc_m.plotscale = False


cfg_orb_XY_m           = GIplotConfig()
cfg_orb_XY_m.ax        = [-1.6e10, 8.0e9 , -9.0e9 , 4.5e9 ]
cfg_orb_XY_m.postPlot  = [plotSimParams, plotDiskMass]
cfg_orb_XY_m.colorFunc = colorBodyAndOrbAndMat
cfg_orb_XY_m.ORzoom    = False

cfg_eng_m           = GIplotConfig()
cfg_eng_m.parm_vc   = [0.6, 0.90]
#cfg_eng_m.dmass_vc  = [0.72, 0.20]
cfg_eng_m.rmax      = 5.01e10
cfg_eng_m.task      = "energy"
cfg_eng_m.colorFunc = colorBodyAndMat
cfg_eng_m.prePlot   = []
cfg_eng_m.postPlot  = [plotSimParams]
cfg_eng_m.plotscale = False


cfg_grd_rho_XY_m             = GIplotConfig()
cfg_grd_rho_XY_m.ax          = [-3.2e10, 1.6e10, -1.8e10, 0.9e10]
cfg_grd_rho_XY_m.prePlot     = [plotGrid, interpolateScalar]
cfg_grd_rho_XY_m.postPlot    = [plotSimParams]
cfg_grd_rho_XY_m.colorFunc   = colorBodyAndMat
cfg_grd_rho_XY_m.gridcmd     = "sph2grid_MSS"
cfg_grd_rho_XY_m.gridcmap    = "gray"
cfg_grd_rho_XY_m.gridvarname = "rho"
cfg_grd_rho_XY_m.gridvmin    = -6.
cfg_grd_rho_XY_m.gridvmax    = 1.
cfg_grd_rho_XY_m.gridlog     = True
cfg_grd_rho_XY_m.plotGas     = False


cfg_grd_p___XY_m             = GIplotConfig()
cfg_grd_p___XY_m.ax          = [-3.2e10, 1.6e10, -1.8e10, 0.9e10]
cfg_grd_p___XY_m.prePlot     = [plotGrid, interpolateScalar]
cfg_grd_p___XY_m.postPlot    = [plotSimParams]
cfg_grd_p___XY_m.colorFunc   = colorBodyAndMat
cfg_grd_p___XY_m.gridcmd     = "sph2grid_MSS"
cfg_grd_p___XY_m.gridcmap    = "gray"
cfg_grd_p___XY_m.gridvarname = "p"
cfg_grd_p___XY_m.gridvmin    =  5.
cfg_grd_p___XY_m.gridvmax    = 10.
cfg_grd_p___XY_m.gridlog     = True
cfg_grd_p___XY_m.plotGas     = False


cfg_T___XY_n = GIplotConfig()
cfg_T___XY_n.ax        = [-3.2e10, 3.2e10, -1.8e10, 1.8e10]
cfg_T___XY_n.scal_min  =  1500. / eVinK
cfg_T___XY_n.scal_max  = 12000. / eVinK
cfg_T___XY_n.scal_name = "T"
cfg_T___XY_n.cbar_plot = True
cfg_T___XY_n.cbar_sc   = eVinK
cfg_T___XY_n.cbar_ut   = "K"
cfg_T___XY_n.cbar_fmt  = '%5.0f'
cfg_T___XY_n.colorFunc = colorScalarLog10
cfg_T___XY_n.postPlot  = [plotSimParams]
cfg_T___XY_n.ORzoom    = False

cfg_mat_XY_n = GIplotConfig()
cfg_mat_XY_n.ax = [-3.2e10, 3.2e10, -1.8e10, 1.8e10]
cfg_mat_XY_n.postPlot  = [plotSimParams]
cfg_mat_XY_n.ORzoom    = False

cfg_fat_XY_Z = GIplotConfig()
cfg_fat_XY_Z.ax = [-2.00e9 , 2.00e9 , -1.50e9 , 0.75e9 ]
cfg_fat_XY_Z.postPlot  = [plotSimParams]
cfg_fat_XY_Z.colorFunc = colorBodyAndMat
cfg_fat_XY_Z.filt = "aux"
cfg_fat_XY_Z.filtFunc = filtLaterOrb
cfg_fat_XY_Z.orbtype = 2



class GIplot(object):
  def __init__(self, cfg=GIplotConfig()):
    self.fig = plt.figure( figsize=(cfg.xinch, cfg.yinch) )
    
    self.physa = plt.axes( [0.0, 0.0, 1.0, 1.0], axisbg=cfg.axisbg)
    self.norma = plt.axes( [0.0, 0.0, 1.0, 1.0], frameon=False)
    self.norma.set_xticks([])
    self.norma.set_yticks([])
    
    self.cfg = cfg
    
    ax = cfg.ax
    if cfg.axreluse:
      ax[0] = cfg.r * cfg.axrelscl * cfg.axrel[0]
      ax[1] = cfg.r * cfg.axrelscl * cfg.axrel[1]
      ax[2] = cfg.r * cfg.axrelscl * cfg.axrel[2]
      ax[3] = cfg.r * cfg.axrelscl * cfg.axrel[3]
      if cfg.verbose:
        print "use relative scaled axis: ",ax
    
    xcen = ( ax[0] + ax[1] ) / 2.
    xscl = ( ax[1] - ax[0] ) / cfg.xinch
    ycen = ( ax[2] + ax[3] ) / 2.
    yscl = ( ax[3] - ax[2] ) / cfg.yinch

    scl = max( xscl, yscl )
    corrax = [ xcen - 0.5*scl*cfg.xinch, xcen + 0.5*scl*cfg.xinch, \
        ycen - 0.5*scl*cfg.yinch, ycen + 0.5*scl*cfg.yinch ]

    self.scl = scl
    self.corrax = corrax
    
    self.cfg.corrax = corrax

  def doPlot(self,pfile,cfile,ifile):
    cfg = self.cfg
    cfg.pfile = pfile
    cfg.cfile = cfile
    cfg.ifile = ifile
    
    self._load(pfile,cfile)
    self._plotPre()
    
    task = cfg.task
    if task == "partpos":
      self._prepPlotParts()
      self._plotParts(self.physa, self.corrax)

      if cfg.ORzoom:
        mtot  = np.sum( self.cdump.m[:,0] )
        mLR   = self.cdump.m[1,0]
      
        posLR = self.cdump.pos[1,:]
        posOR = -posLR*mLR/(mtot - mLR)
        
        relX = 0.5 / cfg.xinch
        relY = 0.5 / cfg.yinch

        ORa = plt.axes( [0.75, 0.15, relX, relY], axisbg=cfg.axisbg)
        ORa.patch.set_edgecolor("darkgrey")
        ORa.patch.set_lw(1.)
        ORa.set_xticks([])
        ORa.set_yticks([])


        dX = 0.5*relX*cfg.xinch*self.scl
        dY = 0.5*relY*cfg.yinch*self.scl

        ORax = [posOR[cfg.X] - dX, posOR[cfg.X] + dX,
                posOR[cfg.Y] - dY, posOR[cfg.Y] + dY]
        ORa.grid(True)

        self._plotParts(ORa, ORax)
      

    if task == "ecc":
      self._plotEcc()
  
    if task == "energy":
      self._plotEnergy()
    
    self._plotMisc()
    self._plotPost()
    self._figSave(ifile)


  def _load(self,pfile,cfile):
    cfg = self.cfg
    if cfg.verbose:
      print "loading particles ..."
    pdumpf = H5PartDump(pfile)
    sname = (pdumpf.getStepNames())[0]

    pdump = pdumpf.getStep(sname)
    self.nop = (pdump.m.shape)[0]

    if cfg.verbose:
      print "load time and G ..."
    self.G     = pdumpf.getAttr(sname,"gravconst")
    self.time  = pdumpf.getAttr(sname,"time")
    
    if cfg.verbose:
      print "loading clumps ..."
    
    cdumpf = H5PartDump(cfile)
    cdump = cdumpf.getStep(sname)
    self.noc = (cdump.m.shape)[0]

    self.pdumpf = pdumpf
    self.pdump  = pdump
    
    self.cdumpf = cdumpf
    self.cdump  = cdump

  def _plotEcc(self):
    cfg = self.cfg
    pax = self.physa
    nax = self.norma
    
    ecca = plt.axes( [0.10, 0.05, 0.86, 0.75], axisbg=cfg.axisbg, frameon=False)
    
    ecca.xaxis.set_major_formatter(expm_formatter)
    ecca.xaxis.set_minor_formatter(expm_formatter)
    ecca.xaxis.tick_top()
    ecca.yaxis.set_major_formatter(math_formatter)
    ecca.yaxis.set_minor_formatter(math_formatter)
    ecca.grid(True, lw=0.1, color='darkgrey', alpha=0.5)
    self.ecca = ecca

    nop = self.nop
    noc = self.noc

    pdump = self.pdump
    cdump = self.cdump
    
    if cfg.verbose:
      print "selecting particles"
    filt = ( pdump.orbit[:,0] != 1 ) & ( pdump.clumpid[:,0] == 1 ) 

    if cfg.verbose:
      print "filtering particles, using ", sum(filt), "/", nop
    pos = pdump.pos[:,:].compress(filt, axis=0)[:,:]
    mat = pdump.mat[:,0].compress(filt)[:]
    bod = np.int32( (pdump.id[:,0].compress(filt)[:] > cfg.idfilt) )
    ecc = pdump.ecc[:,0].compress(filt)[:]
   
    rclmp = cdump.rclmp[1,0]
    rmean = cdump.rc[1,0]
    rhoclmp = cdump.rho[1,0]
    rhosat  = cfg.rhosat
    #rrche = 2.9*rclmp
    rrche = 2.455*rmean*pow( (rhoclmp/rhosat), 1./3. )

    com = cdump.posclmp[1,:]
    rvec = com - pos
    r   = np.sqrt( (rvec*rvec).sum(axis=1) )
    
    if cfg.verbose:
      print "get point colors and size"
    (pcol, pt2size) = cfg.colorFunc(pdump, cdump, filt, cfg)
    
    if cfg.verbose:
      print "plotting points ...   "
    cfg.cf = ecca.scatter( r, ecc, cfg.sctsize, pcol, lw=0, vmin=0., vmax=1., cmap=cfg.cmap)
    
    rmax = cfg.rmax
    rmin = -0.01*rmax
    r    = np.arange( rmin, rmax, 0.01*(rmax-rmin) )
    
    ecca.plot( [rclmp, rclmp], [0., 1.05], lw=0.3, color=cfg.txtc, ls='-' ) 
    ecca.plot( [0, 0], [0., 1.1], lw=0.3, color=cfg.txtc, ls=':' ) 
    ecca.plot( [rrche, rrche], [0., 1.05], lw=0.3, color='g', ls=':' ) 
    ecca.text( rclmp, 1.05, r'$r_{\mathrm{surf} }$', color=cfg.txtc, size=cfg.parm_txts )
    
    ecca.text( rrche, 1.05, r'$r_{\mathrm{Roche} }$', color='g', size=cfg.parm_txts )

    ecca.plot( r, (r-rrche)/(r+rrche), lw=0.15, color='g', ls='-' ) 
    ecca.plot( r, (r-rclmp)/(r+rclmp), lw=0.15, color=cfg.txtc, ls='-' ) 

    ecca.axis( [rmin, rmax, 0., 1.1] )

    for lbl in ecca.get_xticklabels():
      lbl.set_color(cfg.txtc)
      lbl.set_size(cfg.parm_txts)
    
    for lbl in ecca.get_yticklabels():
      lbl.set_color(cfg.txtc)
      lbl.set_size(cfg.parm_txts)

    ecca.set_xlabel('$\mathrm{radius}$', color=cfg.txtc, size=cfg.parm_txts)
    ecca.set_ylabel('$\mathrm{eccentricity}$', color=cfg.txtc, size=cfg.parm_txts)
    ecca.xaxis.set_label_position('top')
    

  def _plotEnergy(self):
    cfg = self.cfg
    pax = self.physa
    nax = self.norma
    
    enga = plt.axes( [0.10, 0.05, 0.86, 0.75], axisbg=cfg.axisbg, frameon=False)
    
    enga.xaxis.set_major_formatter(expm_formatter)
    enga.xaxis.set_minor_formatter(expm_formatter)
    enga.xaxis.tick_top()
    enga.yaxis.set_major_formatter(math_formatter)
    enga.yaxis.set_minor_formatter(math_formatter)
    enga.grid(True, lw=0.1, color='darkgrey', alpha=0.5)
    self.enga = enga

    nop = self.nop
    noc = self.noc

    pdump = self.pdump
    cdump = self.cdump
   
    filt = np.ones([nop])
   
    rclmp = cdump.rclmp[1,0]
    rmean = cdump.rc[1,0]
    rhoclmp = cdump.rho[1,0]

    rcom = cdump.posclmp[1,:]
    vcom = cdump.vel[1,:]
    
    rvec = rcom - pdump.pos[:,:]
    vvec = vcom - pdump.vel[:,:]
    
    #ekin = 0.5*(vvec*vvec).sum(axis=1)
    #epot = pdump.pot[:,0]
    r    = np.sqrt( (rvec*rvec).sum(axis=1) )
    vesc = np.sqrt( -2.*pdump.pot[:,0] )
    vz   = ( ( rvec*vvec ).sum(axis=1) ) / r
    
    for lbl in enga.get_xticklabels():
      lbl.set_color('white')
      lbl.set_size(cfg.parm_txts)
    
    for lbl in enga.get_yticklabels():
      lbl.set_color('white')
      lbl.set_size(cfg.parm_txts)

    
    if cfg.verbose:
      print "get point colors and size"
    (pcol, pt2size) = cfg.colorFunc(pdump, cdump, filt, cfg)
    
    if cfg.verbose:
      print "plotting points ...   "
    cfg.cf = enga.scatter( r, vz / vesc, 0.01, pcol, lw=0, vmin=0., vmax=1., cmap=cfg.cmap)
    
    rmax = cfg.rmax
    rmin = -0.01*rmax
    enga.axis( [rmin, rmax, -1.5, 2.0] )

  def _prepPlotParts(self):
    cfg = self.cfg
    pax = self.physa
    nax = self.norma

    nop = self.nop
    noc = self.noc

    pdump = self.pdump
    cdump = self.cdump
    
    cposz = np.zeros( max( cdump.id[:,0] ) + 1)
    for i in range(noc):
      cposz[ int(cdump.id[i]) ] = cdump.pos[i,cfg.Z]
    
    if cfg.verbose:
      print "determine mean h ..."
    
    if not np.isnan(cfg.hmed):
      self.hmed = cfg.hmed
    elif not cfg.plotGas:
      self.hmed = np.median(pdump.h[:,0].compress( pdump.mat[:,0] > 0 )[:] )
    else:
      self.hmed = np.median(pdump.h[:,0])
    hmed = self.hmed

    if cfg.verbose:
      print "selecting particles, filter is:",cfg.filt
    filt = np.ones(nop, dtype=np.bool)

    if cfg.filt == "negzonly":
      filt = pdump.pos[:,cfg.Z] < 0.
    
    if cfg.filt == "clumpcut":
      for i in range(nop):
        cidi = int( pdump.clumpid[i,0] )
        if (cidi == 0 or cidi >= noc):
          filt[i] = True
        else:
          cidi = int( pdump.clumpid[i,0] )
          zrel = pdump.pos[i,cfg.Z] - cposz[cidi]
          filt[i] = ( ( zrel > -hmed ) & ( zrel < hmed ) )
    
    if cfg.filt == "clumpnegz":
      # if there is no clump (noc = 0 or 1), then fall back to "negzonly"
      if noc < 2:
        filt = pdump.pos[:,cfg.Z] < 0.
      else:
        for i in range(nop):
          cidi = pdump.clumpid[i,0]
	  if (cidi > 0 and cidi < noc ):
            filt[i] = ( pdump.pos[i,cfg.Z] < cposz[cidi] )
          else:
            filt[i] = True

    if cfg.filt == "slice":
      filt = ( pdump.pos[:,cfg.Z] > -hmed ) & ( pdump.pos[:,cfg.Z] <  hmed ) & ( pdump.mat[i,0] > 0)

  
    if cfg.filt == "aux":
      filt = cfg.filtFunc(pdump, cdump, filt, cfg)
    
    if not cfg.plotGas:
      if cfg.verbose:
        print "removing gas particles"
      for i in range(nop):
        if ( pdump.mat[i,0] == 0 ):
          filt[i] = False
        
    if cfg.verbose:
      print "filtering particles, using ", sum(filt), "/", nop
    pos = pdump.pos[:,:].compress(filt, axis=0)[:,:]
    mat = pdump.mat[:,0].compress(filt)[:]
    bod = np.int32( (pdump.id[:,0].compress(filt)[:] > cfg.idfilt) )

    # auto determine the scatter ptsize^2 (no idea why we need a fudge factor)
    cfg.pt2size = 0.008*np.power( hmed*cfg.dpi / self.scl , 2.)

    self.pos  = pos
    self.filt = filt

    if cfg.verbose:
      print "z-sorting points ...   "
    sidx = pos[:,cfg.Z].argsort()
    self.sidx    = sidx
    
    if cfg.verbose:
      print "get point colors and size"
    (pcol, pt2size) = cfg.colorFunc(pdump, cdump, filt, cfg)

    self.pt2size = pt2size
    self.pcol    = pcol
    

  def _plotParts(self,cax,caxis):
    cfg = self.cfg

    noc = self.noc
    
    pos = self.pos
    sidx = self.sidx
    
    pt2size = self.pt2size
    pcol    = self.pcol   

    pdump = self.pdump
    cdump = self.cdump
    
    if cfg.plottraj:
      clumps = FewBodies( cdump, G, cfg.Mmin)
      if cfg.verbose:
        print "integrate clumps ...   "
      clumps.integrate(dt, cfg.ds)
    
      if cfg.verbose:
        print "plotting clumps with trajectories ... "
      for i in range(1,noc):
        curtraj = clumps.traj[i]
        if ( clumps.m[i] > cfg.MminPlt ):
          cax.add_patch( mp.patches.Circle((curtraj[0,cfg.X], curtraj[0,cfg.Y]),\
            radius=clumps.rc[i], ec=cfg.clpc_ec, fc=cfg.clpc_fc, \
            lw=cfg.clpc_lw, alpha=cfg.clpc_al) )
    
          if ( clumps.m[i] > cfg.MminLbl ):
            cax.text( curtraj[0,cfg.X], curtraj[0,cfg.Y], \
                '$\mathrm{'+ '%1.4f' % ( clumps.m[i] / cfg.ME ) +r' M_{\oplus}}$', \
                size=cfg.clpc_txts, color=cfg.txtc)
    
          cax.add_patch( \
              mp.patches.Circle((curtraj[-1,cfg.X], curtraj[-1,cfg.Y]), \
              radius=clumps.rc[i], ec=cfg.clpp_ec, fc=cfg.clpp_fc, \
              lw=cfg.clpp_lw, alpha=cfg.clpp_al) )
          cax.plot( curtraj[:,cfg.X], curtraj[:,cfg.Y], color=cfg.clpt_fc, \
              lw=cfg.clpt_lw, alpha=cfg.clpt_al)
    else:
      if cfg.verbose:
        print "plotting clumps ... "
      for i in range(1,noc):
        if ( cdump.m[i,0] > cfg.MminPlt ):
          
          rclmp = cdump.rclmp[i,0]
          rmean = cdump.rc[i,0]
          rhoclmp = cdump.rho[i,0]
          rhosat  = cfg.rhosat

          rrche   = 2.455*rmean*pow( (rhoclmp/rhosat), 1./3. )

          cax.add_patch( mp.patches.Circle((cdump.pos[i,cfg.X], \
              cdump.pos[i,cfg.Y]), radius=cdump.rc[i,0], \
              ec=cfg.clpp_ec, fc=cfg.clpp_fc, lw=cfg.clpp_lw, \
              alpha=cfg.clpp_al) )
          
          if ( rclmp > 0. ):
            if cfg.clp_plotsurf:
              cax.add_patch( mp.patches.Circle((cdump.posclmp[i,cfg.X], \
                  cdump.posclmp[i,cfg.Y]), radius=rclmp, \
                  ec=cfg.clpp_ec, fc=cfg.clpp_fc, lw=cfg.clpp_lw, \
                  alpha=cfg.clpp_al, ls='dotted') )
            
            if cfg.clp_plotrche:
              cax.add_patch( mp.patches.Circle((cdump.posclmp[i,cfg.X], \
                  cdump.posclmp[i,cfg.Y]), radius=rrche, \
                  ec=cfg.clpr_ec, fc=cfg.clpp_fc, lw=cfg.clpp_lw, \
                  alpha=cfg.clpp_al, ls=cfg.clpr_ls ) )
          
          if ( cdump.m[i] > cfg.MminLbl ):
            # custom function
            cax.text( cdump.pos[i,cfg.X], cdump.pos[i,cfg.Y], \
                '$\mathrm{'+ '%1.4f' % ( cdump.m[i,0] / cfg.ME ) +r' M_{\oplus}}$', \
                size=cfg.clpc_txts, color=cfg.txtc)


    if cfg.verbose:
      print "plotting points ...   "
    cfg.cf = cax.scatter( pos[sidx,cfg.X], pos[sidx,cfg.Y], pt2size, \
        pcol[sidx,:], lw=0, vmin=0., vmax=1., cmap=cfg.cmap)

    cax.axis("scaled")
    cax.axis(caxis)
    if cfg.verbose:
      print "corrected axis: ",self.corrax

  def _plotPre(self):
    cfg = self.cfg
    if cfg.verbose:
      print "pre-plot ... "
    for func in cfg.prePlot:
      func(self.norma, self.physa, plt, self.pdump, self.cdump, cfg)

  def _plotPost(self):
    cfg = self.cfg
    if cfg.verbose:
      print "post-plot ... "
    for func in cfg.postPlot:
      func(self.norma, self.physa, plt, self.pdump, self.cdump, cfg)

  def _figSave(self,ifile):
    cfg = self.cfg
    if cfg.verbose:
      print "save figure     ...   "
    plt.rc('savefig', dpi=self.cfg.dpi)
    plt.savefig(ifile)

  def _plotMisc(self):
    nax = self.norma
    cfg = self.cfg
    fig = self.fig

    # plot a scale
    if cfg.plotscale:
      if cfg.verbose:
        print "plot a scale    ...   "
      nax.plot( [cfg.scale_vc[0], cfg.scale_vc[2]], \
          [cfg.scale_vc[1], cfg.scale_vc[3]], lw=.3, color=cfg.txtc )
      
      scldst = (cfg.scale_vc[2] - cfg.scale_vc[0])*( self.corrax[1] - \
          self.corrax[0] )*cfg.scale_sc
      scltxt = latexExp10( '$%6.1e' % (scldst) )\
          + cfg.scale_ut + '$'
      nax.text( cfg.scale_vc[0], cfg.scale_vc[1] + 0.01, scltxt, \
          color=cfg.txtc, size=cfg.scale_txts)
    
    nax.axis([0., 1., 0., 1.])
    
    if cfg.verbose:
      print "plot time ..."
    timetxt = '$t = %6.2f' % (self.time*cfg.time_sc) + ' ' + cfg.time_ut + '$'
    nax.text( cfg.time_vc[0], cfg.time_vc[1], timetxt, \
        color=cfg.txtc, size=cfg.time_txts)
    
    if cfg.verbose:
      print "parameters and copyright ..."
    nax.text( cfg.parm_vc[0], cfg.parm_vc[1], cfg.parm_txt, \
        color=cfg.txtc, size=cfg.parm_txts )
    
    nax.text( cfg.copy_vc[0], cfg.copy_vc[1], cfg.copy_txt, \
        color=cfg.txtc, size=cfg.copy_txts )
    
    if cfg.cbar_plot:
      if cfg.verbose:
        print "plot colorbar   ...   "
      cbax = plt.axes( cfg.cbar_extent, frameon=False)
      cbar = fig.colorbar(cfg.cf, cax=cbax, orientation=cfg.cbar_orientation, \
          ticks=cfg.cbar_ticksx, drawedges=False)

      if cfg.cbar_orientation == 'horizontal':
        cbax.set_xticklabels(cfg.cbar_tickstxt)
        for t in cbax.get_xticklabels():
          t.set_color(cfg.txtc)
          t.set_size(cfg.cbar_txts)
      else:
        cbax.set_yticklabels(cfg.cbar_tickstxt)
        for t in cbax.get_yticklabels():
          t.set_color(cfg.txtc)
          t.set_size(cfg.cbar_txts)

      cbar.outline.set_lw(cfg.cbar_lw)



    
    

  




