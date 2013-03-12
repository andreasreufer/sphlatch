#!/usr/bin/env ipython
import tables as pt
import numpy as np

ifile = pt.openFile("mod_IG_Mc10_rho2881.h5")
ofile = pt.openFile("profcgs.hdf5", mode="w")

ofile.createArray( ofile.root, "r", ifile.root.R[:]*1.e2)
ofile.createArray( ofile.root, "rho", ifile.root.density[:]*1.e-3 )
ofile.createArray( ofile.root, "u", ifile.root.energy[:]*1.e4 )
ofile.createArray( ofile.root, "p", ifile.root.P[:]*1.e1 )

ofile.root._v_attrs["gamma"] = 1.28571

ofile.close()
ifile.close()

