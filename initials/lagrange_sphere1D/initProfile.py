#!/usr/bin/env ipython
import tables as pt
import numpy as np

me = 5.9742e27
mtot = 2*me
eVinK = 11604.505;
filename = "in.hdf5"

nocells = 1000
dm = mtot / nocells

m   = np.zeros( (nocells,1) )
p   = np.zeros( (nocells,1) )
rho = np.zeros( (nocells,1) )
S   = np.zeros( (nocells,1) )
mat = np.zeros( (nocells,1), dtype=np.int32 )

ironrng = range(  0, 120)
ironmat = 5
ironrho =  8.0
ironS   = 1.60e11

sio2rng = range(120, 500)
sio2mat = 1
sio2rho = 2.00
sio2S   = 4.03e11

watrrng = range(500,1000)
watrmat = 2
watrrho = 0.80
watrS   = 4.64e11

m[:,0] = dm

S[ironrng,0] = ironS
S[sio2rng,0] = sio2S
S[watrrng,0] = watrS

mat[ironrng,0] = ironmat
mat[sio2rng,0] = sio2mat
mat[watrrng,0] = watrmat

rho[ironrng,0] = ironrho
rho[sio2rng,0] = sio2rho
rho[watrrng,0] = watrrho

# chris profile mat 1
#S[:,0] = 7.67966e11
#mat[:,0] = 1
#rho[:,0] = 2.65
#p[nocells-1,0] = 2.5570e10

# chris profile mat 4
#S[:,0] = 7.6534e11
#mat[:,0] = 4
#rho[:,0] = 1.5
#p[nocells-1,0] = 7.031e10


h5file = pt.openFile(filename, mode = "w")

h5file.createArray( h5file.root, "m",   m )
h5file.createArray( h5file.root, "p",   p )
h5file.createArray( h5file.root, "rho", rho )
h5file.createArray( h5file.root, "S",   S )
h5file.createArray( h5file.root, "mat", mat )

h5file.root._v_attrs["gravconst"] = 6.67e-08

h5file.close()


