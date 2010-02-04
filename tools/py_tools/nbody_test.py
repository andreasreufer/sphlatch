import tables
import nbody

fileh = tables.openFile("dump0010494_T03.0000e+04_clumps.h5part", "r")
clumps = fileh.root.current

#nbody.init( clumps )



