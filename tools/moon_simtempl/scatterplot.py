#!/usr/bin/python

execfile("config.sh")

import matplotlib
matplotlib.use('Agg')

from matplotlib.patches import Circle

from pylab import *
filename = sys.argv[1]
clmpname = sys.argv[2]
timesec = float(sys.argv[3])
timestring = '$\\mathrm{ t = ' + '%6.2f' % ( timesec/3600 ) + ' h }$'

plotline1 = '$\\mathrm{mass~ratio~' + '%1.2f' % ( MIMP/MTAR ) + '}$'
plotline2 = '$\\mathrm{total~mass~' + '%1.2f' % ( MIMP+MTAR ) + '~M_e}$'
plotline3 = '$\\mathrm{L_{tot}~=~'  + '%1.2f' % ( LTOT      ) + '~lem}$'

rc('grid', color='grey')
rc('text', usetex=True)
rc('savefig', dpi=300)

figure( figsize=(3, 2) )
a = axes( [0.0, 0.0, 1.0, 1.0], axisbg='k')

slice = loadtxt(filename)
clump = loadtxt(clmpname)

# slice[:,2] is 0 for mantle and 1 for iron
# slice[:,3] is 0 for protoearth and 1 for impactor
mat = slice[:,2].astype(int)
bod = slice[:,3].astype(int)

map = zeros((2,6), float)
map[0,2] = 0.30 # ice
map[0,4] = 0.95 # dun
map[0,5] = 0.05 # iron
map[1,2] = 0.35
map[1,4] = 0.90
map[1,5] = 0.10

color = map[ bod, mat ]

scatter( slice[:,0], slice[:,1], 0.10, color, linewidth=0)
ca().add_patch( Circle( ( clump[0,0], clump[0,1]) , radius=clump[0,2], ec='yellow', fc='yellow', lw=1., ) )


#scatter( slice[:,0], slice[:,1], 0.10, color, linewidth=0)
clim(0., 1.)
axhline(color='w', ls=':', lw=0.10, alpha=0.2)
axvline(color='w', ls=':', lw=0.10, alpha=0.2)
gca().add_patch( Arrow(-6.4e9,-4.1e9, 6.371e8, 0., alpha=1.0, ec='r', fc='g', lw=0.001, width=0.) )
annotate('$\mathrm{r_{earth}}$', [-7.08e9, -4.2e9], color='w', size=4)
gca().add_patch( Arrow(-6.4e9,-4.5e9, 3.841e9, 0., alpha=1.0, ec='r', fc='g', lw=0.001, width=0.) )
annotate('$ 0.1 \mathrm{a_{moon}}$', [-7.5e9, -4.6e9], color='w', size=4)
annotate('$\mathrm{Giant~impact}$', [-7.50e9, 2.30e9], color='w', size=8)
annotate(timestring,                [-7.50e9, 1.80e9], color='w', size=8)
annotate(plotline1, [ 1.40e9, 2.50e9], color='w', size=4)
annotate(plotline2, [ 1.40e9, 2.20e9], color='w', size=4)
annotate(plotline3, [ 1.40e9, 1.90e9], color='w', size=4)
annotate('$\mathrm{Andreas~Reufer}$',       [ 1.40e9,-4.30e9], color='w', size=4)
annotate('$\mathrm{University~of~Bern}$',   [ 1.40e9,-4.60e9], color='w', size=4)

grid(False)
xlabel("$\\mathrm{x [cm]}$")
ylabel("$\\mathrm{y [cm]}$")
#title(timestring)

#axis("scaled")
axis(FULLSCALE)
gca().set_aspect("equal")
savefig(filename + "_full.png")

gca().clear();

scatter( slice[:,0], slice[:,1], 0.20, color, linewidth=0)
clim(0., 1.)
axhline(color='w', ls=':', lw=0.10, alpha=0.2)
axvline(color='w', ls=':', lw=0.10, alpha=0.2)
gca().add_patch( Arrow(-3.2e9,-2.05e9, 6.371e8, 0., alpha=1.0, ec='r', fc='g', lw=0.001, width=0.) )
annotate('$\mathrm{r_{earth}}$', [-3.54e9, -2.1e9], color='w', size=4)
gca().add_patch( Arrow(-3.2e9,-2.25e9, 3.841e9, 0., alpha=1.0, ec='r', fc='g', lw=0.001, width=0.) )
annotate('$ 0.1 \mathrm{a_{moon}}$', [-3.75e9, -2.3e9], color='w', size=4)
annotate('$\mathrm{Giant~impact}$', [-3.75e9, 1.15e9], color='w', size=8)
annotate(timestring,[-3.75e9, 0.90e9], color='w', size=8)
annotate(plotline1, [ 0.70e9, 1.25e9], color='w', size=4)
annotate(plotline2, [ 0.70e9, 1.10e9], color='w', size=4)
annotate(plotline3, [ 0.70e9, 0.95e9], color='w', size=4)
annotate('$\mathrm{Andreas~Reufer}$',       [ 0.70e9,-2.15e9], color='w', size=4)
annotate('$\mathrm{University~of~Bern}$',   [ 0.70e9,-2.30e9], color='w', size=4)
grid(False)
xlabel("$\\mathrm{x [cm]}$")
ylabel("$\\mathrm{y [cm]}$")
#title(timestring)

axis("scaled")
axis(ZOOMSCALE)
savefig(filename + "_zoom.png")
