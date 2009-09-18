#!/usr/bin/env ipython
import matplotlib
matplotlib.use('Agg')

from pylab import *


filename = sys.argv[1]
timeraw = sys.argv[2]
zoom     = float(sys.argv[3])
partsize = float(sys.argv[4])

timems = float(timeraw)*1000.
timestr = '$\\mathrm{ t = ' + '%6.2f' % ( timems ) + ' ms }$'
ptsize = partsize*zoom;

slice = load(filename)
noparts = len( slice )

def retnullstring(x, pos):
  return ''

empty_formatter = FuncFormatter(retnullstring)

rc('grid', color='grey')
rc('text', usetex=True)
rc('savefig', dpi=75)

figure( figsize=(10, 5) )

# dummy scatter plot, so that we can get a colorbar
axes( [0.067, 0.13, 0.40, 0.73], axisbg='w')
scatter( slice[:,0], slice[:,1], ptsize, slice[:,5], linewidth=0)
clim(0.75,1)

xlabel("$x [cm]$")
ylabel("$z [cm]$")
title("$damage$")
axis("scaled")
axis([-600/zoom,600/zoom,-700/zoom,400/zoom])
gca().yaxis.tick_left()
gca().yaxis.set_label_position("left")
annotate(timestr, [-560/zoom, 300/zoom])
annotate('$1m$', [-540/zoom, -440/zoom], color='white')
gca().add_patch( Arrow(-550./zoom,-450./zoom, 100./zoom, 0., ec='w', fc='w', lw=0.001, width=0.) )
gca().add_patch( Arrow(-550./zoom,-450./zoom, 0., 100./zoom, ec='w', fc='w', lw=0.001, width=0.) )

gca().set_aspect("equal")

# plot colorbar
cba = axes( [0.23, 0.75, 0.20, 0.040], axisbg='w')
cb = colorbar( cax = cba , orientation='horizontal', ticks=[0.75, 0.80, 0.85, 0.90, 0.95, 1.00])
cb.ax.set_xticklabels(['$ < 75\%$', '', '', '$90\%$', '', '$100\%$',])# horizontal colorbar

logp = zeros(noparts, 'f')
pcol = [[]]*noparts
for i in range(noparts):
  if slice[i,10] > 0.:
    logp[i] = log10( slice[i,10] / 10. )
  else:
    logp[i] = -100.

  if slice[i,11] > 3.e10:
    pcol[i] = 'r'
  elif slice[i,11] > 1.e10:
    pcol[i] = 'g'
  else:
    pcol[i] = 'b'

  if slice[i,6] > 3.:
    pcol[i] = 'grey'

# dummy scatter plot, so that we can get a colorbar
axes( [0.533, 0.533, 0.40, 0.33], axisbg='w', frameon=True)
scatter( slice[:,0], slice[:,1], ptsize, pcol, linewidth=0)
title("$ice~peak~pressure$")
axis("scaled")
axis([-600/zoom,600/zoom,-500/zoom,0/zoom])
clim(7., 9.5)
annotate('$\mathrm{red:}$',     [150/zoom, -320/zoom], color='white')
annotate('$\mathrm{green:}$',   [150/zoom, -370/zoom], color='white')
annotate('$\mathrm{blue:}$',    [150/zoom, -420/zoom], color='white')
annotate('$\mathrm{grey:}$',    [150/zoom, -470/zoom], color='white')
annotate('$> 3.0~\mathrm{GPa}$', [300/zoom, -320/zoom], color='white')
annotate('$> 1.0~\mathrm{GPa}$', [300/zoom, -370/zoom], color='white')
annotate('$< 1.0~\mathrm{GPa}$', [300/zoom, -420/zoom], color='white')
annotate('$\mathrm{dunite}$',    [300/zoom, -470/zoom], color='white')
annotate('$1m$', [-540/zoom, -440/zoom], color='white')
gca().add_patch( Arrow(-550./zoom,-450./zoom, 100./zoom, 0., ec='w', fc='w', lw=0.001, width=0.) )
gca().add_patch( Arrow(-550./zoom,-450./zoom, 0., 100./zoom, ec='w', fc='w', lw=0.001, width=0.) )
gca().set_xticks([],[])
gca().set_yticks([],[])


axes( [0.533, 0.133, 0.40, 0.33], axisbg='w', frameon=True)
scatter( slice[:,0], slice[:,1], ptsize, logp, linewidth=0)
title("$pressure$")
axis("scaled")
axis([-600/zoom,600/zoom,-500/zoom,0/zoom])
clim(7., 9.5)
annotate('$1m$', [-540/zoom, -440/zoom], color='white')
gca().add_patch( Arrow(-550./zoom,-450./zoom, 100./zoom, 0., ec='w', fc='w', lw=0.001, width=0.) )
gca().add_patch( Arrow(-550./zoom,-450./zoom, 0., 100./zoom, ec='w', fc='w', lw=0.001, width=0.) )
gca().set_xticks([],[])
gca().set_yticks([],[])

cba = axes( [0.70, 0.07, 0.20, 0.04], axisbg='w')
cb = colorbar( cax = cba , orientation='horizontal', ticks=[7., 8., 9., 9.4771])
cb.ax.set_xticklabels(['$<0.01$', '$0.1$', '$1$', '$3\mathrm{~GPa}$ ',])

savefig(filename + ".png")

