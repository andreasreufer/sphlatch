# the simulations name and directory
SIMNAME="can030"
SIMDIR=~/moon_impact/sim_${SIMNAME} 
NOCPU=4

# sources
SRCDIR=~/repos/sphlatch/apps/impact
MAKETARG="impact_GHUA__"
BINARY="impact_GHUA__"
EXTRAFILES="aneos.input"

# bodies
TARGFILE=~/moon_impact/bodies_03rlx/body_A_0.965me_T2000K_0.3iron_0.7dun_018k.h5part
IMPCFILE=~/moon_impact/bodies_03rlx/body_A_0.144me_T2000K_0.3iron_0.7dun_002k.h5part

# initial condition parameters
MTAR=0.965
MIMP=0.144
RTAR=0.985
RIMP=0.549
VINF="4.293e5"
LTOT=1.34
RSEP=5.
UMIN="1.e8"   # minimal spec. energy
MAXR="7.e9"  # minimal radius to proto-earth for removal of part.

# scale of scatter plots and max. binning radius for mass profiles
FULLSCALE="[-10.0e9,4.0e9,-5.0e9,5.0e9]"
ZOOMSCALE="[-2.0e9,2.0e9,-1.5e9,1.5e9]"
PPMAXRAD="7.e9"

# times
SAVETIME=1e2
STOPTIME=1.296e5
