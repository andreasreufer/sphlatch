# the simulations name and directory
#SIMNAME="can030"
#SIMDIR=~/moon_impact/sim_${SIMNAME} 
#NOCPU=8
SIMNAME="fuc04c"
SIMDIR="~/pmc_colls/sim_fuc04c"
NOCPU=8

# sources
#SRCDIR=~/repos/sphlatch/apps/impact
#MAKETARG="impact_GHUA__"
#BINARY="impact_GHUA__"
SRCDIR="~/repos/sphlatch/apps/simple_sph/"
MAKETARG="simple_sph_GHUa__"
BINARY="simple_sph_GHUa__"
EXTRAFILES="aneos.input"

# bodies
TARGFILE="~/moon_impact/bodies_03rlx/body_A_0.965me_T2000K_0.3iron_0.7dun_174k.h5part"
IMPCFILE="~/moon_impact/bodies_03rlx/body_A_0.144me_T2000K_0.3iron_0.7dun_026k.h5part"

# initial condition parameters
MTAR=0.965
MIMP=0.144
RTAR=1.014
RIMP=0.578
VINF=0.000e5
LTOT=1.12
RSEP=3.
UMIN=1.e8   # minimal spec. energy
MAXR=7.e9  # minimal radius to proto-earth for removal of part.

# scale of scatter plots and max. binning radius for mass profiles
FULLSCALE=[-10.0e9,4.0e9,-5.0e9,5.0e9]
ZOOMSCALE=[-2.0e9,2.0e9,-1.5e9,1.5e9]
PPMAXRAD="7.e9"

# times
SAVETIME=1e3
STOPTIME=8.64e4
