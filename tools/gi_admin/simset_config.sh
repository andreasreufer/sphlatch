# definition of 
SIMSETBDIR="./test001"
BODIESDB="./bodies"
LOGFILE="simset.log"

TEMP="1500"

MIMP="0.01 0.1 1.0"
MTAR="0.01 0.1 1.0"
VIMPREL="1.0 1.5 2.0 2.5 3.0"
IMPA="0 30 45 60 75"
SIMCOND="(mimp <= mtar)"

BODTOLERANCE="0.1"

AUXFILES="./aux/aneos.input"
SRCDIR="~/repos/sphlatch/apps/simple_sph"
MAKETARG="simple_sph_GHUa__"
BINARY="simple_sph_GHUa__"

SUBCMD=""
RUNARGS="16"
NOCPUS="8"

RELSEP="3.0"
GRAVCONST="6.67429e-8"
ATTRKEYS="UMIN KEY"
ATTRVALS="1.e9 -4.e2"

