#!/usr/bin/env python

# UBELIX:
#qsub_sgl = "qsub -cwd -M andreas.reufer@space.unibe.ch -l h_cpu=1:00:00 -l h_vmem=512M -N $JOBNAME $JOBCMD"
#qsub_omp = "qsub -M andreas.reufer@space.unibe.ch -b y -cwd -l h_cpu=240:00:00 -l h_vmem=512M -pe smp $NOCPUS -R y -N $SIMNAME ./$BINARY initial.h5part $SAVETIME $STOPTIME $RUNARGS"
#qdel     = "qdel $JOBID"
#basedir  = "/home/ubelix/space/areufer/"
#srcdir   = basedir + "repos/sphlatch/apps/simple_sph"

# ISIS2
#qsub_sgl = "qsubSALL $JOBNAME $JOBCMD"
#qsub_omp = "qsubPSMP $NOCPUS $SIMNAME ./$BINARY initial.h5part $SAVETIME $STOPTIME $RUNARGS"
#qdel     = "qdel $JOBID"
#basedir  = "/home0/areufer/"
#srcdir   = basedir + "repos/sphlatch/apps/simple_sph"

# caliban
qsub_sgl = ""
qsub_omp = ""
qdel     = ""
basedir  = "/Users/areufer/Documents/UniTAPS/gi_sims/"
srcdir   = "/Users/areufer/repos/sphlatch/apps/simple_sph"





