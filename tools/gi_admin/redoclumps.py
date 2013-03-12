import commands
import os
import stat

def redoClumps(sim):
  mclmp = 0.01*sim.mtot
  morbt = 0.10*sim.mtot

  cfile   = sim.dir + "clumps.h5part"
  dskfile = sim.dir + "disk.hdf5"

  print ("rm " + cfile)
  (exstat, out) = commands.getstatusoutput("rm " + cfile)
  jobname = sim.dir + "redoclumps.sh"
  if os.path.exists(jobname):
    (exstat, out) = commands.getstatusoutput("rm " + jobname)
  scriptfile = open(jobname,"w")
  
  print >>scriptfile, "#!/bin/bash"
  
  for dump in sim.dumps:
    (dfilerel, dtime) = dump
    dfile = sim.dir + dfilerel
    cmd = "h5part_writeattr -i " + dfile + " -k rhominclump -v 1.9875 "
    print >>scriptfile, cmd
    
    cmd = "find_clumps_A_ " + dfile + " " + cfile + " " + str(mclmp) + " " + str(morbt) + " 1"
    print >>scriptfile, cmd
    
    cmd = "moon_pp " + dfile + " " + dskfile + " 1 0. 2.e10 100 0.1 1.5 1.e23"
    print >>scriptfile, cmd
  
  print >>scriptfile, "rm " + jobname
  scriptfile.close()

  os.chmod(jobname, stat.S_IRWXU)
  print jobname
  jobsubstr = "qsubSALL $JOBNAME $JOBCMD"
  jobsubstr = "qsub -M andreas.reufer@space.unibe.ch -b n -cwd -l h_cpu=96:00:00 -l h_vmem=512M -N $JOBNAME $JOBCMD"

  jobsubstr = jobsubstr.replace("$JOBNAME", "redoclumps_"+sim.params.key)
  jobsubstr = jobsubstr.replace("$JOBCMD",  "redoclumps.sh")
  
  oldwd = os.getcwd()
  os.chdir(sim.dir)
  (exstat, out) = commands.getstatusoutput(jobsubstr)
  print jobsubstr
  os.chdir(oldwd)
  print out


