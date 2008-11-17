#!/bin/bash

if [ ! -e config.sh ]
then
  echo "config.sh not found!"
  exit 1
fi

source config.sh

TMPLDIR=`pwd`

# create simdir
if [ -e $SIMDIR ]
then
  echo "$SIMDIR already exists!"
  exit 1
fi
mkdir $SIMDIR

# copy everything to new simdir
cp -Rpv bf_nbody.m config.sh scatterplot.sh setup_gi.m setup_gi_show.m \
 $EXTRAFILES $SIMDIR >/dev/null
cd $SIMDIR

# copy body files to simdir/bodies
mkdir bodies
cp -Rpv $TARGFILE ./bodies/target.h5part
cp -Rpv $IMPCFILE ./bodies/impactor.h5part

# setup_gi.sh
PARAMFILE=`tempfile`
echo $MTAR >>$PARAMFILE
echo $MIMP >>$PARAMFILE
echo $RTAR >>$PARAMFILE
echo $RIMP >>$PARAMFILE
echo $VINF >>$PARAMFILE
echo $LTOT >>$PARAMFILE
echo $RSEP >>$PARAMFILE
#./setup_gi_show.m < $PARAMFILE >/dev/null
./setup_gi.m < $PARAMFILE >/dev/null
./setup_gi.sh
sleep 1
h5part_writeattr -i bodies/combined.h5part -k umin -v $UMIN
h5part_writeattr -i bodies/combined.h5part -k maxradius -v $MAXR
ln -s bodies/combined.h5part initial.h5part
rm $PARAMFILE bf_nbody.m setup_gi_show.m setup_gi.m


# prepare binary
ln -s $SRCDIR src
cd src
make $MAKETARG
cd ..
cp -Rpv src/$BINARY .

# prepare postproc script
mkdir gifs_full gifs_zoom postproc
cp -Rpv $TMPLDIR/_postproc.sh postproc.sh
echo "# time         Mpe          Lpe          Mdisk        Ldisk        Mesc         Lesc         Lrem         Ltot       DiskIronFrac" >simSummary

echo 
echo " start simulation with:"
echo "cd $SIMDIR && \\"
echo "qsubPOPTERON $NOCPU $SIMNAME $BINARY initial.h5part $SAVETIME $STOPTIME"
echo

