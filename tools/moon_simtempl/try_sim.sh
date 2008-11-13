#!/bin/bash

if [ ! -e config.sh ]
then
  echo "config.sh not found!"
  exit 1
fi

source config.sh

# setup_gi.sh
# PARAMFILE=`tempfile`
PARAMFILE=params
echo $MTAR >>$PARAMFILE
echo $MIMP >>$PARAMFILE
echo $RTAR >>$PARAMFILE
echo $RIMP >>$PARAMFILE
echo $VINF >>$PARAMFILE
echo $LTOT >>$PARAMFILE
echo $RSEP >>$PARAMFILE
./setup_gi_show.m < $PARAMFILE >/dev/null
#rm $PARAMFILE

cat setup_gi.sh | head -n36 | tail -n35
rm setup_gi.sh

