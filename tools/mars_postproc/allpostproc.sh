#!/bin/bash

SLEEPTIME=60
BASEDIR=`pwd`

while true; do
  DIRS=`find -L . -maxdepth 2 -name postproc.sh|sort|sed -e 's/\/postproc\.sh//g'|sort`

  for CURDIR in $DIRS
  do
    echo
    echo " try $CURDIR"
    echo
    cd $CURDIR
    nice ./postproc.sh
    cd $BASEDIR
    echo
    echo " done $CURDIR"
    echo
  done

  echo
  echo " sleep for ${SLEEPTIME}s ..."
  echo
  sleep $SLEEPTIME
done

