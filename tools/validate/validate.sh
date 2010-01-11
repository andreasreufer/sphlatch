#!/bin/bash

GOOD=$'\e[32;01m'
WARN=$'\e[33;01m'
BAD=$'\e[31;01m'
NORMAL=$'\e[0m'

source $1

BASEDIR=`pwd`

cd $SRCDIR
SRCDIR=`pwd`
make $MAKETARG
    
if [[ $? == 0 ]]; then
  echo "make $MAKETARG $GOOD succeeded $NORMAL"
else
  echo "make $MAKETARG $BAD failed $NORMAL"
fi

cp -Rpv $SRCDIR/$RUNBIN $RUNDIR/

cd $RUNDIR
RUNDIR=`pwd`

$RUNCMD >$RUNOUT
mv $RUNOUT $BASEDIR

if [[ $? == 0 ]]; then
  echo "$RUNCMD $MAKETARG $GOOD succeeded $NORMAL"
else
  echo "$RUNCMD $MAKETARG $BAD failed $NORMAL"
fi

cd $BASEDIR
exit 0

