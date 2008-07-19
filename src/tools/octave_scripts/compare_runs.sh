#!/bin/bash

GOOD=$'\e[32;01m'
WARN=$'\e[33;01m'
BAD=$'\e[31;01m'
NORMAL=$'\e[0m'

REFDIR=$1
CMPDIR=$2

REFFILES=`cd $REFDIR && find . -maxdepth 1 -iname '*.h5part'|sed 's/\.\///'|sort`

OCTSCRIPT=`mktemp /tmp/XXXXXXX`

cat >$OCTSCRIPT<<EOT
arg_list = argv();
exitStat = 0;
diff = h5part_compare( load( arg_list{1}).current, load( arg_list{2}).current );
for [diffVal, diffKey] = diff
  if ( strncmp("r_", diffKey, 2) )
    maxDev = max(max( diffVal ) );
    if ( maxDev > 0. )
      printf(" %s: max relative difference %e\n", diffKey, maxDev)
      exitStat = 1;
    endif
  endif
endfor
exit(exitStat);
EOT

echo "$WARN   $REFDIR vs. $CMPDIR:$NORMAL"
for FILE in $REFFILES
do
  if [ -a $CMPDIR/$FILE ]; then
    octave -q $OCTSCRIPT $REFDIR/$FILE $CMPDIR/$FILE
    if [[ $? == 0 ]]; then
      echo "$GOOD$FILE $NORMAL equal"
    else
      echo "$BAD$FILE $NORMAL differs!"
    fi
    echo
  fi
done

rm $OCTSCRIPT
