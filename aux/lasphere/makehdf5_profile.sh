#!/bin/bash

prepareConfigfile() {
echo "PATH $1                 "  >$CFGFILE
echo "INPUT-CLASS $2          " >>$CFGFILE
echo "INPUT-SIZE 64           " >>$CFGFILE
echo "RANK 1                  " >>$CFGFILE
echo "DIMENSION-SIZES $NOCELLS" >>$CFGFILE
echo "OUTPUT-CLASS FP         " >>$CFGFILE
echo "OUTPUT-SIZE 64          " >>$CFGFILE
echo "OUTPUT-ARCHITECTURE IEEE" >>$CFGFILE
echo "OUTPUT-BYTE-ORDER LE    " >>$CFGFILE
}

ascii2hdf5() {
sed -e 's/D/E/g' $INPUTFILE | cut -c $1 $RAWFILE >$ARRFILE
prepareConfigfile $2 $3
h5import $ARRFILE -c $CFGFILE -o $OUTPUTFILE
}

ARRFILE=`tempfile`
CFGFILE=`tempfile`

INPUTFILE=$1
OUTPUTFILE=$2

NOCELLS=`cat $INPUTFILE|wc -l`

ascii2hdf5  7-18 "r"   "TEXTFP"
ascii2hdf5 19-30 "v"   "TEXTFP"
ascii2hdf5 31-42 "p"   "TEXTFP"
ascii2hdf5 43-54 "q"   "TEXTFP"
ascii2hdf5 55-66 "rho" "TEXTFP"
ascii2hdf5 67-78 "u"   "TEXTFP"
ascii2hdf5 79-81 "mat" "TEXTIN"

rm $ARRFILE $CFGFILE
