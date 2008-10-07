#!/bin/bash

DELAY="5"
POVQUAL="7"
POVOPTS="+A +AM2"

#ANIMFILE="anim-norm.gif"
#SCALE="0.05"

ANIMFILE="anim-zoom.gif"
SCALE="0.15"

HEADER="cut_header.pov"
HDF52POVRAY="./hdf52povray_cut"

FILES=`find . -maxdepth 1 -iname '*.h5part'|sed 's/\.h5part//'|sed 's/\.\///'|sort `

mkdir $ANIMFILE-frames 2>&1 >/dev/null

noTasks=0
for FILE in $FILES
do
  let noTasks=noTasks+1
done

JOBSCRIPT=`tempfile`
chmod 0775 $JOBSCRIPT

echo '#!/bin/bash'        > $JOBSCRIPT
echo '#$ -t 1-'$noTasks  >> $JOBSCRIPT
echo '#$ -cwd          ' >> $JOBSCRIPT
echo '#$ -l qname=all.q' >> $JOBSCRIPT
echo '#$ -N '$ANIMFILE   >> $JOBSCRIPT
echo '#'                 >> $JOBSCRIPT
echo 'export TMPDIR=/tmp/'    >> $JOBSCRIPT

curTask=1  
for FILE in $FILES
do
  echo 'if [[ "$SGE_TASK_ID" -eq '$curTask' ]]'                                                          >> $JOBSCRIPT
  echo 'then'                                                                                            >> $JOBSCRIPT
  echo '   '$HDF52POVRAY -i $FILE.h5part -o '$TMPDIR'/$FILE.pov -s $SCALE -f $HEADER                     >> $JOBSCRIPT
#  echo '   'povray -d +H576 +W768 +Q$POVQUAL $POVOPTS +O'$TMPDIR'/$FILE.png '$TMPDIR'/$FILE.pov          >> $JOBSCRIPT
#  echo '   'convert '$TMPDIR'/$FILE.png $ANIMFILE-frames/$FILE.gif                                       >> $JOBSCRIPT
  echo '   'povray -d +H576 +W768 +Q$POVQUAL $POVOPTS +O$FILE.png '$TMPDIR'/$FILE.pov   >> $JOBSCRIPT
#  echo '   rm $TMPDIR/'$FILE.png' $TMPDIR/'$FILE.pov                                                     >> $JOBSCRIPT
  echo '   rm $TMPDIR/'$FILE.pov                                                                         >> $JOBSCRIPT
  echo 'fi'                                                                                              >> $JOBSCRIPT
  let curTask=curTask+1
done

qsub $JOBSCRIPT
rm $JOBSCRIPT
