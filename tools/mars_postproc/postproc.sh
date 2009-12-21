#!/bin/bash

if [ ! -e config.sh ]
then 
  echo "config.sh not found!"
  exit 1
fi
source config.sh

if [ "$SURFDIM" = "2.0m" ]
then
  ZOOM="4"
  PTSIZE="0.25"
  SMOLEN="1.7"
fi
if [ "$SURFDIM" = "4.0m" ]
then
  ZOOM="2"
  PTSIZE="0.50"
  SMOLEN="3.5"
fi
if [ "$SURFDIM" = "8.0m" ]
then
  ZOOM="1"
  PTSIZE="1.00"
  SMOLEN="7.0"
fi


writeOctScript()
{
cat >$1 << EOT
#!/usr/bin/octave
args = argv;
n = size(args)(1);

evToK = 11605.333;

for i=1:n
  filename = args{i};
  
  slice = load(strcat(filename,"_cut.h5part")).current;
  
  #isice = double( 1 - ( ( slice.mat - 2 ) / 2 ) );

  fid = fopen( strcat(filename,"_cut.txt"), 'w', 'native');
  noparts=size(slice.pos)(2);
  for j=1:noparts
    fprintf(fid, '%e\t%e\t%e\t%e\t%e\t%e\t%u\t%u\t%e\t%u\t%e\t%e\n',
            slice.pos(1,j), slice.pos(3,j),
            slice.vel(1,j), slice.vel(3,j),
            slice.u(j), slice.dam(j),
            slice.mat(j), slice.id(j),
            evToK*slice.T(j), slice.phase(j), slice.p(j), slice.peakp(j));
  endfor
  fclose(fid);
endfor
EOT
chmod 0755 $1
}

FILES=`find . -maxdepth 1 -iname 'out???.xdr' | sort`

for FILE in $FILES
do
    FILEPRFX=`echo $FILE|sed -e "s/\.xdr//"|sed -e "s/\.\///"`

    # convert to h5part
    if [ ! -e ${FILEPRFX}.h5part ] && [ ! -e ${FILEPRFX}.lock ]
    then
      touch ${FILEPRFX}.lock
      echo "${FILEPRFX}.xdr -> ${FILEPRFX}.h5part"
      sleep 30
      cdat2h5part -i ${FILEPRFX}.xdr -o ${FILEPRFX}.h5part
      h5part_add_ANEOSout_table -i ${FILEPRFX}.h5part
      #h5part_get_phase_histogram -i ${FILEPRFX}.h5part -o ice_phases.txt -n 2
      #h5part_get_phase_histogram -i ${FILEPRFX}.h5part -o dun_phases.txt -n 4
      rm ${FILEPRFX}.lock
    fi
    
    # convert to txt
    if [ ! -e ${FILEPRFX}_cut.txt ] && [ -e ${FILEPRFX}.h5part ] && [ ! -e ${FILEPRFX}.lock ]
    then
      touch ${FILEPRFX}.lock
      echo "${FILEPRFX}.h5part -> ${FILEPRFX}_cut.txt"
      OCTSCRIPT=`mktemp ./octXXXXXX`
      writeOctScript $OCTSCRIPT
      h5part_slice --ymin=-${SMOLEN} --ymax=${SMOLEN} -i ${FILEPRFX}.h5part -o ${FILEPRFX}_cut.h5part
      $OCTSCRIPT ${FILEPRFX}
      rm ${FILEPRFX}_cut.h5part $OCTSCRIPT
      rm ${FILEPRFX}.lock
    fi

    # make plot
    if [ ! -e gifs/${FILEPRFX}.gif ] && [ -e ${FILEPRFX}_cut.txt ] && [ ! -e ${FILEPRFX}.lock ]
    then
      touch ${FILEPRFX}.lock
      echo "plot ${FILEPRFX}.gif"
      TIMESTR=`h5part_readattr -i $FILEPRFX.h5part -k TIME|cut -c 7-17`
      ./scatterplot.py ${FILEPRFX}_cut.txt $TIMESTR $ZOOM $PTSIZE
      echo "${FILEPRFX}_cut.txt.png -> gifs/$FILEPRFX.gif"
      convert ${FILEPRFX}_cut.txt.png gifs/$FILEPRFX.gif
      rm ${FILEPRFX}_cut.txt.png
      rm ${FILEPRFX}.lock
    fi
done

gifsicle --delay 5 --colors 256 gifs/out*.gif  >${SIMNAME}.gif

echo "finishing"
exit 0

