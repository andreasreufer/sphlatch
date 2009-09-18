#!/bin/bash

if [ ! -e config.sh ]
then 
  echo "config.sh not found!"
  exit 1
fi
source config.sh


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
    if [ ! -e ${FILEPRFX}.h5part ]
    then
      touch ${FILEPRFX}.h5part
      sleep 30
      cdat2h5part -i ${FILEPRFX}.xdr -o ${FILEPRFX}.h5part
      h5part_add_ANEOSout_table -i ${FILEPRFX}.h5part
      #h5part_get_phase_histogram -i ${FILEPRFX}.h5part -o ice_phases.txt -n 2
      #h5part_get_phase_histogram -i ${FILEPRFX}.h5part -o dun_phases.txt -n 4
    fi
    
    # convert to txt
    if [ ! -e ${FILEPRFX}_cut.txt ]
    then
      touch ${FILEPRFX}_cut.txt
      OCTSCRIPT=`mktemp ./octXXXXXX`
      writeConverterScript $OCTSCRIPT
      h5part_slice --ymin=-3.5 --ymax=3.5 -i ${FILEPRFX}.h5part -o ${FILEPRFX}_cut.h5part
      $OCTSCRIPT ${FILEPRFX}
      rm ${FILEPRFX}_cut.h5part $OCTSCRIPT
    fi

    # make plot
    if [ ! -e gifs/${FILEPRFX}.gif ]
    then
      touch gifs/${FILEPRFX}.gif
      TIMESTR=`h5part_readattr -i $FILEPRFX.h5part -k TIME|cut -c 7-17`
      ./scatterplot.py ${FILEPRFX}_cut.txt $TIMESTR
      convert ${FILEPRFX}_cut.txt.png gifs/$FILEPRFX.gif
      rm ${FILEPRFX}_cut.txt.png
    fi
done

gifsicle --delay 5 --colors 256 gifs/out*.gif  >scatter.gif

echo "finishing"
exit 0

