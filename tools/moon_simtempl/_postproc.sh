#!/bin/bash
if [ ! -e config.sh ]
then 
  echo "config.sh not found!"
  exit 1
fi
source config.sh

wPartConv()
{
cat >$1 << EOT2
#!/usr/bin/octave
args = argv;
n = size(args)(1);

for i=1:n
  filename = args{i};
  
  slice = load(filename).current;
  fid = fopen( strcat(filename, "_slice.txt"), 'w', 'native');
  noparts=size(slice.pos)(2);
  for j=1:noparts
    fprintf(fid, '%e\t%e\t%e\t%e\n',
            slice.pos(3,j), slice.pos(1,j),
            slice.mat(j), ((slice.id(j) / 2e6)) );
  endfor
  fclose(fid);
endfor
EOT2
chmod 0755 $1
}

wClumpConv()
{
cat >$1 << EOT3
#!/usr/bin/octave
args = argv;
n = size(args)(1);

for i=1:n
  filename = args{i};
  
  slice = load(filename).current;
  fid = fopen( strcat(filename, ".txt"), 'w', 'native');
  noparts=size(slice.pos)(2);
  for j=1:noparts
    if ( slice.m(j) > 1.e23 )
      fprintf(fid, '%e\t%e\t%e\t%e\n',
              slice.pos(3,j), slice.pos(1,j),
              slice.rc(j), slice.id(j) );
    end
  endfor
  fclose(fid);
endfor
EOT3
chmod 0755 $1
}


FILES=`find . -maxdepth 1 -iname 'dump*.h5part' | sort`

#cat escDom0* | sort -n >escapees
#echo '1.e99' >>escapees

for FILE in $FILES
do
  FILEPRFX=`echo $FILE| sed -e's/\.h5part//g'`
  if [ ! -e gifs_full/${FILEPRFX}.gif ] && [ ! -e ${FILEPRFX}.lock ]
  then
    touch ${FILEPRFX}.lock
    OCTSCRIPT=`mktemp ./octXXXXXX`
    SLICETMP=`mktemp  ./slcXXXXXX`

    TIME=`h5part_readattr -k time -i ${FILEPRFX}.h5part |cut -c 6-14`
    h5part_slice --ymin=-0.5e8 --ymax=0.5e8 -i ${FILEPRFX}.h5part -o $SLICETMP
    wPartConv $OCTSCRIPT
    $OCTSCRIPT $SLICETMP
    
    wClumpConv $OCTSCRIPT
    $OCTSCRIPT ${FILEPRFX}_clumps.h5part

    rm $OCTSCRIPT $SLICETMP

    ./scatterplot.py ${SLICETMP}_slice.txt ${FILEPRFX}_clumps.h5part.txt $TIME
    #rm ${SLICETMP}_slice.txt ${FILEPRFX}_clumps.h5part.txt
    rm ${SLICETMP}_slice.txt

    #convert ${SLICETMP}_slice.txt_full.png gifs_full/$FILEPRFX.gif
    #convert ${SLICETMP}_slice.txt_zoom.png gifs_zoom/$FILEPRFX.gif
    #rm ${SLICETMP}_slice.txt_????.png

    #moon_pp -i ${FILE} -o postproc/${FILE}_pp.hdf5 -r $PPMAXRAD >>simSummary
    rm ${FILEPRFX}.lock
  fi
done

#gifsicle --delay 5 --colors 256 gifs_full/dump*.gif  >${SIMNAME}_scatter_full.gif
#gifsicle --delay 5 --colors 256 gifs_zoom/dump*.gif  >${SIMNAME}_scatter_zoom.gif

