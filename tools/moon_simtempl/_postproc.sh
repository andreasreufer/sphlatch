#!/bin/bash
if [ ! -e config.sh ]
then 
  echo "config.sh not found!"
  exit 1
fi
source config.sh

FILES=`find . -maxdepth 1 -iname 'dump*.h5part' | sort`

cat escDom0* | sort -n >escapees
echo '1.e99' >>escapees

for FILE in $FILES
do
  if [ ! -e gifs_full/$FILE.gif ]
  then
    touch gifs_full/$FILE.gif
    ./scatterplot.sh ${FILE}
    moon_pp -i ${FILE} -o postproc/${FILE}_pp.hdf5 -r $PPMAXRAD >>simSummary
  fi
done

gifsicle --delay 5 --colors 256 gifs_full/dump*.gif  >${SIMNAME}_scatter_full.gif
gifsicle --delay 5 --colors 256 gifs_zoom/dump*.gif  >${SIMNAME}_scatter_zoom.gif

