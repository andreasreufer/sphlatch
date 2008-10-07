
DELAY=5
POVQUAL=2
POVOPTS="+A +AM2"
SCALE=0.05

ZOOMHEADER=zoom_header.pov
NORMHEADER=default_header.pov

FILES=`find . -maxdepth 1 -iname '*.h5part'|sed 's/\.h5part//'|sed 's/\.\///'|sort `

mkdir gifs 2>&1 >/dev/null

for FILE in $FILES
do
  echo hdf52povray -i $FILE.h5part -o $FILE-zoom.pov -s $SCALE -f $ZOOMHEADER
  hdf52povray -i $FILE.h5part -o $FILE-zoom.pov -s $SCALE -f $ZOOMHEADER
  echo povray -d +H576 +W768 +Q$POVQUAL $POVOPTS $FILE-zoom.pov
  povray -d +H576 +W768 +Q$POVQUAL $POVOPTS $FILE-zoom.pov
  convert $FILE-zoom.png $FILE-zoom.gif
  rm $FILE-zoom.png $FILE-zoom.pov
  mv $FILE-zoom.gif gifs
  
  echo hdf52povray -i $FILE.h5part -o $FILE-zoom.pov -s $SCALE
  hdf52povray -i $FILE.h5part -o $FILE-norm.pov -s $SCALE
  echo povray -d +H576 +W768 +Q$POVQUAL $POVOPTS $FILE-norm.pov
  povray -d +H576 +W768 +Q$POVQUAL $POVOPTS $FILE-norm.pov
  convert $FILE-norm.png $FILE-norm.gif
  rm $FILE-norm.png $FILE-norm.pov
  mv $FILE-norm.gif gifs
done

gifsicle --colors 256 --delay $DELAY -i gifs/*-zoom.gif > anim-zoom.gif
gifsicle --colors 256 --delay $DELAY -i gifs/*-norm.gif > anim-norm.gif



