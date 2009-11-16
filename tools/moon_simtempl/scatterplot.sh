#!/bin/bash

writePythonScript()
{
cat >$1 << EOT
#!/usr/bin/python

from pylab import *
filename = sys.argv[1]
timesec = float(filename[15:26])
timestring = '$\\mathrm{ t = ' + '%6.2f' % ( timesec/3600 ) + ' h }$'

rc('grid', color='grey')
rc('text', usetex=True)

slice = load(filename)

axes(axisbg='black')

# slice[:,2] is 0 for mantle and 1 for iron
# slice[:,3] is 0 for protoearth and 1 for impactor
color = 1-(((slice[:,2]-0.5)/(slice[:,3]+1))+1)

scatter( slice[:,0], slice[:,1], 1, color, linewidth=0)
grid(True)
xlabel('$\\mathrm{x [cm]}$')
ylabel('$\\mathrm{y [cm]}$')
title(timestring)

axis("scaled")
axis($FULLSCALE)
savefig(filename + "_full.png")

axis("scaled")
axis($ZOOMSCALE)
savefig(filename + "_zoom.png")
EOT
chmod 0755 $1
}

writeConverterScript()
{
cat >$1 << EOT
#!/usr/bin/octave
args = argv;
n = size(args)(1);

for i=1:n
  filename = args{i};
  
  dump = load(filename).current;
  slice = h5part_slice( dump, [-inf,-0.5e8,-inf], [inf,0.5e8,inf] );
  fid = fopen( strcat(filename, "_slice.txt"), 'w', 'native');
  noparts=size(slice.pos)(2);
  for j=1:noparts
    fprintf(fid, '%e\t%e\t%e\t%e\n',
            slice.pos(3,j), slice.pos(1,j),
            (slice.mat(j)-4)', ((slice.id(j) / 2e6)) );
  endfor
  fclose(fid);
endfor
EOT
chmod 0755 $1
}

if [ ! -e config.sh ]
then
  echo "config.sh not found!"
  exit 1
fi
source config.sh

FILE=$1

OCTSCRIPT=`tempfile`
writeConverterScript $OCTSCRIPT
$OCTSCRIPT ${FILE}

IPYSCRIPT=`tempfile`
writePythonScript $IPYSCRIPT
$IPYSCRIPT ${FILE}_slice.txt
rm ${FILE}_slice.txt

convert ${FILE}_slice.txt_full.png gifs_full/$FILE.gif
rm ${FILE}_slice.txt_full.png

convert ${FILE}_slice.txt_zoom.png gifs_zoom/$FILE.gif
rm ${FILE}_slice.txt_zoom.png

rm $IPYSCRIPT
rm $OCTSCRIPT

echo "finishing"
exit 0

