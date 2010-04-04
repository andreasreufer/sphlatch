#!/usr/bin/env python

import os
import sys
import shutil
import tempfile

def isimg(str):
  return os.path.splitext(str)[1] == sfx

dir = sys.argv[1]
sfx = sys.argv[2]
out = sys.argv[3]
#frg = sys.argv[4]
frg = '-b 20000k'

filelist = filter( isimg, os.listdir( dir ) )

tmpdir = tempfile.mkdtemp( dir=dir )

count = 0
for file in filelist:
  print "../" + file, "  -->  ", tmpdir + "/img" + '%06d' % count + sfx
  os.symlink( "../" + file, tmpdir + "/img" + '%06d' % count + sfx )
  count += 1

print "ffmpeg -i " + tmpdir + "/img%06d" + sfx + " " + frg + " " + out
os.system( "ffmpeg -i " + tmpdir + "/img%06d" + sfx + " " + frg + " " + out )

count = 0
for file in filelist:
  os.remove(tmpdir + "/img" + '%06d' % count + sfx)
  count += 1
os.removedirs(tmpdir)


