#!/bin/bash

TIMESTR="${2}h"
FILEIN="impact_9_1_povray_back.gif.${1}"
FILEOUT="out${1}.gif"

convert \
  -pointsize 36 -annotate +30+50 "Giant impact\n\nt = $TIMESTR" \
  -pointsize 36 -annotate +750+50 "mass ratio 9:1\ntotal mass 1Me\nL = 1lem"\
  -pointsize 18 -annotate +30+670 '"cold" run, Tillotson EOS,\nSPHLATCH code (SPH), 2M particles' \
  -pointsize 18 -annotate +750+670 '\nAndreas Reufer, University of Bern' \
  -fill white $FILEIN $FILEOUT


