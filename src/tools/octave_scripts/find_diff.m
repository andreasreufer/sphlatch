#!/usr/bin/octave -q

arg_list = argv();
exitStat = 0;

diff = h5part_compare( load( arg_list{1}).current, load( arg_list{2}).current );

for [diffVal, diffKey] = diff

  maxDev = max(max( diffVal ) );
  if ( maxDev > 0. && strncmp("d_", diffKey, 2) )
    printf("non-zero difference: %s (%f)\n", diffKey, maxDev);
    exitStat = 1;
  endif

endfor

exit(exitStat);
