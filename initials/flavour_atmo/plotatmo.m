#!/usr/bin/env octave

dump = load("out.h5part").current;
dump = h5part_projvect( dump, 'pos');

prof = load("profcgs.hdf5");

#figure(1)
#loglog( prof.r, prof.rho, '3-')
#hold on
#loglog( dump.dist, dump.rho, '1.')

figure(2)
loglog( prof.r, prof.p, '3-')
hold on
loglog( dump.dist, dump.p, '1.')
axis([1.0000e+08   1.0000e+11,   1.0000e+02  , 1.0000e+14])

#figure(3)
#semilogx( dump.dist, dump.m, 'x')
#hold on

#figure(4)
#semilogx( dump.dist, dump.h, 'x')
#hold on

figure(1)

