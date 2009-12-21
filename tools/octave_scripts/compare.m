#!/usr/bin/env octave

arg_list = argv();

dumpA = load(arg_list{1}).current;
dumpB = load(arg_list{2}).current;
dumpC = load(arg_list{3}).current;

dumpA.acc_scal = sqrt( dot( dumpA.acc(:,:), dumpA.acc(:,:) ) );
dumpB.acc_scal = sqrt( dot( dumpB.acc(:,:), dumpB.acc(:,:) ) );
dumpC.acc_scal = sqrt( dot( dumpB.acc(:,:), dumpC.acc(:,:) ) );

acc_relB = ( dumpB.acc_scal - dumpA.acc_scal ) ./ dumpA.acc_scal;
acc_relC = ( dumpC.acc_scal - dumpA.acc_scal ) ./ dumpA.acc_scal;

x = 10.^[-9:2];
[nB,x] = hist( abs(acc_relB), x);
[nC,x] = hist( abs(acc_relC), x);
[x',nB',nC']

[xs, ysB] = stairs( x,nB );
[xs, ysC] = stairs( x,nC );
loglog(xs, ysB+1, '1-');
hold on
loglog(xs, ysC+1, '2-');

print( arg_list{4} );

