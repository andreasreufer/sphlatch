iter = load("out_iterate.h5part").current;
tabl = load("out_table.h5part").current;

dp = ( tabl.p - iter.p ) ./ iter.p;
dcs = ( tabl.cs - iter.cs ) ./ iter.cs;

subplot(1,2,1);
plot( tabl.p, iter.p, '1.');

subplot(1,2,2);
plot( tabl.cs, iter.cs, '1.');

