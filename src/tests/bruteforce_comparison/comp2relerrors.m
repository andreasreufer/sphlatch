filename = input("inputfile: ");

comp = load(filename);

a = sqrt( comp(:,8) .* comp(:,8) + comp(:,10) .* comp(:,10) + comp(:,12) .* comp(:,12));
da = sqrt( comp(:,9) .* comp(:,9) + comp(:,11) .* comp(:,11) + comp(:,13) .* comp(:,13));

x = ( 2.1544.^[0:16] ) / 2.1539e+05;
rel = da ./ a;
id = comp(:,1);

[ counts, x ] = hist( rel, x );

output = [ counts' , x' ]

save -text "histogram.txt" output;
