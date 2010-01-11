#!/usr/bin/env octave

output_max_field_width(15);
output_precision(10);

arg_list = argv();

dump = load(arg_list{1}).current;
idlim = str2double( arg_list{2} );

tar = h5part_filter(dump, dump.id, idlim);
tar_mass = sum(tar.m);

tar_pos = zeros(3,1);
tar_pos(1) = sum( tar.pos(1,:) .* tar.m ) / tar_mass;
tar_pos(2) = sum( tar.pos(2,:) .* tar.m ) / tar_mass;
tar_pos(3) = sum( tar.pos(3,:) .* tar.m ) / tar_mass;

tar_vel = zeros(3,1);
tar_vel(1) = sum( tar.vel(1,:) .* tar.m ) / tar_mass;
tar_vel(2) = sum( tar.vel(2,:) .* tar.m ) / tar_mass;
tar_vel(3) = sum( tar.vel(3,:) .* tar.m ) / tar_mass;

tar_acc = zeros(3,1);
tar_acc(1) = sum( tar.acc(1,:) .* tar.m ) / tar_mass;
tar_acc(2) = sum( tar.acc(2,:) .* tar.m ) / tar_mass;
tar_acc(3) = sum( tar.acc(3,:) .* tar.m ) / tar_mass;
[tar_pos, tar_vel, tar_acc]


imp = h5part_filter(dump, idlim, dump.id);
imp_mass = sum(imp.m);

imp_pos = zeros(3,1);
imp_pos(1) = sum( imp.pos(1,:) .* imp.m ) / imp_mass;
imp_pos(2) = sum( imp.pos(2,:) .* imp.m ) / imp_mass;
imp_pos(3) = sum( imp.pos(3,:) .* imp.m ) / imp_mass;

imp_vel = zeros(3,1);
imp_vel(1) = sum( imp.vel(1,:) .* imp.m ) / imp_mass;
imp_vel(2) = sum( imp.vel(2,:) .* imp.m ) / imp_mass;
imp_vel(3) = sum( imp.vel(3,:) .* imp.m ) / imp_mass;

imp_acc = zeros(3,1);
imp_acc(1) = sum( imp.acc(1,:) .* imp.m ) / imp_mass;
imp_acc(2) = sum( imp.acc(2,:) .* imp.m ) / imp_mass;
imp_acc(3) = sum( imp.acc(3,:) .* imp.m ) / imp_mass;
 
[imp_pos, imp_vel, imp_acc]

fid = fopen( 'target', 'a', 'native');
fprintf(fid, '%s\t', arg_list{1} )
fprintf(fid, '%e\t', tar_pos, tar_vel, tar_acc)
fprintf(fid, '\n' )
fclose(fid);

fid = fopen( 'impactor', 'a', 'native');
fprintf(fid, '%s\t', arg_list{1} )
fprintf(fid, '%e\t', imp_pos, imp_vel, imp_acc)
fprintf(fid, '\n' )
fclose(fid);

