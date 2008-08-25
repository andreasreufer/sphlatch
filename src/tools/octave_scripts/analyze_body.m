#!/usr/bin/octave --persist
args = argv;

for i=1:size(args)(1)
  filename = args{i};
  dump = load(filename).current;

  com = sum( dump.pos(:,:) * dump.m' , 2) / sum(dump.m);
  dump.pos(1,:) -= com(1);
  dump.pos(2,:) -= com(2);
  dump.pos(3,:) -= com(3);

  dump = h5part_projvect( dump, 'vel');
  dump = h5part_projvect( dump, 'acc');

  subplot(2,2,1);
  title("rho");
  plot( dump.dist, dump.rho, strcat(num2str(i),".;",filename,";"));
  hold on;
  
  subplot(2,2,2);
  title("p");
  plot( dump.dist, dump.p, strcat(num2str(i),".;",filename,";"));
  hold on;
  
  subplot(2,2,3);
  title("acc");
  plot( dump.dist, dump.acc_proj, strcat(num2str(i),".;",filename,";"));
  hold on;
  
  subplot(2,2,4);
  title("vel");
  plot( dump.dist, dump.vel_proj, strcat(num2str(i),".;",filename,";"));
  hold on;
endfor
