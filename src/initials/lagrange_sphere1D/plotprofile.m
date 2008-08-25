
if (! exist("filename") )
  filename="profile.hdf5";
endif

if (exist("color_idx") )
  color_idx++;
else
  color_idx = 1;
endif

color = strcat(num2str(color_idx),' ');
prof = load(filename);

subplot(2,2,1);
plot( prof.r, prof.rho, color);
title('density');
hold on;

subplot(2,2,2);
plot( prof.r, prof.u, color);
title('spec. energy');
hold on;

subplot(2,2,3);
semilogy( prof.r, prof.p, color);
title('pressure');
hold on;

subplot(2,2,4);
plot( prof.r, prof.mat, strcat('x',color));
title('material');
hold on;
