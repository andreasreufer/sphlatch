
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
hold on;
plot( prof.r, prof.rho, color);
title('density');

subplot(2,2,2);
hold on;
plot( prof.r, prof.u, color);
title('spec. energy');

subplot(2,2,3);
hold on;
semilogy( prof.r, prof.p, color);
title('pressure');

subplot(2,2,4);
hold on;
plot( prof.r, 11604.*prof.T, color);
title('material');

