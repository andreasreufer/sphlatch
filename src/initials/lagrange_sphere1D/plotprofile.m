
color = strcat(num2str(color_idx),' ');
prof = load("profile.hdf5");

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
