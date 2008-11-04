#!/usr/bin/octave --persist

setup_gi

parts.pos(1,1) = rimp0(1);
parts.pos(2,1) = rimp0(2);
parts.vel(1,1) = vimp0(1);
parts.vel(2,1) = vimp0(2);
parts.m(1) = mimp;

parts.pos(1,2) = rtar0(1);
parts.pos(2,2) = rtar0(2);
parts.vel(1,2) = vtar0(1);
parts.vel(2,2) = vtar0(2);
parts.m(2) = mtar;

parts.acc = zeros(2,2);

parts_orig = parts;

tabs = timp - t0;
ttot = 1.10*abs(tabs);


nsteps = 5000;
dt = ttot / nsteps;
pos = zeros(nsteps, 4);
dist = zeros(nsteps, 1);
vrel = zeros(nsteps, 1);
t    = zeros(nsteps, 1);

hit = false;
hitstep = nsteps;

for step = 1:nsteps
  parts = bf_nbody(parts, G, dt);

  pos(step,:) = [ parts.pos(:,1); parts.pos(:,2) ];

  rvect = parts.pos(:,1) - parts.pos(:,2);
  dist(step) = sqrt( sum( rvect .* rvect ) );
  vvect = parts.vel(:,1) - parts.vel(:,2);
  vrel(step) = sqrt( sum( vvect .* vvect ) );
  
  if ( !hit && ( dist(step) < rimp ) )
    hitstep = step;
    hit = true;
  end

  tabs += dt;
  t(step) = tabs;
end

figure(1)
plot( pos(:,1), pos(:,2), '1.');
hold on
plot( pos(nsteps,1), pos(nsteps,2), '1o');
plot( pos(hitstep,1), pos(hitstep,2), '1x');
title( 'trajectories [cm]' );
plot( pos(:,3), pos(:,4), '3.');
plot( pos(nsteps,3), pos(nsteps,4), '3o');
plot( pos(hitstep,3), pos(hitstep,4), '3x');
axis('square');
hold off

figure(2)
plot( t, dist, '1.' );
hold on
plot( t(hitstep), dist(hitstep), '1x');
title( 'distance [cm]' );
hold off

figure(3)
plot( t, vrel, '1.' );
hold on
plot( t(hitstep), vrel(hitstep), '1x');
title( 'relative velocity [cm/s]' );
hold off
