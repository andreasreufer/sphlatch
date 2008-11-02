G   = 6.6742e-8;
Re  = 6.37814e8;
me  = 5.9742e27;
lem = 3.4000e41;

#mtar =  input("m_target   [me] = ");
#mimp =  input("m_impactor [me] = ");
mtar = 0.90;
mimp = 0.10;
mtar *= me;
mimp *= me;

#Rtar =  input("r_target   [Re] = ");
#Rimp =  input("r_impactor [Re] = ");
Rtar = 0.9;
Rimp = 0.6;
Rtar *= Re;
Rimp *= Re;

#vinf =  input("v_inf    [cm/s] = ");
#Ltot  = input("L_tot     [lem] = ");
vinf = 0.0000e4;
Ltot = 1.520000;
Ltot *= lem;

# determine impact parameters
mtot = mtar + mimp;
gamma = mimp / mtot;

# some constants
vesc = sqrt(  2.*G*( mtar + mimp ) / (Rtar + Rimp) );
vimp = sqrt( vesc^2. + vinf^2. );
rimp = Rtar + Rimp;

mu = G*mtot;
k1 = rimp * (vimp^2. ) / mu;

vgraz = sqrt( ( 2.*mu / (Rtar + Rimp) ) + vinf^2 );
Lgraz = vgraz * mimp * (Rtar + Rimp);
bscal = Ltot / Lgraz;
betaimp = acos( bscal );


printf("\n parameters:\n");
printf("gamma = %f\n", gamma);
printf("vgraz = %e cm/s\n", vgraz);
printf("Lgraz = %f lem\n", Lgraz/lem);
printf("bscal = %f \n", bscal);
printf("vesc  = %e cm/s\n", vesc);

# total energy
Etot = ( (vimp^2.) / 2. ) - ( mu / rimp );

# determine impactor orbit eccentricity and true anomaly of impact
e = sqrt( ((k1 - 1)*cos(betaimp))^2. + (sin(betaimp))^2. );

thetaimp = atan2(( k1 * sin(betaimp) * cos(betaimp) ),
( k1 * (cos(betaimp) ^2. ) - 1.) );

printf("e     = %f\n", e);

# determine perigee radius and velocity
rperih = rimp *( 1. + e*cos(thetaimp) ) / ( 1. + e );
vperih = sqrt( vinf^2. + ( 2.*mu / rperih ) );
printf("\n perihelon:\n");
printf("r     = %e cm/s\n", rperih );
printf("v     = %e cm/s\n", vperih);
printf("theta = %f° (true anomaly)\n", 0.);

printf("\n impact:\n");
printf("r     = %e cm\n", rimp);
printf("v     = %e cm/s\n", vimp);
printf("theta = %f° (true anomaly)\n", thetaimp*360/(2*pi));
printf("beta  = %f° (impact angle)\n", betaimp*360/(2*pi));

# determine true anomaly theta0 for r0 = x*rimp
r0 = 5.*rimp;
v0 = sqrt( vinf^2. + ( 2.*mu / r0 ) );
theta0 = acos( ( rperih*(1. + e ) - r0 ) / ( e * r0 ) );
beta0 = acos( Ltot / ( r0*v0*mimp ) );

printf("\n initial setup:\n");
printf("( r0 = %f rimp )\n", r0 / rimp);
printf("r     = %e cm\n", r0);
printf("v     = %e cm/s\n", v0);
printf("theta = %f° (true anomaly)\n", theta0*360/(2*pi));
printf("beta  = %f° \n", beta0*360/(2*pi));

# initial velocity vector
r0vect(1) = -r0* cos( (pi/2.) - thetaimp + theta0 );
r0vect(2) =  r0* sin( (pi/2.) - thetaimp + theta0 );

alpha = theta0 - thetaimp - beta0;
v0vect(1) = -v0 * cos(alpha);
v0vect(2) = -v0 * sin(alpha);

# now place everything, so that center of mass is (0,0)
# and there is no resulting momentum

rimp0(1) =  ( mtar / mtot )*r0vect(1);
rimp0(2) =  ( mtar / mtot )*r0vect(2);
vimp0(1) =  ( mtar / mtot )*v0vect(1);
vimp0(2) = -( mtar / mtot )*v0vect(2);

rtar0(1) = -( mimp / mtot )*r0vect(1);
rtar0(2) = -( mimp / mtot )*r0vect(2);
vtar0(1) = -( mimp / mtot )*v0vect(1);
vtar0(2) =  ( mimp / mtot )*v0vect(2);


# determine time to perihelon for impact and for initial setup

if e > 1.
  a = rperih / ( e - 1. );
  k2 = sqrt( mu / (a^3. ) );
  
  timp = ( (e*sqrt(e^2.-1.)*sin(thetaimp))/(1. + e*cos(thetaimp) ) -
    log( (sqrt(e^2.-1.) + (e - 1.)*tan( thetaimp / 2. ))/
         (sqrt(e^2.-1.) - (e - 1.)*tan( thetaimp / 2. )) ) ) / k2;
  
  t0 = ( (e*sqrt(e^2.-1.)*sin(theta0))/(1. + e*cos(theta0) ) -
    log( (sqrt(e^2.-1.) + (e - 1.)*tan( theta0 / 2. ))/
         (sqrt(e^2.-1.) - (e - 1.)*tan( theta0 / 2. )) ) ) / k2;
else
  k2 = sqrt( mu / ( 8. * rperih^3.) );

  timp = abs(0.5*(tan(thetaimp / 2.) + (1./3.)*((tan(thetaimp / 2.))^3.)) / k2);
  t0   = abs(0.5*(tan(theta0   / 2.) + (1./3.)*((tan(theta0   / 2.))^3.)) / k2);
end

printf("t     = %f s\n", timp - t0);
printf("\n");

printf(" commands to setup the initial file:\n");
printf("h5part_displace -i impactor.h5part --pos [%e,%e,%e] --vel [%e,%e,%e] --id 2000000 \n", rimp0(1), 0., rimp0(2), vimp0(1), 0., vimp0(2));
printf("h5part_displace -i target.h5part --pos [%e,%e,%e] --vel [%e,%e,%e]\n", rtar0(1), 0., rtar0(2), vtar0(1), 0., vtar0(2));
printf("h5part_combine -a target.h5part -b impactor.h5part -o combined.h5part\n");
printf("h5part_writeattr -i combined.h5part -k time -v=%e\n", timp - t0 );
printf("\n");

