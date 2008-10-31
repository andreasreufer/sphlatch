G   = 6.6742e-8;
Re  = 6.37814e8;
me  = 5.9742e27;
lem = 3.4000e41;

#mtar =  input("m_target   [me] = ");
#mimp =  input("m_impactor [me] = ");
mtar = 0.9;
mimp = 0.1;
mtar *= me;
mimp *= me;

#Rtar =  input("r_target   [Re] = ");
#Rimp =  input("r_impactor [Re] = ");
Rtar = 0.8;
Rimp = 0.4;
Rtar *= Re;
Rimp *= Re;

#vinf =  input("v_inf    [cm/s] = ");
#Ltot  = input("L_tot     [lem] = ");
vinf = 0.e6;
Ltot = 1.00000;
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


printf("\n parameters summary:\n\n");
printf("gamma = %f\n", gamma);
printf("vgraz = %e cm/s\n", vgraz);
printf("Lgraz = %f lem\n", Lgraz/lem);
printf("bscal = %f  (impact angle = %f°)\n", bscal, betaimp*360/(2*pi) );
printf("vimp  = %e cm/s\n", vimp);
printf("vesc  = %e cm/s\n", vesc);

# total energy
Etot = ( (vimp^2.) / 2. ) - ( mu / rimp );


# determine impactor orbit eccentricity and true anomaly of impact
e = sqrt( ((k1 - 1)*cos(betaimp))^2. + (sin(betaimp))^2. );
thetaimp = atan(( k1 * sin(betaimp) * cos(betaimp) ) 
  / ( k1 * (cos(betaimp) ^2. ) - 1.) );

printf("e     = %f\n", e);
printf("theta = %f° (impact true anomaly)\n", thetaimp*360/(2*pi));
printf("\n");

# determine 


#

