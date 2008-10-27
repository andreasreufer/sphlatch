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
#vinf = 0.;
vinf = 1.e6;
Ltot = 1.;
Ltot *= lem;

# determine impact parameters

mtot = mtar + mimp;
gamma = mimp / mtot;

vgraz = sqrt( ( 2.*G*mtar / (Rtar + Rimp) ) + vinf^2 );
Lgraz = vgraz * mimp * (Rtar + Rimp);
bscal = Ltot / Lgraz;
impangl = asin( bscal );

vesc = sqrt(  2.*G*( mtar + mimp ) / (Rtar + Rimp) );
vimp = sqrt( vesc^2. + vinf^2. );
rimp = Rtar + Rimp;

printf("\n parameters summary:\n\n");
printf("gamma = %f\n", gamma);
printf("vgraz = %e cm/s\n", vgraz);
printf("Lgraz = %f lem\n", Lgraz/lem);
printf("bscal = %f  (impact angle = %f°)\n", bscal, impangl*360/(2*pi) );
printf("vimp  = %e cm/s\n", vimp);
printf("\n");

# determine impactor orbit parameters

vimpr = vimp*cos(impangl);
vimpt = vimp*sin(impangl);
mu = G*mtar;
p = (rimp*vimpt)^2. / mu;
e = sqrt( (vimpt/vimp - 1)^2.+(vimpr/vimp)^2. );

