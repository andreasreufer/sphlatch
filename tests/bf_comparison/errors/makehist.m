x = [ - ( 2.1544.^[-(-15:0)] ) / 9.9976e4 , 0 , ( 2.1544.^[0:15] ) / 9.9976e4 ]
theta = load("random_monopol_theta_1.0");
[theta_1, xx] = hist( theta, x )
theta = load("random_monopol_theta_0.9");
[theta_09, xx] = hist( theta, x )
theta = load("random_monopol_theta_0.8");
[theta_08, xx] = hist( theta, x );
theta = load("random_monopol_theta_0.7");
[theta_07, xx] = hist( theta, x );
theta = load("random_monopol_theta_0.6");
[theta_06, xx] = hist( theta, x );
theta = load("random_monopol_theta_0.5");
[theta_05, xx] = hist( theta, x );
theta = load("random_monopol_theta_0.4");
[theta_04, xx] = hist( theta, x );
theta = load("random_monopol_theta_0.3");
[theta_03, xx] = hist( theta, x );
theta = load("random_monopol_theta_0.2");
[theta_02, xx] = hist( theta, x );
theta = load("random_monopol_theta_0.1");
[theta_01, xx] = hist( theta, x );
theta = load("random_monopol_theta_0.05");
[theta_005, xx] = hist( theta, x );
theta = load("random_monopol_theta_0.01");
[theta_001, xx] = hist( theta, x );

output = [x', theta_1', theta_09', theta_08',theta_07',theta_06',theta_05',theta_04',theta_03',theta_02',theta_01',theta_005',theta_001']

save -ascii "histograms" output
