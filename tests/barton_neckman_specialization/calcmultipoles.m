cells = load("cells");

noCells = size(cells)(1)
parent = cells(1,:);

cx = cells(2:noCells, 7);
cy = cells(2:noCells, 8);
cz = cells(2:noCells, 9);
cm = cells(2:noCells, 10);
q11 = cells(2:noCells, 11);
q22 = cells(2:noCells, 12);
q33 = cells(2:noCells, 13);
q12 = cells(2:noCells, 14);
q13 = cells(2:noCells, 15);
q23 = cells(2:noCells, 16);
s11 = cells(2:noCells, 17);
s22 = cells(2:noCells, 18);
s33 = cells(2:noCells, 19);
s12 = cells(2:noCells, 20);
s21 = cells(2:noCells, 21);
s13 = cells(2:noCells, 22);
s31 = cells(2:noCells, 23);
s23 = cells(2:noCells, 24);
s32 = cells(2:noCells, 25);
s123 = cells(2:noCells, 26);

cmP = sum( cm );

cxP = sum( cx.*cm ) / cmP;
cyP = sum( cy.*cm ) / cmP;
czP = sum( cz.*cm ) / cmP;

rx = cx - cxP;
ry = cy - cyP;
rz = cz - czP;

rr = rx.*rx + ry.*ry + rz.*rz;

q11P = sum( cm.*( 3*rx.*rx - rr ) + q11 );
q22P = sum( cm.*( 3*ry.*ry - rr ) + q22 );
q33P = sum( cm.*( 3*rz.*rz - rr ) + q33 );

q12P = sum( cm.*( 3*rx.*ry ) + q12 );
q13P = sum( cm.*( 3*rx.*rz ) + q13 );
q23P = sum( cm.*( 3*ry.*rz ) + q23 );


s11P = sum( cm.*( 5*rx.*rx - 3*rr).*rx + 2.5*rx.*q11
            - rx.*q11 - ry.*q12 - rz.*q13 + s11 );

s22P = sum( cm.*( 5*ry.*ry - 3*rr).*ry + 2.5*ry.*q22
            - rx.*q12 - ry.*q22 - rz.*q23 + s22 );

s33P = sum( cm.*( 5*rz.*rz - 3*rr).*rz + 2.5*rz.*q33
            - rx.*q13 - ry.*q23 - rz.*q33 + s33 );

s12P = sum( cm.*( 15*rx.*rx - 3*rr).*ry
            + 5.0*rx.*q12 + 2.5*ry.*q11
            - rx.*q12 - ry.*q22 - rz.*q23 + s12 );
s21P = sum( cm.*( 15*ry.*ry - 3*rr).*rx
            + 5.0*ry.*q12 + 2.5*rx.*q22
            - rx.*q11 - ry.*q12 - rz.*q13 + s21 );

s13P = sum( cm.*( 15*rx.*rx - 3*rr).*rz
            + 5.0*rx.*q13 + 2.5*rz.*q11
            - rx.*q13 - ry.*q23 - rz.*q33 + s13 );
s31P = sum( cm.*( 15*rz.*rz - 3*rr).*rx
            + 5.0*rz.*q13 + 2.5*rx.*q33
            - rx.*q11 - ry.*q12 - rz.*q13 + s31 );

s23P = sum( cm.*( 15*ry.*ry - 3*rr).*rz
            + 5.0*ry.*q23 + 2.5*rz.*q22
            - rx.*q13 - ry.*q23 - rz.*q33 + s23 );
s32P = sum( cm.*( 15*rz.*rz - 3*rr).*ry
            + 5.0*rz.*q23 + 2.5*ry.*q33
            - rx.*q12 - ry.*q22 - rz.*q23 + s32 );

s123P = sum( 15.0*cm.*rx.*ry.*rz + 25.0*(rx.*q23 + ry.*q13 + rz.*q12) + 15.0*s123 );

input = parent(7:26)';
calc = [cxP, cyP, czP, cmP, q11P, q22P, q33P, q12P, q13P, q23P, s11P, s22P, s33P, s12P, s21P, s13P, s31P, s23P, s32P, s123P]';

[input, calc, input - calc]
