function [defgrad] = contAxi3ts(ec, ed)

r1 = ec(1,1);
r2 = ec(1,2);
r3 = ec(1,3);

z1 = ec(2,1);
z2 = ec(2,2);
z3 = ec(2,3);

z0 = (z1+z2+z3)/3;
r0 = (r1+r2+r3)/3;


A = [1,  1,  1;
    r1, r2, r3;
    z1, z2, z3];

A0e = 1/2*det(A);

N1 = (r2*z3-r3*z2+(z2-z3)*r0+(r3-r2)*z0)*(1/(2*A0e));
N2 = (r3*z1-r1*z3+(z3-z1)*r0+(r1-r3)*z0)*(1/(2*A0e));
N3 = (r1*z2-r2*z1+(z1-z2)*r0+(r2-r1)*z0)*(1/(2*A0e));


H0e = 1/(2*A0e)*[z2-z3 , 0  , z3-z1,  0  , z1-z2, 0;
      r3-r2,  0, r1-r3, 0, r2-r1, 0;
      0, z2-z3, 0,z3-z1,  0 ,z1-z2;
      0   ,r3-r2,  0 ,r1-r3, 0,r2-r1;
      2*A0e*N1/r0, 0,2*A0e*N2/r0, 0, 2*A0e*N3/r0, 0];

Ae = H0e*ed';

defgrad = [ 1+Ae(1), Ae(2), Ae(3), 1+Ae(4), 1+Ae(5)]'; % numeriskt fel på tredje raden???


end

