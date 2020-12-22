function [Ke] = cont2D3ue(ec, t, D, stress)   % Compute the element stiffness matrix in an updated Lagrangian formulation

C = [1           1        1;
    ec(1, 1) ec(1, 2) ec(1, 3);
    ec(2, 1) ec(2, 2) ec(2, 3)];

A0e = 1/2*det(C);
v0 = A0e*t;

dN1dx0 = 1/(2*A0e)*(ec(2, 2) - ec(2, 3));
dN2dx0 = 1/(2*A0e)*(ec(2, 3) - ec(2, 1));
dN3dx0 = 1/(2*A0e)*(ec(2, 1) - ec(2, 2));
dN1dy0 = 1/(2*A0e)*(ec(1, 3) - ec(1, 2));
dN2dy0 = 1/(2*A0e)*(ec(1, 1) - ec(1, 3));
dN3dy0 = 1/(2*A0e)*(ec(1, 2) - ec(1, 1));

B = [dN1dx0    0   dN2dx0    0   dN3dx0    0;        % 3x6
    0      dN1dy0   0    dN2dy0   0    dN3dy0;
    dN1dy0 dN1dx0 dN2dy0 dN2dx0 dN3dy0 dN3dx0];

H = [dN1dx0 0 dN2dx0 0 dN3dx0 0;      % 4x6
    dN1dy0 0 dN2dy0 0 dN3dy0 0;
    0 dN1dx0 0 dN2dx0 0 dN3dx0;
    0 dN1dy0 0 dN2dy0 0 dN3dy0];


S = [stress(1), stress(3);
    stress(3), stress(2)];
z = zeros(2);

R = [S z;
    z S];

Ke = (B'*D*B + H'*R*H)*v0;
end

