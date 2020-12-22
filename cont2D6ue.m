function [Ke] = cont2D6ue(ec, t, D, stress)
x1 = ec(1,1);
x2 = ec(1,2);
x3 = ec(1,3);
x4 = ec(1,4);
x5 = ec(1,5);
x6 = ec(1,6);
y1 = ec(2,1);
y2 = ec(2,2);
y3 = ec(2,3);
y4 = ec(2,4);
y5 = ec(2,5);
y6 = ec(2,6);

zmatrix = 1/3*[2    1   1;
               1    2   1;
               1    1   2];

% zmatrix = [ 1   0  0;
%             0   1  0;
%             0   0  1;
%            1/2 1/2 0;
%            0 1/2 1/2;
%            1/2 0 1/2];

Ke = zeros(12);

for n = 1:3
    
z1 = zmatrix(n, 1);
z2 = zmatrix(n, 2);
z3 = zmatrix(n, 3);

% N1 = z1*(2*z1-1);
% N2 = z2*(2*z2-1);
% N3 = z3*(2*z3-1);
% N4 = 4*z1*z2;
% N5 = 4*z2*z3;
% N6 = 4*z3*z1;

% Components for Jacobian, 5.53

% With constraint z1+z2+z3 = 1
dx_dz1 = (4*z1-1)*x1 + 4*z2*x4 + 4*z3*x6;
dx_dz2 = (4*z2-1)*x2 + 4*z1*x4 + 4*z3*x5;
dx_dz3 = (4*z3-1)*x3 + 4*z2*x5 + 4*z1*x6;

dy_dz1 = (4*z1-1)*y1 + 4*z2*y4 + 4*z3*y6;
dy_dz2 = (4*z2-1)*y2 + 4*z1*y4 + 4*z3*y5;
dy_dz3 = (4*z3-1)*y3 + 4*z2*y5 + 4*z1*y6;

%Jacobian and its inverse
J = [   1      1      1;
     dx_dz1 dx_dz2 dx_dz3;
     dy_dz1 dy_dz2 dy_dz3];
Jinv = inv(J);

P = [0 1 0; 0 0 1]*Jinv';               % 2x3

% With constraint z1+z2+z3 = 1
dN = P * [4*z1-1 0 0 4*z2 0 4*z3;       % 2x6
          0 4*z2-1 0 4*z1 4*z3 0;
          0 0 4*z3-1 0 4*z2 4*z1];

%---------------------------------------------------

B0 = [dN(1,1),  0  ,dN(1,2),  0  ,dN(1,3),  0  ,dN(1,4),  0, dN(1,5), 0, dN(1,6), 0;                    % 3x12
          0  ,dN(2, 1),  0  ,dN(2, 2),  0  ,dN(2,3),  0  ,dN(2,4), 0, dN(2,5), 0, dN(2,6);
          dN(2,1),dN(1,1),dN(2,2),dN(1,2),dN(2,3),dN(1,3),dN(2,4),dN(1,4),dN(2,5),dN(1,5),dN(2,6),dN(1,6)];

H0 = [dN(1,1),  0  ,dN(1,2),  0  ,dN(1,3),  0  ,dN(1,4),  0  ,dN(1,5),  0  ,dN(1,6),  0;                % 4x12
      dN(2,1),  0  ,dN(2,2),  0  ,dN(2,3),  0  ,dN(2,4),  0  ,dN(2,5),  0  ,dN(2,6),  0;
        0,    dN(1,1),  0  ,dN(1,2),  0  ,dN(1,3),  0  ,dN(1,4),  0  ,dN(1,5),  0  ,dN(1,6);
        0,    dN(2,1),  0  ,dN(2,2),  0  ,dN(2,3),  0  ,dN(2,4),  0  ,dN(2,5),  0  ,dN(2,6)];

stressvec = stress{n};

R0 = [stressvec(1) stressvec(3)      0           0;
      stressvec(3) stressvec(2)      0           0;
          0             0       stressvec(1) stressvec(3);
          0             0       stressvec(3) stressvec(2)];
        
Ke = Ke + (B0'*D{n}*B0 + H0'*R0*H0)*1/2*det(J)*t*1/3;
end
end