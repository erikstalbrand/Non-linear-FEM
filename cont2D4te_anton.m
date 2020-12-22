%function Ke=cont2D4te(ec,t,D,ed,stress)
x1=ec(1,1);
x2=ec(1,2);
x3=ec(1,3);
x4=ec(1,4);
y1=ec(2,1);
y2=ec(2,2);
y3=ec(2,3);
y4=ec(2,4);

syms eta e;

N1 = 1/4*(1-e)*(1-eta);
N2 = 1/4*(1+e)*(1-eta);
N3 = 1/4*(1+e)*(1+eta);
N4 = 1/4*(1-e)*(1+eta);

x0 = N1*x1 + N2*x2 + N3*x3 + N4*x4;
y0 = N1*y1 + N2*y2 + N3*y3 + N4*y4;

F = [x0 y0];
v = [e eta];

%Jacobian and inverse of it
J = jacobian(F, v);
Jinv = inv(J);

% Area
A = [ 1,  1,  1;
     x1, x2, x3;
     y1, y2, y3];

Ae = 1/2*det(A);

gausscoord = 0.577350269189626;
gaussmatrix = [-gausscoord -gausscoord;  % xi eta 1
                gausscoord -gausscoord;  % xi eta 2
                gausscoord  gausscoord;  % xi eta 3
               -gausscoord  gausscoord]; % xi eta 4
%---------------------------------------------------
dN1 = Jinv'*[diff(N1,e);
             diff(N1,eta)];
dN1 = subs(dN1, v, gaussmatrix(1, :));

dN2 = Jinv'*[diff(N2,e);
             diff(N2,eta)];
dN2 = subs(dN2, v, gaussmatrix(2, :));

dN3 = Jinv'*[diff(N3,e);
             diff(N3,eta)];
dN3 = subs(dN3, v, gaussmatrix(3, :));

dN4 = Jinv'*[diff(N4,e);
             diff(N4,eta)];
dN4 = subs(dN4, v, gaussmatrix(4, :));
%---------------------------------------------------

B0le = [dN1(1),  0  ,dN2(1),  0  ,dN3(1),  0  ,dN4(1),  0  ;            % 3x8
            0  ,dN1(2),  0  ,dN2(2),  0  ,dN3(2),  0  ,dN4(2);
          dN1(1),dN1(2),dN2(1),dN2(2),dN3(1),dN3(2),dN4(1),dN4(2)];

H0 = [dN1(1),  0  ,dN2(1),  0  ,dN3(1),  0  ,dN4(1),  0  ;              % 4x8
      dN1(2),  0  ,dN2(2),  0  ,dN3(2),  0  ,dN4(2),  0  ;
        0  ,dN1(1),  0  ,dN2(1),  0  ,dN3(1),  0  ,dN4(1);
        0  ,dN1(2),  0  ,dN2(2),  0  ,dN3(2),  0  ,dN4(2)];


AES = H0*ed';                           % 4x1

AE = [AES(1),   0  ,  AES(3),   0  ;       % 3x4
       0    , AES(2),  0   ,  AES(4);
      AES(2), AES(1), AES(4), AES(3)];

B0 = B0le + AE * H0;                          % 3x8 + 3x8

B0 = double(B0);
H0 = double(H0);

Ke1 = zeros(8);
for n = 1:4

    J = jacobian(F, v);
    J = subs(J, v, gaussmatrix(n, :));
    
stressvec = stress{n, 1};
    
R0 = [stressvec(1) stressvec(3)     0           0;
      stressvec(3) stressvec(2)     0           0;
          0             0       stressvec(1) stressvec(3);
          0             0       stressvec(3) stressvec(2)];
        
    Ke1 = Ke1 + (B0'*D{n}*B0 + H0'*R0*H0)*det(J);
end
% Ke1 = double(Ke1);
%end